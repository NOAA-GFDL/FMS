!-----------------------------------------------------------------------
!   Domain decomposition and domain update for message-passing codes
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------
#include <os.h>

!shmalloc is used on MPP SGI/Cray systems for shmem
#if defined(use_libSMA) && defined(SGICRAY_MPP)
#define use_shmalloc
#endif

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999
  use mpp_mod
  implicit none
  private
  character(len=128), private :: version= &
       '$Id: mpp_domains.F90,v 6.3 2002/02/22 19:09:12 fms Exp $'
  character(len=128), private :: name= &
       '$Name: galway $'
  character(len=128), private :: version_update_domains2D, version_global_reduce, version_global_sum, version_global_field

!parameters used to define domains: these are passed to the flags argument of mpp_define_domains
!  if data domain is to have global extent, set GLOBAL_DATA_DOMAIN
!  if global domain has periodic boundaries, set CYCLIC_GLOBAL_DOMAIN
!  sum flags together if more than one of the above conditions is to be met.
  integer, parameter, private :: GLOBAL=0, CYCLIC=1
  integer, parameter, private :: WEST=2, EAST=3, SOUTH=4, NORTH=5
  integer, parameter, private :: SEND=1, RECV=2
  integer, parameter, public :: GLOBAL_DATA_DOMAIN=2**GLOBAL, CYCLIC_GLOBAL_DOMAIN=2**CYCLIC
!gridtypes
  integer, parameter, private :: AGRID=0, BGRID=1, CGRID=2
  integer, parameter, public :: BGRID_NE=BGRID+2**NORTH+2**EAST
  integer, parameter, public :: BGRID_SW=BGRID+2**SOUTH+2**WEST
  integer, parameter, public :: CGRID_NE=CGRID+2**NORTH+2**EAST
  integer, parameter, public :: CGRID_SW=CGRID+2**SOUTH+2**WEST
  integer, private :: grid_offset_type=AGRID
!folds
  integer, parameter, public :: FOLD_WEST_EDGE = 2**WEST, FOLD_EAST_EDGE = 2**EAST
  integer, parameter, public :: FOLD_SOUTH_EDGE=2**SOUTH, FOLD_NORTH_EDGE=2**NORTH
!update
  integer, parameter, public :: WUPDATE=2**WEST, EUPDATE=2**EAST, SUPDATE=2**SOUTH, NUPDATE=2**NORTH
  integer, parameter, public :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE
!used by mpp_global_sum
  integer, parameter, public :: BITWISE_EXACT_SUM=1

  type, public :: domain_axis_spec        !type used to specify index limits along an axis of a domain
     private
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec
  type, public :: domain1D
     private
     type(domain_axis_spec) :: compute, data, global
     logical :: cyclic
     type(domain1D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type, private :: rectangle
     integer :: is, ie, js, je
     logical :: overlap, folded
  end type rectangle
  type, public :: domain2D
     private
     type(domain1D) :: x
     type(domain1D) :: y
     type(domain2D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer :: fold, gridtype
     logical :: overlap
     type(rectangle) :: recv_e, recv_se, recv_s, recv_sw, &
                        recv_w, recv_nw, recv_n, recv_ne
     type(rectangle) :: send_e, send_se, send_s, send_sw, &
                        send_w, send_nw, send_n, send_ne
     logical :: remote_domains_initialized
     type(rectangle) :: recv_e_off, recv_se_off, recv_s_off, recv_sw_off, &
                        recv_w_off, recv_nw_off, recv_n_off, recv_ne_off
     type(rectangle) :: send_e_off, send_se_off, send_s_off, send_sw_off, &
                        send_w_off, send_nw_off, send_n_off, send_ne_off
     logical :: remote_off_domains_initialized
  end type domain2D
  type(domain1D), public :: NULL_DOMAIN1D
  type(domain2D), public :: NULL_DOMAIN2D

  integer, private :: pe, npes

  integer, private :: tk
  logical, private :: verbose=.FALSE., debug=.FALSE., domain_clocks_on=.FALSE.
  logical, private :: mpp_domains_initialized=.FALSE.
  integer, parameter, public :: MPP_DOMAIN_TIME=MPP_DEBUG+1
  integer :: send_clock=0, recv_clock=0, unpk_clock=0, wait_clock=0, pack_clock=0, pack_loop_clock=0

!stack for internal buffers
!allocated differently if use_shmalloc
#ifdef use_shmalloc
  real(DOUBLE_KIND), private :: mpp_domains_stack(1)
  pointer( ptr_stack, mpp_domains_stack )
#else
  real(DOUBLE_KIND), private, allocatable :: mpp_domains_stack(:)
#endif
  integer, private :: mpp_domains_stack_size=0, mpp_domains_stack_hwm=0

!used by mpp_define_domains2D_new to transmit data
  integer :: domain_info_buf(16)
#ifdef use_shmalloc
  pointer( ptr_info, domain_info_buf )
#endif

!public interfaces
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
     module procedure mpp_update_domain2D_r8_2dv
     module procedure mpp_update_domain2D_r8_3dv
     module procedure mpp_update_domain2D_r8_4dv
     module procedure mpp_update_domain2D_r8_5dv
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_i8_5d
     module procedure mpp_update_domain2D_l8_2d
     module procedure mpp_update_domain2D_l8_3d
     module procedure mpp_update_domain2D_l8_4d
     module procedure mpp_update_domain2D_l8_5d
#endif
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
     module procedure mpp_update_domain2D_l4_2d
     module procedure mpp_update_domain2D_l4_3d
     module procedure mpp_update_domain2D_l4_4d
     module procedure mpp_update_domain2D_l4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv

!     module procedure mpp_update_domain1D_r8_2d
!     module procedure mpp_update_domain1D_r8_3d
!     module procedure mpp_update_domain1D_r8_4d
!     module procedure mpp_update_domain1D_r8_5d
!     module procedure mpp_update_domain1D_c8_2d
!     module procedure mpp_update_domain1D_c8_3d
!     module procedure mpp_update_domain1D_c8_4d
!     module procedure mpp_update_domain1D_c8_5d
!#ifndef no_8byte_integers
!     module procedure mpp_update_domain1D_i8_2d
!     module procedure mpp_update_domain1D_i8_3d
!     module procedure mpp_update_domain1D_i8_4d
!     module procedure mpp_update_domain1D_i8_5d
!     module procedure mpp_update_domain1D_l8_2d
!     module procedure mpp_update_domain1D_l8_3d
!     module procedure mpp_update_domain1D_l8_4d
!     module procedure mpp_update_domain1D_l8_5d
!#endif
!     module procedure mpp_update_domain1D_r4_2d
!     module procedure mpp_update_domain1D_r4_3d
!     module procedure mpp_update_domain1D_r4_4d
!     module procedure mpp_update_domain1D_r4_5d
!     module procedure mpp_update_domain1D_c4_2d
!     module procedure mpp_update_domain1D_c4_3d
!     module procedure mpp_update_domain1D_c4_4d
!     module procedure mpp_update_domain1D_c4_5d
!     module procedure mpp_update_domain1D_i4_2d
!     module procedure mpp_update_domain1D_i4_3d
!     module procedure mpp_update_domain1D_i4_4d
!     module procedure mpp_update_domain1D_i4_5d
!     module procedure mpp_update_domain1D_l4_2d
!     module procedure mpp_update_domain1D_l4_3d
!     module procedure mpp_update_domain1D_l4_4d
!     module procedure mpp_update_domain1D_l4_5d
  end interface

  interface mpp_redistribute
     module procedure mpp_redistribute_r8_2D
     module procedure mpp_redistribute_r8_3D
     module procedure mpp_redistribute_r8_4D
     module procedure mpp_redistribute_r8_5D
     module procedure mpp_redistribute_c8_2D
     module procedure mpp_redistribute_c8_3D
     module procedure mpp_redistribute_c8_4D
     module procedure mpp_redistribute_c8_5D
#ifndef no_8byte_integers
     module procedure mpp_redistribute_i8_2D
     module procedure mpp_redistribute_i8_3D
     module procedure mpp_redistribute_i8_4D
     module procedure mpp_redistribute_i8_5D
     module procedure mpp_redistribute_l8_2D
     module procedure mpp_redistribute_l8_3D
     module procedure mpp_redistribute_l8_4D
     module procedure mpp_redistribute_l8_5D
#endif
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
     module procedure mpp_redistribute_l4_2D
     module procedure mpp_redistribute_l4_3D
     module procedure mpp_redistribute_l4_4D
     module procedure mpp_redistribute_l4_5D
  end interface

  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
#endif
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
     module procedure mpp_global_field2D_c4_2d
     module procedure mpp_global_field2D_c4_3d
     module procedure mpp_global_field2D_c4_4d
     module procedure mpp_global_field2D_c4_5d
     module procedure mpp_global_field2D_i4_2d
     module procedure mpp_global_field2D_i4_3d
     module procedure mpp_global_field2D_i4_4d
     module procedure mpp_global_field2D_i4_5d
     module procedure mpp_global_field2D_l4_2d
     module procedure mpp_global_field2D_l4_3d
     module procedure mpp_global_field2D_l4_4d
     module procedure mpp_global_field2D_l4_5d

     module procedure mpp_global_field1D_r8_2d
     module procedure mpp_global_field1D_c8_2d
#ifndef no_8byte_integers
     module procedure mpp_global_field1D_i8_2d
     module procedure mpp_global_field1D_l8_2d
#endif
     module procedure mpp_global_field1D_r4_2d
     module procedure mpp_global_field1D_c4_2d
     module procedure mpp_global_field1D_i4_2d
     module procedure mpp_global_field1D_l4_2d
  end interface

  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
#ifndef no_8byte_integers
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
#endif
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
#ifndef no_8byte_integers
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
#endif
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#ifndef no_8byte_integers
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
#endif
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface

  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface

  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface

  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface

  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface

  public :: mpp_define_layout, mpp_define_domains, mpp_domains_init, mpp_domains_set_stack_size, mpp_domains_exit,&
            mpp_get_compute_domain, mpp_get_compute_domains, mpp_get_data_domain, mpp_get_global_domain, &
            mpp_get_domain_components, mpp_get_layout, mpp_get_pelist, mpp_redistribute, mpp_update_domains, &
            mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum, operator(.EQ.), operator(.NE.)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                 MPP_DOMAINS: initialization and termination                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_domains_init(flags)
!initialize domain decomp package
      integer, intent(in), optional :: flags
      integer :: l=0

      if( mpp_domains_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      pe = mpp_pe()
      npes = mpp_npes()
      mpp_domains_initialized = .TRUE.
      if( pe.EQ.mpp_root_pe() )then
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(version)
!          write( stdlog(),'(a)' )trim(version_update_domains2D)
      end if

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
          domain_clocks_on = flags.EQ.MPP_DOMAIN_TIME
      end if

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary
#ifdef use_shmalloc
      call mpp_malloc( ptr_info, 16, l )
#endif

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      NULL_DOMAIN1D%global%begin  = -1; NULL_DOMAIN1D%global%end  = -1; NULL_DOMAIN1D%global%size = 0
      NULL_DOMAIN1D%data%begin    = -1; NULL_DOMAIN1D%data%end    = -1; NULL_DOMAIN1D%data%size = 0
      NULL_DOMAIN1D%compute%begin = -1; NULL_DOMAIN1D%compute%end = -1; NULL_DOMAIN1D%compute%size = 0
      NULL_DOMAIN1D%pe = NULL_PE
      NULL_DOMAIN2D%x = NULL_DOMAIN1D
      NULL_DOMAIN2D%y = NULL_DOMAIN1D
      NULL_DOMAIN2D%pe = NULL_PE

      if( domain_clocks_on )then
          pack_clock = mpp_clock_id( 'Halo pack' )
          pack_loop_clock = mpp_clock_id( 'Halo pack loop' )
          send_clock = mpp_clock_id( 'Halo send' )
          recv_clock = mpp_clock_id( 'Halo recv' )
          unpk_clock = mpp_clock_id( 'Halo unpk' )
          wait_clock = mpp_clock_id( 'Halo wait' )
      end if
      return
    end subroutine mpp_domains_init

    subroutine mpp_domains_set_stack_size(n)
!set the mpp_domains_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=8) :: text

      if( n.LE.mpp_domains_stack_size )return
#ifdef use_shmalloc
      call mpp_malloc( ptr_stack, n, mpp_domains_stack_size )
#else
      if( allocated(mpp_domains_stack) )deallocate(mpp_domains_stack)
      allocate( mpp_domains_stack(n) )
      mpp_domains_stack_size = n
#endif
      write( text,'(i8)' )n
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_DOMAINS_SET_STACK_SIZE: stack size set to '//text//'.' )

      return
    end subroutine mpp_domains_set_stack_size

    subroutine mpp_domains_exit()
!currently does not have much to do, but provides the possibility of re-initialization
      if( .NOT.mpp_domains_initialized )return
      call mpp_max(mpp_domains_stack_hwm)
      if( pe.EQ.mpp_root_pe() )write( stdout(),* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      mpp_domains_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                MPP_DOMAINS: overloaded operators (==, /=)                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function mpp_domain1D_eq( a, b )
      logical :: mpp_domain1D_eq
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_eq = .FALSE.
      if( a%compute%begin.EQ.b%compute%begin .AND. &
          a%compute%end  .EQ.b%compute%end   .AND. &
          a%data%begin   .EQ.b%data%begin    .AND. &
          a%data%end     .EQ.b%data%end      .AND. & 
          a%global%begin .EQ.b%global%begin  .AND. &
          a%global%end   .EQ.b%global%end )mpp_domain1D_eq = .TRUE.
      return
    end function mpp_domain1D_eq

    function mpp_domain1D_ne( a, b )
      logical :: mpp_domain1D_ne
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_ne = .NOT. ( a.EQ.b )
      return
    end function mpp_domain1D_ne

    function mpp_domain2D_eq( a, b )
      logical :: mpp_domain2D_eq
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_eq = a%x.EQ.b%x .AND. a%y.EQ.b%y
      return
    end function mpp_domain2D_eq

    function mpp_domain2D_ne( a, b )
      logical :: mpp_domain2D_ne
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_ne = .NOT. ( a.EQ.b )
      return
    end function mpp_domain2D_ne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_DEFINE_DOMAINS: define layout and decomposition            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap )
!routine to divide global array indices among domains, and assign domains to PEs
!domain is of type domain1D
!ARGUMENTS:
!      global_indices(2)=(isg,ieg) gives the extent of global domain
!      ndivs is number of divisions of domain: even divisions unless extent is present.
!      domain is the returned domain1D
!      pelist (optional) list of PEs to which domains are to be assigned (default 0...npes-1)
!                 size of pelist must correspond to number of mask=.TRUE. divisions
!      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
!      halo (optional) defines halo width (currently the same on both sides)
!      extent (optional) array defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
!      maskmap (optional) a division whose maskmap=.FALSE. is not assigned to any domain
!  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
!  as the input parameters permit
      integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
      integer, intent(in) :: ndivs
      type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: flags, halo
      integer, intent(in), optional :: extent(0:)
      logical, intent(in), optional :: maskmap(0:)

      logical :: compute_domain_is_global, data_domain_is_global
      integer :: ndiv, n, isg, ieg, is, ie, i, off, pos, hs, he
      integer, allocatable :: pes(:)
      logical, allocatable :: mask(:)
      integer :: halosz
!used by symmetry algorithm
      integer :: imax, ndmax, ndmirror
      logical :: symmetrize
!statement functions
      logical :: even, odd
      even(n) = (mod(n,2).EQ.0)
      odd (n) = (mod(n,2).EQ.1)
      
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
!get global indices
      isg = global_indices(1)
      ieg = global_indices(2)
      if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
!get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
      if( PRESENT(pelist) )then
          if( .NOT.any(pelist.EQ.pe) )then
              write( stderr(),* )'pe=', pe, ' pelist=', pelist
              call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
          end if
          allocate( pes(0:size(pelist)-1) )
          pes(:) = pelist(:)
      else
          allocate( pes(0:npes-1) )
          pes(:) = (/ (i,i=0,npes-1) /)
      end if

!get number of real domains: 1 mask domain per PE in pes
      allocate( mask(0:ndivs-1) )
      mask = .TRUE.                 !default mask
      if( PRESENT(maskmap) )then
          if( size(maskmap).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: maskmap array size must equal number of domain divisions.' )
          mask(:) = maskmap(:)
      end if
      if( count(mask).NE.size(pes) ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: number of TRUEs in maskmap array must match PE count.' )
      if( PRESENT(extent) )then
          if( size(extent).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
      end if
!get halosize
      halosz = 0
      if( PRESENT(halo) )halosz = halo

!get flags
      compute_domain_is_global = .FALSE.
      data_domain_is_global    = .FALSE.
      domain%cyclic = .FALSE.
      if( PRESENT(flags) )then
!NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
          compute_domain_is_global = ndivs.EQ.1
!if compute domain is global, data domain must also be
          data_domain_is_global    = BTEST(flags,GLOBAL) .OR. compute_domain_is_global
          domain%cyclic  = BTEST(flags,CYCLIC) .AND. halosz.NE.0
      end if

!set up links list
      allocate( domain%list(0:ndivs-1) )

!set global domain
      domain%list(:)%global%begin     = isg
      domain%list(:)%global%end       = ieg
      domain%list(:)%global%size      = ieg-isg+1
      domain%list(:)%global%max_size  = ieg-isg+1
      domain%list(:)%global%is_global = .TRUE. !always

!get compute domain
      if( compute_domain_is_global )then
          domain%list(:)%compute%begin = isg
          domain%list(:)%compute%end   = ieg
          domain%list(:)%compute%is_global = .TRUE.
          domain%list(:)%pe = pes(:)
          domain%pos = 0
      else
          domain%list(:)%compute%is_global = .FALSE.
          is = isg
          n = 0
          do ndiv=0,ndivs-1
             if( PRESENT(extent) )then
                 ie = is + extent(ndiv) - 1
                 if( ndiv.EQ.ndivs-1 .AND. ie.NE.ieg ) &
                      call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extent array limits do not match global domain.' )
             else
!modified for mirror-symmetry
!original line
!                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

!problem of dividing nx points into n domains maintaining symmetry
!i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
!this will always work for nx even n even or odd
!this will always work for nx odd, n odd
!this will never  work for nx odd, n even: for this case we supersede the mirror calculation
!                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
!nx even n odd fails if n>nx/2
                 symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                              (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                              (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

!mirror domains are stored in the list and retrieved if required.
                 if( ndiv.EQ.0 )then
!initialize max points and max domains
                     imax = ieg
                     ndmax = ndivs
                 end if
!do bottom half of decomposition, going over the midpoint for odd ndivs
                 if( ndiv.LT.(ndivs-1)/2+1 )then
!domain is sized by dividing remaining points by remaining domains
                     ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                     ndmirror = (ndivs-1) - ndiv !mirror domain
                     if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
!mirror extents, the max(,) is to eliminate overlaps
                         domain%list(ndmirror)%compute%begin = max( isg+ieg-ie, ie+1 )
                         domain%list(ndmirror)%compute%end   = max( isg+ieg-is, ie+1 )
                         imax = domain%list(ndmirror)%compute%begin - 1
                         ndmax = ndmax - 1
                     end if
                 else
                     if( symmetrize )then
!do top half of decomposition by retrieving saved values
                         is = domain%list(ndiv)%compute%begin
                         ie = domain%list(ndiv)%compute%end
                     else
                         ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                     end if
                 end if
             end if
             if( ie.LT.is )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents must be positive definite.' )
             domain%list(ndiv)%compute%begin = is
             domain%list(ndiv)%compute%end   = ie
             if( ndiv.GT.0 .AND. is.NE.domain%list(ndiv-1)%compute%end+1 ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
             if( ndiv.EQ.ndivs-1 .AND. domain%list(ndiv)%compute%end.NE.ieg ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
             if( mask(ndiv) )then
                 domain%list(ndiv)%pe = pes(n)
                 if( pe.EQ.pes(n) )domain%pos = ndiv
                 n = n + 1
             end if
             is = ie + 1
          end do
      end if
          
      domain%list(:)%compute%size  = domain%list(:)%compute%end - domain%list(:)%compute%begin + 1

!get data domain
!data domain is at least equal to compute domain
      domain%list(:)%data%begin = domain%list(:)%compute%begin
      domain%list(:)%data%end   = domain%list(:)%compute%end
      domain%list(:)%data%is_global = .FALSE.
!apply global flags
      if( data_domain_is_global )then
          domain%list(:)%data%begin  = isg
          domain%list(:)%data%end    = ieg
          domain%list(:)%data%is_global = .TRUE.
      end if
!apply margins
      domain%list(:)%data%begin = domain%list(:)%data%begin - halosz
      domain%list(:)%data%end   = domain%list(:)%data%end   + halosz  
      domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

!      domain = domain%list(pos) !load domain from domain%list(pos)
      domain%compute = domain%list(domain%pos)%compute
      domain%data = domain%list(domain%pos)%data
      domain%global = domain%list(domain%pos)%global
      domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
      domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
      domain%global%max_size  = domain%global%size

!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask )
      return

      contains

        function if_overlap( hs, he, cs, ce, os, oe )
!function to look if a portion of the halo [hs,he] lies with in the compute region [cs,ce]
!if yes, if_overlap returns true, and the overlap region is returned in [os,oe]
          logical :: if_overlap
          integer, intent(in) :: hs, he, cs, ce
          integer, intent(out) :: os, oe
          os = max(hs,cs)
          oe = min(he,ce)
          if( debug )write( stderr(),'(a,7i4)' ) &
               'MPP_DEFINE_DOMAINS1D: pe, hs, he, cs, ce, os, oe=', pe, hs, he, cs, ce, os, oe
          if_overlap = oe.GE.os
          return
        end function if_overlap
          
    end subroutine mpp_define_domains1D

    subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, &
                                                     xflags, yflags, xhalo, yhalo, xextent, yextent, maskmap, name )
!define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      integer, intent(in) :: layout(2)
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
      integer, intent(in), optional :: xextent(0:), yextent(0:)
      logical, intent(in), optional :: maskmap(0:,0:)
      character(len=*), optional :: name
      integer :: i, j, m, n
      integer :: ipos, jpos, pos
      integer :: ndivx, ndivy, isg, ieg, jsg, jeg, isd, ied, jsd, jed
      
      logical, allocatable :: mask(:,:)
      integer, allocatable :: pes(:), pearray(:,:)
      character(len=8) :: text

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
      ndivx = layout(1); ndivy = layout(2)
      isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

      if( PRESENT(pelist) )then
          allocate( pes(0:size(pelist)-1) )
          pes = pelist
      else
          allocate( pes(0:npes-1) )
          pes = (/ (i,i=0,npes-1) /)
      end if

      allocate( mask(0:ndivx-1,0:ndivy-1) )
      mask = .TRUE.
      if( PRESENT(maskmap) )then
          if( size(maskmap,1).NE.ndivx .OR. size(maskmap,2).NE.ndivy ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: maskmap array does not match layout.' )
          mask(:,:) = maskmap(:,:)
      end if
!number of unmask domains in layout must equal number of PEs assigned
      n = count(mask)
      if( n.NE.size(pes) )then
          write( text,'(i8)' )n
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for this layout and maskmap. Use ' &
               //text//' PEs for this domain decomposition.' )
      end if

!place on PE array; need flag to assign them to j first and then i
      allocate( pearray(0:ndivx-1,0:ndivy-1) )
      pearray(:,:) = NULL_PE
      ipos = -1; jpos = -1; pos = -1
      n = 0
      do j = 0,ndivy-1
         do i = 0,ndivx-1
            if( mask(i,j) )then
                pearray(i,j) = pes(n)
                if( pes(n).EQ.pe )then
                    pos = n
                    ipos = i
                    jpos = j
                end if
                n = n + 1
            end if
         end do
      end do
      if( ipos.EQ.-1 .OR. jpos.EQ.-1 .or. pos.EQ.-1 ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: pelist must include this PE.' )
      if( debug )write( stderr(), * )'pe, ipos, jpos=', pe, ipos, jpos, ' pearray(:,jpos)=', pearray(:,jpos), &
                                                                        ' pearray(ipos,:)=', pearray(ipos,:)
      
!do domain decomposition using 1D versions in X and Y
      call mpp_define_domains( global_indices(1:2), ndivx, domain%x, &
           pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos) )
      call mpp_define_domains( global_indices(3:4), ndivy, domain%y, &
           pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:) )
      if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
      domain%pos = pos

!set up fold
      domain%fold = 0
      if( PRESENT(xflags) )then
          if( BTEST(xflags,WEST) )domain%fold = domain%fold + FOLD_WEST_EDGE
          if( BTEST(xflags,EAST) )domain%fold = domain%fold + FOLD_EAST_EDGE
      end if
      if( PRESENT(yflags) )then
          if( BTEST(yflags,SOUTH) )domain%fold = domain%fold + FOLD_SOUTH_EDGE
          if( BTEST(yflags,NORTH) )domain%fold = domain%fold + FOLD_NORTH_EDGE
      end if
          
      if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,NORTH) )then
          if( domain%y%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
!check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x%list(i)%compute%size.NE.domain%x%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
      end if
      if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
          if( domain%x%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
!check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y%list(i)%compute%size.NE.domain%y%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
      end if

!set up domain%list
      if( debug )write( stderr(),'(a,9i4)' )'pe, domain=', pe, domain_info_buf(1:8)
      if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )then
          write( stdlog(), '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( stdlog(), '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
      end if
      call mpp_sync()
      call mpp_get_compute_domain( domain, domain_info_buf(1), domain_info_buf(2), domain_info_buf(3), domain_info_buf(4) )
      call mpp_get_data_domain   ( domain, domain_info_buf(5), domain_info_buf(6), domain_info_buf(7), domain_info_buf(8) )
      n = size(pes)
      allocate( domain%list(0:n-1) ) !this is only used for storage of remote compute and data domain info
      do i = 0,n-1
         m = mod(pos+i,n)
         domain%list(m)%pe = pes(m)
         call mpp_transmit( domain_info_buf(1), 8, pes(mod(pos+n-i,n)), domain_info_buf(9), 8, pes(m) )
         domain%list(m)%x%compute%begin = domain_info_buf(9)
         domain%list(m)%x%compute%end   = domain_info_buf(10)
         domain%list(m)%y%compute%begin = domain_info_buf(11)
         domain%list(m)%y%compute%end   = domain_info_buf(12)
         domain%list(m)%x%data%begin = domain_info_buf(13)
         domain%list(m)%x%data%end   = domain_info_buf(14)
         domain%list(m)%y%data%begin = domain_info_buf(15)
         domain%list(m)%y%data%end   = domain_info_buf(16)
         if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )write( stdlog(), '(2x,i3,x,4i5,3x,4i5)' )pes(m), domain_info_buf(9:)
      end do
      call mpp_sync_self(pes)
      domain%list(:)%x%compute%size = domain%list(:)%x%compute%end - domain%list(:)%x%compute%begin + 1
      domain%list(:)%y%compute%size = domain%list(:)%y%compute%end - domain%list(:)%y%compute%begin + 1
      domain%list(:)%x%data%size = domain%list(:)%x%data%end - domain%list(:)%x%data%begin + 1
      domain%list(:)%y%data%size = domain%list(:)%y%data%end - domain%list(:)%y%data%begin + 1

      domain%remote_domains_initialized = .FALSE.
      domain%remote_off_domains_initialized = .FALSE.
      call compute_overlaps(domain)
!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask, pearray )
          
      return
    end subroutine mpp_define_domains2D

    subroutine compute_overlaps( domain )
!computes remote domain overlaps
!assumes only one in each direction
      type(domain2D), intent(inout) :: domain
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff
      integer :: list

      if( grid_offset_type.EQ.AGRID .AND. domain%remote_domains_initialized     )return
      if( grid_offset_type.NE.AGRID .AND. domain%remote_off_domains_initialized )return
      domain%gridtype = grid_offset_type
      n = size(domain%list)
!send
      call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ioff, ysize=joff ) !cyclic offsets
      domain%list(:)%overlap = .FALSE.
      do list = 0,n-1
         m = mod( domain%pos+list, n )
!to_pe's eastern halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_w_off%overlap = .TRUE.
                 domain%list(m)%send_w_off%is = is
                 domain%list(m)%send_w_off%ie = ie
                 domain%list(m)%send_w_off%js = js
                 domain%list(m)%send_w_off%je = je
             else
                 domain%list(m)%send_w%overlap = .TRUE.
                 domain%list(m)%send_w%is = is
                 domain%list(m)%send_w%ie = ie
                 domain%list(m)%send_w%js = js
                 domain%list(m)%send_w%je = je
             end if
         else
             domain%list(m)%send_w%overlap = .FALSE.
             domain%list(m)%send_w_off%overlap = .FALSE.
         end if
!to_pe's SE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_nw_off%overlap = .TRUE.
                 domain%list(m)%send_nw_off%is = is
                 domain%list(m)%send_nw_off%ie = ie
                 domain%list(m)%send_nw_off%js = js
                 domain%list(m)%send_nw_off%je = je
             else
                 domain%list(m)%send_nw%overlap = .TRUE.
                 domain%list(m)%send_nw%is = is
                 domain%list(m)%send_nw%ie = ie
                 domain%list(m)%send_nw%js = js
                 domain%list(m)%send_nw%je = je
             end if
         else
             domain%list(m)%send_nw%overlap = .FALSE.
             domain%list(m)%send_nw_off%overlap = .FALSE.
         end if
!to_pe's southern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_n_off%overlap = .TRUE.
                 domain%list(m)%send_n_off%is = is
                 domain%list(m)%send_n_off%ie = ie
                 domain%list(m)%send_n_off%js = js
                 domain%list(m)%send_n_off%je = je
             else
                 domain%list(m)%send_n%overlap = .TRUE.
                 domain%list(m)%send_n%is = is
                 domain%list(m)%send_n%ie = ie
                 domain%list(m)%send_n%js = js
                 domain%list(m)%send_n%je = je
             end if
         else
             domain%list(m)%send_n%overlap = .FALSE.
             domain%list(m)%send_n_off%overlap = .FALSE.
         end if
!to_pe's SW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_ne_off%overlap = .TRUE.
                 domain%list(m)%send_ne_off%is = is
                 domain%list(m)%send_ne_off%ie = ie
                 domain%list(m)%send_ne_off%js = js
                 domain%list(m)%send_ne_off%je = je
             else
                 domain%list(m)%send_ne%overlap = .TRUE.
                 domain%list(m)%send_ne%is = is
                 domain%list(m)%send_ne%ie = ie
                 domain%list(m)%send_ne%js = js
                 domain%list(m)%send_ne%je = je
             end if
         else
             domain%list(m)%send_ne%overlap = .FALSE.
             domain%list(m)%send_ne_off%overlap = .FALSE.
         end if
!to_pe's western halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_e_off%overlap = .TRUE.
                 domain%list(m)%send_e_off%is = is
                 domain%list(m)%send_e_off%ie = ie
                 domain%list(m)%send_e_off%js = js
                 domain%list(m)%send_e_off%je = je
             else
                 domain%list(m)%send_e%overlap = .TRUE.
                 domain%list(m)%send_e%is = is
                 domain%list(m)%send_e%ie = ie
                 domain%list(m)%send_e%js = js
                 domain%list(m)%send_e%je = je
             end if
         else
             domain%list(m)%send_e%overlap = .FALSE.
             domain%list(m)%send_e_off%overlap = .FALSE.
         end if
!to_pe's NW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_se_off%overlap = .TRUE.
                 domain%list(m)%send_se_off%is = is
                 domain%list(m)%send_se_off%ie = ie
                 domain%list(m)%send_se_off%js = js
                 domain%list(m)%send_se_off%je = je
             else
                 domain%list(m)%send_se%overlap = .TRUE.
                 domain%list(m)%send_se%is = is
                 domain%list(m)%send_se%ie = ie
                 domain%list(m)%send_se%js = js
                 domain%list(m)%send_se%je = je
             end if
         else
             domain%list(m)%send_se%overlap = .FALSE.
             domain%list(m)%send_se_off%overlap = .FALSE.
         end if
!to_pe's northern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_s_off%overlap = .TRUE.
                 domain%list(m)%send_s_off%is = is
                 domain%list(m)%send_s_off%ie = ie
                 domain%list(m)%send_s_off%js = js
                 domain%list(m)%send_s_off%je = je
             else
                 domain%list(m)%send_s%overlap = .TRUE.
                 domain%list(m)%send_s%is = is
                 domain%list(m)%send_s%ie = ie
                 domain%list(m)%send_s%js = js
                 domain%list(m)%send_s%je = je
             end if
         else
             domain%list(m)%send_s%overlap = .FALSE.
             domain%list(m)%send_s_off%overlap = .FALSE.
         end if
!to_pe's NE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if
         end if
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_sw_off%overlap = .TRUE.
                 domain%list(m)%send_sw_off%is = is
                 domain%list(m)%send_sw_off%ie = ie
                 domain%list(m)%send_sw_off%js = js
                 domain%list(m)%send_sw_off%je = je
             else
                 domain%list(m)%send_sw%overlap = .TRUE.
                 domain%list(m)%send_sw%is = is
                 domain%list(m)%send_sw%ie = ie
                 domain%list(m)%send_sw%js = js
                 domain%list(m)%send_sw%je = je
             end if
         else
             domain%list(m)%send_sw%overlap = .FALSE.
             domain%list(m)%send_sw_off%overlap = .FALSE.
         end if
      end do
             
!recv
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
!recv_e
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_e%folded = .FALSE.
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 domain%list(m)%recv_e%folded = .TRUE.
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_e_off%overlap = .TRUE.
                 domain%list(m)%recv_e_off%is = is
                 domain%list(m)%recv_e_off%ie = ie
                 domain%list(m)%recv_e_off%js = js
                 domain%list(m)%recv_e_off%je = je
             else
                 domain%list(m)%recv_e%overlap = .TRUE.
                 domain%list(m)%recv_e%is = is
                 domain%list(m)%recv_e%ie = ie
                 domain%list(m)%recv_e%js = js
                 domain%list(m)%recv_e%je = je
             endif
         else
             domain%list(m)%recv_e%overlap = .FALSE.
             domain%list(m)%recv_e_off%overlap = .FALSE.
         end if
!recv_se
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_se%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_se_off%overlap = .TRUE.
                 domain%list(m)%recv_se_off%is = is
                 domain%list(m)%recv_se_off%ie = ie
                 domain%list(m)%recv_se_off%js = js
                 domain%list(m)%recv_se_off%je = je
             else
                 domain%list(m)%recv_se%overlap = .TRUE.
                 domain%list(m)%recv_se%is = is
                 domain%list(m)%recv_se%ie = ie
                 domain%list(m)%recv_se%js = js
                 domain%list(m)%recv_se%je = je
             endif
         else
             domain%list(m)%recv_se%overlap = .FALSE.
             domain%list(m)%recv_se_off%overlap = .FALSE.
         end if
!recv_s
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_s%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_s%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_s_off%overlap = .TRUE.
                 domain%list(m)%recv_s_off%is = is
                 domain%list(m)%recv_s_off%ie = ie
                 domain%list(m)%recv_s_off%js = js
                 domain%list(m)%recv_s_off%je = je
             else
                 domain%list(m)%recv_s%overlap = .TRUE.
                 domain%list(m)%recv_s%is = is
                 domain%list(m)%recv_s%ie = ie
                 domain%list(m)%recv_s%js = js
                 domain%list(m)%recv_s%je = je
             endif
         else
             domain%list(m)%recv_s%overlap = .FALSE.
             domain%list(m)%recv_s_off%overlap = .FALSE.

         end if
!recv_sw
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_sw%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_sw_off%overlap = .TRUE.
                 domain%list(m)%recv_sw_off%is = is
                 domain%list(m)%recv_sw_off%ie = ie
                 domain%list(m)%recv_sw_off%js = js
                 domain%list(m)%recv_sw_off%je = je
             else
                 domain%list(m)%recv_sw%overlap = .TRUE.
                 domain%list(m)%recv_sw%is = is
                 domain%list(m)%recv_sw%ie = ie
                 domain%list(m)%recv_sw%js = js
                 domain%list(m)%recv_sw%je = je
             endif
         else
             domain%list(m)%recv_sw%overlap = .FALSE.
             domain%list(m)%recv_sw_off%overlap = .FALSE.
         end if
!recv_w
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_w%folded = .FALSE.
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_w%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_w_off%overlap = .TRUE.
                 domain%list(m)%recv_w_off%is = is
                 domain%list(m)%recv_w_off%ie = ie
                 domain%list(m)%recv_w_off%js = js
                 domain%list(m)%recv_w_off%je = je
             else
                 domain%list(m)%recv_w%overlap = .TRUE.
                 domain%list(m)%recv_w%is = is
                 domain%list(m)%recv_w%ie = ie
                 domain%list(m)%recv_w%js = js
                 domain%list(m)%recv_w%je = je
             endif
         else
             domain%list(m)%recv_w%overlap = .FALSE.
             domain%list(m)%recv_w_off%overlap = .FALSE.
         end if
!recv_nw
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_nw%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_nw_off%overlap = .TRUE.
                 domain%list(m)%recv_nw_off%is = is
                 domain%list(m)%recv_nw_off%ie = ie
                 domain%list(m)%recv_nw_off%js = js
                 domain%list(m)%recv_nw_off%je = je
             else
                 domain%list(m)%recv_nw%overlap = .TRUE.
                 domain%list(m)%recv_nw%is = is
                 domain%list(m)%recv_nw%ie = ie
                 domain%list(m)%recv_nw%js = js
                 domain%list(m)%recv_nw%je = je
             endif
         else
             domain%list(m)%recv_nw%overlap = .FALSE.
             domain%list(m)%recv_nw_off%overlap = .FALSE.
         end if
!recv_n
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_n%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_n%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_n_off%overlap = .TRUE.
                 domain%list(m)%recv_n_off%is = is
                 domain%list(m)%recv_n_off%ie = ie
                 domain%list(m)%recv_n_off%js = js
                 domain%list(m)%recv_n_off%je = je
             else
                 domain%list(m)%recv_n%overlap = .TRUE.
                 domain%list(m)%recv_n%is = is
                 domain%list(m)%recv_n%ie = ie
                 domain%list(m)%recv_n%js = js
                 domain%list(m)%recv_n%je = je
             end if
         else
             domain%list(m)%recv_n%overlap = .FALSE.
             domain%list(m)%recv_n_off%overlap = .FALSE.
         end if
!recv_ne
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_ne%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_ne_off%overlap = .TRUE.
                 domain%list(m)%recv_ne_off%is = is
                 domain%list(m)%recv_ne_off%ie = ie
                 domain%list(m)%recv_ne_off%js = js
                 domain%list(m)%recv_ne_off%je = je
             else
                 domain%list(m)%recv_ne%overlap = .TRUE.
                 domain%list(m)%recv_ne%is = is
                 domain%list(m)%recv_ne%ie = ie
                 domain%list(m)%recv_ne%js = js
                 domain%list(m)%recv_ne%je = je
             end if
         else
             domain%list(m)%recv_ne%overlap = .FALSE.
             domain%list(m)%recv_ne_off%overlap = .FALSE.
         end if
      end do
      if( grid_offset_type.EQ.AGRID )domain%remote_domains_initialized = .TRUE.
      if( grid_offset_type.NE.AGRID )domain%remote_off_domains_initialized = .TRUE.
      return
    end subroutine compute_overlaps

    subroutine mpp_define_layout2D( global_indices, ndivs, layout )
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      integer, intent(in) :: ndivs !number of divisions to divide global domain
      integer, intent(out) :: layout(2)

      integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

      isg = global_indices(1)
      ieg = global_indices(2)
      jsg = global_indices(3)
      jeg = global_indices(4)

      isz = ieg - isg + 1
      jsz = jeg - jsg + 1
!first try to divide ndivs in the domain aspect ratio: if imperfect aspect, reduce idiv till it divides ndivs
      idiv = nint( sqrt(float(ndivs*isz)/jsz) )
      idiv = max(idiv,1) !for isz=1 line above can give 0
      do while( mod(ndivs,idiv).NE.0 )
         idiv = idiv - 1
      end do                 !will terminate at idiv=1 if not before
      jdiv = ndivs/idiv

      layout = (/ idiv, jdiv /)
      return
    end subroutine mpp_define_layout2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!     MPP_GET and SET routiness: retrieve various components of domains       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_get_compute_domain1D( domain, begin, end, size, max_size, is_global )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size
      logical, intent(out), optional :: is_global

      if( PRESENT(begin)     )begin     = domain%compute%begin
      if( PRESENT(end)       )end       = domain%compute%end
      if( PRESENT(size)      )size      = domain%compute%size
      if( PRESENT(max_size)  )max_size  = domain%compute%max_size
      if( PRESENT(is_global) )is_global = domain%compute%is_global
      return
    end subroutine mpp_get_compute_domain1D

    subroutine mpp_get_data_domain1D( domain, begin, end, size, max_size, is_global )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size
      logical, intent(out), optional :: is_global

      if( PRESENT(begin)     )begin     = domain%data%begin
      if( PRESENT(end)       )end       = domain%data%end
      if( PRESENT(size)      )size      = domain%data%size
      if( PRESENT(max_size)  )max_size  = domain%data%max_size
      if( PRESENT(is_global) )is_global = domain%data%is_global
      return
    end subroutine mpp_get_data_domain1D

    subroutine mpp_get_global_domain1D( domain, begin, end, size, max_size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size

      if( PRESENT(begin)    )begin    = domain%global%begin
      if( PRESENT(end)      )end      = domain%global%end
      if( PRESENT(size)     )size     = domain%global%size
      if( PRESENT(max_size) )max_size = domain%global%max_size
      return
    end subroutine mpp_get_global_domain1D

    subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
         x_is_global, y_is_global )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      logical, intent(out), optional :: x_is_global, y_is_global
      call mpp_get_compute_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
      call mpp_get_compute_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
      return
    end subroutine mpp_get_compute_domain2D

    subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
         x_is_global, y_is_global )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      logical, intent(out), optional :: x_is_global, y_is_global
      call mpp_get_data_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
      call mpp_get_data_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
      return
    end subroutine mpp_get_data_domain2D

    subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_global_domain( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_global_domain( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_global_domain2D

    subroutine mpp_get_domain_components( domain, x, y )
      type(domain2D), intent(in) :: domain
      type(domain1D), intent(out), optional :: x, y
      if( PRESENT(x) )x = domain%x
      if( PRESENT(y) )y = domain%y
      return
    end subroutine mpp_get_domain_components

    subroutine mpp_get_compute_domains1D( domain, begin, end, size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional, dimension(:) :: begin, end, size 

      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )
!we use shape instead of size for error checks because size is used as an argument
      if( PRESENT(begin) )then
          if( any(shape(begin).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: begin array size does not match domain.' )
          begin(:) = domain%list(:)%compute%begin
      end if
      if( PRESENT(end) )then
          if( any(shape(end).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: end array size does not match domain.' )
          end(:) = domain%list(:)%compute%end
      end if
      if( PRESENT(size) )then
          if( any(shape(size).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: size array size does not match domain.' )
          size(:) = domain%list(:)%compute%size
      end if
      return
    end subroutine mpp_get_compute_domains1D

    subroutine mpp_get_compute_domains2D( domain, xbegin, xend, xsize, ybegin, yend, ysize )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional, dimension(:) :: xbegin, xend, xsize, ybegin, yend, ysize

      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )

      if( PRESENT(xbegin) )then
          if( size(xbegin).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xbegin array size does not match domain.' )
          xbegin(:) = domain%list(:)%x%compute%begin
      end if
      if( PRESENT(xend) )then
          if( size(xend).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xend array size does not match domain.' )
          xend(:) = domain%list(:)%x%compute%end
      end if
      if( PRESENT(xsize) )then
          if( size(xsize).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xsize array size does not match domain.' )
          xsize(:) = domain%list(:)%x%compute%size
      end if
      if( PRESENT(ybegin) )then
          if( size(ybegin).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ybegin array size does not match domain.' )
          ybegin(:) = domain%list(:)%y%compute%begin
      end if
      if( PRESENT(yend) )then
          if( size(yend).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: yend array size does not match domain.' )
          yend(:) = domain%list(:)%y%compute%end
      end if
      if( PRESENT(ysize) )then
          if( size(ysize).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ysize array size does not match domain.' )
          ysize(:) = domain%list(:)%y%compute%size
      end if
      return
    end subroutine mpp_get_compute_domains2D

    subroutine mpp_get_pelist1D( domain, pelist, pos )
      type(domain1D), intent(in) :: domain
      integer, intent(out) :: pelist(:)
      integer, intent(out), optional :: pos
      integer :: ndivs

      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
      ndivs = size(domain%list)
      
      if( size(pelist).NE.ndivs ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

      pelist(:) = domain%list(0:ndivs-1)%pe
      if( PRESENT(pos) )pos = domain%pos
      return
    end subroutine mpp_get_pelist1D

    subroutine mpp_get_pelist2D( domain, pelist, pos )
      type(domain2D), intent(in) :: domain
      integer, intent(out) :: pelist(:)
      integer, intent(out), optional :: pos

      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
      if( size(pelist).NE.size(domain%list) ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

      pelist(:) = domain%list(:)%pe
      if( PRESENT(pos) )pos = domain%pos
      return
    end subroutine mpp_get_pelist2D

    subroutine mpp_get_layout1D( domain, layout )
      type(domain1D), intent(in) :: domain
      integer, intent(out) :: layout
      
      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

      layout = size(domain%list)
      return
    end subroutine mpp_get_layout1D

    subroutine mpp_get_layout2D( domain, layout )
      type(domain2D), intent(in) :: domain
      integer, intent(out) :: layout(2)
      
      if( .NOT.mpp_domains_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

      layout(1) = size(domain%x%list)
      layout(2) = size(domain%y%list)
      return
    end subroutine mpp_get_layout2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_UPDATE_DOMAINS: fill halos for 2D decomposition            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTOR_FIELD_
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r8_5D
#ifdef  VECTOR_FIELD_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r8_2Dv
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r8_3Dv
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r8_4Dv
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r8_5Dv
#endif
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r8_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c8_5D
#include <mpp_update_domains2D.h>

#ifndef no_8byte_integers
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i8_5D
#include <mpp_update_domains2D.h>

#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_l8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_l8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_l8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_l8_5D
#include <mpp_update_domains2D.h>
#endif

#define VECTOR_FIELD_
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D
#ifdef  VECTOR_FIELD_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r4_2Dv
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r4_3Dv
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r4_4Dv
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r4_5Dv
#endif
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r4_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c4_5D
#include <mpp_update_domains2D.h>

#define MPP_TYPE_ integer(INT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i4_5D
#include <mpp_update_domains2D.h>

#define MPP_TYPE_ logical(INT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_l4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_l4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_l4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_l4_5D
#include <mpp_update_domains2D.h>

!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r8_5D
!#define MPP_TYPE_ real(DOUBLE_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c8_5D
!#define MPP_TYPE_ complex(DOUBLE_KIND)
!#include <mpp_update_domains1D.h>
!
!#ifndef no_8byte_integers
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i8_5D
!#define MPP_TYPE_ integer(LONG_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l8_5D
!#define MPP_TYPE_ logical(LONG_KIND)
!#include <mpp_update_domains1D.h>
!#endif
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r4_5D
!#define MPP_TYPE_ real(FLOAT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c4_5D
!#define MPP_TYPE_ complex(FLOAT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i4_5D
!#define MPP_TYPE_ integer(INT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l4_5D
!#define MPP_TYPE_ logical(INT_KIND)
!#include <mpp_update_domains1D.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_REDUCE: get global max/min of field                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#ifndef no_8byte_integers
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                   MPP_GLOBAL_SUM: global sum of field                       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#ifndef no_8byte_integers
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>
#endif

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_FIELD: get global field from domain field           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_r8_2d
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_c8_2d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_field.h>

#ifndef no_8byte_integers
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_i8_2d
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_l8_2d
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_global_field.h>
#endif

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_r4_2d
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_c4_2d
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_i4_2d
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_l4_2d
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_global_field.h>

end module mpp_domains_mod

#ifdef test_mpp_domains
program mpp_domains_test
  use mpp_mod
  use mpp_domains_mod
  implicit none
  integer :: pe, npes
  type(domain2D) :: domain
  integer :: nx=128, ny=128, nz=40, halo=2, stackmax=32768
  real(DOUBLE_KIND), dimension(:,:,:), allocatable :: local, localy, global, gglobal
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: tk, tk0, tks_per_sec
  integer :: i,j,k, unit=7, layout(2)
  integer :: id
  real :: t
  real :: lsum, gsum
  logical :: debug=.FALSE., opened
  namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, debug

  call mpp_init()
  call mpp_set_warn_level(FATAL)
!possibly open a file called mpp_domains.nml
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='mpp_domains.nml', err=10 )
  read( unit,mpp_domains_nml )
  close(unit)
10 continue

  pe = mpp_pe()
  npes = mpp_npes()

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)

  if( pe.EQ.mpp_root_pe() )then
      print '(a,5i4)', 'npes, nx, ny, nz, halo=', npes, nx, ny, nz, halo
      print *, 'Using NEW domaintypes and calls...'
  end if

  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name='Simple halo update' )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

  allocate( global(1-halo:nx+halo,1-halo:ny+halo,nz) )
  allocate( gglobal(nx,ny,nz) )
  allocate( local(isd:ied,jsd:jed,nz) )
  allocate( localy(isd:ied,jsd:jed,nz) )

!fill in global array: with k.iiijjj
  global = 0.
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           global(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)

!fill in gglobal array
  id = mpp_clock_id( 'Global field', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_global_field( domain, local, gglobal )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( global(1:nx,1:ny,:), gglobal, 'mpp_global_field' )

!fill in halos
  id = mpp_clock_id( 'Halo update', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_update_domains( local, domain )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Halo update' )

!test cyclic boundary conditions
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
       xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name='Cyclic halo update' )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
!fill in cyclic global array
  global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
  global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
  global(1-halo:nx+halo,    1-halo:0,:) = global(1-halo:nx+halo,ny-halo+1:ny,:)
  global(1-halo:nx+halo,ny+1:ny+halo,:) = global(1-halo:nx+halo,1:halo,      :)
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
  id = mpp_clock_id( 'Cyclic halo update', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_update_domains( local, domain )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Cyclic halo update' )

!test folded boundary conditions
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
       xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name='Folded halo update' )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
!fill in folded north edge, cyclic east and west edge
  global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
  global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
  global(1-halo:nx+halo,ny+1:ny+halo,:) = global(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
  global(1-halo:nx+halo,1-halo:0,:) = 0
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
  id = mpp_clock_id( 'Folded halo update', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_update_domains( local, domain )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( local(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), 'Folded halo update' )
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos: partial update
  id = mpp_clock_id( 'Folded halo N+E update', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_update_domains( local, domain, NUPDATE+EUPDATE )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( local(is:ied,js:jed,:), global(is:ied,js:jed,:), 'Folded halo N+E update' )
!folded, with sign flip at fold (vector component)
!fill in folded north edge, cyclic east and west edge
  global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
  global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
  global(1-halo:nx+halo-1,ny+1:ny+halo,:) = -global(nx+halo-1:1-halo:-1,ny-1:ny-halo:-1,:)
  global(nx+halo,ny+1:ny+halo,:) = -global(nx-halo,ny-1:ny-halo:-1,:)
  global(1-halo:nx+halo,1-halo:0,:) = 0
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
  localy = local
!fill in halos
  id = mpp_clock_id( 'Folded halo update, vector field', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
  call mpp_clock_begin(id)
  call mpp_update_domains( local, localy, domain, gridtype=BGRID_NE )
  call mpp_clock_end  (id)
!compare checksums between global and local arrays
  call compare_checksums( local(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), 'Folded halo update, vector field' )

!timing tests
  if( pe.EQ.mpp_root_pe() )print '(a)', 'TIMING TESTS:'
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
       xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name='Timing' )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
!fill in cyclic global array
  global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
  global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
  global(1-halo:nx+halo,    1-halo:0,:) = global(1-halo:nx+halo,ny-halo+1:ny,:)
  global(1-halo:nx+halo,ny+1:ny+halo,:) = global(1-halo:nx+halo,1:halo,      :)
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
  call mpp_sync()          !this ensures you time only the update_domains call
  call SYSTEM_CLOCK(tk0)
  call mpp_update_domains( local, domain )
  call SYSTEM_CLOCK(tk)
  t = float(tk-tk0)/tks_per_sec
!compare checksums between global and local arrays
  call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Timed cyclic halo update' )
!words transferred
  j = ( (ied-isd+1)*(jed-jsd+1) - (ie-is+1)*(je-js+1) )*nz
  call mpp_max(j)
  if( pe.EQ.mpp_root_pe() ) &
       print '(a,i4,i8,es12.4,f10.3)', 'Halo width, words, time (sec), bandwidth (MB/sec)=', halo, j, t, j*8e-6/t
!test and time mpp_global_sum
  gsum = sum( global(1:nx,1:ny,:) )
  call mpp_sync()          !this ensures you time only the update_domains call
  call SYSTEM_CLOCK(tk0)
  lsum = mpp_global_sum( domain, local )
  call SYSTEM_CLOCK(tk)
  t = float(tk-tk0)/tks_per_sec
  if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Fast sum=', lsum, gsum, ' Time (sec)=', t
  call mpp_sync()          !this ensures you time only the update_domains call
  call SYSTEM_CLOCK(tk0)
  lsum = mpp_global_sum( domain, local, BITWISE_EXACT_SUM )
  call SYSTEM_CLOCK(tk)
  t = float(tk-tk0)/tks_per_sec
  if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Bitwise-exact sum=', lsum, gsum, ' Time (sec)=', t

  call mpp_domains_exit()
  call mpp_exit()

contains
  subroutine compare_checksums( a, b, string )
    real(DOUBLE_KIND), intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer :: i, j

    call mpp_sync()
    i = mpp_chksum( a, (/pe/) )
    j = mpp_chksum( b, (/pe/) )
    if( i.EQ.j )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': chksums are correct.' )
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums
end program mpp_domains_test
#endif
