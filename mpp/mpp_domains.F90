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
       '$Id: mpp_domains.F90,v 6.2 2001/10/25 17:55:11 fms Exp $'
  character(len=128), private :: name= &
       '$Name: fez $'

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5,   stdout=6,   stderr=0
#endif
  integer, private :: pe, npes, root

  type, public :: domain_axis_spec        !type used to specify index limits along an axis of a domain
#ifndef MPP_DOMAINS_PRE_FEZ
     private
#endif
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
#ifdef MPP_DOMAINS_PRE_FEZ
!next two items for backward compatibility
     integer :: start_index, end_index
#endif
  end type domain_axis_spec
  type, public :: domain1D
#ifndef MPP_DOMAINS_PRE_FEZ
     private
#endif
     type(domain_axis_spec) :: compute, data, global, active, putb, getb, putf, getf
     logical :: mustputb, mustgetb, mustputf, mustgetf, folded
     type(domain1D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer, pointer :: pemap(:)
#ifdef MPP_DOMAINS_PRE_FEZ
!next two items for backward compatibility
     integer :: ndomains        !number of domains over which data has been distributed
     integer, pointer :: pelist(:), sizelist(:)
#endif
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type, public :: domain2D
#ifndef MPP_DOMAINS_PRE_FEZ
     private
#endif
     type(domain1D) :: x
     type(domain1D) :: y
     type(domain2D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
!     integer, pointer :: pemap(:,:)
  end type domain2D
  type(domain1D), public :: NULL_DOMAIN1D
  type(domain2D), public :: NULL_DOMAIN2D

!parameters used to define domains: these are passed to the flags argument of mpp_define_domains
!  if compute domain is to have global extent, set GLOBAL_COMPUTE_DOMAIN
!  if data domain is to have global extent, set GLOBAL_DATA_DOMAIN
!  if global domain has periodic boundaries, set CYCLIC_GLOBAL_DOMAIN
!  sum flags together if more than one of the above conditions is to be met.
  integer, parameter, public :: GLOBAL_COMPUTE_DOMAIN=1, GLOBAL_DATA_DOMAIN=2, CYCLIC_GLOBAL_DOMAIN=4
!the next 3 pairs must use the same values
  integer, parameter, public :: FOLD_WEST_EDGE =8, FOLD_EAST_EDGE =16
  integer, parameter, public :: FOLD_SOUTH_EDGE=8, FOLD_NORTH_EDGE=16
  integer, parameter, public :: FOLD_LOWER_EDGE=8, FOLD_UPPER_EDGE=16
  integer, parameter, public :: WUPDATE=1, EUPDATE=2, SUPDATE=4, NUPDATE=8, TRANSPOSE=16
  integer, parameter, public :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE
  integer, parameter, public :: VECTOR_COMPONENT=1
  integer, parameter, public :: BITWISE_EXACT_SUM=1

  integer, private :: tk
  logical, private :: verbose=.FALSE., debug=.FALSE.
  logical, private :: mpp_domains_initialized=.FALSE.

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
  integer, allocatable, dimension(:) :: local_dom, remote_dom

!public interfaces
  interface mpp_define_domains
#ifdef MPP_DOMAINS_PRE_FEZ
     module procedure mpp_define_domains1D_old
     module procedure mpp_define_domains2D_old
#endif
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
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

     module procedure mpp_update_domain1D_r8_2d
     module procedure mpp_update_domain1D_r8_3d
     module procedure mpp_update_domain1D_r8_4d
     module procedure mpp_update_domain1D_r8_5d
     module procedure mpp_update_domain1D_c8_2d
     module procedure mpp_update_domain1D_c8_3d
     module procedure mpp_update_domain1D_c8_4d
     module procedure mpp_update_domain1D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_update_domain1D_i8_2d
     module procedure mpp_update_domain1D_i8_3d
     module procedure mpp_update_domain1D_i8_4d
     module procedure mpp_update_domain1D_i8_5d
     module procedure mpp_update_domain1D_l8_2d
     module procedure mpp_update_domain1D_l8_3d
     module procedure mpp_update_domain1D_l8_4d
     module procedure mpp_update_domain1D_l8_5d
#endif
     module procedure mpp_update_domain1D_r4_2d
     module procedure mpp_update_domain1D_r4_3d
     module procedure mpp_update_domain1D_r4_4d
     module procedure mpp_update_domain1D_r4_5d
     module procedure mpp_update_domain1D_c4_2d
     module procedure mpp_update_domain1D_c4_3d
     module procedure mpp_update_domain1D_c4_4d
     module procedure mpp_update_domain1D_c4_5d
     module procedure mpp_update_domain1D_i4_2d
     module procedure mpp_update_domain1D_i4_3d
     module procedure mpp_update_domain1D_i4_4d
     module procedure mpp_update_domain1D_i4_5d
     module procedure mpp_update_domain1D_l4_2d
     module procedure mpp_update_domain1D_l4_3d
     module procedure mpp_update_domain1D_l4_4d
     module procedure mpp_update_domain1D_l4_5d
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

  interface mpp_get_active_domain
     module procedure mpp_get_active_domain1D
     module procedure mpp_get_active_domain2D
  end interface

  interface mpp_set_active_domain
     module procedure mpp_set_active_domain1D
     module procedure mpp_set_active_domain2D
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
            mpp_get_active_domain, mpp_get_compute_domain, mpp_get_compute_domains, mpp_get_data_domain, mpp_get_global_domain, &
            mpp_get_domain_components, mpp_get_layout, mpp_get_pelist, mpp_set_active_domain, mpp_update_domains, &
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

      if( mpp_domains_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      pe = mpp_pe()
      npes = mpp_npes()
      root = mpp_root_pe()
      mpp_domains_initialized = .TRUE.
      if( pe.EQ.root )write( stdout,'(/a)' )'MPP_DOMAINS module '//trim(version)

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      NULL_DOMAIN1D%global%begin  = -1; NULL_DOMAIN1D%global%end  = -1; NULL_DOMAIN1D%global%size = 0
      NULL_DOMAIN1D%data%begin    = -1; NULL_DOMAIN1D%data%end    = -1; NULL_DOMAIN1D%data%size = 0
      NULL_DOMAIN1D%compute%begin = -1; NULL_DOMAIN1D%compute%end = -1; NULL_DOMAIN1D%compute%size = 0
      NULL_DOMAIN1D%pe = NULL_PE
      NULL_DOMAIN2D%x = NULL_DOMAIN1D
      NULL_DOMAIN2D%y = NULL_DOMAIN1D
      NULL_DOMAIN2D%pe = NULL_PE
#ifdef MPP_DOMAINS_PRE_FEZ
!backward-compatibility
      NULL_DOMAIN1D%global%start_index  = -1; NULL_DOMAIN1D%global%end_index  = -1
      NULL_DOMAIN1D%data%start_index    = -1; NULL_DOMAIN1D%data%end_index    = -1
      NULL_DOMAIN1D%compute%start_index = -1; NULL_DOMAIN1D%compute%end_index = -1
#endif
      
      return
    end subroutine mpp_domains_init

    subroutine mpp_domains_set_stack_size(n)
!set the mpp_domains_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=8) :: text
#ifdef use_shmalloc
      call mpp_malloc( ptr_stack, n, mpp_domains_stack_size )
#else
      if( n.GT.mpp_domains_stack_size .AND. allocated(mpp_domains_stack) )deallocate(mpp_domains_stack)
      if( .NOT.allocated(mpp_domains_stack) )then
          allocate( mpp_domains_stack(n) )
          mpp_domains_stack_size = n
      end if
#endif
      write( text,'(i8)' )n
      if( pe.EQ.root )call mpp_error( NOTE, 'MPP_DOMAINS_SET_STACK_SIZE: stack size set to '//text//'.' )

      return
    end subroutine mpp_domains_set_stack_size

    subroutine mpp_domains_exit()
!currently does not have much to do, but provides the possibility of re-initialization
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DOMAINS_EXIT: You must first call mpp_domains_init.' )
      call mpp_max(mpp_domains_stack_hwm)
      if( pe.EQ.root )then
          write( stdout,'(/a)' )'Exiting MPP_DOMAINS module...'
          write( stdout,* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      end if
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

      logical :: compute_domain_is_global, data_domain_is_global, global_domain_is_cyclic
      logical :: upper_edge_is_folded, lower_edge_is_folded
      integer :: ndiv, n, isg, ieg, is, ie, i, off, pos, hs, he
      integer, allocatable :: pes(:)
      logical, allocatable :: mask(:)
      integer :: halosz
      
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
!get global indices
      isg = global_indices(1)
      ieg = global_indices(2)
      if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
!get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
      if( PRESENT(pelist) )then
          if( .NOT.any(pelist.EQ.pe) )then
              write( stderr,* )'pe=', pe, ' pelist=', pelist
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
      global_domain_is_cyclic  = .FALSE.
      upper_edge_is_folded = .FALSE.
      lower_edge_is_folded = .FALSE.
      if( PRESENT(flags) )then
!NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
          compute_domain_is_global = ndivs.EQ.1
!if compute domain is global, data domain must also be
          data_domain_is_global    = BTEST(flags,1) .OR. compute_domain_is_global
          global_domain_is_cyclic  = BTEST(flags,2)
          if( global_domain_is_cyclic .AND. ( BTEST(flags,3) .OR. BTEST(flags,4) ) ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: An axis cannot be both cyclic and folded.' )
          lower_edge_is_folded = BTEST(flags,3) .AND. halosz.NE.0
          upper_edge_is_folded = BTEST(flags,4) .AND. halosz.NE.0
      end if

!set up links list
      if( ASSOCIATED(domain%list) )NULLIFY(domain%list)
      if( lower_edge_is_folded .OR. upper_edge_is_folded )then
          domain%folded = .TRUE.
!folded domains requires twice as many links
!masks and PEs for domains over the fold will be assigned by the calling routine (usually mpp_define_domains2D)
          allocate( domain%list(0:2*ndivs-1) )
          domain%list(ndivs:2*ndivs-1)%folded = .TRUE.
      else
          domain%folded = .FALSE.
          allocate( domain%list(0:  ndivs-1) )
      end if
      domain%list(0:ndivs-1)%folded = .FALSE.

!set global domain
      domain%list(:)%global%begin     = isg
      domain%list(:)%global%end       = ieg
      domain%list(:)%global%size      = ieg-isg+1
      domain%list(:)%global%max_size  = ieg-isg+1
      domain%list(:)%global%is_global = .TRUE. !always

!get compute domain
!      nullify(domain%pemap)
!      allocate( domain%pemap(isg:ieg) )
!      domain%pemap(isg:ieg) = NULL_PE
      if( compute_domain_is_global )then
          domain%list(:)%compute%begin = isg
          domain%list(:)%compute%end   = ieg
          domain%list(:)%compute%is_global = .TRUE.
          domain%list(:)%pe = pes(:)
          domain%pos = 0
!          domain%pemap(isg:ieg) = pes(0)
      else
          domain%list(:)%compute%is_global = .FALSE.
          is = isg
          n = 0
          do ndiv=0,ndivs-1
             if( PRESENT(extent) )then
                 if( extent(ndiv).LE.0 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extents must be positive definite.' )
                 ie = is + extent(ndiv) - 1
                 if( ndiv.EQ.ndivs-1 .AND. ie.NE.ieg ) &
                      call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extent array limits do not match global domain.' )
             else
                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1
             end if
             domain%list(ndiv)%compute%begin = is
             domain%list(ndiv)%compute%end   = ie
             if( lower_edge_is_folded .OR. upper_edge_is_folded )then
                 domain%list(2*ndivs-1-ndiv)%compute%begin = is
                 domain%list(2*ndivs-1-ndiv)%compute%end   = ie
             end if
             if( mask(ndiv) )then
                 domain%list(ndiv)%pe = pes(n)
!                 domain%pemap(is:ie) = pes(n)
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
!apply global flags
      if( data_domain_is_global )then
          domain%list(:)%data%begin  = isg
          domain%list(:)%data%end    = ieg
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
      domain%active = domain%data !active domain is initialized to data domain

!create list of put and get domains
      domain%list(:)%mustputb = .FALSE.
      domain%list(:)%mustgetb = .FALSE.
      domain%list(:)%mustputf = .FALSE.
      domain%list(:)%mustgetf = .FALSE.
      if( halosz.NE.0 )then
!overlaps with your own domain (only for cyclic domains)
          pos = domain%pos
          if( global_domain_is_cyclic )then
              off = ieg-isg+1 !use this offset for repositioning
!putb: repositioned domain upper halo must lie in my compute domain
              domain%list(pos)%mustputb = if_overlap( domain%compute%end+1-off, domain%data%end-off, &
                                                      domain%compute%begin, domain%compute%end, &
                                                      domain%list(pos)%putb%begin, domain%list(pos)%putb%end )
!getb: my lower halo must lie in repositioned domain compute domain
              domain%list(pos)%mustgetb = if_overlap( domain%data%begin, domain%compute%begin-1, &
                                                      domain%compute%begin-off, domain%compute%end-off, &
                                                      domain%list(pos)%getb%begin, domain%list(pos)%getb%end )
!putf: repositioned domain lower halo must lie in my compute domain
              domain%list(pos)%mustputf = if_overlap( domain%data%begin+off, domain%compute%begin-1+off, &
                                                      domain%compute%begin, domain%compute%end, &
                                                      domain%list(pos)%putf%begin, domain%list(pos)%putf%end )
!getf: my upper halo must lie in domain compute domain
              domain%list(pos)%mustgetf = if_overlap( domain%compute%end+1, domain%data%end, &
                                                      domain%compute%begin+off, domain%compute%end+off, &
                                                      domain%list(pos)%getf%begin, domain%list(pos)%getf%end )
          end if

!overlaps with other domains (if any)
          do ndiv = 1,size(domain%list)-1
             pos = domain%pos-ndiv
             if( pos.GE.0 )then
!putb: domain%list(pos) upper halo must lie in my compute domain
                 domain%list(pos)%mustputb = if_overlap( domain%list(pos)%compute%end+1, domain%list(pos)%data%end, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putb%begin, domain%list(pos)%putb%end ) .AND. mask(pos)
!getb: my lower halo must lie in domain%list(pos) compute domain
                 domain%list(pos)%mustgetb = if_overlap( domain%data%begin, domain%compute%begin-1, &
                                                         domain%list(pos)%compute%begin, domain%list(pos)%compute%end, &
                                                         domain%list(pos)%getb%begin, domain%list(pos)%getb%end ) .AND. mask(pos)
             else if( global_domain_is_cyclic )then
                 pos = pos + ndivs
                 off = ieg-isg+1 !use this offset for repositioning
!putb: repositioned domain%list(pos) upper halo must lie in my compute domain
                 domain%list(pos)%mustputb = if_overlap( domain%list(pos)%compute%end+1-off, domain%list(pos)%data%end-off, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putb%begin, domain%list(pos)%putb%end ) .AND. mask(pos)
!getb: my lower halo must lie in repositioned domain%list(pos) compute domain
                 domain%list(pos)%mustgetb = if_overlap( domain%data%begin, domain%compute%begin-1, &
                                                         domain%list(pos)%compute%begin-off, domain%list(pos)%compute%end-off, &
                                                         domain%list(pos)%getb%begin, domain%list(pos)%getb%end ) .AND. mask(pos)
             else if( lower_edge_is_folded )then
!use putf rather than putb!
                 pos = pos + 2*ndivs
!putf: repositioned and reversed domain%list(pos) LOWER halo must lie in my compute domain
                 off = 2*isg - 1
                 domain%list(pos)%mustputf = if_overlap( off-domain%list(pos)%compute%begin-1, off-domain%list(pos)%data%begin, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putf%begin, domain%list(pos)%putf%end )
!getb: my lower halo must lie in repositioned and reversed domain%list(pos) compute domain
                 domain%list(pos)%mustgetb = if_overlap( domain%data%begin, domain%compute%begin-1, &
                                                         off-domain%list(pos)%compute%end, off-domain%list(pos)%compute%begin, &
                                                         domain%list(pos)%getb%begin, domain%list(pos)%getb%end )
             end if

             pos = domain%pos+ndiv
             if( pos.LT.ndivs )then
!putf: domain%list(pos) lower halo must lie in my compute domain
                 domain%list(pos)%mustputf = if_overlap( domain%list(pos)%data%begin, domain%list(pos)%compute%begin-1, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putf%begin, domain%list(pos)%putf%end ) .AND. mask(pos)
!getf: my upper halo must lie in domain%list(pos) compute domain
                 domain%list(pos)%mustgetf = if_overlap( domain%compute%end+1, domain%data%end, &
                                                         domain%list(pos)%compute%begin, domain%list(pos)%compute%end, &
                                                         domain%list(pos)%getf%begin, domain%list(pos)%getf%end ) .AND. mask(pos)
             else if( global_domain_is_cyclic )then
                 pos = pos - ndivs
                 off = ieg-isg+1
!putf: repositioned domain%list(pos) lower halo must lie in my compute domain
                 domain%list(pos)%mustputf = if_overlap( domain%list(pos)%data%begin+off, domain%list(pos)%compute%begin-1+off, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putf%begin, domain%list(pos)%putf%end ) .AND. mask(pos)
!getf: my upper halo must lie in domain%list(pos) compute domain
                 domain%list(pos)%mustgetf = if_overlap( domain%compute%end+1, domain%data%end, &
                                                         domain%list(pos)%compute%begin+off, domain%list(pos)%compute%end+off, &
                                                         domain%list(pos)%getf%begin, domain%list(pos)%getf%end ) .AND. mask(pos)
             else if( upper_edge_is_folded )then
!use putb rather than putf!
                 off = 2*ieg + 1
!putb: repositioned and reversed domain%list(pos) UPPER halo must lie in my compute domain
                 domain%list(pos)%mustputb = if_overlap( off-domain%list(pos)%data%end, off-domain%list(pos)%compute%end+1, &
                                                         domain%compute%begin, domain%compute%end, &
                                                         domain%list(pos)%putb%begin, domain%list(pos)%putb%end )
!getf: my upper halo must lie in repositioned and reversed domain%list(pos) compute domain
                 domain%list(pos)%mustgetf = if_overlap( domain%compute%end+1, domain%data%end, &
                                                         off-domain%list(pos)%compute%end, off-domain%list(pos)%compute%begin, &
                                                         domain%list(pos)%getf%begin, domain%list(pos)%getf%end )
             end if
          end do
          domain%list(:)%putb%size = domain%list(:)%putb%end - domain%list(:)%putb%begin + 1
          domain%list(:)%getb%size = domain%list(:)%getb%end - domain%list(:)%getb%begin + 1
          domain%list(:)%putf%size = domain%list(:)%putf%end - domain%list(:)%putf%begin + 1
          domain%list(:)%getf%size = domain%list(:)%getf%end - domain%list(:)%getf%begin + 1
      end if
!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask )
#ifdef MPP_DOMAINS_PRE_FEZ
!backward compatibility
      domain%compute%start_index = domain%compute%begin
      domain%compute%end_index = domain%compute%end
      domain%data%start_index = domain%data%begin
      domain%data%end_index = domain%data%end
      domain%global%start_index = domain%global%begin
      domain%global%end_index = domain%global%end
      domain%list(:)%compute%start_index = domain%list(:)%compute%begin
      domain%list(:)%compute%end_index = domain%list(:)%compute%end
      domain%list(:)%data%start_index = domain%list(:)%data%begin
      domain%list(:)%data%end_index = domain%list(:)%data%end
      domain%list(:)%global%start_index = domain%list(:)%global%begin
      domain%list(:)%global%end_index = domain%list(:)%global%end
      domain%ndomains = ndivs
      allocate( domain%sizelist(ndivs) )
      allocate( domain%pelist(ndivs) )
      domain%sizelist(:) = domain%list(0:ndivs-1)%compute%size
      domain%pelist(:) = domain%list(0:ndivs-1)%pe
#endif
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
          if( debug )write( stderr,'(a,7i4)' ) &
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
      integer :: i, j, n, ipos, jpos, pos, remote_pos, ndivx, ndivy, isg, ieg, jsg, jeg
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
      if( debug )write( stderr, * )'pe, ipos, jpos=', pe, ipos, jpos, ' pearray(:,jpos)=', pearray(:,jpos), &
                                                                      ' pearray(ipos,:)=', pearray(ipos,:)

!do domain decomposition using 1D versions in X and Y
      call mpp_define_domains( global_indices(1:2), ndivx, domain%x, &
           pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos) )
      call mpp_define_domains( global_indices(3:4), ndivy, domain%y, &
           pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:) )
      if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
      domain%pos = pos

      if( size(domain%y%list).EQ.2*ndivy )then ! there is a fold in Y
!check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x%list(i)%compute%size.NE.domain%x%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
          domain%y%list(ndivy:2*ndivy-1)%pe = pearray(n-ipos,ndivy-1:0:-1)
      end if
      if( size(domain%x%list).EQ.2*ndivx )then ! there is a fold in X
!check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y%list(i)%compute%size.NE.domain%y%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
          domain%x%list(ndivx:2*ndivx-1)%pe = pearray(ndivx-1:0:-1,n-jpos)
      end if

      if( .NOT.allocated(local_dom) )allocate( local_dom(8) )
      if( .NOT.allocated(remote_dom) )allocate( remote_dom(8) )
      call mpp_get_compute_domain( domain, local_dom(1), local_dom(2), local_dom(3), local_dom(4) )
      call mpp_get_data_domain   ( domain, local_dom(5), local_dom(6), local_dom(7), local_dom(8) )
      if( debug )write( stderr,'(a,9i4)' )'pe, domain=', pe, local_dom
      n = size(pes)
      if( ASSOCIATED(domain%list) )NULLIFY(domain%list)
      allocate( domain%list(0:n-1) ) !this is only used for storage of remote compute and data domain info
      if( pe.EQ.root .AND. PRESENT(name) )then
          write( stdout, '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( stdout, '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
      end if
      call mpp_sync()
!      if( ASSOCIATED(domain%pemap) )NULLIFY(domain%pemap)
!      ALLOCATE( domain%pemap(isg:ieg,jsg:jeg) )
!      domain%pemap(:,:) = NULL_PE
      do i = 0,n-1
         remote_pos = mod(pos+i,n)
         domain%list(remote_pos)%pe = pes(remote_pos)
         call mpp_transmit( local_dom, 8, pes(mod(pos+n-i,n)), remote_dom, 8, pes(remote_pos) )
         domain%list(remote_pos)%x%compute%begin = remote_dom(1)
         domain%list(remote_pos)%x%compute%end   = remote_dom(2)
         domain%list(remote_pos)%y%compute%begin = remote_dom(3)
         domain%list(remote_pos)%y%compute%end   = remote_dom(4)
         domain%list(remote_pos)%x%data%begin    = remote_dom(5)
         domain%list(remote_pos)%x%data%end      = remote_dom(6)
         domain%list(remote_pos)%y%data%begin    = remote_dom(7)
         domain%list(remote_pos)%y%data%end      = remote_dom(8)
!         domain%pemap(remote_dom(1):remote_dom(2),remote_dom(3):remote_dom(4)) = pes(remote_pos)
         if( pe.EQ.root .AND. PRESENT(name) )write( stdout, '(2x,i3,x,4i5,3x,4i5)' )pes(remote_pos), remote_dom
      end do
      call mpp_sync_self(pes)
      domain%list(:)%x%compute%size = domain%list(:)%x%compute%end - domain%list(:)%x%compute%begin + 1
      domain%list(:)%y%compute%size = domain%list(:)%y%compute%end - domain%list(:)%y%compute%begin + 1
      domain%list(:)%x%data%size = domain%list(:)%x%data%end - domain%list(:)%x%data%begin + 1
      domain%list(:)%y%data%size = domain%list(:)%y%data%end - domain%list(:)%y%data%begin + 1
!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask, pearray )
#ifdef MPP_DOMAINS_PRE_FEZ
      domain%list(:)%x%compute%start_index = domain%list(:)%x%compute%begin
      domain%list(:)%x%compute%end_index = domain%list(:)%x%compute%end
      domain%list(:)%y%compute%start_index = domain%list(:)%y%compute%begin
      domain%list(:)%y%compute%end_index = domain%list(:)%y%compute%end
      domain%list(:)%x%global%start_index = domain%list(:)%x%global%begin
      domain%list(:)%x%global%end_index = domain%list(:)%x%global%end
      domain%list(:)%y%global%start_index = domain%list(:)%y%global%begin
      domain%list(:)%y%global%end_index = domain%list(:)%y%global%end
      domain%list(:)%x%data%start_index = domain%list(:)%x%data%begin
      domain%list(:)%x%data%end_index = domain%list(:)%x%data%end
      domain%list(:)%y%data%start_index = domain%list(:)%y%data%begin
      domain%list(:)%y%data%end_index = domain%list(:)%y%data%end
#endif
          
      return
    end subroutine mpp_define_domains2D

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

    subroutine mpp_get_active_domain1D( domain, begin, end, size, max_size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size

      if( PRESENT(begin)    )begin    = domain%active%begin
      if( PRESENT(end)      )end      = domain%active%end
      if( PRESENT(size)     )size     = domain%active%size
      if( PRESENT(max_size) )max_size = domain%active%max_size
      return
    end subroutine mpp_get_active_domain1D

    subroutine mpp_set_active_domain1D( domain, begin, end )
      type(domain1D), intent(inout) :: domain
      integer, intent(in), optional :: begin, end

      if( PRESENT(begin) )then
          if( begin.GT.domain%compute%begin ) &
               call mpp_error( FATAL, 'MPP_SET_ACTIVE_DOMAIN1D: active domain must contain compute domain.' )
          if( begin.LT.domain%data%begin ) &
               call mpp_error( FATAL, 'MPP_SET_ACTIVE_DOMAIN1D: active domain must lie within data domain.' )
          domain%active%begin = begin
      end if
      if( PRESENT(end) )then
          if( end.LT.domain%compute%end ) &
               call mpp_error( FATAL, 'MPP_SET_ACTIVE_DOMAIN1D: active domain must contain compute domain.' )
          if( end.GT.domain%data%end ) &
               call mpp_error( FATAL, 'MPP_SET_ACTIVE_DOMAIN1D: active domain must lie within data domain.' )
          domain%active%end   = end
      end if
      domain%active%size = domain%active%end - domain%active%begin + 1
      return
    end subroutine mpp_set_active_domain1D

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

    subroutine mpp_get_active_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_active_domain( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_active_domain( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_active_domain2D

    subroutine mpp_set_active_domain2D( domain, xbegin, xend, ybegin, yend )
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: xbegin, xend, ybegin, yend
      call mpp_set_active_domain( domain%x, xbegin, xend )
      call mpp_set_active_domain( domain%y, ybegin, yend )
      return
    end subroutine mpp_set_active_domain2D

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
      if( domain%folded )ndivs = ndivs/2
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

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r8_5D
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D
#define MPP_TYPE_ complex(DOUBLE_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#ifndef no_8byte_integers
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D
#define MPP_TYPE_ integer(LONG_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l8_5D
#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>
#endif

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D
#define MPP_TYPE_ real(FLOAT_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D
#define MPP_TYPE_ complex(FLOAT_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D
#define MPP_TYPE_ integer(INT_KIND)
#undef MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l4_5D
#define MPP_TYPE_ logical(INT_KIND)
#define MPP_TYPE_IS_LOGICAL_
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r8_5D
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_update_domains1D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c8_5D
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_update_domains1D.h>

#ifndef no_8byte_integers
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i8_5D
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_update_domains1D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l8_5D
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_update_domains1D.h>
#endif

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r4_5D
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_update_domains1D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c4_5D
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_update_domains1D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i4_5D
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_update_domains1D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l4_5D
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_update_domains1D.h>

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

#ifdef MPP_DOMAINS_PRE_FEZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!      OLD ROUTINES DUE TO BE RETIRED WITH EUGENE RELEASE OF FMS              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_define_domains1D_old( global_indices, domains, pelist, flags, halo, extent )
!routine to divide global array indices among domains, and assign domains to PEs
!domain are an array of type(domain1D) of required size
!ARGUMENTS:
!      global_indices(2)=(isg,ieg) gives the extent of global domain
!      domain are an array of type(domain1D) of required size
!      (optional) pelist list of PEs to which domains are to be assigned (default 0...npes-1)
!      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
!      (optional) halo defines halo width (currently the same on both sides)
!      (optional) array extent defines width of each domain (used for non-uniform domain decomp, for e.g load-balancing)
!  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
!  as the input parameters permit
      integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
      type(domain1D), intent(out), target :: domains(:)
      integer, intent(in), optional :: pelist(:)
      integer, intent(in), optional :: flags, halo
      integer, intent(in), optional :: extent(:)

      type(domain1D) :: domain

      call mpp_define_domains( global_indices, size(domains), domain, pelist, flags, halo, extent )
      domains(:) = domain%list(:)

      return
    end subroutine mpp_define_domains1D_old

    subroutine mpp_define_domains2D_old( global_indices, domains, pelist, xflags, yflags, xhalo, yhalo, xextent, yextent, &
                                     domain_layout, pe_layout )
!define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      type(domain2D), intent(out), target :: domains(0:)
      integer, intent(in), optional :: pelist(:)
      integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
      integer, intent(in), optional :: xextent(:), yextent(:)
      integer, intent(in), optional :: domain_layout(2), pe_layout(2)

      type(domain2D) :: domain
      integer :: layout(2)
      if( PRESENT(pe_layout) )then
          layout = pe_layout
      else
          call mpp_define_layout( global_indices, size(domains), layout )
      end if

      call mpp_define_domains( global_indices, layout, domain, pelist, &
                                               xflags, yflags, xhalo, yhalo, xextent, yextent )
      domains(:) = domain%list(:)
!call again to fill in all the values at domain%pos
      call mpp_define_domains( global_indices, layout, domains(domain%pos), pelist, &
                                               xflags, yflags, xhalo, yhalo, xextent, yextent )

      return
    end subroutine mpp_define_domains2D_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GET_GLOBAL: get global field from domain field             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define mpp_get_global mpp_global_field
#define MPP_GET_GLOBAL mpp_global_field
#endif MPP_DOMAINS_PRE_FEZ
end module mpp_domains_mod

#ifdef test_mpp_domains
program mpp_domains_test
  use mpp_mod
  use mpp_domains_mod
  implicit none
  integer :: pe, npes, root
  type(domain2D) :: domain
  integer :: nx=128, ny=128, nz=40, halo=2, stackmax=32768
  real(DOUBLE_KIND), dimension(:,:,:), allocatable :: local, global, gglobal
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: tk, tk0, tks_per_sec
  integer :: i,j,k, unit=7, layout(2)
  real :: t
  real :: lsum, gsum
  logical :: debug=.FALSE., opened
  namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, debug

  call mpp_init()
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
  root = mpp_root_pe()

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init()
  end if
  call mpp_domains_set_stack_size(stackmax)

  if( pe.EQ.root )then
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
!fill in halos
  call mpp_update_domains( local, domain )
!compare checksums between global and local arrays
  call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Halo update' )

!fill in gglobal array
  call mpp_global_field( domain, local, gglobal )
!compare checksums between global and local arrays
  call compare_checksums( global(1:nx,1:ny,:), gglobal, 'mpp_global_field' )

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
  call mpp_update_domains( local, domain )
!compare checksums between global and local arrays
  call compare_checksums( local(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), 'Folded halo update' )
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos: partial update
  call mpp_update_domains( local, domain, NUPDATE+EUPDATE )
!compare checksums between global and local arrays
  call compare_checksums( local(is:ied,js:jed,:), global(is:ied,js:jed,:), 'Folded halo N+E update' )
!folded, with sign flip at fold (vector component)
!fill in folded north edge, cyclic east and west edge
  global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
  global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
  global(1-halo:nx+halo,ny+1:ny+halo,:) = -global(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
  global(1-halo:nx+halo,1-halo:0,:) = 0
!fill in local array
  local = 0.
  local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
  call mpp_update_domains( local, domain, type=VECTOR_COMPONENT )
!compare checksums between global and local arrays
  call compare_checksums( local(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), 'Folded halo update, vector field' )

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
  call mpp_update_domains( local, domain )
!compare checksums between global and local arrays
  call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Cyclic halo update' )

!timing tests
  if( pe.EQ.root )print '(a)', 'TIMING TESTS:'
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
       xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name='Timing' )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
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
  if( pe.EQ.root ) &
       print '(a,i4,i8,es12.4,f10.3)', 'Halo width, words, time (sec), bandwidth (MB/sec)=', halo, j, t, j*8e-6/t
!test and time mpp_global_sum
  gsum = sum( global(1:nx,1:ny,:) )
  call mpp_sync()          !this ensures you time only the update_domains call
  call SYSTEM_CLOCK(tk0)
  lsum = mpp_global_sum( domain, local )
  call SYSTEM_CLOCK(tk)
  t = float(tk-tk0)/tks_per_sec
  if( pe.EQ.root )print '(a,2es15.8,a,es12.4)', 'Fast sum=', lsum, gsum, ' Time (sec)=', t
  call mpp_sync()          !this ensures you time only the update_domains call
  call SYSTEM_CLOCK(tk0)
  lsum = mpp_global_sum( domain, local, BITWISE_EXACT_SUM )
  call SYSTEM_CLOCK(tk)
  t = float(tk-tk0)/tks_per_sec
  if( pe.EQ.root )print '(a,2es15.8,a,es12.4)', 'Bitwise-exact sum=', lsum, gsum, ' Time (sec)=', t

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
        if( pe.EQ.root )call mpp_error( NOTE, trim(string)//': chksums are correct.' )
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums
end program mpp_domains_test
#endif
