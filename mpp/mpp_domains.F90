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
       '$Id: mpp_domains.F90,v 6.0 2001/03/06 20:26:45 fms Exp $'
  character(len=128), private :: name= &
       '$Name: eugene $'

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5,   stdout=6,   stderr=0
#endif
  integer, private :: pe, npes, root

  type, public :: domain_axis_spec        !type used to specify index limits along an axis of a domain
!     private
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     integer :: start_index, end_index !OBSOLETE , size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec
  type, public :: domain1D
!     private
     type(domain_axis_spec) :: compute, data, global, active, putb, getb, putf, getf
     integer :: ndomains!OBSOLETE: number of domains over which data has been distributed
     integer :: pe              !PE to which this domain is assigned
     integer, pointer :: pelist(:), sizelist(:) !OBSOLETE
     integer :: lhalo, rhalo!OBSOLETE can be used to keep track of active (i.e usable) halo points
     logical :: mustputb, mustgetb, mustputf, mustgetf
     type(domain1D), pointer :: prev, next !neighbours toward decreasing and increasing indices
     type(domain1D), pointer, dimension(:) :: list
     integer :: pos             !position of this PE within list, i.e domain%list(pos)%pe = pe
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type, public :: domain2D
!     private
     type(domain1D) :: x
     type(domain1D) :: y
     type(domain2D), pointer, dimension(:) :: list
     integer :: pos             !position of this PE within list, i.e domain%list(pos)%pe = pe
     integer :: pe!OBSOLETE PE to which this domain is assigned
     integer :: whalo, ehalo, shalo, nhalo!OBSOLETE active halo widths
     type(domain2D), pointer :: west, east, south, north!OBSOLETE neighbours: west and south are toward decreasing indices
  end type domain2D
  type(domain1D), target, public :: NULL_DOMAIN1D !OBSOLETE
  type(domain2D), target, public :: NULL_DOMAIN2D !OBSOLETE

!parameters used to define domains: these are passed to the flags argument of mpp_define_domains
!  if compute domain is to have global extent, set GLOBAL_COMPUTE_DOMAIN
!  if data domain is to have global extent, set GLOBAL_DATA_DOMAIN
!  if global domain has periodic boundaries, set CYCLIC_GLOBAL_DOMAIN
!  sum flags together if more than one of the above conditions is to be met.
  integer, parameter, public :: GLOBAL_COMPUTE_DOMAIN=1, GLOBAL_DATA_DOMAIN=2, CYCLIC_GLOBAL_DOMAIN=4
  integer, parameter, public :: FOLD_WEST_EDGE =8, FOLD_EAST_EDGE =16
  integer, parameter, public :: FOLD_SOUTH_EDGE=8, FOLD_NORTH_EDGE=16
  integer, parameter, public :: WUPDATE=1, EUPDATE=2, SUPDATE=4, NUPDATE=8, TRANSPOSE=16
  integer, parameter, public :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE
  integer, parameter, public :: FLIP_AT_FOLD=32

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
  integer, private :: mpp_domains_stack_size=0

!used by mpp_define_domains2D_new to transmit data
  integer, allocatable, dimension(:) :: local_dom, remote_dom

!public interfaces
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
     module procedure mpp_define_domains1D_new
     module procedure mpp_define_domains2D_new
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

  interface mpp_update_domains_new
     module procedure mpp_update_domain2D_r8_2d_new
     module procedure mpp_update_domain2D_r8_3d_new
     module procedure mpp_update_domain2D_r8_4d_new
     module procedure mpp_update_domain2D_r8_5d_new
     module procedure mpp_update_domain2D_c8_2d_new
     module procedure mpp_update_domain2D_c8_3d_new
     module procedure mpp_update_domain2D_c8_4d_new
     module procedure mpp_update_domain2D_c8_5d_new
#ifndef no_8byte_integers
     module procedure mpp_update_domain2D_i8_2d_new
     module procedure mpp_update_domain2D_i8_3d_new
     module procedure mpp_update_domain2D_i8_4d_new
     module procedure mpp_update_domain2D_i8_5d_new
     module procedure mpp_update_domain2D_l8_2d_new
     module procedure mpp_update_domain2D_l8_3d_new
     module procedure mpp_update_domain2D_l8_4d_new
     module procedure mpp_update_domain2D_l8_5d_new
#endif
     module procedure mpp_update_domain2D_r4_2d_new
     module procedure mpp_update_domain2D_r4_3d_new
     module procedure mpp_update_domain2D_r4_4d_new
     module procedure mpp_update_domain2D_r4_5d_new
     module procedure mpp_update_domain2D_c4_2d_new
     module procedure mpp_update_domain2D_c4_3d_new
     module procedure mpp_update_domain2D_c4_4d_new
     module procedure mpp_update_domain2D_c4_5d_new
     module procedure mpp_update_domain2D_i4_2d_new
     module procedure mpp_update_domain2D_i4_3d_new
     module procedure mpp_update_domain2D_i4_4d_new
     module procedure mpp_update_domain2D_i4_5d_new
     module procedure mpp_update_domain2D_l4_2d_new
     module procedure mpp_update_domain2D_l4_3d_new
     module procedure mpp_update_domain2D_l4_4d_new
     module procedure mpp_update_domain2D_l4_5d_new

     module procedure mpp_update_domain1D_r8_2d_new
     module procedure mpp_update_domain1D_r8_3d_new
     module procedure mpp_update_domain1D_r8_4d_new
     module procedure mpp_update_domain1D_r8_5d_new
     module procedure mpp_update_domain1D_c8_2d_new
     module procedure mpp_update_domain1D_c8_3d_new
     module procedure mpp_update_domain1D_c8_4d_new
     module procedure mpp_update_domain1D_c8_5d_new
#ifndef no_8byte_integers
     module procedure mpp_update_domain1D_i8_2d_new
     module procedure mpp_update_domain1D_i8_3d_new
     module procedure mpp_update_domain1D_i8_4d_new
     module procedure mpp_update_domain1D_i8_5d_new
     module procedure mpp_update_domain1D_l8_2d_new
     module procedure mpp_update_domain1D_l8_3d_new
     module procedure mpp_update_domain1D_l8_4d_new
     module procedure mpp_update_domain1D_l8_5d_new
#endif
     module procedure mpp_update_domain1D_r4_2d_new
     module procedure mpp_update_domain1D_r4_3d_new
     module procedure mpp_update_domain1D_r4_4d_new
     module procedure mpp_update_domain1D_r4_5d_new
     module procedure mpp_update_domain1D_c4_2d_new
     module procedure mpp_update_domain1D_c4_3d_new
     module procedure mpp_update_domain1D_c4_4d_new
     module procedure mpp_update_domain1D_c4_5d_new
     module procedure mpp_update_domain1D_i4_2d_new
     module procedure mpp_update_domain1D_i4_3d_new
     module procedure mpp_update_domain1D_i4_4d_new
     module procedure mpp_update_domain1D_i4_5d_new
     module procedure mpp_update_domain1D_l4_2d_new
     module procedure mpp_update_domain1D_l4_3d_new
     module procedure mpp_update_domain1D_l4_4d_new
     module procedure mpp_update_domain1D_l4_5d_new
  end interface
  interface mpp_get_global
     module procedure mpp_get_global2d_r8_2d
     module procedure mpp_get_global2d_r8_3d
     module procedure mpp_get_global2d_r8_4d
     module procedure mpp_get_global2d_r8_5d
     module procedure mpp_get_global2d_c8_2d
     module procedure mpp_get_global2d_c8_3d
     module procedure mpp_get_global2d_c8_4d
     module procedure mpp_get_global2d_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_get_global2d_i8_2d
     module procedure mpp_get_global2d_i8_3d
     module procedure mpp_get_global2d_i8_4d
     module procedure mpp_get_global2d_i8_5d
#endif
     module procedure mpp_get_global2d_r4_2d
     module procedure mpp_get_global2d_r4_3d
     module procedure mpp_get_global2d_r4_4d
     module procedure mpp_get_global2d_r4_5d
     module procedure mpp_get_global2d_c4_2d
     module procedure mpp_get_global2d_c4_3d
     module procedure mpp_get_global2d_c4_4d
     module procedure mpp_get_global2d_c4_5d
     module procedure mpp_get_global2d_i4_2d
     module procedure mpp_get_global2d_i4_3d
     module procedure mpp_get_global2d_i4_4d
     module procedure mpp_get_global2d_i4_5d
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
  interface mpp_get_layout
     module procedure mpp_get_layout2D
  end interface
  public :: mpp_define_domains, mpp_domains_init, mpp_domains_set_stack_size, mpp_domains_exit, mpp_get_global, &
            mpp_update_domains, mpp_update_domains_new, operator(.EQ.), operator(.NE.), &
            mpp_get_active_domain, mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain, mpp_get_layout, &
            mpp_set_active_domain

  contains

    subroutine mpp_domains_init(flags,stackmax)
!initialize domain decomp package
      integer, intent(in), optional :: flags, stackmax

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

      if( PRESENT(stackmax) )then
          call mpp_domains_set_stack_size(stackmax)
      else
          call mpp_domains_set_stack_size(32768) !default, pretty arbitrary
      end if

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      NULL_DOMAIN1D%global%start_index  = -1; NULL_DOMAIN1D%global%end_index  = -1
      NULL_DOMAIN1D%data%start_index    = -1; NULL_DOMAIN1D%data%end_index    = -1
      NULL_DOMAIN1D%compute%start_index = -1; NULL_DOMAIN1D%compute%end_index = -1
      NULL_DOMAIN1D%pe = NULL_PE
      NULL_DOMAIN2D%x = NULL_DOMAIN1D
      NULL_DOMAIN2D%y = NULL_DOMAIN1D
      NULL_DOMAIN2D%pe = NULL_PE
      
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
      mpp_domains_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

    function mpp_domain1D_eq( a, b )
      logical :: mpp_domain1D_eq
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_eq = .FALSE.
      if( a%compute%start_index.EQ.b%compute%start_index .AND. &
          a%compute%end_index  .EQ.b%compute%end_index   .AND. &
          a%data%start_index   .EQ.b%data%start_index    .AND. &
          a%data%end_index     .EQ.b%data%end_index      .AND. & 
          a%global%start_index .EQ.b%global%start_index  .AND. &
          a%global%end_index   .EQ.b%global%end_index )mpp_domain1D_eq = .TRUE.
      return
    end function mpp_domain1D_eq

    function mpp_domain1D_ne( a, b )
      logical :: mpp_domain1D_ne
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_ne = .NOT. mpp_domain1D_eq( a, b )
      return
    end function mpp_domain1D_ne

    function mpp_domain2D_eq( a, b )
      logical :: mpp_domain2D_eq
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_eq = mpp_domain1D_eq( a%x, b%x ) .AND. mpp_domain1D_eq( a%y, b%y )
      return
    end function mpp_domain2D_eq

    function mpp_domain2D_ne( a, b )
      logical :: mpp_domain2D_ne
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_ne = .NOT. mpp_domain2D_eq( a, b )
      return
    end function mpp_domain2D_ne

    subroutine mpp_define_domains1D( global_indices, domain, pelist, flags, halo, extent )
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
      type(domain1D), intent(out), target :: domain(0:)
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: flags, halo
      integer, intent(in), optional :: extent(0:)

      logical :: compute_domain_is_global, data_domain_is_global, global_domain_is_cyclic
      integer :: numdomains, numpes, n, isg, ieg, is, ie, i
      
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )

!get global indices
      isg = global_indices(1)
      ieg = global_indices(2)
      if( size(domain).GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more domains requested than rows available.' )

!get flags
      compute_domain_is_global = .FALSE.
      data_domain_is_global    = .FALSE.
      global_domain_is_cyclic  = .FALSE.
      if( PRESENT(flags) )then
          compute_domain_is_global = BTEST(flags,0)
          data_domain_is_global    = BTEST(flags,1)
          global_domain_is_cyclic  = BTEST(flags,2)
      end if
!if compute domain is global, data domain must also be
      if( compute_domain_is_global )data_domain_is_global = .TRUE.

      numdomains = size(domain)
      numpes = npes
      if( PRESENT(pelist) )numpes = size(pelist)
      if( numdomains.EQ.1 )then
          compute_domain_is_global = .TRUE.
          data_domain_is_global    = .TRUE.
      end if

      domain(:)%global%start_index = isg
      domain(:)%global%end_index   = ieg
      domain(:)%global%size        = ieg-isg+1
      domain(:)%global%max_size    = ieg-isg+1
      domain(:)%global%is_global = .TRUE. !always
      domain(:)%data%is_global = data_domain_is_global
      domain(:)%compute%is_global = compute_domain_is_global
      domain(:)%ndomains = numdomains

!get compute domain
      if( compute_domain_is_global )then
          domain(:)%compute%start_index = isg
          domain(:)%compute%end_index   = ieg
      else
          is = isg
          do n=0,numdomains-1
             if( PRESENT(extent) )then
                 ie = is + extent(n) - 1
             else
                 ie = is + CEILING( float(ieg-is+1)/(numdomains-n) ) - 1
             end if
             domain(n)%compute%start_index = is
             domain(n)%compute%end_index   = ie
             is = ie + 1
          end do
      end if
      if( domain(numdomains-1)%compute%end_index.NE.ieg ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: problem with domain decomposition.' )
      domain(:)%compute%size = domain(:)%compute%end_index - domain(:)%compute%start_index + 1
!find max_size
      domain(:)%compute%max_size = MAXVAL( domain(:)%compute%size )
!create and fill in the sizelist array
      do n = 0,numdomains-1
         allocate( domain(n)%sizelist(numdomains) )
         domain(n)%sizelist(:) = domain(:)%compute%size
      end do

!get data domain
!data domain is at least equal to compute domain
      domain(:)%data%start_index = domain(:)%compute%start_index
      domain(:)%data%end_index   = domain(:)%compute%end_index
!apply global flags
      if( data_domain_is_global )then
          domain(:)%data%start_index  = isg
          domain(:)%data%end_index    = ieg
      end if
!apply margins
      if( PRESENT(halo) )then
          domain(:)%data%start_index = domain(:)%data%start_index - halo
          domain(:)%data%end_index   = domain(:)%data%end_index   + halo  
      end if
      domain(:)%data%size = domain(:)%data%end_index - domain(:)%data%start_index + 1
!find max_size
      domain(:)%data%max_size = MAXVAL( domain(:)%data%size )

!assign to PEs: currently does round-robin allocation
      do n = 0,numdomains-1
         allocate( domain(n)%pelist(numpes) )
         if( PRESENT(pelist) )then
             domain(n)%pe = pelist( mod(n,numpes) )
             domain(n)%pelist = pelist
             domain(n)%pos = mod(n,numpes) + 1
         else
             domain(n)%pe = mod(n,npes)
             domain(n)%pelist = (/ (i,i=0,npes-1) /)
             domain(n)%pos = mod(n,npes) + 1
         end if
      end do
!create linked list of neighbours
      do n = 0,numdomains-1
         if( n.NE.0            )domain(n)%prev => domain(n-1)
         if( n.NE.numdomains-1 )domain(n)%next => domain(n+1)
      end do
      if( global_domain_is_cyclic )then
          domain(0)%prev => domain(numdomains-1)
          domain(numdomains-1)%next => domain(0)
      end if

      if( verbose )then
          do n = 0,numdomains-1
             write(stdout,'(a,i2,a,3(x,4i3),a,i2)')'DD on PE', pe, ' compute, data, global=',&
                  domain(n)%compute%start_index, domain(n)%compute%end_index, domain(n)%compute%size, domain(n)%compute%max_size, &
                  domain(n)%data%start_index,    domain(n)%data%end_index,    domain(n)%data%size,    domain(n)%data%max_size,    &
                  domain(n)%global%start_index,  domain(n)%global%end_index,  domain(n)%global%size,  domain(n)%global%max_size,  &
                  ' pe(n)=', domain(n)%pe
          end do
      end if
      
      return
    end subroutine mpp_define_domains1D

    subroutine mpp_define_domains2D( global_indices, domain, pelist, xflags, yflags, xhalo, yhalo, xextent, yextent, &
                                     domain_layout, pe_layout )
!define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      type(domain2D), intent(out), target :: domain(0:)
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
      integer, intent(in), optional :: xextent(:), yextent(:)
      integer, intent(in), optional :: domain_layout(2), pe_layout(2)
      integer :: i, j, is, js, ie, je, isd, jsd, ied, jed, isg, jsg, ieg, jeg, npe, isz, jsz
      integer :: n, numdomains, ndomx, ndomy, numpes, npex, npey, pe_start
      logical :: x_compute_domain_is_global, y_compute_domain_is_global
      logical :: x_global_domain_is_cyclic, y_global_domain_is_cyclic

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )

      isg = global_indices(1)
      ieg = global_indices(2)
      jsg = global_indices(3)
      jeg = global_indices(4)

      isz = ieg - isg + 1
      jsz = jeg - jsg + 1

      numdomains = size(domain)

!get flags: no need to examine data_domain_is_global, that is passed directly to mpp_define_domains1d
      x_compute_domain_is_global = .FALSE.
      x_global_domain_is_cyclic  = .FALSE.
      if( PRESENT(xflags) )then
          x_compute_domain_is_global = BTEST(xflags,0)
          x_global_domain_is_cyclic  = BTEST(xflags,2)
      end if
      y_compute_domain_is_global = .FALSE.
      y_global_domain_is_cyclic  = .FALSE.
      if( PRESENT(yflags) )then
          y_compute_domain_is_global = BTEST(yflags,0)
          y_global_domain_is_cyclic  = BTEST(yflags,2)
      end if

!determine domain layout in X and Y: guarantee ndomx*ndomy=numdomains
      if( PRESENT(domain_layout) )then
          if( .NOT.PRESENT(pe_layout) ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: if you specify domain_layout, you must specify pe_layout as well.' )
          ndomx = domain_layout(1)
          ndomy = domain_layout(2)
          if( ndomx*ndomy.NE.numdomains ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain_layout(1)*domain_layout(2) must equal size(domain).'  )
      else                      !if not passed via domain_layout then
!first try to divide numdomains in the domain aspect ratio: if imperfect aspect, reduce ndomx till it divides numdomains
          ndomx = nint( sqrt(float(numdomains*isz)/jsz) )
          ndomx = max( ndomx,1 ) !for isz=1 line above can give 0
          do while( mod(numdomains,ndomx).NE.0 )
             ndomx = ndomx - 1
          end do                 !will terminate at ndomx=1 if not before
          ndomy = numdomains/ndomx
!then apply global flags
          if( x_compute_domain_is_global )then
              ndomx = 1
              ndomy = numdomains
!if numdomains>1 and compute domain is both x_ and y_global, we assume full redundancy
              if( y_compute_domain_is_global )ndomx = numdomains
          else if( y_compute_domain_is_global )then     !y is global but not x
              ndomy = 1
              ndomx = numdomains
          end if
      end if
!do the same as above for PE count
      if( PRESENT(pe_layout) )then
          npex = pe_layout(1)
          npey = pe_layout(2)
          numpes = npex*npey
      else
          numpes = npes
          if( PRESENT(pelist) )numpes = size(pelist)
          npex = nint( sqrt(float(numpes*isz)/jsz) )
          npex = max(npex,1)
          do while( mod(numpes,npex).NE.0 )
             npex = npex - 1
          end do
          npey = numpes/npex
          if( x_compute_domain_is_global )then
              npex = 1
              npey = numpes
              if( y_compute_domain_is_global )npex = numpes
          else if( y_compute_domain_is_global )then
              npey = 1
              npex = numpes
          end if
      end if
!we now have ndomx, ndomy, npex, npey
      if( PRESENT(pelist) )then
          do n = 0,ndomy-1
             pe_start = pelist( mod(n,npey)*npex )
             call mpp_define_domains1D( (/isg,ieg/), domain(ndomx*n:ndomx*(n+1)-1)%x, &
                                                     pelist(pe_start:pe_start+npex-1), xflags, xhalo, xextent )
          end do
          do n = 0,ndomx-1
             pe_start = pelist( mod(n,npex) )
             call mpp_define_domains1D( (/jsg,jeg/), domain(n:numdomains-1:ndomx)%y, &
                                                     pelist(pe_start:numpes-1:npex),   yflags, yhalo, yextent )
          end do
      else
          do n = 0,ndomy-1
             pe_start = mod(n,npey)*npex
             call mpp_define_domains1D( (/isg,ieg/), domain(ndomx*n:ndomx*(n+1)-1)%x, &
                                                     (/ (i,i=pe_start,pe_start+npex-1) /), xflags, xhalo, xextent )
          end do
          do n = 0,ndomx-1
             pe_start = mod(n,npex)
             call mpp_define_domains1D( (/jsg,jeg/), domain(n:numdomains-1:ndomx)%y, &
                                                     (/ (i,i=pe_start,numpes-1,npex) /),   yflags, yhalo, yextent )
          end do
      end if
!assign domains to PEs
      do n = 0,numdomains-1
         if( domain(n)%x%pe.NE.domain(n)%y%pe ) &
              call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: Problem with PE assignment to domain, x%pe.NE.y%pe.' )
         domain(n)%pe = domain(n)%x%pe
      end do
!create linked list of domains
      do n = 0,numdomains-1
         if( n.GE.ndomx            )domain(n)%south => domain(n-ndomx)
         if( n.LT.numdomains-ndomx )domain(n)%north => domain(n+ndomx)
         if( mod(n,ndomx).NE.0       )domain(n)%west => domain(n-1)
         if( mod(n,ndomx).NE.ndomx-1 )domain(n)%east => domain(n+1)
         if( x_global_domain_is_cyclic )then
             if( mod(n,ndomx).EQ.0       )domain(n)%west => domain(n+ndomx-1)
             if( mod(n,ndomx).EQ.ndomx-1 )domain(n)%east => domain(n-ndomx+1)
         end if
         if( y_global_domain_is_cyclic )then
             if( n.LT.ndomx            )domain(n)%south => domain(n+ndomx*(ndomy-1))
             if( n.GE.numdomains-ndomx )domain(n)%north => domain(n-ndomx*(ndomy-1))
         end if
      end do
          
      return
    end subroutine mpp_define_domains2D

    subroutine mpp_define_domains1D_new( global_indices, ndivs, domain, pelist, flags, halo, extent, mask )
!routine to divide global array indices among domains, and assign domains to PEs
!domain is of type domain1D
!ARGUMENTS:
!      global_indices(2)=(isg,ieg) gives the extent of global domain
!      ndivs is number of divisions of domain: even divisions unless extent is present.
!      domain is the returned domain1D
!      (optional) pelist list of PEs to which domains are to be assigned (default 0...npes-1)
!                 size of pelist must correspond to number of unmasked divisions
!      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
!      (optional) halo defines halo width (currently the same on both sides)
!      (optional) array extent defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
!      (optional) a division all whose points are masked are not assigned to any domain
!  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
!  as the input parameters permit
      integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
      integer, intent(in) :: ndivs
      type(domain1D), intent(inout) :: domain !declared inout so that existing list, if any, can be nullified
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: flags, halo
      integer, intent(in), optional :: extent(0:)
      logical, intent(in), optional :: mask(0:)

      logical :: compute_domain_is_global, data_domain_is_global, global_domain_is_cyclic
      integer :: ndomains, ndiv, n, isg, ieg, is, ie, i, off, pos
      integer, allocatable :: pes(:)
      logical, allocatable :: masked(:)
      
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
!get global indices
      isg = global_indices(1)
      ieg = global_indices(2)
      if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
!get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
      domain%pos = -1
      if( PRESENT(pelist) )then
          allocate( pes(0:size(pelist)-1) )
          pes = pelist
          do i = 0,size(pes)-1
             if( pe.EQ.pelist(i) )then
                 domain%pos = i
                 exit
             end if
          end do
          if( domain%pos.EQ.-1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
      else
          allocate( pes(0:npes-1) )
          pes = (/ (i,i=0,npes-1) /)
          domain%pos = pe
      end if
      pos = domain%pos
!the size of the pes array is the number of divisions that will be assigned to domains on PEs
      ndomains = size(pes)
!get number of real domains: 1 unmasked domain per PE in pes
      allocate( masked(0:ndivs-1) )
      masked = .FALSE.                 !default unmasked
      if( PRESENT(mask) )then
          if( size(mask).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: mask array size must equal number of domain divisions.' )
          masked(:) = mask(:)
      end if
      if( count(.NOT.masked).NE.ndomains ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: mask arrays contains the wrong number of unmasked domains.' )
      if( PRESENT(extent) )then
          if( size(extent).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
      end if

      if( ASSOCIATED(domain%list) )NULLIFY(domain%list)
      allocate( domain%list(0:ndomains-1) )

!get flags
      compute_domain_is_global = .FALSE.
      data_domain_is_global    = .FALSE.
      global_domain_is_cyclic  = .FALSE.
      if( PRESENT(flags) )then
          compute_domain_is_global = BTEST(flags,0) .OR. ndomains.EQ.1
!if compute domain is global, data domain must also be
          data_domain_is_global    = BTEST(flags,1) .OR. compute_domain_is_global
          global_domain_is_cyclic  = BTEST(flags,2)
      end if

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
      else
          is = isg
          n = 0
          do ndiv=0,ndivs-1
             if( PRESENT(extent) )then
                 ie = is + extent(n) - 1
             else
                 ie = is + CEILING( float(ieg-is+1)/(ndomains-n) ) - 1
             end if
             if( .NOT.masked(ndiv) )then
                 domain%list(n)%compute%begin = is
                 domain%list(n)%compute%end   = ie
                 domain%list(n)%compute%is_global = .FALSE.
                 n = n + 1
             end if
             is = ie + 1
          end do
          if( n.NE.ndomains )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: problem with domain decomposition.' )
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
      if( PRESENT(halo) )then
          domain%list(:)%data%begin = domain%list(:)%data%begin - halo
          domain%list(:)%data%end   = domain%list(:)%data%end   + halo  
      end if
      domain%list(:)%data%size = domain%list(:)%data%end - domain%list(:)%data%begin + 1

!      domain = domain%list(pos) !load domain from domain%list(pos)
      domain%compute = domain%list(pos)%compute
      domain%data = domain%list(pos)%data
      domain%global = domain%list(pos)%global
      domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
      domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
      domain%global%max_size  = domain%global%size
      domain%active = domain%data !active domain is initialized to data domain

!assign to PEs
      domain%list(:)%pe = pes(:)

!create list of put and get domains
      domain%list(:)%mustputb = .FALSE.
      domain%list(:)%mustgetb = .FALSE.
      domain%list(:)%mustputf = .FALSE.
      domain%list(:)%mustgetf = .FALSE.
      do n = 0,ndomains-1

         pos = domain%pos-n
         if( global_domain_is_cyclic )pos = mod(pos+ndomains,ndomains)
         if( pos.GE.0 .AND. pos.LT.ndomains )then
!putb: domain%list(pos) upper halo must lie in my compute domain
             off = 0
             if( global_domain_is_cyclic .AND. domain%list(pos)%compute%end.EQ.domain%list(pos)%global%end )off = ieg-isg+1
             domain%list(pos)%mustputb = if_overlap( domain%list(pos)%compute%end+1-off, domain%list(pos)%data%end-off, &
                  domain%compute%begin, domain%compute%end, &
                  domain%list(pos)%putb%begin, domain%list(pos)%putb%end )
             domain%list(pos)%putb%size = domain%list(pos)%putb%end - domain%list(pos)%putb%begin + 1
!getb: my lower halo must lie in domain%list(pos) compute domain
             off = 0
             if( global_domain_is_cyclic .AND. domain%compute%begin.EQ.domain%global%begin )off = ieg-isg+1
             domain%list(pos)%mustgetb = if_overlap( domain%data%begin, domain%compute%begin-1, &
                  domain%list(pos)%compute%begin-off, domain%list(pos)%compute%end-off, &
                  domain%list(pos)%getb%begin, domain%list(pos)%getb%end )
             domain%list(pos)%getb%size = domain%list(pos)%getb%end - domain%list(pos)%getb%begin + 1
             if( debug )write( stderr,* )'MPP_DEFINE_DOMAINS1D_NEW: pe, domain%pos, cyclic, n, pos, off=', &
                  pe, domain%pos, global_domain_is_cyclic, n, pos, off, &
                  ' putb=', domain%list(pos)%mustputb, domain%list(pos)%putb%begin, domain%list(pos)%putb%end, &
                  ' getb=', domain%list(pos)%mustgetb, domain%list(pos)%getb%begin, domain%list(pos)%getb%end
         end if

         pos = domain%pos+n
         if( global_domain_is_cyclic )pos = mod(pos,ndomains)
         if( pos.GE.0 .AND. pos.LT.ndomains )then
!putf: domain%list(pos) lower halo must lie in my compute domain
             off = 0
             if( global_domain_is_cyclic .AND. domain%list(pos)%compute%begin.EQ.domain%list(pos)%global%begin )off = ieg-isg+1
             domain%list(pos)%mustputf = if_overlap( domain%list(pos)%data%begin+off, domain%list(pos)%compute%begin-1+off, &
                  domain%compute%begin, domain%compute%end, &
                  domain%list(pos)%putf%begin, domain%list(pos)%putf%end )
             domain%list(pos)%putf%size = domain%list(pos)%putf%end - domain%list(pos)%putf%begin + 1
!getf: my upper halo must lie in domain%list(pos) compute domain
             off = 0
             if( global_domain_is_cyclic .AND. domain%compute%end.EQ.domain%global%end )off = ieg-isg+1
             domain%list(pos)%mustgetf = if_overlap( domain%compute%end+1, domain%data%end, &
                  domain%list(pos)%compute%begin+off, domain%list(pos)%compute%end+off, &
                  domain%list(pos)%getf%begin, domain%list(pos)%getf%end )
             domain%list(pos)%getf%size = domain%list(pos)%getf%end - domain%list(pos)%getf%begin + 1
             if( debug )write( stderr,* )'MPP_DEFINE_DOMAINS1D_NEW: pe, domain%pos, cyclic, n, pos, off=', &
                  pe, domain%pos, global_domain_is_cyclic, n, pos, off, &
                  ' putf=', domain%list(pos)%mustputf, domain%list(pos)%putf%begin, domain%list(pos)%putf%end, &
                  ' getf=', domain%list(pos)%mustgetf, domain%list(pos)%getf%begin, domain%list(pos)%getf%end
         end if
      end do
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
               'MPP_DEFINE_DOMAINS1D_NEW: pe, hs, he, cs, ce, os, oe=', pe, hs, he, cs, ce, os, oe
          if_overlap = oe.GE.os
          return
        end function if_overlap
          
    end subroutine mpp_define_domains1D_new

    subroutine mpp_define_domains2D_new( global_indices, layout, domain, pelist, &
         xflags, yflags, xhalo, yhalo, xextent, yextent, mask )
!define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      integer, intent(in) :: layout(2)
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
      integer, intent(in), optional :: xextent(0:), yextent(0:)
      logical, intent(in), optional :: mask(0:,0:)
      integer :: i, j, n, ipos, jpos, pos, remote_pos
      logical, allocatable :: masked(:,:)
      integer, allocatable :: pes(:), pearray(:,:)
      character(len=8) :: text
      logical :: w_edge_is_folded, e_edge_is_folded, s_edge_is_folded, n_edge_is_folded

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )

      if( PRESENT(pelist) )then
          allocate( pes(0:size(pelist)-1) )
          pes = pelist
      else
          allocate( pes(0:npes-1) )
          pes = (/ (i,i=0,npes-1) /)
      end if

      allocate( masked(0:layout(1)-1,0:layout(2)-1) )
      masked = .FALSE.
      if( PRESENT(mask) )then
          if( size(mask,1).NE.layout(1) .OR. size(mask,2).NE.layout(2) ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: mask array does not match layout.' )
          masked(:,:) = mask(:,:)
      end if
!number of unmasked domains in layout must equal number of PEs assigned
      n = count(.NOT.masked)
      if( n.NE.size(pes) )then
          write( text,'(i8)' )n
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for this layout and mask. Use ' &
               //text//' PEs for this domain decomposition.' )
      end if

!place on PE array; need flag to assign them to j first and then i
      allocate( pearray(0:layout(1)-1,0:layout(2)-1) )
      ipos = -1; jpos = -1; pos = -1
      n = 0
      do j = 0,layout(2)-1
         do i = 0,layout(1)-1
            if( .NOT.masked(i,j) )then
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
      call mpp_define_domains( global_indices(1:2), layout(1), domain%x, &
           pack(pearray(:,jpos),.NOT.masked(:,jpos)), xflags, xhalo, xextent, masked(:,jpos) )
      call mpp_define_domains( global_indices(3:4), layout(2), domain%y, &
           pack(pearray(ipos,:),.NOT.masked(ipos,:)), yflags, yhalo, yextent, masked(ipos,:) )
      if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
      domain%pos = pos
      if( debug )then
          write( stderr,'(a,4i4,a,(8i4))' )'MPP_DEFINE_DOMAINS2D_NEW: pe, pos, ipos, jpos=', pe, pos, ipos, jpos, &
               'pes=', pes
          write( stderr,'(a,5i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, bounds(xlist), bounds(ylist)=', pe, &
               lbound(domain%x%list), ubound(domain%x%list), lbound(domain%y%list), ubound(domain%y%list)
          do i = 0,size(domain%x%list(:))-1
             if( domain%x%list(i)%mustputb )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, i, putb=', &
                  pe, i, domain%x%list(i)%putb%begin, domain%x%list(i)%putb%end
             if( domain%x%list(i)%mustputf )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, i, putf=', &
                  pe, i, domain%x%list(i)%putf%begin, domain%x%list(i)%putf%end
             if( domain%x%list(i)%mustgetb )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, i, getb=', &
                  pe, i, domain%x%list(i)%getb%begin, domain%x%list(i)%getb%end
             if( domain%x%list(i)%mustgetf )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, i, getf=', &
                  pe, i, domain%x%list(i)%getf%begin, domain%x%list(i)%getf%end
          end do
          do i = 0,size(domain%y%list(:))-1
             if( domain%y%list(i)%mustputb )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, j, putb=', &
                  pe, i, domain%y%list(i)%putb%begin, domain%y%list(i)%putb%end
             if( domain%y%list(i)%mustputf )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, j, putf=', &
                  pe, i, domain%y%list(i)%putf%begin, domain%y%list(i)%putf%end
             if( domain%y%list(i)%mustgetb )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, j, getb=', &
                  pe, i, domain%y%list(i)%getb%begin, domain%y%list(i)%getb%end
             if( domain%y%list(i)%mustgetf )write( stderr,'(a,4i4)' )'MPP_DEFINE_DOMAINS2D_NEW: pe, j, getf=', &
                  pe, i, domain%y%list(i)%getf%begin, domain%y%list(i)%getf%end
          end do
      end if

      if( .NOT.allocated(local_dom) )allocate( local_dom(8) )
      if( .NOT.allocated(remote_dom) )allocate( remote_dom(8) )
      call mpp_get_compute_domain( domain, local_dom(1), local_dom(2), local_dom(3), local_dom(4) )
      call mpp_get_data_domain   ( domain, local_dom(5), local_dom(6), local_dom(7), local_dom(8) )
      if( debug )write( stderr,'(a,9i4)' )'pe, domain=', pe, local_dom
      n = size(pes)
      if( ASSOCIATED(domain%list) )NULLIFY(domain%list)
      allocate( domain%list(0:n-1) )
      if( pe.EQ.root ) &
           print *, 'Domain decomposition: pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
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
         if( pe.EQ.root )print '(22x,i3,x,4i5,3x,4i5)', pes(remote_pos), remote_dom
      end do
      call mpp_sync_self(pes)
          
      return
    end subroutine mpp_define_domains2D_new

    subroutine mpp_get_layout2D( global_indices, ndivs, layout )
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
    end subroutine mpp_get_layout2D

    subroutine mpp_get_compute_domain1D( domain, begin, end, size, max_size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size

      if( PRESENT(begin)    )begin    = domain%compute%begin
      if( PRESENT(end)      )end      = domain%compute%end
      if( PRESENT(size)     )size     = domain%compute%size
      if( PRESENT(max_size) )max_size = domain%compute%max_size
      return
    end subroutine mpp_get_compute_domain1D

    subroutine mpp_get_data_domain1D( domain, begin, end, size, max_size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size

      if( PRESENT(begin)    )begin    = domain%data%begin
      if( PRESENT(end)      )end      = domain%data%end
      if( PRESENT(size)     )size     = domain%data%size
      if( PRESENT(max_size) )max_size = domain%data%max_size
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

    subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_compute_domain1D( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_compute_domain1D( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_compute_domain2D

    subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_data_domain1D( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_data_domain1D( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_data_domain2D

    subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_global_domain1D( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_global_domain1D( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_global_domain2D

    subroutine mpp_get_active_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_active_domain1D( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_active_domain1D( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_active_domain2D

    subroutine mpp_set_active_domain2D( domain, xbegin, xend, ybegin, yend )
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: xbegin, xend, ybegin, yend
      call mpp_set_active_domain1D( domain%x, xbegin, xend )
      call mpp_set_active_domain1D( domain%y, ybegin, yend )
      return
    end subroutine mpp_set_active_domain2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GET_GLOBAL: get global field from domain field             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GET_GLOBAL_ mpp_get_global2D_r8_2d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r8_3d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r8_4d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c8_2d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c8_3d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c8_4d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>

#ifndef no_8byte_integers
#define MPP_GET_GLOBAL_ mpp_get_global2D_i8_2d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i8_3d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i8_4d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>
#endif

#define MPP_GET_GLOBAL_ mpp_get_global2D_r4_2d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r4_3d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r4_4d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c4_2d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c4_3d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c4_4d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_c4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i4_2d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_EXTRA_INDICES_
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i4_3d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_EXTRA_INDICES_ ,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i4_4d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:
#include <mpp_get_global.h>

#define MPP_GET_GLOBAL_ mpp_get_global2D_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_EXTRA_INDICES_ ,:,:,:
#include <mpp_get_global.h>

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
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_update_domains2D.h>

#ifndef no_8byte_integers
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l8_5D
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_update_domains2D.h>
#endif

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_update_domains2D.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l4_5D
#define MPP_TYPE_ logical(INT_KIND)
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

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r8_5D_new
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D_new
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_update_domains2D_new.h>

#ifndef no_8byte_integers
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D_new
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l8_5D_new
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_update_domains2D_new.h>
#endif

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D_new
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D_new
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D_new
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l4_5D_new
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_update_domains2D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r8_5D_new
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_update_domains1D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c8_5D_new
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_update_domains1D_new.h>

#ifndef no_8byte_integers
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i8_5D_new
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_update_domains1D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l8_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l8_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l8_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l8_5D_new
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_update_domains1D_new.h>
#endif

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r4_5D_new
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_update_domains1D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c4_5D_new
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_update_domains1D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i4_5D_new
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_update_domains1D_new.h>

#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l4_2D_new
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l4_3D_new
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l4_4D_new
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l4_5D_new
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_update_domains1D_new.h>

    subroutine get_halos_1D( domain, put_domain, get_domain, sense, put_pe, get_pe, put_start, put_end, get_start, get_end )
!compute halo zones from given domain1D types
!sense=-1 for left halo, +1 for right halo
!returns put_pe, put_start and put_end within domain%compute, to be sent to put_domain's halo
!    and get_pe, get_start and get_end within domain's halo, to be got from get_domain%compute
      type(domain1D), intent(in) :: domain, put_domain, get_domain
      integer, intent(in) :: sense
      integer, intent(out) :: put_pe, get_pe, put_start, put_end, get_start, get_end
      logical :: wrapped

      if( is_halo(domain,put_domain,sense,put_start,put_end,wrapped) )then
          put_pe = put_domain%pe
      else
          put_pe = NULL_PE
      end if
      if( is_halo(get_domain,domain,sense,get_start,get_end,wrapped) )then
          get_pe = get_domain%pe
!is_halo returns numbers relative to get_domain: if wrapped we need to shift them to get numbers relative to <domain>
          if( wrapped )then
              get_start = get_start + sense*domain%global%size
              get_end   = get_end   + sense*domain%global%size
          end if
      else
          get_pe = NULL_PE
      end if
              
      return
    end subroutine get_halos_1D

    function is_halo( donor, receiver, sense, halo_start, halo_end, wrapped )
!returns true if a portion of receiver's halo (returned in [halo_start,halo_end]) lies within donor's compute domain
!else returns false
!this routine is aware of cyclic domain wraparound, and will make appropriate adjustments
![halo_start,halo_end] always lies within donor's compute zone
!if wrapped returns true, the numbers relative to receiver domain will have to be shifted by sense*global%size
      logical :: is_halo
      type(domain1D), intent(in) :: donor, receiver
      integer, intent(in) :: sense
      integer, intent(out) :: halo_start, halo_end
      logical, intent(out) :: wrapped
      integer :: h_start, h_end, d_start, d_end, g_start, g_end, g_size

      if( sense.NE.1 .AND. sense.NE.-1 )call mpp_error( FATAL, 'IS_HALO: sense must be +1 or -1.' )
      g_start = donor%global%start_index; g_end = donor%global%end_index; g_size = donor%global%size !shorthand
      d_start = donor%compute%start_index; d_end = donor%compute%end_index !shorthand
      wrapped = .FALSE.
      if( receiver.EQ.NULL_DOMAIN1D )then
          is_halo = .FALSE.
      else if( sense.EQ.-1 )then
          h_start = receiver%data%start_index; h_end = receiver%compute%start_index-1 !left halo of receiver
          if( h_start.LT.g_start )then !run off left edge of cyclic domain
              if( h_end.GE.g_start )then !straddling left edge, first check inside, then outside
                  is_halo = is_overlap(d_start,d_end,g_start,h_end,halo_start,halo_end)
                  if( .NOT.is_halo )then
                      is_halo = is_overlap(d_start,d_end,h_start+g_size,g_start-1+g_size,halo_start,halo_end)
                      wrapped = .TRUE.
                  end if
              else              !entirely off left edge
                  is_halo = is_overlap(d_start,d_end,h_start+g_size,h_end+g_size,halo_start,halo_end)
                  wrapped = .TRUE.
              end if
          else                  !entirely within global_domain
              is_halo = is_overlap(d_start,d_end,h_start,h_end,halo_start,halo_end)
          end if
      else                      !sense=+1
          h_start = receiver%compute%end_index+1; h_end = receiver%data%end_index !right halo of receiver
          if( h_end.GT.g_end )then !run off right edge of cyclic domain
              if( h_start.LE.g_end )then !straddling right edge, first check inside, then outside
                  is_halo = is_overlap(d_start,d_end,h_start,g_end,halo_start,halo_end)
                  if( .NOT.is_halo )then
                      is_halo = is_overlap(d_start,d_end,g_end+1-g_size,h_end-g_size,halo_start,halo_end)
                      wrapped = .TRUE.
                  end if
              else              !entirely off right edge
                  is_halo = is_overlap(d_start,d_end,h_start-g_size,h_end-g_size,halo_start,halo_end)
                  wrapped = .TRUE.
              end if
          else                  !entirely within global_domain
              is_halo = is_overlap(d_start,d_end,h_start,h_end,halo_start,halo_end)
          end if
      end if

      return
    end function is_halo

    function is_overlap( ds, de, hs, he, is, ie )
!if a portion of [hs,he] lies within [ds,de], return true, AND return the intersection in [is,ie]
!if there is no intersection, return false
!assumes hs.LE.he and ds.LE.de: possible problem otherwise
      logical :: is_overlap
      integer, intent(in) :: ds, de, hs, he
      integer, intent(out) :: is, ie

      is = max(hs,ds)
      ie = min(he,de)
      is_overlap = is.LE.de .AND. ie.GE.ds
      return
    end function is_overlap

  end module mpp_domains_mod

#ifdef test_mpp_domains
  program mpp_domains_test
    use mpp_mod
    use mpp_domains_mod
    implicit none
    integer :: pe, npes, root
    type(domain2D), allocatable, target :: domains(:)
    integer :: nx=128, ny=128, nz=40, halo=2, stackmax=32768
    real(DOUBLE_KIND), dimension(:,:,:), allocatable :: local, global, gglobal
    integer :: is, ie, js, je, isd, ied, jsd, jed
    integer :: tk, tk0, tks_per_sec
    integer :: i,j,k, unit=7, layout(2)
    real :: t
    logical :: new=.FALSE., debug=.FALSE., opened
    namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, new, debug

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
    call mpp_set_stack_size(stackmax)

    call SYSTEM_CLOCK( count_rate=tks_per_sec )
    if( debug )then
        call mpp_domains_init(MPP_DEBUG)
    else
        call mpp_domains_init()
    end if

    allocate( domains(0:npes-1) )
    if( new )then
        call mpp_get_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domains(pe), xhalo=halo, yhalo=halo )
        call mpp_get_compute_domain( domains(pe), is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domains(pe), isd, ied, jsd, jed )
    else
        call mpp_define_domains( (/1,nx,1,ny/), domains, xhalo=halo, yhalo=halo )
        if( pe.EQ.root )then
            print '(a,5i4)', 'ndomains, nx, ny, nz, halo=', size(domains), nx, ny, nz, halo
            print '(a)', 'DOMAIN DECOMPOSITION:'
            print *, ' pe,  is,  ie,  js,  je,  isd,  ied,  jsd,  jed'
            do i = 0,size(domains)-1
               print '(i4,4i5,4i6)', i, &
                    domains(i)%x%compute%start_index, &
                    domains(i)%x%compute%end_index,   &
                    domains(i)%y%compute%start_index, &
                    domains(i)%y%compute%end_index,   &
                    domains(i)%x%data%start_index, &
                    domains(i)%x%data%end_index,   &
                    domains(i)%y%data%start_index, &
                    domains(i)%y%data%end_index
            end do
        end if
        is = domains(pe)%x%compute%start_index
        ie = domains(pe)%x%compute%end_index
        js = domains(pe)%y%compute%start_index
        je = domains(pe)%y%compute%end_index
        isd = domains(pe)%x%data%start_index
        ied = domains(pe)%x%data%end_index
        jsd = domains(pe)%y%data%start_index
        jed = domains(pe)%y%data%end_index
    end if
       
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
    if( new )then
        call mpp_update_domains_new( local, domains(pe) )
    else
        call mpp_update_domains( local, domains(pe) )
    end if
!compare checksums between global and local arrays
    call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Halo update' )

!fill in gglobal array
    if( new )then
    else
        call mpp_get_global( domains(pe), local, gglobal )
!compare checksums between global and local arrays
        call compare_checksums( global(1:nx,1:ny,:), gglobal, 'mpp_get_global' )
    end if

!test cyclic boundary conditions
    if( new )then
        call mpp_define_domains( (/1,nx,1,ny/), layout, domains(pe), xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
        call mpp_get_compute_domain( domains(pe), is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domains(pe), isd, ied, jsd, jed )
    else
        call mpp_define_domains( (/1,nx,1,ny/), domains, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
        is = domains(pe)%x%compute%start_index
        ie = domains(pe)%x%compute%end_index
        js = domains(pe)%y%compute%start_index
        je = domains(pe)%y%compute%end_index
        isd = domains(pe)%x%data%start_index
        ied = domains(pe)%x%data%end_index
        jsd = domains(pe)%y%data%start_index
        jed = domains(pe)%y%data%end_index
    end if
!fill in cyclic global array
    global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
    global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
    global(1-halo:nx+halo,    1-halo:0,:) = global(1-halo:nx+halo,ny-halo+1:ny,:)
    global(1-halo:nx+halo,ny+1:ny+halo,:) = global(1-halo:nx+halo,1:halo,      :)
!fill in local array
    local = 0.
    local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
    if( new )then
        call mpp_update_domains_new( local, domains(pe) )
    else
        call mpp_update_domains( local, domains(pe) )
    end if
!compare checksums between global and local arrays
    call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Cyclic halo update' )

!timing tests
    if( pe.EQ.root )print '(a)', 'TIMING TESTS:'
    if( new )then
        call mpp_define_domains( (/1,nx,1,ny/), layout, domains(pe), xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
        call mpp_get_compute_domain( domains(pe), is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domains(pe), isd, ied, jsd, jed )
    else
        call mpp_define_domains( (/1,nx,1,ny/), domains, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
        is = domains(pe)%x%compute%start_index
        ie = domains(pe)%x%compute%end_index
        js = domains(pe)%y%compute%start_index
        je = domains(pe)%y%compute%end_index
        isd = domains(pe)%x%data%start_index
        ied = domains(pe)%x%data%end_index
        jsd = domains(pe)%y%data%start_index
        jed = domains(pe)%y%data%end_index
    end if
    if( ALLOCATED(local) )deallocate(local)
    allocate( local(isd:ied,jsd:jed,nz) )
!fill in local array
    local = 0.
    local(is:ie,js:je,:) = global(is:ie,js:je,:)
!fill in halos
    call mpp_sync()          !this ensures you time only the update_domains call
    call SYSTEM_CLOCK(tk0)
    if( new )then
        call mpp_update_domains_new( local, domains(pe) )
    else
        call mpp_update_domains( local, domains(pe) )
    end if
    call SYSTEM_CLOCK(tk)
    t = float(tk-tk0)/tks_per_sec
!compare checksums between global and local arrays
    call compare_checksums( local, global(isd:ied,jsd:jed,:), 'Timed cyclic halo update' )
!words transferred
    j = ( (ied-isd+1)*(jed-jsd+1) - (ie-is+1)*(je-js+1) )*nz
    call mpp_max(j)
    if( pe.EQ.root ) &
         print '(a,i4,i8,es12.4,f10.3)', 'Halo width, words, time(s), bandwidth (MB/s)=', halo, j, t, j*8e-6/t

    call mpp_domains_exit()
    call mpp_exit()

    contains
      subroutine compare_checksums( a, b, string )
        real(DOUBLE_KIND), intent(in), dimension(:,:,:) :: a, b
        character(len=*), intent(in) :: string
        integer :: i, j
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
