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
!these are used to determine hardware/OS/compiler

#ifdef __sgi
#ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!compilers that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)
#define use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#endif

#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints
!at any price
#define LONG_KIND 4
#endif

!if using shmem calls on Origin, you may need to use shmalloc
#ifdef use_libSMA
#ifdef sgi_mipspro
#define use_shmalloc
#endif
#endif

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999
  use mpp_mod
  implicit none
  private
  character(len=256), private :: version='$Id: mpp_domains.F90,v 5.3 2000/05/30 19:46:56 vb Exp $'

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5,   stdout=6,   stderr=0
#endif
  integer, private :: pe, npes

  type, public :: domain_axis_spec        !type used to specify index limits along an axis of a domain
     sequence
     integer :: start_index, end_index, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec
  type, public :: domain1D
     sequence
     type(domain_axis_spec) :: compute, data, global
     integer :: ndomains        !number of domains over which data has been distributed
     integer :: pe              !PE to which this domain is assigned
     integer, pointer :: pelist(:), sizelist(:)
     integer :: pos             !position of this PE within pelist, i.e pelist(pos) = pe
     integer :: lhalo, rhalo    !can be used to keep track of active (i.e usable) halo points
     type(domain1D), pointer :: prev, next !neighbours toward decreasing and increasing indices
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type, public :: domain2D
     sequence
     type(domain1D) :: x
     type(domain1D) :: y
     integer :: pe              !PE to which this domain is assigned
     integer :: whalo, ehalo, shalo, nhalo !active halo widths
     type(domain2D), pointer :: west, east, south, north !neighbours: west and south are toward decreasing indices
  end type domain2D
  type(domain1D), target, public :: NULL_DOMAIN1D
  type(domain2D), target, public :: NULL_DOMAIN2D

!parameters used to define domains: these are passed to the flags argument of mpp_define_domains
!  if compute domain is to have global extent, set GLOBAL_COMPUTE_DOMAIN
!  if data domain is to have global extent, set GLOBAL_DATA_DOMAIN
!  if global domain has periodic boundaries, set CYCLIC_GLOBAL_DOMAIN
!  sum flags together if more than one of the above conditions is to be met.
  integer, parameter, public :: GLOBAL_COMPUTE_DOMAIN=1, GLOBAL_DATA_DOMAIN=2, CYCLIC_GLOBAL_DOMAIN=4
  integer, parameter, public :: WUPDATE=1, EUPDATE=2, SUPDATE=4, NUPDATE=8, TRANSPOSE=16
  integer, parameter, public :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE

!buffer area for halo regions
  integer, private :: max_halo_size=0
#ifdef use_shmalloc
!put_r8, etc will be pointed into shared_heap
  real :: shared_heap(1)
  integer :: len_heap=0
  pointer( ptr_heap, shared_heap )
#else
  real, allocatable :: put_r8(:), get_r8(:)
  complex, allocatable :: put_c8(:), get_c8(:)
#endif

  integer, private :: tick
  logical, private :: verbose=.FALSE.
  logical, private :: mpp_domains_initialized=.FALSE.

!public interfaces
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface
  interface mpp_update_domains
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain1D_i8_1d
     module procedure mpp_update_domain1D_i8_2d
     module procedure mpp_update_domain1D_i8_3d
     module procedure mpp_update_domain1D_r8_1d
     module procedure mpp_update_domain1D_r8_2d
     module procedure mpp_update_domain1D_r8_3d
  end interface
  interface mpp_get_global
     module procedure mpp_get_global2d_r8_3d
     module procedure mpp_get_global2d_r8_2d
     module procedure mpp_get_global2d_i8_3d
     module procedure mpp_get_global2d_i8_2d
     module procedure mpp_get_global2d_c8_3d
     module procedure mpp_get_global2d_c8_2d
  end interface
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface
  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface
  public :: mpp_define_domains, mpp_domains_init, mpp_domains_exit, mpp_get_global, mpp_get_halo_size, mpp_set_halo_size, &
            mpp_update_domains, operator(.EQ.), operator(.NE.)

  contains

    subroutine mpp_domains_init(flags,halosize)
!initialize domain decomp package
      integer, intent(in), optional :: flags, halosize

      if( mpp_domains_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      pe = mpp_pe()
      npes = mpp_npes()
      mpp_domains_initialized = .TRUE.
      if( pe.EQ.0 )write( stdout,'(/a)' )'MPP_DOMAINS module '//trim(version)

      if( PRESENT(flags) )verbose = flags.EQ.MPP_VERBOSE

      if( PRESENT(halosize) )then
          call mpp_set_halo_size(halosize)
      else
          call mpp_set_halo_size(32768) !default, pretty arbitrary
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

    subroutine mpp_set_halo_size(max_halo_size_in)
      integer, intent(in) :: max_halo_size_in

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_SET_HALO_SIZE: You must first call mpp_domains_init.' )
      if( max_halo_size_in.LE.max_halo_size )return
      max_halo_size = max_halo_size_in
#ifdef use_shmalloc
      call mpp_malloc( ptr_heap, 4*max_halo_size, len_heap ) !must hold 2 complex arrays of size max_halo_size
#else
      if( ALLOCATED(put_r8) )then
          deallocate(put_r8)
          deallocate(get_r8)
          deallocate(put_c8)
          deallocate(get_c8)
      end if
      allocate( put_r8(max_halo_size) )
      allocate( get_r8(max_halo_size) )
      allocate( put_c8(max_halo_size) )
      allocate( get_c8(max_halo_size) )
#endif
      if( pe.EQ.0 )write( stdout,* )'MPP_DOMAINS: size of buffer arrays for halo regions=', max_halo_size
      return
    end subroutine mpp_set_halo_size

    subroutine mpp_get_halo_size(max_halo_size_out)
      integer, intent(out) :: max_halo_size_out

      max_halo_size_out = max_halo_size
      return
    end subroutine mpp_get_halo_size

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!              MPP_GET_GLOBAL: get global field from domain field              !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_get_global2D_r8_3d( domain, local_field, global_field )
!given a 3D array defined on the data domain <domain>, constructs a global field
!USE WITH CARE! a global 3D array could occupy a lot of memory!
      type(domain2D), intent(in) :: domain
      real(DOUBLE_KIND), intent(in)  ::  local_field(domain%x%data%start_index:,domain%y%data%start_index:,:)
      real(DOUBLE_KIND), intent(out) :: global_field(domain%x%global%start_index:,domain%y%global%start_index:,:)
      type(domain2D), allocatable :: global_domain(:)

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_GET_GLOBAL2D_R8_3D: You must first call mpp_domains_init.' )
      if( size(local_field,1).NE.domain%x%data%size .OR. &
          size(local_field,2).NE.domain%y%data%size .OR. &
          size(global_field,1).NE.domain%x%global%size .OR. &
          size(global_field,2).NE.domain%y%global%size .OR. &
          size(global_field,3).NE.size(local_field,3) ) &
           call mpp_error( FATAL, 'MPP_GET_GLOBAL2D_R8_3D: argument mismatch, check domain and array sizes.' )
!copy compute domain from local to global array
      global_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index,:) = &
       local_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index,:)
      if( npes.EQ.1 )return

!update global domain
      call mpp_set_halo_size(domain%x%global%max_size*domain%y%compute%max_size*size(local_field,3))
      if( domain%x%data%is_global .AND. domain%y%data%is_global )then
          call mpp_update_domains( global_field, domain )
      else
!construct global domains
          allocate( global_domain(0:domain%x%ndomains*domain%y%ndomains-1) )
          call mpp_define_domains( (/domain%x%global%start_index,domain%x%global%end_index,   &
                                     domain%y%global%start_index,domain%y%global%end_index/), &
                                     global_domain, xflags=GLOBAL_DATA_DOMAIN, yflags=GLOBAL_DATA_DOMAIN, &
                                     domain_layout=(/domain%x%ndomains,domain%y%ndomains/),     &
                                     pe_layout=(/size(domain%x%pelist),size(domain%y%pelist)/), &
                                     xextent=domain%x%sizelist, yextent=domain%y%sizelist )
          call mpp_update_domains( global_field, global_domain(pe) )
      end if

      return
    end subroutine mpp_get_global2D_r8_3d

    subroutine mpp_get_global2D_c8_3d( domain, local_field, global_field )
!given a 3D array defined on the data domain <domain>, constructs a global field
!USE WITH CARE! a global 3D array could occupy a lot of memory!
      type(domain2D), intent(in) :: domain
      complex(DOUBLE_KIND), intent(in)  ::  local_field(domain%x%data%start_index:,domain%y%data%start_index:,:)
      complex(DOUBLE_KIND), intent(out) :: global_field(domain%x%global%start_index:,domain%y%global%start_index:,:)
      type(domain2D), allocatable :: global_domain(:)

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_GET_GLOBAL2D_C8_3D: You must first call mpp_domains_init.' )
      if( size(local_field,1).NE.domain%x%data%size .OR. &
          size(local_field,2).NE.domain%y%data%size .OR. &
          size(global_field,1).NE.domain%x%global%size .OR. &
          size(global_field,2).NE.domain%y%global%size .OR. &
          size(global_field,3).NE.size(local_field,3) ) &
           call mpp_error( FATAL, 'MPP_GET_GLOBAL2D_C8_3D: argument mismatch, check domain and array sizes.' )
!copy compute domain from local to global array
      global_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index,:) = &
       local_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index,:)
      if( npes.EQ.1 )return

!update global domain
      call mpp_set_halo_size(domain%x%global%max_size*domain%y%compute%max_size*size(local_field,3))
      if( domain%x%data%is_global .AND. domain%y%data%is_global )then
          call mpp_update_domains( global_field, domain )
      else
!construct global domains
          allocate( global_domain(0:domain%x%ndomains*domain%y%ndomains-1) )
          call mpp_define_domains( (/domain%x%global%start_index,domain%x%global%end_index,   &
                                     domain%y%global%start_index,domain%y%global%end_index/), &
                                     global_domain, xflags=GLOBAL_DATA_DOMAIN, yflags=GLOBAL_DATA_DOMAIN, &
                                     domain_layout=(/domain%x%ndomains,domain%y%ndomains/),     &
                                     pe_layout=(/size(domain%x%pelist),size(domain%y%pelist)/), &
                                     xextent=domain%x%sizelist, yextent=domain%y%sizelist )
          call mpp_update_domains( global_field, global_domain(pe) )
      end if

      return
    end subroutine mpp_get_global2D_c8_3d

    subroutine mpp_update_domain1D_r8_2d( field, domain )
!updates data domain of 2D field whose computational domains have been computed
      type(domain1D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(domain%data%start_index:,:)

      type(domain1D), pointer :: put_domain, get_domain
!limits of computation, remote put and get domains
      integer :: isc, iec,  isp, iep,  isg, ieg,  put_pe, get_pe

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      isc = domain%compute%start_index; iec = domain%compute%end_index
!update left ("prev") boundary
      put_domain => domain; get_domain => domain
      do while( put_domain.NE.NULL_DOMAIN1D .OR. get_domain.NE.NULL_DOMAIN1D )
         if( ASSOCIATED(put_domain%next) )then
             put_domain => put_domain%next
         else
             put_domain => NULL_DOMAIN1D
         end if
         if( ASSOCIATED(get_domain%prev) )then
             get_domain => get_domain%prev
         else
             get_domain => NULL_DOMAIN1D
         end if
         call get_halos_1D( domain, put_domain, get_domain, -1, put_pe, get_pe, isp, iep, isg, ieg )
         call buffer_and_transmit
         if( isg.EQ.domain%data%start_index )get_domain => NULL_DOMAIN1D
         if( put_domain.EQ.domain )put_domain => NULL_DOMAIN1D
      end do
      call mpp_sync()
!update right ("next") boundary
      put_domain => domain; get_domain => domain
      do while( put_domain.NE.NULL_DOMAIN1D .OR. get_domain.NE.NULL_DOMAIN1D )
         if( ASSOCIATED(put_domain%prev) )then
             put_domain => put_domain%prev
         else
             put_domain => NULL_DOMAIN1D
         end if
         if( ASSOCIATED(get_domain%next) )then
             get_domain => get_domain%next
         else
             get_domain => NULL_DOMAIN1D
         end if
         call get_halos_1D( domain, put_domain, get_domain, +1, put_pe, get_pe, isp, iep, isg, ieg )
         call buffer_and_transmit
         if( ieg.EQ.domain%data%end_index )get_domain => NULL_DOMAIN1D
         if( put_domain.EQ.domain )put_domain => NULL_DOMAIN1D
      end do
      call mpp_sync()

      return

      contains
        subroutine buffer_and_transmit
!buffer the halo region and send
!isp, iep, jsp, jep: limits of put domain
!put_pe, get_pe: pe of put and get domains
          integer :: i, k
          integer :: offset, put_len, get_len
          character(len=8) :: text
#ifdef use_shmalloc
          real, dimension(max_halo_size) :: put_r8, get_r8
          pointer( ptr$put, put_r8 )
          pointer( ptr$get, get_r8 )
          ptr$put = LOC(shared_heap(              1))
          ptr$get = LOC(shared_heap(max_halo_size+1))
#endif

          if( verbose )then
              call SYSTEM_CLOCK(tick)
              write( stdout,'(a,i18,a,i5,a,2i2,4i4)' ) &
                   'T=',tick, ' PE=',pe, ' BUFFER_AND_TRANSMIT: ', put_pe, get_pe, isp, iep, isg, ieg
          end if
          if( put_pe.NE.NULL_PE )then  !put is to be done: buffer input
              put_len = (iep-isp+1)*size(field,2)
              if( put_len.GT.max_halo_size )then
                  write( text,'(i8)' )put_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              call mpp_sync_self()  !check if put_r8 is still in use
              do k = 1,size(field,2)
                 offset = (iep-isp+1)*(k-1) - isp + 1
#ifdef _CRAYT3E
                 call SHMEM_GET( put_r8(isp+offset), field(isp,k), iep-isp+1, pe )
#else
                 do i = isp,iep
                    put_r8(i+offset) = field(i,k)
                 end do
#endif
              end do
          else
              put_len = 1
          end if
          get_len = (ieg-isg+1)*size(field,2)
          if( get_pe.EQ.NULL_PE )get_len = 1
          call mpp_transmit( put_r8, put_len, put_pe, get_r8, get_len, get_pe )
          if( get_pe.NE.NULL_PE )then  !get was done: unbuffer output
              if( get_len.GT.max_halo_size )then
                  write( text,'(i8)' )get_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              do k = 1,size(field,2)
                 offset = (ieg-isg+1)*(k-1) - isg + 1
#ifdef _CRAYT3E
                 call SHMEM_GET( field(isg,k), get_r8(isg+offset), ieg-isg+1, pe )
#else
                 do i = isg,ieg
                    field(i,k) = get_r8(i+offset)
                 end do
#endif
              end do
          end if

          return
        end subroutine buffer_and_transmit
    end subroutine mpp_update_domain1D_r8_2d

    subroutine mpp_update_domain2D_r8_3d( field, domain, flags_in )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(domain%x%data%start_index:,domain%y%data%start_index:,:)
      integer, intent(in), optional :: flags_in
      integer :: flags

      type(domain2D), pointer :: put_domain, get_domain
!limits of computation, remote put and get domains
      integer :: isc, iec, jsc, jec,  isp, iep, jsp, jep,  isg, ieg, jsg, jeg,  put_pe, get_pe

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags_in) )flags = flags_in
      if( flags.EQ.TRANSPOSE ) &
           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: if TRANSPOSE is set, an update direction must be set also.' )

      isc = domain%x%compute%start_index; iec = domain%x%compute%end_index
      jsc = domain%y%compute%start_index; jec = domain%y%compute%end_index
      call mpp_sync()
      if( BTEST(flags,0) )then  !WUPDATE: update western halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%east) )then
                 put_domain => put_domain%east
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%west) )then
                 get_domain => get_domain%west
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, -1, put_pe, get_pe, isp, iep, isg, ieg )
             call buffer_and_transmit
             if( isg.EQ.domain%x%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          isc = domain%x%data%start_index !reset compute domain left edge so that Y update will do corners
      end if
      call mpp_sync()
      if( BTEST(flags,1) )then  !EUPDATE: update eastern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%west) )then
                 put_domain => put_domain%west
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%east) )then
                 get_domain => get_domain%east
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, +1, put_pe, get_pe, isp, iep, isg, ieg )
             call buffer_and_transmit
             if( ieg.EQ.domain%x%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          iec = domain%x%data%end_index !reset compute domain right edge so that Y update will do corners
      end if
      call mpp_sync()
      if( BTEST(flags,2) )then  !SUPDATE: update southern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%north) )then
                 put_domain => put_domain%north
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%south) )then
                 get_domain => get_domain%south
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, -1, put_pe, get_pe, jsp, jep, jsg, jeg )
             call buffer_and_transmit
             if( jsg.EQ.domain%y%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
      call mpp_sync()
      if( BTEST(flags,3) )then  !NUPDATE: update northern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%south) )then
                 put_domain => put_domain%south
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%north) )then
                 get_domain => get_domain%north
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, +1, put_pe, get_pe, jsp, jep, jsg, jeg )
             call buffer_and_transmit
             if( jeg.EQ.domain%y%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
      call mpp_sync()

      return

      contains
        subroutine buffer_and_transmit
!buffer the halo region and send
!isp, iep, jsp, jep: limits of put domain
!isg, ieg, jsg, jeg: limits of get domain
!put_pe, get_pe: pe of put and get domains
          integer :: i, j, k
          integer :: offset, put_len, get_len
          character(len=8) :: text
#ifdef use_shmalloc
          real, dimension(max_halo_size) :: put_r8, get_r8
          pointer( ptr$put, put_r8 )
          pointer( ptr$get, get_r8 )
          ptr$put = LOC(shared_heap(              1))
          ptr$get = LOC(shared_heap(max_halo_size+1))
#endif

          if( verbose )then
              call SYSTEM_CLOCK(tick)
              write( stdout,'(a,i18,a,i5,a,2i2,8i4)' ) &
                   'T=',tick, ' PE=',pe, ' BUFFER_AND_TRANSMIT: ', put_pe, get_pe, isp, iep, jsp, jep, isg, ieg, jsg, jeg
          end if
          if( put_pe.NE.NULL_PE )then  !put is to be done: buffer input
              put_len = (iep-isp+1)*(jep-jsp+1)*size(field,3)
              if( put_len.GT.max_halo_size )then
                  write( text,'(i8)' )put_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              call mpp_sync_self()  !check if put_r8 is still in use
              do k = 1,size(field,3)
                 do j = jsp,jep
                    offset = (iep-isp+1)*( (jep-jsp+1)*(k-1) + j - jsp ) - isp + 1
#ifdef _CRAYT3E
                    call SHMEM_GET( put_r8(isp+offset), field(isp,j,k), iep-isp+1, pe )
#else
                    do i = isp,iep
                       put_r8(i+offset) = field(i,j,k)
                    end do
#endif
                 end do
              end do
          else
              put_len = 1
          end if
          get_len = (ieg-isg+1)*(jeg-jsg+1)*size(field,3)
          if( get_pe.EQ.NULL_PE )get_len = 1
          call mpp_transmit( put_r8, put_len, put_pe, get_r8, get_len, get_pe )
          if( get_pe.NE.NULL_PE )then  !get was done: unbuffer output
              if( get_len.GT.max_halo_size )then
                  write( text,'(i8)' )get_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              do k = 1,size(field,3)
                 do j = jsg,jeg
                    offset = (ieg-isg+1)*( (jeg-jsg+1)*(k-1) + j - jsg ) - isg + 1
#ifdef _CRAYT3E
                    call SHMEM_GET( field(isg,j,k), get_r8(isg+offset), ieg-isg+1, pe )
#else
                    do i = isg,ieg
                       field(i,j,k) = get_r8(i+offset)
                    end do
#endif
                 end do
              end do
          end if

          return
        end subroutine buffer_and_transmit
    end subroutine mpp_update_domain2D_r8_3d

    subroutine mpp_update_domain2D_c8_3d( field, domain, flags_in )
!updates data domain of 3D field whose computational domains have been computed
      type(domain2D), intent(in), target :: domain
      complex(DOUBLE_KIND), intent(inout) :: field(domain%x%data%start_index:,domain%y%data%start_index:,:)
      integer, intent(in), optional :: flags_in
      integer :: flags

      type(domain2D), pointer :: put_domain, get_domain
!limits of computation, remote put and get domains
      integer :: isc, iec, jsc, jec,  isp, iep, jsp, jep,  isg, ieg, jsg, jeg,  put_pe, get_pe

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: You must first call mpp_domains_init.' )

      flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags_in) )flags = flags_in
      if( flags.EQ.TRANSPOSE ) &
           call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: if TRANSPOSE is set, an update direction must be set also.' )

      isc = domain%x%compute%start_index; iec = domain%x%compute%end_index
      jsc = domain%y%compute%start_index; jec = domain%y%compute%end_index
      call mpp_sync()
      if( BTEST(flags,0) )then  !WUPDATE: update western halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%east) )then
                 put_domain => put_domain%east
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%west) )then
                 get_domain => get_domain%west
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, -1, put_pe, get_pe, isp, iep, isg, ieg )
             call buffer_and_transmit
             if( isg.EQ.domain%x%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          isc = domain%x%data%start_index !reset compute domain left edge so that Y update will do corners
      end if
      call mpp_sync()
      if( BTEST(flags,1) )then  !EUPDATE: update eastern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%west) )then
                 put_domain => put_domain%west
                 if( BTEST(flags,4) )then
                     jsp = put_domain%y%compute%start_index; jep = put_domain%y%compute%end_index
                 else
                     jsp = jsc; jep = jec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%east) )then
                 get_domain => get_domain%east
                 jsg = jsc; jeg = jec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%x, put_domain%x, get_domain%x, +1, put_pe, get_pe, isp, iep, isg, ieg )
             call buffer_and_transmit
             if( ieg.EQ.domain%x%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
          iec = domain%x%data%end_index !reset compute domain right edge so that Y update will do corners
      end if
      call mpp_sync()
      if( BTEST(flags,2) )then  !SUPDATE: update southern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%north) )then
                 put_domain => put_domain%north
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%south) )then
                 get_domain => get_domain%south
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, -1, put_pe, get_pe, jsp, jep, jsg, jeg )
             call buffer_and_transmit
             if( jsg.EQ.domain%y%data%start_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
      call mpp_sync()
      if( BTEST(flags,3) )then  !NUPDATE: update northern halo
          put_domain => domain; get_domain => domain
          do while( put_domain.NE.NULL_DOMAIN2D .OR. get_domain.NE.NULL_DOMAIN2D )
             if( ASSOCIATED(put_domain%south) )then
                 put_domain => put_domain%south
                 if( BTEST(flags,4) )then
                     isp = put_domain%x%compute%start_index; iep = put_domain%x%compute%end_index
                 else
                     isp = isc; iep = iec
                 end if
             else
                 put_domain => NULL_DOMAIN2D
             end if
             if( ASSOCIATED(get_domain%north) )then
                 get_domain => get_domain%north
                 isg = isc; ieg = iec
             else
                 get_domain => NULL_DOMAIN2D
             end if
             call get_halos_1D( domain%y, put_domain%y, get_domain%y, +1, put_pe, get_pe, jsp, jep, jsg, jeg )
             call buffer_and_transmit
             if( jeg.EQ.domain%y%data%end_index )get_domain => NULL_DOMAIN2D
             if( put_domain.EQ.domain )put_domain => NULL_DOMAIN2D
          end do
      end if
      call mpp_sync()

      return

      contains
        subroutine buffer_and_transmit
!buffer the halo region and send
!isp, iep, jsp, jep: limits of put domain
!isg, ieg, jsg, jeg: limits of get domain
!put_pe, get_pe: pe of put and get domains
          integer :: i, j, k
          integer :: offset, put_len, get_len
          character(len=8) :: text
#ifdef use_shmalloc
          complex, dimension(max_halo_size) :: put_c8, get_c8
          pointer( ptr$put, put_c8 )
          pointer( ptr$get, get_c8 )
          ptr$put = LOC(shared_heap(              1))
          ptr$get = LOC(shared_heap(2*max_halo_size+1))
#endif

          if( verbose )then
              call SYSTEM_CLOCK(tick)
              write( stdout,'(a,i18,a,i5,a,2i2,8i4)' ) &
                   'T=',tick, ' PE=',pe, ' BUFFER_AND_TRANSMIT: ', put_pe, get_pe, isp, iep, jsp, jep, isg, ieg, jsg, jeg
          end if
          if( put_pe.NE.NULL_PE )then  !put is to be done: buffer input
              put_len = (iep-isp+1)*(jep-jsp+1)*size(field,3)
              if( put_len.GT.max_halo_size )then
                  write( text,'(i8)' )put_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              call mpp_sync_self()  !check if put_c8 is still in use
              do k = 1,size(field,3)
                 do j = jsp,jep
                    offset = (iep-isp+1)*( (jep-jsp+1)*(k-1) + j - jsp ) - isp + 1
#ifdef _CRAYT3E
                    call SHMEM_GET( put_c8(isp+offset), field(isp,j,k), 2*(iep-isp+1), pe )
#else
                    do i = isp,iep
                       put_c8(i+offset) = field(i,j,k)
                    end do
#endif
                 end do
              end do
          else
              put_len = 1
          end if
          get_len = (ieg-isg+1)*(jeg-jsg+1)*size(field,3)
          if( get_pe.EQ.NULL_PE )get_len = 1
          call mpp_transmit( put_c8, put_len, put_pe, get_c8, get_len, get_pe )
          if( get_pe.NE.NULL_PE )then  !get was done: unbuffer output
              if( get_len.GT.max_halo_size )then
                  write( text,'(i8)' )get_len
                  call mpp_error( FATAL, 'BUFFER_AND_TRANSMIT: halosize too small: you require at least '//text//'.' )
              end if
              do k = 1,size(field,3)
                 do j = jsg,jeg
                    offset = (ieg-isg+1)*( (jeg-jsg+1)*(k-1) + j - jsg ) - isg + 1
#ifdef _CRAYT3E
                    call SHMEM_GET( field(isg,j,k), get_c8(isg+offset), 2*(ieg-isg+1), pe )
#else
                    do i = isg,ieg
                       field(i,j,k) = get_c8(i+offset)
                    end do
#endif
                 end do
              end do
          end if

          return
        end subroutine buffer_and_transmit
    end subroutine mpp_update_domain2D_c8_3d

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!       OVERLOADED FUNCTIONS CONGRUENT TO VARIOUS ROUTINES ABOVE               !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_update_domain2D_r8_2d( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
!converts to 3D and calls 3D version
      type(domain2D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(:,:)
      integer, intent(in), optional :: flags
      real(DOUBLE_KIND) :: field3D(size(field,1),size(field,2),1)

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'You must first call mpp_domains_init.' )
      field3D(:,:,1) = field(:,:)
      call mpp_update_domain2D_r8_3d( field3D, domain, flags )
      field(:,:) = field3D(:,:,1)
      return
    end subroutine mpp_update_domain2D_r8_2d

    subroutine mpp_update_domain2D_r8_4d( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
!converts to 3D and calls 3D version
!NOTE: assumes whole field is being passed!
      type(domain2D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(:,:,:,:)
      integer, intent(in), optional :: flags
      real(DOUBLE_KIND) :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
      integer :: i, j, k
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
#endif

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'You must first call mpp_domains_init.' )

#ifdef use_CRI_pointers
      ptr = LOC(field)
      call mpp_update_domain2D_r8_3d( field3D, domain, flags )
#else !use_CRI_pointers
      i = 0
      do k = 1,size(field,4)
         do j = 1,size(field,3)
            i = i + 1
            field3D(:,:,i) = field(:,:,j,k)
         end do
      end do
      call mpp_update_domain2D_r8_3d( field3D, domain, flags )
      i = 0
      do k = 1,size(field,4)
         do j = 1,size(field,3)
            i = i + 1
            field(:,:,j,k) = field3D(:,:,i)
         end do
      end do
#endif use_CRI_pointers
      return
    end subroutine mpp_update_domain2D_r8_4d

    subroutine mpp_update_domain2D_c8_2d( field, domain, flags )
!updates data domain of 2D field whose computational domains have been computed
!converts to 3D and calls 3D version
      type(domain2D), intent(in), target :: domain
      complex(DOUBLE_KIND), intent(inout) :: field(:,:)
      integer, intent(in), optional :: flags
      complex(DOUBLE_KIND) :: field3D(size(field,1),size(field,2),1)

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'You must first call mpp_domains_init.' )
      field3D(:,:,1) = field(:,:)
      call mpp_update_domain2D_c8_3d( field3D, domain, flags )
      field(:,:) = field3D(:,:,1)
      return
    end subroutine mpp_update_domain2D_c8_2d

    subroutine mpp_update_domain2D_c8_4d( field, domain, flags )
!updates data domain of 4D field whose computational domains have been computed
!converts to 3D and calls 3D version
!NOTE: assumes whole field is being passed!
      type(domain2D), intent(in), target :: domain
      complex(DOUBLE_KIND), intent(inout) :: field(:,:,:,:)
      integer, intent(in), optional :: flags
      complex(DOUBLE_KIND) :: field3D(size(field,1),size(field,2),size(field,3)*size(field,4))
      integer :: i, j, k
#ifdef use_CRI_pointers
      pointer( ptr, field3D )
#endif

      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'You must first call mpp_domains_init.' )

#ifdef use_CRI_pointers
      ptr = LOC(field)
      call mpp_update_domain2D_c8_3d( field3D, domain, flags )
#else !use_CRI_pointers
      i = 0
      do k = 1,size(field,4)
         do j = 1,size(field,3)
            i = i + 1
            field3D(:,:,i) = field(:,:,j,k)
         end do
      end do
      call mpp_update_domain2D_c8_3d( field3D, domain, flags )
      i = 0
      do k = 1,size(field,4)
         do j = 1,size(field,3)
            i = i + 1
            field(:,:,j,k) = field3D(:,:,i)
         end do
      end do
#endif use_CRI_pointers
      return
    end subroutine mpp_update_domain2D_c8_4d

    subroutine mpp_update_domain2D_i8_2d( field, domain, flags )
      type(domain2D), intent(in) :: domain
      integer(LONG_KIND), intent(inout) :: field(:,:)
      integer, intent(in), optional :: flags
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain, flags )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN2D_I8_2D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain2D_i8_2d

    subroutine mpp_update_domain2D_i8_3d( field, domain, flags )
      type(domain2D), intent(in) :: domain
      integer(LONG_KIND), intent(inout) :: field(:,:,:)
      integer, intent(in), optional :: flags
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2),size(field,3))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain, flags )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN2D_I8_3D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain2D_i8_3d

    subroutine mpp_update_domain2D_i8_4d( field, domain, flags )
      type(domain2D), intent(in) :: domain
      integer(LONG_KIND), intent(inout) :: field(:,:,:,:)
      integer, intent(in), optional :: flags
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2),size(field,3),size(field,4))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain, flags )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN2D_I8_4D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain2D_i8_4d

    subroutine mpp_update_domain1D_i8_2d( field, domain )
      type(domain1D), intent(in), target :: domain
      integer(LONG_KIND), intent(inout) :: field(:,:)
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN1D_I8_2D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain1D_i8_2d

    subroutine mpp_update_domain1D_i8_1d( field, domain )
      type(domain1D), intent(in), target :: domain
      integer(LONG_KIND), intent(inout) :: field(:)
      real(DOUBLE_KIND) :: field_r8(size(field,1),1)
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN1D_I8_1D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain1D_i8_1d

    subroutine mpp_update_domain1D_r8_1d( field, domain )
      type(domain1D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(:)
      real(DOUBLE_KIND) :: field_r8(size(field,1),1)
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN1D_R8_1D: currently requires CRI pointers.' )
#endif
      return
    end subroutine mpp_update_domain1D_r8_1d

    subroutine mpp_update_domain1D_r8_3d( field, domain )
      type(domain1D), intent(in), target :: domain
      real(DOUBLE_KIND), intent(inout) :: field(:,:,:)
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2)*size(field,3))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN1D_R8_3D: currently requires CRI pointers.' )
#endif
    end subroutine mpp_update_domain1D_r8_3d

    subroutine mpp_update_domain1D_i8_3d( field, domain )
      type(domain1D), intent(in), target :: domain
      integer(LONG_KIND), intent(inout) :: field(:,:,:)
      real(DOUBLE_KIND) :: field_r8(size(field,1),size(field,2)*size(field,3))
#ifdef use_CRI_pointers
      pointer( ptr, field_r8 )
      ptr = LOC(field)
      call mpp_update_domains( field_r8, domain )
#else
      call mpp_error( FATAL, 'MPP_UPDATE_DOMAIN1D_I8_3D: currently requires CRI pointers.' )
#endif
    end subroutine mpp_update_domain1D_i8_3d
              
    subroutine mpp_get_global2D_r8_2d( domain, local_field, global_field )
      type(domain2D), intent(in) :: domain
      real(DOUBLE_KIND), intent(in)  ::  local_field(:,:)
      real(DOUBLE_KIND), intent(out) :: global_field(:,:)
      real(DOUBLE_KIND) ::  local_3d(size( local_field,1),size( local_field,2),1)
      real(DOUBLE_KIND) :: global_3d(size(global_field,1),size(global_field,2),1)
      local_3d(:,:,1) = local_field(:,:)
      call mpp_get_global( domain,  local_3d, global_3d )
      global_field(:,:) = global_3d(:,:,1)
      return
    end subroutine mpp_get_global2D_r8_2d

    subroutine mpp_get_global2D_c8_2d( domain, local_field, global_field )
      type(domain2D), intent(in) :: domain
      complex(DOUBLE_KIND), intent(in)  ::  local_field(:,:)
      complex(DOUBLE_KIND), intent(out) :: global_field(:,:)
      complex(DOUBLE_KIND) ::  local_3d(size( local_field,1),size( local_field,2),1)
      complex(DOUBLE_KIND) :: global_3d(size(global_field,1),size(global_field,2),1)
      local_3d(:,:,1) = local_field(:,:)
      call mpp_get_global( domain,  local_3d, global_3d )
      global_field(:,:) = global_3d(:,:,1)
      return
    end subroutine mpp_get_global2D_c8_2d

    subroutine mpp_get_global2D_i8_3d( domain, local_field, global_field )
      type(domain2D), intent(in) :: domain
      integer(LONG_KIND), intent(in)  ::  local_field(:,:,:)
      integer(LONG_KIND), intent(out) :: global_field(:,:,:)
      real(DOUBLE_KIND) ::  local_r8(size( local_field,1),size( local_field,2),size( local_field,3))
      real(DOUBLE_KIND) :: global_r8(size(global_field,1),size(global_field,2),size(global_field,3))
#ifdef use_CRI_pointers
      pointer(  local_ptr,  local_r8 )
      pointer( global_ptr, global_r8 )
       local_ptr = LOC( local_field)
      global_ptr = LOC(global_field)
      call mpp_get_global( domain, local_r8, global_r8 )
#else
      call mpp_error( FATAL, 'MPP_GET_GLOBAL2D_I8_3D: currently requires the use of CRI pointers.' )
#endif
      return
    end subroutine mpp_get_global2D_i8_3d
      
    subroutine mpp_get_global2D_i8_2d( domain, local_field, global_field )
      type(domain2D), intent(in) :: domain
      integer(LONG_KIND), intent(in)  ::  local_field(:,:)
      integer(LONG_KIND), intent(out) :: global_field(:,:)
      integer(LONG_KIND) ::  local_3d(size( local_field,1),size( local_field,2),1)
      integer(LONG_KIND) :: global_3d(size(global_field,1),size(global_field,2),1)
      local_3d(:,:,1) = local_field(:,:)
      call mpp_get_global( domain, local_3d, global_3d )
      global_field(:,:) = global_3d(:,:,1)
      return
    end subroutine mpp_get_global2D_i8_2d
      
  end module mpp_domains_mod
