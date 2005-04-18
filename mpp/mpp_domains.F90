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

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_domains.html"/>

! <OVERVIEW>
!   <TT>mpp_domains_mod</TT> is a set of simple calls for domain
!   decomposition and domain updates on rectilinear grids. It requires the
!   module <LINK SRC="mpp.html">mpp_mod</LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   Scalable implementations of finite-difference codes are generally
!   based on decomposing the model domain into subdomains that are
!   distributed among processors. These domains will then be obliged to
!   exchange data at their boundaries if data dependencies are merely
!   nearest-neighbour, or may need to acquire information from the global
!   domain if there are extended data dependencies, as in the spectral
!   transform. The domain decomposition is a key operation in the
!   development of parallel codes.
!   
!   <TT>mpp_domains_mod</TT> provides a domain decomposition and domain
!   update API for <I>rectilinear</I> grids, built on top of the <LINK
!   SRC="mpp.html">mpp_mod</LINK> API for message passing. Features
!   of <TT>mpp_domains_mod</TT> include:
! 
!   Simple, minimal API, with free access to underlying API for more complicated stuff.
!
!   Design toward typical use in climate/weather CFD codes.
!  
!   <H4>Domains</H4>
! 
!   I have assumed that domain decomposition will mainly be in 2
!   horizontal dimensions, which will in general be the two
!   fastest-varying indices. There is a separate implementation of 1D
!   decomposition on the fastest-varying index, and 1D decomposition on
!   the second index, treated as a special case of 2D decomposition, is
!   also possible. We define <I>domain</I> as the grid associated with a <I>task</I>.
!   We define the <I>compute domain</I> as the set of gridpoints that are
!   computed by a task, and the <I>data domain</I> as the set of points
!   that are required by the task for the calculation. There can in
!   general be more than 1 task per PE, though often
!   the number of domains is the same as the processor count. We define
!   the <I>global domain</I> as the global computational domain of the
!   entire model (i.e, the same as the computational domain if run on a
!   single processor). 2D domains are defined using a derived type <TT>domain2D</TT>,
!   constructed as follows (see comments in code for more details):
!   
!   <PRE>
!     type, public :: domain_axis_spec
!        private
!        integer :: begin, end, size, max_size
!        logical :: is_global
!     end type domain_axis_spec
!     type, public :: domain1D
!        private
!        type(domain_axis_spec) :: compute, data, global, active
!        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
!        type(domain1D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain
!     type, public :: domain2D
!        private
!        type(domain1D) :: x
!        type(domain1D) :: y
!        type(domain2D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain2D
!     type(domain1D), public :: NULL_DOMAIN1D
!     type(domain2D), public :: NULL_DOMAIN2D
!   </PRE>

!   The <TT>domain2D</TT> type contains all the necessary information to
!   define the global, compute and data domains of each task, as well as the PE
!   associated with the task. The PEs from which remote data may be
!   acquired to update the data domain are also contained in a linked list
!   of neighbours.
! </DESCRIPTION>

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999

  use mpp_parameter_mod,      only : GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE
  use mpp_parameter_mod,      only : BGRID_SW, CGRID_NE, CGRID_SW, FOLD_WEST_EDGE
  use mpp_parameter_mod,      only : FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  use mpp_parameter_mod,      only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  use mpp_parameter_mod,      only : NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME
  use mpp_datatype_mod,       only : domain_axis_spec, domain1D, domain2D, DomainCommunicator2D
  use mpp_data_mod,           only : NULL_DOMAIN1D, NULL_DOMAIN2D

  use mpp_domains_util_mod,   only : mpp_domains_set_stack_size, mpp_get_compute_domain
  use mpp_domains_util_mod,   only : mpp_get_compute_domains, mpp_get_data_domain
  use mpp_domains_util_mod,   only : mpp_get_global_domain, mpp_get_domain_components
  use mpp_domains_util_mod,   only : mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.)
  use mpp_domains_reduce_mod, only : mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum
  use mpp_domains_misc_mod,   only : mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit
  use mpp_domains_misc_mod,   only : mpp_redistribute, mpp_update_domains, mpp_check_field
  use mpp_domains_define_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain

  implicit none
  private

  !--- public paramters imported from mpp_domains_parameter_mod
  public :: GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, BGRID_SW, CGRID_NE
  public :: CGRID_SW, FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE
  public :: WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
  public :: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, MPP_DOMAIN_TIME

  !--- public data type imported from mpp_datatype_mod
  public :: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D

  !--- public data imported from mpp_data_mod
  public :: NULL_DOMAIN1D, NULL_DOMAIN2D

  !--- public interface imported from mpp_domains_util_mod
  public :: mpp_domains_set_stack_size, mpp_get_compute_domain, mpp_get_compute_domains
  public :: mpp_get_data_domain, mpp_get_global_domain, mpp_get_domain_components
  public :: mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.) 

  !--- public interface imported from mpp_domains_reduce_mod
  public :: mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum

  !--- public interface imported from mpp_domains_misc_mod
  public :: mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit, mpp_redistribute
  public ::  mpp_update_domains, mpp_check_field

  !--- public interface imported from mpp_domains_define_mod
  public :: mpp_define_layout, mpp_define_domains, mpp_modify_domain

  !--- version information variables
  character(len=128), public :: version= &
       '$Id: mpp_domains.F90,v 12.0 2005/04/14 17:57:52 fms Exp $'
  character(len=128), public :: tagname= &
       '$Name: lima $'

end module mpp_domains_mod

#ifdef test_mpp_domains
program mpp_domains_test
  use mpp_mod,         only : FATAL, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, FOLD_NORTH_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain

  implicit none
#include <fms_platform.h>
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, halo=2, stackmax=4000000
  real, dimension(:,:,:), allocatable :: global, gcheck
  integer :: unit=7
  logical :: debug=.FALSE., opened
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
                                       ! mpes should be equal to npes     
  integer :: mpes = 0
  integer :: xhalo =0, yhalo =0
  namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, debug, mpes, check_parallel, xhalo, yhalo
  integer :: i, j, k
  integer :: layout(2)
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: id

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
  
  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)
  
  if( pe.EQ.mpp_root_pe() )print '(a,6i4)', 'npes, mpes, nx, ny, nz, halo=', npes, mpes, nx, ny, nz, halo
  
  allocate( global(1-halo:nx+halo,1-halo:ny+halo,nz) )
  allocate( gcheck(nx,ny,nz) )
  
!fill in global array: with k.iiijjj
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           global(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do
     
  if( .not. check_parallel) then
      call test_modify_domain()
      call test_halo_update( 'Simple' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic' )
      call test_halo_update( 'Folded' ) !includes vector field test
      
      call test_redistribute( 'Complete pelist' )
      call test_redistribute( 'Overlap  pelist' )
      call test_redistribute( 'Disjoint pelist' )
  else
      call test_parallel( )
  endif

!Balaji adding openMP tests
  call test_openmp()
 
  call mpp_domains_exit()
  call mpp_exit()
  
contains
  subroutine test_openmp()
#ifdef _OPENMP
    integer :: omp_get_num_thread, omp_get_max_threads, omp_get_thread_num
    real, allocatable :: a(:,:,:)
    type(domain2D) :: domain
    integer :: layout(2)
    integer :: i,j,k, jthr
    integer :: thrnum, maxthr
    integer(LONG_KIND) :: sum1, sum2
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( a(isd:ied,jsd:jed,nz) )
    maxthr = omp_get_max_threads()
    write( stdout(),'(a,4i4)' )'pe,js,je,maxthr=', pe, js, je, maxthr
!    write( stderr(),'(a,2i4)' )'pe,mldid=', pe, mld_id()
    if( mod(je-js+1,maxthr).NE.0 ) &
         call mpp_error( FATAL, 'maxthr must divide domain (TEMPORARY).' )
    jthr = (je-js+1)/maxthr
!$OMP PARALLEL PRIVATE(i,j,k,thrnum)
    thrnum = omp_get_thread_num()
    write( stdout(),'(a,4i4)' )'pe,thrnum,js,je=', &
         pe, thrnum, js+thrnum*jthr,js+(thrnum+1)*jthr-1
    write( stdout(),'(a,3i4)' )'pe,thrnum,node=', pe, thrnum, mpp_node()
!!$OMP DO
    do k = 1,nz
!when omp DO is commented out, user must compute j loop limits
!with omp DO, let OMP figure it out
       do j = js+thrnum*jthr,js+(thrnum+1)*jthr-1
!       do j = js,je
          do i = is,ie
             a(i,j,k) = global(i,j,k)
          end do
       end do
    end do
!!$OMP END DO
!$OMP END PARALLEL
    sum1 = mpp_chksum( a(is:ie,js:je,:) )
    sum2 = mpp_chksum( global(is:ie,js:je,:) )
    if( sum1.EQ.sum2 )then
        call mpp_error( NOTE, 'OMP parallel test OK.' )
    else
        if( mpp_pe().EQ.mpp_root_pe() )write( stderr(),'(a,2z18)' )'OMP checksums: ', sum1, sum2
        call mpp_error( FATAL, 'OMP parallel test failed.' )
    end if
#endif
    return
  end subroutine test_openmp

  subroutine test_redistribute( type )
!test redistribute between two domains
    character(len=*), intent(in) :: type
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real, allocatable, dimension(:,:,:), save :: x, y
    real, allocatable, dimension(:,:,:), save :: x2, y2
    real, allocatable, dimension(:,:,:), save :: x3, y3
    real, allocatable, dimension(:,:,:), save :: x4, y4
    real, allocatable, dimension(:,:,:), save :: x5, y5
    real, allocatable, dimension(:,:,:), save :: x6, y6
    integer, allocatable :: pelist(:)
    integer :: pemax
    
    pemax = npes/2              !the partial pelist will run from 0...pemax
    domainx%list=>NULL() !otherwise it retains memory between calls
    domainy%list=>NULL()
    
!select pelists
    select case(type)
    case( 'Complete pelist' )
!both pelists run from 0...npes-1
        allocate( pelist(0:npes-1) )
        pelist = (/ (i,i=0,npes-1) /)
        call mpp_declare_pelist( pelist )
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
        if( pemax+1.GE.npes )return
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
        call mpp_declare_pelist( (/ (i,i=pemax+1,npes-1) /))
    case default
        call mpp_error( FATAL, 'TEST_REDISTRIBUTE: no such test: '//type )
    end select
        
!set up x and y arrays
    select case(type)
    case( 'Complete pelist' )
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy, name=type )
        call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
        allocate( y(isd:ied,jsd:jed,nz) )
        allocate( y2(isd:ied,jsd:jed,nz) )
        allocate( y3(isd:ied,jsd:jed,nz) )
        allocate( y4(isd:ied,jsd:jed,nz) )
        allocate( y5(isd:ied,jsd:jed,nz) )
        allocate( y6(isd:ied,jsd:jed,nz) )
        y = 0.
        y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        end if
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
    
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        else
!set up x array
            call mpp_set_current_pelist( (/ (i,i=pemax+1,npes-1) /) )
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
            call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
            allocate( x(isd:ied,jsd:jed,nz) )
            allocate( x2(isd:ied,jsd:jed,nz) )
            allocate( x3(isd:ied,jsd:jed,nz) )
            allocate( x4(isd:ied,jsd:jed,nz) )
            allocate( x5(isd:ied,jsd:jed,nz) )
            allocate( x6(isd:ied,jsd:jed,nz) )
            x = 0.
            x(is:ie,js:je,:) = global(is:ie,js:je,:)
            x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
         end if
    end select
         
!go global and redistribute
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)
    
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_redistribute( domainx, x, domainy, y )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
        
    call mpp_set_current_pelist()

    call mpp_clock_begin(id)
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()

    if(type == 'Complete pelist')then
      write(stdout(),*) 'Use domain communicator handle'
      call mpp_clock_begin(id)
      if(ALLOCATED(y))then
         y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
      endif
      call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
      call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
      call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
      call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
      call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
      call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
      call mpp_clock_end  (id)
    
!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
    endif
    dch =>NULL()
        
    call mpp_set_current_pelist()

    call mpp_sync()
    
    if(ALLOCATED(x))then
      call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
      deallocate(x,x2,x3,x4,x5,x6)
    endif
    if(ALLOCATED(y))deallocate(y,y2,y3,y4,y5,y6)
  end subroutine test_redistribute
    
  subroutine test_halo_update( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y2, y3, y4
    type(domain2D) :: domain
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real :: lsum, gsum
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type )
        global(1-halo:0,    :,:) = 0
        global(nx+1:nx+halo,:,:) = 0
        global(:,    1-halo:0,:) = 0
        global(:,ny+1:ny+halo,:) = 0
    case( 'Cyclic' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name=type )
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo,    1-halo:0,:) = global(1-halo:nx+halo,ny-halo+1:ny,:)
        global(1-halo:nx+halo,ny+1:ny+halo,:) = global(1-halo:nx+halo,1:halo,      :)
    case( 'Folded' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type )
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo,ny+1:ny+halo,:) = global(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
        global(1-halo:nx+halo,1-halo:0,:) = 0
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select
        
!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.; x2 = 0.; x3 = 0.; x4 = 0.
    x(is:ie,js:je,:) = global(is:ie,js:je,:)
    x2(is:ie,js:je,:) = global(is:ie,js:je,:)
    x3(is:ie,js:je,:) = global(is:ie,js:je,:)
    x4(is:ie,js:je,:) = global(is:ie,js:je,:)
    
!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain, NUPDATE+EUPDATE,  complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial' )
    call compare_checksums( x2(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial' )
    call compare_checksums( x3(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial' )
    call compare_checksums( x4(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial' )
    
!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), type )
    
    select case(type)           !extra tests
    case( 'Simple' )
    
!test mpp_global_field
        id = mpp_clock_id( type//' global field', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        call mpp_global_field( domain, x, gcheck )
        call mpp_clock_end  (id)
!compare checksums between global and x arrays
        call compare_checksums( global(1:nx,1:ny,:), gcheck, 'mpp_global_field' )
        
        gcheck=0.d0
        call mpp_clock_begin(id)
        call mpp_global_field( domain, x, gcheck, new=.true., dc_handle=dch )
        call mpp_clock_end  (id)                 
!compare checksums between global and x arrays
        call compare_checksums( global(1:nx,1:ny,:), gcheck, 'mpp_global_field_new' )

        gcheck=0.d0
        call mpp_clock_begin(id)
        call mpp_global_field( domain, x, gcheck, new=.true., dc_handle=dch )
        dch =>NULL()
        call mpp_clock_end  (id)                                          
!compare checksums between global and x arrays  
        call compare_checksums( global(1:nx,1:ny,:), gcheck, 'mpp glob fld new w/ dom comm handle' )

        
!test mpp_global_sum
        gsum = sum( global(1:nx,1:ny,:) )
        id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        lsum = mpp_global_sum( domain, x )
        call mpp_clock_end  (id)
        if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Fast sum=', lsum, gsum
!test exact mpp_global_sum
        id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM )
        call mpp_clock_end  (id)
        if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Bitwise-exact sum=', lsum, gsum
        call mpp_clock_begin(id)
        lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM )
        call mpp_clock_end  (id)          
        if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'New Bitwise-exact sum=', lsum, gsum
    case( 'Folded' )!test vector update: folded, with sign flip at fold
!fill in folded north edge, cyclic east and west edge
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo-1,ny+1:ny+halo,:) = -global(nx+halo-1:1-halo:-1,ny-1:ny-halo:-1,:)
        global(nx+halo,ny+1:ny+halo,:) = -global(nx-halo,ny-1:ny-halo:-1,:)
        global(1-halo:nx+halo,1-halo:0,:) = 0
        
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
!set up y array
        allocate( y(isd:ied,jsd:jed,nz) )
        allocate( y2(isd:ied,jsd:jed,nz) )
        allocate( y3(isd:ied,jsd:jed,nz) )
        allocate( y4(isd:ied,jsd:jed,nz) )
        y = x; x2 = x; x3 = x; x4 = x
        y = x; y2 = x; y3 = x; y4 = x
        id = mpp_clock_id( type//' vector', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        call mpp_update_domains( x, y, domain, gridtype=BGRID_NE, complete=.false. )
        call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
        call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
        call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true., dc_handle=dch )
        call mpp_clock_end  (id)
!redundant points must be equal and opposite
        global(nx/2,ny,:) = 0.  !pole points must have 0 velocity
        global(nx  ,ny,:) = 0.  !pole points must have 0 velocity
        global(nx/2+1:nx-1, ny,:) = -global(nx/2-1:1:-1, ny,:)
        global(1-halo:0,    ny,:) = -global(nx-halo+1:nx,ny,:)
        global(nx+1:nx+halo,ny,:) = -global(1:halo,      ny,:)
        call compare_checksums( x, global(isd:ied,jsd:jed,:), type//' X' )
        call compare_checksums( x2, global(isd:ied,jsd:jed,:), type//' X2' )
        call compare_checksums( x3, global(isd:ied,jsd:jed,:), type//' X3' )
        call compare_checksums( x4, global(isd:ied,jsd:jed,:), type//' X4' )
        call compare_checksums( y, global(isd:ied,jsd:jed,:), type//' Y' )
        call compare_checksums( y2, global(isd:ied,jsd:jed,:), type//' Y2' )
        call compare_checksums( y3, global(isd:ied,jsd:jed,:), type//' Y3' )
        call compare_checksums( y4, global(isd:ied,jsd:jed,:), type//' Y4' )
    end select
    
    if(ALLOCATED(x))then
      call mpp_update_domains( x, domain, flags=NUPDATE+EUPDATE, free=.true., list_size=4 )
      deallocate(x,x2,x3,x4)
    endif
    if(ALLOCATED(y))then
      deallocate(y,y2,y3,y4)
    endif
  end subroutine test_halo_update
    
  subroutine test_parallel ( )
  
    integer :: npes, layout(2), i, j, k,is, ie, js, je, isd, ied, jsd, jed
    real, dimension(:,:), allocatable :: field, lfield
    real, dimension(:,:,:), allocatable :: field3d, lfield3d
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: pelist1 , pelist2
    logical :: group1, group2
    character(len=128)  :: mesg
    
    npes = mpp_npes()
    allocate(pelist1(npes-mpes), pelist2(mpes))
    pelist1 = (/(i, i = 0, npes-mpes -1)/)
    pelist2 = (/(i, i = npes-mpes, npes - 1)/)
    call mpp_declare_pelist(pelist1)
    call mpp_declare_pelist(pelist2)
    group1 = .FALSE. ; group2 = .FALSE.
    if(any(pelist1==pe)) group1 = .TRUE.
    if(any(pelist2==pe)) group2 = .TRUE.
    mesg = 'parallel checking'
    
  if(group1) then
     call mpp_set_current_pelist(pelist1)
     call mpp_define_layout( (/1,nx,1,ny/), npes-mpes, layout )
  else if(group2) then
     call mpp_set_current_pelist(pelist2)
     call mpp_define_layout( (/1,nx,1,ny/), mpes, layout )
  endif
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain,&
                                   xhalo=xhalo, yhalo=yhalo)
                                   
     call mpp_set_current_pelist() 
     
     call mpp_get_compute_domain(domain, is, ie, js, je)
     call mpp_get_data_domain(domain, isd, ied, jsd, jed)
     allocate(lfield(is:ie,js:je),field(isd:ied,jsd:jed))
     allocate(lfield3d(is:ie,js:je,nz),field3d(isd:ied,jsd:jed,nz))
     
     do i = is, ie
     do j = js, je
        lfield(i,j) = real(i)+real(j)*0.001
     enddo
     enddo
     do i = is, ie
     do j = js, je
     do k = 1, nz
        lfield3d(i,j,k) = real(i)+real(j)*0.001+real(k)*0.00001
     enddo
     enddo
     enddo
     field = 0.0
     field3d = 0.0
     field(is:ie,js:je)= lfield(is:ie,js:je)
     field3d(is:ie,js:je,:) = lfield3d(is:ie,js:je,:)
     call mpp_update_domains(field,domain)
     call mpp_update_domains(field3d,domain)
     
    call mpp_check_field(field, pelist1, pelist2,domain, '2D '//mesg, w_halo = xhalo, &
                            s_halo = yhalo, e_halo = xhalo, n_halo = yhalo)
    call mpp_check_field(field3d, pelist1, pelist2,domain, '3D '//mesg, w_halo = xhalo, &
                            s_halo = yhalo, e_halo = xhalo, n_halo = yhalo)
                            
  end subroutine test_parallel
  
  subroutine test_modify_domain( )
  
    type(domain2D) :: domain2d_no_halo, domain2d_with_halo
    type(domain1D) :: domain1d_no_halo, domain1d_with_halo
    integer :: is,ie,js,je,isd,ied,jsd,jed
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain2d_no_halo,   &
                            yflags=CYCLIC_GLOBAL_DOMAIN, xhalo=0, yhalo=0)
    call mpp_get_compute_domain(domain2d_no_halo, is, ie, js, je)
    call mpp_get_data_domain(domain2d_no_halo, isd, ied, jsd, jed)
    print*, "at pe ", mpp_pe(), "the compute domain decomposition of domain without halo", is, ie, js, je
    print*, "at pe ", mpp_pe(), "the data domain decomposition of domain without halo", isd, ied, jsd, jed
    call mpp_modify_domain(domain2d_no_halo, domain2d_with_halo, xhalo=xhalo,yhalo=yhalo)
    call mpp_get_compute_domain(domain2d_with_halo, is, ie, js, je)
    call mpp_get_data_domain(domain2d_with_halo, isd, ied, jsd, jed)
    print*, "at pe ", mpp_pe(), "the compute domain decomposition of domain with halo", is, ie, js, je
    print*, "at pe ", mpp_pe(), "the data domain decomposition of domain with halo", isd, ied, jsd, jed
    
    return
    
end subroutine test_modify_domain

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: i, j
    
    call mpp_sync()
    i = mpp_chksum( a, (/pe/) )
    j = mpp_chksum( b, (/pe/) )
    if( i.EQ.j )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums
end program mpp_domains_test
#endif

! <INFO>

!   <COMPILER NAME="">     
!     Any module or program unit using <TT>mpp_domains_mod</TT>
!     must contain the line

!     <PRE>
!     use mpp_domains_mod
!     </PRE>

!     <TT>mpp_domains_mod</TT> <TT>use</TT>s <LINK
!     SRC="mpp.html">mpp_mod</LINK>, and therefore is subject to the <LINK
!     SRC="mpp.html#COMPILING AND LINKING SOURCE">compiling and linking requirements of that module.</LINK>
!   </COMPILER>
!   <PRECOMP FLAG="">      
!     <TT>mpp_domains_mod</TT> uses standard f90, and has no special
!     requirements. There are some OS-dependent
!     pre-processor directives that you might need to modify on
!     non-SGI/Cray systems and compilers. The <LINK
!     SRC="mpp.html#PORTABILITY">portability of mpp_mod</LINK>
!     obviously is a constraint, since this module is built on top of
!     it. Contact me, Balaji, SGI/GFDL, with questions.
!   </PRECOMP> 
!   <LOADER FLAG="">       
!     The <TT>mpp_domains</TT> source consists of the main source file
!     <TT>mpp_domains.F90</TT> and also requires the following include files:
!    <PRE>
!     <TT>fms_platform.h</TT>
!     <TT>mpp_update_domains2D.h</TT>
!     <TT>mpp_global_reduce.h</TT>
!     <TT>mpp_global_sum.h</TT>
!     <TT>mpp_global_field.h</TT>
!    </PRE>
!    GFDL users can check it out of the main CVS repository as part of
!    the <TT>mpp</TT> CVS module. The current public tag is <TT>galway</TT>.
!    External users can download the latest <TT>mpp</TT> package <LINK SRC=
!    "ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.tar.Z">here</LINK>. Public access
!    to the GFDL CVS repository will soon be made available.

!   </LOADER>

! </INFO>
