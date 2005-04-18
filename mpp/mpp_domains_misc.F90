module mpp_domains_misc_mod
#include <fms_platform.h>
use mpp_mod,              only : FATAL, ANY_PE, NULL_PE, MPP_DEBUG, MPP_VERBOSE
use mpp_mod,              only : mpp_pe, mpp_npes, mpp_root_pe, stdlog, stdout, stderr
use mpp_mod,              only : mpp_error, mpp_sync, mpp_sync_self, mpp_get_current_pelist
use mpp_mod,              only : mpp_clock_id, mpp_clock_end, mpp_clock_begin, mpp_malloc
use mpp_mod,              only : mpp_init, mpp_max, mpp_send, mpp_recv, mpp_broadcast, mpp_chksum
use mpp_mod,              only : mpp_version=>version, mpp_tagname=>tagname

use mpp_parameter_mod,    only : MPP_DOMAIN_TIME, XUPDATE, YUPDATE, BGRID_NE, BGRID_SW
use mpp_parameter_mod,    only : CGRID_NE, CGRID_SW, WEST, EAST, SOUTH, NORTH, SCALAR_BIT, AGRID
use mpp_parameter_mod,    only : MAX_DOMAIN_FIELDS
use mpp_datatype_mod,     only : domain1D, domain2D, DomainCommunicator2D
use mpp_datatype_mod,     only : mpp_datatype_version=>version, mpp_datatype_tagname=>tagname

use mpp_data_mod,         only : NULL_DOMAIN1D, NULL_DOMAIN2D, pe, debug=>debug_mpp_domains
use mpp_data_mod,         only : grid_offset_type, mpp_domains_stack_size, mpp_domains_stack
use mpp_data_mod,         only : PTR_DOMAINS_STACK, ptr_info, mpp_domains_stack_hwm
use mpp_data_mod,         only : debug_gsm, verbose, module_is_initialized=>mpp_domains_is_initialized
use mpp_domains_util_mod, only : compute_overlaps, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_util_mod, only : mpp_get_compute_domain, mpp_domains_set_stack_size
use mpp_domains_util_mod, only : mpp_domains_util_version=>version, mpp_domains_util_tagname=>tagname
use mpp_domains_comm_mod, only : mpp_update_init_comm, mpp_update_free_comm
use mpp_domains_comm_mod, only : mpp_redistribute_init_comm, mpp_redistribute_free_comm
#ifdef use_CAF
use mpp_domains_comm_mod, only : mpp_associate_caf_field
#endif

implicit none
private

public :: mpp_broadcast_domain, mpp_domains_init, mpp_domains_exit 
public :: mpp_redistribute, mpp_update_domains, mpp_check_field

!--- private variables 
  logical :: domain_clocks_on=.FALSE.
  integer :: send_clock=0, recv_clock=0, unpk_clock=0
  integer :: wait_clock=0, pack_clock=0, pack_loop_clock=0

!--- version information
  character(len=128), private :: version= &
       '$Id: mpp_domains_misc.F90,v 12.0 2005/04/14 17:58:08 fms Exp $'
  character(len=128), private :: tagname= &
       '$Name: lima $'

!--- public interface

! <INTERFACE NAME="mpp_update_domains">
!  <OVERVIEW>
!     Halo updates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_update_domains</TT> is used to perform a halo update of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
!    two input data fields) is only present for <TT>real</TT> types.
!    
!    For 2D domain updates, if there are halos present along both
!    <TT>x</TT> and <TT>y</TT>, we can choose to update one only, by
!    specifying <TT>flags=XUPDATE</TT> or <TT>flags=YUPDATE</TT>. In
!    addition, one-sided updates can be performed by setting <TT>flags</TT>
!    to any combination of <TT>WUPDATE</TT>, <TT>EUPDATE</TT>,
!    <TT>SUPDATE</TT> and <TT>NUPDATE</TT>, to update the west, east, north
!    and south halos respectively. Any combination of halos may be used by
!    adding the requisite flags, e.g: <TT>flags=XUPDATE+SUPDATE</TT> or
!    <TT>flags=EUPDATE+WUPDATE+SUPDATE</TT> will update the east, west and
!    south halos.
!    
!    If a call to <TT>mpp_update_domains</TT> involves at least one E-W
!    halo and one N-S halo, the corners involved will also be updated, i.e,
!    in the example above, the SE and SW corners will be updated.
!    
!    If <TT>flags</TT> is not supplied, that is
!    equivalent to <TT>flags=XUPDATE+YUPDATE</TT>.
!    
!    The vector version is passed the <TT>x</TT> and <TT>y</TT>
!    components of a vector field in tandem, and both are updated upon
!    return. They are passed together to treat parity issues on various
!    grids. For example, on a cubic sphere projection, the <TT>x</TT> and
!    <TT>y</TT> components may be interchanged when passing from an
!    equatorial cube face to a polar face. For grids with folds, vector
!    components change sign on crossing the fold.  Paired scalar quantities
!    can also be passed with the vector version if flags=SCALAR_PAIR, in which
!    case components are appropriately interchanged, but signs are not.
!    
!    Special treatment at boundaries such as folds is also required for
!    staggered grids. The following types of staggered grids are
!    recognized:
!    
!    1) <TT>AGRID</TT>: values are at grid centers.<BR/>
!    2) <TT>BGRID_NE</TT>: vector fields are at the NE vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i+&#189;,j+&#189;) with respect to the grid centers.<BR/>
!    3) <TT>BGRID_SW</TT>: vector fields are at the SW vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i-&#189;,j-&#189;) with respect to the grid centers.<BR/>
!    4) <TT>CGRID_NE</TT>: vector fields are at the N and E faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i+&#189;,j) and (i,j+&#189;) with respect to the
!    grid centers.<BR/>
!    5) <TT>CGRID_SW</TT>: vector fields are at the S and W faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i-&#189;,j) and (i,j-&#189;) with respect to the
!    grid centers.
!
!    The gridtypes listed above are all available by use association as
!    integer parameters. The scalar version of <TT>mpp_update_domains</TT>
!    assumes that the values of a scalar field are always at <TT>AGRID</TT>
!    locations, and no special boundary treatment is required. If vector
!    fields are at staggered locations, the optional argument
!    <TT>gridtype</TT> must be appropriately set for correct treatment at
!    boundaries.
!    
!    It is safe to apply vector field updates to the appropriate arrays
!    irrespective of the domain topology: if the topology requires no
!    special treatment of vector fields, specifying <TT>gridtype</TT> will
!    do no harm.
!
!    <TT>mpp_update_domains</TT> internally buffers the date being sent
!    and received into single messages for efficiency. A turnable internal
!    buffer area in memory is provided for this purpose by
!    <TT>mpp_domains_mod</TT>. The size of this buffer area can be set by
!    the user by calling <LINK SRC="mpp_domains.html#mpp_domains_set_stack_size">
!    <TT>mpp_domains_set_stack_size</TT></LINK>.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_update_domains( field, domain, flags )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_update_domains( fieldx, fieldy, domain, flags, gridtype )
!  </TEMPLATE>
! </INTERFACE>
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
#ifndef no_4byte_reals
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv
#endif
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
     module procedure mpp_update_domain2D_l4_2d
     module procedure mpp_update_domain2D_l4_3d
     module procedure mpp_update_domain2D_l4_4d
     module procedure mpp_update_domain2D_l4_5d
  end interface


  interface mpp_do_update
     module procedure mpp_do_update_new_r8_3d
     module procedure mpp_do_update_old_r8_3d
     module procedure mpp_do_update_new_r8_3dv
     module procedure mpp_do_update_old_r8_3dv
     module procedure mpp_do_update_new_c8_3d
     module procedure mpp_do_update_old_c8_3d
#ifndef no_8byte_integers
     module procedure mpp_do_update_new_i8_3d
     module procedure mpp_do_update_old_i8_3d
     module procedure mpp_do_update_new_l8_3d
     module procedure mpp_do_update_old_l8_3d
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_update_new_r4_3d
     module procedure mpp_do_update_old_r4_3d
     module procedure mpp_do_update_new_r4_3dv
     module procedure mpp_do_update_old_r4_3dv
     module procedure mpp_do_update_new_c4_3d
     module procedure mpp_do_update_old_c4_3d
#endif
     module procedure mpp_do_update_new_i4_3d
     module procedure mpp_do_update_old_i4_3d
     module procedure mpp_do_update_new_l4_3d
     module procedure mpp_do_update_old_l4_3d
  end interface

! <INTERFACE NAME="mpp_redistribute">
!  <OVERVIEW>
!    Reorganization of distributed global arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_redistribute</TT> is used to reorganize a distributed
!    array.  <TT>MPP_TYPE_</TT> can be of type <TT>integer</TT>,
!    <TT>complex</TT>, or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_redistribute( domain_in, field_in, domain_out, field_out )
!  </TEMPLATE>
!  <IN NAME="field_in" TYPE="MPP_TYPE_">
!    <TT>field_in</TT> is dimensioned on the data domain of <TT>domain_in</TT>.
!  </IN>
!  <OUT NAME="field_out" TYPE="MPP_TYPE_">
!    <TT>field_out</TT> on the data domain of <TT>domain_out</TT>.
!  </OUT>
! </INTERFACE>
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
#ifndef no_4byte_reals
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
#endif
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
     module procedure mpp_redistribute_l4_2D
     module procedure mpp_redistribute_l4_3D
     module procedure mpp_redistribute_l4_4D
     module procedure mpp_redistribute_l4_5D
  end interface

  interface mpp_do_redistribute
     module procedure mpp_do_redistribute_new_r8_3D
     module procedure mpp_do_redistribute_old_r8_3D
     module procedure mpp_do_redistribute_new_c8_3D
     module procedure mpp_do_redistribute_old_c8_3D
#ifndef no_8byte_integers
     module procedure mpp_do_redistribute_new_i8_3D
     module procedure mpp_do_redistribute_old_i8_3D
     module procedure mpp_do_redistribute_new_l8_3D
     module procedure mpp_do_redistribute_old_l8_3D
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_redistribute_new_r4_3D
     module procedure mpp_do_redistribute_old_r4_3D
     module procedure mpp_do_redistribute_new_c4_3D
     module procedure mpp_do_redistribute_old_c4_3D
#endif
     module procedure mpp_do_redistribute_new_i4_3D
     module procedure mpp_do_redistribute_old_i4_3D
     module procedure mpp_do_redistribute_new_l4_3D
     module procedure mpp_do_redistribute_old_l4_3D
  end interface


! <INTERFACE NAME="mpp_check_field">
!   <OVERVIEW>
!     Parallel checking between two ensembles which run
!     on different set pes at the same time.
!   </OVERVIEW>
!   <DESCRIPTION>
!     There are two forms for the <TT>mpp_check_field</TT> call. The 2D
!     version is generally to be used and 3D version is  built by repeated calls to the
!     2D version.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_check_field(field_in, pelist1, pelist2, domain, mesg, &
!                                w_halo, s_halo, e_halo, n_halo, force_abort  )
!   </TEMPLATE>
!   <IN NAME="field_in" >
!     Field to be checked
!   </IN>
!   <IN NAME="pelist1, pelist2">
!     Pelist of the two ensembles to be compared
!   </IN>
!   <IN NAME="domain">
!     Domain of current pe
!   </IN>
!   <IN NAME="mesg" >
!     Message to be printed out
!   </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo">
!     Halo size to be checked. Default value is 0.
!   </IN>
!   <IN NAME="force_abort">
!     When true, abort program when any difference found. Default value is false.
!   </IN>
! </INTERFACE>

  interface mpp_check_field
     module procedure mpp_check_field_2D
     module procedure mpp_check_field_3D
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                 MPP_DOMAINS: initialization and termination                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_domains_init">
!  <OVERVIEW>
!    Initialize domain decomp package.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Called to initialize the <TT>mpp_domains_mod</TT> package.
!    
!    <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to have
!    <TT>mpp_domains_mod</TT> keep you informed of what it's up
!    to. <TT>MPP_DEBUG</TT> returns even more information for debugging.
!    
!    <TT>mpp_domains_init</TT> will call <TT>mpp_init</TT>, to make sure
!    <LINK SRC="mpp.html"><TT>mpp_mod</TT></LINK> is initialized. (Repeated
!    calls to <TT>mpp_init</TT> do no harm, so don't worry if you already
!    called it).
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_init(flags)
!  </TEMPLATE>
!  <IN NAME="flags" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_domains_init(flags)
      integer, intent(in), optional :: flags
      integer                       :: l=0

      if( module_is_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      module_is_initialized = .TRUE.
      if( pe.EQ.mpp_root_pe() )then
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(version)//trim(tagname)
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(mpp_version)//trim(mpp_tagname)
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(mpp_datatype_version)//trim(mpp_datatype_tagname)
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(mpp_domains_util_version)//trim(mpp_domains_util_tagname)
!          write( stdlog(),'(a)' )trim(version_update_domains2D)
      end if

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
          domain_clocks_on = flags.EQ.MPP_DOMAIN_TIME
      end if

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary
#ifdef use_libSMA
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
          pack_clock      = mpp_clock_id( 'Halo pack' )
          pack_loop_clock = mpp_clock_id( 'Halo pack loop' )
          send_clock      = mpp_clock_id( 'Halo send' )
          recv_clock      = mpp_clock_id( 'Halo recv' )
          unpk_clock      = mpp_clock_id( 'Halo unpk' )
          wait_clock      = mpp_clock_id( 'Halo wait' )
      end if
      return
    end subroutine mpp_domains_init

!#####################################################################
! <SUBROUTINE NAME="mpp_domains_exit">
!  <OVERVIEW>
!    Exit <TT>mpp_domains_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Serves no particular purpose, but is provided should you require to
!    re-initialize <TT>mpp_domains_mod</TT>, for some odd reason.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_exit()
!  </TEMPLATE>
! </SUBROUTINE>
    subroutine mpp_domains_exit()
      if( .NOT.module_is_initialized )return
      call mpp_max(mpp_domains_stack_hwm)
      if( pe.EQ.mpp_root_pe() )write( stdout(),* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      module_is_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_CHECK_FIELD: Check parallel                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_check_field_3D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:,:)" > </IN>
!   <IN NAME="pelist1, pelist2" TYPE="integer, dimension(:)" > </IN>
!   <IN NAME="domain" TYPE="type(domain2d)" > </IN>
!   <IN NAME="mesg" TYPE="character(len=*)" > </IN>
!   <IN NAME="w_halo, s_halo, e_halo, n_halo" TYPE="integer, optional" > </IN>
!   <IN NAME="force_abort" TYPE="logical,optional" > </IN>
! </SUBROUTINE>

  subroutine mpp_check_field_3D(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort  )
!  This routine is used to do parallel checking for 3d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when true, call mpp_error if any difference
                                                    ! found. default value is false.              

    integer :: k
    character(len=256) :: temp_mesg


    do k = 1, size(field_in,3)
       write(temp_mesg, '(a, i3)') trim(mesg)//" at level " , k
       call mpp_check_field_2d(field_in(:,:,k), pelist1, pelist2, domain, temp_mesg, &
                                 w_halo, s_halo, e_halo, n_halo, force_abort )
    enddo
    
  end subroutine mpp_check_field_3D
  
  
!#####################################################################################
! <SUBROUTINE NAME="mpp_check_field_2D" INTERFACE="mpp_check_field">
!   <IN NAME="field_in" TYPE="real, dimension(:,:)" > </IN>
! </SUBROUTINE>
  subroutine mpp_check_field_2d(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort  )
!  This routine is used to do parallel checking for 2d data between n and m pe. The comparison is
!  is done on pelist2. When size of pelist2 is 1, we can check the halo; otherwise,
!  halo can not be checked.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
                                                    
     if(size(pelist2(:)) == 1) then                 
        call mpp_check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) == 1) then
        call mpp_check_field_2d_type1(field_in, pelist2, pelist1, domain, mesg, &
                                      w_halo, s_halo, e_halo, n_halo, force_abort )
     else if(size(pelist1(:)) .gt. 1 .and. size(pelist2(:)) .gt. 1) then
        call mpp_check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg, force_abort )
     else
        call mpp_error(FATAL, 'mpp_check_field: size of both pelists should be greater than 0')
     endif

  end subroutine mpp_check_field_2D


!####################################################################################

  subroutine mpp_check_field_2d_type1(field_in, pelist1, pelist2, domain, mesg, &
                                 w_halo, s_halo, e_halo, n_halo,force_abort  )
!  This routine is used to check field between running on 1 pe (pelist2) and
!  n pe(pelist1). The need_to_be_checked data is sent to the pelist2 and All the
!  comparison is done on pelist2.

  real, dimension(:,:), intent(in)  :: field_in ! field to be checked
  integer, dimension(:), intent(in) :: pelist1, pelist2 ! pe list for the two groups
  type(domain2d),        intent(in) :: domain   ! domain for each pe
  character(len=*),     intent(in)  :: mesg     ! message to be printed out
                                               ! if differences found      
  integer, intent(in), optional     :: w_halo,  s_halo, e_halo, n_halo
                                       ! halo size for west, south, east and north
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.         
! some local data

  integer                :: pe,npes, p
  integer                :: hwest, hsouth, heast, hnorth, isg, ieg, jsg, jeg, xhalo, yhalo
  integer                :: i,j,im,jm,l,is,ie,js,je,isc,iec,jsc,jec,isd,ied,jsd,jed
  real,dimension(:,:), allocatable :: field1,field2
  real,dimension(:),   allocatable :: send_buffer
  integer, dimension(4)  ::  ibounds
  logical                :: check_success,  error_exit

  check_success = .TRUE.
  error_exit    = .FALSE.
  if(present(force_abort)) error_exit = force_abort
  hwest  = 0; if(present(w_halo)) hwest  = w_halo
  heast  = 0; if(present(e_halo)) heast  = e_halo
  hsouth = 0; if(present(s_halo)) hsouth = s_halo
  hnorth = 0; if(present(n_halo)) hnorth = n_halo

  pe = mpp_pe ()
  npes = mpp_npes()

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  xhalo = isc - isd
  yhalo = jsc - jsd
  !--- need to checked halo size should not be bigger than x_halo or y_halo
  if(hwest .gt. xhalo .or. heast .gt. xhalo .or. hsouth .gt. yhalo .or. hnorth .gt. yhalo) &
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': The halo size is not correct')

  is = isc - hwest; ie = iec + heast; js = jsc - hsouth; je = jec + hnorth
  allocate(field2(is:ie,js:je))

  ! check if the field_in is on compute domain or data domain
  if((size(field_in,1) .eq. iec-isc+1) .and. (size(field_in,2) .eq. jec-jsc+1)) then
     !if field_in on compute domain, you can not check halo points
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(:,:) = field_in(:,:)
  else if((size(field_in,1) .eq. ied-isd+1) .and. (size(field_in,2) .eq. jed-jsd+1)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else if((size(field_in,1) .eq. ieg-isg+1) .and. (size(field_in,2) .eq. jeg-jsg+1)) then
     if( hwest .ne. 0 .or. heast .ne. 0 .or. hsouth .ne. 0 .or. hnorth .ne. 0 ) &
         call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//': field is on compute domain, can not check halo')
     field2(is:ie,js:je) = field_in(1:ie-is+1,1:je-js+1)
  else if((size(field_in,1) .eq. ieg-isg+1+2*xhalo) .and. (size(field_in,2) .eq. jeg-jsg+1+2*yhalo)) then
     field2(is:ie,js:je) = field_in(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1)
  else
     print*, 'on pe ', pe, 'domain: ', isc, iec, jsc, jec, isd, ied, jsd, jed, 'size of field: ', size(field_in,1), size(field_in,2)
     call mpp_error(FATAL,'mpp_check_field: '//trim(mesg)//':field is not on compute, data or global domain')
  endif
     
  call mpp_sync_self()
  
  if(any(pelist1 == pe)) then  ! send data to root pe
  
     im = ie-is+1; jm=je-js+1
     allocate(send_buffer(im*jm))
     
     ibounds(1) = is; ibounds(2) = ie; ibounds(3) = js; ibounds(4) = je
     l = 0
     do i = is,ie
     do j = js,je
        l = l+1
        send_buffer(l) = field2(i,j)
     enddo
     enddo
!  send the check bounds and data to the root pe
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_send(ibounds(1), plen=4, to_pe=pelist2(1))
     call mpp_send(send_buffer(1),plen=im*jm, to_pe=pelist2(1))
     deallocate(send_buffer)
     
   else if(pelist2(1) == pe) then        ! receive data and compare
     do p = pelist1(1), pelist1(size(pelist1(:)))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(ibounds(1), glen=4,from_pe=p)
        is = ibounds(1); ie = ibounds(2); js=ibounds(3); je=ibounds(4)
        im = ie-is+1; jm=je-js+1
        if(allocated(field1)) deallocate(field1)
        if(allocated(send_buffer)) deallocate(send_buffer)
        allocate(field1(is:ie,js:je),send_buffer(im*jm))
     ! Force use of "scalar", integer pointer mpp interface
        call mpp_recv(send_buffer(1),glen=im*jm,from_pe=p)
        l = 0
        
!  compare here, the comparison criteria can be changed according to need
        do i = is,ie
        do j = js,je
           l = l+1
           field1(i,j) = send_buffer(l)
           if(field1(i,j) .ne. field2(i,j)) then
         !   write to standard output
             print*,trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, pass_field(i,j), field_check(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
           endif
        enddo
        enddo
      enddo
        
      if(check_success) then
         print*, trim(mesg)//": ", 'comparison between 1 pe and ', npes-1, ' pes is ok'
      endif
  ! release memery
      deallocate(field1, send_buffer)
    endif
      
    deallocate(field2)
    
    call mpp_sync()
    
  end subroutine mpp_check_field_2d_type1
  
!####################################################################

  subroutine mpp_check_field_2d_type2(field_in, pelist1, pelist2, domain, mesg,force_abort)
!  This routine is used to check field between running on m pe (root pe) and
!  n pe. This routine can not check halo.
 
  real, dimension(:,:),  intent(in) :: field_in
  type(domain2d),        intent(in) :: domain
  integer, dimension(:), intent(in) :: pelist1
  integer, dimension(:), intent(in) :: pelist2
  character(len=*),      intent(in) :: mesg
  logical, intent(in), optional     :: force_abort   ! when, call mpp_error if any difference
                                                    ! found. default value is false.
! some local variables                              
  logical                :: check_success, error_exit
  real, dimension(:,:), allocatable :: field1, field2
  integer :: i, j, pe, npes, isd,ied,jsd,jed, isc, iec, jsc, jec, is, ie, js, je
  type(domain2d) :: domain1, domain2
  
    check_success = .TRUE.
    error_exit    = .FALSE.
    if(present(force_abort)) error_exit = force_abort
    pe = mpp_pe()
    npes = mpp_npes()
    call mpp_sync_self()
    if(any(pelist1 == pe)) domain1 = domain
    if(any(pelist2 == pe)) domain2 = domain
    
!  Comparison is made on pelist2.
    if(any(pelist2 == pe)) then
       call mpp_get_data_domain(domain2, isd, ied, jsd, jed)
       call mpp_get_compute_domain(domain2, is, ie, js, je)
       allocate(field1(isd:ied, jsd:jed),field2(isd:ied, jsd:jed))
      if((size(field_in,1) .ne. ied-isd+1) .or. (size(field_in,2) .ne. jed-jsd+1)) &
         call mpp_error(FATAL,'mpp_check_field: input field is not on the data domain')
      field2(isd:ied, jsd:jed) = field_in(:,:)
    endif
      
!  broadcast domain
    call mpp_broadcast_domain(domain1)
    call mpp_broadcast_domain(domain2)
    
    call mpp_redistribute(domain1,field_in,domain2,field1)
    
    if(any(pelist2 == pe)) then
        do i =is,ie
        do j =js,je
          if(field1(i,j) .ne. field2(i,j)) then
             print*, trim(mesg)//": ", i, j, field1(i,j), field2(i,j), field1(i,j) - field2(i,j)
!             write(stdout(),'(a,2i,2f)') trim(mesg), i, j, field_check(i,j), field_out(i,j)
             check_success = .FALSE.
             if(error_exit) call mpp_error(FATAL,"mpp_check_field: can not reproduce at this point")
          endif
        enddo
        enddo
        if(check_success) &
             print*, trim(mesg)//": ", 'comparison between ', size(pelist1(:)), ' pes and ', &
                   size(pelist2(:)), ' pe on', pe, ' pes is ok'
    endif          
                   
    if(any(pelist2 == pe))    deallocate(field1, field2)
    
    call mpp_sync()
    
    return
    
  end subroutine mpp_check_field_2d_type2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_BROADCAST_DOMAIN                                           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_domain( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D), intent(inout) :: domain
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(5) :: msg, info         !pe and compute domain of each item in list
      
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN: You must first call mpp_domains_init.' ) 
      
!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)
      
!am I part of this domain?
      native = ASSOCIATED(domain%list)
      
!set local list size
      if( native )then
          listsize = size(domain%list(:))
      else
          listsize = 0
      end if
      call mpp_max(listsize)
      
      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          domain%x%compute%begin =  1
          domain%x%compute%end   = -1
          domain%y%compute%begin =  1
          domain%y%compute%end   = -1
      end if
!initialize values in info
      info(1) = domain%pe
      call mpp_get_compute_domain( domain, info(2), info(3), info(4), info(5) )
      
!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( pe.EQ.pes(n) .AND. debug )write( stderr(),* )'PE ', pe, 'broadcasting msg ', msg
         call mpp_broadcast( msg, 5, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe = msg(1)
             domain%list(listpos)%x%compute%begin = msg(2)
             domain%list(listpos)%x%compute%end   = msg(3)
             domain%list(listpos)%y%compute%begin = msg(4)
             domain%list(listpos)%y%compute%end   = msg(5)
             listpos = listpos + 1
             if( debug )write( stderr(),* )'PE ', pe, 'received domain from PE ', msg(1), 'is,ie,js,je=', msg(2:5)
         end if
      end do 
    end subroutine mpp_broadcast_domain

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

#ifndef no_4byte_reals
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
#endif

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

!*******************************************************
#define VECTOR_FIELD_
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_r8_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_r8_3d
#ifdef  VECTOR_FIELD_
#define MPP_DOMAINS_DO_UPDATE_3Dnew_V_ mpp_do_update_new_r8_3dv
#define MPP_DOMAINS_DO_UPDATE_3Dold_V_ mpp_do_update_old_r8_3dv
#endif
#include <mpp_do_update_old.h>
#include <mpp_do_updateV_old.h>
#if defined(use_GSM)
#include <mpp_do_update_gsm.h>
#include <mpp_do_updateV_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#include <mpp_do_update_caf.h>
#include <mpp_do_updateV_caf.h>
#else
#include <mpp_do_update_new.h>
#include <mpp_do_updateV_new.h>
#endif

#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_c8_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_c8_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif


#ifndef no_8byte_integers
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_i8_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_i8_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif


#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_l8_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_l8_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif

#endif

#ifndef no_4byte_reals
#define VECTOR_FIELD_
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_r4_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_r4_3d
#ifdef  VECTOR_FIELD_
#define MPP_DOMAINS_DO_UPDATE_3Dnew_V_ mpp_do_update_new_r4_3dv
#define MPP_DOMAINS_DO_UPDATE_3Dold_V_ mpp_do_update_old_r4_3dv
#endif
#include <mpp_do_update_old.h>
#include <mpp_do_updateV_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#include <mpp_do_updateV_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#include <mpp_do_update_caf.h>
#include <mpp_do_updateV_caf.h>
#else
#include <mpp_do_update_new.h>
#include <mpp_do_updateV_new.h>
#endif

#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_c4_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_c4_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif

#endif

#define MPP_TYPE_ integer(INT_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_i4_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_i4_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif

#define MPP_TYPE_ logical(INT_KIND)
#define MPP_DOMAINS_DO_UPDATE_3Dnew_ mpp_do_update_new_l4_3d
#define MPP_DOMAINS_DO_UPDATE_3Dold_ mpp_do_update_old_l4_3d
#include <mpp_do_update_old.h>
#if defined(use_GSM)              
#include <mpp_do_update_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#include <mpp_do_update_caf.h>
#else
#include <mpp_do_update_new.h>
#endif


!********************************************************
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_r8_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_r8_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_c8_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_c8_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif


#ifndef no_8byte_integers

#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_i8_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_i8_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif


#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_l8_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_l8_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif

#endif

#ifndef no_4byte_reals

#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_r4_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_r4_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_c4_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_c4_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif

#endif

#define MPP_TYPE_ integer(INT_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_i4_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_i4_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif


#define MPP_TYPE_ logical(INT_KIND)
#define MPP_DO_REDISTRIBUTE_3Dnew_ mpp_do_redistribute_new_l4_3D
#define MPP_DO_REDISTRIBUTE_3Dold_ mpp_do_redistribute_old_l4_3D
#include <mpp_do_redistribute_old.h>
#if defined(use_GSM)
#include <mpp_do_redistribute_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#include <mpp_do_redistribute_caf.h>
#else
#include <mpp_do_redistribute_new.h>
#endif

end module mpp_domains_misc_mod
