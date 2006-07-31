#include <fms_platform.h>

 module sat_vapor_pres_k_mod

! This module is what I (pjp) think a kernel should be.
! There have been many proposals as to what a kernel should look like.
! If fact, so many different ideas have been expressed that the lack
! of agreement has greatly hampered progress.
! The only way to move forward is to limit the requirments for a kernel
! to only what is widely agreeded upon.
! I believe that there are only two things widely agreeded upon.

! 1) A kernel should be independent of the rest of FMS so that it can
!    easily be ported into another programming system.
!    This requires that a kernel does not access anything by use association.
!    The one exception is this kernel, because it is not practical for physics
!    modules to avoid using a module that computes the saturation vapor
!    pressure of water vapor.

! 2) A kernel should be stateless.
!    To facilitate this it is convenient for the kernel to define a derived
!    type which is public and which allows information associated with each
!    user to be available to every interface.
!    Typically this would be data that previously was saved as module globals
!    or was available by use association.

! 3) A kernel should not read from an external file.

! One of the things that was not widely agreeded upon is that a kernel should
! not be a fortran module. This complicates things greatly for questionable
! benefit and could be done as a second step anyway, if necessary.

 implicit none
 private

 character(len=128), parameter :: version = '$Id: sat_vapor_pres_k.F90,v 13.0.2.2 2006/05/26 01:13:36 pjp Exp $'
 character(len=128), parameter :: tagname = '$Name: memphis_2006_07 $'

 public :: sat_vapor_pres_type
 public :: sat_vapor_pres_init_k
 public :: lookup_es_k
 public :: lookup_des_k

 interface lookup_es_k
   module procedure lookup_es_k_0d
   module procedure lookup_es_k_1d
   module procedure lookup_es_k_2d
   module procedure lookup_es_k_3d
 end interface

 interface lookup_des_k
   module procedure lookup_des_k_0d
   module procedure lookup_des_k_1d
   module procedure lookup_des_k_2d
   module procedure lookup_des_k_3d
 end interface

 type sat_vapor_pres_type
   real :: dtres, teps, tmin, dtinv
   integer :: table_size
   real, dimension(:), _ALLOCATABLE :: TABLE _NULL  !  sat vapor pres (es)
   real, dimension(:), _ALLOCATABLE :: DTABLE _NULL !  first derivative of es
   real, dimension(:), _ALLOCATABLE :: D2TABLE _NULL! second derivative of es
 end type

 contains

 subroutine sat_vapor_pres_init_k(svp_data, table_size, tcmin, tcmax, TFREEZE, err_msg)

! This routine has been generalized to return tables for any temperature range and resolution

 type(sat_vapor_pres_type), intent(inout) :: svp_data
 integer, intent(in) :: table_size
 real, intent(in) :: tcmin ! TABLE(1)          = sat vapor pressure at temperature tcmin (deg C)
 real, intent(in) :: tcmax ! TABLE(table_size) = sat vapor pressure at temperature tcmax (deg C)
 real, intent(in) :: TFREEZE
 character(len=*), intent(out) :: err_msg

! increment used to generate derivative table
  real, dimension(3) :: tem(3), es(3)
  real :: hdtinv, tinrc, tfact
  integer :: i

      err_msg = ''

      if(_ALLOCATED(svp_data%TABLE) .or. _ALLOCATED(svp_data%DTABLE) .or. _ALLOCATED(svp_data%D2TABLE)) then
        err_msg = 'Attempt to allocate sat vapor pressure tables when already allocated'
        return
      else
        allocate(svp_data%TABLE(table_size), svp_data%DTABLE(table_size), svp_data%D2TABLE(table_size))
      endif

      svp_data%table_size = table_size
      svp_data%dtres = (tcmax - tcmin)/(table_size-1)
      svp_data%tmin = real(tcmin)+TFREEZE  ! minimum valid temp in table
      svp_data%dtinv = 1./svp_data%dtres
      svp_data%teps = .5*svp_data%dtres
      tinrc = .1*svp_data%dtres

! To be able to compute tables for any temperature range and resolution,
! and at the same time exactly reproduce answers from memphis revision,
! it is necessary to compute ftact differently than it is in memphis.
      tfact = 5*svp_data%dtinv

      hdtinv = svp_data%dtinv*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      do i = 1, table_size
         tem(1) = svp_data%tmin + svp_data%dtres*real(i-1)
         tem(2) = tem(1)-tinrc
         tem(3) = tem(1)+tinrc
         es = compute_es_k (tem, TFREEZE)
         svp_data%TABLE(i) = es(1)
         svp_data%DTABLE(i) = (es(3)-es(2))*tfact
      enddo

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, table_size-1
         svp_data%D2TABLE(i) = 0.25*svp_data%dtinv*(svp_data%DTABLE(i+1)-svp_data%DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries

         svp_data%D2TABLE(1) = 0.50*svp_data%dtinv*(svp_data%DTABLE(2)-svp_data%DTABLE(1))

         svp_data%D2TABLE(table_size) = 0.50*svp_data%dtinv*&
              (svp_data%DTABLE(table_size)-svp_data%DTABLE(table_size-1))

 end subroutine sat_vapor_pres_init_k

!#######################################################################

 function compute_es_k(tem, TFREEZE) result (es)
 real, intent(in) :: tem(:), TFREEZE
 real :: es(size(tem,1))
         
 real    :: x, esice, esh2o, TBASW, TBASI
 integer :: i
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

   TBASW = TFREEZE+100.
   TBASI = TFREEZE

   do i = 1, size(tem)

!  compute es over ice 

     if (tem(i) < TBASI) then
         x = -9.09718*(TBASI/tem(i)-1.0) - 3.56654*log10(TBASI/tem(i)) &
             +0.876793*(1.0-tem(i)/TBASI) + log10(ESBASI)
         esice =10.**(x)
     else
         esice = 0.
     endif

!  compute es over water greater than -20 c.
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

     if (tem(i) > -20.+TBASI) then
         x = -7.90298*(TBASW/tem(i)-1) + 5.02808*log10(TBASW/tem(i)) &
             -1.3816e-07*(10**((1-tem(i)/TBASW)*11.344)-1)        &
             +8.1328e-03*(10**((TBASW/tem(i)-1)*(-3.49149))-1)    &
             +log10(ESBASW)
         esh2o = 10.**(x)
     else
         esh2o = 0.
     endif

!  derive blended es over ice and supercooled water between -20c and 0c

     if (tem(i) <= -20.+TBASI) then
         es(i) = esice
     else if (tem(i) >= TBASI) then
         es(i) = esh2o
     else
         es(i) = 0.05*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.)*esh2o)
     endif

   enddo

 end function compute_es_k

!#######################################################################

 subroutine lookup_es_k_3d(svp_data, temp, esat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-svp_data%tmin
     ind = int(svp_data%dtinv*(tmp+svp_data%teps))
     del = tmp-svp_data%dtres*real(ind)
     esat(i,j,k) = svp_data%TABLE(ind+1) + &
     del*(svp_data%DTABLE(ind+1) + del*svp_data%D2TABLE(ind+1))
     if (ind < 0 .or. ind >= svp_data%table_size) nbad = nbad+1
   enddo
   enddo
   enddo

 end subroutine lookup_es_k_3d

!#######################################################################

 subroutine lookup_des_k_3d(svp_data, temp, desat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-svp_data%tmin
     ind = int(svp_data%dtinv*(tmp+svp_data%teps))
     del = tmp-svp_data%dtres*real(ind)
     desat(i,j,k) = svp_data%DTABLE(ind+1) + 2.*del*svp_data%D2TABLE(ind+1)
     if (ind < 0 .or. ind >= svp_data%table_size) nbad = nbad+1
   enddo
   enddo
   enddo

 end subroutine lookup_des_k_3d

!#######################################################################
 subroutine lookup_des_k_2d(svp_data, temp, desat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: desat
 integer, intent(out) :: nbad
 real, dimension(size(temp,1), size(temp,2), 1) :: temp_3d, desat_3d

 temp_3d(:,:,1) = temp
 call lookup_des_k_3d(svp_data, temp_3d, desat_3d, nbad)
 desat = desat_3d(:,:,1)

 end subroutine lookup_des_k_2d
!#######################################################################
 subroutine lookup_es_k_2d(svp_data, temp, esat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat
 integer, intent(out) :: nbad
 real, dimension(size(temp,1), size(temp,2), 1) :: temp_3d, esat_3d

 temp_3d(:,:,1) = temp
 call lookup_es_k_3d(svp_data, temp_3d, esat_3d, nbad)
 esat = esat_3d(:,:,1)

 end subroutine lookup_es_k_2d
!#######################################################################
 subroutine lookup_des_k_1d(svp_data, temp, desat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: desat
 integer, intent(out) :: nbad
 real, dimension(size(temp),1,1) :: temp_3d, desat_3d

 temp_3d(:,1,1) = temp
 call lookup_des_k_3d(svp_data, temp_3d, desat_3d, nbad)
 desat = desat_3d(:,1,1)

 end subroutine lookup_des_k_1d
!#######################################################################
 subroutine lookup_es_k_1d(svp_data, temp, esat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat
 integer, intent(out) :: nbad
 real, dimension(size(temp),1,1) :: temp_3d, esat_3d

 temp_3d(:,1,1) = temp
 call lookup_es_k_3d(svp_data, temp_3d, esat_3d, nbad)
 esat = esat_3d(:,1,1)

 end subroutine lookup_es_k_1d
!#######################################################################
 subroutine lookup_des_k_0d(svp_data, temp, desat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in)     :: temp
 real, intent(out)    :: desat
 integer, intent(out) :: nbad
 real, dimension(1,1,1) :: temp_3d, desat_3d

 temp_3d(1,1,1) = temp
 call lookup_des_k_3d(svp_data, temp_3d, desat_3d, nbad)
 desat = desat_3d(1,1,1)

 end subroutine lookup_des_k_0d
!#######################################################################
 subroutine lookup_es_k_0d(svp_data, temp, esat, nbad)
 type(sat_vapor_pres_type), intent(in) :: svp_data
 real, intent(in)     :: temp
 real, intent(out)    :: esat
 integer, intent(out) :: nbad
 real, dimension(1,1,1) :: temp_3d, esat_3d

 temp_3d(1,1,1) = temp
 call lookup_es_k_3d(svp_data, temp_3d, esat_3d, nbad)
 esat = esat_3d(1,1,1)

 end subroutine lookup_es_k_0d
!#######################################################################

 end module sat_vapor_pres_k_mod
