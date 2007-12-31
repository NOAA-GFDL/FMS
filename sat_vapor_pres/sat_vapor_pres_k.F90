
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

! 2) For the sake of thread safety, module globals should be written only at initialization.
!    In this case, the module globals are the tables and a handful of scalars.

! 3) A kernel should not read from an external file.

! One of the things that was not widely agreeded upon is that a kernel should
! not be a fortran module. This complicates things greatly for questionable
! benefit and could be done as a second step anyway, if necessary.

 implicit none
 private

 character(len=128), parameter :: version = '$Id: sat_vapor_pres_k.F90,v 15.0 2007/08/14 04:15:43 fms Exp $'
 character(len=128), parameter :: tagname = '$Name: omsk_2007_12 $'

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

 real :: dtres, tepsl, tminl, dtinvl
 integer :: table_siz
 real, dimension(:), allocatable :: TABLE   !  sat vapor pres (es)
 real, dimension(:), allocatable :: DTABLE  !  first derivative of es
 real, dimension(:), allocatable :: D2TABLE ! second derivative of es

 contains

 subroutine sat_vapor_pres_init_k(table_size, tcmin, tcmax, TFREEZE, err_msg, teps, tmin, dtinv)

! This routine has been generalized to return tables for any temperature range and resolution

 integer, intent(in) :: table_size
 real, intent(in) :: tcmin ! TABLE(1)          = sat vapor pressure at temperature tcmin (deg C)
 real, intent(in) :: tcmax ! TABLE(table_size) = sat vapor pressure at temperature tcmax (deg C)
 real, intent(in) :: TFREEZE
 character(len=*), intent(out) :: err_msg
 real, intent(out), optional :: teps, tmin, dtinv

! increment used to generate derivative table
  real, dimension(3) :: tem(3), es(3)
  real :: hdtinv, tinrc, tfact
  integer :: i

      err_msg = ''

      if(allocated(TABLE) .or. allocated(DTABLE) .or. allocated(D2TABLE)) then
        err_msg = 'Attempt to allocate sat vapor pressure tables when already allocated'
        return
      else
        allocate(TABLE(table_size), DTABLE(table_size), D2TABLE(table_size))
      endif

      table_siz = table_size
      dtres = (tcmax - tcmin)/(table_size-1)
      tminl = real(tcmin)+TFREEZE  ! minimum valid temp in table
      dtinvl = 1./dtres
      tepsl = .5*dtres
      tinrc = .1*dtres
      if(present(teps )) teps =tepsl
      if(present(tmin )) tmin =tminl
      if(present(dtinv)) dtinv=dtinvl

! To be able to compute tables for any temperature range and resolution,
! and at the same time exactly reproduce answers from memphis revision,
! it is necessary to compute ftact differently than it is in memphis.
      tfact = 5*dtinvl

      hdtinv = dtinvl*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      do i = 1, table_size
         tem(1) = tminl + dtres*real(i-1)
         tem(2) = tem(1)-tinrc
         tem(3) = tem(1)+tinrc
         es = compute_es_k (tem, TFREEZE)
         TABLE(i) = es(1)
         DTABLE(i) = (es(3)-es(2))*tfact
      enddo

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, table_size-1
         D2TABLE(i) = 0.25*dtinvl*(DTABLE(i+1)-DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries

         D2TABLE(1) = 0.50*dtinvl*(DTABLE(2)-DTABLE(1))

         D2TABLE(table_size) = 0.50*dtinvl*&
              (DTABLE(table_size)-DTABLE(table_size-1))

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

 subroutine lookup_es_k_3d(temp, esat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     esat(i,j,k) = TABLE(ind+1) + &
     del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo
   enddo
   enddo

 end subroutine lookup_es_k_3d

!#######################################################################

 subroutine lookup_des_k_3d(temp, desat, nbad)
 real, intent(in),  dimension(:,:,:)  :: temp
 real, intent(out), dimension(:,:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j, k

   nbad = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     desat(i,j,k) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo
   enddo
   enddo

 end subroutine lookup_des_k_3d

!#######################################################################
 subroutine lookup_des_k_2d(temp, desat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     desat(i,j) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo
   enddo

 end subroutine lookup_des_k_2d
!#######################################################################
 subroutine lookup_es_k_2d(temp, esat, nbad)
 real, intent(in),  dimension(:,:)  :: temp
 real, intent(out), dimension(:,:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i, j

   nbad = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     esat(i,j) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo
   enddo

 end subroutine lookup_es_k_2d
!#######################################################################
 subroutine lookup_des_k_1d(temp, desat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     desat(i) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo

 end subroutine lookup_des_k_1d
!#######################################################################
 subroutine lookup_es_k_1d(temp, esat, nbad)
 real, intent(in),  dimension(:)  :: temp
 real, intent(out), dimension(:)  :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind, i

   nbad = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tminl
     ind = int(dtinvl*(tmp+tepsl))
     del = tmp-dtres*real(ind)
     esat(i) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
     if (ind < 0 .or. ind >= table_siz) nbad = nbad+1
   enddo

 end subroutine lookup_es_k_1d
!#######################################################################
 subroutine lookup_des_k_0d(temp, desat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: desat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   del = tmp-dtres*real(ind)
   desat = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
   if (ind < 0 .or. ind >= table_siz) nbad = nbad+1

 end subroutine lookup_des_k_0d
!#######################################################################
 subroutine lookup_es_k_0d(temp, esat, nbad)
 real, intent(in)     :: temp
 real, intent(out)    :: esat
 integer, intent(out) :: nbad
 real    :: tmp, del
 integer :: ind

   nbad = 0
   tmp = temp-tminl
   ind = int(dtinvl*(tmp+tepsl))
   del = tmp-dtres*real(ind)
   esat = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
   if (ind < 0 .or. ind >= table_siz) nbad = nbad+1

 end subroutine lookup_es_k_0d
!#######################################################################

 end module sat_vapor_pres_k_mod
