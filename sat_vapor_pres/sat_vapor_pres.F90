
module sat_vapor_pres_mod

!-----------------------------------------------------------------------
!
!                 saturation vapor pressure lookup
!
!      routines for computing the saturation vapor pressure (es),
!      the derivation of es with respect to temperature, and
!      initialization of the look-up table.
!
!-----------------------------------------------------------------------
!
!                               usage
!                               -----
!
!                       call lookup_es  (temp,es)
!
!                       call lookup_des (temp,des)
!
!    arguments
!    ---------
!      temp    intent in       temperature in degrees kelvin
!      es      intent out      saturation vapor pressure in Pascals
!      des     intent out      derivative of saturation vapor pressure
!                              with respect to temperature (Pascals/degree)
!
!-----------------------------------------------------------------------

 use  constants_mod, only:  TFREEZE
 use        fms_mod, only:  write_version_number,   &
                            error_mesg, FATAL

implicit none
private

 public :: lookup_es, lookup_des, sat_vapor_pres_init
 public :: compute_es
 public :: escomp, descomp ! for backward compatibility
                           ! use lookup_es, lookup_des instead

!-----------------------------------------------------------------------
 interface lookup_es
   module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
 end interface
! for backward compatibility (to be removed soon)
 interface escomp
   module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
 end interface
!-----------------------------------------------------------------------
 interface lookup_des
   module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
 end interface
! for backward compatibility (to be removed soon)
 interface descomp
   module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
 end interface
!-----------------------------------------------------------------------
 interface compute_es
   module procedure compute_es_0d, compute_es_1d, compute_es_2d, compute_es_3d
 end interface
!-----------------------------------------------------------------------
 interface temp_check
   module procedure temp_check_0d, temp_check_1d, temp_check_2d, temp_check_3d
 end interface
!-----------------------------------------------------------------------
!  cvs version and tag name

character(len=128) :: version = '$Id: sat_vapor_pres.F90,v 1.6 2002/07/16 22:56:58 fms Exp $'
character(len=128) :: tagname = '$Name: havana $'

!-----------------------------------------------------------------------
!  parameters for table size and resolution

 integer, parameter :: tcmin = -160  ! minimum temperature (degC) in lookup table
 integer, parameter :: tcmax =  100  ! maximum temperature (degC) in lookup table
 integer, parameter :: esres =  10   ! table resolution (increments per degree)
 integer, parameter :: nsize = (tcmax-tcmin)*esres+1    !  lookup table size
 integer, parameter :: nlim  = nsize-1

 real    :: tmin, tmax          !  lookup table limits in degK
 real    :: dtres, dtinv, teps

 real ::   TABLE(nsize)    !  sat vapor pres (es)
 real ::  DTABLE(nsize)    !  first derivative of es
 real :: D2TABLE(nsize)    ! second derivative of es

 logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine lookup_es_0d ( temp, esat )

 real, intent(in)  :: temp
 real, intent(out) :: esat

 real    :: tmp, del
 integer :: ind
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   tmp = temp-tmin
   ind = int(dtinv*(tmp+teps))
   del = tmp-dtres*real(ind)
   esat = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!esat = TABLE(ind+1) + del*DTABLE(ind+1)

     if (ind < 0 .or. ind > nlim) call temp_check ( 1, temp )

!-----------------------------------------------

 end subroutine lookup_es_0d

!#######################################################################

 subroutine lookup_es_1d ( temp, esat )

 real, intent(in)  :: temp(:)
 real, intent(out) :: esat(:)

 real    :: tmp, del
 integer :: ind, i, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_1d

!#######################################################################

 subroutine lookup_es_2d ( temp, esat )

 real, intent(in)  :: temp(:,:)
 real, intent(out) :: esat(:,:)

 real    :: tmp, del
 integer :: ind, i, j, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i,j) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i,j) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo

   if ( n > 0 ) call temp_check_2d ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_2d

!#######################################################################

 subroutine lookup_es_3d ( temp, esat )

 real, intent(in)  :: temp(:,:,:)
 real, intent(out) :: esat(:,:,:)

 real    :: tmp, del
 integer :: ind, i, j, k, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i,j,k) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i,j,k) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_3d

!#######################################################################
!  routines for computing derivative of es
!#######################################################################

 subroutine lookup_des_0d ( temp, desat )

 real, intent(in)  :: temp
 real, intent(out) :: desat

 real    :: tmp, del
 integer :: ind
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   tmp = temp-tmin
   ind = int(dtinv*(tmp+teps))
   del = tmp-dtres*real(ind)
   desat = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)

   if (ind < 0 .or. ind > nlim) call temp_check ( 1, temp )

!-----------------------------------------------

 end subroutine lookup_des_0d

!#######################################################################

 subroutine lookup_des_1d ( temp, desat )

 real, intent(in)  :: temp (:)
 real, intent(out) :: desat(:)

 real    :: tmp, del
 integer :: ind, i, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_des_1d

!#######################################################################

 subroutine lookup_des_2d ( temp, desat )

 real, intent(in)  :: temp (:,:)
 real, intent(out) :: desat(:,:)

 real    :: tmp, del
 integer :: ind, i, j, n
!-----------------------------------------------
   
   if (.not.module_is_initialized) call sat_vapor_pres_init
   
   n = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i,j) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_des_2d

!#######################################################################

 subroutine lookup_des_3d ( temp, desat )

 real, intent(in)  :: temp (:,:,:)
 real, intent(out) :: desat(:,:,:)

 real    :: tmp, del
 integer :: ind, i, j, k, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i,j,k) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )
   
!-----------------------------------------------

 end subroutine lookup_des_3d

!#######################################################################
!#######################################################################

 subroutine sat_vapor_pres_init

!  =================================================================
!  +                                                               +
!  +             construction of the es table                      +
!  +                                                               +
!  + this table is constructed from es equations from the          +
!  + smithsonian tables.  the es input is computed from values     +
!  + (in one-tenth of a degree increments) of es over ice          +
!  + from -153c to 0c and values of es over water from 0c to 102c. +
!  + output table contains these data interleaved with their       +
!  + derivatives with respect to temperature except between -20c   +
!  + and 0c where blended (over water and over ice) es values and  +
!  + derivatives are calculated.                                   +
!  +   note: all es computation is done in pascals                 +
!  =================================================================

 real    :: tem(3), es(3), hdtinv
 integer :: i, n

! increment used to generate derivative table
  real, parameter :: tinrc = .01           
  real, parameter :: tfact = 1./(2.*tinrc)

! return silently if this routine has already been called
      if (module_is_initialized) return

! write version number to log file
      call write_version_number (version, tagname)

! global variables
      tmin = real(tcmin)+TFREEZE   ! minimum valid temp in table
      tmax = real(tcmax)+TFREEZE   ! maximum valid temp in table
      dtinv = real(esres)
      dtres = 1./dtinv
      teps = 1./real(2*esres)
! local variables
      hdtinv = dtinv*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      do i = 1, nsize
         tem(1) = tmin + dtres*real(i-1)
         tem(2) = tem(1)-tinrc
         tem(3) = tem(1)+tinrc
         es = compute_es (tem)
          TABLE(i) = es(1)
         DTABLE(i) = (es(3)-es(2))*tfact
      enddo

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, nsize-1
         D2TABLE(i) = 0.25*dtinv*(DTABLE(i+1)-DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries
         D2TABLE(1)     = 0.50*dtinv*(DTABLE(2)    -DTABLE(1))
         D2TABLE(nsize) = 0.50*dtinv*(DTABLE(nsize)-DTABLE(nsize-1))

    module_is_initialized = .true.

 end subroutine sat_vapor_pres_init

!#######################################################################
!#######################################################################
!-------------------------------------------------------------------
!                Computation of the es values
!
!   Saturation vapor pressure (es) values are computed from
!   equations in the Smithsonian meteorological tables page 350.
!   For temperatures < 0C, sat vapor pres is computed over ice.
!   For temperatures > -20C, sat vapor pres is computed over water.
!   Between -20C and 0C the returned value is blended (over water
!   and over ice).  All sat vapor pres values are returned in pascals.
!
!   Reference:  Smithsonian meteorological tables, page 350.
!-------------------------------------------------------------------

 function compute_es_1d (tem) result (es)
 real, intent(in) :: tem(:)
 real :: es(size(tem))

 real, parameter :: TBASW = TFREEZE+100.
 real, parameter :: TBASI = TFREEZE
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

 real    :: x, esice, esh2o
 integer :: i

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
   
 end function compute_es_1d

!--------------------------------------------------------
! overloaded functions (may noy be efficient)

 function compute_es_0d (tem) result (es)
 real, intent(in) :: tem
 real :: es
 real, dimension(1) :: tem1, es1

   tem1(1) = tem
   es1 = compute_es_1d (tem1)
   es = es1(1)

 end function compute_es_0d

!--------------------------

 function compute_es_2d (tem) result (es)
 real, intent(in) :: tem(:,:)
 real, dimension(size(tem,1),size(tem,2)) :: es
 integer :: j

    do j = 1, size(tem,2)
      es(:,j) = compute_es_1d (tem(:,j))
    enddo

 end function compute_es_2d

!--------------------------

 function compute_es_3d (tem) result (es)
 real, intent(in) :: tem(:,:,:)
 real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: es
 integer :: j, k

    do k = 1, size(tem,3)
    do j = 1, size(tem,2)
      es(:,j,k) = compute_es_1d (tem(:,j,k))
    enddo
    enddo

 end function compute_es_3d

!#######################################################################
!#######################################################################

 subroutine error_handler ( n )
 integer, intent(in) :: n
 character(len=28) :: mesg

   write (mesg,'(a21,i7)') 'table overflow, nbad=', n

   call error_mesg ('sat_vapor_pres_mod', mesg, FATAL)

!  print *, 'ERROR: ' // mesg
!  stop 111

 end subroutine error_handler

!#######################################################################

 function check_1d ( temp ) result ( nbad )
 real   , intent(in)  :: temp(:)
 integer :: nbad, ind, i
 real    :: tmp

   nbad = 0
   do i = 1, size(temp,1)
     ind = int(dtinv*(temp(i)-tmin+teps))
     if (ind < 0 .or. ind > nlim) nbad = nbad+1
   enddo

 end function check_1d

!------------------------------------------------

 function check_2d ( temp ) result ( nbad )
 real   , intent(in)  :: temp(:,:)
 integer :: nbad
 integer :: j, ind
 real    :: tmp
    nbad = 0
    do j = 1, size(temp,2)
      nbad = nbad + check_1d ( temp(:,j) )
    enddo
 end function check_2d

!#######################################################################

 subroutine temp_check_0d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp

   call error_handler (nbad)

 end subroutine temp_check_0d

!--------------------------------------------------------------

 subroutine temp_check_1d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp(:)
 integer :: i

   print *, 'Bad temperatures (dimension 1): ', (check_1d(temp(i:i)),i=1,size(temp,1))
   call error_handler (nbad)

 end subroutine temp_check_1d

!--------------------------------------------------------------

 subroutine temp_check_2d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp(:,:)
 integer :: i, j

   print *, 'Bad temperatures (dimension 1): ', (check_1d(temp(i,:)),i=1,size(temp,1))
   print *, 'Bad temperatures (dimension 2): ', (check_1d(temp(:,j)),j=1,size(temp,2))
   call error_handler (nbad)

 end subroutine temp_check_2d

!--------------------------------------------------------------

 subroutine temp_check_3d ( nbad, temp )
 integer, intent(in) :: nbad
 real, intent(in)  :: temp(:,:,:)
 integer :: i, j, k

   print *, 'Bad temperatures (dimension 1): ', (check_2d(temp(i,:,:)),i=1,size(temp,1))
   print *, 'Bad temperatures (dimension 2): ', (check_2d(temp(:,j,:)),j=1,size(temp,2))
   print *, 'Bad temperatures (dimension 3): ', (check_2d(temp(:,:,k)),k=1,size(temp,3))
   call error_handler (nbad)

 end subroutine temp_check_3d

!#######################################################################

end module sat_vapor_pres_mod

