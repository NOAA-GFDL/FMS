
module sat_vapor_pres_mod

!-----------------------------------------------------------------------
!
!                 SATURATION VAPOR PRESSURE LOOKUP
!
!      Routines for computing the saturation vapor pressure (es),
!      the derivation of es with respect to temperature, and
!      initialization of the look-up table.
!
!-----------------------------------------------------------------------
!
!                               USAGE
!                               -----
!
!                       call tcheck  (temp,nbad)
!
!                       call escomp  (temp,es)
!
!                       call descomp (temp,des)
!
!    ARUMENTS
!    --------
!      temp    intent IN       Temperature in degrees Kelvin
!      es      intent OUT      Saturation vapor pressure
!      des     intent OUT      Derivative of saturation vapor pressure
!                              with respect to temperature
!      nbad    intent OUT      Number of temperatures that are outside
!                              the range of the lookup table
!
!-----------------------------------------------------------------------

use  constants_mod, only:  tfreeze
use  utilities_mod, only:  get_my_pe

implicit none
private

public :: escomp, descomp, tcheck

!-----------------------------------------------------------------------

Interface escomp
  module procedure escomp_0d, escomp_1d, escomp_2d, escomp_3d, &
		   escomp_2d_flag
end Interface

!-----------------------------------------------------------------------

Interface descomp
  module procedure descomp_0d, descomp_1d, descomp_2d, descomp_3d
end Interface

!-----------------------------------------------------------------------

Interface tcheck
  module procedure tcheck_0d, tcheck_1d, tcheck_2d, tcheck_3d, &
		   tcheck_0d_reset, tcheck_1d_reset, tcheck_2d_reset, &
		   tcheck_3d_reset
end Interface

!-----------------------------------------------------------------------

!real, parameter :: Tmin=-153.16+tfreeze, Tmax=102.00+tfreeze
real, parameter :: Tmin=-153.15+tfreeze, Tmax=102.00+tfreeze

!-----------------------------------------------------------------------

!  Freezing point -- may want to use global value ???

real, parameter :: tx00=153.00-tfreeze,tx05=153.05-tfreeze

!-----------------------------------------------------------------------

real    :: ETABL(5102)
Logical :: ESSET_DONE = .FALSE.

!-----------------------------------------------------------------------

!    ------- MKS units ------
real  :: UNITS = 0.10

!    ------- CGS units ------
!real :: UNITS = 1.00

!-----------------------------------------------------------------------

CONTAINS

!#######################################################################
!#######################################################################

subroutine escomp_0d (temp,esat)

!------------------- public interface arrays ---------------------------
real, intent(IN)  :: temp
real, intent(OUT) :: esat

!------------------- private local storage -----------------------------
   real ::  work
integer :: iwork
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

       If (.NOT.ESSET_DONE) call ESSET

         iwork = 10.*(temp+tx05)
          work = temp+tx00-0.10*Float(iwork)
         iwork = 2*iwork+1
          esat = es(iwork,work)

!-----------------------------------------------------------------------

end subroutine escomp_0d

!#######################################################################

subroutine escomp_1d (temp,esat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:) :: temp
real, intent(OUT), dimension(:) :: esat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1)) ::  work
integer, dimension(size(temp,1)) :: iwork
integer  i
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.NOT.ESSET_DONE) call ESSET

      do i=1,size(temp,1)
         iwork(i) = 10.*(temp(i)+tx05)
          work(i) = temp(i)+tx00-0.10*Float(iwork(i))
         iwork(i) = 2*iwork(i)+1
          esat(i) = es(iwork(i),work(i))
      enddo

!-----------------------------------------------------------------------

end subroutine escomp_1d

!#######################################################################

subroutine escomp_2d (temp,esat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:,:) :: temp
real, intent(OUT), dimension(:,:) :: esat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1),size(temp,2)) ::  work
integer, dimension(size(temp,1),size(temp,2)) :: iwork
integer  i,j
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do j=1,size(temp,2)
      do i=1,size(temp,1)
         iwork(i,j) = 10.*(temp(i,j)+tx05)
          work(i,j) = temp(i,j)+tx00-0.10*Float(iwork(i,j))
         iwork(i,j) = 2*iwork(i,j)+1
          esat(i,j) = es(iwork(i,j),work(i,j))
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine escomp_2d

!#######################################################################

subroutine escomp_2d_flag (temp, esat, flag)

!------------------- public interface arrays ---------------------------
real, intent(IN),    dimension(:,:) :: temp
logical, intent(in), dimension(:,:) :: flag
real, intent(OUT),   dimension(:,:) :: esat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1),size(temp,2)) ::  work
integer, dimension(size(temp,1),size(temp,2)) :: iwork
integer  i,j
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do j=1,size(temp,2)
        do i=1,size(temp,1)
	  if (.not. flag(i,j)) then
            iwork(i,j) = 10.*(temp(i,j)+tx05)
            work(i,j) = temp(i,j)+tx00-0.10*Float(iwork(i,j))
            iwork(i,j) = 2*iwork(i,j)+1
            esat(i,j) = es(iwork(i,j),work(i,j))
          endif
        enddo
      enddo

!-----------------------------------------------------------------------

end subroutine escomp_2d_flag

!#######################################################################

subroutine escomp_3d (temp,esat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:,:,:) :: temp
real, intent(OUT), dimension(:,:,:) :: esat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1),size(temp,2),size(temp,3)) ::  work
integer, dimension(size(temp,1),size(temp,2),size(temp,3)) :: iwork
integer  i,j,k
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do k=1,size(temp,3)
      do j=1,size(temp,2)
      do i=1,size(temp,1)
         iwork(i,j,k) = 10.*(temp(i,j,k)+tx05)
          work(i,j,k) = temp(i,j,k)+tx00-0.10*Float(iwork(i,j,k))
         iwork(i,j,k) = 2*iwork(i,j,k)+1
          esat(i,j,k) = es(iwork(i,j,k),work(i,j,k))
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine escomp_3d

!#######################################################################
!#######################################################################

subroutine descomp_0d (temp,desat)

!------------------- public interface arrays ---------------------------
real, intent(IN)  :: temp
real, intent(OUT) :: desat

!------------------- private local storage -----------------------------
   real ::  work
integer :: iwork
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

         iwork=10.*(temp+0.50+tx05)
          work=temp+0.50+tx00-0.10*Float(iwork)
         iwork=2*iwork+1
         desat=es(iwork,work)

         iwork=10.*(temp-0.50+tx05)
          work=temp-0.50+tx00-0.10*Float(iwork)
         iwork=2*iwork+1
         desat=desat-es(iwork,work)

!-----------------------------------------------------------------------

end subroutine descomp_0d

!#######################################################################

subroutine descomp_1d (temp,desat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:) :: temp
real, intent(OUT), dimension(:) :: desat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1)) ::  work
integer, dimension(size(temp,1)) :: iwork
integer  i
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do i=1,size(temp,1)
         iwork(i)=10.*(temp(i)+0.50+tx05)
          work(i)=temp(i)+0.50+tx00-0.10*Float(iwork(i))
         iwork(i)=2*iwork(i)+1
         desat(i)=es(iwork(i),work(i))

         iwork(i)=10.*(temp(i)-0.50+tx05)
          work(i)=temp(i)-0.50+tx00-0.10*Float(iwork(i))
         iwork(i)=2*iwork(i)+1
         desat(i)=desat(i)-es(iwork(i),work(i))
      enddo

!-----------------------------------------------------------------------

end subroutine descomp_1d

!#######################################################################

subroutine descomp_2d (temp,desat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:,:) :: temp
real, intent(OUT), dimension(:,:) :: desat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1),size(temp,2)) :: work
integer, dimension(size(temp,1),size(temp,2)) :: iwork
integer  i,j
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do j=1,size(temp,2)
      do i=1,size(temp,1)
         iwork(i,j)=10.*(temp(i,j)+0.50+tx05)
          work(i,j)=temp(i,j)+0.50+tx00-0.10*Float(iwork(i,j))
         iwork(i,j)=2*iwork(i,j)+1
         desat(i,j)=es(iwork(i,j),work(i,j))

         iwork(i,j)=10.*(temp(i,j)-0.50+tx05)
          work(i,j)=temp(i,j)-0.50+tx00-0.10*Float(iwork(i,j))
         iwork(i,j)=2*iwork(i,j)+1
         desat(i,j)=desat(i,j)-es(iwork(i,j),work(i,j))
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine descomp_2d

!#######################################################################

subroutine descomp_3d (temp,desat)

!------------------- public interface arrays ---------------------------
real, intent(IN),  dimension(:,:,:) :: temp
real, intent(OUT), dimension(:,:,:) :: desat

!------------------- private local storage -----------------------------
   real, dimension(size(temp,1),size(temp,2),size(temp,3)) :: work
integer, dimension(size(temp,1),size(temp,2),size(temp,3)) :: iwork
integer  i,j,k
!-----------------------------------------------------------------------
    integer :: ind
    real    :: dt,es
    es(ind,dt)=ETABL(ind)+dt*ETABL(ind+1)
!-----------------------------------------------------------------------

      If (.not.ESSET_DONE) call ESSET

      do k=1,size(temp,3)
      do j=1,size(temp,2)
      do i=1,size(temp,1)
         iwork(i,j,k)=10.*(temp(i,j,k)+0.50+tx05)
          work(i,j,k)=temp(i,j,k)+0.50+tx00-0.10*Float(iwork(i,j,k))
         iwork(i,j,k)=2*iwork(i,j,k)+1
         desat(i,j,k)=es(iwork(i,j,k),work(i,j,k))

         iwork(i,j,k)=10.*(temp(i,j,k)-0.50+tx05)
          work(i,j,k)=temp(i,j,k)-0.50+tx00-0.10*Float(iwork(i,j,k))
         iwork(i,j,k)=2*iwork(i,j,k)+1
         desat(i,j,k)=desat(i,j,k)-es(iwork(i,j,k),work(i,j,k))
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

end subroutine descomp_3d

!#######################################################################
!#######################################################################

subroutine ESSET

!  =================================================================
!  +                                                               +
!  +             CONSTRUCTION OF THE ES TABLE                      +
!  +                                                               +
!  + THIS TABLE IS CONSTRUCTED FROM ES EQUATIONS FROM THE          +
!  + SMITHSONIAN TABLES.  THE ES INPUT IS COMPUTED FROM VALUES     +
!  + (IN ONE-TENTH OF A DEGREE INCREMENTS) OF ES OVER ICE          +
!  + FROM -153C TO 0C AND VALUES OF ES OVER WATER FROM 0C TO 102C. +
!  + OUTPUT TABLE CONTAINS THESE DATA INTERLEAVED WITH THEIR       +
!  + DERIVATIVES WITH RESPECT TO TEMPERATURE EXCEPT BETWEEN -20C   +
!  + AND 0C WHERE BLENDED (OVER WATER AND OVER ICE) ES VALUES AND  +
!  + DERIVATIVES ARE CALCULATED.                                   +
!  +   NOTE: ALL ES COMPUTATION IS DONE IN MICROBARS AND NOT       +
!  +         MILLIBARS AND THEN CONVERTED TO THE APPROPRIATE       +
!  +         UNITS WHEN PUT INTO ESTABL.                           +
!  =================================================================

real   TABLE(5102),ESNEW(2551),ESUPC(200)
real   ESBASW,TBASW,ESBASI,TBASI,TEM,AA,B,C,D,E,ESH2O,WICE,WH2O
integer  i,ii

      ESBASW = 1013246.0
       TBASW =     373.16
      ESBASI =    6107.1
       TBASI =     273.16

!  COMPUTE ES OVER ICE BETWEEN -153 C AND 0 C.
!  SEE SMITHSONIAN METEOROLOGICAL TABLES PAGE 350.

      Do i=1,1530
         TEM = 120.16+0.1*FLOAT(i-1)
         AA  = -9.09718 *(TBASI/TEM-1.0)
         B   = -3.56654 *LOG10(TBASI/TEM)
         C   =  0.876793*(1.0-TEM/TBASI)
         E   = LOG10(ESBASI)
         TABLE(i)=10**(AA+B+C+E)
      EndDo

!  COMPUTE ES OVER WATER BETWEEN -20C AND FREEZING.
!  SEE SMITHSONIAN METEOROLOGICAL TABLES PAGE 350.

      Do i=1,1221
         TEM = 253.16+0.1*FLOAT(i-1)
         AA  = -7.90298*(TBASW/TEM-1)
         B   =  5.02808*LOG10(TBASW/TEM)
         C   = -1.3816E-07*(10**((1-TEM/TBASW)*11.344)-1)
         D   =  8.1328E-03*(10**((TBASW/TEM-1)*(-3.49149))-1)
         E   = LOG10(ESBASW)
         ESH2O  = 10**(AA+B+C+D+E)
         If (i <= 200) Then
           ESUPC(i) = ESH2O
         Else
           TABLE(i+1330) = ESH2O
         EndIf
      EndDo

!  DERIVE BLENDED ES OVER ICE AND SUPERCOOLED WATER BETWEEN -20C AND 0C

      Do i=1,200
         TEM  = 253.16+0.1*FLOAT(i-1)
         WICE = 0.05*(273.16-TEM)
         WH2O = 0.05*(TEM-253.16)
         TABLE(i+1330) = WICE*TABLE(i+1330)+WH2O*ESUPC(i)
      EndDo

!  ESTIMATE ES DERIVATIVE IN MICROBARS PER DEGREE CELSIUS BY
!    DIFFERENCING ES VALUES IN THE TABLE

      Do i=2,2550
         ESNEW(i) = 5.0*(TABLE(i+1)-TABLE(i-1))
      EndDo
         ESNEW(   1) = ESNEW(   2)
         ESNEW(2551) = ESNEW(2550)

!  INTERLEAVE ES AND ITS DERIVATIVE INTO ONE TABLE OF ALTERNATING
!    VARIABLES (ES ENTRIES ARE ODD NUMBERED, DERIVATIVES ARE EVEN)

      Do i=1,2551
         ii = (i-1)*2
         TABLE(5101-ii) = TABLE(2552-i)
      EndDo

      Do i=1,2551
         ii    = (i-1)*2
         TABLE(ii+2) = ESNEW(i)
      EndDo

! CONVERT TO PASCALS AND PUT INTO COMMON

      Do i=1,5102
         ETABL(i)=TABLE(i)*UNITS
      EndDo

      ESSET_DONE = .TRUE.

end subroutine ESSET

!#######################################################################
!#######################################################################

subroutine tcheck_0d (T,nbad)

!-----------------------------------------------------------------------
   real, intent(IN)  :: T
integer, intent(OUT) :: nbad
!-----------------------------------------------------------------------
          If (T <= Tmin .or. T >= Tmax) nbad=1
!-----------------------------------------------------------------------

end subroutine tcheck_0d

!#######################################################################

subroutine tcheck_1d (T,nbad)

!-----------------------------------------------------------------------
   real, intent(IN), dimension(:) :: T
integer, intent(OUT)              :: nbad
!-----------------------------------------------------------------------
integer  i
!-----------------------------------------------------------------------
      nbad=0
      do i=1,size(T,1)
        If (T(i) <= Tmin .or. T(i) >= Tmax) nbad=nbad+1
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_1d

!#######################################################################

subroutine tcheck_2d (T,nbad)

!-----------------------------------------------------------------------
   real, intent(IN), dimension(:,:) :: T
integer, intent(OUT)                :: nbad
!-----------------------------------------------------------------------
integer  i,j
!-----------------------------------------------------------------------
      nbad=0
      do j=1,size(T,2)
      do i=1,size(T,1)
         If (T(i,j) <= Tmin .or. T(i,j) >= Tmax) nbad=nbad+1
      enddo
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_2d

!#######################################################################

subroutine tcheck_3d (T,nbad)

!-----------------------------------------------------------------------
   real, intent(IN), dimension(:,:,:) :: T
integer, intent(OUT)                  :: nbad
!-----------------------------------------------------------------------
integer  i,j,k
!-----------------------------------------------------------------------
      nbad=0
      do k=1,size(T,3)
      do j=1,size(T,2)
      do i=1,size(T,1)
          If (T(i,j,k) <= Tmin .or. T(i,j,k) >= Tmax) nbad=nbad+1
      enddo
      enddo
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_3d

!#######################################################################

subroutine tcheck_0d_reset (T,nbad, deltat)

!-----------------------------------------------------------------------
real, intent(INOUT)  :: T
real, intent(in)     :: deltat
integer, intent(OUT) :: nbad
!-----------------------------------------------------------------------
          If (T <= Tmin )  then
            nbad=1
            print *, 'bad temp pe,temp', get_my_pe(), T, Tmin+deltat
            T = Tmin + deltat
          else If (T >= Tmax)  then
            nbad=1
            print *, 'bad temp pe,temp', get_my_pe(), T, Tmax-deltat
            T = Tmax - deltat
          endif
!-----------------------------------------------------------------------

end subroutine tcheck_0d_reset

!#######################################################################

subroutine tcheck_1d_reset (T,nbad, deltat)

!-----------------------------------------------------------------------
real, intent(INOUT), dimension(:) :: T
real, intent(in)     :: deltat
integer, intent(OUT)              :: nbad
!-----------------------------------------------------------------------
integer  i
!-----------------------------------------------------------------------
      nbad=0
      do i=1,size(T,1)
	If (T(i) <= Tmin ) then
	  nbad=nbad+1
          print *, 'bad temp pe,i,temp', get_my_pe(), &
                                         i, T(i), Tmin+deltat
	  T(i) = Tmin + deltat 
	else If (T(i) >= Tmax) then
	  nbad=nbad+1
          print *, 'bad temp pe,i,temp', get_my_pe(), &
                                         i, T(i), Tmax-deltat
	  T(i) = Tmax - deltat
	endif
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_1d_reset

!#######################################################################

subroutine tcheck_2d_reset (T,nbad, deltat)

!-----------------------------------------------------------------------
real, intent(INOUT), dimension(:,:) :: T
real, intent(in)     :: deltat
integer, intent(OUT)                :: nbad
!-----------------------------------------------------------------------
integer  i,j
!-----------------------------------------------------------------------
      nbad=0
      do j=1,size(T,2)
      do i=1,size(T,1)
         If (T(i,j) <= Tmin ) then
           nbad=nbad+1
           print *, 'bad temp pe,i,j,temp', get_my_pe(), &
                                            i,j, T(i,j), Tmin+deltat
           T(i,j) = Tmin + deltat
         else If ( T(i,j) >= Tmax) then
           nbad=nbad+1
           print *, 'bad temp pe,i,j,temp', get_my_pe(), &
                                            i,j, T(i,j), Tmax-deltat
           T(i,j) = Tmax - deltat
         endif
      enddo
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_2d_reset

!#######################################################################

subroutine tcheck_3d_reset (T,nbad, deltat)

!-----------------------------------------------------------------------
real, intent(in)     :: deltat
real, intent(INOUT), dimension(:,:,:) :: T
integer, intent(OUT)                  :: nbad
!-----------------------------------------------------------------------
integer  i,j,k
!-----------------------------------------------------------------------
      nbad=0
      do k=1,size(T,3)
      do j=1,size(T,2)
      do i=1,size(T,1)
        If (T(i,j,k) <= Tmin)  then
          nbad=nbad+1
          print *, 'bad temp pe,i,j,k,temp', get_my_pe(), &
                                        i,j,k, T(i,j,k), Tmin+deltat
          T(i,j,k) = Tmin + deltat
        else If ( T(i,j,k) >= Tmax) then
          nbad=nbad+1
          print *, 'bad temp pe,i,j,k,temp', get_my_pe(), &
                                        i,j,k, T(i,j,k), Tmax-deltat
          T(i,j,k) = Tmax - deltat
        endif
      enddo
      enddo
      enddo
!-----------------------------------------------------------------------

end subroutine tcheck_3d_reset

!#######################################################################
!#######################################################################


end module sat_vapor_pres_mod

