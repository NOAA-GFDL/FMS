
module amip_interp_mod

!-----------------------------------------------------------------------

use  time_interp_mod, only: time_interp, fraction_of_year

use time_manager_mod, only: time_type, operator(+), operator(>), &
                            get_date, set_time

use  horiz_interp_mod, only: horiz_interp_init, horiz_interp,  &
                             horiz_interp_end, horiz_interp_type

use    utilities_mod, only: file_exist, error_mesg,      &
                            NOTE, WARNING, FATAL,        &
                            check_nml_error, open_file,  &
                            get_my_pe, close_file

use    constants_mod, only: tfreeze

implicit none
private

!-----------------------------------------------------------------------
!----------------- public interfaces -----------------------------------

public amip_interp_init, get_amip_sst, get_amip_ice,  &
       get_sst_grid_boundary, get_sst_grid_size, amip_interp_end,  &
       amip_interp_type

!-----------------------------------------------------------------------
!--------------------- private below here ------------------------------

!  ---- version number -----

character(len=128) :: version = '$Id: amip_interp.F90,v 1.2 2000/08/04 19:53:49 fms Exp $'
character(len=128) :: tag = '$Name: calgary $'

!-----------------------------------------------------------------------
!------ private defined data type --------

type date_type
   sequence
   integer :: year, month, day
end type

interface operator (==)
   module procedure date_equals
end interface

interface operator (/=)
   module procedure date_not_equals
end interface

!-----------------------------------------------------------------------
!----- public data type ------

type amip_interp_type
   private
   type (horiz_interp_type) :: Hintrp
   real, pointer            ::    data1(:,:), data2(:,:)
   type (date_type)         ::    Date1,       Date2
   logical                  :: use_climo, use_annual
end type

!-----------------------------------------------------------------------
!  ---- resolution/grid variables ----

   integer :: mobs, nobs
   real, allocatable :: lon_bnd(:), lat_bnd(:)
   logical, allocatable :: mask_sst(:,:)

!  ---- global unit & date ----

   integer, parameter :: maxc = 128
   integer :: unit
   character(len=maxc) :: file_name

   type (date_type) :: Curr_date = date_type( -99, -99, -99 )

   real    :: tice_crit_k

   logical :: do_init_once = .true.

!-----------------------------------------------------------------------
!---- namelist ----

!  tice_crit    = freezing point of sea water
!  use_no_climo = flag to prevent monthly or annual mean climatologies
!                 from being used 
!  use_amip1    = flag to select the sst data set:  amip1 or reynolds-oi
!  verbose      = controls printed output
!  use_zonal    = flag to selected zonal sst or data set
!
!  parameters for zonal prescribed sst
!
!  teq  = sst at equator
!  tdif = equator to pole sst difference
!  tann = amplitude of annual cycle
!  tlag = offset for time of year (for annual cycle)


 real    :: tice_crit    = -1.80       !  in degC or degK
 logical :: use_no_climo = .true.
 logical :: use_amip1    = .true.
 logical :: use_zonal    = .false.
 integer :: verbose      = 0           !  0 <= verbose <= 3

 real :: teq  = 305.
 real :: tdif = 50.
 real :: tann = 20.
 real :: tlag = 0.875


 namelist /amip_interp_nml/ tice_crit, use_no_climo, use_amip1,  &
                            use_zonal, verbose, teq, tdif, tann, tlag

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine get_amip_sst (Time, Interp, sst)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   ::  sst(:,:)

    real, dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------

if (use_zonal) then
   call zonal_sst (Time, sice, temp)
   call horiz_interp ( Interp%Hintrp, temp, sst )
else

!-----------------------------------------------------------------------
!---------- get new observed sea surface temperature -------------------

! ---- time interpolation for months -----

   call time_interp (Time, fmonth, year1, year2, month1, month2)
   
! ---- force climatology ----
   if (Interp % use_climo) then
       year1=0; year2=0
   endif
   if (Interp % use_annual) then
        year1=0;  year2=0
       month1=0; month2=0
   endif
! ---------------------------

   Date1 = date_type( year1, month1, 0 )
   Date2 = date_type( year2, month2, 0 )

   unit = -1
!-----------------------------------------------------------------------

    if (Date1 /= Interp % Date1) then
!       ---- use Date2 for Date1 ----
        if (Date1 == Interp % Date2) then
            Interp % Date1 = Interp % Date2
            Interp % data1 = Interp % data2
        else
            call get_sst_record (Date1, Udate1, sice, temp)
            call horiz_interp ( Interp%Hintrp, temp, Interp%data1 )
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

        call get_sst_record (Date2, Udate2, sice, temp)
        call horiz_interp ( Interp%Hintrp, temp, Interp%data2 )
        Interp % Date2 = Date2

    endif

!   ---- if the unit was opened, close it and print dates ----

      if (unit /= -1) then
         call close_file (unit)
         if (verbose > 0 .and. get_my_pe() == 0)         &
                               call print_dates (Time,   &
                                Interp % Date1, Udate1,  &
                                Interp % Date2, Udate2, fmonth)
      endif

!-----------------------------------------------------------------------
!---------- time interpolation (between months) of sst's ---------------
!-----------------------------------------------------------------------

   sst = Interp % data1 + fmonth * (Interp % data2 - Interp % data1)

endif

!-----------------------------------------------------------------------

 end subroutine get_amip_sst

!#######################################################################

subroutine get_amip_ice (Time, Interp, ice)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   :: ice(:,:)

    real, dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------

if (use_zonal) then
   call zonal_sst (Time, sice, temp)
   call horiz_interp ( Interp%Hintrp, sice, ice )
else

!-----------------------------------------------------------------------
!---------- get new observed sea surface temperature -------------------

! ---- time interpolation for months -----

   call time_interp (Time, fmonth, year1, year2, month1, month2)
   
! ---- force climatology ----
   if (Interp % use_climo) then
       year1=0; year2=0
   endif
   if (Interp % use_annual) then
        year1=0;  year2=0
       month1=0; month2=0
   endif
! ---------------------------

   Date1 = date_type( year1, month1, 0 )
   Date2 = date_type( year2, month2, 0 )

   unit = -1
!-----------------------------------------------------------------------

    if (Date1 /= Interp % Date1) then
!       ---- use Date2 for Date1 ----
        if (Date1 == Interp % Date2) then
            Interp % Date1 = Interp % Date2
            Interp % data1 = Interp % data2
        else
            call get_sst_record (Date1, Udate1, sice, temp)
            call horiz_interp ( Interp%Hintrp, sice, Interp%data1 )
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

        call get_sst_record (Date2, Udate2, sice, temp)
        call horiz_interp ( Interp%Hintrp, sice, Interp%data2 )
        Interp % Date2 = Date2

    endif

!   ---- if the unit was opened, close it and print dates ----

      if (unit /= -1) then
         call close_file (unit)
         if (verbose > 0 .and. get_my_pe() == 0)         &
                               call print_dates (Time,   &
                                Interp % Date1, Udate1,  &
                                Interp % Date2, Udate2, fmonth)
      endif

!-----------------------------------------------------------------------
!---------- time interpolation (between months) ------------------------
!-----------------------------------------------------------------------

   ice = Interp % data1 + fmonth * (Interp % data2 - Interp % data1)

endif

!-----------------------------------------------------------------------

 end subroutine get_amip_ice

!#######################################################################

 function amip_interp_init ( blon , blat , mask ,  &
                             use_climo, use_annual ) result (Interp)

 real,    intent(in), dimension(:)   :: blon, blat
 logical, intent(in), dimension(:,:) :: mask
 logical, intent(in), optional       :: use_climo, use_annual

   type (amip_interp_type) :: Interp

   integer :: unit,io,ierr

!-----------------------------------------------------------------------
!----- initialization done once ------

 if (do_init_once) then

!   ---- read namelist ----

    if ( file_exist('input.nml')) then
       unit = open_file (file='input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=amip_interp_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'amip_interp_nml')
       enddo
  10   call close_file (unit)
    endif

!  ----- write namelist/version info -----

    unit = open_file ('logfile.out', action='append')
    if (get_my_pe() == 0) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=amip_interp_nml)
    endif
    call close_file (unit)

!   ---- freezing point of sea water in deg K ---

    tice_crit_k = tice_crit
    if ( tice_crit_k < 200. ) tice_crit_k = tice_crit_k + tfreeze

!   ---- set up file dependent variable ----
!   ----   global file name   ----
!   ----   grid box edges     ----
!   ---- initialize zero size grid if not pe 0 ------

    if (use_amip1) then
        file_name = 'INPUT/' // 'amip1_sst.data'
        mobs = 180;  nobs = 91
        call set_sst_grid_edges_amip1
        if (get_my_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 1 sst', NOTE)
    else
        file_name = 'INPUT/' // 'reyoi_sst.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
        if (get_my_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using Reynolds OI sst', &
                                                                NOTE)
    endif


!  --- check existence of sst data file ??? ---

    if (.not.file_exist(trim(file_name))) then
      call error_mesg ('amip_interp_init',  &
                       'requested input data set does not exist', FATAL)
    endif

    do_init_once = .false.

 endif

!-----------------------------------------------------------------------
!  --- optional arguments ? -----

   Interp % use_climo  = .false.
   if (present(use_climo)) Interp % use_climo  = use_climo
   Interp % use_annual = .false.
   if (present(use_annual)) Interp % use_annual  = use_annual

   if ( use_no_climo .and. Interp%use_climo )  call error_mesg &
                ('amip_interp_init', 'use_climo mismatch', FATAL)

   if ( use_no_climo .and. Interp%use_annual )  call error_mesg &
          ('amip_interp_init', 'use_annual(climo) mismatch', FATAL)

!-----------------------------------------------------------------------
!   ---- initialization of horizontal interpolation ----

    call horiz_interp_init ( Interp%Hintrp, lon_bnd, lat_bnd, &
                             blon, blat)


    Interp % Date1 = date_type( -99, -99, -99 )
    Interp % Date2 = date_type( -99, -99, -99 )


    allocate ( Interp % data1 (size(blon)-1,size(blat)-1), &
               Interp % data2 (size(blon)-1,size(blat)-1)  )


   end function amip_interp_init

!#######################################################################

   subroutine amip_interp_end (Interp)
   type (amip_interp_type), intent(inout) :: Interp

     call horiz_interp_end ( Interp%Hintrp )

   end subroutine amip_interp_end

!#######################################################################

   subroutine set_sst_grid_edges_amip1

   integer :: i, j
   real    :: hpie, dlon, dlat, wb, sb

      allocate ( lon_bnd(mobs+1), lat_bnd(nobs+1) )

! ---- compute grid edges (do only once) -----

      hpie = acos(0.0)

      dlon = 4.*hpie/float(mobs);  wb = -0.5*dlon
      do i = 1, mobs+1
          lon_bnd(i) = wb + dlon * float(i-1)
      enddo
          lon_bnd(mobs+1) = lon_bnd(1) + 4.*hpie

      dlat = 2.*hpie/float(nobs-1);  sb = -hpie + 0.5*dlat
      lat_bnd(1) = -hpie;  lat_bnd(nobs+1) = hpie
      do j = 2, nobs
          lat_bnd(j) = sb + dlat * float(j-2)
      enddo

   end subroutine set_sst_grid_edges_amip1

!#######################################################################

   subroutine set_sst_grid_edges_oi

   integer :: i, j
   real    :: hpie, dlon, dlat, wb, sb

      allocate ( lon_bnd(mobs+1), lat_bnd(nobs+1) )

! ---- compute grid edges (do only once) -----

      hpie = acos(0.0)

      dlon = 4.*hpie/float(mobs);  wb = 0.0
          lon_bnd(1) = wb
      do i = 2, mobs+1
          lon_bnd(i) = wb + dlon * float(i-1)
      enddo
          lon_bnd(mobs+1) = lon_bnd(1) + 4.*hpie

      dlat = 2.*hpie/float(nobs);  sb = -hpie
      lat_bnd(1) = sb;  lat_bnd(nobs+1) = hpie
      do j = 2, nobs
          lat_bnd(j) = sb + dlat * float(j-1)
      enddo

   end subroutine set_sst_grid_edges_oi

!#######################################################################

   subroutine get_sst_grid_size (nlon, nlat)
   integer, intent(out) :: nlon, nlat

      if ( do_init_once ) call error_mesg ('get_sst_grid_size',  &
                         'have not called amip_interp_init', FATAL)

      nlon = mobs;  nlat = nobs

   end subroutine get_sst_grid_size

!#######################################################################

   subroutine get_sst_grid_boundary (blon, blat, mask)

   real,    intent(out) :: blon(:), blat(:)
   logical, intent(out) :: mask(:,:)

      if ( do_init_once ) call error_mesg ('get_sst_grid_boundary',  &
                         'have not called amip_interp_init', FATAL)

! ---- check size of argument(s) ----

      if (size(blon) /= mobs+1 .or. size(blat) /= nobs+1)   &
      call error_mesg ('get_sst_grid_boundary in amip_interp_mod',  &
                       'invalid argument dimensions', FATAL)

! ---- return grid box edges -----

      blon = lon_bnd
      blat = lat_bnd

! ---- masking (data exists at all points) ----

      mask = .true.


   end subroutine get_sst_grid_boundary

!#######################################################################

   subroutine get_sst_record (Date, Adate, ice, sst)

   type (date_type), intent(in)  :: Date
   type (date_type), intent(out) :: Adate
   real,             intent(out) :: ice(mobs,nobs), sst(mobs,nobs)

   integer(2) :: isst (mobs,nobs)
   integer :: yr, mo, dy
   character(len=38) :: mesg

!  ---- rewind condition (if unit is open) ----

    if (unit /= -1 .and. Curr_date % year == 0 .and.   &
         date % month <= Curr_date % month ) then
             if (verbose > 1 .and. get_my_pe() == 0)  &
                             print *, ' rewinding ', trim(file_name)
             rewind unit
    endif

    unit = open_file (file_name, form='ieee32', action='read')

    if (verbose > 2 .and. get_my_pe() == 0)  &
                    print *, 'looking for date = ', Date
    do
      read (unit, end=10)  yr, mo, dy, isst
      Adate = date_type( yr, mo, dy )
      Curr_date = Adate
    if (verbose > 2 .and. get_my_pe() == 0)  &
                    print *, '....... checking   ', Adate

!     --- found date ---
      if (Date == Adate)                    exit

!     --- otherwise use monthly climo ---
      if (yr == 0 .and. Date % month == mo) exit
    enddo

!     --- check if climo used when not wanted ---

      if ( use_no_climo ) then
          if (yr == 0 .or. mo == 0) call error_mesg   &
                ('get_sst_record in amip_interp_mod', &
                 'climo data read when NO climo data requested', FATAL)
      endif

!     --- create (real) ice mask ---

      where ( isst == -200 )
         ice = 1.0
      elsewhere
         ice = 0.0
      endwhere

!     --- unpack sst ---

      sst = float(isst)*0.01 + tfreeze

!     --- set ice mask for sst < critical ----

      where ( sst <= tice_crit ) 
         ice = 1.0
         sst = tice_crit
      endwhere

      return

  10  write (mesg, 20) unit
      call error_mesg ('get_sst_record in amip_interp_mod', mesg, FATAL)
                       
  20  format ('end of file reading unit ',i2,' (sst data)')

   end subroutine get_sst_record

!#######################################################################

function date_equals (Left, Right) result (answer)
type (date_type), intent(in) :: Left, Right
logical :: answer

   if (Left % year  == Right % year  .and.  &
       Left % month == Right % month .and.  &
       Left % day   == Right % day ) then
           answer = .true. 
   else
           answer = .false.
   endif

end function date_equals

!#######################################################################

function date_not_equals (Left, Right) result (answer)
type (date_type), intent(in) :: Left, Right
logical :: answer

   if (Left % year  == Right % year  .and.  &
       Left % month == Right % month .and.  &
       Left % day   == Right % day ) then
           answer = .false.
   else
           answer = .true. 
   endif

end function date_not_equals

!#######################################################################

subroutine print_dates (Time, Date1, Udate1,  &
                              Date2, Udate2, fmonth)

   type (time_type), intent(in) :: Time
   type (date_type), intent(in) :: Date1, Udate1, Date2, Udate2
   real,             intent(in) :: fmonth

   integer :: year, month, day, hour, minute, second

   call get_date (Time, year, month, day, hour, minute, second)

   write (*,10) year,month,day, hour,minute,second
   write (*,20) fmonth
   write (*,30) Date1, Udate1
   write (*,40) Date2, Udate2

10 format (/,' date(y/m/d h:m:s) = ',i4,2('/',i2.2),1x,2(i2.2,':'),i2.2)
20 format (' fmonth = ',f9.7)
30 format (' date1(y/m/d) = ',i4,2('/',i2.2),6x, &
                    'used = ',i4,2('/',i2.2),6x  )
40 format (' date2(y/m/d) = ',i4,2('/',i2.2),6x, &
                    'used = ',i4,2('/',i2.2),6x  )
     
end subroutine print_dates

!#######################################################################

subroutine zonal_sst (Time, ice, sst)

   type (time_type), intent(in)  :: Time
   real,             intent(out) :: ice(mobs,nobs), sst(mobs,nobs)

   real    :: tpi, fdate, eps, ph, sph, sph2, ts
   integer :: j

! namelist needed
!
!  teq  = sst at equator
!  tdif = equator to pole sst difference
!  tann = amplitude of annual cycle
!  tlag = offset for time of year (for annual cycle)
!

    tpi = 4.0*acos(0.0)

    fdate = fraction_of_year (Time)

    eps = sin( tpi*(fdate-tlag) ) * tann

    do j = 1, nobs

        ph  = 0.5*(lat_bnd(j)+lat_bnd(j+1))
       sph  = sin(ph)
       sph2 = sph*sph

       ts = teq - tdif*sph2 - eps*sph

       sst(:,j) = ts

    enddo

    where ( sst < tice_crit )
       ice = 1.0
       sst = tice_crit
    elsewhere
       ice  = 0.0
    endwhere


end subroutine zonal_sst

!#######################################################################

end module amip_interp_mod

