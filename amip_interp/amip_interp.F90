
module amip_interp_mod


! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Provides observed SST and ice mask data sets that have been
!   interpolated onto your model's grid.
! </OVERVIEW>

! <DESCRIPTION>
! Three possible data sets are available:
!
!     1)  <LINK SRC="http://www-pcmdi.llnl.gov/amip">AMIP 1</LINK>        from Jan 1979 to Jan 1989 (2 deg x 2 deg)<BR/>
!     2)  <LINK SRC="amip_interp.rey_oi.txt">Reynolds OI</LINK>   from Nov 1981 to Jan 1999 (1 deg x 1 deg)<BR/>
!     3)  <LINK SRC="ftp://podaac.jpl.nasa.gov/pub/sea_surface_temperature/reynolds/rsst/doc/rsst.html">Reynolds EOF</LINK>  from Jan 1950 to Dec 1998 (2 deg x 2 deg)<BR/><BR/>
!
!     All original data are observed monthly means. This module
!     interpolates linearly in time between pairs of monthly means.
!     Horizontal interpolation is done using the horiz_interp module.
!
!     When a requested date falls outside the range of dates available
!     a namelist option allows for use of the climatological monthly
!     mean values which are computed from all of the data in a particular
!     data set.
! </DESCRIPTION>

! <DATASET NAME="AMIP 1">
!   from Jan 1979 to Jan 1989 (2 deg x 2 deg).
! </DATASET>
! <DATASET NAME="Reynolds OI">
!   from Nov 1981 to Jan 1999 (1 deg x 1 deg)
!             The analysis uses in situ and satellite SST's plus
!             SST's simulated by sea-ice cover.
! </DATASET>
! <DATASET NAME="Reynolds EOF">
!   from Jan 1950 to Dec 1998 (2 deg x 2 deg)
!             NCEP Reynolds Historical Reconstructed Sea Surface Temperature
!             The analysis uses both in-situ SSTs and satellite derived SSTs
!             from the NOAA Advanced Very High Resolution Radiometer.
!             In-situ data is used from 1950 to 1981, while both AVHRR derived
!             satellite SSTs and in-situ data are used from 1981 to the
!             end of 1998.
!
! Note: The data set used by this module have been reformatted as 32-bit IEEE.
!   The data values are packed into 16-bit integers.
!
!   The data sets are read from the following files:
!
!         amip1           INPUT/amip1_sst.data
!         reynolds_io     INPUT/reyoi_sst.data
!         reynolds_eof    INPUT/reynolds_sst.data
! </DATASET>
!-----------------------------------------------------------------------

use  time_interp_mod, only: time_interp, fraction_of_year

use time_manager_mod, only: time_type, operator(+), operator(>), &
                             get_date, set_time, set_date

use  horiz_interp_mod, only: horiz_interp_init, horiz_interp,  &
                             horiz_interp_end, horiz_interp_type

use           fms_mod, only: file_exist, error_mesg, write_version_number,  &
                             NOTE, WARNING, FATAL, stdlog, check_nml_error, &
                             open_namelist_file, open_ieee32_file,          &
                             mpp_pe, close_file, lowercase

use    constants_mod, only: tfreeze, pi

implicit none
private

!-----------------------------------------------------------------------
!----------------- Public interfaces -----------------------------------

public amip_interp_init, get_amip_sst, get_amip_ice,  &
       get_sst_grid_boundary, get_sst_grid_size, amip_interp_end,  &
       amip_interp_type

!-----------------------------------------------------------------------
!--------------------- private below here ------------------------------

!  ---- version number -----

character(len=128) :: version = '$Id: amip_interp.F90,v 1.6 2003/04/09 21:15:25 fms Exp $'
character(len=128) :: tagname = '$Name: inchon $'

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

interface operator (>)
   module procedure date_gt
end interface

! <INTERFACE NAME="amip_interp_init">
!   <OVERVIEW>
!     Function that initializes data needed for the horizontal
!         interpolation between the sst grid and model grid. The 
!         returned variable of type amip_interp_type is needed when
!         calling get_amip_sst and get_amip_ice.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function that initializes data needed for the horizontal
!         interpolation between the sst grid and model grid. The 
!         returned variable of type amip_interp_type is needed when
!         calling get_amip_sst and get_amip_ice.
!   </DESCRIPTION>
!   <IN NAME="blon">
!     Longitude in radians of the model's grid box edges.
!   </IN>
!   <IN NAME="blat">
!     Latitude in radians of the model's grid box edges.
!   </IN>
!   <IN NAME="mask">
!     A mask for the model grid.
!   </IN>
!   <IN NAME="use_climo">
!     Flag the specifies that monthly mean climatological values will be used.
!   </IN>
!   <IN NAME="use_annual">
!     Flag the specifies that the annual mean climatological
!              will be used.  If both use_annual = use_climo = true,
!              then use_annual = true will be used.
!   </IN>
!   <IN NAME="interp_method">
!     specify the horiz_interp scheme. = "conservative" means conservative scheme, 
!     = "bilinear" means  bilinear interpolation.
!   </IN>
!   <OUT NAME="Interp">
!     A defined data type variable needed when calling get_amip_sst and get_amip_ice.
!   </OUT>
!   <TEMPLATE>
!     Interp = amip_interp_init ( blon, blat, mask, use_climo, use_annual, interp_method )
!   </TEMPLATE>

!   <NOTE>
!     This function may be called to initialize multiple variables
!     of type amip_interp_type.  However, there currently is no
!     call to release the storage used by this variable.
!   </NOTE>
!   <NOTE>
!     The size of input augment mask must be a function of the size
!     of input augments blon and blat. The first and second dimensions
!     of mask must equal (size(blon)-1, size(blon)-1).
!   </NOTE>

!   <ERROR MSG="obsolete namelist option USE_AMIP1 should be replaced" STATUS="FATAL">
!     Use namelist variable DATA_SET instead of USE_AMIP1.
!   </ERROR>
!   <ERROR MSG="obsolete namelist option USE_NO_CLIMO should be replaced by option DATE_OUT_OF_RANGE - check documentation " STATUS="FATAL">
!     Use namelist variable DATE_OUT_OF_RANGE instead of USE_NO_CLIMO.
!   </ERROR>
!   <ERROR MSG="the value of the namelist parameter DATA_SET being used is not allowed" STATUS="FATAL">
!     Check the value of namelist variable DATA_SET.
!   </ERROR>
!   <ERROR MSG="requested input data set does not exist" STATUS="FATAL">
!     The data set requested is valid but the data does not exist in
!      the INPUT subdirectory. You may have requested amip2 data which
!      has not been officially set up.
!      See the section on DATA SETS to properly set the data up.
!   </ERROR>
!   <ERROR MSG="use_climo mismatch" STATUS="FATAL">
!     The namelist variable date_out_of_range = 'fail' and the amip_interp_init 
!     argument use_climo = true.  This combination is not allowed.
!   </ERROR>
!   <ERROR MSG="use_annual(climo) mismatch" STATUS="FATAL">
!     The namelist variable date_out_of_range = 'fail' and the amip_interp_init 
!     argument use_annual = true.  This combination is not allowed.
!   </ERROR>
interface amip_interp_init
   module procedure amip_interp_init_1d
   module procedure amip_interp_init_2d
end interface
! </INTERFACE>


!-----------------------------------------------------------------------
!----- public data type ------
! <DATA NAME="amip_interp_type"  TYPE="type (horiz_interp_type)"  >
!   All variables in this data type are PRIVATE. It contains information
!   needed by the interpolation module (exchange_mod) and buffers data.
! </DATA>
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

!  ---- global unit & date ----

   integer, parameter :: maxc = 128
   integer :: unit
   character(len=maxc) :: file_name_sst, file_name_ice

   type (date_type) :: Curr_date = date_type( -99, -99, -99 )
   type (date_type) :: Date_end  = date_type( -99, -99, -99 )

   real            :: tice_crit_k
   integer(kind=2) ::  ice_crit

   logical :: do_init_once = .true.

!-----------------------------------------------------------------------
!---- namelist ----

! <NAMELIST NAME="amip_interp_nml">
!   <DATA NAME="data_set" TYPE="character(len=24)" DEFAULT="data_set = 'amip1'">
!     Name/type of SST data that will be used.
!  <BR/>
!        Possible values (case-insensitive) are: <BR/>
!                          1) amip1<BR/>
!                          2) reynolds_eof<BR/>
!                          3) reynolds_oi<BR/>
!        See the <LINK SRC="amip_interp.html#DATA SETS">data set </LINK>section for more on these data.
!   </DATA>

!   <DATA NAME="date_out_of_range" TYPE="character(len=16)" DEFAULT="date_out_of_range = 'fail'">
!     Controls the use of climatological monthly mean data when
!     the requested date falls outside the range of the data set.<BR/>
!     Possible values are:
!     <PRE>
!   fail      - program will fail if requested date is prior
!               to or after the data set period.
!   initclimo - program uses climatological requested data is
!               prior to data set period and will fail if
!               requested date is after data set period.
!   climo     - program uses climatological data anytime.
!    </PRE>
!   </DATA>

!   <DATA NAME="tice_crit" TYPE="real" DEFAULT="tice_crit = -1.80">
!     Freezing point of sea water in degC or degK.
!   </DATA>
!   <DATA NAME="verbose" TYPE="integer" DEFAULT="verbose = 0">
!     Controls printed output, 0 <= verbose <= 3
!   </DATA>
!   <DATA NAME="use_no_climo" TYPE="logical">
!     Flag to prevent monthly or annual mean climatologies from being used.
!   </DATA>
!   <DATA NAME="use_amip1" TYPE="logical">
!     Flag to select the sst data set:  amip1 or reynolds-oi.
!   </DATA>

!---- additional parameters for controlling zonal prescribed sst ----
!---- these parameters only have an effect when use_zonal=.true. ----
!   <DATA NAME="use_zonal" TYPE="logical" DEFAULT=".false.">
!     Flag to selected zonal sst or data set.
!   </DATA>
!   <DATA NAME="teq" TYPE="real" DEFAULT="teq=305.">
!     sst at the equator.
!   </DATA>
!   <DATA NAME="tdif" TYPE="real" DEFAULT="tdif=50.">
!     Equator to pole sst difference.
!   </DATA>
!   <DATA NAME="tann" TYPE="real" DEFAULT="tann=20.">
!     Amplitude of annual cycle.
!   </DATA>
!   <DATA NAME="tlag" TYPE="real" DEFAULT="tlag=0.875">
!     Offset for time of year (for annual cycle).
!   </DATA>

!   <DATA NAME="amip_date" TYPE="integer(3)" DEFAULT="/-1,-1,-1/">
!     Single calendar date in integer "(year,month,day)" format
!     that is used only if set with year>0, month>0, day>0. 
!     If used, model calendar date is replaced by this date, 
!     but model time of day is still used to determine ice/sst.
!     Used for repeating-single-day (rsd) experiments.
!   </DATA>

!   <DATA NAME="sst_pert" TYPE="real" DEFAULT="sst_pert=0.">
!     Temperature perturbation in degrees Kelvin added onto the SST.
!                The perturbation is globally-uniform (even near sea-ice).
!                It is only used when abs(sst_pert) > 1.e-4.  SST perturbation runs
!                may be useful in accessing model sensitivities.
!   </DATA>
 character(len=24) :: data_set = 'amip1'   !  use 'amip1'
                                           !      'amip2'
                                           !      'reynolds_eof'
                                           !      'reynolds_oi'
					   !      'hurrell'

 character(len=16) :: date_out_of_range = 'fail'  !  use 'fail'
                                                  !      'initclimo'
                                                  !      'climo' 

 real    :: tice_crit    = -1.80       !  in degC or degK
 logical :: use_no_climo = .true.      !  obsolete option
 logical :: use_amip1    = .true.      !  obsolete option
 integer :: verbose      = 0           !  0 <= verbose <= 3

!parameters for prescribed zonal sst option
 logical :: use_zonal    = .false.
 real :: teq  = 305.
 real :: tdif = 50.
 real :: tann = 20.
 real :: tlag = 0.875


!amip date for repeating single day (rsd) option
 integer :: amip_date(3)=(/-1,-1,-1/)

!global temperature perturbation used for sensitivity experiments
 real :: sst_pert = 0.

 namelist /amip_interp_nml/ tice_crit, use_no_climo, use_amip1,  &
                            data_set, date_out_of_range,         &
                            use_zonal, teq, tdif, tann, tlag,    &
                            amip_date, sst_pert, verbose
! </NAMELIST>


!-----------------------------------------------------------------------

contains

!#######################################################################
! <SUBROUTINE NAME="get_amip_sst" INTERFACE="get_amip_sst">
!   <IN NAME="Time" TYPE="time_type" ></IN>
!   <OUT NAME="sst" TYPE="real" DIM="(:,:)"> </OUT>                  
!   <INOUT NAME="Interp" TYPE="amip_interp_type"> </INOUT>
! </SUBROUTINE>

subroutine get_amip_sst (Time, Interp, sst)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   ::  sst(:,:)

    real, dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum


!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------



    if (all(amip_date>0)) then
       call get_date(Time,dum,dum,dum,tod(1),tod(2),tod(3))
       Amip_Time = set_date(amip_date(1),amip_date(2),amip_date(3),tod(1),tod(2),tod(3))
    else
       Amip_Time = Time
    endif
       
       
if (use_zonal) then
   call zonal_sst (Amip_Time, sice, temp)
   call horiz_interp ( Interp%Hintrp, temp, sst )
else

!-----------------------------------------------------------------------
!---------- get new observed sea surface temperature -------------------

! ---- time interpolation for months -----
   call time_interp (Amip_Time, fmonth, year1, year2, month1, month2)
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

!  -- open/rewind file --
   unit = -1
!-----------------------------------------------------------------------

    if (Date1 /= Interp % Date1) then
!       ---- use Date2 for Date1 ----
        if (Date1 == Interp % Date2) then
            Interp % Date1 = Interp % Date2
            Interp % data1 = Interp % data2
        else
            call read_record ('sst', Date1, Udate1, temp)
            call horiz_interp ( Interp%Hintrp, temp, Interp%data1 )
            call clip_data ('sst', Interp%data1)
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

        call read_record ('sst', Date2, Udate2, temp)
        call horiz_interp ( Interp%Hintrp, temp, Interp%data2 )
        call clip_data ('sst', Interp%data2)
        Interp % Date2 = Date2

    endif

!   ---- if the unit was opened, close it and print dates ----

      if (unit /= -1) then
         call close_file (unit)
         if (verbose > 0 .and. mpp_pe() == 0)         &
                               call print_dates (Amip_Time,   &
                                Interp % Date1, Udate1,  &
                                Interp % Date2, Udate2, fmonth)
      endif

!-----------------------------------------------------------------------
!---------- time interpolation (between months) of sst's ---------------
!-----------------------------------------------------------------------

   sst = Interp % data1 + fmonth * (Interp % data2 - Interp % data1)

endif

! add on non-zero sea surface temperature perturbation (namelist option)
! this perturbation may be useful in accessing model sensitivities

   if ( abs(sst_pert) > 0.0001 ) then
      sst = sst + sst_pert
   endif

!-----------------------------------------------------------------------

 end subroutine get_amip_sst



!#######################################################################
! <SUBROUTINE NAME="get_amip_ice" INTERFACE="get_amip_ice">
!   <IN NAME="Time"  TYPE="time_type"  > </IN>
!   <OUT NAME="ice" TYPE="real" DIM="(:,:)"> </OUT>                  
!   <INOUT NAME="Interp" TYPE="amip_interp_type"> </INOUT>
! </SUBROUTINE>

subroutine get_amip_ice (Time, Interp, ice)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   :: ice(:,:)

    real, dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum

!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------


    if (any(amip_date>0)) then

       call get_date(Time,dum,dum,dum,tod(1),tod(2),tod(3))

       Amip_Time = set_date(amip_date(1),amip_date(2),amip_date(3),tod(1),tod(2),tod(3))

    else

       Amip_Time = Time
       
    endif
       

if (use_zonal) then
   call zonal_sst (Amip_Time, sice, temp)
   call horiz_interp ( Interp%Hintrp, sice, ice )
else

!-----------------------------------------------------------------------
!---------- get new observed sea surface temperature -------------------

! ---- time interpolation for months -----

   call time_interp (Amip_Time, fmonth, year1, year2, month1, month2)
   
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
            call read_record ('ice', Date1, Udate1, sice)
            call horiz_interp ( Interp%Hintrp, sice, Interp%data1 )
            call clip_data ('ice', Interp%data1)
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

        call read_record ('ice', Date2, Udate2, sice)
        call horiz_interp ( Interp%Hintrp, sice, Interp%data2 )
        call clip_data ('ice', Interp%data2)
        Interp % Date2 = Date2

    endif

!   ---- if the unit was opened, close it and print dates ----

      if (unit /= -1) then
         call close_file (unit)
         if (verbose > 0 .and. mpp_pe() == 0)         &
                               call print_dates (Amip_Time,   &
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

! <FUNCTION NAME="amip_interp_init_1d" INTERFACE="amip_interp_init">

!   <IN NAME="blon" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="blat" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="mask" TYPE="logical" DIM="(:,:)"> </IN>
!   <IN NAME="use_climo" TYPE="logical" DEFAULT="use_climo = .false."> </IN>
!   <IN NAME="use_annual" TYPE="logical" DEFAULT="use_annual = .false."> </IN>
!   <IN NAME="interp_method" TYPE="character(len=*), optional" DEFAULT="interp_method = conservative"></IN>
!   <OUT NAME="Interp" TYPE="amip_interp_type"> </OUT>

 function amip_interp_init_1d ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)

 real,    intent(in), dimension(:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional       :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (amip_interp_type) :: Interp

   integer :: unit,io,ierr

!-----------------------------------------------------------------------
!----- initialization done once ------

 if (do_init_once) then

!   ---- read namelist ----

    if ( file_exist('input.nml')) then
       unit = open_namelist_file( )
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=amip_interp_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'amip_interp_nml')
       enddo
  10   call close_file (unit)
    endif

!  ----- write namelist/version info -----
    call write_version_number (version, tagname)

    unit = stdlog ( )
    if (mpp_pe() == 0) then
        write (unit,nml=amip_interp_nml)
        write (unit,11)
    endif
    call close_file (unit)

!   print error message if default value of use_amip1 has changed

    if (.not.use_amip1) then
       call error_mesg ('amip_interp_mod',  &
                'obsolete namelist option USE_AMIP1 should be replaced &
                 &by option DATA_SET - check documentation', FATAL)
    endif

!   also for use_no_climo

    if (.not.use_no_climo) then
       call error_mesg ('amip_interp_mod',  &
            'obsolete namelist option USE_NO_CLIMO should be replaced &
            &by option DATE_OUT_OF_RANGE - check documentation', FATAL)
    endif

!   ---- freezing point of sea water in deg K ---

    tice_crit_k = tice_crit
    if ( tice_crit_k < 200. ) tice_crit_k = tice_crit_k + tfreeze
    ice_crit = nint((tice_crit_k-tfreeze)*100.)

!   ---- set up file dependent variable ----
!   ----   global file name   ----
!   ----   grid box edges     ----
!   ---- initialize zero size grid if not pe 0 ------

    if (lowercase(trim(data_set)) == 'amip1') then
        file_name_sst = 'INPUT/' // 'amip1_sst.data'
        file_name_ice = 'INPUT/' // 'amip1_sst.data'
        mobs = 180;  nobs = 91
        call set_sst_grid_edges_amip1
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 1 sst', NOTE)
        Date_end = date_type( 1989, 1, 0 )
    else if (lowercase(trim(data_set)) == 'amip2') then
        file_name_sst = 'INPUT/' // 'amip2_sst.data'
        file_name_ice = 'INPUT/' // 'amip2_ice.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
!       --- specfied min for amip2 ---
        tice_crit_k = 271.38
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 2 sst', NOTE)
        Date_end = date_type( 1996, 3, 0 )
    else if (lowercase(trim(data_set)) == 'hurrell') then
        file_name_sst = 'INPUT/' // 'hurrell_sst.data'
        file_name_ice = 'INPUT/' // 'hurrell_ice.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
!       --- specfied min for hurrell ---
        tice_crit_k = 271.38
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using HURRELL sst', NOTE)
        Date_end = date_type( 2001, 12, 0 )
    else if (lowercase(trim(data_set)) == 'reynolds_eof') then
        file_name_sst = 'INPUT/' // 'reynolds_sst.data'
        file_name_ice = 'INPUT/' // 'reynolds_sst.data'
        mobs = 180;  nobs = 90
        call set_sst_grid_edges_oi
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init',  &
             'using NCEP Reynolds Historical Reconstructed SST', NOTE)
        Date_end = date_type( 1998, 12, 0 )
    else if (lowercase(trim(data_set)) == 'reynolds_oi') then
        file_name_sst = 'INPUT/' // 'reyoi_sst.data'
        file_name_ice = 'INPUT/' // 'reyoi_sst.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using Reynolds OI SST', &
                                                                NOTE)
        Date_end = date_type( 1999, 1, 0 )
    else
        call error_mesg ('amip_interp_init', 'the value of the &
        &namelist parameter DATA_SET being used is not allowed', FATAL)
    endif

    if (verbose > 1 .and. mpp_pe() == 0) &
              print *, 'ice_crit,tice_crit_k=',ice_crit,tice_crit_k


!  --- check existence of sst data file ??? ---

    if (.not.file_exist(trim(file_name_sst)) .or. &
        .not.file_exist(trim(file_name_ice))) then
      call error_mesg ('amip_interp_init',  &
                       'requested input data set does not exist', FATAL)
    endif

    do_init_once = .false.

 11 format (/,'namelist option USE_AMIP1 is no longer supported', &
           /,'use option DATA_SET instead')

 endif

!-----------------------------------------------------------------------
!  --- optional arguments ? -----

   Interp % use_climo  = .false.
   if (present(use_climo)) Interp % use_climo  = use_climo
   Interp % use_annual = .false.
   if (present(use_annual)) Interp % use_annual  = use_annual

   if ( date_out_of_range == 'fail' .and. Interp%use_climo ) &
             call error_mesg ('amip_interp_init',  &
                              'use_climo mismatch', FATAL)

   if ( date_out_of_range == 'fail' .and. Interp%use_annual ) &
             call error_mesg ('amip_interp_init', &
                              'use_annual(climo) mismatch', FATAL)

!-----------------------------------------------------------------------
!   ---- initialization of horizontal interpolation ----

    call horiz_interp_init ( Interp%Hintrp, lon_bnd, lat_bnd, &
                             lon, lat, interp_method= interp_method )


    Interp % Date1 = date_type( -99, -99, -99 )
    Interp % Date2 = date_type( -99, -99, -99 )


    allocate ( Interp % data1 (size(lon)-1,size(lat)-1), &
               Interp % data2 (size(lon)-1,size(lat)-1)  )


   end function amip_interp_init_1d
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="amip_interp_init_2d" INTERFACE="amip_interp_init">
!   <IN NAME="blon" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="blat" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="mask" TYPE="logical" DIM="(:,:)"> </IN>
!   <IN NAME="use_climo" TYPE="logical" DEFAULT="use_climo = .false."> </IN>
!   <IN NAME="use_annual" TYPE="logical" DEFAULT="use_annual = .false."> </IN>
!   <IN NAME="interp_method" TYPE="character(len=*), optional" DEFAULT="interp_method = conservative "></IN>
!   <OUT NAME="Interp" TYPE="amip_interp_type"> </OUT>

 function amip_interp_init_2d ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)

 real,    intent(in), dimension(:,:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (amip_interp_type) :: Interp

   integer :: unit,io,ierr

!-----------------------------------------------------------------------
!----- initialization done once ------

 if (do_init_once) then

!   ---- read namelist ----

    if ( file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=amip_interp_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'amip_interp_nml')
       enddo
  10   call close_file (unit)
    endif

!  ----- write namelist/version info -----
    call write_version_number(version,tagname)
    unit = stdlog ( )
    if (mpp_pe() == 0) then
        write (unit,nml=amip_interp_nml)
        write (unit,11)
    endif
    call close_file (unit)

!   print error message if default value of use_amip1 has changed

    if (.not.use_amip1) then
       call error_mesg ('amip_interp_mod',  &
                'obsolete namelist option USE_AMIP1 should be replaced &
                 &by option DATA_SET - check documentation', FATAL)
    endif

!   also for use_no_climo

    if (.not.use_no_climo) then
       call error_mesg ('amip_interp_mod',  &
            'obsolete namelist option USE_NO_CLIMO should be replaced &
            &by option DATE_OUT_OF_RANGE - check documentation', FATAL)
    endif

!   ---- freezing point of sea water in deg K ---

    tice_crit_k = tice_crit
    if ( tice_crit_k < 200. ) tice_crit_k = tice_crit_k + tfreeze
    ice_crit = nint((tice_crit_k-tfreeze)*100.)

!   ---- set up file dependent variable ----
!   ----   global file name   ----
!   ----   grid box edges     ----
!   ---- initialize zero size grid if not pe 0 ------

    if (lowercase(trim(data_set)) == 'amip1') then
        file_name_sst = 'INPUT/' // 'amip1_sst.data'
        file_name_ice = 'INPUT/' // 'amip1_sst.data'
        mobs = 180;  nobs = 91
        call set_sst_grid_edges_amip1
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 1 sst', NOTE)
        Date_end = date_type( 1989, 1, 0 )
    else if (lowercase(trim(data_set)) == 'amip2') then
        file_name_sst = 'INPUT/' // 'amip2_sst.data'
        file_name_ice = 'INPUT/' // 'amip2_ice.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
!       --- specfied min for amip2 ---
        tice_crit_k = 271.38
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 2 sst', NOTE)
        Date_end = date_type( 1996, 3, 0 )
    else if (lowercase(trim(data_set)) == 'hurrell') then
        file_name_sst = 'INPUT/' // 'hurrell_sst.data'
        file_name_ice = 'INPUT/' // 'hurrell_ice.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
!       --- specfied min for hurrell ---
        tice_crit_k = 271.38
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using HURRELL sst', NOTE)
        Date_end = date_type( 2001, 12, 0 )
    else if (lowercase(trim(data_set)) == 'reynolds_eof') then
        file_name_sst = 'INPUT/' // 'reynolds_sst.data'
        file_name_ice = 'INPUT/' // 'reynolds_sst.data'
        mobs = 180;  nobs = 90
        call set_sst_grid_edges_oi
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init',  &
             'using NCEP Reynolds Historical Reconstructed SST', NOTE)
        Date_end = date_type( 1998, 12, 0 )
    else if (lowercase(trim(data_set)) == 'reynolds_oi') then
        file_name_sst = 'INPUT/' // 'reyoi_sst.data'
        file_name_ice = 'INPUT/' // 'reyoi_sst.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using Reynolds OI SST', &
                                                                NOTE)
        Date_end = date_type( 1999, 1, 0 )
    else
        call error_mesg ('amip_interp_init', 'the value of the &
        &namelist parameter DATA_SET being used is not allowed', FATAL)
    endif

    if (verbose > 1 .and. mpp_pe() == 0) &
              print *, 'ice_crit,tice_crit_k=',ice_crit,tice_crit_k


!  --- check existence of sst data file ??? ---

    if (.not.file_exist(trim(file_name_sst)) .or. &
        .not.file_exist(trim(file_name_ice))) then
      call error_mesg ('amip_interp_init',  &
                       'requested input data set does not exist', FATAL)
    endif

    do_init_once = .false.

 11 format (/,'namelist option USE_AMIP1 is no longer supported', &
           /,'use option DATA_SET instead')

 endif

!-----------------------------------------------------------------------
!  --- optional arguments ? -----

   Interp % use_climo  = .false.
   if (present(use_climo)) Interp % use_climo  = use_climo
   Interp % use_annual = .false.
   if (present(use_annual)) Interp % use_annual  = use_annual

   if ( date_out_of_range == 'fail' .and. Interp%use_climo ) &
             call error_mesg ('amip_interp_init',  &
                              'use_climo mismatch', FATAL)

   if ( date_out_of_range == 'fail' .and. Interp%use_annual ) &
             call error_mesg ('amip_interp_init', &
                              'use_annual(climo) mismatch', FATAL)

!-----------------------------------------------------------------------
!   ---- initialization of horizontal interpolation ----

    call horiz_interp_init ( Interp%Hintrp, lon_bnd, lat_bnd, &
                             lon, lat, interp_method = interp_method)


    Interp % Date1 = date_type( -99, -99, -99 )
    Interp % Date2 = date_type( -99, -99, -99 )


    allocate ( Interp % data1 (size(lon,1),size(lat,2)), &
               Interp % data2 (size(lon,1),size(lat,2))) 


   end function amip_interp_init_2d
! </FUNCTION>

!#######################################################################

! <SUBROUTINE NAME="amip_interp_end">

!   <OVERVIEW>
!     Call this routine for all amip_interp_type variables created by amip_interp_init.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Call this routine for all amip_interp_type variables created by amip_interp_init.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call amip_interp_end (Interp)
!   </TEMPLATE>
!   <INOUT NAME="Interp" TYPE="amip_interp_type">
!     A defined data type variable initialized by amip_interp_init
!            and used when calling get_amip_sst and get_amip_ice.
!   </INOUT>

   subroutine amip_interp_end (Interp)
   type (amip_interp_type), intent(inout) :: Interp

     deallocate(Interp%data1, Interp%data2, lon_bnd, lat_bnd)
     call horiz_interp_end ( Interp%Hintrp )

   end subroutine amip_interp_end
! </SUBROUTINE>

!#######################################################################

   subroutine set_sst_grid_edges_amip1

   integer :: i, j
   real    :: hpie, dlon, dlat, wb, sb

      allocate ( lon_bnd(mobs+1), lat_bnd(nobs+1) )

! ---- compute grid edges (do only once) -----

      hpie = 0.5*pi

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

      hpie = 0.5*pi

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

! <SUBROUTINE NAME="get_sst_grid_size">

!   <OVERVIEW>
!     Returns the size (i.e., number of longitude and latitude
!         points) of the observed data grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns the size (i.e., number of longitude and latitude
!         points) of the observed data grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_sst_grid_size (nlon, nlat)
!   </TEMPLATE>
!   <OUT NAME="nlon" TYPE="integer">
!     The number of longitude points (first dimension) in the
!        observed data grid.  For AMIP 1 nlon = 180, and the Reynolds nlon = 360.
!   </OUT>
!   <OUT NAME="nlat" TYPE="integer">
!     The number of latitude points (second dimension) in the
!        observed data grid.  For AMIP 1 nlon = 91, and the Reynolds nlon = 180.
!   </OUT>
!   <ERROR MSG="have not called amip_interp_init" STATUS="FATAL">
!     Must call amip_interp_init before get_sst_grid_size.
!   </ERROR>

   subroutine get_sst_grid_size (nlon, nlat)

   integer, intent(out) :: nlon, nlat

      if ( do_init_once ) call error_mesg ('get_sst_grid_size',  &
                         'have not called amip_interp_init', FATAL)

      nlon = mobs;  nlat = nobs

   end subroutine get_sst_grid_size
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="get_sst_grid_boundary">

!   <OVERVIEW>
!     Returns the grid box boundaries of the observed data grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Returns the grid box boundaries of the observed data grid.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_sst_grid_boundary (blon, blat, mask)
!   </TEMPLATE>
!   <OUT NAME="blon" TYPE="real" DIM="(:)">
!     The grid box edges (in radians) for longitude points of the
!        observed data grid. The size of this argument must be nlon+1.
!   </OUT>
!   <OUT NAME="blat" TYPE="real" DIM="(:)">
!     The grid box edges (in radians) for latitude points of the
!        observed data grid. The size of this argument must be nlat+1.
!   </OUT>
!   <ERROR MSG="have not called amip_interp_init" STATUS="FATAL">
!     Must call amip_interp_init before get_sst_grid_boundary.
!   </ERROR>
!   <ERROR MSG="invalid argument dimensions" STATUS="FATAL">
!     The size of the output argument arrays do not agree with
!      the size of the observed data. See the documentation for
!      interfaces get_sst_grid_size and get_sst_grid_boundary.
!   </ERROR>

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
! </SUBROUTINE>

!#######################################################################

   subroutine read_record (type, Date, Adate, dat)

   character(len=*), intent(in)  :: type
   type (date_type), intent(in)  :: Date
   type (date_type), intent(out) :: Adate
   real,             intent(out) :: dat(mobs,nobs)

   integer(2) :: idat (mobs,nobs)
   integer :: yr, mo, dy, ierr
   character(len=38) :: mesg

!  ---- rewind condition (if unit is open) ----

    if (unit /= -1 .and. Curr_date % year == 0 .and.   &
         date % month <= Curr_date % month ) then
             if (verbose > 1 .and. mpp_pe() == 0)  &
                             print *, ' rewinding unit = ', unit
             rewind unit
    endif

    if (unit == -1) then
       if (type(1:3) == 'sst') then
          unit = open_ieee32_file (file_name_sst, 'read')
       else if (type(1:3) == 'ice') then
          unit = open_ieee32_file (file_name_ice, 'read')
       endif
    endif
    dy = 0

    if (verbose > 2 .and. mpp_pe() == 0)  &
                    print *, 'looking for date = ', Date
    do
      if (lowercase(trim(data_set)) == 'amip2' .or. lowercase(trim(data_set)) == 'hurrell') then
         read (unit, end=10)  yr, mo,      dat
      else
         read (unit, end=10)  yr, mo, dy, idat
!new     read (unit, end=10)  yr, mo, dy
      endif
      Adate = date_type( yr, mo, dy )
      Curr_date = Adate
    if (verbose > 2 .and. mpp_pe() == 0)  &
                    print *, '....... checking   ', Adate

!     --- found date ---
      if (Date == Adate)                    exit

!     --- otherwise use monthly climo ---
      if (yr == 0 .and. Date % month == mo) exit

!     --- skip this data record ---
!new  if (lowercase(trim(data_set)) /= 'amip2') read (unit)
    enddo

!     --- read data ---
!new  if (lowercase(trim(data_set)) /= 'amip2') read (unit) idat

!   --- check if climo used when not wanted ---

    if (yr == 0 .or. mo == 0) then
        ierr = 0
        if (date_out_of_range == 'fail' )               ierr = 1
        if (date_out_of_range == 'initclimo' .and.  &
                                    Date > Date_end )   ierr = 1
        if (ierr /= 0) call error_mesg &
               ('read_record in amip_interp_mod', &
                'climo data read when NO climo data requested', FATAL)
    endif

!     --- create (real) ice mask ---

      if (type(1:3) == 'ice') then
          if (lowercase(trim(data_set)) /= 'amip2' .and. lowercase(trim(data_set)) /= 'hurrell') then
              where ( idat <= ice_crit )
                 dat = 1.
              elsewhere
                 dat = 0.
              endwhere
          else
              dat = dat*0.01
          endif
      else if (type(1:3) == 'sst') then
!         --- unpack sst ---
          if (lowercase(trim(data_set)) /= 'amip2' .and. lowercase(trim(data_set)) /= 'hurrell') &
              dat = float(idat)*0.01 + tfreeze
      endif


      return

  10  write (mesg, 20) unit
      call error_mesg ('read_record in amip_interp_mod', mesg, FATAL)

  20  format ('end of file reading unit ',i2,' (sst data)')

   end subroutine read_record

!#######################################################################

   subroutine clip_data (type, dat)

   character(len=*), intent(in)    :: type
   real,             intent(inout) :: dat(:,:)

   if (type(1:3) == 'ice') then
       dat = min(max(dat,0.0),1.0)
   else if (type(1:3) == 'sst') then
       dat = max(tice_crit_k,dat)
   endif

   end subroutine clip_data

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

function date_gt (Left, Right) result (answer)
type (date_type), intent(in) :: Left, Right
logical :: answer
integer :: i, dif(3)

   dif(1) = Left%year  - Right%year
   dif(2) = Left%month - Right%month
   dif(3) = Left%day   - Right%day
   answer = .false.
   do i = 1, 3
     if (dif(i) == 0) cycle
     if (dif(i)  < 0) exit
     if (dif(i)  > 0) then
         answer = .true.
         exit
     endif
   enddo

end function date_gt

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

    tpi = 2.0*pi

    fdate = fraction_of_year (Time)

    eps = sin( tpi*(fdate-tlag) ) * tann

    do j = 1, nobs

        ph  = 0.5*(lat_bnd(j)+lat_bnd(j+1))
       sph  = sin(ph)
       sph2 = sph*sph

       ts = teq - tdif*sph2 - eps*sph

       sst(:,j) = ts

    enddo

    where ( sst < tice_crit_k )
       ice = 1.0
       sst = tice_crit_k
    elsewhere
       ice  = 0.0
    endwhere


end subroutine zonal_sst

!#######################################################################

end module amip_interp_mod
! <INFO>

!   <FUTURE> 
!     Add AMIP 2 data set.
!
!     Other data sets (or extend current data sets).
!   </FUTURE>

! </INFO>
