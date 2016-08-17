
module amip_interp_mod


! <CONTACT EMAIL="Bruce.Wyman@noaa.gov">
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

! add by JHC
use get_cal_time_mod, only: get_cal_time
use mpp_io_mod, only : mpp_open, mpp_read, MPP_RDONLY, MPP_NETCDF, & 
                       MPP_MULTI, MPP_SINGLE, mpp_close, mpp_get_times
! end add by JHC

use  horiz_interp_mod, only: horiz_interp_init, horiz_interp,  &
                             horiz_interp_new, horiz_interp_del, &
                             horiz_interp_type, assignment(=)

use           fms_mod, only: file_exist, error_mesg, write_version_number,  &
                             NOTE, WARNING, FATAL, stdlog, check_nml_error, &
                             open_namelist_file, open_ieee32_file,          &
                             mpp_pe, close_file, lowercase, mpp_root_pe,    &
                             NOTE, mpp_error, fms_error_handler
use        fms_io_mod, only: read_data, field_size ! add by JHC     
use     constants_mod, only: TFREEZE, pi
use      platform_mod, only: R4_KIND, I2_KIND
use mpp_mod,           only: input_nml_file

implicit none
private

!-----------------------------------------------------------------------
!----------------- Public interfaces -----------------------------------

public amip_interp_init, get_amip_sst, get_amip_ice, amip_interp_new, &
       amip_interp_del, amip_interp_type, assignment(=)

!-----------------------------------------------------------------------
!----------------- Public Data -----------------------------------
integer :: i_sst = 1200
integer :: j_sst = 600
real, parameter:: big_number = 1.E30
logical :: forecast_mode = .false.
real, allocatable, dimension(:,:) ::  sst_ncep, sst_anom

public i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst

!-----------------------------------------------------------------------
!--------------------- private below here ------------------------------

!  ---- version number -----

! Include variable "version" to be written to log file.
#include<file_version.h>

   real, allocatable:: temp1(:,:), temp2(:,:)
! add by JHC
   real, allocatable, dimension(:,:) :: tempamip
! end add by JHC
!-----------------------------------------------------------------------
!------ private defined data type --------

type date_type
   sequence
   integer :: year, month, day
end type

interface assignment(=)
  module procedure  amip_interp_type_eq
end interface

interface operator (==)
   module procedure date_equals
end interface

interface operator (/=)
   module procedure date_not_equals
end interface

interface operator (>)
   module procedure date_gt
end interface

! <INTERFACE NAME="amip_interp_new">
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
!   <IN NAME="lon">
!     Longitude in radians of the model's grid box edges (1d lat/lon grid case)
!     or at grid box mid-point (2d case for arbitrary grids).
!   </IN>
!   <IN NAME="lat">
!     Latitude in radians of the model's grid box edges (1d lat/lon grid case)
!     or at grid box mid-point (2d case for arbitrary grids).
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
!     Interp = amip_interp_new ( lon, lat, mask, use_climo, use_annual, interp_method )
!   </TEMPLATE>

!   <NOTE>
!     This function may be called to initialize multiple variables
!     of type amip_interp_type.  However, there currently is no
!     call to release the storage used by this variable.
!   </NOTE>
!   <NOTE>
!     The size of input augment mask must be a function of the size
!     of input augments lon and lat. The first and second dimensions
!     of mask must equal (size(lon,1)-1, size(lat,2)-1).
!   </NOTE>

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
!     The namelist variable date_out_of_range = 'fail' and the amip_interp_new
!     argument use_climo = true.  This combination is not allowed.
!   </ERROR>
!   <ERROR MSG="use_annual(climo) mismatch" STATUS="FATAL">
!     The namelist variable date_out_of_range = 'fail' and the amip_interp_new
!     argument use_annual = true.  This combination is not allowed.
!   </ERROR>
interface amip_interp_new
   module procedure amip_interp_new_1d
   module procedure amip_interp_new_2d
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
   type (horiz_interp_type) :: Hintrp, Hintrp2 ! add by JHC
   real, pointer            ::    data1(:,:) =>NULL(), &
                                  data2(:,:) =>NULL()
   type (date_type)         ::    Date1,       Date2
   logical                  :: use_climo, use_annual
   logical                  :: I_am_initialized=.false.
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

   real             :: tice_crit_k
   integer(I2_KIND) ::  ice_crit

   logical :: module_is_initialized = .false.

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
                                           ! add by JHC:  
                                           !      'daily', when "use_daily=.T."

 character(len=16) :: date_out_of_range = 'fail'  !  use 'fail'
                                                  !      'initclimo'
                                                  !      'climo' 

 real    :: tice_crit    = -1.80       !  in degC or degK
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

! add by JHC
 character(len=6) :: sst_pert_type = 'fixed'  ! use 'random' or 'fixed'
 logical :: do_sst_pert = .false. 
 logical :: use_daily = .false. ! if '.true.', give 'data_set = 'daily''
! end add by JHC

! SJL: During nudging:   use_ncep_sst = .T.;  no_anom_sst = .T.
!      during forecast:  use_ncep_sst = .T.;  no_anom_sst = .F.
! For seasonal forecast: use_ncep_ice = .F.

 logical :: use_ncep_sst = .false.
 logical ::  no_anom_sst = .true.
 logical :: use_ncep_ice = .false.
 logical :: interp_oi_sst = .false.       ! changed to false for regular runs

 namelist /amip_interp_nml/ use_ncep_sst, no_anom_sst, use_ncep_ice,  tice_crit, &
                            interp_oi_sst, data_set, date_out_of_range,          &
                            use_zonal, teq, tdif, tann, tlag, amip_date,         &
                            ! add by JHC
                            sst_pert, sst_pert_type, do_sst_pert,                &
                            use_daily,                                           &
                            ! end add by JHC
                            verbose, i_sst, j_sst, forecast_mode
! </NAMELIST>


!-----------------------------------------------------------------------

contains

!#######################################################################
! <SUBROUTINE NAME="get_amip_sst" INTERFACE="get_amip_sst">
!   <IN NAME="Time" TYPE="time_type" ></IN>
!   <OUT NAME="sst" TYPE="real" DIM="(:,:)"> </OUT>                  
!   <INOUT NAME="Interp" TYPE="amip_interp_type"> </INOUT>
! </SUBROUTINE>

! modified by JHC
subroutine get_amip_sst (Time, Interp, sst, err_msg, lon_model, lat_model)
!subroutine get_amip_sst (Time, Interp, sst, err_msg)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   ::  sst(:,:)
   character(len=*), optional, intent(out) :: err_msg

   real, dimension(mobs,nobs) :: sice

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum(3)

! add by JHC
    real,    intent(in), dimension(:,:), optional :: lon_model, lat_model
    real :: pert
    integer :: i, j, mobs_sst, nobs_sst 
    integer :: jhctod(6)
    type (time_type) :: Udate
    character(len=4) :: yyyy
    integer :: nrecords, ierr, k, yr, mo, dy
    integer :: siz(4)
    integer, dimension(:), allocatable :: ryr, rmo, rdy
    character(len=30) :: time_unit
    real, dimension(:), allocatable :: timeval
    character(len=maxc) :: ncfilename
! end add by JHC


    if(present(err_msg)) err_msg = ''
    if(.not.Interp%I_am_initialized) then
      if(fms_error_handler('get_amip_sst','The amip_interp_type variable is not initialized',err_msg)) return
    endif

!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------

    if ( use_ncep_sst .and. forecast_mode ) no_anom_sst = .false.

    if (all(amip_date>0)) then
       call get_date(Time,dum(1),dum(2),dum(3),tod(1),tod(2),tod(3))
       Amip_Time = set_date(amip_date(1),amip_date(2),amip_date(3),tod(1),tod(2),tod(3))
    else
       Amip_Time = Time
    endif

! add by JHC
if ( .not.use_daily ) then 
! end add by JHC

   if ( .not. allocated(temp1) ) allocate (temp1(mobs,nobs))
   if ( .not. allocated(temp2) ) allocate (temp2(mobs,nobs))

   if (use_zonal) then
      call zonal_sst (Amip_Time, sice, temp1)
      call horiz_interp ( Interp%Hintrp, temp1, sst )
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
              temp1(:,:) = temp2(:,:)   ! SJL BUG fix: June 24, 2011
          else
              call read_record ('sst', Date1, Udate1, temp1)
              if ( use_ncep_sst .and. (.not. no_anom_sst) ) then
                   temp1(:,:) = temp1(:,:) + sst_anom(:,:)
              endif
              call horiz_interp ( Interp%Hintrp, temp1, Interp%data1 )
              call clip_data ('sst', Interp%data1)
             Interp % Date1 = Date1
          endif
      endif

!-----------------------------------------------------------------------

      if (Date2 /= Interp % Date2) then
          call read_record ('sst', Date2, Udate2, temp2)
          if ( use_ncep_sst .and. (.not. no_anom_sst) ) then
               temp2(:,:) = temp2(:,:) + sst_anom(:,:)
          endif
          call horiz_interp ( Interp%Hintrp, temp2, Interp%data2 )
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

!-------------------------------------------------------------------------------
! SJL mods for NWP and TCSF ---
!      Nudging runs: (Note: NCEP SST updated only every 6-hr)
!      Compute SST anomaly from global SST datasets for subsequent forecast runs
!-------------------------------------------------------------------------------
    if ( use_ncep_sst .and. no_anom_sst ) then
         sst_anom(:,:) = sst_ncep(:,:) - (temp1(:,:) + fmonth*(temp2(:,:) - temp1(:,:)) )
         call horiz_interp ( Interp%Hintrp, sst_ncep, sst )
         call clip_data ('sst', sst)
    endif

!!! DEBUG CODE
!          call get_date(Amip_Time,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
!          if (mpp_pe() == 0) then 
!             write (*,200) 'JHC: use_daily = F, AMIP_Time: ',jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6)
!             write (*,300) 'JHC: use_daily = F, interped SST: ', sst(1,1),sst(5,5),sst(10,10)
!          endif
!!! END DEBUG CODE


  endif

! add by JHC
else
    call get_date(Amip_Time,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
     if (mpp_pe() == mpp_root_pe()) write(*,200) 'amip_interp_mod: use_daily = T, Amip_Time = ',jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6) 
 
    yr = jhctod(1); mo = jhctod(2); dy = jhctod(3)
   
    write (yyyy,'(i4)') jhctod(1)

    file_name_sst = 'INPUT/' // 'sst.day.mean.'//yyyy//'.v2.nc'
    ncfilename = trim(file_name_sst)
    time_unit = 'days since 1978-01-01 00:00:00'

    mobs_sst = 1440;  nobs_sst = 720 
 
    call set_sst_grid_edges_daily(mobs_sst, nobs_sst)
    call horiz_interp_new ( Interp%Hintrp2, lon_bnd, lat_bnd, &
                             lon_model, lat_model, interp_method="bilinear" )


    if ( (.NOT. file_exist(ncfilename))  ) call mpp_error ('amip_interp_mod', &
             'cannot find daily SST input data file: '//trim(ncfilename), NOTE)

    if (file_exist(ncfilename)) then
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('amip_interp_mod', &
             'Reading NetCDF formatted daily SST from: '//trim(ncfilename), NOTE)

        call field_size(ncfilename, 'TIME', siz)
        nrecords = siz (1)
        if (nrecords < 1) call mpp_error('amip_interp_mod', &
                           'Invalid number of SST records in daily SST data file: '//trim(ncfilename), FATAL)
        allocate(timeval(nrecords), ryr(nrecords), rmo(nrecords), rdy(nrecords))

        call mpp_open( unit, ncfilename, MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE )
        call mpp_get_times(unit, timeval)
        call mpp_close(unit)

!!! DEBUG CODE
!          if (mpp_pe() == 0) then 
!             print *, 'JHC: nrecords = ', nrecords
!             print *, 'JHC: TIME = ', timeval
!          endif
!!! END DEBUG CODE

        ierr = 1
        do k = 1, nrecords

          Udate = get_cal_time (timeval(k), time_unit, 'julian')
          call get_date(Udate,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
          ryr(k) = jhctod(1); rmo(k) = jhctod(2); rdy(k) = jhctod(3)

          if ( yr == ryr(k) .and. mo == rmo(k) .and. dy == rdy (k) ) ierr = 0
          if (ierr==0) exit

        enddo
!!! DEBUG CODE
             if (mpp_pe() == 0) then 
             print *, 'JHC: k =', k 
             print *, 'JHC: ryr(k) rmo(k) rdy(k)',ryr(k), rmo(k), rdy(k) 
             print *, 'JHC:  yr     mo     dy   ',yr, mo, dy 
          endif
!!! END DEBUG CODE
        if (ierr .ne. 0) call mpp_error('amip_interp_mod', &
                         'Model time is out of range not in SST data: '//trim(ncfilename), FATAL)
    endif ! if(file_exist(ncfilename))
     

   !---- read NETCDF data ----
     if ( .not. allocated(tempamip) ) allocate (tempamip(mobs_sst,nobs_sst))

     if (file_exist(ncfilename)) then
          call read_data(ncfilename, 'SST', tempamip, timelevel=k, no_domain=.true.)
          tempamip = tempamip + TFREEZE

!!! DEBUG CODE
!          if (mpp_pe() == 0) then 
!             print*, 'JHC: TFREEZE = ', TFREEZE
!             print*, lbound(sst)
!             print*, ubound(sst)
!             print*, lbound(tempamip)
!             print*, ubound(tempamip)
!             write(*,300) 'JHC: tempamip : ', tempamip(100,100), tempamip(200,200), tempamip(300,300)
!          endif
!!! END DEBUG CODE

          call horiz_interp ( Interp%Hintrp2, tempamip, sst )
          call clip_data ('sst', sst)

     endif 

!!! DEBUG CODE
!          if (mpp_pe() == 400) then 
!             write(*,300)'JHC: use_daily = T, daily SST: ', sst(1,1),sst(5,5),sst(10,10)
!             print *,'JHC: use_daily = T, daily SST: ', sst
!          endif
!!! END DEBUG CODE

200 format(a35, 6(i5,1x))
300 format(a35, 3(f7.3,2x))

endif
! end add by JHC

! add by JHC: add on non-zero sea surface temperature perturbation (namelist option)
!             This perturbation may be useful in accessing model sensitivities

 if ( do_sst_pert ) then
      
      if ( trim(sst_pert_type) == 'fixed' ) then
          sst = sst + sst_pert
      else if ( trim(sst_pert_type) == 'random' ) then
          call random_seed()
!!! DEBUG CODE
!    	  if (mpp_pe() == 0) then
!             print*, 'mobs = ', mobs	
!             print*, 'nobs = ', nobs	
!             print*, lbound(sst)
!             print*, ubound(sst)
!          endif
!!! END DEBUG CODE
          do i = 1, size(sst,1)
          do j = 1, size(sst,2)
             call random_number(pert)
             sst (i,j) = sst (i,j) + sst_pert*((pert-0.5)*2)
          end do
          end do
      endif

  endif
! end add by JHC
   
!-----------------------------------------------------------------------

 end subroutine get_amip_sst


!#######################################################################
! <SUBROUTINE NAME="get_amip_ice" INTERFACE="get_amip_ice">
!   <IN NAME="Time"  TYPE="time_type"  > </IN>
!   <OUT NAME="ice" TYPE="real" DIM="(:,:)"> </OUT>                  
!   <INOUT NAME="Interp" TYPE="amip_interp_type"> </INOUT>
! </SUBROUTINE>

subroutine get_amip_ice (Time, Interp, ice, err_msg)

   type (time_type),         intent(in)    :: Time
   type (amip_interp_type),  intent(inout) :: Interp
   real,                     intent(out)   :: ice(:,:)
   character(len=*), optional, intent(out) :: err_msg

    real, dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum(3)

    if(present(err_msg)) err_msg = ''
    if(.not.Interp%I_am_initialized) then
      if(fms_error_handler('get_amip_ice','The amip_interp_type variable is not initialized',err_msg)) return
    endif

!-----------------------------------------------------------------------
!----- compute zonally symetric sst ---------------


    if (any(amip_date>0)) then

       call get_date(Time,dum(1),dum(2),dum(3),tod(1),tod(2),tod(3))

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
!-- SJL -------------------------------------------------------------
! Can NOT use ncep_sst to determine sea_ice For seasonal forecast
! Use climo sea ice for seasonal runs
            if ( use_ncep_sst .and. use_ncep_ice ) then
               where ( sst_ncep <= (TFREEZE+tice_crit) )
                   sice = 1.
               elsewhere
                   sice = 0.
               endwhere
            else
               call read_record ('ice', Date1, Udate1, sice)
            endif
!--------------------------------------------------------------------
            call horiz_interp ( Interp%Hintrp, sice, Interp%data1 )
            call clip_data ('ice', Interp%data1)
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

!-- SJL -------------------------------------------------------------
            if ( use_ncep_sst .and. use_ncep_ice ) then
               where ( sst_ncep <= (TFREEZE+tice_crit) )
                   sice = 1.
               elsewhere
                   sice = 0.
               endwhere
            else
               call read_record ('ice', Date2, Udate2, sice)
            endif
!--------------------------------------------------------------------
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

! <FUNCTION NAME="amip_interp_new_1d" INTERFACE="amip_interp_new">

!   <IN NAME="lon" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:)"> </IN>
!   <IN NAME="mask" TYPE="logical" DIM="(:,:)"> </IN>
!   <IN NAME="use_climo" TYPE="logical" DEFAULT="use_climo = .false."> </IN>
!   <IN NAME="use_annual" TYPE="logical" DEFAULT="use_annual = .false."> </IN>
!   <IN NAME="interp_method" TYPE="character(len=*), optional" DEFAULT="interp_method = conservative"></IN>
!   <OUT NAME="Interp" TYPE="amip_interp_type"> </OUT>

 function amip_interp_new_1d ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)

 real,    intent(in), dimension(:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional       :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (amip_interp_type) :: Interp

   if(.not.module_is_initialized) call amip_interp_init

   Interp % use_climo  = .false.
   if (present(use_climo)) Interp % use_climo  = use_climo
   Interp % use_annual = .false.
   if (present(use_annual)) Interp % use_annual  = use_annual

   if ( date_out_of_range == 'fail' .and. Interp%use_climo ) &
      call error_mesg ('amip_interp_new_1d', 'use_climo mismatch', FATAL)

   if ( date_out_of_range == 'fail' .and. Interp%use_annual ) &
      call error_mesg ('amip_interp_new_1d', 'use_annual(climo) mismatch', FATAL)

   Interp % Date1 = date_type( -99, -99, -99 )
   Interp % Date2 = date_type( -99, -99, -99 )

!-----------------------------------------------------------------------
!   ---- initialization of horizontal interpolation ----

    call horiz_interp_new ( Interp%Hintrp, lon_bnd, lat_bnd, &
                             lon, lat, interp_method= interp_method )

    allocate ( Interp % data1 (size(lon(:))-1,size(lat(:))-1), &
               Interp % data2 (size(lon(:))-1,size(lat(:))-1)  )

    Interp%I_am_initialized = .true.

   end function amip_interp_new_1d
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME="amip_interp_new_2d" INTERFACE="amip_interp_new">
!   <IN NAME="lon" TYPE="real" DIM="(:,:)"> </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)"> </IN>
!   <IN NAME="mask" TYPE="logical" DIM="(:,:)"> </IN>
!   <IN NAME="use_climo" TYPE="logical" DEFAULT="use_climo = .false."> </IN>
!   <IN NAME="use_annual" TYPE="logical" DEFAULT="use_annual = .false."> </IN>
!   <IN NAME="interp_method" TYPE="character(len=*), optional" DEFAULT="interp_method = conservative "></IN>
!   <OUT NAME="Interp" TYPE="amip_interp_type"> </OUT>

 function amip_interp_new_2d ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)

 real,    intent(in), dimension(:,:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (amip_interp_type) :: Interp

   if(.not.module_is_initialized) call amip_interp_init

   Interp % use_climo  = .false.
   if (present(use_climo)) Interp % use_climo  = use_climo
   Interp % use_annual = .false.
   if (present(use_annual)) Interp % use_annual  = use_annual

   if ( date_out_of_range == 'fail' .and. Interp%use_climo ) &
      call error_mesg ('amip_interp_new_2d', 'use_climo mismatch', FATAL)

   if ( date_out_of_range == 'fail' .and. Interp%use_annual ) &
      call error_mesg ('amip_interp_new_2d', 'use_annual(climo) mismatch', FATAL)

   Interp % Date1 = date_type( -99, -99, -99 )
   Interp % Date2 = date_type( -99, -99, -99 )

!-----------------------------------------------------------------------
!   ---- initialization of horizontal interpolation ----

   call horiz_interp_new ( Interp%Hintrp, lon_bnd, lat_bnd, &
                           lon, lat, interp_method = interp_method)

   allocate ( Interp % data1 (size(lon,1),size(lat,2)), &
              Interp % data2 (size(lon,1),size(lat,2))) 

   Interp%I_am_initialized = .true.

   end function amip_interp_new_2d
! </FUNCTION>

!#######################################################################

 subroutine amip_interp_init()

   integer :: unit,io,ierr

!-----------------------------------------------------------------------

    call horiz_interp_init

!   ---- read namelist ----

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, amip_interp_nml, iostat=io)
    ierr = check_nml_error(io,'amip_interp_nml')
#else
    if ( file_exist('input.nml')) then
       unit = open_namelist_file( )
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=amip_interp_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'amip_interp_nml')
       enddo
  10   call close_file (unit)
    endif
#endif

!  ----- write namelist/version info -----
    call write_version_number("AMIP_INTERP_MOD", version)

    unit = stdlog ( )
    if (mpp_pe() == 0) then
        write (unit,nml=amip_interp_nml)
    endif
    call close_file (unit)

    if ( .not. use_ncep_sst ) interp_oi_sst = .false.

!   ---- freezing point of sea water in deg K ---

    tice_crit_k = tice_crit
    if ( tice_crit_k < 200. ) tice_crit_k = tice_crit_k + TFREEZE
    ice_crit = nint((tice_crit_k-TFREEZE)*100.)

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
        Date_end = date_type( 2011, 8, 16 ) ! updated by JHC
! add by JHC
    else if (lowercase(trim(data_set)) == 'daily') then
        file_name_sst = 'INPUT/' // 'hurrell_sst.data'
        file_name_ice = 'INPUT/' // 'hurrell_ice.data'
        mobs = 360;  nobs = 180 
        call set_sst_grid_edges_oi
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AVHRR daily sst', NOTE)
        Date_end = date_type( 2011, 8, 16 )
! end add by JHC
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
!--- Added by SJL ---------------------------------------------- 
        if ( use_ncep_sst ) then
             mobs = i_sst;  nobs = j_sst
            if (.not. allocated (sst_ncep)) then
                allocate (sst_ncep(i_sst,j_sst))
                sst_ncep(:,:) = big_number
            endif
            if (.not. allocated (sst_anom)) then
                allocate (sst_anom(i_sst,j_sst))
                sst_anom(:,:) = big_number
            endif
        else
             mobs = 360;    nobs = 180
        endif
!--- Added by SJL ---------------------------------------------- 
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

    if (.not.file_exist(trim(file_name_sst)) .and. .not.file_exist(trim(file_name_sst)//'.nc')) then
      call error_mesg ('amip_interp_init', &
            'Neither '//trim(file_name_sst)//' or '//trim(file_name_sst)//'.nc exists', FATAL)
    endif
    if (.not.file_exist(trim(file_name_ice)) .and. .not.file_exist(trim(file_name_ice)//'.nc')) then
      call error_mesg ('amip_interp_init', &
            'Neither '//trim(file_name_ice)//' or '//trim(file_name_ice)//'.nc exists', FATAL)
    endif

    module_is_initialized = .true.

 end subroutine amip_interp_init

!#######################################################################

! <SUBROUTINE NAME="amip_interp_del">

!   <OVERVIEW>
!     Call this routine for all amip_interp_type variables created by amip_interp_new.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Call this routine for all amip_interp_type variables created by amip_interp_new.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call amip_interp_del (Interp)
!   </TEMPLATE>
!   <INOUT NAME="Interp" TYPE="amip_interp_type">
!     A defined data type variable initialized by amip_interp_new
!            and used when calling get_amip_sst and get_amip_ice.
!   </INOUT>

   subroutine amip_interp_del (Interp)
   type (amip_interp_type), intent(inout) :: Interp

     if(associated(Interp%data1)) deallocate(Interp%data1)
     if(associated(Interp%data2)) deallocate(Interp%data2)
     if(allocated(lon_bnd))       deallocate(lon_bnd)
     if(allocated(lat_bnd))       deallocate(lat_bnd)
     call horiz_interp_del ( Interp%Hintrp )

     Interp%I_am_initialized = .false.

   end subroutine amip_interp_del
!#######################################################################

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

! add by JHC
      if(allocated(lon_bnd))       deallocate(lon_bnd)
      if(allocated(lat_bnd))       deallocate(lat_bnd)
! end add by JHC
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
! add by JHC
   subroutine set_sst_grid_edges_daily(mobs_sst, nobs_sst)

   integer :: i, j, mobs_sst, nobs_sst
   real    :: hpie, dlon, dlat, wb, sb

      if(allocated(lon_bnd))       deallocate(lon_bnd)
      if(allocated(lat_bnd))       deallocate(lat_bnd)
      allocate ( lon_bnd(mobs_sst+1), lat_bnd(nobs_sst+1) )

! ---- compute grid edges (do only once) -----

      hpie = 0.5*pi

      dlon = 4.*hpie/float(mobs_sst);  wb = 0.0
          lon_bnd(1) = wb
      do i = 2, mobs_sst+1
          lon_bnd(i) = wb + dlon * float(i-1)
      enddo
          lon_bnd(mobs_sst+1) = lon_bnd(1) + 4.*hpie

      dlat = 2.*hpie/float(nobs_sst);  sb = -hpie
      lat_bnd(1) = sb;  lat_bnd(nobs_sst+1) = hpie
      do j = 2, nobs_sst
          lat_bnd(j) = sb + dlat * float(j-1)
      enddo

   end subroutine set_sst_grid_edges_daily
! end add by JHC
!#######################################################################


   subroutine a2a_bilinear(nx, ny, dat1, n1, n2, dat2)
   integer, intent(in):: nx, ny
   integer, intent(in):: n1, n2
   real, intent(in) :: dat1(nx,ny)
   real, intent(out):: dat2(n1,n2)      ! output interpolated data

! local:
  real:: lon1(nx), lat1(ny)
  real:: lon2(n1), lat2(n2)
  real:: dx1, dy1, dx2, dy2
  real:: xc, yc
  real:: a1, b1, c1, c2, c3, c4
  integer i1, i2, jc, i0, j0, it, jt
  integer i,j


!-----------------------------------------------------------
! * Interpolate from "FMS" 1x1 SST data grid to a finer grid
!                     lon: 0.5, 1.5, ..., 359.5
!                     lat: -89.5, -88.5, ... , 88.5, 89.5
!-----------------------------------------------------------

! INput Grid
  dx1 = 360./real(nx)
  dy1 = 180./real(ny)

  do i=1,nx
     lon1(i) = 0.5*dx1 + real(i-1)*dx1
  enddo
  do j=1,ny
     lat1(j) = -90. + 0.5*dy1 + real(j-1)*dy1
  enddo

! OutPut Grid:
  dx2 = 360./real(n1)
  dy2 = 180./real(n2)

  do i=1,n1
     lon2(i) = 0.5*dx2 + real(i-1)*dx2
  enddo
  do j=1,n2
     lat2(j) = -90. + 0.5*dy2 + real(j-1)*dy2
  enddo

  jt = 1
  do 5000 j=1,n2

     yc = lat2(j)
     if ( yc<lat1(1) ) then
            jc = 1
            b1 = 0.
     elseif ( yc>lat1(ny) ) then
            jc = ny-1
            b1 = 1.
     else
          do j0=jt,ny-1
          if ( yc>=lat1(j0) .and. yc<=lat1(j0+1) ) then
               jc = j0
               jt = j0
               b1 = (yc-lat1(jc)) / dy1
               go to 222
          endif
          enddo
     endif
222  continue

     it = 1
     do i=1,n1
        xc = lon2(i)
       if ( xc>lon1(nx) ) then
            i1 = nx;     i2 = 1
            a1 = (xc-lon1(nx)) / dx1
       elseif ( xc<lon1(1) ) then
            i1 = nx;     i2 = 1
            a1 = (xc+360.-lon1(nx)) / dx1
       else
            do i0=it,nx-1
            if ( xc>=lon1(i0) .and. xc<=lon1(i0+1) ) then
               i1 = i0;  i2 = i0+1
               it = i0
               a1 = (xc-lon1(i1)) / dx1
               go to 111
            endif
            enddo
       endif
111    continue

! Debug code:
       if ( a1<-0.001 .or. a1>1.001 .or.  b1<-0.001 .or. b1>1.001 ) then
            write(*,*) i,j,a1, b1
            call mpp_error(FATAL,'a2a bilinear interpolation')
       endif

       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1

! Bilinear interpolation:
       dat2(i,j) = c1*dat1(i1,jc) + c2*dat1(i2,jc) + c3*dat1(i2,jc+1) + c4*dat1(i1,jc+1)

     enddo   !i-loop

5000 continue   ! j-loop

   end subroutine a2a_bilinear

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
!   <ERROR MSG="have not called amip_interp_new" STATUS="FATAL">
!     Must call amip_interp_new before get_sst_grid_size.
!   </ERROR>

   subroutine get_sst_grid_size (nlon, nlat)

   integer, intent(out) :: nlon, nlat

      if ( .not.module_is_initialized ) call amip_interp_init

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
!   <ERROR MSG="have not called amip_interp_new" STATUS="FATAL">
!     Must call amip_interp_new before get_sst_grid_boundary.
!   </ERROR>
!   <ERROR MSG="invalid argument dimensions" STATUS="FATAL">
!     The size of the output argument arrays do not agree with
!      the size of the observed data. See the documentation for
!      interfaces get_sst_grid_size and get_sst_grid_boundary.
!   </ERROR>

   subroutine get_sst_grid_boundary (blon, blat, mask)

   real,    intent(out) :: blon(:), blat(:)
   logical, intent(out) :: mask(:,:)

      if ( .not.module_is_initialized ) call amip_interp_init

! ---- check size of argument(s) ----

      if (size(blon(:)) /= mobs+1 .or. size(blat(:)) /= nobs+1)   &
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
     type (date_type), intent(inout) :: Adate
     real,             intent(out) :: dat(mobs,nobs)
     real :: tmp_dat(360,180)

     real   (R4_KIND) :: dat4(mobs,nobs)
     integer(I2_KIND) :: idat(mobs,nobs)
     integer :: nrecords, yr, mo, dy, ierr, k
     integer, dimension(:), allocatable :: ryr, rmo, rdy
     character(len=38)   :: mesg
     character(len=maxc) :: ncfilename, ncfieldname

    !---- set file and field name for NETCDF data sets ----

        ncfieldname = 'sst'
     if(type(1:3) == 'sst') then
        ncfilename = trim(file_name_sst)//'.nc'
     else if(type(1:3) == 'ice') then
        ncfilename = trim(file_name_ice)//'.nc'
        if (lowercase(trim(data_set)) == 'amip2' .or. &
            lowercase(trim(data_set)) == 'hurrell' .or. & 
            lowercase(trim(data_set)) == 'daily') ncfieldname = 'ice' ! modified by JHC
     endif

    !---- make sure IEEE format file is open ----

     if ( (.NOT. file_exist(ncfilename))  ) then

       ! rewind condition (if unit is open)
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

     endif

     dy = 0 ! only processing monthly data

     if (verbose > 2 .and. mpp_pe() == 0)  &
          print *, 'looking for date = ', Date

    !---- check dates in NETCDF file -----

     ! This code can handle amip1, reynolds, or reyoi type SST data files in netCDF format
     if (file_exist(ncfilename)) then
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('amip_interp_mod', &
             'Reading NetCDF formatted input data file: '//trim(ncfilename), NOTE)
        call read_data (ncfilename, 'nrecords', nrecords, no_domain=.true.)
        if (nrecords < 1) call mpp_error('amip_interp_mod', &
                           'Invalid number of SST records in SST datafile: '//trim(ncfilename), FATAL)
        allocate(ryr(nrecords), rmo(nrecords), rdy(nrecords))
        call read_data(ncfilename, 'yr', ryr, no_domain=.true.)
        call read_data(ncfilename, 'mo', rmo, no_domain=.true.)
        call read_data(ncfilename, 'dy', rdy, no_domain=.true.)
        ierr = 1
        do k = 1, nrecords
          yr = ryr(k);  mo = rmo(k)
          Adate = date_type( yr, mo, 0)
          Curr_date = Adate
          if (verbose > 2 .and. mpp_pe() == 0)  &
                print *, '....... checking   ', Adate
          if (Date == Adate) ierr = 0 
          if (yr == 0 .and. mo == Date%month) ierr = 0
          if (ierr == 0) exit
        enddo
        if (ierr .ne. 0) call mpp_error('amip_interp_mod', &
                         'Model time is out of range not in SST data: '//trim(ncfilename), FATAL)
        deallocate(ryr, rmo, rdy)
       !PRINT *, 'New SST data: ', k, yr, mo, dy, Date%year, Date%month, Date%day, ryr(1), rmo(1)

    !---- check dates in IEEE file -----

     else
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('amip_interp_mod', &
             'Reading native formatted input data file: '//trim(data_set), NOTE)
        k = 0
        do
           k = k + 1
           if (lowercase(trim(data_set)) == 'amip2' .or. lowercase(trim(data_set)) == 'hurrell') then
              read (unit, end=10)  yr, mo,     dat4
              dat=dat4
           else
              read (unit, end=10)  yr, mo, dy, idat
           endif
           !new     read (unit, end=10)  yr, mo, dy
           Adate = date_type( yr, mo, dy )
           Curr_date = Adate
           if (verbose > 2 .and. mpp_pe() == 0)  &
                print *, '....... checking   ', Adate
           
           !     --- found date ---
           if (Date == Adate)                    exit
           if (Date%month == mo .and. Date%day == dy .and. Date%year == yr ) exit
           !     --- otherwise use monthly climo ---
           if (yr == 0 .and. Date % month == mo) exit
           
           !     --- skip this data record ---
           !new  if (lowercase(trim(data_set)) /= 'amip2') read (unit)
        enddo
        
        !     --- read data ---
        !new  if (lowercase(trim(data_set)) /= 'amip2') read (unit) idat
        
        !   --- check if climo used when not wanted ---
        
     endif ! if(file_exist(ncfilename))
     
   !---- check if climatological data should be used ----

     if (yr == 0 .or. mo == 0) then
        ierr = 0
        if (date_out_of_range == 'fail' )               ierr = 1
        if (date_out_of_range == 'initclimo' .and.  &
             Date > Date_end )   ierr = 1
        if (ierr /= 0) call error_mesg &
             ('read_record in amip_interp_mod', &
             'climo data read when NO climo data requested', FATAL)
     endif

   !---- read NETCDF data ----

     if (file_exist(ncfilename)) then
         if ( interp_oi_sst ) then
              call read_data(ncfilename, ncfieldname, tmp_dat, timelevel=k, no_domain=.true.)
!     interpolate tmp_dat(360, 180) ---> dat(mobs,nobs) (to enable SST anom computation)
              if ( mobs/=360 .or. nobs/=180 ) then
                   call a2a_bilinear(360, 180, tmp_dat, mobs, nobs, dat)
              else
                   dat(:,:) = tmp_dat(:,:)
              endif
         else
              call read_data(ncfilename, ncfieldname, dat, timelevel=k, no_domain=.true.)
         endif
        idat =  nint(dat*100.) ! reconstruct packed data for reproducibity
     endif

   !---- unpacking of data ----

     if (type(1:3) == 'ice') then
        !---- create fractional [0,1] ice mask
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
        !---- unpack sst ----
        if (lowercase(trim(data_set)) /= 'amip2' .and. lowercase(trim(data_set)) /= 'hurrell') then
               dat = real(idat)*0.01 + TFREEZE
        endif
     endif


     return

10   write (mesg, 20) unit
     call error_mesg ('read_record in amip_interp_mod', mesg, FATAL)

20   format ('end of file reading unit ',i2,' (sst data)')

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

subroutine amip_interp_type_eq(amip_interp_out, amip_interp_in)
    type(amip_interp_type), intent(inout) :: amip_interp_out
    type(amip_interp_type), intent(in)    :: amip_interp_in

    if(.not.amip_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'amip_interp_type_eq: amip_interp_type variable on right hand side is unassigned')
    endif

    amip_interp_out%Hintrp     =  amip_interp_in%Hintrp
    amip_interp_out%data1      => amip_interp_in%data1
    amip_interp_out%data2      => amip_interp_in%data2
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%use_climo  =  amip_interp_in%use_climo
    amip_interp_out%use_annual =  amip_interp_in%use_annual
    amip_interp_out%I_am_initialized = .true.

end subroutine amip_interp_type_eq

!#######################################################################

end module amip_interp_mod
! <INFO>

!   <FUTURE> 
!     Add AMIP 2 data set.
!
!     Other data sets (or extend current data sets).
!   </FUTURE>

! </INFO>
