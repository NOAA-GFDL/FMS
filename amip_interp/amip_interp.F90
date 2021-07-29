!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!
!> @defgroup amip_interp_mod amip_interp_mod
!> @ingroup amip_interp
!> @brief Provides observed sea surface temperature and ice mask data sets that have been
!! interpolated onto your model's grid.
!!
!> @author Bruce Wyman
!!
!> When using these routines three possible data sets are available:
!!
!! 1. AMIP @link http://www-pcmdi.llnl.gov/amip @endlink from Jan 1979 to Jan 1989 (2 deg x 2 deg)
!! 2. Reynolds OI @link amip_interp.rey_oi.txt @endlink from Nov 1981 to Jan 1999 (1 deg x 1 deg)
!! 3. Reynolds EOF @link ftp://podaac.jpl.nasa.gov/pub/sea_surface_temperature/reynolds/rsst/doc/rsst.html @endlink from Jan 1950 to Dec 1998 (2 deg x 2 deg)
!!
!! All original data are observed monthly means. This module
!! interpolates linearly in time between pairs of monthly means.
!! Horizontal interpolation is done using the horiz_interp module.
!!
!! When a requested date falls outside the range of dates available
!! a namelist option allows for use of the climatological monthly
!! mean values which are computed from all of the data in a particular
!! data set. \n
!! \n AMIP 1:\n
!!   from Jan 1979 to Jan 1989 (2 deg x 2 deg).\n\n
!! Reynolds OI:\n
!!   from Nov 1981 to Jan 1999 (1 deg x 1 deg)\n
!!   The analysis uses in situ and satellite SST's plus
!!   SST's simulated by sea-ice cover.\n\n
!! Reynold's EOF:\n
!!   from Jan 1950 to Dec 1998 (2 deg x 2 deg)\n
!!   NCEP Reynolds Historical Reconstructed Sea Surface Temperature
!!   The analysis uses both in-situ SSTs and satellite derived SSTs
!!   from the NOAA Advanced Very High Resolution Radiometer.
!!   In-situ data is used from 1950 to 1981, while both AVHRR derived
!!   satellite SSTs and in-situ data are used from 1981 to the
!!   end of 1998.
!!
!> @note The data set used by this module have been reformatted as 32-bit IEEE.
!!   The data values are packed into 16-bit integers.
!!
!!   The data sets are read from the following files:
!!
!!         amip1           INPUT/amip1_sst.data
!!         reynolds_io     INPUT/reyoi_sst.data
!!         reynolds_eof    INPUT/reynolds_sst.data


!> @file
!> File for amip_interp_mod

!> @addtogroup amip_interp_mod
!> @{
module amip_interp_mod

use  time_interp_mod, only: time_interp, fraction_of_year

use time_manager_mod, only: time_type, operator(+), operator(>), &
                             get_date, set_time, set_date

! add by JHC
use get_cal_time_mod, only: get_cal_time

! end add by JHC

use  horiz_interp_mod, only: horiz_interp_init, horiz_interp,  &
                             horiz_interp_new, horiz_interp_del, &
                             horiz_interp_type, assignment(=)

use           fms_mod, only: error_mesg, write_version_number,  &
                             NOTE, WARNING, FATAL, stdlog, check_nml_error, &
                             mpp_pe, lowercase, mpp_root_pe,    &
                             NOTE, mpp_error, fms_error_handler

use     constants_mod, only: TFREEZE, pi
use      platform_mod, only: R4_KIND, I2_KIND
use mpp_mod,           only: input_nml_file
use fms2_io_mod,       only: FmsNetcdfFile_t, fms2_io_file_exists=>file_exists, open_file, close_file, &
                             get_dimension_size, fms2_io_read_data=>read_data

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

!> @}

!> @brief Private data type for representing a calendar date
!> @ingroup amip_interp_mod
type date_type
   sequence
   integer :: year, month, day
end type

!> Assignment overload to allow native assignment between amip_interp_type variables.
!> @ingroup amip_interp_mod
interface assignment(=)
  module procedure  amip_interp_type_eq
end interface

!> Private logical equality overload for amip_interp_type
!> @ingroup amip_interp_mod
interface operator (==)
   module procedure date_equals
end interface

!> Private logical inequality overload for amip_interp_type
!> @ingroup amip_interp_mod
interface operator (/=)
   module procedure date_not_equals
end interface

!> Private logical greater than overload for amip_interp_type
!> @ingroup amip_interp_mod
interface operator (>)
   module procedure date_gt
end interface


!> Initializes data needed for the horizontal
!! interpolation between the sst data and model grid.
!!
!> The returned variable of type amip_interp_type is needed when
!! calling get_amip_sst and get_amip_ice.
!!
!> @param lon
!!     Longitude in radians of the model's grid box edges (1d lat/lon grid case)
!!     or at grid box mid-point (2d case for arbitrary grids).
!> @param lat
!!     Latitude in radians of the model's grid box edges (1d lat/lon grid case)
!!     or at grid box mid-point (2d case for arbitrary grids).
!> @param mask
!!     A mask for the model grid.
!> @param use_climo
!!     Flag the specifies that monthly mean climatological values will be used.
!> @param use_annual
!!     Flag the specifies that the annual mean climatological
!!              will be used.  If both use_annual = use_climo = true,
!!              then use_annual = true will be used.
!> @param interp_method
!!     specify the horiz_interp scheme. = "conservative" means conservative scheme,
!!     = "bilinear" means  bilinear interpolation.
!!
!> @return interp, a defined data type variable needed when calling get_amip_sst and get_amip_ice.
!!
!! \n Example usage:
!!
!!                 Interp = amip_interp_new ( lon, lat, mask, use_climo, use_annual, interp_method )
!!
!! This function may be called to initialize multiple variables
!! of type amip_interp_type.  However, there currently is no
!! call to release the storage used by this variable.
!!
!! The size of input augment mask must be a function of the size
!! of input augments lon and lat. The first and second dimensions
!! of mask must equal (size(lon,1)-1, size(lat,2)-1).
!!
!> @throws "FATAL: the value of the namelist parameter DATA_SET being used is not allowed"
!! Check the value of namelist variable DATA_SET.
!!
!> @throws "FATAL: requested input data set does not exist"
!! The data set requested is valid but the data does not exist in
!! the INPUT subdirectory. You may have requested amip2 data which
!! has not been officially set up.
!! See the section on DATA SETS to properly set the data up.
!!
!> @throws "FATAL: use_climo mismatch"
!! The namelist variable date_out_of_range = 'fail' and the amip_interp_new
!! argument use_climo = true.  This combination is not allowed.
!!
!> @throws "FATAL: use_annual(climo) mismatch"
!! The namelist variable date_out_of_range = 'fail' and the amip_interp_new
!! argument use_annual = true.  This combination is not allowed.
!!
!> @ingroup amip_interp_mod
interface amip_interp_new
   module procedure amip_interp_new_1d
   module procedure amip_interp_new_2d
end interface


!----- public data type ------

!> @brief Contains information needed by the interpolation module (exchange_mod) and buffers data.
!> @ingroup amip_interp_mod
type amip_interp_type
   private
   type (horiz_interp_type) :: Hintrp, Hintrp2 ! add by JHC
   real, pointer            ::    data1(:,:) =>NULL(), &
                                  data2(:,:) =>NULL()
   type (date_type)         ::    Date1,       Date2
   logical                  :: use_climo, use_annual
   logical                  :: I_am_initialized=.false.
end type

!> @addtogroup amip_interp_mod
!> @{
!-----------------------------------------------------------------------
!  ---- resolution/grid variables ----

   integer :: mobs, nobs
   real, allocatable :: lon_bnd(:), lat_bnd(:)

!  ---- global unit & date ----

   integer, parameter :: maxc = 128
   integer :: unit
   character(len=maxc) :: file_name_sst, file_name_ice
   type(FmsNetcdfFile_t), target :: fileobj_sst, fileobj_ice

   type (date_type) :: Curr_date = date_type( -99, -99, -99 )
   type (date_type) :: Date_end  = date_type( -99, -99, -99 )

   real             :: tice_crit_k
   integer(I2_KIND) ::  ice_crit

   logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------
!---- namelist ----

 character(len=24) :: data_set = 'amip1' !< use 'amip1', 'amip2', 'reynolds_eof'
                                         !! 'reynolds_oi', 'hurrell', or 'daily',
                                         !! when "use_daily=.T."
                                         ! add by JHC

 character(len=16) :: date_out_of_range = 'fail' !< use 'fail', 'initclimo', or 'climo'

 real    :: tice_crit    = -1.80 !<  in degC or degK
 integer :: verbose      = 0     !<  0 <= verbose <= 3

 logical :: use_zonal    = .false. !< parameters for prescribed zonal sst option
 real :: teq  = 305. !< parameters for prescribed zonal sst option
 real :: tdif = 50. !< parameters for prescribed zonal sst option
 real :: tann = 20. !< parameters for prescribed zonal sst option
 real :: tlag = 0.875 !< parameters for prescribed zonal sst option


 integer :: amip_date(3)=(/-1,-1,-1/) !< amip date for repeating single day (rsd) option

 real :: sst_pert = 0. !< global temperature perturbation used for sensitivity experiments

 character(len=6) :: sst_pert_type = 'fixed'  !< use 'random' or 'fixed'
 logical :: do_sst_pert = .false.
 logical :: use_daily = .false. !< if '.true.', give 'data_set = 'daily''

 logical :: use_ncep_sst = .false. !< SJL: During nudging:   use_ncep_sst = .T.;  no_anom_sst = .T.
                                   !!      during forecast:  use_ncep_sst = .T.;  no_anom_sst = .F.
 logical ::  no_anom_sst = .true.  !< SJL: During nudging:   use_ncep_sst = .T.;  no_anom_sst = .T.
                                   !!      during forecast:  use_ncep_sst = .T.;  no_anom_sst = .F.
 logical :: use_ncep_ice = .false. !< For seasonal forecast: use_ncep_ice = .F.
 logical :: interp_oi_sst = .false. !< changed to false for regular runs
 logical :: use_mpp_io = .false. !< Set to .true. to use mpp_io, otherwise fms2io is used

!> @page amip_interp_nml @ref amip_interp_mod Namelist
!!
!> @var character(len=24) data_set
!! Name/type of SST data that will be used.
!!        Possible values (case-insensitive) are:
!!                          1) amip1
!!                          2) reynolds_eof
!!                          3) reynolds_oi
!!        See the @ref amip_interp_oi page for more information
!! @var character(len=16) date_out_of_range
!!     Controls the use of climatological monthly mean data when
!!     the requested date falls outside the range of the data set.<BR/>
!!     Possible values are:
!!     <PRE>
!!   fail      - program will fail if requested date is prior
!!               to or after the data set period.
!!   initclimo - program uses climatological requested data is
!!               prior to data set period and will fail if
!!               requested date is after data set period.
!!   climo     - program uses climatological data anytime.
!!    </PRE>
!! @var real tice_crit
!!     Freezing point of sea water in degC or degK. Defaults to -1.80
!! @var integer verbose
!!     Controls printed output, 0 <= verbose <= 3, default=0
!!     additional parameters for controlling zonal prescribed sst ----
!!     these parameters only have an effect when use_zonal=.true. ----
!! @var logical use_zonal
!!     Flag to selected zonal sst or data set. Default=.false.
!! @var real teq
!!     sst at the equator. Default=305
!! @var real tdif
!!     Equator to pole sst difference. Default=50
!! @var real tann
!!     Amplitude of annual cycle. Default=20
!! @var real tlag
!!     Offset for time of year (for annual cycle). Default=0.875
!! @var integer amip_date
!!     Single calendar date in integer "(year,month,day)" format
!!     that is used only if set with year>0, month>0, day>0.
!!     If used, model calendar date is replaced by this date,
!!     but model time of day is still used to determine ice/sst.
!!     Used for repeating-single-day (rsd) experiments.
!!     Default=/-1,-1,-1/
!! @var real sst_pert
!!     Temperature perturbation in degrees Kelvin added onto the SST.
!!                The perturbation is globally-uniform (even near sea-ice).
!!                It is only used when abs(sst_pert) > 1.e-4.  SST perturbation runs
!!                may be useful in accessing model sensitivities.
!!     Default=0.
 namelist /amip_interp_nml/ use_ncep_sst, no_anom_sst, use_ncep_ice,  tice_crit, &
                            interp_oi_sst, data_set, date_out_of_range,          &
                            use_zonal, teq, tdif, tann, tlag, amip_date,         &
                            ! add by JHC
                            sst_pert, sst_pert_type, do_sst_pert,                &
                            use_daily,                                           &
                            ! end add by JHC
                            verbose, i_sst, j_sst, forecast_mode,                &
                            use_mpp_io

!-----------------------------------------------------------------------

contains

! modified by JHC
!> Retrieve sea surface temperature data and interpolated grid
subroutine get_amip_sst (Time, Interp, sst, err_msg, lon_model, lat_model)

   type (time_type),         intent(in)    :: Time !< Time to interpolate
   type (amip_interp_type),  intent(inout) :: Interp !< Holds data for interpolation
   real,                     intent(out)   ::  sst(:,:) !< Sea surface temperature data
   character(len=*), optional, intent(out) :: err_msg !< Holds error message string if present

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
    type(FmsNetcdfFile_t) :: fileobj
    logical :: the_file_exists
! end add by JHC
    logical, parameter :: DEBUG = .false. !> switch for debugging output
    !> These are fms_io specific
    integer :: unit

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

!! DEBUG CODE
    if (DEBUG) then
          call get_date(Amip_Time,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
          if (mpp_pe() == 0) then
             write (*,200) 'JHC: use_daily = F, AMIP_Time: ',jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6)
             write (*,300) 'JHC: use_daily = F, interped SST: ', sst(1,1),sst(5,5),sst(10,10)
          endif
    endif


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

    the_file_exists = fms2_io_file_exists(ncfilename)

    if ( (.NOT. the_file_exists)  ) then
        call mpp_error ('amip_interp_mod', &
             'cannot find daily SST input data file: '//trim(ncfilename), NOTE)
    else
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('amip_interp_mod', &
             'Reading NetCDF formatted daily SST from: '//trim(ncfilename), NOTE)

            if(.not. open_file(fileobj, trim(ncfilename), 'read')) &
                call error_mesg ('get_amip_sst', 'Error in opening file '//trim(ncfilename), FATAL)

            call get_dimension_size(fileobj, 'TIME', nrecords)
            if (nrecords < 1) call mpp_error('amip_interp_mod', &
                           'Invalid number of SST records in daily SST data file: '//trim(ncfilename), FATAL)
            allocate(timeval(nrecords), ryr(nrecords), rmo(nrecords), rdy(nrecords))
            call fms2_io_read_data(fileobj, 'TIME', timeval)
!!! DEBUG CODE
        if(DEBUG) then
          if (mpp_pe() == 0) then
             print *, 'JHC: nrecords = ', nrecords
             print *, 'JHC: TIME = ', timeval
          endif
        endif

        ierr = 1
        do k = 1, nrecords

          Udate = get_cal_time (timeval(k), time_unit, 'julian')
          call get_date(Udate,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
          ryr(k) = jhctod(1); rmo(k) = jhctod(2); rdy(k) = jhctod(3)

          if ( yr == ryr(k) .and. mo == rmo(k) .and. dy == rdy (k) ) ierr = 0
          if (ierr==0) exit

        enddo

        if(DEBUG) then
          if (mpp_pe() == 0) then
            print *, 'JHC: k =', k
            print *, 'JHC: ryr(k) rmo(k) rdy(k)',ryr(k), rmo(k), rdy(k)
            print *, 'JHC:  yr     mo     dy   ',yr, mo, dy
          endif
        endif

        if (ierr .ne. 0) call mpp_error('amip_interp_mod', &
                         'Model time is out of range not in SST data: '//trim(ncfilename), FATAL)
    endif ! if(file_exist(ncfilename))


   !---- read NETCDF data ----
     if ( .not. allocated(tempamip) ) allocate (tempamip(mobs_sst,nobs_sst))

     if (the_file_exists) then
          call fms2_io_read_data(fileobj, 'SST', tempamip, unlim_dim_level=k)
          call close_file(fileobj)
          tempamip = tempamip + TFREEZE

!!! DEBUG CODE
          if(DEBUG) then
            if (mpp_pe() == 0) then
              print*, 'JHC: TFREEZE = ', TFREEZE
              print*, lbound(sst)
              print*, ubound(sst)
              print*, lbound(tempamip)
              print*, ubound(tempamip)
              write(*,300) 'JHC: tempamip : ', tempamip(100,100), tempamip(200,200), tempamip(300,300)
            endif
          endif

          call horiz_interp ( Interp%Hintrp2, tempamip, sst )
          call clip_data ('sst', sst)

     endif

    if(DEBUG) then
      if (mpp_pe() == 400) then
        write(*,300)'JHC: use_daily = T, daily SST: ', sst(1,1),sst(5,5),sst(10,10)
        print *,'JHC: use_daily = T, daily SST: ', sst
      endif
    endif

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

       if(DEBUG) then
         if (mpp_pe() == 0) then
             print*, 'mobs = ', mobs
             print*, 'nobs = ', nobs
             print*, lbound(sst)
             print*, ubound(sst)
          endif
       endif

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

!> AMIP interpolation for ice
subroutine get_amip_ice (Time, Interp, ice, err_msg)

   type (time_type),         intent(in)    :: Time !< Time to interpolate
   type (amip_interp_type),  intent(inout) :: Interp !< Holds data for interpolation
   real,                     intent(out)   :: ice(:,:) !< ice data
   character(len=*), optional, intent(out) :: err_msg !< Holds error message string if present

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

!-----------------------------------------------------------------------
!---------- time interpolation (between months) ------------------------
!-----------------------------------------------------------------------

   ice = Interp % data1 + fmonth * (Interp % data2 - Interp % data1)

endif

!-----------------------------------------------------------------------

 end subroutine get_amip_ice

!#######################################################################

 !> @return A newly created @ref amip_interp_type
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

 !> @return A newly created @ref amip_interp_type
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

!#######################################################################

 !> initialize @ref amip_interp_mod for use
 subroutine amip_interp_init()

   integer :: unit,io,ierr

!-----------------------------------------------------------------------

    call horiz_interp_init

!   ---- read namelist ----

    read (input_nml_file, amip_interp_nml, iostat=io)
    ierr = check_nml_error(io,'amip_interp_nml')

!  ----- write namelist/version info -----
    call write_version_number("AMIP_INTERP_MOD", version)

    unit = stdlog ( )
    if (mpp_pe() == 0) then
        write (unit,nml=amip_interp_nml)
    endif

    if (use_mpp_io) then
            !! USE_MPP_IO_WARNING
            call mpp_error ('amip_interp_mod', &
             'MPP_IO is no longer supported.  Please remove use_mpp_io from amip_interp_nml',&
              FATAL)
    endif
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
    file_name_sst = trim(file_name_sst)//'.nc'
    file_name_ice = trim(file_name_ice)//'.nc'

    if (.not. fms2_io_file_exists(trim(file_name_sst)) ) then
        call error_mesg ('amip_interp_init', &
             'file '//trim(file_name_sst)//' does not exist', FATAL)
    endif
    if (.not. fms2_io_file_exists(trim(file_name_ice)) ) then
        call error_mesg ('amip_interp_init', &
             'file '//trim(file_name_ice)//' does not exist', FATAL)
    endif

    if (.not. open_file(fileobj_sst, trim(file_name_sst), 'read')) &
        call error_mesg ('amip_interp_init', 'Error in opening file '//trim(file_name_sst), FATAL)
    if (.not. open_file(fileobj_ice, trim(file_name_ice), 'read')) &
        call error_mesg ('amip_interp_init', 'Error in opening file '//trim(file_name_ice), FATAL)
    module_is_initialized = .true.

 end subroutine amip_interp_init

!#######################################################################

!> Frees data associated with a amip_interp_type variable. Should be used for any
!! variables initialized via @ref amip_interp_new.
!> @param[inout] Interp A defined data type variable initialized by amip_interp_new and used
!! when calling get_amip_sst and get_amip_ice.
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
   real, intent(out):: dat2(n1,n2)      !> output interpolated data

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

  dx1 = 360./real(nx) !> INput Grid
  dy1 = 180./real(ny) !> INput Grid

  do i=1,nx
     lon1(i) = 0.5*dx1 + real(i-1)*dx1
  enddo
  do j=1,ny
     lat1(j) = -90. + 0.5*dy1 + real(j-1)*dy1
  enddo

  dx2 = 360./real(n1) !> OutPut Grid:
  dy2 = 180./real(n2) !> OutPut Grid:

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

!> @brief Returns the size (i.e., number of longitude and latitude
!!         points) of the observed data grid.
!! @throws FATAL have not called amip_interp_new
!!     Must call amip_interp_new before get_sst_grid_size.
   subroutine get_sst_grid_size (nlon, nlat)

   integer, intent(out) :: nlon !> The number of longitude points (first dimension) in the
                                !! observed data grid.  For AMIP 1 nlon = 180, and the Reynolds nlon = 360.
   integer, intent(out) :: nlat !> The number of latitude points (second dimension) in the
                                !! observed data grid.  For AMIP 1 nlon = 91, and the Reynolds nlon = 180.

      if ( .not.module_is_initialized ) call amip_interp_init

      nlon = mobs;  nlat = nobs

   end subroutine get_sst_grid_size

!#######################################################################

!> @brief Returns the grid box boundaries of the observed data grid.
!!
!! @throws FATAL, have not called amip_interp_new
!!     Must call amip_interp_new before get_sst_grid_boundary.
!!
!! @throws FATAL, invalid argument dimensions
!!     The size of the output argument arrays do not agree with
!!     the size of the observed data. See the documentation for
!!     interfaces get_sst_grid_size and get_sst_grid_boundary.
   subroutine get_sst_grid_boundary (blon, blat, mask)

   real,    intent(out) :: blon(:) !> The grid box edges (in radians) for longitude points of the
                                   !! observed data grid. The size of this argument must be nlon+1.
   real,    intent(out) :: blat(:) !> The grid box edges (in radians) for latitude points of the
                                   !! observed data grid. The size of this argument must be nlat+1.
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

!#######################################################################

   subroutine read_record (type, Date, Adate, dat)

     character(len=*), intent(in)  :: type
     type (date_type), intent(in)  :: Date
     type (date_type), intent(inout) :: Adate
     real,             intent(out) :: dat(mobs,nobs)
     real :: tmp_dat(360,180)

     integer(I2_KIND) :: idat(mobs,nobs)
     integer :: nrecords, yr, mo, dy, ierr, k
     integer, dimension(:), allocatable :: ryr, rmo, rdy
     character(len=38)   :: mesg
     character(len=maxc) :: ncfilename, ncfieldname
     type(FmsNetcdfFile_t), pointer :: fileobj

    !---- set file and field name for NETCDF data sets ----

        ncfieldname = 'sst'
     if(type(1:3) == 'sst') then
        ncfilename = trim(file_name_sst)
        fileobj => fileobj_sst
     else if(type(1:3) == 'ice') then
        ncfilename = trim(file_name_ice)
        fileobj => fileobj_ice
        if (lowercase(trim(data_set)) == 'amip2' .or. &
            lowercase(trim(data_set)) == 'hurrell' .or. &
            lowercase(trim(data_set)) == 'daily') ncfieldname = 'ice' ! modified by JHC
     endif

     dy = 0 ! only processing monthly data

     if (verbose > 2 .and. mpp_pe() == 0)  &
          print *, 'looking for date = ', Date

     ! This code can handle amip1, reynolds, or reyoi type SST data files in netCDF format
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('amip_interp_mod', &
          'Reading NetCDF formatted input data file: '//trim(ncfilename), NOTE)

        call fms2_io_read_data (fileobj, 'nrecords', nrecords)
        if (nrecords < 1) call mpp_error('amip_interp_mod', &
                           'Invalid number of SST records in SST datafile: '//trim(ncfilename), FATAL)
        allocate(ryr(nrecords), rmo(nrecords), rdy(nrecords))
        call fms2_io_read_data(fileobj, 'yr', ryr)
        call fms2_io_read_data(fileobj, 'mo', rmo)
        call fms2_io_read_data(fileobj, 'dy', rdy)

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

     if ( interp_oi_sst ) then
          call fms2_io_read_data(fileobj, ncfieldname, tmp_dat, unlim_dim_level=k)
!     interpolate tmp_dat(360, 180) ---> dat(mobs,nobs) (to enable SST anom computation)
          if ( mobs/=360 .or. nobs/=180 ) then
               call a2a_bilinear(360, 180, tmp_dat, mobs, nobs, dat)
          else
               dat(:,:) = tmp_dat(:,:)
          endif
     else
          call fms2_io_read_data(fileobj, ncfieldname, dat, unlim_dim_level=k)
     endif
     idat =  nint(dat) ! reconstruct packed data for reproducibility

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

!> @return logical answer
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

!> @return logical answer
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

!> @return logical answer
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
!> @}
! <INFO>

!   <FUTURE>
!     Add AMIP 2 data set.
!
!     Other data sets (or extend current data sets).
!   </FUTURE>

! </INFO>
