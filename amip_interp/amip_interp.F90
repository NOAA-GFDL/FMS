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
!! 1. AMIP http://www.pcmdi.github.io/mips/amip from Jan 1979 to Jan 1989 (2 deg x 2 deg)
!! 2. Reynolds OI @ref amip_interp.rey_oi.txt from Nov 1981 to Jan 1999 (1 deg x 1 deg)
!! 3. Reynolds EOF podaac.jpl.nasa.gov/ from Jan 1950 to Dec 1998 (2 deg x 2 deg)
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
use      platform_mod, only: r4_kind, r8_kind, i2_kind, FMS_FILE_LEN
use mpp_mod,           only: input_nml_file
use fms2_io_mod,       only: FmsNetcdfFile_t, fms2_io_file_exists=>file_exists, open_file, close_file, &
                             get_dimension_size, fms2_io_read_data=>read_data
use netcdf,            only: NF90_MAX_NAME

implicit none
private

!-----------------------------------------------------------------------
!----------------- Public interfaces -----------------------------------

public amip_interp_init, get_amip_sst, get_amip_ice, amip_interp_new, &
     & amip_interp_del, amip_interp_type, assignment(=)

!-----------------------------------------------------------------------
!----------------- Public Data -----------------------------------
integer :: i_sst = 1200
integer :: j_sst = 600
real(r8_kind), parameter:: big_number = 1.E30_r8_kind
logical :: forecast_mode = .false.
real(r8_kind), allocatable, dimension(:,:) ::  sst_ncep, sst_anom

public i_sst, j_sst, sst_ncep, sst_anom, forecast_mode, use_ncep_sst

!-----------------------------------------------------------------------
!--------------------- private below here ------------------------------

!  ---- version number -----

! Include variable "version" to be written to log file.
#include<file_version.h>

! add by JHC
   real(r8_kind), allocatable, dimension(:,:) :: tempamip
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
  module procedure amip_interp_type_eq
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

!> Retrieve sea surface temperature data and interpolated grid
interface get_amip_sst
  module procedure get_amip_sst_r4, get_amip_sst_r8
end interface

!> AMIP interpolation for ice
interface get_amip_ice
  module procedure get_amip_ice_r4, get_amip_ice_r8
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
   module procedure amip_interp_new_1d_r4, amip_interp_new_1d_r8
   module procedure amip_interp_new_2d_r4, amip_interp_new_2d_r8
end interface

!----- public data type ------

!> @brief Contains information needed by the interpolation module (exchange_mod) and buffers
!! data (r4_kind flavor).
!> @ingroup amip_interp_mod
type amip_interp_type
   private
   type (horiz_interp_type)              :: Hintrp, Hintrp2 ! add by JHC
   real(r4_kind), dimension(:,:), allocatable :: data1_r4, data2_r4
   real(r8_kind), dimension(:,:), allocatable :: data1_r8, data2_r8
   type (date_type)                      :: Date1, Date2
   logical                               :: use_climo, use_annual
   logical                               :: I_am_initialized=.false.
end type amip_interp_type

!> @addtogroup amip_interp_mod
!> @{
!-----------------------------------------------------------------------
!  ---- resolution/grid variables ----

   integer :: mobs, nobs
   real(r8_kind), allocatable, dimension(:) :: lon_bnd, lat_bnd

!  ---- global unit & date ----

   integer :: iunit
   character(len=FMS_FILE_LEN) :: file_name_sst, file_name_ice
   type(FmsNetcdfFile_t), target :: fileobj_sst, fileobj_ice

   type (date_type) :: Curr_date = date_type( -99, -99, -99 )
   type (date_type) :: Date_end  = date_type( -99, -99, -99 )

   real(r8_kind)    :: tice_crit_k
   integer(i2_kind) ::  ice_crit

   logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------
!---- namelist ----

 character(len=24) :: data_set = 'amip1' !< use 'amip1', 'amip2', 'reynolds_eof'
                                         !! 'reynolds_oi', 'hurrell', or 'daily',
                                         !! when "use_daily=.T."
                                         ! add by JHC

 character(len=16) :: date_out_of_range = 'fail' !< use 'fail', 'initclimo', or 'climo'

 real(r8_kind)    :: tice_crit    = -1.80_r8_kind !<  in degC or degK
 integer :: verbose      = 0     !<  0 <= verbose <= 3

 logical :: use_zonal    = .false. !< parameters for prescribed zonal sst option
 real(r8_kind) :: teq  = 305._r8_kind !< parameters for prescribed zonal sst option
 real(r8_kind) :: tdif = 50._r8_kind !< parameters for prescribed zonal sst option
 real(r8_kind) :: tann = 20._r8_kind !< parameters for prescribed zonal sst option
 real(r8_kind) :: tlag = 0.875_r8_kind !< parameters for prescribed zonal sst option


 integer :: amip_date(3)=(/-1,-1,-1/) !< amip date for repeating single day (rsd) option

 real(r8_kind) :: sst_pert = 0._r8_kind !< global temperature perturbation used for sensitivity experiments

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

 !> initialize @ref amip_interp_mod for use
 subroutine amip_interp_init
   integer :: iunit,io,ierr

!-----------------------------------------------------------------------

    call horiz_interp_init

!   ---- read namelist ----

    read (input_nml_file, amip_interp_nml, iostat=io)
    ierr = check_nml_error(io,'amip_interp_nml')

!  ----- write namelist/version info -----
    call write_version_number("AMIP_INTERP_MOD", version)

    iunit = stdlog ( )
    if (mpp_pe() == 0) then
        write (iunit,nml=amip_interp_nml)
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
    if ( tice_crit_k < 200._r8_kind ) then
      tice_crit_k = tice_crit_k + TFREEZE
    endif
    ice_crit = nint((tice_crit_k-TFREEZE)*100._r8_kind, I2_KIND)

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
        tice_crit_k = 271.38_r8_kind
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 2 sst', NOTE)
        Date_end = date_type( 1996, 3, 0 )
    else if (lowercase(trim(data_set)) == 'hurrell') then
        file_name_sst = 'INPUT/' // 'hurrell_sst.data'
        file_name_ice = 'INPUT/' // 'hurrell_ice.data'
        mobs = 360;  nobs = 180
        call set_sst_grid_edges_oi
!       --- specfied min for hurrell ---
        tice_crit_k = 271.38_r8_kind
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

   subroutine set_sst_grid_edges_amip1
   integer :: i, j
   real(r8_kind) :: hpie, dlon, dlat, wb, sb

      allocate(lon_bnd(mobs+1))
      allocate(lat_bnd(nobs+1))

! ---- compute grid edges (do only once) -----

      hpie = pi / 2._r8_kind

      dlon = 4._r8_kind*hpie/real(mobs, r8_kind)
      wb = -0.5_r8_kind*dlon

      do i = 1, mobs+1
          lon_bnd(i) = wb + dlon*real(i-1, r8_kind)
      enddo
      lon_bnd(mobs+1) = lon_bnd(1) + 4._r8_kind*hpie

      dlat = 2._r8_kind*hpie/real(nobs-1, r8_kind)
      sb = -hpie + 0.5_r8_kind*dlat

      lat_bnd(1) = -hpie
      lat_bnd(nobs+1) = hpie
      do j = 2, nobs
          lat_bnd(j) = sb + dlat * real(j-2, r8_kind)
      enddo
   end subroutine set_sst_grid_edges_amip1

   subroutine set_sst_grid_edges_oi
   integer :: i, j
   real(r8_kind) :: hpie, dlon, dlat, wb, sb

! add by JHC
      if(allocated(lon_bnd)) deallocate(lon_bnd)
      if(allocated(lat_bnd)) deallocate(lat_bnd)
! end add by JHC

      allocate(lon_bnd(mobs+1))
      allocate(lat_bnd(nobs+1))

! ---- compute grid edges (do only once) -----

      hpie = pi / 2._r8_kind
      dlon = 4._r8_kind*hpie/real(mobs, r8_kind)
      wb = 0.0_r8_kind

      lon_bnd(1) = wb
      do i = 2, mobs+1
          lon_bnd(i) = wb + dlon * real(i-1, r8_kind)
      enddo
      lon_bnd(mobs+1) = lon_bnd(1) + 4._r8_kind*hpie

      dlat = 2._r8_kind*hpie/real(nobs, r8_kind)
      sb = -hpie

      lat_bnd(1) = sb
      lat_bnd(nobs+1) = hpie
      do j = 2, nobs
          lat_bnd(j) = sb + dlat * real(j-1, r8_kind)
      enddo
   end subroutine set_sst_grid_edges_oi

!> Frees data associated with a amip_interp_type variable. Should be used for any
!! variables initialized via @ref amip_interp_new.
!> @param[inout] Interp A defined data type variable initialized by amip_interp_new and used
!! when calling get_amip_sst and get_amip_ice.
   subroutine amip_interp_del (Interp)
     type (amip_interp_type), intent(inout) :: Interp

     if(allocated(Interp%data1_r4)) deallocate(Interp%data1_r4)
     if(allocated(Interp%data1_r8)) deallocate(Interp%data1_r8)
     if(allocated(Interp%data2_r4)) deallocate(Interp%data2_r4)
     if(allocated(Interp%data2_r8)) deallocate(Interp%data2_r8)

     if(allocated(lon_bnd))   deallocate(lon_bnd)
     if(allocated(lat_bnd))   deallocate(lat_bnd)

     call horiz_interp_del ( Interp%Hintrp )

     Interp%I_am_initialized = .false.
   end subroutine amip_interp_del

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

subroutine amip_interp_type_eq (amip_interp_out, amip_interp_in)
    type(amip_interp_type), intent(inout) :: amip_interp_out
    type(amip_interp_type), intent(in)    :: amip_interp_in

    if(.not.amip_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'amip_interp_type_eq: amip_interp_type variable on right hand side is unassigned')
    endif

    amip_interp_out%Hintrp     =  amip_interp_in%Hintrp
    amip_interp_out%Hintrp2    =  amip_interp_in%Hintrp2 !< missing assignment statement; added by GPP
    amip_interp_out%data1_r4   =  amip_interp_in%data1_r4
    amip_interp_out%data1_r8   =  amip_interp_in%data1_r8
    amip_interp_out%data2_r4   =  amip_interp_in%data2_r4
    amip_interp_out%data2_r8   =  amip_interp_in%data2_r8
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%use_climo  =  amip_interp_in%use_climo
    amip_interp_out%use_annual =  amip_interp_in%use_annual
    amip_interp_out%I_am_initialized = .true.
end subroutine amip_interp_type_eq

#include "amip_interp_r4.fh"
#include "amip_interp_r8.fh"

end module amip_interp_mod
!> @}
! <INFO>

!   <FUTURE>
!     Add AMIP 2 data set.
!
!     Other data sets (or extend current data sets).
!   </FUTURE>

! </INFO>
