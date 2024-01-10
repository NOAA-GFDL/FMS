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
!> @defgroup interpolator_mod interpolator_mod
!> @ingroup interpolator
!> @brief A module to interpolate climatology data to model the grid.
!> @author William Cooke <William.Cooke@noaa.gov>
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_init receives various data as input in order
!!        to initialize interpolating.
!!
!! @param [inout] <clim_type> An interpolate type containing the necessary file and field
!!                                     data to be passed to the interpolator routine
!! @param [in] <file_name> Climatology filename
!! @param [in] <lonb_mod> The corners of the model grid-box longitudes
!! @param [in] <latb_mod> The corners of the model grid_box latitudes
!! @param [in] <data_names> OPTIONAL: A list of the names of components within the climatology
!!                                               file which you wish to read
!! @param [in] <data_out_of_bounds> A list of the flags that are to be used in determining
!!                                             what to do if the pressure levels in the model go out of
!!                                      bounds from those of the climatology
!! @param [in] <vert_interp> OPTIONAL: Flag to determine type of vertical interpolation
!! @param [out] <clim_units> OPTIONAL: A list of the units for the components listed in data_names
!! @param [out] <single_year_file> OPTIONAL: No description
subroutine INTERPOLATOR_INIT_( clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units, single_year_file)
type(interpolate_type), intent(inout)   :: clim_type
character(len=*), intent(in)            :: file_name
real(FMS_INTP_KIND_), intent(in)        :: lonb_mod(:,:), latb_mod(:,:)
character(len=*), intent(in) , optional :: data_names(:)
!++lwh
integer         , intent(in)            :: data_out_of_bounds(:)
integer         , intent(in), optional  :: vert_interp(:)
!--lwh
character(len=*), intent(out), optional :: clim_units(:)
logical,          intent(out), optional :: single_year_file
!
! INTENT IN
!  file_name  :: Climatology filename
!  lonb_mod   :: The corners of the model grid-box longitudes.
!  latb_mod   :: The corners of the model grid_box latitudes.
!  data_names :: A list of the names of components within the climatology file which you wish to read.
!  data_out_of_bounds :: A list of the flags that are to be used in determining what to do
!  if the pressure levels in the model
!                        go out of bounds from those of the climatology.
!  vert_interp:: Flag to determine type of vertical interpolation
!
! INTENT OUT
!  clim_type  :: An interpolate type containing the necessary file and field data to be passed
!  to the interpolator routine.
!  clim_units :: A list of the units for the components listed in data_names.
!
integer :: io, ierr

if (.not. module_is_initialized) then
  call fms_init
  call diag_manager_init
  call horiz_interp_init

!--------------------------------------------------------------------
! namelist input
!--------------------------------------------------------------------

  read (input_nml_file, nml=interpolator_nml, iostat=io)
  ierr = check_nml_error(io,'interpolator_nml')

  ! retain_cm3_bug is no longer supported.
  if (retain_cm3_bug) then
    call mpp_error(FATAL, "interpolator_init: You have overridden the default value of " // &
       "retain_cm3_bug and set it to .true. in interpolator_nml. This was a temporary workaround " // &
       "that is no longer supported. Please remove this namelist variable.")
  endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
call write_version_number("INTERPOLATOR_MOD", version)

if (mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=interpolator_nml)

module_is_initialized = .true.

endif !> if (module_is_initilized)
if (use_mpp_io) then
   call mpp_error(FATAL, "Interpolator::nml=interpolator_nml " //&
           'MPP_IO is no longer supported.  Please remove from use_mpp_io from interpolator_nml')
endif

clim_type%FMS_INTP_TYPE_%is_allocated = .true.

call fms2io_interpolator_init(clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units, single_year_file)

end subroutine INTERPOLATOR_INIT_

subroutine FMS2IO_INTERPOLATOR_INIT_(clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units, single_year_file)

type(interpolate_type), intent(inout)   :: clim_type
character(len=*), intent(in)            :: file_name
real(FMS_INTP_KIND_) , intent(in)       :: lonb_mod(:,:), latb_mod(:,:)
character(len=*), intent(in) , optional :: data_names(:)
!++lwh
integer         , intent(in)            :: data_out_of_bounds(:)
integer         , intent(in), optional  :: vert_interp(:)
!--lwh
character(len=*), intent(out), optional :: clim_units(:)
logical,          intent(out), optional :: single_year_file

character(len=64)            :: src_file
!++lwh
real(FMS_INTP_KIND_)         :: dlat, dlon
!--lwh
type(time_type)              :: base_time
integer                      :: fileday, filemon, fileyr, filehr, filemin,filesec, m,m1
character(len= 20)           :: fileunits
character(len=128)           :: var_dimname(6)
character(len=128), allocatable :: var_names(:)
integer,            allocatable :: var_ndims(:)
integer   :: j, i
logical :: non_monthly
character(len=24) :: file_calendar
character(len=256) :: error_mesg
integer :: model_calendar
integer :: yr, mo, dy, hr, mn, sc
integer :: n
type(time_type) :: Julian_time, Noleap_time
real(r8_kind), allocatable :: time_in(:)
real(FMS_INTP_KIND_), allocatable, save :: agrid_mod(:,:,:)
integer :: nx, ny
type(FmsNetcdfFile_t) :: fileobj

integer, parameter :: lkind=FMS_INTP_KIND_
real(FMS_INTP_KIND_), parameter :: lPI=real(PI,FMS_INTP_KIND_)

!> variables used to set time
logical :: yearly !< flags to indicate if time data is in units of months or years
integer :: num_years !< number of years
integer :: base_year, base_month, base_day, base_hour, base_minute, base_second
integer :: nn !< counter
logical :: noleap_file_calendar !< is the file calendar noleap or julian
real(r8_kind) :: num_days, frac_year !< variables for yearly=.true.
type(time_type) :: n_time !< temporary time

clim_type%separate_time_vary_calc = .false.

num_fields = 0

!--------------------------------------------------------------------
! open source file containing fields to be interpolated
!--------------------------------------------------------------------
src_file = 'INPUT/'//trim(file_name)

if(fms2_io_file_exist(trim(src_file))) then
   if(.not. open_file(clim_type%fileobj, trim(src_file), 'read')) &
        call mpp_error(FATAL, 'Interpolator_init: Error in opening file '//trim(src_file))
   fileobj = clim_type%fileobj
else
!Climatology file doesn't exist, so exit
   call mpp_error(FATAL,'Interpolator_init : Data file '//trim(src_file)//' does not exist')
endif

!Find the number of variables (nvar) in this file
clim_type%file_name = trim(file_name)
nvar = get_num_variables(fileobj)
num_fields = nvar
if(present(data_names)) then
   num_fields= size(data_names(:))
else
   allocate(var_names(nvar), var_ndims(nvar))
   call get_variable_names(fileobj, var_names)
   !--- loop through all the vars to exclude scalar or 1-D array
   do i=1,nvar
      var_ndims(i) = get_variable_num_dimensions(fileobj, var_names(i))
   enddo
   num_fields = count(var_ndims>1)
endif

! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------

nlon=0 ! Number of longitudes (center-points) in the climatology.
nlat=0 ! Number of latitudes (center-points) in the climatology.
nlev=0 ! Number of levels (center-points) in the climatology.
nlatb=0 ! Number of longitudes (boundaries) in the climatology.
nlonb=0 ! Number of latitudes (boundaries) in the climatology.
nlevh=0 ! Number of levels (boundaries) in the climatology.

clim_type%level_type = 0 ! Default value

!++lwh
! -------------------------------------------------------------------
! For 2-D fields, set a default value of nlev=nlevh=1
! -------------------------------------------------------------------
nlev = 1
nlevh = 1
!--lwh
        clim_type%vertical_indices = 0  ! initial value

!--- get clim_type%lat
if(dimension_exists(fileobj, "lat")) then
   call get_dimension_size(fileobj, "lat", nlat)
else
   call mpp_error(FATAL,'Interpolator_init : dimension lat does not exist in file '//trim(src_file) )
endif
allocate(clim_type%FMS_INTP_TYPE_%lat(nlat))
call get_axis_latlon_data(fileobj, 'lat', clim_type%FMS_INTP_TYPE_%lat)

!--- get clim_type%lon
if(dimension_exists(fileobj, "lon")) then
   call get_dimension_size(fileobj, "lon", nlon)
else
   call mpp_error(FATAL,'Interpolator_init : dimension lon does not exist in file '//trim(src_file) )
endif
allocate(clim_type%FMS_INTP_TYPE_%lon(nlon))
call get_axis_latlon_data(fileobj, 'lon', clim_type%FMS_INTP_TYPE_%lon)

!--- get clim_type%latb
if(dimension_exists(fileobj, "latb")) then
   call get_dimension_size(fileobj, "latb", nlatb)
   allocate(clim_type%FMS_INTP_TYPE_%latb(nlatb))
   call get_axis_latlon_data(fileobj, 'latb', clim_type%FMS_INTP_TYPE_%latb)
else
   if(nlat == 1) call mpp_error(FATAL,'Interpolator_init : nlat is 1')
   ! In the case where only the grid midpoints of the latitudes are defined we force the
   ! definition of the boundaries to be half-way between the midpoints.
   allocate(clim_type%FMS_INTP_TYPE_%latb(nlat+1))
   dlat = (clim_type%FMS_INTP_TYPE_%lat(2)-clim_type%FMS_INTP_TYPE_%lat(1)) * 0.5_lkind
!  clim_type%latb(1) = min( 90., max(-90., clim_type%lat(1) - dlat) )
   clim_type%FMS_INTP_TYPE_%latb(1) = min( lPI/2._lkind, max(-lPI/2._lkind, clim_type%FMS_INTP_TYPE_%lat(1) - dlat) )
   clim_type%FMS_INTP_TYPE_%latb(2:nlat) = ( clim_type%FMS_INTP_TYPE_%lat(1:nlat-1) + &
                                             clim_type%FMS_INTP_TYPE_%lat(2:nlat) ) * 0.5_lkind
   dlat = ( clim_type%FMS_INTP_TYPE_%lat(nlat) - clim_type%FMS_INTP_TYPE_%lat(nlat-1) ) * 0.5_lkind
!  clim_type%latb(nlat+1) = min( 90., max(-90., clim_type%lat(nlat) + dlat) )
   clim_type%FMS_INTP_TYPE_%latb(nlat+1) = min( lPI/2._lkind, max(-lPI/2._lkind, &
                                                clim_type%FMS_INTP_TYPE_%lat(nlat) + dlat) )
endif

!--- get clim_type%lonb
if(dimension_exists(fileobj, "lonb")) then
   call get_dimension_size(fileobj, "lonb", nlonb)
   allocate(clim_type%FMS_INTP_TYPE_%lonb(nlonb))
   call get_axis_latlon_data(fileobj, 'lonb', clim_type%FMS_INTP_TYPE_%lonb)
else
! In the case where only the midpoints of the longitudes are defined we force the definition
! of the boundaries to be half-way between the midpoints.
  if (size(clim_type%FMS_INTP_TYPE_%lon(:)) /= 1) then
    allocate(clim_type%FMS_INTP_TYPE_%lonb(size(clim_type%FMS_INTP_TYPE_%lon(:))+1))
    dlon = ( clim_type%FMS_INTP_TYPE_%lon(2)-clim_type%FMS_INTP_TYPE_%lon(1) )/2.0_lkind
    clim_type%FMS_INTP_TYPE_%lonb(1)  = clim_type%FMS_INTP_TYPE_%lon(1) - dlon
    clim_type%FMS_INTP_TYPE_%lonb(2:) = clim_type%FMS_INTP_TYPE_%lon(1:) + dlon
  else
!! this is the case for zonal mean data, lon = 1, lonb not present
!! in file.
    allocate(clim_type%FMS_INTP_TYPE_%lonb(2))
    clim_type%FMS_INTP_TYPE_%lonb(1) = -360.0_lkind*real(dtr,FMS_INTP_KIND_)
    clim_type%FMS_INTP_TYPE_%lonb(2) =  360.0_lkind*real(dtr,FMS_INTP_KIND_)
    clim_type%FMS_INTP_TYPE_%lon(1) = 0.0_lkind
  endif
endif

!--- get clim_type%levs
clim_type%level_type = 0
clim_type%vertical_indices = 0
if(dimension_exists(fileobj, "pfull")) then
   call get_dimension_size(fileobj, "pfull", nlev)
   allocate(clim_type%FMS_INTP_TYPE_%levs(nlev))
   call get_axis_level_data(fileobj, 'pfull', clim_type%FMS_INTP_TYPE_%levs, &
                            clim_type%level_type, clim_type%vertical_indices)
else if(dimension_exists(fileobj, "sigma_full")) then
   call get_dimension_size(fileobj, "sigma_full", nlev)
   allocate(clim_type%FMS_INTP_TYPE_%levs(nlev))
   call fms2_io_read_data(fileobj, "sigma_full", clim_type%FMS_INTP_TYPE_%levs)
   clim_type%level_type = SIGMA
else
   ! -------------------------------------------------------------------
   ! For 2-D fields, allocate levs and halflevs here
   !  code is still needed for case when only halflevs are in data file.
   ! -------------------------------------------------------------------
   nlev = 1
   allocate( clim_type%FMS_INTP_TYPE_%levs(nlev) )
   clim_type%FMS_INTP_TYPE_%levs = 0.0_lkind
endif

!--- get clim_type%halflevs
if(dimension_exists(fileobj, "phalf")) then
   call get_dimension_size(fileobj, "phalf", nlevh)
   allocate(clim_type%FMS_INTP_TYPE_%halflevs(nlevh))
   call get_axis_level_data(fileobj, 'phalf', clim_type%FMS_INTP_TYPE_%halflevs, &
                            clim_type%level_type, clim_type%vertical_indices)
else if(dimension_exists(fileobj, "sigma_half")) then
   call get_dimension_size(fileobj, "sigma_half", nlevh)
   allocate(clim_type%FMS_INTP_TYPE_%halflevs(nlevh))
   call fms2_io_read_data(fileobj, "sigma_half", clim_type%FMS_INTP_TYPE_%halflevs)
   clim_type%level_type = SIGMA
else
   allocate( clim_type%FMS_INTP_TYPE_%halflevs(nlev+1) )
   clim_type%FMS_INTP_TYPE_%halflevs(1) = 0.0_lkind
   if (clim_type%level_type == PRESSURE) then
      clim_type%FMS_INTP_TYPE_%halflevs(nlev+1) = 1013.25_lkind* 100.0_lkind   ! MKS
   else if (clim_type%level_type == SIGMA   ) then
      clim_type%FMS_INTP_TYPE_%halflevs(nlev+1) = 1.0_lkind
   endif
      do n=2,nlev
         clim_type%FMS_INTP_TYPE_%halflevs(n) = 0.5_lkind*(clim_type%FMS_INTP_TYPE_%levs(n) + &
                                                           clim_type%FMS_INTP_TYPE_%levs(n-1))
      end do
endif

!get time informaiton
if(dimension_exists(fileobj, "time")) then
   call get_dimension_size(fileobj, "time", ntime)

   call get_variable_units(fileobj, "time", units)
   call get_time_calendar(fileobj, "time", file_calendar)
   model_calendar = get_calendar_type()
   fileday = 0
   filemon = 0
   fileyr  = 0
   filehr  = 0
   filemin = 0
   filesec = 0
   yearly  = .false.
   select case(units(:3))
     case('day')
       fileunits = units(12:) !Assuming "days since YYYY-MM-DD HH:MM:SS"
       if ( len_trim(fileunits) < 19 ) then
         write(error_mesg, '(A49,A,A49,A)' ) &
           'Interpolator_init : Incorrect time units in file ', &
           trim(file_name), '. Expecting days since YYYY-MM-DD HH:MM:SS, found', &
           trim(units)
         call mpp_error(FATAL,error_mesg)
       endif
       read(fileunits(1:4)  , *)  fileyr
       read(fileunits(6:7)  , *)  filemon
       read(fileunits(9:10) , *)  fileday
       read(fileunits(12:13), *)  filehr
       read(fileunits(15:16), *)  filemin
       read(fileunits(18:19), *)  filesec
     case('mon')
       fileunits = units(14:) !Assuming "months since YYYY-MM-DD HH:MM:SS"
       if ( len_trim(fileunits) < 19 ) then
         write(error_mesg, '(A49,A,A51,A)' ) &
           'Interpolator_init : Incorrect time units in file ', &
           trim(file_name), '. Expecting months since YYYY-MM-DD HH:MM:SS, found', &
           trim(units)
         call mpp_error(FATAL,error_mesg)
       endif
       read(fileunits(1:4)  , *)  fileyr
       read(fileunits(6:7)  , *)  filemon
       read(fileunits(9:10) , *)  fileday
       read(fileunits(12:13), *)  filehr
       read(fileunits(15:16), *)  filemin
       read(fileunits(18:19), *)  filesec
    case( 'yea')
       fileunits = units(13:) !Assuming "years since YYYY-MM-DD HH:MM:SS"
       if ( len_trim(fileunits) < 19 ) then
         write(error_mesg, '(A49,A,A51,A)' ) &
           'Interpolator_init : Incorrect time units in file ', &
           trim(file_name), '. Expecting years since YYYY-MM-DD HH:MM:SS, found', &
           trim(units)
         call mpp_error(FATAL,error_mesg)
       endif
       read(fileunits(1:4)  , *)  fileyr
       read(fileunits(6:7)  , *)  filemon
       read(fileunits(9:10) , *)  fileday
       read(fileunits(12:13), *)  filehr
       read(fileunits(15:16), *)  filemin
       read(fileunits(18:19), *)  filesec
       yearly = .true.
     case default
       call mpp_error(FATAL,'Interpolator_init : Time units not recognized in file '//file_name)
   end select

   clim_type%climatological_year = (fileyr == 0)

   if (.not. clim_type%climatological_year) then

      !----------------------------------------------------------------------
      !    if file date has a non-zero year in the base time, determine that
      !    base_time based on the netcdf info.
      !----------------------------------------------------------------------
      noleap_file_calendar=.false.
      if ( (model_calendar == JULIAN .and. trim(adjustl(lowercase(file_calendar))) == 'julian')) then
         call mpp_error (NOTE, 'interpolator_mod: Model and file&
              & calendars are the same for file ' //   &
              & trim(file_name) // '; no calendar conversion  &
              &needed')
         base_time = set_date (fileyr, filemon, fileday, filehr,filemin,filesec)
         noleap_file_calendar=.false.
      else if((model_calendar == NOLEAP .and. trim(adjustl(lowercase(file_calendar))) == 'noleap') )  then
         call mpp_error (NOTE, 'interpolator_mod: Model and file&
              & calendars are the same for file ' //   &
              & trim(file_name) // '; no calendar conversion  &
              &needed')
         base_time = set_date (fileyr, filemon, fileday, filehr,filemin,filesec)
         noleap_file_calendar=.true.
      else if ( (model_calendar == JULIAN .and. trim(adjustl(lowercase(file_calendar))) == 'noleap')) then
         call mpp_error (NOTE, 'interpolator_mod: Using julian &
                          &model calendar and noleap file calendar&
                          & for file ' // trim(file_name) //   &
                          &'; calendar conversion needed')
        base_time = set_date_no_leap (fileyr, filemon, fileday,filehr, filemin, filesec)
        noleap_file_calendar=.true.
      else if ( (model_calendar == NOLEAP .and. trim(adjustl(lowercase(file_calendar))) == 'julian')) then
        call mpp_error (NOTE, 'interpolator_mod: Using noleap &
                          &model calendar and julian file calendar&
                          & for file ' // trim(file_name) //  &
                          &'; calendar conversion needed')
        base_time = set_date_julian (fileyr, filemon, fileday,filehr, filemin, filesec)
        noleap_file_calendar=.false.
      else
        call mpp_error (FATAL , 'interpolator_mod: Model and file&
             & calendars ( ' // trim(file_calendar) // ' ) differ  &
             &for file ' // trim(file_name) // ';  this calendar  &
             &conversion not currently available')
      endif

   else

       !! if the year is specified as '0000', then the file is intended to
       !! apply to all years -- the time variables within the file refer to
       !! the displacement from the start of each year to the time of the
       !! associated data. Time interpolation is to be done with interface
       !! time_interp_list, with the optional argument modtime=YEAR. base_time
       !! is set to an arbitrary value here; it's only use will be as a
       !! timestamp for optionally generated diagnostics.
       base_time = get_base_time ()
   endif


   ntime_in = 1
   if (ntime > 0) then
      allocate(time_in(ntime), clim_type%time_slice(ntime))
      allocate(clim_type%clim_times(12,(ntime+11)/12))
      time_in = 0.0_r8_kind
      clim_type%time_slice = set_time(0,0) + base_time
      clim_type%clim_times = set_time(0,0) + base_time
      call fms2_io_read_data(fileobj, "time", time_in)
      ntime_in = ntime
      !> convert the number of years passed to days if yearly=.true.
      if(yearly) then
         do n = 1, ntime
            if(.not.noleap_file_calendar) then
               ! Julian calendar
               num_years  = int(time_in(n))
               frac_year = time_in(n) - real(num_years, r8_kind)
               call get_date_julian(base_time, base_year, base_month, base_day, base_hour, base_minute, base_second)
               num_days = 0.0_r8_kind
               do nn=1, num_years
                  if( mod(base_year+nn-1,4)==0) num_days = num_days + 366._r8_kind
                  if( mod(base_year+nn-1,4)/=0) num_days = num_days + 365._r8_kind
               end do
               if( mod(base_year+num_years,4)==0) num_days = num_days + 366._r8_kind*frac_year
               if( mod(base_year+num_years,4)/=0) num_days = num_days + 365._r8_kind*frac_year
            else
               num_days = time_in(n)*365._r8_kind
            end if
            time_in(n)=num_days
         end do
      end if
      ! determine whether the data is a continuous set of monthly values or
      ! a series of annual cycles spread throughout the period of data
      non_monthly = .false.
      do n = 1, ntime-1
         !  Assume that the times in the data file correspond to days only.
         if (time_in(n+1) > (time_in(n) + 32._r8_kind)) then
           non_monthly = .true.
           exit
         endif
      end do
      if (clim_type%climatological_year) then
         call mpp_error (NOTE, 'interpolator_mod :'  // &
                      & trim(file_name) // ' is a year-independent climatology file')
      else
         call mpp_error (NOTE, 'interpolator_mod :' // &
                      & trim(file_name) // ' is a timeseries file')
      endif

      do n = 1, ntime
         !Assume that the times in the data file correspond to days only.


         if (clim_type%climatological_year) then
           !! RSH NOTE:
           !! for this case, do not add base_time. time_slice will be sent to
           !! time_interp_list with the optional argument modtime=YEAR, so that
           !! the time that is needed in time_slice is the displacement into the
           !! year, not the displacement from a base_time.
            clim_type%time_slice(n) = set_time( INT( (time_in(n)-real(INT(time_in(n)),r8_kind))*SECONDS_PER_DAY), &
                 INT(time_in(n)))
         else

           !--------------------------------------------------------------------
           !    if fileyr /= 0 (i.e., climatological_year=F),
           !    then define the times associated with each time-
           !    slice. if calendar conversion between data file and model calendar
           !    is needed, do it so that data from the file is associated with the
           !    same calendar time in the model. here the time_slice needs to
           !    include the base_time; values will be generated relative to the
           !    "real" time.
           !--------------------------------------------------------------------
           if ( (model_calendar == JULIAN .and. trim(adjustl(lowercase(file_calendar))) == 'julian')  .or. &
                & (model_calendar == NOLEAP .and. trim(adjustl(lowercase(file_calendar))) == 'noleap') )  then

             !---------------------------------------------------------------------
             !    no calendar conversion needed.
             !---------------------------------------------------------------------
             clim_type%time_slice(n) = &
               set_time( INT( (time_in(n)-real(INT(time_in(n)),r8_kind))*SECONDS_PER_DAY), &
               INT(time_in(n))) + base_time

           !---------------------------------------------------------------------
           !    convert file times from noleap to julian.
           !---------------------------------------------------------------------
           else if ( (model_calendar == JULIAN .and. trim(adjustl(lowercase(file_calendar))) == 'noleap')) then
             Noleap_time = set_time (0, INT(time_in(n))) + base_time
             call get_date_no_leap (Noleap_time, yr, mo, dy, hr, mn, sc)
             clim_type%time_slice(n) = set_date_julian (yr, mo, dy, hr, mn, sc)
             if (n == 1) then
               call print_date (clim_type%time_slice(1), &
                       str= 'for file ' // trim(file_name) // ', the first time slice is mapped to :')
             endif
             if (n == ntime) then
               call print_date (clim_type%time_slice(ntime), &
                        str= 'for file ' // trim(file_name) // ', the last time slice is mapped to:')
             endif


           !---------------------------------------------------------------------
           !    convert file times from julian to noleap.
           !---------------------------------------------------------------------
           else if ( (model_calendar == NOLEAP .and. trim(adjustl(lowercase(file_calendar))) == 'julian')) then
             Julian_time = set_time (0, INT(time_in(n))) + base_time
             call get_date_julian (Julian_time, yr, mo, dy, hr, mn, sc)
             clim_type%time_slice(n) = set_date_no_leap (yr, mo, dy,hr, mn, sc)
             if (n == 1) then
               call print_date (clim_type%time_slice(1), &
                                str= 'for file ' // trim(file_name) // ', the first time slice is mapped to :')
             endif
             if (n == ntime) then
               call print_date (clim_type%time_slice(ntime), &
                                str= 'for file ' // trim(file_name) // ', the last time slice is mapped to:')
             endif

           !---------------------------------------------------------------------
           !    any other calendar combinations would have caused a fatal error
           !    above.
           !---------------------------------------------------------------------
           endif
         endif

         m = (n-1)/12 +1 ; m1 = n- (m-1)*12
         clim_type%clim_times(m1,m) = clim_type%time_slice(n)
      enddo
   else
      allocate(time_in(1), clim_type%time_slice(1))
      allocate(clim_type%clim_times(1,1))
      time_in = 0.0_r8_kind
      clim_type%time_slice = set_time(0,0) + base_time
      clim_type%clim_times(1,1) = set_time(0,0) + base_time
   endif
   deallocate(time_in)
else
   call mpp_error(FATAL, 'Interpolator_init: time axis does not exist in file '//trim(src_file))
endif

!Assume that the horizontal interpolation within a file is the same for each variable.

 if (conservative_interp) then
    call horiz_interp_new (clim_type%interph, &
                           clim_type%FMS_INTP_TYPE_%lonb, clim_type%FMS_INTP_TYPE_%latb, &
                           lonb_mod, latb_mod)
 else

    call mpp_error(NOTE, "Using Bilinear interpolation")

    !!! DEBUG CODE
    if (.not. allocated(agrid_mod)) then
       nx = size(lonb_mod,1)-1
       ny = size(latb_mod,2)-1
       allocate(agrid_mod(nx,ny,2))
       do j=1,ny
       do i=1,nx
          call cell_center2((/lonb_mod(i,j),latb_mod(i,j)/),     &
                            (/lonb_mod(i+1,j),latb_mod(i+1,j)/), &
                            (/lonb_mod(i,j+1),latb_mod(i,j+1)/), &
                            (/lonb_mod(i+1,j+1),latb_mod(i+1,j+1)/), agrid_mod(i,j,:))
       enddo
       enddo
    endif

    !!! END DEBUG CODE

    call horiz_interp_new (clim_type%interph, &
                           clim_type%FMS_INTP_TYPE_%lonb, clim_type%FMS_INTP_TYPE_%latb, &
                           agrid_mod(:,:,1), agrid_mod(:,:,2), interp_method="bilinear")
 endif

!--------------------------------------------------------------------
!  allocate the variable clim_type%data . This will be the climatology
!  data horizontally interpolated, so it will be on the model horizontal
!  grid, but it will still be on the climatology vertical grid.
!--------------------------------------------------------------------

select case(ntime)
 case (13:)
! This may  be data that does not have a continous time-line
! i.e. IPCC data where decadal data is present but we wish to retain
! the seasonal nature of the data.
!! RSH: the following test will not always work; instead use the
!! RSH: non-monthly variable to test on.
!RSHlast_time = clim_type%time_slice(1) + ( ntime -1 ) * &
!RSH        ( clim_type%time_slice(2) - clim_type%time_slice(1) )

!RSHif ( last_time < clim_type%time_slice(ntime)) then

 if (non_monthly) then
! We have a broken time-line. e.g. We have monthly data but only for years ending in 0. 1960,1970 etc.
!   allocate(clim_type%data5d(size(lonb_mod(:))-1, size(latb_mod(:))-1, nlev, 2, num_fields))
   allocate(clim_type%FMS_INTP_TYPE_%pmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%FMS_INTP_TYPE_%pmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%FMS_INTP_TYPE_%nmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%FMS_INTP_TYPE_%nmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   clim_type%FMS_INTP_TYPE_%pmon_pyear = 0.0_lkind
   clim_type%FMS_INTP_TYPE_%pmon_nyear = 0.0_lkind
   clim_type%FMS_INTP_TYPE_%nmon_nyear = 0.0_lkind
   clim_type%FMS_INTP_TYPE_%nmon_pyear = 0.0_lkind
   clim_type%TIME_FLAG = BILINEAR
else
! We have a continuous time-line so treat as for 5-12 timelevels as below.
   if ( .not. read_all_on_init) then
   allocate(clim_type%FMS_INTP_TYPE_%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%FMS_INTP_TYPE_%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
            ntime, num_fields))
   endif
   clim_type%FMS_INTP_TYPE_%data5d = 0.0_lkind
   clim_type%TIME_FLAG = LINEAR
endif


!++lwh
 case (1:12)
!--lwh
! We have more than 4 timelevels
! Assume we have monthly or higher time resolution datasets (climatology or time series)
! So we only need to read 2 datasets and apply linear temporal interpolation.
   if ( .not. read_all_on_init) then
   allocate(clim_type%FMS_INTP_TYPE_%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%FMS_INTP_TYPE_%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
            ntime, num_fields))
   endif
   clim_type%FMS_INTP_TYPE_%data5d = 0.0_lkind
   clim_type%TIME_FLAG = LINEAR
!++lwh
!case (1:4)
! Assume we have seasonal data and read in all the data.
! We can apply sine curves to these data.

!  allocate(clim_type%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, ntime, num_fields))
!  clim_type%data5d = 0.0
!  clim_type%TIME_FLAG = SEASONAL
!--lwh
! case (default)
 case(:0)
   clim_type%TIME_FLAG = NOTIME
   allocate(clim_type%FMS_INTP_TYPE_%data5d(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 1, num_fields))
end select


!------------------------------------------------------------------
!    Allocate space for the single time level of the climatology on its
!    grid size.
!----------------------------------------------------------------------

   if(clim_type%TIME_FLAG .eq. LINEAR ) then
   allocate(clim_type%time_init(num_fields,2))
   else
   allocate(clim_type%time_init(num_fields,ntime))
   endif
   allocate (clim_type%indexm(num_fields),   &
             clim_type%indexp(num_fields),   &
             clim_type%climatology(num_fields))
   clim_type%time_init(:,:) = 0
   clim_type%indexm(:)      = 0
   clim_type%indexp(:)      = 0
   clim_type%climatology(:) = 0

allocate(clim_type%has_level(num_fields))
allocate(clim_type%field_name(num_fields))
allocate(clim_type%mr(num_fields))
allocate(clim_type%out_of_bounds(num_fields))
clim_type%out_of_bounds(:)=0
allocate(clim_type%vert_interp(num_fields))
clim_type%vert_interp(:)=0

if(present(data_names)) then

!++lwh
   if ( size(data_out_of_bounds(:)) /= size(data_names(:)) .and. size(data_out_of_bounds(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the data_out_of_bounds array must be 1&
                            & or size(data_names)')
   if (present(vert_interp)) then
      if( size(vert_interp(:)) /= size(data_names(:)) .and. size(vert_interp(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the vert_interp array must be 1&
                            & or size(data_names)')
   endif
! Only read the fields named in data_names
   do j=1,size(data_names(:))
      if(variable_exists(fileobj, data_names(j)) ) then
         call get_variable_units(fileobj, data_names(j), units)
         ndim = get_variable_num_dimensions(fileobj, data_names(j))
         clim_type%has_level(j) = .false.
         if(ndim > 2) then
            call get_variable_dimension_names(fileobj, data_names(j), var_dimname(1:ndim))
            if(trim(var_dimname(3)) == "pfull" .OR. trim(var_dimname(3)) == "sigma_full") &
               & clim_type%has_level(j) = .true.
         endif

         units=chomp(units)
         if (mpp_pe() == 0 ) write(*,*) 'Initializing src field : ',trim(data_names(j))
         clim_type%field_name(j) = data_names(j)
         clim_type%mr(j)         = check_climo_units(units)
         if (present(clim_units)) clim_units(j) = units
         clim_type%out_of_bounds(j) = data_out_of_bounds( MIN(j,SIZE(data_out_of_bounds(:))) )
         if( clim_type%out_of_bounds(j) /= CONSTANT .and. clim_type%out_of_bounds(j) /= ZERO ) &
             call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                                    & set to ZERO or CONSTANT")
         if( present(vert_interp) ) then
            clim_type%vert_interp(j) = vert_interp( MIN(j,SIZE(vert_interp(:))) )
            if( clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. clim_type%vert_interp(j) /= INTERP_LINEAR_P ) &
               call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                    & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
         else
            clim_type%vert_interp(j) = INTERP_WEIGHTED_P
         end if
      else
         call mpp_error(FATAL,'interpolator_init : Check names of fields being passed. ' &
                              &//trim(data_names(j))//' does not exist.')
      endif
   enddo
else
   if ( size(data_out_of_bounds(:)) /= nvar .and. size(data_out_of_bounds(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the out of bounds array must be 1&
                           & or the number of fields in the climatology dataset')
   if ( present(vert_interp) ) then
      if (size(vert_interp(:)) /= nvar .and. size(vert_interp(:)) /= 1 ) &
      call mpp_error(FATAL,'interpolator_init : The size of the vert_interp array must be 1&
                           & or the number of fields in the climatology dataset')
   endif

! Read all the fields within the climatology data file.
   j = 0
   do i=1,nvar
      ndim = var_ndims(i)
      if(ndim .LE. 1) cycle
      j = j + 1
      clim_type%has_level(j) = .false.
      if(ndim > 2) then
         call get_variable_dimension_names(fileobj, var_names(i), var_dimname(1:ndim))
         if(trim(var_dimname(3)) == "pfull" .OR. trim(var_dimname(3)) == "sigma_full") clim_type%has_level(j) = .true.
      endif

      call get_variable_units(fileobj, var_names(i), units)
      if (mpp_pe() ==0 ) write(*,*) 'Initializing src field : ',trim(var_names(i))
      clim_type%field_name(j) = trim(var_names(i))
      clim_type%mr(j)         = check_climo_units(units)
      if (present(clim_units)) clim_units(j) = units
      clim_type%out_of_bounds(j) = data_out_of_bounds( MIN(i,SIZE(data_out_of_bounds(:))) )
      if( clim_type%out_of_bounds(j) /= CONSTANT .and. &
          clim_type%out_of_bounds(j) /= ZERO ) &
         call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                              & set to ZERO or CONSTANT")
      if( present(vert_interp) ) then
         clim_type%vert_interp(j) = vert_interp( MIN(i,SIZE(vert_interp(:))) )
         if( clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. &
             clim_type%vert_interp(j) /= INTERP_LINEAR_P ) &
            call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                 & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
      else
         clim_type%vert_interp(j) = INTERP_WEIGHTED_P
      end if
   end do
   if(j .NE. num_fields) call mpp_error(FATAL,"Interpolator_init: j does not equal to num_fields")
   deallocate(var_names, var_ndims)
!--lwh
endif

if( clim_type%TIME_FLAG .eq. SEASONAL ) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call read_data( clim_type, clim_type%field_name(i), &
                         clim_type%FMS_INTP_TYPE_%data5d(:,:,:,n,i), n, i, base_time )
      enddo
   enddo
endif

if( clim_type%TIME_FLAG .eq. LINEAR  .and. read_all_on_init) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call read_data( clim_type, clim_type%field_name(i), &
                         clim_type%FMS_INTP_TYPE_%data5d(:,:,:,n,i), n, i, base_time )
      enddo
   enddo

   call close_file (fileobj)
endif

if( clim_type%TIME_FLAG .eq. NOTIME ) then
! Read all the data at this point.
   do i=1,num_fields
     call read_data_no_time_axis( clim_type, clim_type%field_name(i), &
                                  clim_type%FMS_INTP_TYPE_%data5d(:,:,:,1,i), i )
   enddo
   call close_file (fileobj)
endif

if (present (single_year_file)) then
  single_year_file = clim_type%climatological_year
endif

end subroutine FMS2IO_INTERPOLATOR_INIT_

subroutine GET_AXIS_LATLON_DATA_(fileobj, name, latlon_data)
   type(FmsNetcdfFile_t), intent(in) :: fileobj
   character(len=*),      intent(in) :: name
   real(FMS_INTP_KIND_), dimension(:),   intent(out) :: latlon_data


   if(variable_exists(fileobj, name)) then
      call fms2_io_read_data(fileobj, name, latlon_data)
   else
      call mpp_error(FATAL,'get_axis_latlon_data: variable '// &
                     & trim(name)//' does not exist in file '//trim(fileobj%path) )
   endif
   call get_variable_units(fileobj, name, units)
   select case(units(1:6))
   case('degree')
      latlon_data = latlon_data*real(DTR,FMS_INTP_KIND_)
   case('radian')
   case default
      call mpp_error(FATAL, "get_axis_latlon_data : Units for '// &
                     & trim(name)//' not recognised in file "//trim(fileobj%path))
   end select

 end subroutine GET_AXIS_LATLON_DATA_


subroutine GET_AXIS_LEVEL_DATA_(fileobj, name, level_data, level_type, vertical_indices)
   type(FmsNetcdfFile_t), intent(in) :: fileobj
   character(len=*),      intent(in) :: name
   real(FMS_INTP_KIND_), dimension(:),   intent(out) :: level_data
   integer,              intent(out) :: level_type, vertical_indices
   real(FMS_INTP_KIND_), dimension(:), allocatable   :: alpha
   integer                           :: n, nlev

   integer, parameter :: lkind=FMS_INTP_KIND_

   if(variable_exists(fileobj, name)) then
      call fms2_io_read_data(fileobj, name, level_data)
   else
      call mpp_error(FATAL,'get_axis_level_data: variable '// &
                     & trim(name)//' does not exist in file '//trim(fileobj%path) )
   endif
   call get_variable_units(fileobj, name, units)
   level_type = PRESSURE
   ! Convert to Pa
   if( trim(adjustl(lowercase(chomp(units)))) == "mb" .or. trim(adjustl(lowercase(chomp(units)))) == "hpa") then
      level_data = level_data * 100._lkind
   endif
   nlev = size(level_data(:))
   sense = get_variable_sense(fileobj, name)
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
      if( sense == 1 ) then
        vertical_indices = INCREASING_UPWARD
        allocate (alpha(nlev))
        do n = 1, nlev
          alpha(n) = level_data(nlev-n+1)
        end do
        do n = 1, nlev
          level_data(n) = alpha(n)
        end do
        deallocate (alpha)
      else
        vertical_indices = INCREASING_DOWNWARD
      endif


end subroutine GET_AXIS_LEVEL_DATA_


!
!---------------------------------------------------------------------
!
!> @brief cell_center2 receives the variables q1, q2, q3, and q4
!!        as inputs and returns e2.
!!
!! @param [in] <q1> No description
!! @param [in] <q2> No description
!! @param [in] <q3> No description
!! @param [in] <q4> No description
!! @param [out] <e2> No description
 subroutine CELL_CENTER2_(q1, q2, q3, q4, e2)
      real(FMS_INTP_KIND_) , intent(in ) :: q1(2), q2(2), q3(2), q4(2)
      real(FMS_INTP_KIND_) , intent(out) :: e2(2)
! Local
      real(FMS_INTP_KIND_) :: p1(3), p2(3), p3(3), p4(3)
      real(FMS_INTP_KIND_) :: ec(3,1)
      real(FMS_INTP_KIND_) :: dd
      integer k

      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

      do k=1,3
         ec(k,1) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1,1)**2 + ec(2,1)**2 + ec(3,1)**2 )

      do k=1,3
         ec(k,1) = ec(k,1) / dd
      enddo

      call cart_to_latlon(1, ec, e2(1:1), e2(2:2))

 end subroutine CELL_CENTER2_
!
!---------------------------------------------------------------------
!> @brief car_to_latlon receives the variables np, q, xs, and ys
!!        as inputs and returns q, xs, and ys.
!!
!! @param [in] <np> No description
!! @param [inout] <q> No description
!! @param [inout] <xs> No description
!! @param [inout] <ys> No description
 subroutine CART_TO_LATLON_(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real(FMS_INTP_KIND_), intent(inout):: q(3,np)
  real(FMS_INTP_KIND_), intent(inout):: xs(np), ys(np)
! local
  real(f_p), parameter:: esl=1.e-10_f_p
  real(f_p):: p(3)
  real(f_p):: dist, lat, lon
  integer i,k

  do i=1,np
     do k=1,3
        p(k) = real(q(k,i), f_p)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0._f_p
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0._f_p ) lon = 2._f_p*real(PI,f_p) + lon
     lat = asin(p(3))

     xs(i) = real(lon,FMS_INTP_KIND_)
     ys(i) = real(lat,FMS_INTP_KIND_)
! q Normalized:
     do k=1,3
        q(k,i) = real(p(k),FMS_INTP_KIND_)
     enddo
  enddo

 end  subroutine CART_TO_LATLON_
!
!---------------------------------------------------------------------
!> @brief latlon2xyz receives the variable p as input and returns e
!!        as output in order to map (lon, lat) to (x,y,z).
!!
!! @param [in] <p> No description
!! @param [out] <e> No description
 subroutine LATLON2XYZ_(p, e)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real(FMS_INTP_KIND_), intent(in) :: p(2)
 real(FMS_INTP_KIND_), intent(out):: e(3)

 integer :: n
 real (f_p):: q(2)
 real (f_p):: e1, e2, e3

    do n=1,2
       q(n) = real(p(n),f_p)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = real(e1,FMS_INTP_KIND_)
    e(2) = real(e2,FMS_INTP_KIND_)
    e(3) = real(e3,FMS_INTP_KIND_)

 end subroutine LATLON2XYZ_

!
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief init_clim_diag is a routine to register diagnostic fields
!!        for the climatology file. This routine calculates the domain
!!        decompostion of the climatology fields for later export
!!        through send_data. The ids created here are for column
!!        burdens that will diagnose the vertical interpolation
!!        routine.
!!
!! @param [inout] <clim_type> The interpolate type containing the
!!                      names of the fields in the climatology file
!! @param [in] <mod_axes> The axes of the model
!! @param [in] <init_time> The model initialization time
!!
!! @throw FATAL, "init_clim_diag : You must call interpolator_init before calling init_clim_diag"
!! @throw FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data"
subroutine INIT_CLIM_DIAG_(clim_type, mod_axes, init_time)
!
! Routine to register diagnostic fields for the climatology file.
! This routine calculates the domain decompostion of the climatology fields
! for later export through send_data.
! The ids created here are for column burdens that will diagnose the vertical interpolation routine.
! climo_diag_id : 'module_name = climo' is intended for use with the model vertical resolution.
! hinterp_id    : 'module_name = 'hinterp' is intended for use with the climatology vertical resolution.

! INTENT INOUT :
!    clim_type : The interpolate type containing the names of the fields in the climatology file.
!
! INTENT IN    :
!   mod_axes   : The axes of the model.
!   init_time  : The model initialization time.
!
type(interpolate_type), intent(inout)  :: clim_type
integer               , intent(in)     :: mod_axes(:)
type(time_type)       , intent(in)     :: init_time

integer :: axes(2),nxd,nyd,ndivs,i
type(domain2d) :: domain
integer :: domain_layout(2), iscomp, iecomp,jscomp,jecomp


if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "init_clim_diag : You must call interpolator_init before calling init_clim_diag")


ndivs = mpp_npes()
nxd = size(clim_type%FMS_INTP_TYPE_%lon(:))
nyd = size(clim_type%FMS_INTP_TYPE_%lat(:))

! Define the domain decomposition of the climatology file. This may be (probably is) different from the model domain.
call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)
call mpp_get_compute_domain(domain, iscomp, iecomp, jscomp, jecomp)
axes(1) = diag_axis_init(clim_type%file_name(1:5)//'x',clim_type%FMS_INTP_TYPE_%lon,units='degrees',&
                         cart_name='x',domain2=domain)
axes(2) = diag_axis_init(clim_type%file_name(1:5)//'y',clim_type%FMS_INTP_TYPE_%lat,units='degrees',&
                         cart_name='y',domain2=domain)
clim_type%is = iscomp
clim_type%ie = iecomp
clim_type%js = jscomp
clim_type%je = jecomp

!init_time = set_date(1980,1,1,0,0,0)

if ((num_clim_diag + size(clim_type%field_name(:))) .gt. max_diag_fields )  &
    call mpp_error(FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data")
do i=1,size(clim_type%field_name(:))
climo_diag_name(i+num_clim_diag) = clim_type%field_name(i)
climo_diag_id(i+num_clim_diag) =  register_diag_field('climo',clim_type%field_name(i),axes(1:2),init_time,&
                                                      'climo_'//clim_type%field_name(i), 'kg/kg', missing_value)
hinterp_id(i+num_clim_diag) =  register_diag_field('hinterp',clim_type%field_name(i),mod_axes(1:2),init_time,&
                                                   'interp_'//clim_type%field_name(i),'kg/kg' , missing_value)
enddo
! Total number of climatology diagnostics (num_clim_diag). This can be from multiple climatology
! fields with different spatial axes.
! It is simply a holder for the diagnostic indices.
num_clim_diag = num_clim_diag+size(clim_type%field_name(:))

clim_diag_initialized = .true.

end subroutine INIT_CLIM_DIAG_



!
!---------------------------------------------------------------------
!> @brief obtain_interpolator_time_slices makes sure that the
!!        appropriate time slices are available for interpolation on
!!        this time step.
!!
!! @param [inout] <clim_type> The interpolate type previously defined
!!                      by a call to interpolator_init
!! @param [in] <Time> The model time that you wish to interpolate to
!!
!! @throw FATAL "interpolator_timeslice 1:  file="
!! @throw FATAL "interpolator_timeslice 2:  file="
!! @throw FATAL "interpolator_timeslice 3:  file="
!! @throw FATAL "interpolator_timeslice 4:  file="
!! @throw FATAL "interpolator_timeslice 5:  file="
!! @throw FATAL "interpolator_timeslice : No data from the previous
!!                    climatology time but we have the next time. How did
!!                    this happen?"
subroutine OBTAIN_INTERPOLATOR_TIME_SLICES_ (clim_type, Time)

!  Makes sure that appropriate time slices are available for interpolation
!  on this time step
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   Time        : The model time that you wish to interpolate to.

type(interpolate_type), intent(inout)  :: clim_type
type(time_type)       , intent(in)  :: Time

integer :: taum, taup
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, n
character(len=256) :: err_msg

integer, parameter :: lkind=FMS_INTP_KIND_

    if (clim_type%climatological_year) then
       !++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, &
                           taum, taup, modtime=YEAR, err_msg=err_msg )
          if(trim(err_msg) /= '') then
             call mpp_error(FATAL,'interpolator_timeslice 1: '// &
                            & trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
          endif
       else
          taum = 1
          taup = 1
          clim_type%FMS_INTP_TYPE_%tweight = 0._lkind
       end if
       !--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, taum, taup, err_msg=err_msg )
       if(trim(err_msg) /= '') then
          call mpp_error(FATAL,'interpolator_timeslice 2: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
       endif
    endif

    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
      ! the climatology year into the appropriate place.


      ! We need to get the previous months data for the climatology year before
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then
            if ( modday <  day ) then
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif

        enddo
        if ( indexm == 0 ) then
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
          call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
          call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

        call time_interp(Time, month, clim_type%FMS_INTP_TYPE_%tweight3, taum, taup, err_msg=err_msg ) ! tweight3 is
                                                                                 !! the time weight between the months.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight3 = 1._lkind
        if(trim(err_msg) /= '') then
          call mpp_error(FATAL,'interpolator_timeslice 3: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%FMS_INTP_TYPE_%tweight1, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !>  protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight1 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_timeslice 4: '// &
                          & trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%FMS_INTP_TYPE_%tweight2, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight2 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_timeslice 5: '// &
                          & trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        if (indexm == clim_type%indexm(1) .and.  &
            indexp == clim_type%indexp(1) .and. &
            climatology == clim_type%climatology(1)) then
        else
          clim_type%indexm(:) = indexm
          clim_type%indexp(:) = indexp
          clim_type%climatology(:) = climatology
          do i=1, size(clim_type%field_name(:))
               call read_data(clim_type,clim_type%field_name(i),  &
                    clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i),   &
                    clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
               call read_data(clim_type,clim_type%field_name(i),  &
                    clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i),   &
                    clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
               call read_data(clim_type,clim_type%field_name(i),  &
                    clim_type%FMS_INTP_TYPE_%pmon_nyear(:,:,:,i),  &
                    clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
               call read_data(clim_type,clim_type%field_name(i),  &
                    clim_type%FMS_INTP_TYPE_%nmon_nyear(:,:,:,i),  &
                    clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
           end do
         endif


      else ! We are within a climatology data set

          do i=1, size(clim_type%field_name(:))
             if (taum /= clim_type%time_init(i,1) .or. &
                taup /= clim_type%time_init(i,2) ) then


              call read_data(clim_type,clim_type%field_name(i),   &
                           clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
              call read_data(clim_type,clim_type%field_name(i),   &
                           clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i), taup,i,Time)
              clim_type%time_init(i,1) = taum
              clim_type%time_init(i,2) = taup
             endif
           end do
!       clim_type%pmon_nyear = 0.0
!       clim_type%nmon_nyear = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(:) = 0
          clim_type%indexp(:) = 0
          clim_type%climatology(:) = 0


!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%FMS_INTP_TYPE_%tweight1 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight2 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight3 = clim_type%FMS_INTP_TYPE_%tweight
      endif
    endif   !(BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(1,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(1,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
      do i=1, size(clim_type%field_name(:))
          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,1,i), taum,i,Time)
               clim_type%time_init(i,1) = taum
               clim_type%itaum = 1
          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,2,i), taup,i,Time)
               clim_type%time_init(i,2) = taup
               clim_type%itaup = 2
      end do
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_timeslice : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
        do i=1, size(clim_type%field_name(:))
           call read_data(clim_type,clim_type%field_name(i), &
                          clim_type%FMS_INTP_TYPE_%data5d(:,:,:,clim_type%itaup,i), taup,i, Time)
           clim_type%time_init(i,clim_type%itaup)=taup
        end do
      endif


    endif! TIME_FLAG

    clim_type%separate_time_vary_calc = .true.

!-------------------------------------------------------------------


 end subroutine OBTAIN_INTERPOLATOR_TIME_SLICES_


!#####################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_4D receives a field name as input and
!!        interpolates the field to model a 4D grid and time axis.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <Time> The model time that you wish to interpolate to
!! @param [in] <phalf> The half level model pressure field
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_4D : You must call interpolator_init
!!                    before calling interpolator"
!! @throw FATAL "interpolator_mod: cannot use 4D interface to interpolator for this file"
!! @throw FATAL "interpolator_4D 1:  file="
!! @throw FATAL "interpolator_4D 2:  file="
!! @throw FATAL "interpolator_4D 3:  file="
!! @throw FATAL "interpolator_4D 4:  file="
!! @throw FATAL "interpolator_4D 5:  file="
!! @throw FATAL "interpolator_3D : No data from the previous climatology
!!                    time but we have the next time. How did this happen?"
!! @throw NOTE "Interpolator: model surface pressure is greater than
!!              climatology surface pressure for "
!! @throw NOTE "Interpolator: model top pressure is less than
!!                    climatology top pressure for "
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_4D_(clim_type, Time, phalf, interp_data,  &
                           field_name, is,js, clim_units)
!
! Return 4-D field interpolated to model grid and time
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of a field that you wish to interpolate.
!                 all variables within this interpolate_type variable
!                 will be interpolated on this call. field_name may
!                 be any one of the variables.
!   Time        : The model time that you wish to interpolate to.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model fields with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
type(time_type)       , intent(in)  :: Time
real(FMS_INTP_KIND_), dimension(:,:,:), intent(in)  :: phalf
real(FMS_INTP_KIND_), dimension(:,:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup, ilon
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),size(interp_data,2),&
                        size(clim_type%FMS_INTP_TYPE_%levs(:)),size(clim_type%field_name(:)))
real(FMS_INTP_KIND_) :: p_fact(size(interp_data,1),size(interp_data,2))
real(FMS_INTP_KIND_) :: col_data(size(interp_data,1),size(interp_data,2),   &
                           size(clim_type%field_name(:)))
real(FMS_INTP_KIND_) :: pclim(size(clim_type%FMS_INTP_TYPE_%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, j, k, n
character(len=256) :: err_msg

integer, parameter :: lkind=FMS_INTP_KIND_

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_4D : You must call interpolator_init before calling interpolator")

   do n=2,size(clim_type%field_name(:))
     if (clim_type%vert_interp(n) /= clim_type%vert_interp(n-1) .or. &
      clim_type%out_of_bounds(n) /= clim_type%out_of_bounds(n-1)) then
       if (mpp_pe() == mpp_root_pe() ) then
         print *, 'processing file ' // trim(clim_type%file_name)
       endif
       call mpp_error (FATAL, 'interpolator_mod: &
               &cannot use 4D interface to interpolator for this file')
     endif
   end do




istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

  do i= 1,size(clim_type%field_name(:))
!!++lwh
   if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then
!--lwh
    found_field=.true.
    exit
 endif
end do
   i = 1

    if(present(clim_units)) then
         call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
    endif


!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------

if ( .not. clim_type%separate_time_vary_calc) then
!   print *, 'TIME INTERPOLATION NOT SEPARATED 4d--',  &
!                                trim(clim_type%file_name), mpp_pe()

    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, &
                           taum, taup, modtime=YEAR, err_msg=err_msg )
          if(trim(err_msg) /= '') then
             call mpp_error(FATAL,'interpolator_4D 1: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
          endif
       else
          taum = 1
          taup = 1
          clim_type%FMS_INTP_TYPE_%tweight = 0._lkind
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, taum, taup, err_msg=err_msg )
       if(trim(err_msg) /= '') then
          call mpp_error(FATAL,'interpolator_4D 2: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
       endif
    endif

    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
      ! the climatology year into the appropriate place.


      ! We need to get the previous months data for the climatology year before
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then
            if ( modday <  day ) then
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif

        enddo
        if ( indexm == 0 ) then
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
          call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
          call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                        climyear, climmonth, climday, climhour, climminute, climsecond)
          month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

        call time_interp(Time, month, clim_type%FMS_INTP_TYPE_%tweight3, taum, taup, err_msg=err_msg ) ! tweight3 is
                                                                                 !! the time weight between the months.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight3 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_4D 3: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%FMS_INTP_TYPE_%tweight1, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight1 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_4D 4: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif
        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%FMS_INTP_TYPE_%tweight2, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight2 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_4D 5: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        if (indexm == clim_type%indexm(1) .and.  &
            indexp == clim_type%indexp(1) .and. &
            climatology == clim_type%climatology(1)) then
        else
          clim_type%indexm(:) = indexm
          clim_type%indexp(:) = indexp
          clim_type%climatology(:) = climatology
          do i=1, size(clim_type%field_name(:))
            call read_data(clim_type,clim_type%field_name(i),  &
                 clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i),   &
             clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_name(i),  &
             clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i),   &
             clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
            call read_data(clim_type,clim_type%field_name(i),  &
              clim_type%FMS_INTP_TYPE_%pmon_nyear(:,:,:,i),  &
              clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
            call read_data(clim_type,clim_type%field_name(i),  &
              clim_type%FMS_INTP_TYPE_%nmon_nyear(:,:,:,i),  &
              clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
          end do
        endif



      else ! We are within a climatology data set

        do i=1, size(clim_type%field_name(:))
          if (taum /= clim_type%time_init(i,1) .or. &
              taup /= clim_type%time_init(i,2) ) then


            call read_data(clim_type,clim_type%field_name(i),   &
                           clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type,clim_type%field_name(i),   &
                           clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i), taup,i,Time)
            clim_type%time_init(i,1) = taum
            clim_type%time_init(i,2) = taup
          endif
        end do
!       clim_type%pmon_nyear = 0.0
!       clim_type%nmon_nyear = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(:) = 0
          clim_type%indexp(:) = 0
          clim_type%climatology(:) = 0


!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%FMS_INTP_TYPE_%tweight1 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight2 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight3 = clim_type%FMS_INTP_TYPE_%tweight
      endif
    endif   !(BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(1,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(1,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    do i=1, size(clim_type%field_name(:))
    call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
    call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
    end do
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_3D : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
    do i=1, size(clim_type%field_name(:))
        call read_data(clim_type,clim_type%field_name(i), &
                       clim_type%FMS_INTP_TYPE_%data5d(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
     end do
      endif


    endif! TIME_FLAG


endif ! (.not. separate_time_vary_calc)


select case(clim_type%TIME_FLAG)
  case (LINEAR)
    do n=1, size(clim_type%field_name(:))
      hinterp_data(:,:,:,n) = (1._lkind-clim_type%FMS_INTP_TYPE_%tweight)*  &
                clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,clim_type%itaum,n)  +  &
                                 clim_type%FMS_INTP_TYPE_%tweight*   &
                clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,clim_type%itaup,n)
    end do
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    do n=1, size(clim_type%field_name(:))
      hinterp_data(:,:,:,n) = (1._lkind-clim_type%FMS_INTP_TYPE_%tweight1)*&
                              (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3)*   &
                              clim_type%FMS_INTP_TYPE_%pmon_pyear(istart:iend,jstart:jend,:,n) + &
                              (1._lkind-clim_type%FMS_INTP_TYPE_%tweight2)*clim_type%FMS_INTP_TYPE_%tweight3* &
                              clim_type%FMS_INTP_TYPE_%nmon_pyear(istart:iend,jstart:jend,:,n) + &
                              clim_type%FMS_INTP_TYPE_%tweight1* (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3)*&
                              clim_type%FMS_INTP_TYPE_%pmon_nyear(istart:iend,jstart:jend,:,n) + &
                              clim_type%FMS_INTP_TYPE_%tweight2* clim_type%FMS_INTP_TYPE_%tweight3*   &
                              clim_type%FMS_INTP_TYPE_%nmon_nyear(istart:iend,jstart:jend,:,n)

    end do

end select

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0_lkind
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

col_data(:,:,:)=0.0_lkind
     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(NO_CONV)
    do k = 1,size(hinterp_data,3)
   col_data(:,:,i) = col_data(:,:,i) + hinterp_data(:,:,k,i)* &
      (clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))/real(grav,FMS_INTP_KIND_)
    enddo

  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       col_data(:,:,i) = col_data(:,:,i) + hinterp_data(:,:,k,i)
       hinterp_data(:,:,k,i) = hinterp_data(:,:,k,i)/ &
         ((clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))*p_fact)
    enddo
end select
    enddo

     do i= 1, size(clim_type%field_name(:))
found = .false.
do j = 1,size(climo_diag_name(:))
  if ( trim(adjustl(lowercase(climo_diag_name(j)))) .eq. trim(adjustl(lowercase(clim_type%field_name(i))))) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),col_data(:,:,i),Time,is_in=istart,js_in=jstart)
  endif
endif

  end do

   i = 1

!++lwh
do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%FMS_INTP_TYPE_%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             & // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,:),interp_data(ilon,j,:,:))
         case(INTERP_LINEAR_P)
          do n=1, size(clim_type%field_name(:))
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,n),interp_data(ilon,j,:,n))
          end do
!        case(INTERP_LOG)
      end select
   enddo
enddo

!--lwh
     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k,i) = interp_data(:,:,k,i)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

     end do

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine INTERPOLATOR_4D_
!
!#######################################################################
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_3D receives a field name as input and
!!        interpolates the field to model a 3D grid and time axis.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <Time> The model time that you wish to interpolate to
!! @param [in] <phalf> The half level model pressure field
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_3D : You must call interpolator_init
!!                    before calling interpolator"
!! @throw FATAL "interpolator_3D 1:  file="
!! @throw FATAL "interpolator_3D 2:  file="
!! @throw FATAL "interpolator_3D 3:  file="
!! @throw FATAL "interpolator_3D 4:  file="
!! @throw FATAL "interpolator_3D 5:  file="
!! @throw FATAL "interpolator_3D : No data from the previous climatology
!!                    time but we have the next time. How did this happen?"
!! @throw NOTE "Interpolator: model surface pressure is greater than
!!              climatology surface pressure for "
!! @throw NOTE "Interpolator: model top pressure is less than
!!                    climatology top pressure for "
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_3D_(clim_type, Time, phalf, interp_data,field_name, is,js, clim_units)
!
! Return 3-D field interpolated to model grid and time
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   Time        : The model time that you wish to interpolate to.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
type(time_type)       , intent(in)  :: Time
real(FMS_INTP_KIND_), dimension(:,:,:), intent(in)  :: phalf
real(FMS_INTP_KIND_), dimension(:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup, ilon
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),&
                        size(interp_data,2),size(clim_type%FMS_INTP_TYPE_%levs(:)))
real(FMS_INTP_KIND_) :: p_fact(size(interp_data,1),size(interp_data,2))
real(FMS_INTP_KIND_) :: col_data(size(interp_data,1),size(interp_data,2))
real(FMS_INTP_KIND_) :: pclim(size(clim_type%FMS_INTP_TYPE_%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: i, j, k, n
character(len=256) :: err_msg

integer, parameter :: lkind=FMS_INTP_KIND_

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_3D : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
!++lwh
  if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then
!--lwh
    found_field=.true.
    if(present(clim_units)) then
      call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
      clim_units = chomp(clim_units)
    endif

!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------

if ( .not. clim_type%separate_time_vary_calc) then
!   print *, 'TIME INTERPOLATION NOT SEPARATED 3d--',  &
!                                trim(clim_type%file_name), mpp_pe()
    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, &
                           taum, taup, modtime=YEAR, err_msg=err_msg )
          if(trim(err_msg) /= '') then
             call mpp_error(FATAL,'interpolator_3D 1: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
          endif
       else
          taum = 1
          taup = 1
          clim_type%FMS_INTP_TYPE_%tweight = 0._lkind
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, taum, taup, err_msg=err_msg )
       if(trim(err_msg) /= '') then
          call mpp_error(FATAL,'interpolator_3D 2: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
       endif
    endif

!   if(clim_type%TIME_FLAG .ne. LINEAR ) then
    if(clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init ) then
      clim_type%itaum=taum
      clim_type%itaup=taup
    endif

    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
      ! the climatology year into the appropriate place.


      ! We need to get the previous months data for the climatology year before
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then
            if ( modday <  day ) then
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif

        enddo
        if ( indexm == 0 ) then
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
        call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
        call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

        call time_interp(Time, month, clim_type%FMS_INTP_TYPE_%tweight3, taum, taup, err_msg=err_msg ) ! tweight3 is
                                                                                 !! the time weight between the months.
        if (taum==2 .and. taup==2) &
             clim_type%FMS_INTP_TYPE_%tweight3 = 1._lkind ! protect against post-perth time_interp behavior
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_3D 3: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%FMS_INTP_TYPE_%tweight1, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        if (taum==2 .and. taup==2) &
             clim_type%FMS_INTP_TYPE_%tweight1 = 1._lkind ! protect against post-perth time_interp behavior
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_3D 4: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%FMS_INTP_TYPE_%tweight2, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        if (taum==2 .and. taup==2) &
             clim_type%FMS_INTP_TYPE_%tweight2 = 1._lkind ! protect against post-perth time_interp behavior
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_3D 5: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif


        if (indexm == clim_type%indexm(i) .and.  &
          indexp == clim_type%indexp(i) .and. &
          climatology == clim_type%climatology(i)) then
        else
          clim_type%indexm(i) = indexm
          clim_type%indexp(i) = indexp
          clim_type%climatology(i) = climatology
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i),  &
            clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i),   &
            clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
          call read_data(clim_type,clim_type%field_name(i),   &
            clim_type%FMS_INTP_TYPE_%pmon_nyear(:,:,:,i),  &
            clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%nmon_nyear(:,:,:,i),  &
            clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
        endif




      else ! We are within a climatology data set

        if (taum /= clim_type%time_init(i,1) .or. &
            taup /= clim_type%time_init(i,2) ) then

          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i), taup,i,Time)
!RSHbug   clim_type%pmon_nyear = 0.0
!RSHbug   clim_type%nmon_nyear = 0.0

!         clim_type%pmon_nyear(:,:,:,i) = 0.0
!         clim_type%nmon_nyear(:,:,:,i) = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(i) = 0
          clim_type%indexp(i) = 0
          clim_type%climatology(i) = 0


          clim_type%time_init(i,1) = taum
          clim_type%time_init(i,2) = taup
        endif
!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%FMS_INTP_TYPE_%tweight1 = 0.0_lkind ; clim_type%FMS_INTP_TYPE_%tweight2 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight3 = clim_type%FMS_INTP_TYPE_%tweight
      endif

    endif ! (BILINEAR)


    if(clim_type%TIME_FLAG .eq. LINEAR  .and.   &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(i,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(i,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
    call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolator_3D : No data from the previous climatology time &
                         & but we have the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
!We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
        call read_data(clim_type,clim_type%field_name(i), &
                       clim_type%FMS_INTP_TYPE_%data5d(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
      endif


    endif! TIME_FLAG

    endif   !( .not. clim_type%separate_time_vary_calc)
select case(clim_type%TIME_FLAG)
  case (LINEAR)
    hinterp_data = (1._lkind-clim_type%FMS_INTP_TYPE_%tweight) * &
                    clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,clim_type%itaum,i) + &
                    clim_type%FMS_INTP_TYPE_%tweight * &
                    clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,clim_type%itaup,i)
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    hinterp_data = &
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight1)*&
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3)*clim_type%FMS_INTP_TYPE_%pmon_pyear(istart:iend,jstart:jend,:,i) + &
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight2) * clim_type%FMS_INTP_TYPE_%tweight3  &
    * clim_type%FMS_INTP_TYPE_%nmon_pyear(istart:iend,jstart:jend,:,i) + &
    clim_type%FMS_INTP_TYPE_%tweight1 * (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3) * &
    clim_type%FMS_INTP_TYPE_%pmon_nyear(istart:iend,jstart:jend,:,i) + &
    clim_type%FMS_INTP_TYPE_%tweight2 * clim_type%FMS_INTP_TYPE_%tweight3 * &
    clim_type%FMS_INTP_TYPE_%nmon_nyear(istart:iend,jstart:jend,:,i)



end select

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0_lkind
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select



col_data(:,:)=0.0_lkind
select case(clim_type%mr(i))
  case(NO_CONV)
    do k = 1,size(hinterp_data,3)
   col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)* &
      (clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))/real(grav,FMS_INTP_KIND_)
    enddo

  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)
       hinterp_data(:,:,k) = hinterp_data(:,:,k)/ &
         ((clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))*p_fact)
    enddo
end select

found = .false.
do j = 1,size(climo_diag_name(:))
  if ( trim(adjustl(lowercase(climo_diag_name(j)))) .eq. trim(adjustl(lowercase(clim_type%field_name(i))))) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),col_data,Time,is_in=istart,js_in=jstart)
  endif
endif


!++lwh
do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%FMS_INTP_TYPE_%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
         case(INTERP_LINEAR_P)
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
!        case(INTERP_LOG)
      end select
   enddo
enddo
!--lwh

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k) = interp_data(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select
  endif !field_name
enddo !End of i loop
if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine INTERPOLATOR_3D_
!
!#######################################################################
!
!++lwh
!---------------------------------------------------------------------
!> @brief interpolator_2D receives a field name as input and
!!        interpolates the field to model a 2D grid and time axis.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <Time> The model time that you wish to interpolate to
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_2D : You must call interpolator_init
!!                    before calling interpolator"
!! @throw FATAL "interpolator_2D 1:  file="
!! @throw FATAL "interpolator_2D 2:  file="
!! @throw FATAL "interpolator_2D 3:  file="
!! @throw FATAL "interpolator_2D 4:  file="
!! @throw FATAL "interpolator_2D 5:  file="
!! @throw FATAL "interpolator_2D : No data from the previous climatology
!!                    time but we have the next time. How did this happen?"
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_2D_(clim_type, Time, interp_data, field_name, is, js, clim_units)
!
! Return 2-D field interpolated to model grid and time
!
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init
!
! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   Time        : The model time that you wish to interpolate to.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)     :: field_name
type(time_type)       , intent(in)     :: Time
real(FMS_INTP_KIND_), dimension(:,:),   intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: taum, taup
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%FMS_INTP_TYPE_%levs(:)))
integer :: istart,iend,jstart,jend
logical :: result, found
logical :: found_field=.false.
integer :: modyear, modmonth, modday, modhour, modminute, modsecond
integer :: climyear, climmonth, climday, climhour, climminute, climsecond
integer :: year1, month1, day, hour, minute, second
integer :: climatology, m
type(time_type) :: t_prev, t_next
type(time_type), dimension(2) :: month
integer :: indexm, indexp, yearm, yearp
integer :: j, i, n
character(len=256) :: err_msg

integer, parameter :: lkind=FMS_INTP_KIND_

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_2D : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
!++lwh
  if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then
!--lwh

    found_field=.true.

    if(present(clim_units)) then
      call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
      clim_units = chomp(clim_units)
    endif

!----------------------------------------------------------------------
!   skip the time interpolation portion of this routine if subroutine
!   obtain_interpolator_time_slices has already been called on this
!   stewp for this interpolate_type variable.
!----------------------------------------------------------------------

if ( .not. clim_type%separate_time_vary_calc) then
!   print *, 'TIME INTERPOLATION NOT SEPARATED 2d--',  &
!                                   trim(clim_type%file_name), mpp_pe()
    if (clim_type%climatological_year) then
!++lwh
       if (size(clim_type%time_slice) > 1) then
          call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, &
                           taum, taup, modtime=YEAR, err_msg=err_msg )
          if(trim(err_msg) /= '') then
             call mpp_error(FATAL,'interpolator_2D 1: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
          endif
       else
          taum = 1
          taup = 1
          clim_type%FMS_INTP_TYPE_%tweight = 0._lkind
       end if
!--lwh
    else
       call time_interp(Time, clim_type%time_slice, clim_type%FMS_INTP_TYPE_%tweight, taum, taup, err_msg=err_msg )
       if(trim(err_msg) /= '') then
          call mpp_error(FATAL,'interpolator_2D 2: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
       endif
    endif

! If the climatology file has seasonal, a split time-line or has all the data
! read in then enter this loop.
!
    if(clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init) then
      clim_type%itaum=taum
      clim_type%itaup=taup
    endif

!    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
!      ! Check if delta-time is greater than delta of first two climatology time-slices.
!      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
!           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
!      ! The difference between the model time and the last climatology time-slice previous to the model time.
!      ! We need 2 time levels. Check we have the correct data.
!        itaum=0
!        itaup=0
!      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
!      ! the climatology year into the appropriate place.
!
!        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
!        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datem = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!        call time_interp(clim_datem, clim_type%time_slice, tweight1, taum1, taup1 )
!
!
!        call get_date(clim_type%time_slice(taup), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datep = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!
!
!      endif
!
!    endif
    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
      ! The difference between the model time and the last climatology time-slice previous to the model time.
      ! We need 2 time levels.
        clim_type%itaum=0
        clim_type%itaup=0
      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
      ! the climatology year into the appropriate place.


      ! We need to get the previous months data for the climatology year before
      ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:,:),2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1,m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m
        enddo
        do m = 1,12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m,climatology), year1, month1, day, hour, minute, second)
          if ( month1 == modmonth ) then
!RSHBUGFX   if ( modday <= day ) then
            if ( modday <  day ) then
              indexm = m-1 ; indexp = m
            else
              indexm = m ; indexp = m+1
            endif
          endif

        enddo
        if ( indexm == 0 ) then
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        endif
        call get_date(clim_type%time_slice(indexm+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if ( indexp == 13 ) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        endif
        call get_date(clim_type%time_slice(indexp+(climatology-1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

        call time_interp(Time, month, clim_type%FMS_INTP_TYPE_%tweight3, taum, taup, err_msg=err_msg ) ! tweight3 is
                                                                                 !! the time weight between the months.
        if (taum==2 .and. taup==2) &
             clim_type%FMS_INTP_TYPE_%tweight3 = 1._lkind ! protect against post-perth time_interp behavior
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_2D 3: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexm+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexm+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, clim_type%FMS_INTP_TYPE_%tweight1, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight1 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_2D 4: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif

        month(1) = clim_type%time_slice(indexp+(climatology-1)*12)
        month(2) = clim_type%time_slice(indexp+climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, clim_type%FMS_INTP_TYPE_%tweight2, taum, taup, err_msg=err_msg ) !tweight1 is
                                                                      !! the time weight between the climatology years.
        !> protect against post-perth time_interp behavior below
        if (taum==2 .and. taup==2) clim_type%FMS_INTP_TYPE_%tweight2 = 1._lkind
        if(trim(err_msg) /= '') then
           call mpp_error(FATAL,'interpolator_2D 5: '//trim(err_msg)//' file='//trim(clim_type%file_name), FATAL)
        endif


        if (indexm == clim_type%indexm(i) .and.  &
          indexp == clim_type%indexp(i) .and. &
          climatology == clim_type%climatology(i)) then
        else
          clim_type%indexm(i) = indexm
          clim_type%indexp(i) = indexp
          clim_type%climatology(i) = climatology
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i),  &
            clim_type%indexm(i)+(clim_type%climatology(i)-1)*12,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i),   &
            clim_type%indexp(i)+(clim_type%climatology(i)-1)*12,i,Time)
          call read_data(clim_type,clim_type%field_name(i),   &
            clim_type%FMS_INTP_TYPE_%pmon_nyear(:,:,:,i),  &
            clim_type%indexm(i)+clim_type%climatology(i)*12,i,Time)
          call read_data(clim_type,clim_type%field_name(i),  &
            clim_type%FMS_INTP_TYPE_%nmon_nyear(:,:,:,i),  &
            clim_type%indexp(i)+clim_type%climatology(i)*12,i,Time)
        endif




      else ! We are within a climatology data set

        if (taum /= clim_type%time_init(i,1) .or. &
            taup /= clim_type%time_init(i,2) ) then

          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%pmon_pyear(:,:,:,i), taum,i,Time)
! Read the data for the next month in the previous climatology.
          call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%nmon_pyear(:,:,:,i), taup,i,Time)
!RSHbug   clim_type%pmon_nyear = 0.0
!RSHbug   clim_type%nmon_nyear = 0.0

!         clim_type%pmon_nyear(:,:,:,i) = 0.0
!         clim_type%nmon_nyear(:,:,:,i) = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
          clim_type%indexm(i) = 0
          clim_type%indexp(i) = 0
          clim_type%climatology(i) = 0


          clim_type%time_init(i,1) = taum
          clim_type%time_init(i,2) = taup
        endif
!       clim_type%tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        clim_type%FMS_INTP_TYPE_%tweight1 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight2 = 0.0_lkind
        clim_type%FMS_INTP_TYPE_%tweight3 = clim_type%FMS_INTP_TYPE_%tweight
      endif

    endif ! (BILINEAR)

    if(clim_type%TIME_FLAG .eq. LINEAR .and. &
        (.not. read_all_on_init) ) then
! We need 2 time levels. Check we have the correct data.
      clim_type%itaum=0
      clim_type%itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(i,n) .eq. taum ) clim_type%itaum = n
        if (clim_type%time_init(i,n) .eq. taup ) clim_type%itaup = n
      enddo

      if (clim_type%itaum.eq.0 .and. clim_type%itaup.eq.0) then
      !Neither time is set so we need to read 2 time slices.
      !Set up
      ! field(:,:,:,1) as the previous time slice.
      ! field(:,:,:,2) as the next time slice.
        call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          clim_type%itaum = 1
        call read_data(clim_type,clim_type%field_name(i), clim_type%FMS_INTP_TYPE_%data5d(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          clim_type%itaup = 2
      endif ! clim_type%itaum.eq.clim_type%itaup.eq.0
      if (clim_type%itaum.eq.0 .and. clim_type%itaup.ne.0) then
      ! Can't think of a situation where we would have the next time level but not the previous.
        call mpp_error(FATAL,'interpolator_2D : No data from the previous climatology time but we have&
                            & the next time. How did this happen?')
      endif
      if (clim_type%itaum.ne.0 .and. clim_type%itaup.eq.0) then
      !We have the previous time step but not the next time step data
        clim_type%itaup = 1
        if (clim_type%itaum .eq. 1 ) clim_type%itaup = 2
        call read_data(clim_type,clim_type%field_name(i), &
                       clim_type%FMS_INTP_TYPE_%data5d(:,:,:,clim_type%itaup,i), taup,i, Time)
        clim_type%time_init(i,clim_type%itaup)=taup
      endif
    endif! TIME_FLAG .eq. LINEAR .and. (.not. read_all_on_init)

  endif ! (.not. separate_time_vary_calc)



select case(clim_type%TIME_FLAG)
  case (LINEAR)
    hinterp_data = (1._lkind-clim_type%FMS_INTP_TYPE_%tweight)*&
         clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,clim_type%itaum,i) &
         + clim_type%FMS_INTP_TYPE_%tweight*clim_type%FMS_INTP_TYPE_%data5d &
         (istart:iend,jstart:jend,:,clim_type%itaup,i)
! case (SEASONAL)
! Do sine fit to data at this point
  case (BILINEAR)
    hinterp_data = &
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight1) * &
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3) * &
    clim_type%FMS_INTP_TYPE_%pmon_pyear(istart:iend,jstart:jend,:,i) + &
    (1._lkind-clim_type%FMS_INTP_TYPE_%tweight2) * &
    clim_type%FMS_INTP_TYPE_%tweight3  * clim_type%FMS_INTP_TYPE_%nmon_pyear(istart:iend,jstart:jend,:,i) + &
    clim_type%FMS_INTP_TYPE_%tweight1 * (1._lkind-clim_type%FMS_INTP_TYPE_%tweight3) * &
    clim_type%FMS_INTP_TYPE_%pmon_nyear(istart:iend,jstart:jend,:,i) + &
    clim_type%FMS_INTP_TYPE_%tweight2 * clim_type%FMS_INTP_TYPE_%tweight3 * &
    clim_type%FMS_INTP_TYPE_%nmon_nyear(istart:iend,jstart:jend,:,i)

end select

found = .false.
do j = 1,size(climo_diag_name(:))
  if (trim(adjustl(lowercase(climo_diag_name(j)))) .eq. trim(adjustl(lowercase(clim_type%field_name(i))))) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),hinterp_data,Time,is_in=istart,js_in=jstart)
  endif
endif

  interp_data(:,:) = hinterp_data(:,:,1)

  endif !field_name
enddo !End of i loop

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine INTERPOLATOR_2D_
!--lwh
!
!#######################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_4D_no_time_axis receives a field name as input and
!!        interpolates the field to model a 4D grid.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <phalf> The half level model pressure field
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_4D_no_time_axis : You must call
!!                    interpolator_init before calling interpolator"
!! @throw FATAL "interpolator_mod: cannot use 4D interface to
!!                    interpolator for this file"
!! @throw NOTE "Interpolator: model surface pressure is greater than
!!                    surface pressure of input data for "
!! @throw NOTE "Interpolator: model top pressure is less than surface
!!                    pressure of input data for "
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_4D_NO_TIME_AXIS_(clim_type, phalf, interp_data, field_name, is,js, clim_units)

! Return 4-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of a field that you wish to interpolate.
!                 all variables within this interpolate_type variable
!                 will be interpolated on this call. field_name may
!                 be any one of the variables.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model fields
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
real(FMS_INTP_KIND_), dimension(:,:,:),   intent(in)  :: phalf
real(FMS_INTP_KIND_), dimension(:,:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: ilon
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),&
                        size(interp_data,2),size(clim_type%FMS_INTP_TYPE_%levs(:)),size(clim_type%field_name(:)))
real(FMS_INTP_KIND_) :: p_fact(size(interp_data,1),size(interp_data,2))
real(FMS_INTP_KIND_) :: pclim(size(clim_type%FMS_INTP_TYPE_%halflevs(:)))
integer :: istart,iend,jstart,jend
logical :: found_field=.false.
integer :: i, j, k, n

integer, parameter :: lkind=FMS_INTP_KIND_

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_4D_no_time_axis : You must call interpolator_init before calling interpolator")

do n=2,size(clim_type%field_name(:))
  if (clim_type%vert_interp(n) /= clim_type%vert_interp(n-1) .or. &
   clim_type%out_of_bounds(n) /= clim_type%out_of_bounds(n-1)) then
    if (mpp_pe() == mpp_root_pe() ) then
      print *, 'processing file ' // trim(clim_type%file_name)
    endif
    call mpp_error (FATAL, 'interpolator_mod: &
            &cannot use 4D interface to interpolator for this file')
  endif
end do

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then
    found_field=.true.
    exit
  endif
end do
i = 1

if(present(clim_units)) then
  call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
  clim_units = chomp(clim_units)
endif

do n=1, size(clim_type%field_name(:))
  hinterp_data(:,:,:,n) = clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,1,n)
end do

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0_lkind
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

    do i= 1, size(clim_type%field_name(:))
      select case(clim_type%mr(i))
      case(KG_M2)
        do k = 1,size(hinterp_data,3)
          hinterp_data(:,:,k,i) = hinterp_data(:,:,k,i)/&
               ((clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))*p_fact)
        enddo
      end select
    enddo

   i = 1

do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%FMS_INTP_TYPE_%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than surface pressure of input data for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than top pressure of input data for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,:),interp_data(ilon,j,:,:))
         case(INTERP_LINEAR_P)
          do n=1, size(clim_type%field_name(:))
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:,n),interp_data(ilon,j,:,n))
          end do
      end select
   enddo
enddo

     do i= 1, size(clim_type%field_name(:))

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k,i) = interp_data(:,:,k,i)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

     end do

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine INTERPOLATOR_4D_NO_TIME_AXIS_

!#######################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_3D_no_time_axis receives a field name as input and
!!        interpolates the field to model a 3D grid.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <phalf> The half level model pressure field
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_3D_no_time_axis : You must call
!!                    interpolator_init before calling interpolator"
!! @throw FATAL "interpolator_mod: cannot use 4D interface to
!!                    interpolator for this file"
!! @throw NOTE "Interpolator: model surface pressure is greater than
!!                    climatology surface pressure for "
!! @throw NOTE "Interpolator: model top pressure is less than
!!                    climatology top pressure for "
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_3D_NO_TIME_AXIS_(clim_type, phalf, interp_data, field_name, is,js, clim_units)

! Return 3-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
real(FMS_INTP_KIND_), dimension(:,:,:), intent(in)  :: phalf
real(FMS_INTP_KIND_), dimension(:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
integer :: ilon          !< No description
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),&
                        size(interp_data,2),size(clim_type%FMS_INTP_TYPE_%levs(:))) !< No description
real(FMS_INTP_KIND_) :: p_fact(size(interp_data,1),size(interp_data,2))   !< No description
real(FMS_INTP_KIND_) :: pclim(size(clim_type%FMS_INTP_TYPE_%halflevs(:))) !< No description
integer :: istart,iend,jstart,jend                                        !< No description
logical :: found_field=.false.          !< No description
integer :: i, j, k                   !< No description

integer, parameter :: lkind=FMS_INTP_KIND_

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_3D_no_time_axis : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then
    found_field=.true.
    if(present(clim_units)) then
      call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
      clim_units = chomp(clim_units)
    endif

    hinterp_data = clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,1,i)

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0_lkind
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       hinterp_data(:,:,k) = hinterp_data(:,:,k)/&
            ((clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))*p_fact)
    enddo
end select

do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%FMS_INTP_TYPE_%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         if (verbose > 3) then
         call mpp_error(NOTE,"Interpolator: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         endif
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
         case(INTERP_LINEAR_P)
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
      end select
   enddo
enddo

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k) = interp_data(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

  endif !field_name
enddo !End of i loop
if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif
end subroutine INTERPOLATOR_3D_NO_TIME_AXIS_

!#######################################################################
!
!---------------------------------------------------------------------
!> @brief interpolator_2D_no_time_axis receives a field name as input and
!!        interpolates the field to model a 2D grid.
!!
!! @param [inout] <clim_type> The interpolate type previously defined by a call to interpolator_init
!! @param [in] <field_name> The name of a field that you wish to interpolate
!! @param [in] <is> OPTIONAL: Index for the physics window
!! @param [in] <js> OPTIONAL: Index for the physics window
!! @param [out] <interp_data> The model fields with the interpolated climatology data
!! @param [out] <clim_units> OPTIONAL: The units of field_name
!!
!! @throw FATAL "interpolator_2D_no_time_axis : You must call
!!                    interpolator_init before calling interpolator"
!! @throw FATAL "Interpolator: the field name is not contained in this
!!                    intepolate_type: "
subroutine INTERPOLATOR_2D_NO_TIME_AXIS_(clim_type, interp_data, field_name, is, js, clim_units)

! Return 2-D field interpolated to model grid

! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolator_init

! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   is, js      : The indices of the physics window.

! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name

type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)     :: field_name
real(FMS_INTP_KIND_), dimension(:,:),   intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
real(FMS_INTP_KIND_) :: hinterp_data(size(interp_data,1),&
                                     size(interp_data,2),size(clim_type%FMS_INTP_TYPE_%levs(:)))
integer :: istart,iend,jstart,jend
logical :: found_field=.false.
integer :: i

if (.not. module_is_initialized .or. .not. allocated(clim_type%FMS_INTP_TYPE_%lon)) &
   call mpp_error(FATAL, "interpolator_2D_no_time_axis : You must call interpolator_init before calling interpolator")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name(:))
  if ( trim(adjustl(lowercase(field_name))) == trim(adjustl(lowercase(clim_type%field_name(i)))) ) then

    found_field=.true.

    if(present(clim_units)) then
      call get_variable_units(clim_type%fileobj, clim_type%field_name(i), clim_units)
      clim_units = chomp(clim_units)
    endif

    hinterp_data = clim_type%FMS_INTP_TYPE_%data5d(istart:iend,jstart:jend,:,1,i)

    interp_data(:,:) = hinterp_data(:,:,1)

  endif !field_name
enddo !End of i loop

if( .not. found_field) then !field name is not in interpolator file.ERROR.
  call mpp_error(FATAL,"Interpolator: the field name is not contained in this &
                   &intepolate_type: "//trim(field_name))
endif

end subroutine INTERPOLATOR_2D_NO_TIME_AXIS_

!#######################################################################
!
!---------------------------------------------------------------------
!> @brief read_data receives various climate data as inputs and
!!        returns a horizontally interpolated climatology field.
!!
!! @param [in] <clim_type> The interpolate type which contains the data
!! @param [in] <src_field> The field type
!! @param [in] <nt> The index of the time slice of the climatology that you wish to read
!! @param [in] <i> OPTIONAL: The index of the field name that you are trying to read
!! @param [in] <Time> OPTIONAL: The model time. Used for diagnostic purposes only
!! @param [out] <hdata> The horizontally interpolated climatology field. This
!                       field will still be on the climatology vertical grid
subroutine READ_DATA_(clim_type,field_name, hdata, nt,i, Time)
!
!  INTENT IN
!    clim_type : The interpolate type which contains the data
!    src_field : The field type
!    nt        : The index of the time slice of the climatology that you wish to read.
!    i         : The index of the field name that you are trying to read. (optional)
!    Time      : The model time. Used for diagnostic purposes only. (optional)
!
!  INTENT OUT
!
!    hdata     : The horizontally interpolated climatology field. This
!                field will still be on the climatology vertical grid.
!
type(interpolate_type)   , intent(in)  :: clim_type
character(len=*)         , intent(in)  :: field_name
integer                  , intent(in)  :: nt
real(FMS_INTP_KIND_)     , intent(out) :: hdata(:,:,:)
integer                  , intent(in)  :: i
type(time_type), optional, intent(in)  :: Time

integer   :: k, km
! sjs
real(FMS_INTP_KIND_), allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%FMS_INTP_TYPE_%lon(:)),size(clim_type%FMS_INTP_TYPE_%lat(:)), &
                        size(clim_type%FMS_INTP_TYPE_%levs(:))))
      if(clim_type%has_level(i)) then ! has vertical level
          call fms2_io_read_data(clim_type%fileobj,field_name, climdata,nt)
      else  ! no vertical level
          call fms2_io_read_data(clim_type%fileobj,field_name, climdata(:,:,1),nt)
      endif

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%FMS_INTP_TYPE_%lon(:)),   &
                           size(clim_type%FMS_INTP_TYPE_%lat(:)), &
                           size(clim_type%FMS_INTP_TYPE_%levs(:))))
        km = size(clim_type%FMS_INTP_TYPE_%levs(:))
        do k=1, km
          climdata2(:,:,k) = climdata(:,:,km+1-k)
        end do
        climdata = climdata2
        deallocate (climdata2)
      endif
      call horiz_interp(clim_type%interph, climdata, hdata)
      if (clim_diag_initialized) &
        call diag_read_data(clim_type,climdata,i, Time)
      deallocate(climdata)


 end subroutine READ_DATA_

!#######################################################################
!
!---------------------------------------------------------------------
!> @brief read_data_no_time_axis receives various climate data as inputs and
!!        returns a horizontally interpolated climatology field without the
!!            time axis.
!!
!! @param [in] <clim_type> The interpolate type which contains the data
!! @param [in] <src_field> The field type
!! @param [in] <i> OPTIONAL: The index of the field name that you are trying to read
!! @param [out] <hdata> The horizontally interpolated climatology field. This
!                       field will still be on the climatology vertical grid
subroutine READ_DATA_NO_TIME_AXIS_(clim_type,field_name, hdata, i)

!  INTENT IN
!    clim_type : The interpolate type which contains the data
!    src_field : The field type
!    i         : The index of the field name that you are trying to read. (optional)

!  INTENT OUT

!    hdata     : The horizontally interpolated climatology field. This
!                field will still be on the climatology vertical grid.

type(interpolate_type)   , intent(in)  :: clim_type
character(len=*)         , intent(in)  :: field_name
real(FMS_INTP_KIND_)     , intent(out) :: hdata(:,:,:)
integer        , intent(in)  :: i

integer   :: k, km
! sjs
real(FMS_INTP_KIND_), allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%FMS_INTP_TYPE_%lon(:)),&
                        size(clim_type%FMS_INTP_TYPE_%lat(:)),&
                        size(clim_type%FMS_INTP_TYPE_%levs(:))))

      if(clim_type%has_level(i)) then ! has vertical level
         call fms2_io_read_data(clim_type%fileobj,field_name, climdata)
      else  ! no vertical level
         call fms2_io_read_data(clim_type%fileobj,field_name, climdata(:,:,1))
      endif
!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%FMS_INTP_TYPE_%lon(:)),   &
                           size(clim_type%FMS_INTP_TYPE_%lat(:)), &
                           size(clim_type%FMS_INTP_TYPE_%levs(:))))
        km = size(clim_type%FMS_INTP_TYPE_%levs(:))
        do k=1, km
          climdata2(:,:,k) = climdata(:,:,km+1-k)
        end do
        climdata = climdata2
        deallocate (climdata2)
      endif

      call horiz_interp(clim_type%interph, climdata, hdata)
      deallocate(climdata)

end subroutine READ_DATA_NO_TIME_AXIS_

!#######################################################################
!
!---------------------------------------------------------------------
!> @brief diag_read_data receives the data read in by read_data as
!!            inputs and runs a diagnosis.
!!
!! @param [in] <clim_type> The interpolate type which contains the data
!! @param [in] <model_data> The data read in from file that is being diagnosed
!! @param [in] <i> The index of the field name that you are diagnosing
!! @param [in] <Time> The model time.
subroutine DIAG_READ_DATA_(clim_type,model_data, i, Time)
!
! A routine to diagnose the data read in by read_data
!
!  INTENT IN
!    clim_type  : The interpolate type.
!    model_data : The data read in from file that is being diagnosed.
!    i          : The index of the field name that you are diagnosing.
!    Time       : The model time
!
type(interpolate_type), intent(in) :: clim_type
real(FMS_INTP_KIND_)  , intent(in) :: model_data(:,:,:)
integer               , intent(in) :: i
type(time_type)       , intent(in) :: Time

integer :: j,k
real(FMS_INTP_KIND_) :: col_data(size(model_data,1),size(model_data,2))
logical :: result, found

integer, parameter :: lkind=FMS_INTP_KIND_

found = .false.
do j = 1,size(climo_diag_name(:))
  if (trim(adjustl(lowercase(climo_diag_name(j)))) .eq. trim(adjustl(lowercase(clim_type%field_name(i))))) then
      found = .true.
      exit
  endif
enddo

if(found) then
  if(climo_diag_id(j)>0) then
  col_data(:,:)=0.0_lkind
    do k=1,size(model_data,3)
      col_data(:,:) = col_data(:,:) + &
        model_data(:,:,k)* &
        (clim_type%FMS_INTP_TYPE_%halflevs(k+1)-clim_type%FMS_INTP_TYPE_%halflevs(k))/real(grav,FMS_INTP_KIND_)
    enddo
    result = send_data(climo_diag_id(j),col_data(clim_type%is:clim_type%ie,clim_type%js:clim_type%je),Time)
  endif
endif

end subroutine DIAG_READ_DATA_
!
!#################################################################
!
!
!---------------------------------------------------------------------
!> @brief interp_weighted_scalar_2D receives the variables grdin,
!!            grdout, and datin as inputs and returns datout.
!!
!! @param [in] <grdin> No description
!! @param [in] <grdout> No description
!! @param [in] <datin> No description
!! @param [out] <datout> No description
subroutine INTERP_WEIGHTED_SCALAR_2D_ (grdin, grdout, datin, datout )
real(FMS_INTP_KIND_), intent(in),  dimension(:) :: grdin, grdout
real(FMS_INTP_KIND_), intent(in),  dimension(:,:) :: datin
real(FMS_INTP_KIND_), intent(out), dimension(:,:) :: datout

integer :: j, k, n
integer, parameter :: lkind=FMS_INTP_KIND_

if (size(grdin(:)).ne. (size(datin,1)+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout,1 )+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : output data and pressure do not have the same number of levels')

  do k = 1, size(datout,1 )
   datout(k,:) = 0.0_lkind

     do j = 1, size(datin,1 )

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdout(k))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdin(j))
          end do

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

          do n= 1, size(datin,2)
           datout(k,n) = datout(k,n) + datin(j,n)*(grdin(j+1)-grdin(j))
          end do

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

          do n= 1, size(datin,2)
          datout(k,n) = datout(k,n) + datin(j,n)*(grdout(k+1)-grdout(k))

          end do
        endif

     enddo

     do n= 1, size(datin,2)
       datout(k,n) = datout(k,n)/(grdout(k+1)-grdout(k))
     end do

  enddo

end subroutine INTERP_WEIGHTED_SCALAR_2D_

!---------------------------------------------------------------------
!> @brief interp_weighted_scalar_1D receives the variables grdin,
!!        grdout, and datin as inputs and returns datout.
!!
!! @param [in] <grdin> No description
!! @param [in] <grdout> No description
!! @param [in] <datin> No description
!! @param [out] <datout> No description
subroutine INTERP_WEIGHTED_SCALAR_1D_ (grdin, grdout, datin, datout )
real(FMS_INTP_KIND_), intent(in),  dimension(:) :: grdin, grdout, datin
real(FMS_INTP_KIND_), intent(out), dimension(:) :: datout

integer :: j, k
integer, parameter :: lkind=FMS_INTP_KIND_

if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call  mpp_error(FATAL,'interp_weighted_scalar : output data and pressure do not have the same number of levels')

  do k = 1, size(datout(:))
   datout(k) = 0.0_lkind

     do j = 1, size(datin(:))

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdout(k))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdin(j))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdin(j))

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdout(k))

        endif

     enddo

     datout(k) = datout(k)/(grdout(k+1)-grdout(k))

  enddo

end subroutine INTERP_WEIGHTED_SCALAR_1D_
!
!#################################################################
!
!---------------------------------------------------------------------
!> @brief interp_linear receives the variables grdin,
!!            grdout, and datin as inputs and returns a linear
!!            interpolation.
!!
!! @param [in] <grdin> No description
!! @param [in] <grdout> No description
!! @param [in] <datin> No description
!! @param [out] <datout> No description
subroutine INTERP_LINEAR_ ( grdin, grdout, datin, datout )
real(FMS_INTP_KIND_), intent(in),  dimension(:) :: grdin, grdout, datin
real(FMS_INTP_KIND_), intent(out), dimension(:) :: datout

integer :: j, k, n
real(FMS_INTP_KIND_) :: wt

integer, parameter :: lkind=FMS_INTP_KIND_

if (size(grdin(:)).ne. (size(datin(:))+1)) &
 call mpp_error(FATAL,'interp_linear : input data and pressure do not have the same number of levels')
if (size(grdout(:)).ne. (size(datout(:))+1)) &
 call mpp_error(FATAL,'interp_linear : output data and pressure do not have the same number of levels')


  n = size(grdin(:))

  do k= 1, size(datout(:))

   ! ascending grid values
     if (grdin(1) < grdin(n)) then
         do j = 2, size(grdin(:))-1
           if (grdout(k) <= grdin(j)) exit
         enddo
   ! descending grid values
     else
         do j = size(grdin(:)), 3, -1
           if (grdout(k) <= grdin(j-1)) exit
         enddo
     endif

   ! linear interpolation
     wt = (grdout(k)-grdin(j-1)) / (grdin(j)-grdin(j-1))
!print '(a,2i3,4f6.1)', 'k,j=',k,j,grdout(k),grdin(j-1),grdin(j),wt
   ! constant value extrapolation
   ! wt = min(max(wt,0.),1.)

     datout(k) = (1._lkind-wt)*datin(j-1) + wt*datin(j)

  enddo

end subroutine INTERP_LINEAR_
!
!########################################################################

!> @}
! close documentation grouping
!
!#######################################################################
