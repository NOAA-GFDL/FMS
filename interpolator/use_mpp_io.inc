subroutine fmsio_interpolator_init(clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units, single_year_file)

type(interpolate_type), intent(inout) :: clim_type
character(len=*), intent(in)            :: file_name
real            , intent(in)            :: lonb_mod(:,:), latb_mod(:,:)
character(len=*), intent(in) , optional :: data_names(:)
!++lwh
integer         , intent(in)            :: data_out_of_bounds(:)
integer         , intent(in), optional  :: vert_interp(:)
!--lwh
character(len=*), intent(out), optional :: clim_units(:)
logical,          intent(out), optional :: single_year_file

integer                      :: unit
character(len=64)            :: src_file
!++lwh
real                         :: dlat, dlon
!--lwh
type(time_type)              :: base_time
logical                      :: NAME_PRESENT
real                         :: dtr,tpi
integer                      :: fileday, filemon, fileyr, filehr, filemin,filesec, m,m1
character(len= 20)           :: fileunits
real, dimension(:), allocatable  :: alpha
integer   :: j, i
logical :: non_monthly
character(len=24) :: file_calendar
character(len=256) :: error_mesg
integer :: model_calendar
integer :: yr, mo, dy, hr, mn, sc
integer :: n
type(time_type) :: Julian_time, Noleap_time
real, allocatable :: time_in(:)
real, allocatable, save :: agrid_mod(:,:,:)
integer :: nx, ny

clim_type%separate_time_vary_calc = .false.

tpi = 2.0*PI ! 4.*acos(0.)
dtr = tpi/360.

num_fields = 0

!--------------------------------------------------------------------
! open source file containing fields to be interpolated
!--------------------------------------------------------------------
src_file = 'INPUT/'//trim(file_name)

if(fms_io_file_exist(trim(src_file))) then
   call mpp_open( unit, trim(src_file), action=MPP_RDONLY, &
                  form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
else
!Climatology file doesn't exist, so exit
   call mpp_error(FATAL,'Interpolator_init : Data file '//trim(src_file)//' does not exist')
endif

!Find the number of variables (nvar) in this file
call mpp_get_info(unit, ndim, nvar, natt, ntime)
clim_type%unit      = unit
clim_type%file_name = trim(file_name)

num_fields = nvar
if(present(data_names)) num_fields= size(data_names(:))

! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------
allocate(axes(ndim))
call mpp_get_axes(unit, axes, time_axis)

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

do i = 1, ndim
  call mpp_get_atts(axes(i), name=name,len=len,units=units,  &
                    calendar=file_calendar, sense=sense)
  select case(name)
    case('lat')
      nlat=len
      allocate(clim_type%lat(nlat))
      call mpp_get_axis_data(axes(i),clim_type%lat)
      select case(units(1:6))
        case('degree')
          clim_type%lat = clim_type%lat*dtr
        case('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lat not recognised in file "//file_name)
      end select
    case('lon')
      nlon=len
      allocate(clim_type%lon(nlon))
      call mpp_get_axis_data(axes(i),clim_type%lon)
      select case(units(1:6))
        case('degree')
          clim_type%lon = clim_type%lon*dtr
        case('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lon not recognised in file "//file_name)
      end select
    case('latb')
      nlatb=len
      allocate(clim_type%latb(nlatb))
      call mpp_get_axis_data(axes(i),clim_type%latb)
      select case(units(1:6))
        case('degree')
          clim_type%latb = clim_type%latb*dtr
        case('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for latb not recognised in file "//file_name)
      end select
    case('lonb')
      nlonb=len
      allocate(clim_type%lonb(nlonb))
      call mpp_get_axis_data(axes(i),clim_type%lonb)
      select case(units(1:6))
        case('degree')
          clim_type%lonb = clim_type%lonb*dtr
        case('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lonb not recognised in file "//file_name)
      end select
    case('pfull')
      nlev=len
      allocate(clim_type%levs(nlev))
      call mpp_get_axis_data(axes(i),clim_type%levs)
      clim_type%level_type = PRESSURE
  ! Convert to Pa
      if( trim(adjustl(lowercase(chomp(units)))) == "mb" .or. trim(adjustl(lowercase(chomp(units)))) == "hpa") then
         clim_type%levs = clim_type%levs * 100.
      end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
      if( sense == 1 ) then
        clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlev))
          do n = 1, nlev
          alpha(n) = clim_type%levs(nlev-n+1)
          end do
          do n = 1, nlev
          clim_type%levs(n) = alpha(n)
          end do
          deallocate (alpha)
      else
        clim_type%vertical_indices = INCREASING_DOWNWARD
      endif

    case('phalf')
      nlevh=len
      allocate(clim_type%halflevs(nlevh))
      call mpp_get_axis_data(axes(i),clim_type%halflevs)
      clim_type%level_type = PRESSURE
  ! Convert to Pa
      if( trim(adjustl(lowercase(chomp(units)))) == "mb" .or. trim(adjustl(lowercase(chomp(units)))) == "hpa") then
         clim_type%halflevs = clim_type%halflevs * 100.
      end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
      if( sense == 1 ) then
        clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlevh))
          do n = 1, nlevh
          alpha(n) = clim_type%halflevs(nlevh-n+1)
          end do
          do n = 1, nlevh
          clim_type%halflevs(n) = alpha(n)
          end do
          deallocate (alpha)
      else
        clim_type%vertical_indices = INCREASING_DOWNWARD
      endif
    case('sigma_full')
      nlev=len
      allocate(clim_type%levs(nlev))
      call mpp_get_axis_data(axes(i),clim_type%levs)
      clim_type%level_type = SIGMA
    case('sigma_half')
      nlevh=len
      allocate(clim_type%halflevs(nlevh))
      call mpp_get_axis_data(axes(i),clim_type%halflevs)
      clim_type%level_type = SIGMA

    case('time')
      model_calendar = get_calendar_type()
      fileday = 0
      filemon = 0
      fileyr = 0
      filehr = 0
      filemin= 0
      filesec = 0
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
        case default
          call mpp_error(FATAL,'Interpolator_init : Time units not recognised in file '//file_name)
      end select

       clim_type%climatological_year = (fileyr == 0)

      if (.not. clim_type%climatological_year) then

!----------------------------------------------------------------------
!    if file date has a non-zero year in the base time, determine that
!    base_time based on the netcdf info.
!----------------------------------------------------------------------
        if ( (model_calendar == JULIAN .and.   &
             & trim(adjustl(lowercase(file_calendar))) == 'julian')  .or. &
             & (model_calendar == NOLEAP .and.   &
             & trim(adjustl(lowercase(file_calendar))) == 'noleap') )  then
          call mpp_error (NOTE, 'interpolator[1]_mod: Model and file&
                    & calendars are the same for file ' //   &
                    & trim(file_name) // '; no calendar conversion  &
                    &needed')
          base_time = set_date (fileyr, filemon, fileday, filehr, &
                                filemin,filesec)
        else if ( (model_calendar == JULIAN .and.   &
             & trim(adjustl(lowercase(file_calendar))) == 'noleap')) then
          call mpp_error (NOTE, 'interpolator[1]_mod: Using julian &
                            &model calendar and noleap file calendar&
                            & for file ' // trim(file_name) //   &
                            &'; calendar conversion needed')
          base_time = set_date_no_leap (fileyr, filemon, fileday,  &
               & filehr, filemin, filesec)
        else if ( (model_calendar == NOLEAP .and.   &
             & trim(adjustl(lowercase(file_calendar))) == 'julian')) then
          call mpp_error (NOTE, 'interpolator[1]_mod: Using noleap &
                            &model calendar and julian file calendar&
                            & for file ' // trim(file_name) //  &
                            &'; calendar conversion needed')
          base_time = set_date_julian (fileyr, filemon, fileday,  &
               & filehr, filemin, filesec)
        else
          call mpp_error (FATAL , 'interpolator[1]_mod: Model and file&
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
        time_in = 0.0
        clim_type%time_slice = set_time(0,0) + base_time
        clim_type%clim_times = set_time(0,0) + base_time
        call mpp_get_times(clim_type%unit, time_in)
        ntime_in = ntime
! determine whether the data is a continuous set of monthly values or
! a series of annual cycles spread throughout the period of data
        non_monthly = .false.
        do n = 1, ntime-1
!  Assume that the times in the data file correspond to days only.
          if (time_in(n+1) > (time_in(n) + 32.)) then
            non_monthly = .true.
            exit
          endif
        end do
        if (clim_type%climatological_year) then
          call mpp_error (NOTE, 'interpolator[1]_mod :'  // &
          trim(file_name) // ' is a year-independent climatology file')
        else
          call mpp_error (NOTE, 'interpolator[1]_mod :' // &
            trim(file_name) // ' is a timeseries file')
        endif

        do n = 1, ntime
!Assume that the times in the data file correspond to days only.


          if (clim_type%climatological_year) then
!! RSH NOTE:
!! for this case, do not add base_time. time_slice will be sent to
!! time_interp_list with the optional argument modtime=YEAR, so that
!! the time that is needed in time_slice is the displacement into the
!! year, not the displacement from a base_time.
            clim_type%time_slice(n) = &
                set_time(INT( ( time_in(n) - INT(time_in(n)) ) * SECONDS_PER_DAY ), &
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
            if ( (model_calendar == JULIAN .and.   &
                 & trim(adjustl(lowercase(file_calendar))) == 'julian')  .or. &
                 & (model_calendar == NOLEAP .and.   &
                 & trim(adjustl(lowercase(file_calendar))) == 'noleap') )  then

!---------------------------------------------------------------------
!    no calendar conversion needed.
!---------------------------------------------------------------------
              clim_type%time_slice(n) = &
                 set_time(INT( ( time_in(n) - INT(time_in(n)) ) * SECONDS_PER_DAY ),&
                                 INT(time_in(n)))  &
                  + base_time

!---------------------------------------------------------------------
!    convert file times from noleap to julian.
!---------------------------------------------------------------------
            else if ( (model_calendar == JULIAN .and.   &
                 & trim(adjustl(lowercase(file_calendar))) == 'noleap')) then
              Noleap_time = set_time (0, INT(time_in(n))) + base_time
              call get_date_no_leap (Noleap_time, yr, mo, dy, hr,  &
                                     mn, sc)
              clim_type%time_slice(n) = set_date_julian (yr, mo, dy,  &
                                                         hr, mn, sc)
              if (n == 1) then
                call print_date (clim_type%time_slice(1), &
                        str= 'for file ' // trim(file_name) // ', the &
                              &first time slice is mapped to :')
              endif
              if (n == ntime) then
                call print_date (clim_type%time_slice(ntime), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &last time slice is mapped to:')
              endif


!---------------------------------------------------------------------
!    convert file times from julian to noleap.
!---------------------------------------------------------------------
            else if ( (model_calendar == NOLEAP .and.   &
                 & trim(adjustl(lowercase(file_calendar))) == 'julian')) then
              Julian_time = set_time (0, INT(time_in(n))) + base_time
              call get_date_julian (Julian_time, yr, mo, dy, hr, mn, sc)
              clim_type%time_slice(n) = set_date_no_leap (yr, mo, dy, &
                                                          hr, mn, sc)
              if (n == 1) then
                call print_date (clim_type%time_slice(1), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &first time slice is mapped to :')
              endif
              if (n == ntime) then
                call print_date (clim_type%time_slice(ntime), &
                         str= 'for file ' // trim(file_name) // ', the &
                               &last time slice is mapped to:')
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
        time_in = 0.0
        clim_type%time_slice = set_time(0,0) + base_time
        clim_type%clim_times(1,1) = set_time(0,0) + base_time
      endif
      deallocate(time_in)
  end select ! case(name)
enddo


! -------------------------------------------------------------------
! For 2-D fields, allocate levs and halflevs here
!  code is still needed for case when only halflevs are in data file.
! -------------------------------------------------------------------
    if( .not. associated(clim_type%levs) ) then
        allocate( clim_type%levs(nlev) )
        clim_type%levs = 0.0
    endif
    if( .not. associated(clim_type%halflevs) )  then
        allocate( clim_type%halflevs(nlev+1) )
        clim_type%halflevs(1) = 0.0
        if (clim_type%level_type == PRESSURE) then
          clim_type%halflevs(nlev+1) = 1013.25* 100.0   ! MKS
        else if (clim_type%level_type == SIGMA   ) then
          clim_type%halflevs(nlev+1) = 1.0
        endif
        do n=2,nlev
           clim_type%halflevs(n) = 0.5*(clim_type%levs(n) + &
                                         clim_type%levs(n-1))
        end do
    endif
deallocate(axes)


! In the case where only the midpoints of the longitudes are defined we force the definition
! of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%lon) .and. .not. associated(clim_type%lonb)) &
   call mpp_error(FATAL,'Interpolator_init : There appears to be no longitude axis in file '//file_name)

if (.not. associated(clim_type%lonb) ) then

  if (size(clim_type%lon(:)) /= 1) then
    allocate(clim_type%lonb(size(clim_type%lon(:))+1))
    dlon = (clim_type%lon(2)-clim_type%lon(1))/2.0
    clim_type%lonb(1) = clim_type%lon(1) - dlon
    clim_type%lonb(2:) = clim_type%lon(1:) + dlon
  else

!! this is the case for zonal mean data, lon = 1, lonb not present
!! in file.

    allocate(clim_type%lonb(2))
    clim_type%lonb(1) = -360.*dtr
    clim_type%lonb(2) = 360.0*dtr
    clim_type%lon(1) = 0.0
  endif
endif

!clim_type%lonb=clim_type%lonb*dtr
! This assumes the lonb are in degrees in the NetCDF file!

if (.not. associated(clim_type%lat) .and. .not. associated(clim_type%latb)) &
   call mpp_error(FATAL,'Interpolator_init : There appears to be no latitude axis in file '//file_name)
! In the case where only the grid midpoints of the latitudes are defined we force the
! definition of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%latb) ) then
   allocate(clim_type%latb(nlat+1))
   dlat = (clim_type%lat(2)-clim_type%lat(1)) * 0.5
!  clim_type%latb(1) = min( 90., max(-90., clim_type%lat(1) - dlat) )
   clim_type%latb(1) = min( PI/2., max(-PI/2., clim_type%lat(1) - dlat) )
   clim_type%latb(2:nlat) = ( clim_type%lat(1:nlat-1) + clim_type%lat(2:nlat) ) * 0.5
   dlat = ( clim_type%lat(nlat) - clim_type%lat(nlat-1) ) * 0.5
!  clim_type%latb(nlat+1) = min( 90., max(-90., clim_type%lat(nlat) + dlat) )
   clim_type%latb(nlat+1) = min( PI/2., max(-PI/2., clim_type%lat(nlat) + dlat) )
endif
!clim_type%latb=clim_type%latb*dtr

!Assume that the horizontal interpolation within a file is the same for each variable.

 if (conservative_interp) then
    call horiz_interp_new (clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
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
          call cell_center2((/lonb_mod(i,j),latb_mod(i,j)/), &
               (/lonb_mod(i+1,j),latb_mod(i+1,j)/), &
               (/lonb_mod(i,j+1),latb_mod(i,j+1)/), &
               (/lonb_mod(i+1,j+1),latb_mod(i+1,j+1)/),  agrid_mod(i,j,:))
       enddo
       enddo
    endif

    !!! END DEBUG CODE

    call horiz_interp_new (clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
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
!   allocate(clim_type%data(size(lonb_mod(:))-1, size(latb_mod(:))-1, nlev, 2, num_fields))
   allocate(clim_type%pmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%pmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%nmon_nyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   allocate(clim_type%nmon_pyear(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, num_fields))
   clim_type%pmon_pyear = 0.0
   clim_type%pmon_nyear = 0.0
   clim_type%nmon_nyear = 0.0
   clim_type%nmon_pyear = 0.0
   clim_type%TIME_FLAG = BILINEAR
else
! We have a continuous time-line so treat as for 5-12 timelevels as below.
   if ( .not. read_all_on_init) then
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
               ntime, num_fields))
   endif
   clim_type%data = 0.0
   clim_type%TIME_FLAG = LINEAR
endif


!++lwh
 case (1:12)
!--lwh
! We have more than 4 timelevels
! Assume we have monthly or higher time resolution datasets (climatology or time series)
! So we only need to read 2 datasets and apply linear temporal interpolation.
   if ( .not. read_all_on_init) then
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 2, num_fields))
   else
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, &
               ntime, num_fields))
   endif
   clim_type%data = 0.0
   clim_type%TIME_FLAG = LINEAR
!++lwh
!case (1:4)
! Assume we have seasonal data and read in all the data.
! We can apply sine curves to these data.

!  allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, ntime, num_fields))
!  clim_type%data = 0.0
!  clim_type%TIME_FLAG = SEASONAL
!--lwh
! case (default)
 case(:0)
   clim_type%TIME_FLAG = NOTIME
   allocate(clim_type%data(size(lonb_mod,1)-1, size(latb_mod,2)-1, nlev, 1, num_fields))
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


allocate(clim_type%field_name(num_fields))
allocate(clim_type%field_type(num_fields))
allocate(clim_type%mr(num_fields))
allocate(clim_type%out_of_bounds(num_fields))
clim_type%out_of_bounds(:)=0
allocate(clim_type%vert_interp(num_fields))
clim_type%vert_interp(:)=0
!--------------------------------------------------------------------
!Allocate the space for the fields within the climatology data file.
allocate(varfields(nvar))
!--------------------------------------------------------------------
! Get the variable names out of the file.
call mpp_get_fields(clim_type%unit, varfields)

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
      NAME_PRESENT = .FALSE.
      do i=1,nvar
         call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if( trim(adjustl(lowercase(name))) == trim(adjustl(lowercase(data_names(j)))) ) then
            units=chomp(units)
            if (mpp_pe() == 0 ) write(*,*) 'Initializing src field : ',trim(name)
            clim_type%field_name(j) = name
            clim_type%field_type(j) = varfields(i)
            clim_type%mr(j)         = check_climo_units(units)
            NAME_PRESENT = .TRUE.
            if (present(clim_units)) clim_units(j) = units
            clim_type%out_of_bounds(j) = data_out_of_bounds( MIN(j,SIZE(data_out_of_bounds(:))) )
            if( clim_type%out_of_bounds(j) /= CONSTANT .and. &
                clim_type%out_of_bounds(j) /= ZERO ) &
               call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                                    & set to ZERO or CONSTANT")
            if( present(vert_interp) ) then
               clim_type%vert_interp(j) = vert_interp( MIN(j,SIZE(vert_interp(:))) )
               if( clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. &
                   clim_type%vert_interp(j) /= INTERP_LINEAR_P ) &
                  call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                       & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
            else
               clim_type%vert_interp(j) = INTERP_WEIGHTED_P
            end if
         endif
      enddo
      if(.not. NAME_PRESENT) &
         call mpp_error(FATAL,'interpolator_init : Check names of fields being passed. ' &
                              //trim(data_names(j))//' does not exist.')
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
   do i=1,nvar
      call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if (mpp_pe() ==0 ) write(*,*) 'Initializing src field : ',trim(name)
         clim_type%field_name(i) = lowercase(trim(name))
         clim_type%field_type(i) = varfields(i)
         clim_type%mr(i)         = check_climo_units(units)
         if (present(clim_units)) clim_units(i) = units
         clim_type%out_of_bounds(i) = data_out_of_bounds( MIN(i,SIZE(data_out_of_bounds(:))) )
         if( clim_type%out_of_bounds(i) /= CONSTANT .and. &
             clim_type%out_of_bounds(i) /= ZERO ) &
            call mpp_error(FATAL,"Interpolator_init: data_out_of_bounds must be&
                                 & set to ZERO or CONSTANT")
         if( present(vert_interp) ) then
            clim_type%vert_interp(i) = vert_interp( MIN(i,SIZE(vert_interp(:))) )
            if( clim_type%vert_interp(i) /= INTERP_WEIGHTED_P .and. &
                clim_type%vert_interp(i) /= INTERP_LINEAR_P ) &
               call mpp_error(FATAL,"Interpolator_init: vert_interp must be&
                                    & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
         else
            clim_type%vert_interp(i) = INTERP_WEIGHTED_P
         end if
   end do
!--lwh
endif

deallocate(varfields)


if( clim_type%TIME_FLAG .eq. SEASONAL ) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call interp_read_data_io( clim_type, clim_type%field_type(i), &
                         clim_type%data(:,:,:,n,i), n, i, base_time )
      enddo
   enddo
endif

if( clim_type%TIME_FLAG .eq. LINEAR  .and. read_all_on_init) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call interp_read_data_io( clim_type, clim_type%field_type(i), &
                         clim_type%data(:,:,:,n,i), n, i, base_time )
      enddo
   enddo

   call mpp_close (unit)
endif

if( clim_type%TIME_FLAG .eq. NOTIME ) then
! Read all the data at this point.
   do i=1,num_fields
     call interp_read_data_io_no_time_axis( clim_type, clim_type%field_type(i), &
                                  clim_type%data(:,:,:,1,i), i )
   enddo
   call mpp_close (unit)
endif

if (present (single_year_file)) then
  single_year_file = clim_type%climatological_year
endif

end subroutine fmsio_interpolator_init

!> \brief interp_read_data_io receives various climate data as inputs and
!!        returns a horizontally interpolated climatology field.
!!
!! \param [in] <clim_type> The interpolate type which contains the data
!! \param [in] <src_field> The field type
!! \param [in] <nt> The index of the time slice of the climatology that you wish to read
!! \param [in] <i> OPTIONAL: The index of the field name that you are trying to read
!! \param [in] <Time> OPTIONAL: The model time. Used for diagnostic purposes only
!! \param [out] <hdata> The horizontally interpolated climatology field. This
!                       field will still be on the climatology vertical grid
subroutine interp_read_data_io(clim_type,src_field, hdata, nt,i, Time)
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
type(fieldtype)          , intent(in)  :: src_field
integer                  , intent(in)  :: nt
real                     , intent(out) :: hdata(:,:,:)
integer        , optional, intent(in)  :: i
type(time_type), optional, intent(in)  :: Time

integer   :: k, km
! sjs
real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%lon(:)),size(clim_type%lat(:)), &
                        size(clim_type%levs(:))))

      call mpp_read(clim_type%unit,src_field, climdata,nt)

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%lon(:)),   &
                           size(clim_type%lat(:)), &
                           size(clim_type%levs(:))))
        km = size(clim_type%levs(:))
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


end subroutine interp_read_data_io

!> \brief interp_read_data_io_no_time_axis receives various climate data as inputs and
!!        returns a horizontally interpolated climatology field without the
!!            time axis.
!!
!! \param [in] <clim_type> The interpolate type which contains the data
!! \param [in] <src_field> The field type
!! \param [in] <i> OPTIONAL: The index of the field name that you are trying to read
!! \param [out] <hdata> The horizontally interpolated climatology field. This
!                       field will still be on the climatology vertical grid
subroutine interp_read_data_io_no_time_axis(clim_type,src_field, hdata, i)

!  INTENT IN
!    clim_type : The interpolate type which contains the data
!    src_field : The field type
!    i         : The index of the field name that you are trying to read. (optional)

!  INTENT OUT

!    hdata     : The horizontally interpolated climatology field. This
!                field will still be on the climatology vertical grid.

type(interpolate_type)   , intent(in)  :: clim_type
type(fieldtype)          , intent(in)  :: src_field
real                     , intent(out) :: hdata(:,:,:)
integer        , optional, intent(in)  :: i

integer   :: k, km
! sjs
real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

      allocate(climdata(size(clim_type%lon(:)),size(clim_type%lat(:)), size(clim_type%levs(:))))

      call mpp_read(clim_type%unit,src_field, climdata)

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
      if (clim_type%vertical_indices == INCREASING_UPWARD) then
        allocate(climdata2(size(clim_type%lon(:)),   &
                           size(clim_type%lat(:)), &
                           size(clim_type%levs(:))))
        km = size(clim_type%levs(:))
        do k=1, km
          climdata2(:,:,k) = climdata(:,:,km+1-k)
        end do
        climdata = climdata2
        deallocate (climdata2)
      endif

      call horiz_interp(clim_type%interph, climdata, hdata)
      deallocate(climdata)

end subroutine interp_read_data_io_no_time_axis
