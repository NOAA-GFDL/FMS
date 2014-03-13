module oda_core_mod
!
!<CONTACT EMAIL="matthew.harrison@noaa.gov"> Matthew Harrison
!</CONTACT>
!
!<OVERVIEW>
! core ocean data assimilation subroutines for ocean models. This module
! includes interfaces to :
!  
! (i)   initialize the ODA core with observations and model
! (ii)  request ocean state increments using observed or idealized data 
! (iii) calculate model guess differences and output to file
! (iv)  terminate the DA core.  
!  
! NOTE: Profile files conform to v0.1a profile metadata standard.
! This metadata standard will evolve in time and will maintain
! backward compability.  
!  
! NOTE: This code is still under development. ODA data structures should be 
! opaque in future releases.
!  
!</OVERVIEW>
!
!<NAMELIST NAME="oda_core_nml">
! <DATA NAME="max_misfit" TYPE="real">
!  flag measurement if abs(omf) exceeds max_misfit
!</DATA>  
! <DATA NAME="min_obs_err_t" TYPE="real">
! minumum rms temperature observation error (degC) for
! variable obs error, else nominal temp error  
!</DATA>    
! <DATA NAME="min_obs_err_s" TYPE="real">
! minimum rms salinity observation error (g/kg) for
! variable obs error,else nominal salt error  
!</DATA>    
! <DATA NAME="eta_tide_const" TYPE="real">
! Tidal internal displacement amplitude (meters) for observational error estimation.  
!</DATA>    
! <DATA NAME="min_prof_depth" TYPE="real">
! Minimum profile depth (meters)  
!</DATA>    
! <DATA NAME="max_prof_spacing" TYPE="real">
! Data must be contiguous to this resolution (meters). Otherwise flag is activated.
! </DATA>
! <DATA NAME="data_window" TYPE="integer">
! Half-width of profile time window (days)
! </DATA>
! <DATA NAME="add_tidal_aliasing" TYPE="logical">
! Add tidal aliasing to observation error. Use eta_tide_const * dT/dz to estimate
! observation.  
! </DATA>
! <DATA NAME="max_profiles" TYPE="integer">
! Allocation size of profile array
! </DATA>
! <DATA NAME="max_sfc_obs" TYPE="integer">
! Allocation size of surface data array
! </DATA>
! <DATA NAME="temp_obs_rmse" TYPE="real">
! nominal temperature error rms error
! </DATA>
! <DATA NAME="salt_obs_rmse" TYPE="real">
! nominal salinity error rms error
! </DATA>      
!</NAMELIST>  
  use fms_mod, only : file_exist,read_data
  use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_sum, stdout,&
                      mpp_sync_self, mpp_pe,mpp_npes,mpp_root_pe,&
                      mpp_broadcast, input_nml_file
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
       domain2d, mpp_get_global_domain, mpp_update_domains
  use time_manager_mod, only : time_type, operator( <= ), operator( - ), &
       operator( > ), operator ( < ),  set_time, set_date, &
       get_date, get_time
  use get_cal_time_mod, only : get_cal_time
  use axis_utils_mod, only : frac_index  
  use constants_mod, only : radian, pi
  use oda_types_mod
  use write_ocean_data_mod, only : write_ocean_data_init
  
  implicit none

  private
  
  real, private :: max_misfit = 5.0 ! reject fg diff gt max_misfit


  real, parameter, private :: depth_min=0.0, depth_max=10000.
  real, parameter, private :: temp_min=-3.0, temp_max=40.
  real, parameter, private :: salt_min=0.0, salt_max=45.
  
  integer :: max_profiles = 250000, max_sfc_obs = 1 
  integer, parameter, private :: max_files=100 
  integer, parameter, private :: PROFILE_FILE = 1,SFC_FILE= 2,&
                                 IDEALIZED_PROFILES=3
  
  
! parameters for obs error inflation due to internal tidal displacements
  
  real, private :: min_obs_err_t = 0.5, min_obs_err_s=0.1, eta_tide_const = 7.0
  
  type(ocean_profile_type), target, save, private, allocatable  :: profiles(:)
  type(ocean_surface_type), target, save, private, allocatable  :: sfc_obs(:)  

  integer, private, save :: num_profiles, num_sfc_obs ! total number of observations 

  integer, private, save :: isc, iec, jsc, jec, isd, ied, jsd, jed  ! indices for local and global domain
  integer, private, save :: isg, ieg, jsg, jeg
  integer, private, save :: nk

! parameters used to restrict profiles
  real, private, save :: min_prof_depth = 200.0 ! profile data ending
                                                   ! above this level are
                                                   ! discarded.  
  real, private, save :: max_prof_spacing = 1.e5 ! reject profile data
                                                    ! if not contiguous
                                                    ! at this resolution

! internal arrays for forward and backward observations
  real, dimension(:,:,:), allocatable, private, save :: sum_wgt, nobs

  type(time_type) , dimension(0:100), public :: time_window

! the DA module maintains a unique grid,domain association
  
  type(grid_type), pointer :: Grd

! DA grid (1-d) coordinates.  Generalized horizontal coordinates are not supporte.
  
  real, allocatable, dimension(:) :: x_grid, y_grid

  real :: temp_obs_rmse = 0.7071
  real :: salt_obs_rmse = 0.1
  
  logical :: add_tidal_aliasing=.false.

  logical :: localize_data = .true.

  logical :: debug=.false.

  integer :: ndebug=10

  integer, allocatable, dimension(:) :: nprof_in_comp_domain
  
  namelist /oda_core_nml/ max_misfit,add_tidal_aliasing,min_obs_err_t,&
                          min_obs_err_s, eta_tide_const, debug, max_profiles,&
                          max_sfc_obs, temp_obs_rmse, salt_obs_rmse, ndebug
                          

! NOTE: Surface observation type is not yet implemented.

! transform from Grd => Obs
  
  interface forward_obs
     module procedure forward_obs_profile
     module procedure forward_obs_sfc
  end interface

! transform from Obs => Grd
  
  interface backward_obs
     module procedure backward_obs_profile
     module procedure backward_obs_sfc
  end interface


! one-time association between profile and grid, used for storing
! interpolation weights
  
  interface assign_forward_model
     module procedure assign_forward_model_profile
     module procedure assign_forward_model_sfc
  end interface

! duplicate observation array
  
  interface copy_obs
     module procedure copy_obs_prof
     module procedure copy_obs_sfc
  end interface

! multiply observation data by inverse error variance
 
  interface mult_obs_I_mse
     module procedure mult_obs_I_mse_profile
     module procedure mult_obs_I_mse_sfc
  end interface

! difference between two observations (e.g. model first guess and observations)
  
  interface diff_obs
     module procedure diff_obs_profile
     module procedure diff_obs_sfc
  end interface

! inflate observational error based on first guess misfit 
! and time window.
  
  interface adjust_obs_error
     module procedure adjust_obs_error_profile
     module procedure adjust_obs_error_sfc
  end interface


  interface nullify_obs
     module procedure nullify_obs_prof
  end interface
  
  public :: forward_obs, backward_obs, &
            copy_obs, adjust_obs_error, &
            oda_core_init, open_profile_dataset, get_obs, &
            assign_forward_model, diff_obs, mult_obs_I_mse, &
            purge_obs, nullify_obs
  
  contains


  subroutine open_profile_dataset(filename, localize)  
!    <DESCRIPTION>
!    open dataset containing profile information in fms station format.
!    store internally.
!    </DESCRIPTION>    
    
    use mpp_io_mod, only : mpp_open, mpp_get_atts, mpp_get_info, &
         mpp_get_fields, mpp_read, MPP_SINGLE, MPP_MULTI, MPP_NETCDF,&
         axistype, atttype, fieldtype, MPP_RDONLY, mpp_get_axes, mpp_close,&
         mpp_get_att_char
    
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: localize
    logical :: found_neighbor,continue_looking
    integer, parameter :: max_levels=1000
    integer :: unit, ndim, nvar, natt, nstation
    character(len=32) :: fldname, axisname
    type(atttype), allocatable, dimension(:), target :: global_atts
    type(atttype), pointer :: version => NULL()
    type(axistype), pointer :: depth_axis => NULL(), station_axis => NULL()
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    
    type(fieldtype), pointer :: field_lon => NULL(), field_lat => NULL(), field_probe => NULL(),&
         field_time => NULL(), field_depth => NULL()
    type(fieldtype), pointer :: field_temp => NULL(), field_salt => NULL(), &
                                field_error => NULL(), field_link => NULL()
    ! NOTE: fields are restricted to be in separate files
    real :: lon, lat, time, depth(max_levels), temp(max_levels), salt(max_levels),&
         error(max_levels), rprobe, profile_error, rlink
    integer :: nlev, probe, yr, mon, day, hr, minu, sec, kl, outunit
    integer :: num_levs, num_levs_temp, num_levs_salt,&
               k, kk, ll, i, i0, j0, k0, nlevs, a, nn, ii, nlinks
    real :: ri0, rj0, rk0, dx1, dx2, dy1, dy2
    character(len=128) :: time_units, attname, catt
    type(time_type) :: profile_time
    integer :: flag_t(max_levels), flag_s(max_levels), cpe
    logical :: data_is_local, &
               continue
    logical :: found_temp=.false., found_salt=.false.
    real, dimension(max_links,max_levels) :: temp_bfr, salt_bfr, depth_bfr
    integer, dimension(max_links,max_levels) :: flag_t_bfr, flag_s_bfr
    real :: temp_missing=missing_value,salt_missing=missing_value,&
         depth_missing=missing_value
    real :: max_prof_depth, zdist
    
    
    ! read timeseries of local observational data from NetCDF files
    ! and allocate ocean_obs arrays.  

    ! File structure:
    !    dimensions:
    !       depth_index;station=UNLIMITED;
    !    variables:
    !       depth_index(depth_index);
    !       station(station);
    !       longitude(station);
    !       latitude(station);
    !       time(station); 
    !       data(station,depth_index);
    !       depth(station,depth_index);
    !       probe(station);
    !       err(station, depth_index);

    cpe = mpp_pe()

    dx1 = (x_grid(isc)-x_grid(isc-1))/2.0
    dx2 = (x_grid(iec+1)-x_grid(iec))/2.0
    dy1 = (y_grid(jsc)-y_grid(jsc-1))/2.0
    dy2 = (y_grid(jec+1)-y_grid(jec))/2.0
    
    localize_data = .true.

    if (PRESENT(localize)) localize_data = localize

    call mpp_open(unit,filename,form=MPP_NETCDF,fileset=MPP_SINGLE,&
         threading=MPP_MULTI,action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, nstation)

    outunit = stdout()
    write(outunit,*) 'Opened profile dataset :',trim(filename)

    ! get version number of profiles

    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)

    do i=1,natt
       catt = mpp_get_att_char(global_atts(i))
       select case (lowercase(trim(catt)))
       case ('version')
          version =>  global_atts(i)
       end select
    end do

    if (.NOT.ASSOCIATED(version)) then
        call mpp_error(NOTE,'no version number available for profile file, assuming v0.1a ')
    else
        write(outunit,*) 'Reading profile dataset version = ',trim(catt)
    endif
    
    
    ! get axis information

    allocate (axes(ndim))
    call mpp_get_axes(unit,axes)
    do i=1,ndim
       call mpp_get_atts(axes(i),name=axisname)
       select case (lowercase(trim(axisname)))
       case ('depth_index')
          depth_axis => axes(i)
       case ('station_index')
          station_axis => axes(i)
       end select
    end do

    if (.NOT.ASSOCIATED(depth_axis) .or. .NOT.ASSOCIATED(station_axis)) then
        call mpp_error(FATAL,'depth and/or station axes do not exist in input file')
    endif

    
! get selected field information.
! NOTE: not checking for all variables here.    

    allocate(fields(nvar))
    call mpp_get_fields(unit,fields)
    do i=1,nvar
       call mpp_get_atts(fields(i),name=fldname)
       select case (lowercase(trim(fldname)))
       case ('longitude')
           field_lon => fields(i)
       case ('latitude')
           field_lat => fields(i)
       case ('probe') 
           field_probe => fields(i)
       case ('time')
           field_time => fields(i)
       case ('temp')
           field_temp => fields(i)
       case ('salt')
           field_salt => fields(i)           
       case ('depth')
           field_depth => fields(i)
       case ('link')
           field_link => fields(i)
       case ('error')
           field_error => fields(i)
       end select
   enddo

    call mpp_get_atts(depth_axis,len=nlevs)

    if (nlevs > max_levels) call mpp_error(FATAL,'increase parameter max_levels ')

    if (nlevs < 1) call mpp_error(FATAL)

    outunit = stdout()
    write(outunit,*) 'There are ', nstation, ' records in this dataset'
    write(outunit,*) 'Searching for profiles matching selection criteria ...'

    if (ASSOCIATED(field_temp)) found_temp=.true.
    if (ASSOCIATED(field_salt)) found_salt=.true.

    if (.not. found_temp .and. .not. found_salt) then
        write(outunit,*) 'temp or salt not found in profile file'
        call mpp_error(FATAL)
    endif
    
    call mpp_get_atts(field_time,units=time_units)
    if (found_temp) call mpp_get_atts(field_temp,missing=temp_missing)
    if (found_salt) call mpp_get_atts(field_salt,missing=salt_missing)
    
    call mpp_get_atts(field_depth,missing=depth_missing)        

    if (found_salt) then
        write(outunit,*) 'temperature and salinity where available'
    else
        write(outunit,*) 'temperature only records'
    endif

    
    i=1
    continue=.true.
    
    do while (continue)

       depth(:) = missing_value
       temp(:)  = missing_value
       salt(:)  = missing_value       
! required fields       
       call mpp_read(unit,field_lon,lon,tindex=i)
       call mpp_read(unit,field_lat,lat, tindex=i)
       call mpp_read(unit,field_time,time,tindex=i)
       call mpp_read(unit,field_depth,depth(1:nlevs),tindex=i)
       if (found_temp) call mpp_read(unit,field_temp,temp(1:nlevs),tindex=i)
       if (found_salt) call mpp_read(unit,field_salt,salt(1:nlevs),tindex=i)       
! not required fields
       if (ASSOCIATED(field_error)) then
           call mpp_read(unit,field_error,profile_error,tindex=i)
       endif
       if (ASSOCIATED(field_probe)) then
           call mpp_read(unit, field_probe, rprobe,tindex=i)
       endif
       if (ASSOCIATED(field_link)) then
           call mpp_read(unit,field_link,rlink,tindex=i)
       else
           rlink = 0.0
       endif
       probe=rprobe
       data_is_local = .false.
! NOTE: assuming grid is modulo 360 here. This needs to be generalized.
       
       if (lon .lt. x_grid(isg-1) ) lon = lon + 360.0
       if (lon .gt. x_grid(ieg+1) ) lon = lon - 360.0

! localized data is within region bounded by halo points
! (halo size = 1) adjacent to boundary points of computational domain
       
       if (lon >= x_grid(isc-1) .and. &
           lon <  x_grid(iec+1) .and. &
           lat >= y_grid(jsc-1) .and. &
           lat <  y_grid(jec+1)) data_is_local = .true.

       
       profile_time = get_cal_time(time,time_units,'julian')

       
       if ( data_is_local .OR. .NOT.localize_data) then
           
           num_profiles=num_profiles+1
           if (num_profiles > max_profiles) then
               call mpp_error(FATAL,'maximum number of profiles exceeded.&
                    &Resize parameter max_profiles in ocean_obs_mod')
               
           endif

           call nullify_obs(Profiles(num_profiles))
           
           num_levs_temp = 0
           num_levs_salt = 0           
           do k = 1, nlevs
              
! flag=0 denotes a valid profile level, anything else
! is invalid. See NODC codes.
!================================================================
!0 -     accepted station
!1 -     failed annual standard deviation check
!2 -     two or more density inversions (Levitus, 1982 criteria)
!3 -     flagged cruise
!4 -     failed seasonal standard deviation check
!5 -     failed monthly standard deviation check
!6 -     flag 1 and flag 4
!7 -     bullseye from standard level data or failed annual and monthly
!        standard deviation check
!8 -     failed seasonal and monthly standard deviation check
!9 -     failed annual, seasonal, and monthly standard deviation check
!================================================================
              
              flag_t(k) = 0;flag_s(k) = 0
              
              if (.not.found_salt) then
                  flag_s(k) = 1
              endif
              
              if (depth(k) .eq. depth_missing .or. depth(k) .lt.depth_min&
                   .or. depth(k) .gt. depth_max) then
                  depth(k) = missing_value
                  flag_t(k)=1
                  flag_s(k)=1
              endif
              
              if (found_temp .and. flag_t(k) .eq. 0) then
                  if (temp(k) .eq. temp_missing .or. temp(k) .lt. temp_min&
                       .or. temp(k) .gt. temp_max) then
                      temp(k) = missing_value
                      flag_t(k) = 1
                      flag_s(k) = 1 ! flag salt if no temperature data
                  else 
                      num_levs_temp = num_levs_temp+1
                  endif
              endif
              if (found_salt .and. flag_s(k) .eq. 0) then
                  if (salt(k) .eq. salt_missing .or. salt(k) .lt. salt_min&
                       .or. salt(k) .gt. salt_max) then
                      salt(k) = missing_value
                      flag_s(k) = 1
                  else 
                      num_levs_salt = num_levs_salt+1
                  endif
              endif
              
           enddo

! large profile are stored externally in separate records
! follow the links to get complete profile
           
           ii=i+1
           nlinks = 0
           do while (rlink > 0.0 .and. nlinks .lt. max_links)
              nlinks=nlinks+1
              depth_bfr(nlinks,:) = missing_value
              temp_bfr(nlinks,:) = missing_value
              salt_bfr(nlinks,:) = missing_value              
              call mpp_read(unit,field_depth,depth_bfr(nlinks,1:nlevs),tindex=ii)
              if (found_temp) call mpp_read(unit,field_temp,temp_bfr(nlinks,1:nlevs),tindex=ii)
              if (found_salt) call mpp_read(unit,field_salt,salt_bfr(nlinks,1:nlevs),tindex=ii)
              call mpp_read(unit,field_link,rlink,tindex=ii)
              ii=ii+1
           enddo
           i=ii ! set record counter to start of next profile

           if (nlinks > 0) then
               do nn = 1, nlinks
                  do k=1, nlevs
                     flag_t_bfr(nn,k) = 0
                     flag_s_bfr(nn,k) = 0
                     if (depth_bfr(nn,k) .eq. depth_missing .or.&
                          depth_bfr(nn,k) .lt. depth_min .or. &
                          depth_bfr(nn,k) .gt. depth_max) then
                         depth_bfr(nn,k) = missing_value
                         flag_t_bfr(nn,k)  = 1
                         flag_s_bfr(nn,k)  = 1
                     endif
                     if (found_temp .and. flag_t_bfr(nn,k) .eq. 0) then
                         if (temp_bfr(nn,k) .eq. temp_missing .or.&
                              temp_bfr(nn,k) .lt. temp_min .or.&
                              temp_bfr(nn,k) .gt. temp_max) then
                             temp_bfr(nn,k) = missing_value
                             flag_t_bfr(nn,k) = 1
                             flag_s_bfr(nn,k) = 1                         
                         else 
                             num_levs_temp = num_levs_temp+1
                         endif
                     endif
                     if (found_salt .and. flag_s_bfr(nn,k) .eq. 0) then
                         if (salt_bfr(nn,k) .eq. salt_missing  .or.&
                              salt_bfr(nn,k) .lt. salt_min .or.&
                              salt_bfr(nn,k) .gt. salt_max) then
                             salt_bfr(nn,k) = missing_value
                             flag_t_bfr(nn,k) = 0                         
                             flag_s_bfr(nn,k) = 1
                         else
                             num_levs_salt = num_levs_salt+1
                         endif
                     endif
                  enddo
               enddo
           endif

           num_levs = max(num_levs_temp,num_levs_salt)
           
           if (num_levs == 0) then
               if (i .gt. nstation) continue = .false.
               cycle
           endif

           allocate(profiles(num_profiles)%depth(num_levs))
           profiles(num_profiles)%depth=missing_value
           if (num_levs_temp .gt. 0) then
               allocate(profiles(num_profiles)%data_t(num_levs))
               profiles(num_profiles)%data_t=missing_value
               allocate(profiles(num_profiles)%flag_t(num_levs))
               profiles(num_profiles)%flag_t= 1
               profiles(num_profiles)%nvar=1
           endif

           if (num_levs_salt .gt. 0) then
               allocate(profiles(num_profiles)%data_s(num_levs))
               profiles(num_profiles)%data_s=missing_value
               allocate(profiles(num_profiles)%flag_s(num_levs))
               profiles(num_profiles)%flag_s= 1
               profiles(num_profiles)%nvar=profiles(num_profiles)%nvar + 1
           endif
           

           if (probe < 1 )   probe = 0
           profiles(num_profiles)%probe = probe
           profiles(num_profiles)%levels = num_levs
           profiles(num_profiles)%lat = lat
           profiles(num_profiles)%lon = lon
           allocate(profiles(num_profiles)%ms_t(num_levs))
           profiles(num_profiles)%ms_t(:) = temp_obs_rmse**2.0 ! default error variance for temperature

           if(num_levs_salt .gt. 0) then
               allocate(profiles(num_profiles)%ms_s(num_levs))
               profiles(num_profiles)%ms_s(:) = salt_obs_rmse**2.0  ! default error variance for salinity
           endif
           
           kk= 1
           do k = 1, nlevs
              if (flag_t(k) .eq. 0) then
                  if (kk > profiles(num_profiles)%levels) then
                      call mpp_error(FATAL)
                  endif
                  profiles(num_profiles)%depth(kk) = depth(k)
                  profiles(num_profiles)%data_t(kk) = temp(k)
                  profiles(num_profiles)%flag_t(kk) = 0                  
                  if (found_salt .and. flag_s(k) .eq. 0) then
                      profiles(num_profiles)%data_s(kk) = salt(k)
                      profiles(num_profiles)%flag_s(kk) = 0
                  endif
                  kk=kk+1
              endif
           enddo

           do nn = 1, nlinks
              do k = 1, nlevs
                 if (flag_t_bfr(nn,k) .eq. 0) then
                     if (kk > profiles(num_profiles)%levels) then
                         call mpp_error(FATAL)
                     endif
                     profiles(num_profiles)%depth(kk) = depth_bfr(nn,k)
                     profiles(num_profiles)%data_t(kk) = temp_bfr(nn,k)
                     profiles(num_profiles)%flag_t(kk) = 0
                     if (found_salt .and. flag_s_bfr(nn,k) .eq. 0) then
                         profiles(num_profiles)%data_s(kk) = salt_bfr(nn,k)
                         profiles(num_profiles)%flag_s(kk) = 0
                     endif
                     kk=kk+1
                 endif
              enddo
           enddo
           
           profiles(num_profiles)%time = profile_time
           
! calculate interpolation coefficients (make sure to account for grid offsets here!)
! NOTE: this only works for lat/lon grids. Lower left indices.           
!       
           
           ri0 = frac_index(lon, x_grid(isg-1:ieg+1)) - 1.
           rj0 = frac_index(lat, y_grid(jsg-1:jeg+1)) - 1.
           i0 = floor(ri0)
           j0 = floor(rj0)
           Profiles(num_profiles)%i_index = ri0
           Profiles(num_profiles)%j_index = rj0           
           Profiles(num_profiles)%accepted = .true.
           if (i0 < 0 .or. j0 < 0) then
               Profiles(num_profiles)%accepted = .false.
           endif
           if (i0 > ieg .or. j0 > jeg) then  
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
           if ((i0 < isc-1 .or. i0 > iec) .and. localize_data) then
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
           if ((j0 < jsc-1 .or. j0 > jec) .and. localize_data) then
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
!
! flag the profile if it sits on a model land point
!           
           if (Profiles(num_profiles)%accepted ) then
               if (Grd%mask(i0,j0,1) == 0.0 .or. &
                    Grd%mask(i0+1,j0,1) == 0.0 .or. &
                    Grd%mask(i0,j0+1,1) == 0.0 .or. &
                    Grd%mask(i0+1,j0+1,1) == 0.0) then
                   Profiles(num_profiles)%accepted = .false.
               endif
           endif

           if (Profiles(num_profiles)%accepted) then
               allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
               max_prof_depth=0.0
               do k=1, Profiles(num_profiles)%levels
                  k0=0
                  if (Profiles(num_profiles)%flag_t(k).eq.0) then
                      rk0 = frac_index(Profiles(num_profiles)%depth(k), Grd%z(:))
                      k0 = floor(rk0)
                      if ( k0 == -1) then
                          if (Profiles(num_profiles)%depth(k) .lt. Grd%z(1)) then
                              k0 = 1
                              rk0 = 1.0
                          else if (Profiles(num_profiles)%depth(k) .gt. Grd%z(Grd%nk)) then
                              Profiles(num_profiles)%flag_t(k) = 1
                          endif
                      endif
                  else
                      cycle
                  endif

                  if (k0 .gt. size(Grd%z)-1 ) then
                      write(*,*) 'k0 out of bounds, rk0,k0= ',rk0,k0
                     write(*,*) 'Z_bound= ',Grd%z_bound
                     write(*,*) 'Profile%depth= ',Profiles(num_profiles)%depth
                     
                      call mpp_error(FATAL)
                  endif
                  
                  Profiles(num_profiles)%k_index(k) = rk0

! flag depth level if adjacent model grid points are land

                  if (Profiles(num_profiles)%flag_t(k) .eq. 0) then
                      if (i0 .lt. 0 .or. j0 .lt. 0 .or. k0 .lt. 0) then
                          write(*,*) 'profile index out of bounds'
                          write(*,*) 'i0,j0,k0=',i0,j0,k0
                          write(*,*) 'lon,lat,depth=',Profiles(num_profiles)%lon,&
                               Profiles(num_profiles)%lat,Profiles(num_profiles)%depth(k)
                          call mpp_error(FATAL)
                      endif
                      
                      if (Grd%mask(i0,j0,k0) == 0.0 .or. &
                          Grd%mask(i0+1,j0,k0) == 0.0 .or. &
                          Grd%mask(i0,j0+1,k0) == 0.0 .or. &
                          Grd%mask(i0+1,j0+1,k0) == 0.0) then
                          Profiles(num_profiles)%flag_t(k) = 1
                      endif
                      if (Grd%mask(i0,j0,k0+1) == 0.0 .or. &
                          Grd%mask(i0+1,j0,k0+1) == 0.0 .or. &
                          Grd%mask(i0,j0+1,k0+1) == 0.0 .or. &
                          Grd%mask(i0+1,j0+1,k0+1) == 0.0) then
                          Profiles(num_profiles)%flag_t(k) = 1
                      endif
                      if (Profiles(num_profiles)%flag_t(k) .eq. 0) then
                          max_prof_depth = Profiles(num_profiles)%depth(k)
                      endif
                  endif

               enddo ! Prof%levels loop

! Flag profile if it is too shallow.
               
               if (max_prof_depth .lt. min_prof_depth) then
                   Profiles(num_profiles)%accepted = .false.
               endif

               found_neighbor=.false.
               
               do k=2,Profiles(num_profiles)%levels - 1
                  if (Profiles(num_profiles)%flag_t(k) .eq. 0) then
                      kk = k-1
                      found_neighbor = .false.
                      continue_looking = .true.
                      do while (continue_looking .and. kk .ge. 1)
                         if (Profiles(num_profiles)%flag_t(kk) .eq. 0) then
                             zdist = Profiles(num_profiles)%depth(k) - Profiles(num_profiles)%depth(kk)
                             if (zdist .gt. max_prof_spacing) then
                                 Profiles(num_profiles)%accepted = .false.
                                 goto 199
                             else 
                                 continue_looking = .false.
                                 found_neighbor = .true.
                             endif
                         else
                             kk = kk - 1
                         endif
                      end do
                      kk = k+1
                      continue_looking = .true.
                      do while (continue_looking .and. kk .le. Profiles(num_profiles)%levels)
                         if (Profiles(num_profiles)%flag_t(kk).eq. 0) then
                             zdist = Profiles(num_profiles)%depth(kk) - Profiles(num_profiles)%depth(k)
                             if (zdist .gt. max_prof_spacing) then
                                 Profiles(num_profiles)%accepted = .false.
                                 goto 199
                             else
                                 continue_looking = .false.
                                 found_neighbor = .true.
                             endif
                         else
                             kk = kk+1
                         endif
                      enddo
                  endif
               enddo

               if (.not. found_neighbor) Profiles(num_profiles)%accepted = .false.

199            continue
               
           endif ! if Prof%accept
       else ! data is not local
           i = i+1
       endif ! if data_is_local

       
       if (i .gt. nstation) continue = .false.

    enddo

!    a = nprof_in_comp_domain(cpe)


    
!    call mpp_broadcast(nprof_in_comp_domain(cpe),cpe)
    
!    call mpp_sum(a)

!    write(stdout(),*) 'A grand total of ',int(a),' profiles satisify acceptance criteria'

!    do i=0,mpp_npes()
!       write(stdout(),*) 'pe=',i,'profile count=',nprof_in_comp_domain(i)
!    enddo
    
    call mpp_sync_self()
    call mpp_close(unit)

  end subroutine open_profile_dataset

  subroutine get_obs(model_time, Prof, Sfc, nprof, nsfc)


    ! get profiles and sfc
    ! obs relevant to current analysis interval

    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:) :: Prof
    type(ocean_surface_type), dimension(:) :: Sfc
    integer, intent(inout) :: nprof, nsfc

    integer :: i,k,yr,mon,day,hr,minu,sec,a,mon_obs,yr_obs, outunit
    type(time_type) :: tdiff
    character(len=1) :: cchar

    nprof=0
    nsfc=0

    outunit = stdout()
    write(outunit,*) 'Gathering profiles for current analysis time'
    call get_date(model_time,yr,mon,day,hr,minu,sec)
    write(outunit,'(a,i4,a,i2,a,i2)') 'Current yyyy/mm/dd= ',yr,'/',mon,'/',day
    write(outunit,*) 'num_profiles=',num_profiles

    
    do i=1,num_profiles

       if (debug .and. i.le.ndebug) then
           call get_date(Profiles(i)%time,yr,mon,day,hr,minu,sec)
           write(*,*) 'in get_obs prof time: yy/mm/dd= ',yr,mon,day
       endif
       
       if (Profiles(i)%time <= model_time) then
           tdiff = model_time - Profiles(i)%time
       else
           tdiff = Profiles(i)%time - model_time
       endif

       if (debug .and. i .le. ndebug) then
           write(*,*) 'Prof%accepted=',Profiles(i)%accepted
       endif
       
       if (tdiff <= time_window(0) .and. &
            Profiles(i)%accepted) then
           nprof=nprof+1
           if (nprof > size(Prof,1)) &
                call mpp_error(FATAL,'increase size of Prof array before call to get_obs')
           call copy_obs(Profiles(i:i),Prof(nprof:nprof))
           Prof(nprof)%tdiff = tdiff
           if (debug .and. nprof .le. ndebug) then
              call get_time(tdiff,sec,day)
              write(*,'(a,i3,a,2f10.5)') 'Accepted profile #',i,' : lon,lat= ',Prof(nprof)%lon,Prof(nprof)%lat
              do k=1,Prof(nprof)%levels
                 if (Prof(nprof)%nvar .eq. 2) then
                     write(*,'(a,i3,a,2f10.5,2i2,2f8.5)') 'Profile #',i,' : temp,salt,flag_t,flag_s,ms_t,ms_s= ',&
                          Prof(nprof)%data_t(k),Prof(nprof)%data_s(k),Prof(nprof)%flag_t(k),Prof(nprof)%flag_s(k),&
                          Prof(nprof)%ms_t(k),Prof(nprof)%ms_s(k)
                 else
                     write(*,'(a,i3,a,2f10.5)') 'Profile #',i,' : temp,flag_t= ',Prof(nprof)%data_t(k),Prof(nprof)%flag_t(k)
                 endif
              enddo
          endif
          
      else
          if (debug .and. i .le. ndebug) then
              call get_time(tdiff,sec,day)
              write(*,'(a,i3,a,2f10.5)') 'Rejected profile #',i,' : lon,lat= ',Prof(i)%lon,Prof(i)%lat
              do k=1,Prof(i)%levels
                 if (Prof(i)%nvar .eq. 2) then
                     write(*,'(a,i3,a,2f10.5,2i2)') 'Profile #',i,' : temp,salt,flag_t,flag_s= ',Prof(i)%data_t(k),Prof(i)%data_s(k),Prof(i)%flag_t(k),Prof(i)%flag_s(k)
                 else
                     write(*,'(a,i3,a,2f10.5)') 'Profile #',i,' : temp,flag_t= ',Prof(i)%data_t(k),Prof(i)%flag_t(k)
                 endif
              enddo
          endif
      endif

    enddo

    a=nprof
    call mpp_sum(a)
    write(outunit,*) 'A total of ',a,'  profiles are being used for the current analysis step'

    return

  end subroutine get_obs

  subroutine oda_core_init(Domain, Grid, localize)

    use fms_mod, only : open_namelist_file, check_nml_error, close_file
    
    type(domain2d), intent(in) :: Domain
    type(grid_type), target, intent(in) :: Grid
    logical, intent(in), optional :: localize

      
    integer :: ioun, ierr, io_status
    
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, oda_core_nml, iostat=io_status)
      ierr = check_nml_error(io_status,'oda_core_nml')
#else
    ioun = open_namelist_file()
    read(ioun,nml=oda_core_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'oda_core_nml')
    call close_file(ioun)      
#endif

!    allocate(nprof_in_comp_domain(0:mpp_npes()-1))
!    nprof_in_comp_domain = 0
    
    Grd => Grid
    
    call mpp_get_compute_domain(Domain,isc,iec,jsc,jec)
    call mpp_get_data_domain(Domain,isd,ied,jsd,jed)
    call mpp_get_global_domain(Domain,isg,ieg,jsg,jeg)
    nk = size(Grid%z)

    allocate(sum_wgt(isd:ied,jsd:jed,nk))
    allocate(nobs(isd:ied,jsd:jed,nk))

    if (PRESENT(localize)) localize_data = localize
    
    call init_observations(localize_data)

    call write_ocean_data_init()
    
  end subroutine oda_core_init

  subroutine purge_obs()

    num_profiles=0
    num_sfc_obs=0

  end subroutine purge_obs

  subroutine forward_obs_profile(Model_obs, fg_t, fg_s)

! map first guess to observation locations
! note that forward operator only becomes associated after
! this call

    type(ocean_profile_type), dimension(:), intent(inout) ::Model_obs    
    real, dimension(isd:ied,jsd:jed,nk) , intent(in) :: fg_t ! first guess for temperature
    real, dimension(isd:ied,jsd:jed,nk) , intent(in), optional :: fg_s ! first guess for salinity

    integer :: n, i0, j0, k, num_prof, k0
    real :: a,b,c

    character(len=128) :: mesg
    
    num_prof = size(Model_obs)

    sum_wgt = 0.0

    do n = 1, num_prof
       i0 = floor(Model_obs(n)%i_index)
       j0 = floor(Model_obs(n)%j_index)
       a = Model_obs(n)%i_index - i0 
       b = Model_obs(n)%j_index - j0 

       if (a >= 1.0 .or. a < 0.0 ) call mpp_error(FATAL)
       if (b >= 1.0 .or. b < 0.0 ) call mpp_error(FATAL)


       if (i0 < isc-1 .or. i0 > iec) then
           write(mesg,'(a,i4,2x,i4)') 'out of bounds in forward_obs: i0,j0= ',i0,j0
           call mpp_error(FATAL,trim(mesg))
       endif

       if (j0 < jsc-1 .or. j0 > jec) then
           write(mesg,*) 'out of bounds in forward_obs: i0,j0= ',i0,j0
           call mpp_error(FATAL,trim(mesg))
       endif

       if (ASSOCIATED(Model_obs(n)%data_t) .and. Model_obs(n)%accepted) then

           if (ASSOCIATED(Model_obs(n)%Forward_model%wgt))&
                Model_obs(n)%Forward_model%wgt => NULL()
           allocate(Model_obs(n)%Forward_model%wgt(2,2,2,Model_obs(n)%levels))

           Model_obs(n)%Forward_model%wgt = 0.0
           
           do k = 1, Model_obs(n)%levels
              if (Model_obs(n)%flag_t(k) .eq. 0) then
                  k0 = floor(Model_obs(n)%k_index(k))
                  if (k0 < 1 .or. k0 > Grd%nk-1) then
                      write(mesg,*) 'out of bounds in forward_obs: k0= ',k0
                      call mpp_error(FATAL,trim(mesg))
                  endif          
                  c =	 Model_obs(n)%k_index(k) - k0

                  if (c >= 1.0 .or. c < 0.0 ) call mpp_error(FATAL)

                  Model_obs(n)%Forward_model%wgt(1,1,1,k) = (1.0-a)*(1.0-b)*(1.0-c)
                  Model_obs(n)%Forward_model%wgt(2,1,1,k) = a*(1.0-b)*(1.0-c)
                  Model_obs(n)%Forward_model%wgt(1,2,1,k) = (1.0-a)*b*(1.0-c)
                  Model_obs(n)%Forward_model%wgt(2,2,1,k) = a*b*(1.0-c)
                  Model_obs(n)%Forward_model%wgt(1,1,2,k) = (1.0-a)*(1.0-b)*c
                  Model_obs(n)%Forward_model%wgt(2,1,2,k) = a*(1.0-b)*c
                  Model_obs(n)%Forward_model%wgt(1,2,2,k) = (1.0-a)*b*c
                  Model_obs(n)%Forward_model%wgt(2,2,2,k) = a*b*c

                  sum_wgt(i0,j0,k0) = sum_wgt(i0,j0,k0)+&
                       Model_obs(n)%Forward_model%wgt(1,1,1,k)
                  sum_wgt(i0+1,j0,k0) = sum_wgt(i0+1,j0,k0)+&
                       Model_obs(n)%Forward_model%wgt(2,1,1,k)
                  sum_wgt(i0,j0+1,k0) = sum_wgt(i0,j0+1,k0)+&
                       Model_obs(n)%Forward_model%wgt(1,2,1,k)
                  sum_wgt(i0+1,j0+1,k0) = sum_wgt(i0+1,j0+1,k0)+&
                       Model_obs(n)%Forward_model%wgt(2,2,1,k)
                  sum_wgt(i0,j0,k0+1) = sum_wgt(i0,j0,k0+1)+&
                       Model_obs(n)%Forward_model%wgt(1,1,2,k)
                  sum_wgt(i0+1,j0,k0+1) = sum_wgt(i0+1,j0,k0+1)+&
                       Model_obs(n)%Forward_model%wgt(2,1,2,k)
                  sum_wgt(i0,j0+1,k0+1) = sum_wgt(i0,j0+1,k0+1)+&
                       Model_obs(n)%Forward_model%wgt(1,2,2,k)
                  sum_wgt(i0+1,j0+1,k0+1) = sum_wgt(i0+1,j0+1,k0+1)+&
                       Model_obs(n)%Forward_model%wgt(2,2,2,k)
                  
                      
                  Model_obs(n)%data_t(k) = &
                       fg_t(i0,j0,k0)*Model_obs(n)%Forward_model%wgt(1,1,1,k) &
                     + fg_t(i0+1,j0,k0)*Model_obs(n)%Forward_model%wgt(2,1,1,k) &
                     + fg_t(i0,j0+1,k0)*Model_obs(n)%Forward_model%wgt(1,2,1,k) &
                     + fg_t(i0+1,j0+1,k0)*Model_obs(n)%Forward_model%wgt(2,2,1,k) &
                     + fg_t(i0,j0,k0+1)*Model_obs(n)%Forward_model%wgt(1,1,2,k) &
                     + fg_t(i0+1,j0,k0+1)*Model_obs(n)%Forward_model%wgt(2,1,2,k) &
                     + fg_t(i0,j0+1,k0+1)*Model_obs(n)%Forward_model%wgt(1,2,2,k) &
                     + fg_t(i0+1,j0+1,k0+1)*Model_obs(n)%Forward_model%wgt(2,2,2,k)

                  if (ASSOCIATED(Model_obs(n)%data_s) .and. PRESENT(fg_s)) then
                      Model_obs(n)%data_s(k) = &
                       fg_s(i0,j0,k0)*Model_obs(n)%Forward_model%wgt(1,1,1,k) &
                     + fg_s(i0+1,j0,k0)*Model_obs(n)%Forward_model%wgt(2,1,1,k) &
                     + fg_s(i0,j0+1,k0)*Model_obs(n)%Forward_model%wgt(1,2,1,k) &
                     + fg_s(i0+1,j0+1,k0)*Model_obs(n)%Forward_model%wgt(2,2,1,k) &
                     + fg_s(i0,j0,k0+1)*Model_obs(n)%Forward_model%wgt(1,1,2,k) &
                     + fg_s(i0+1,j0,k0+1)*Model_obs(n)%Forward_model%wgt(2,1,2,k) &
                     + fg_s(i0,j0+1,k0+1)*Model_obs(n)%Forward_model%wgt(1,2,2,k) &
                     + fg_s(i0+1,j0+1,k0+1)*Model_obs(n)%Forward_model%wgt(2,2,2,k) 
                  endif
              else
                  if (ASSOCIATED(Model_obs(n)%data_t)) then
                      Model_obs(n)%data_t(k) = missing_value
                  endif
                  if (ASSOCIATED(Model_obs(n)%data_s)) then
                      Model_obs(n)%data_s(k) = missing_value
                  endif                  
              endif
           enddo
       endif
    enddo

    return

  end subroutine forward_obs_profile





  subroutine forward_obs_sfc(Sfc, Guess, Diff)
    type(ocean_surface_type), intent(in) :: Sfc
    type(ocean_dist_type), intent(in) :: Guess
    type(ocean_surface_type), intent(inout) :: Diff


    return

  end subroutine forward_obs_sfc

  subroutine backward_obs_profile(Obs,model_t,model_s)
!
! map observations back to model locations
!
    type(ocean_profile_type), dimension(:), intent(in) :: Obs
    real, dimension(isd:ied,jsd:jed,nk), intent(inout) :: model_t
    real, dimension(isd:ied,jsd:jed,nk), intent(inout), optional :: model_s    

    real :: tmp
    integer :: i0,j0,k0,k,n

    model_t = 0.0
    if (PRESENT(model_s)) model_s = 0.0
    
    do n=1,size(Obs)
       if (ASSOCIATED(Obs(n)%data_t) .and. Obs(n)%accepted) then
           i0 = floor(Obs(n)%i_index)
           j0 = floor(Obs(n)%j_index) 
! profiles are assumed to lie           
! in domain bounded by first halo points

           if (i0 < isd .or. i0 > ied-1) cycle
           if (j0 < jsd .or. j0 > jed-1) cycle
           
           if (.not.ASSOCIATED(Obs(n)%Forward_model%wgt)) call mpp_error(FATAL,'forward operator not associated with obs')
           
           do k=1, Obs(n)%levels
              if (Obs(n)%flag_t(k) .eq. 0) then
                  k0 = floor(Obs(n)%k_index(k))
                  if (k0 < 1 .or. k0 > Grd%nk-1) call mpp_error(FATAL,'profile k indx out of bnds')
                  
                  tmp = Obs(n)%data_t(k) 
                  model_t(i0,j0,k0)       = model_t(i0,j0,k0) + tmp*Obs(n)%Forward_model%wgt(1,1,1,k)
                  model_t(i0+1,j0,k0)     = model_t(i0+1,j0,k0) + tmp*Obs(n)%Forward_model%wgt(2,1,1,k)
                  model_t(i0,j0+1,k0)     = model_t(i0,j0+1,k0) + tmp*Obs(n)%Forward_model%wgt(1,2,1,k)
                  model_t(i0+1,j0+1,k0)   = model_t(i0+1,j0+1,k0) + tmp*Obs(n)%Forward_model%wgt(2,2,1,k)
                  model_t(i0,j0,k0+1)     = model_t(i0,j0,k0+1) + tmp*Obs(n)%Forward_model%wgt(1,1,2,k)
                  model_t(i0+1,j0,k0+1)   = model_t(i0+1,j0,k0+1) + tmp*Obs(n)%Forward_model%wgt(2,1,2,k)
                  model_t(i0,j0+1,k0+1)   = model_t(i0,j0+1,k0+1) + tmp*Obs(n)%Forward_model%wgt(1,2,2,k)
                  model_t(i0+1,j0+1,k0+1) = model_t(i0+1,j0+1,k0+1) + tmp*Obs(n)%Forward_model%wgt(2,2,2,k)

                  if (PRESENT(model_s)) then

                      tmp = Obs(n)%data_s(k) 
                      model_s(i0,j0,k0)       = model_s(i0,j0,k0) + tmp*Obs(n)%Forward_model%wgt(1,1,1,k)
                      model_s(i0+1,j0,k0)     = model_s(i0+1,j0,k0) + tmp*Obs(n)%Forward_model%wgt(2,1,1,k)
                      model_s(i0,j0+1,k0)     = model_s(i0,j0+1,k0) + tmp*Obs(n)%Forward_model%wgt(1,2,1,k)
                      model_s(i0+1,j0+1,k0)   = model_s(i0+1,j0+1,k0) + tmp*Obs(n)%Forward_model%wgt(2,2,1,k)
                      model_s(i0,j0,k0+1)     = model_s(i0,j0,k0+1) + tmp*Obs(n)%Forward_model%wgt(1,1,2,k)
                      model_s(i0+1,j0,k0+1)   = model_s(i0+1,j0,k0+1) + tmp*Obs(n)%Forward_model%wgt(2,1,2,k)
                      model_s(i0,j0+1,k0+1)   = model_s(i0,j0+1,k0+1) + tmp*Obs(n)%Forward_model%wgt(1,2,2,k)
                      model_s(i0+1,j0+1,k0+1) = model_s(i0+1,j0+1,k0+1) + tmp*Obs(n)%Forward_model%wgt(2,2,2,k)                      
                  endif

              end if

                  
           end do
       end if
    end do


    where(sum_wgt > 0.0)
        model_t = model_t /sum_wgt
    elsewhere
        model_t = 0.0
    end where

    if (PRESENT(model_s)) then
        where(sum_wgt > 0.0)
            model_s = model_s /sum_wgt
        elsewhere
            model_s = 0.0
        end where
    endif
    
  end subroutine backward_obs_profile


  subroutine backward_obs_sfc(Obs,model)

    type(ocean_surface_type), dimension(:), intent(in) :: Obs
    real, dimension(isd:ied,jsd:jed,nk), intent(inout) :: model


  end subroutine backward_obs_sfc
  
  subroutine assign_forward_model_profile(Obs1,Obs2)

    type(ocean_profile_type), dimension(:), target, intent(in) :: Obs1
    type(ocean_profile_type), dimension(:), intent(inout) :: Obs2    


    integer :: n

    if (size(Obs1) .ne. size(Obs2)) call mpp_error(FATAL)

    do n=1,size(Obs1)

       Obs2(n)%Forward_model%wgt => Obs1(n)%Forward_model%wgt

    enddo

  end subroutine assign_forward_model_profile


  subroutine assign_forward_model_sfc(Obs1,Obs2)

    type(ocean_surface_type), target, intent(in) :: Obs1
    type(ocean_surface_type), intent(inout) :: Obs2    
    

    return
  end subroutine assign_forward_model_sfc
  
  subroutine diff_obs_profile(prof1, prof2, Diff)

    type(ocean_profile_type), dimension(:), intent(in) :: prof1
    type(ocean_profile_type), dimension(:), intent(in) :: prof2    
    type(ocean_profile_type), dimension(:), intent(inout) :: Diff

    integer :: n,k

    if (size(prof1) .ne. size(prof2) ) call mpp_error(FATAL)
    
    if (size(prof1) .ne. size(Diff) ) call mpp_error(FATAL)    


    do n=1,size(prof1)
       do k=1,prof1(n)%levels
          if (prof1(n)%flag_t(k) .eq. 0) then
              Diff(n)%data_t(k) = prof2(n)%data_t(k)-prof1(n)%data_t(k)
          else
              Diff(n)%data_t(k) = missing_value
          endif
          if (abs(Diff(n)%data_t(k)) .gt. max_misfit) then
              Diff(n)%flag_t(k) = 1
          endif
         if (ASSOCIATED(prof1(n)%data_s)) then
             if (prof1(n)%flag_s(k) .eq. 0) then
                 Diff(n)%data_s(k) = prof2(n)%data_s(k)-prof1(n)%data_s(k)
             else
                 Diff(n)%data_s(k) = missing_value
             endif
             if (abs(Diff(n)%data_s(k)) .gt. max_misfit) then
                 Diff(n)%flag_s(k) = 1
             endif
         endif
       enddo
    enddo


  end subroutine diff_obs_profile

  subroutine diff_obs_sfc(prof1,prof2,Diff)

    type(ocean_surface_type), dimension(:), intent(in) :: prof1, prof2
    type(ocean_surface_type), dimension(:), intent(inout) :: Diff


  end subroutine diff_obs_sfc
  
  subroutine copy_obs_prof(obs_in, obs_out)

    type(ocean_profile_type), dimension(:), intent(in) :: obs_in
    type(ocean_profile_type), dimension(:), intent(inout) :: obs_out


    integer :: n

    if (size(obs_in) .ne. size(obs_out)) call mpp_error(FATAL,&
         'size mismatch in call to copy_obs_prof')


    do n=1,size(obs_in)
       call nullify_obs(obs_out(n))
       Obs_out(n)%nvar = Obs_in(n)%nvar
       Obs_out(n)%project = Obs_in(n)%project
       Obs_out(n)%probe = Obs_in(n)%probe
       Obs_out(n)%ref_inst = Obs_in(n)%ref_inst
       Obs_out(n)%wod_cast_num = Obs_in(n)%wod_cast_num
       Obs_out(n)%fix_depth = Obs_in(n)%fix_depth
       Obs_out(n)%ocn_vehicle = Obs_in(n)%ocn_vehicle
       Obs_out(n)%database_id = Obs_in(n)%database_id
       Obs_out(n)%levels = Obs_in(n)%levels
       Obs_out(n)%profile_flag = Obs_in(n)%profile_flag
       Obs_out(n)%profile_flag_s = Obs_in(n)%profile_flag_s       
       Obs_out(n)%lon = Obs_in(n)%lon
       Obs_out(n)%lat = Obs_in(n)%lat
       Obs_out(n)%accepted = Obs_in(n)%accepted
       ALLOCATE(Obs_out(n)%depth(Obs_in(n)%levels))
       Obs_out(n)%depth(:) = Obs_in(n)%depth(:)
       ALLOCATE(Obs_out(n)%data_t(Obs_in(n)%levels))
       Obs_out(n)%data_t(:) = Obs_in(n)%data_t(:)
       ALLOCATE(Obs_out(n)%flag_t(Obs_in(n)%levels))
       Obs_out(n)%flag_t(:) = Obs_in(n)%flag_t(:)
       if (ASSOCIATED(Obs_in(n)%data_s)) then
           ALLOCATE(Obs_out(n)%data_s(Obs_in(n)%levels))
           Obs_out(n)%data_s(:) = Obs_in(n)%data_s(:)
           ALLOCATE(Obs_out(n)%flag_s(Obs_in(n)%levels))
           Obs_out(n)%flag_s(:) = Obs_in(n)%flag_s(:)          
       endif
       
       Obs_out(n)%time = Obs_in(n)%time
       Obs_out(n)%yyyy = Obs_in(n)%yyyy
       Obs_out(n)%mmdd = Obs_in(n)%mmdd
       Obs_out(n)%i_index = Obs_in(n)%i_index
       Obs_out(n)%j_index = Obs_in(n)%j_index
       ALLOCATE(Obs_out(n)%k_index(Obs_in(n)%levels))          
       Obs_out(n)%k_index = Obs_in(n)%k_index
       ALLOCATE(Obs_out(n)%ms_t(Obs_in(n)%levels))          
       Obs_out(n)%ms_t = Obs_in(n)%ms_t
       if (ASSOCIATED(Obs_in(n)%ms_s)) then       
           ALLOCATE(Obs_out(n)%ms_s(Obs_in(n)%levels))          
           Obs_out(n)%ms_s = Obs_in(n)%ms_s
       endif
       


       Obs_out(n)%tdiff = Obs_in(n)%tdiff
       Obs_out(n)%nbr_index = Obs_in(n)%nbr_index
       Obs_out(n)%nbr_dist  = Obs_in(n)%nbr_dist
       if (ASSOCIATED(Obs_in(n)%Model_grid)) &
            Obs_out(n)%Model_grid => Obs_in(n)%Model_Grid
    enddo
 
end subroutine copy_obs_prof

  subroutine copy_obs_sfc(Obs_in, Obs_out)
    type(ocean_surface_type), dimension(:), intent(in) :: Obs_in
    type(ocean_surface_type), dimension(:), intent(inout) :: Obs_out


    return

  end subroutine copy_obs_sfc

  subroutine adjust_obs_error_profile(Prof)

    use time_manager_mod, only : get_time
    
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer :: secs, days, n, k, secs_w, days_w, m
    real :: tfac, Ims
    
    do n=1,size(Prof)
       call get_time(Prof(n)%tdiff,secs, days)
       m=Prof(n)%probe
       call get_time(time_window(m),secs_w,days_w)
       tfac = (days + secs/86400.) / days_w
       tfac = 1. - min(1.,tfac)
       if (tfac > 1.0 ) call mpp_error(FATAL)
       if (tfac < 0.0 ) call mpp_error(FATAL)
       do k=1,Prof(n)%levels
           Prof(n)%ms_t(k) = 1.0/ max(1.e-1,tfac) * Prof(n)%ms_t(k)
           if (ASSOCIATED(Prof(n)%data_s)) then
               Prof(n)%ms_s(k) = 1.0/ max(1.e-1,tfac) * Prof(n)%ms_s(k)
          endif
       end do
    end do

  end subroutine adjust_obs_error_profile

  subroutine adjust_obs_error_sfc(Diff)

    type(ocean_surface_type), intent(inout) :: Diff

    return
    
  end subroutine adjust_obs_error_sfc


  subroutine mult_obs_I_mse_profile(Obs)

    type(ocean_profile_type), dimension(:), intent(inout) :: Obs

    integer :: n,k
    real :: Ims
    
    do n=1,size(Obs)
       do k = 1, Obs(n)%levels
          Ims = 1./Obs(n)%ms_t(k)
          if (Obs(n)%flag_t(k) .eq. 0) Obs(n)%data_t(k) = Ims*Obs(n)%data_t(k)
          if (ASSOCIATED(Obs(n)%data_s)) then
              Ims = 1/Obs(n)%ms_s(k)
              if (Obs(n)%flag_s(k) .eq. 0) Obs(n)%data_s(k) = Ims*Obs(n)%data_s(k)
          endif
       end do
    end do
    
  end subroutine mult_obs_I_mse_profile

  subroutine mult_obs_I_mse_sfc(a, Obs)

    real, dimension(:), intent(in) :: a
    type(ocean_surface_type), intent(inout) :: Obs

  end subroutine mult_obs_I_mse_sfc


       

 ! </SUBROUTINE>
! <FUNCTION NAME="lowercase">
!   <DESCRIPTION>
! Turn a string from uppercase to lowercase, do nothing if the
! string is already in lowercase.
!   </DESCRIPTION>
 function lowercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: lowercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    ca = transfer(cs,"x",len(cs)) 
    where (ca >= "A" .and. ca <= "Z") ca = achar(iachar(ca)+co) 
    lowercase = transfer(ca,cs) 
    
  end function lowercase

  subroutine init_observations(localize)  

    use fms_mod, only : open_namelist_file,close_file,check_nml_error
    use mpp_io_mod, only : mpp_open, MPP_ASCII, MPP_RDONLY, MPP_MULTI, MPP_SINGLE
    use mpp_domains_mod, only : mpp_global_field
    
    logical, intent(in), optional :: localize

    integer :: data_window = 15 ! default data half-window is 15 days

    integer :: i,j
    
    
    type obs_entry_type
       character(len=128) :: filename
       character(len=16)  :: file_type
    end type obs_entry_type

    namelist /ocean_obs_nml/ data_window, max_prof_spacing, min_prof_depth
    
    character(len=128) :: input_files(max_files) = ''
    integer :: nfiles, filetype(max_files), ioun, io_status, ierr,&
                unit, nrecs, n
    character(len=256) :: record
    type(obs_entry_type) :: tbl_entry

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, ocean_obs_nml, iostat=io_status)
      ierr = check_nml_error(io_status,'ocean_obs_nml')
#else
    ioun = open_namelist_file()
    read(ioun,nml=ocean_obs_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'ocean_obs_nml')
    call close_file(ioun)    
#endif

    time_window(:) = set_time(0,data_window)

    nfiles=0
    
    if (file_exist('ocean_obs_table') ) then
        call mpp_open(unit,'ocean_obs_table',action=MPP_RDONLY)
        nfiles = 0;nrecs=0
        do while (nfiles <= max_files)
           read(unit,'(a)',end=99,err=98) record
           nrecs=nrecs+1
           if (record(1:1) == '#') cycle
           read(record,*,err=98,end=98) tbl_entry
           nfiles=nfiles+1       
           input_files(nfiles) = tbl_entry%filename
           select case (trim(tbl_entry%file_type))
           case ('profiles')
               filetype(nfiles)    = PROFILE_FILE
           case ('sfc')
               filetype(nfiles)    = SFC_FILE
           case ('idealized')
               filetype(nfiles)    = IDEALIZED_PROFILES
           case default
               call mpp_error(FATAL,'error in obs_table entry format')
           end select
98         continue
        enddo
        call mpp_error(FATAL,' number of obs files exceeds max_files parameter')
99      continue

    endif
    num_profiles=0
    num_sfc_obs=0
    
! get local indices for Model grid
! Since we currently only support regular grids, the
! input 2-d grid array is converted to 1-d
! halo points are added

!    xhalo=isc-isd
!    yhalo=jsc-jsd
!    if (xhalo.ne.ied-iec) call mpp_error(FATAL)
!    if (yhalo.ne.jed-jec) call mpp_error(FATAL)    

    allocate(x_grid(isg-1:ieg+1))
    allocate(y_grid(jsg-1:jeg+1))


    x_grid(isg:ieg) = Grd%x(isg:ieg,jsg)
    y_grid(jsg:jeg) = Grd%y(ieg/4,jsg:jeg)

    allocate(Profiles(max_profiles))

     if (Grd%cyclic) then
         x_grid(isg-1) = x_grid(ieg) - 360. ! assume grid is modulo 360 which is reasonable for data assimilation
         x_grid(ieg+1) = x_grid(isg) + 360.
     else
         x_grid(isg-1) = x_grid(isg) - 1.e-10
         x_grid(ieg+1) = x_grid(ieg) + 1.e-10
     endif

     y_grid(jsg-1) = y_grid(jsg) - 1.e-10
     y_grid(jeg+1) = y_grid(jeg) + 1.e-10
     
    
    do n=1, nfiles
       select case (filetype(n))
       case (PROFILE_FILE)
           call open_profile_dataset(trim(input_files(n)), localize)
       case (IDEALIZED_PROFILES)
           call create_ideal_profiles(localize)
       case default
          call mpp_error(FATAL,'filetype not currently supported')
       end select
    enddo


    return

  end subroutine init_observations
  
   subroutine add_tidal_error(Prof)
! NOT IMPLEMENTED YET !!!
     type(ocean_profile_type), intent(inout) :: Prof

     integer :: k
     real :: dtdz, err, a1, a2
    
     if (.not.ASSOCIATED(prof%ms_t)) then
         allocate(prof%ms_t(prof%levels))
         prof%ms_t(:) = min_obs_err_t
     endif
     
     do k=2,prof%levels - 1
        if (prof%flag_t(k-1) .eq. 0 .and. prof%flag_t(k+1) .eq. 0) then 
            dtdz = (prof%data_t(k+1)-prof%data_t(k-1))/(prof%depth(k+1)-prof%depth(k-1))
            a1 = abs(dtdz) * eta_tide_const
            err = max(a1,min_obs_err_t)
            prof%ms_t(k) = err*err
            if (ASSOCIATED(prof%data_s)) then
                dtdz = (prof%data_s(k+1)-prof%data_s(k-1))/(prof%depth(k+1)-prof%depth(k-1))
                a1 = abs(dtdz) * eta_tide_const
                err = max(a1,min_obs_err_s)
                prof%ms_s(k) = err*err
            endif
        endif
     enddo

   end subroutine add_tidal_error

  subroutine create_ideal_profiles(localize)
!
    use field_manager_mod, only: MODEL_OCEAN, parse, find_field_index, get_field_methods, method_type, get_field_info
    
    logical, intent(in), optional :: localize
    logical :: localize_data = .true.
    integer, parameter :: nlevels = 100 ! number of vertical levels for idealized profiles
    real, parameter :: width_trans = 250.0 ! with over which to transition from sfc value to bottom value
    real, parameter :: bot_depth = 2000.0 ! bottom depth for idealized profiles
    real, allocatable, dimension(:) :: lon,lat, sst, sss, bot_temp, bot_salt, depth
    real, allocatable, dimension(:,:) :: temp, salt, temp_error, salt_error
    integer, allocatable, dimension(:) :: yr, mon, day, hr, mm, ss
    integer :: nstation, unit, n, noobs, i, k
    real :: ri0,rj0,rk0, mid_depth, dtdf, temp_cent, depthC_I_trans, dsdf, salt_cent
    type(time_type) :: profile_time
    logical :: data_is_local
    integer :: model, parse_ok, cpe
    integer :: i0,j0,k0
    real :: dz, a, dx1, dx2, dy1, dy2
    
    real :: temp_missing=missing_value,salt_missing=missing_value,depth_missing=missing_value
    character(len=32) :: fld_type, fld_name
    type(method_type), allocatable, dimension(:) :: ocean_obs_methods
    
    if (PRESENT(localize)) localize_data = localize


    cpe = mpp_pe()

    dx1 = (x_grid(isc)-x_grid(isc-1))/2.0
    dx2 = (x_grid(iec+1)-x_grid(iec))/2.0
    dy1 = (y_grid(jsc)-y_grid(jsc-1))/2.0
    dy2 = (y_grid(jec+1)-y_grid(jec))/2.0
    
    model = model_ocean
    n = find_field_index(model,'ideal_profiles')
    call get_field_info(n,fld_type,fld_name,model,noobs)

    allocate(ocean_obs_methods(noobs))
    allocate(lon(noobs),lat(noobs), yr(noobs), mon(noobs), day(noobs), &
         hr(noobs), mm(noobs), ss(noobs), &
         sst(noobs), sss(noobs), bot_temp(noobs), bot_salt(noobs))
    allocate(temp(noobs,nlevels), salt(noobs,nlevels), temp_error(noobs,nlevels), salt_error(noobs,nlevels))
    allocate(depth(nlevels))
    
    call get_field_methods(n,ocean_obs_methods)
    do i=1,noobs
       parse_ok = parse(ocean_obs_methods(i)%method_control,'lon',lon(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       if (lon(i) .lt. x_grid(isg) ) lon(i) = lon(i) + 360.0
       if (lon(i) .gt. x_grid(ieg) ) lon(i) = lon(i) - 360.0       
       parse_ok = parse(ocean_obs_methods(i)%method_control,'lat',lat(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'yr',yr(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')              
       parse_ok = parse(ocean_obs_methods(i)%method_control,'mon',mon(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'day',day(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'hr',hr(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'mm',mm(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'ss',ss(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'sst',sst(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'sss',sss(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'bot_temp',bot_temp(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
       parse_ok = parse(ocean_obs_methods(i)%method_control,'bot_salt',bot_salt(i))
       if (parse_ok == 0) call mpp_error(FATAL,'==>Error oda_core_mod: idealized_ocean_profiles table entry error')
    enddo

    if (noobs == 0 ) then
        call mpp_error(FATAL,'==> NOTE from oda_core_mod: no idealized profiles given in field table')
        return
    endif

    dz = bot_depth/(nlevels-1)
    
    do k=1,nlevels
       depth(k) = (k-1)*dz
    enddo
    
    mid_depth = bot_depth/2.0
    depthC_I_trans = mid_depth / width_trans
    do i=1,noobs
       dtdf = (bot_temp(i) - sst(i)) / (2.0*atan(1.0) + atan(depthC_I_trans))
       temp_cent = sst(i) + dtdf * atan(depthC_I_trans)
       temp(i,1) = sst(i)
       do k=2,nlevels-1
          temp(i,k) = temp_cent + dtdf * atan((depth(k)  - mid_depth)/width_trans)
       enddo
       temp(i,nlevels) = bot_temp(i)

       dsdf = (bot_salt(i) - sss(i)) / (2.0*atan(1.0) + atan(depthC_I_trans))
       salt_cent = sss(i) + dsdf * atan(depthC_I_trans)
       salt(i,1) = sss(i)
       do k=2,nlevels-1
          salt(i,k) = salt_cent + dsdf * atan((depth(k)  - mid_depth)/width_trans)
       enddo
       salt(i,nlevels) = bot_salt(i)
    enddo


    num_profiles=0
    do i=1,noobs

       data_is_local = .false.


! localized data is within region bounded by halo points
! (halo size = 1) adjacent to boundary points of computational domain
       
       if (lon(i) >= x_grid(isc-1) .and. &
            lon(i) <  x_grid(iec+1) .and. &
            lat(i) >= y_grid(jsc-1) .and. &
            lat(i) <  y_grid(jec+1)) data_is_local = .true.

      
       profile_time = set_date(yr(i),mon(i),day(i),hr(i),mm(i),ss(i))
       
       if ( data_is_local .OR. .NOT.localize_data) then
           if (lon(i) >= x_grid(isc)-dx1 .and. &
                lon(i) < x_grid(iec)+dx2 .and. &
                lat(i) >= y_grid(jsc)-dy1 .and. &
                lat(i) <  y_grid(jec)+dy2) then
!               nprof_in_comp_domain(cpe) = nprof_in_comp_domain(cpe)+1
           endif
           num_profiles=num_profiles+1
           if (num_profiles > max_profiles) then
               call mpp_error(FATAL,'maximum number of profiles exceeded.  Resize parameter max_profiles in ocean_obs_mod')
           endif


           profiles(num_profiles)%Model_Grid => Grd
           profiles(num_profiles)%nvar = 2
           profiles(num_profiles)%profile_flag = 0
           profiles(num_profiles)%profile_flag_s = 0
           profiles(num_profiles)%accepted = .true.
           allocate(profiles(num_profiles)%depth(nlevels))
           profiles(num_profiles)%depth=depth(1:nlevels)
           allocate(profiles(num_profiles)%data_t(nlevels))
           profiles(num_profiles)%data_t=temp(i,:)
           allocate(profiles(num_profiles)%flag_t(nlevels))
           profiles(num_profiles)%flag_t= 0               

           allocate(profiles(num_profiles)%data_s(nlevels))
           profiles(num_profiles)%data_s=salt(i,:)
           allocate(profiles(num_profiles)%flag_s(nlevels))
           profiles(num_profiles)%flag_s= 0

           profiles(num_profiles)%probe = 0.
           profiles(num_profiles)%levels = nlevels
           profiles(num_profiles)%lat = lat(i)
           profiles(num_profiles)%lon = lon(i)
           allocate(profiles(num_profiles)%ms_t(nlevels))
           profiles(num_profiles)%ms_t(:) = min_obs_err_t ! default error variance for temperature
           allocate(profiles(num_profiles)%ms_s(nlevels))
           profiles(num_profiles)%ms_s(:) = min_obs_err_s  ! default error variance for salinity

           profiles(num_profiles)%time = profile_time
           
! calculate interpolation coefficients (make sure to account for grid offsets here!)
! note that this only works for lat/lon grids
           
           ri0 = frac_index(lon(i), x_grid(isg-1:ieg+1)) - 1.
           rj0 = frac_index(lat(i), y_grid(jsg-1:jeg+1)) - 1.
           i0 = floor(ri0)
           j0 = floor(rj0)
           Profiles(num_profiles)%i_index = ri0
           Profiles(num_profiles)%j_index = rj0           
           Profiles(num_profiles)%accepted = .true.
           if (i0 < 0 .or. j0 < 0) then
               Profiles(num_profiles)%accepted = .false.
           endif
           if (i0 > ieg+1 .or. j0 > jeg+1) then
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
           if (i0 < isc-1 .or. i0 > iec) then
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
           if (j0 < jsc-1 .or. j0 > jec) then
               call mpp_error(FATAL,'grid lat/lon index is out of bounds ')
           endif
           if (Profiles(num_profiles)%accepted ) then
               if (Grd%mask(i0,j0,1) == 0.0 .or. &
                    Grd%mask(i0+1,j0,1) == 0.0 .or. &
                    Grd%mask(i0,j0+1,1) == 0.0 .or. &
                    Grd%mask(i0+1,j0+1,1) == 0.0) then
                   Profiles(num_profiles)%accepted = .false.
               endif
           endif


           
           if (Profiles(num_profiles)%accepted) then
               allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
               do k=1, Profiles(num_profiles)%levels
                  rk0 = frac_index(Profiles(num_profiles)%depth(k), Grd%z(:))
                  k0 = floor(rk0)
                  if ( k0 == -1) then
                      if (Profiles(num_profiles)%depth(k) .ne. missing_value .and. &
                           Profiles(num_profiles)%depth(k) .lt. Grd%z(1)) then
                           k0 = 1
                           rk0 = 1.0
                       endif
                   endif


                  if (k0 .gt. size(Grd%z)-1 ) then
                      write(*,*) 'k0 out of bounds, rk0,k0= ',rk0,k0
                     write(*,*) 'Z_bound= ',Grd%z_bound
                      write(*,*) 'Profile%depth= ',Profiles(num_profiles)%depth
                      call mpp_error(FATAL)
                  endif
                  
                  Profiles(num_profiles)%k_index(k) = rk0
                  
                  if (Profiles(num_profiles)%flag_t(k) .eq. 0) then
                      if (Grd%mask(i0,j0,k0) == 0.0 .or. &
                          Grd%mask(i0+1,j0,k0) == 0.0 .or. &
                          Grd%mask(i0,j0+1,k0) == 0.0 .or. &
                          Grd%mask(i0+1,j0+1,k0) == 0.0) then
                          Profiles(num_profiles)%flag_t(k) = 1
                      endif
                      if (Grd%mask(i0,j0,k0+1) == 0.0 .or. &
                          Grd%mask(i0+1,j0,k0+1) == 0.0 .or. &
                          Grd%mask(i0,j0+1,k0+1) == 0.0 .or. &
                          Grd%mask(i0+1,j0+1,k0+1) == 0.0) then
                          Profiles(num_profiles)%flag_t(k) = 1
                      endif
                      if (Profiles(num_profiles)%data_t(k) == missing_value &
                         .or. Profiles(num_profiles)%depth(k) == missing_value) then
                          Profiles(num_profiles)%flag_t(k) = 1
                      endif
                  endif
                  
               enddo
           endif          
       endif

    enddo

!    a = nprof_in_comp_domain(cpe)

!    call mpp_broadcast(nprof_in_comp_domain(cpe),cpe)
    
!    call mpp_sum(a)

!    write(stdout(),*) 'A grand total of ',a,' profiles satisify acceptance criteria'

!    do i=0,mpp_npes()-1
!       write(stdout(),*) 'pe=',i,'profile count=',nprof_in_comp_domain(i)
!    enddo
    
  end subroutine create_ideal_profiles


  subroutine nullify_obs_prof(profile)

    type(ocean_profile_type), intent(inout) :: profile


    profile%nvar = 0
    profile%project=0.
    profile%probe=0.
    profile%ref_inst=0.
    profile%wod_cast_num=0
    profile%fix_depth=0.
    profile%ocn_vehicle=0.
    profile%database_id=0.
    profile%levels=0
    profile%profile_flag=-1
    profile%profile_flag_s=-1
    profile%lon=-1.0e10
    profile%lat=-1.0e10
    profile%accepted=.false.
    profile%nlinks=0
    if (ASSOCIATED(profile%next)) profile%next=>NULL()
    if (ASSOCIATED(profile%depth)) profile%depth=>NULL()
    if (ASSOCIATED(profile%data_t)) profile%data_t=>NULL()
    if (ASSOCIATED(profile%data_s)) profile%data_s=>NULL()
    if (ASSOCIATED(profile%flag_t)) profile%flag_t=>NULL()
    if (ASSOCIATED(profile%flag_s)) profile%flag_s=>NULL()
    profile%temp_err=0.0
    profile%salt_err=0.0
    if (ASSOCIATED(profile%ms_t)) profile%ms_t=>NULL()
    if (ASSOCIATED(profile%ms_s)) profile%ms_s=>NULL()    
    profile%time = set_time(0,0)
    profile%yyyy = 0
    profile%mmdd = 0
    if (ASSOCIATED(profile%model_time)) profile%model_time=>NULL()
    if (ASSOCIATED(profile%model_grid)) profile%model_grid=>NULL()
    if (ASSOCIATED(profile%k_index)) profile%k_index=>NULL()
    profile%i_index=-1.0
    profile%j_index=-1.0
    profile%tdiff = set_time(0,0)
    
  end subroutine nullify_obs_prof
  
end module oda_core_mod
