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

! modified by JHC
!> Retrieve sea surface temperature data and interpolated grid
subroutine GET_AMIP_SST_ (Time, Interp, sst, err_msg, lon_model, lat_model)
   type (time_type),         intent(in)    :: Time !< Time to interpolate
   type (AMIP_INTERP_TYPE_),  intent(inout) :: Interp !< Holds data for interpolation
   real(AMIP_INTERP_KIND_),                     intent(out)   ::  sst(:,:) !< Sea surface temperature data
   character(len=*), optional, intent(out) :: err_msg !< Holds error message string if present

   real(AMIP_INTERP_KIND_), dimension(mobs,nobs) :: sice
   real(AMIP_INTERP_KIND_), allocatable:: temp1(:,:), temp2(:,:)

    integer :: year1, year2, month1, month2
    real(AMIP_INTERP_KIND_)    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum(3)

! add by JHC
    real(AMIP_INTERP_KIND_),    intent(in), dimension(:,:), optional :: lon_model, lat_model
    real(AMIP_INTERP_KIND_) :: pert
    integer :: i, j, mobs_sst, nobs_sst
    integer :: jhctod(6)
    type (time_type) :: Udate
    character(len=4) :: yyyy
    integer :: nrecords, ierr, k, yr, mo, dy
    integer, dimension(:), allocatable :: ryr, rmo, rdy
    character(len=30) :: time_unit
    real(AMIP_INTERP_KIND_), dimension(:), allocatable :: timeval
    character(len=maxc) :: ncfilename
    type(FmsNetcdfFile_t) :: fileobj
    logical :: the_file_exists
! end add by JHC
    logical, parameter :: DEBUG = .false. !> switch for debugging output
    !> These are fms_io specific
    integer :: unit
    integer, parameter :: lkind = AMIP_INTERP_KIND_

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
      call ZONAL_SST_ (Amip_Time, sice, temp1)
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
              call READ_RECORD_ ('sst', Date1, Udate1, temp1)
              if ( use_ncep_sst .and. (.not. no_anom_sst) ) then
                   temp1(:,:) = temp1(:,:) + SST_ANOM_ (:,:)
              endif
              call horiz_interp ( Interp%Hintrp, temp1, Interp%data1 )
              call CLIP_DATA_ ('sst', Interp%data1)
             Interp % Date1 = Date1
          endif
      endif

!-----------------------------------------------------------------------

      if (Date2 /= Interp % Date2) then
          call READ_RECORD_ ('sst', Date2, Udate2, temp2)
          if ( use_ncep_sst .and. (.not. no_anom_sst) ) then
               temp2(:,:) = temp2(:,:) + SST_ANOM_ (:,:)
          endif
          call horiz_interp ( Interp%Hintrp, temp2, Interp%data2 )
          call CLIP_DATA_ ('sst', Interp%data2)
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
         SST_ANOM_ (:,:) = SST_NCEP_ (:,:) - (temp1(:,:) + fmonth*(temp2(:,:) - temp1(:,:)) )
         call horiz_interp ( Interp%Hintrp, SST_NCEP_, sst )
         call CLIP_DATA_ ('sst', sst)
    endif

!! DEBUG CODE
    if (DEBUG) then
          call get_date(Amip_Time,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
          if (mpp_pe() == 0) then
             write (*,200) 'JHC: use_daily = F, AMIP_Time: ',jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5), &
                   & jhctod(6)
             write (*,300) 'JHC: use_daily = F, interped SST: ', sst(1,1),sst(5,5),sst(10,10)
          endif
    endif


  endif

! add by JHC
else
    call get_date(Amip_Time,jhctod(1),jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6))
     if (mpp_pe() == mpp_root_pe()) write(*,200) 'amip_interp_mod: use_daily = T, Amip_Time = ',jhctod(1), &
        & jhctod(2),jhctod(3),jhctod(4),jhctod(5),jhctod(6)

    yr = jhctod(1); mo = jhctod(2); dy = jhctod(3)

    write (yyyy,'(i4)') jhctod(1)

    file_name_sst = 'INPUT/' // 'sst.day.mean.'//yyyy//'.v2.nc'
    ncfilename = trim(file_name_sst)
    time_unit = 'days since 1978-01-01 00:00:00'

    mobs_sst = 1440;  nobs_sst = 720

    call SET_SST_GRID_EDGES_DAILY_ (mobs_sst, nobs_sst)
    call horiz_interp_new ( Interp%Hintrp2, LON_BND_, LAT_BND_, &
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
     if ( .not. allocated(TEMPAMIP_) ) &
     & allocate (TEMPAMIP_ (mobs_sst,nobs_sst))

     if (the_file_exists) then
          call fms2_io_read_data(fileobj, 'SST', TEMPAMIP_, unlim_dim_level=k)
          call close_file(fileobj)
          TEMPAMIP_ = TEMPAMIP_ + real(TFREEZE, AMIP_INTERP_KIND_)

!!! DEBUG CODE
          if(DEBUG) then
            if (mpp_pe() == 0) then
              print*, 'JHC: TFREEZE = ', real(TFREEZE, AMIP_INTERP_KIND_)
              print*, lbound(sst)
              print*, ubound(sst)
              print*, lbound(TEMPAMIP_)
              print*, ubound(TEMPAMIP_)
              write(*,300) 'JHC: tempamip : ', TEMPAMIP_ (100,100), TEMPAMIP_ (200,200), TEMPAMIP_ (300,300)
            endif
          endif

          call horiz_interp ( Interp%Hintrp2, TEMPAMIP_, sst )
          call CLIP_DATA_ ('sst', sst)

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
          sst = sst + real(sst_pert, AMIP_INTERP_KIND_)
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
             sst (i,j) = sst (i,j) + real(sst_pert, AMIP_INTERP_KIND_)*((pert-0.5_lkind)*2)
          end do
          end do
      endif

  endif
! end add by JHC
 end subroutine GET_AMIP_SST_

!> AMIP interpolation for ice
subroutine GET_AMIP_ICE_ (Time, Interp, ice, err_msg)
   type (time_type),         intent(in)    :: Time !< Time to interpolate
   type (AMIP_INTERP_TYPE_),  intent(inout) :: Interp !< Holds data for interpolation
   real(AMIP_INTERP_KIND_),                     intent(out)   :: ice(:,:) !< ice data
   character(len=*), optional, intent(out) :: err_msg !< Holds error message string if present

    real(AMIP_INTERP_KIND_), dimension(mobs,nobs) :: sice, temp

    integer :: year1, year2, month1, month2
    real(AMIP_INTERP_KIND_)    :: fmonth
    type (date_type) :: Date1, Date2, Udate1, Udate2

    type(time_type) :: Amip_Time
    integer :: tod(3),dum(3)
    integer, parameter :: lkind = AMIP_INTERP_KIND_

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
   call ZONAL_SST_ (Amip_Time, sice, temp)
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
               where ( SST_NCEP_ <= (real(TFREEZE, AMIP_INTERP_KIND_)+real(tice_crit, AMIP_INTERP_KIND_)) )
                   sice = 1._lkind
               elsewhere
                   sice = 0._lkind
               endwhere
            else
               call READ_RECORD_ ('ice', Date1, Udate1, sice)
            endif
!--------------------------------------------------------------------
            call horiz_interp ( Interp%Hintrp, sice, Interp%data1 )
            call CLIP_DATA_ ('ice', Interp%data1)
            Interp % Date1 = Date1
        endif
    endif

!-----------------------------------------------------------------------

    if (Date2 /= Interp % Date2) then

!-- SJL -------------------------------------------------------------
            if ( use_ncep_sst .and. use_ncep_ice ) then
               where ( SST_NCEP_ <= (real(TFREEZE, AMIP_INTERP_KIND_)+real(tice_crit, AMIP_INTERP_KIND_)) )
                   sice = 1._lkind
               elsewhere
                   sice = 0._lkind
               endwhere
            else
               call READ_RECORD_ ('ice', Date2, Udate2, sice)
            endif
!--------------------------------------------------------------------
        call horiz_interp ( Interp%Hintrp, sice, Interp%data2 )
        call CLIP_DATA_ ('ice', Interp%data2)
        Interp % Date2 = Date2

    endif

!-----------------------------------------------------------------------
!---------- time interpolation (between months) ------------------------
!-----------------------------------------------------------------------

   ice = Interp % data1 + fmonth * (Interp % data2 - Interp % data1)

endif
 end subroutine GET_AMIP_ICE_

 !> @return A newly created @ref amip_interp_type
 function AMIP_INTERP_NEW_1D_ ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)
 real(AMIP_INTERP_KIND_),    intent(in), dimension(:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional       :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (AMIP_INTERP_TYPE_) :: Interp

   if(.not.module_is_initialized) call AMIP_INTERP_INIT_

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

    call horiz_interp_new ( Interp%Hintrp, LON_BND_, LAT_BND_, &
                             lon, lat, interp_method= interp_method )

    allocate(Interp % data1 (size(lon(:))-1,size(lat(:))-1))
    allocate(Interp % data2 (size(lon(:))-1,size(lat(:))-1))

    Interp%I_am_initialized = .true.
   end function AMIP_INTERP_NEW_1D_

 !> @return A newly created @ref amip_interp_type
 function AMIP_INTERP_NEW_2D_ ( lon , lat , mask , use_climo, use_annual, &
                                interp_method ) result (Interp)
 real(AMIP_INTERP_KIND_),    intent(in), dimension(:,:)   :: lon, lat
 logical, intent(in), dimension(:,:) :: mask
 character(len=*), intent(in), optional :: interp_method
 logical, intent(in), optional       :: use_climo, use_annual

   type (AMIP_INTERP_TYPE_) :: Interp

   if(.not.module_is_initialized) call AMIP_INTERP_INIT_

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

   call horiz_interp_new ( Interp%Hintrp, LON_BND_, LAT_BND_, &
                           lon, lat, interp_method = interp_method)

   allocate(Interp % data1 (size(lon,1),size(lat,2)))
   allocate(Interp % data2 (size(lon,1),size(lat,2)))

   Interp%I_am_initialized = .true.
   end function AMIP_INTERP_NEW_2D_

 !> initialize @ref amip_interp_mod for use
 subroutine AMIP_INTERP_INIT_
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
        call SET_SST_GRID_EDGES_AMIP1_
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 1 sst', NOTE)
        Date_end = date_type( 1989, 1, 0 )
    else if (lowercase(trim(data_set)) == 'amip2') then
        file_name_sst = 'INPUT/' // 'amip2_sst.data'
        file_name_ice = 'INPUT/' // 'amip2_ice.data'
        mobs = 360;  nobs = 180
        call SET_SST_GRID_EDGES_OI_
!       --- specfied min for amip2 ---
        tice_crit_k = 271.38_r8_kind
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AMIP 2 sst', NOTE)
        Date_end = date_type( 1996, 3, 0 )
    else if (lowercase(trim(data_set)) == 'hurrell') then
        file_name_sst = 'INPUT/' // 'hurrell_sst.data'
        file_name_ice = 'INPUT/' // 'hurrell_ice.data'
        mobs = 360;  nobs = 180
        call SET_SST_GRID_EDGES_OI_
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
        call SET_SST_GRID_EDGES_OI_
        if (mpp_pe() == 0) &
        call error_mesg ('amip_interp_init', 'using AVHRR daily sst', NOTE)
        Date_end = date_type( 2011, 8, 16 )
! end add by JHC
    else if (lowercase(trim(data_set)) == 'reynolds_eof') then
        file_name_sst = 'INPUT/' // 'reynolds_sst.data'
        file_name_ice = 'INPUT/' // 'reynolds_sst.data'
        mobs = 180;  nobs = 90
        call SET_SST_GRID_EDGES_OI_
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
            if (.not. allocated (SST_NCEP_)) then
                allocate (SST_NCEP_ (i_sst,j_sst))
                SST_NCEP_ (:,:) = real(big_number, AMIP_INTERP_KIND_)
            endif
            if (.not. allocated (SST_ANOM_)) then
                allocate (SST_ANOM_ (i_sst,j_sst))
                SST_ANOM_ (:,:) = real(big_number, AMIP_INTERP_KIND_)
            endif
        else
             mobs = 360;    nobs = 180
        endif
!--- Added by SJL ----------------------------------------------
        call SET_SST_GRID_EDGES_OI_
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
 end subroutine AMIP_INTERP_INIT_

!> Frees data associated with a amip_interp_type variable. Should be used for any
!! variables initialized via @ref amip_interp_new.
!> @param[inout] Interp A defined data type variable initialized by amip_interp_new and used
!! when calling get_amip_sst and get_amip_ice.
   subroutine AMIP_INTERP_DEL_ (Interp)
   type (AMIP_INTERP_TYPE_), intent(inout) :: Interp
     if(allocated(Interp%data1)) deallocate(Interp%data1)
     if(allocated(Interp%data2)) deallocate(Interp%data2)
     if(allocated(LON_BND_))   deallocate(LON_BND_)
     if(allocated(LAT_BND_))   deallocate(LAT_BND_)
     call horiz_interp_del ( Interp%Hintrp )

     Interp%I_am_initialized = .false.
   end subroutine AMIP_INTERP_DEL_

   subroutine SET_SST_GRID_EDGES_AMIP1_
   integer :: i, j
   real(AMIP_INTERP_KIND_) :: hpie, dlon, dlat, wb, sb
   integer, parameter :: lkind = AMIP_INTERP_KIND_

      allocate(LON_BND_ (mobs+1))
      allocate(LAT_BND_ (nobs+1))

! ---- compute grid edges (do only once) -----

      hpie = 0.5_lkind*real(pi, AMIP_INTERP_KIND_)

      dlon = 4._lkind*hpie/real(mobs, AMIP_INTERP_KIND_);  wb = -0.5_lkind*dlon
      do i = 1, mobs+1
          LON_BND_ (i) = wb + dlon * real(i-1, AMIP_INTERP_KIND_)
      enddo
          LON_BND_ (mobs+1) = LON_BND_ (1) + 4._lkind*hpie

      dlat = 2._lkind*hpie/real(nobs-1, AMIP_INTERP_KIND_);  sb = -hpie + 0.5_lkind*dlat
      LAT_BND_ (1) = -hpie;  LAT_BND_ (nobs+1) = hpie
      do j = 2, nobs
          LAT_BND_ (j) = sb + dlat * real(j-2, AMIP_INTERP_KIND_)
      enddo
   end subroutine SET_SST_GRID_EDGES_AMIP1_

   subroutine SET_SST_GRID_EDGES_OI_
   integer :: i, j
   real(AMIP_INTERP_KIND_)    :: hpie, dlon, dlat, wb, sb
   integer, parameter :: lkind = AMIP_INTERP_KIND_

! add by JHC
      if(allocated(LON_BND_))       deallocate(LON_BND_)
      if(allocated(LAT_BND_))       deallocate(LAT_BND_)
! end add by JHC
      allocate(LON_BND_ (mobs+1))
      allocate(LAT_BND_ (nobs+1))

! ---- compute grid edges (do only once) -----

      hpie = 0.5_lkind*real(pi, AMIP_INTERP_KIND_)

      dlon = 4._lkind*hpie/real(mobs, AMIP_INTERP_KIND_);  wb = 0.0_lkind
          LON_BND_ (1) = wb
      do i = 2, mobs+1
          LON_BND_ (i) = wb + dlon * real(i-1, AMIP_INTERP_KIND_)
      enddo
          LON_BND_ (mobs+1) = LON_BND_ (1) + 4._lkind*hpie

      dlat = 2._lkind*hpie/real(nobs, AMIP_INTERP_KIND_);  sb = -hpie
      LAT_BND_ (1) = sb;  LAT_BND_ (nobs+1) = hpie
      do j = 2, nobs
          LAT_BND_ (j) = sb + dlat * real(j-1, AMIP_INTERP_KIND_)
      enddo
   end subroutine SET_SST_GRID_EDGES_OI_

! add by JHC
   subroutine SET_SST_GRID_EDGES_DAILY_ (mobs_sst, nobs_sst)
   integer :: i, j, mobs_sst, nobs_sst
   real(AMIP_INTERP_KIND_)    :: hpie, dlon, dlat, wb, sb
   integer, parameter :: lkind = AMIP_INTERP_KIND_

      if(allocated(LON_BND_))       deallocate(LON_BND_)
      if(allocated(LAT_BND_))       deallocate(LAT_BND_)
      allocate(LON_BND_ (mobs_sst+1))
      allocate(LAT_BND_ (nobs_sst+1))

! ---- compute grid edges (do only once) -----

      hpie = 0.5_lkind*real(pi, AMIP_INTERP_KIND_)

      dlon = 4._lkind*hpie/real(mobs_sst, AMIP_INTERP_KIND_);  wb = 0.0_lkind
          LON_BND_ (1) = wb
      do i = 2, mobs_sst+1
          LON_BND_ (i) = wb + dlon * real(i-1, AMIP_INTERP_KIND_)
      enddo
          LON_BND_ (mobs_sst+1) = LON_BND_ (1) + 4._lkind*hpie

      dlat = 2._lkind*hpie/real(nobs_sst, AMIP_INTERP_KIND_);  sb = -hpie
      LAT_BND_ (1) = sb;  LAT_BND_ (nobs_sst+1) = hpie
      do j = 2, nobs_sst
          LAT_BND_ (j) = sb + dlat * real(j-1, AMIP_INTERP_KIND_)
      enddo
   end subroutine SET_SST_GRID_EDGES_DAILY_
! end add by JHC

   subroutine A2A_BILINEAR_ (nx, ny, dat1, n1, n2, dat2)
   integer, intent(in):: nx, ny
   integer, intent(in):: n1, n2
   real(AMIP_INTERP_KIND_), intent(in) :: dat1(nx,ny)
   real(AMIP_INTERP_KIND_), intent(out):: dat2(n1,n2)      !> output interpolated data

! local:
  real(AMIP_INTERP_KIND_):: lon1(nx), lat1(ny)
  real(AMIP_INTERP_KIND_):: lon2(n1), lat2(n2)
  real(AMIP_INTERP_KIND_):: dx1, dy1, dx2, dy2
  real(AMIP_INTERP_KIND_):: xc, yc
  real(AMIP_INTERP_KIND_):: a1, b1, c1, c2, c3, c4
  integer i1, i2, jc, i0, j0, it, jt
  integer i,j
  integer, parameter :: lkind = AMIP_INTERP_KIND_


!-----------------------------------------------------------
! * Interpolate from "FMS" 1x1 SST data grid to a finer grid
!                     lon: 0.5, 1.5, ..., 359.5
!                     lat: -89.5, -88.5, ... , 88.5, 89.5
!-----------------------------------------------------------

  dx1 = 360._lkind/real(nx, AMIP_INTERP_KIND_) !> INput Grid
  dy1 = 180._lkind/real(ny, AMIP_INTERP_KIND_) !> INput Grid

  do i=1,nx
     lon1(i) = 0.5_lkind*dx1 + real(i-1, AMIP_INTERP_KIND_)*dx1
  enddo
  do j=1,ny
     lat1(j) = -90._lkind + 0.5_lkind*dy1 + real(j-1, AMIP_INTERP_KIND_)*dy1
  enddo

  dx2 = 360._lkind/real(n1, AMIP_INTERP_KIND_) !> OutPut Grid:
  dy2 = 180._lkind/real(n2, AMIP_INTERP_KIND_) !> OutPut Grid:

  do i=1,n1
     lon2(i) = 0.5_lkind*dx2 + real(i-1, AMIP_INTERP_KIND_)*dx2
  enddo
  do j=1,n2
     lat2(j) = -90._lkind + 0.5_lkind*dy2 + real(j-1, AMIP_INTERP_KIND_)*dy2
  enddo

  jt = 1
  do 5000 j=1,n2

     yc = lat2(j)
     if ( yc<lat1(1) ) then
            jc = 1
            b1 = 0._lkind
     elseif ( yc>lat1(ny) ) then
            jc = ny-1
            b1 = 1._lkind
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
            a1 = (xc+360._lkind-lon1(nx)) / dx1
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
       if ( a1<-0.001_lkind .or. a1>1.001_lkind .or.  b1<-0.001_lkind .or. b1>1.001_lkind ) then
            write(*,*) i,j,a1, b1
            call mpp_error(FATAL,'a2a bilinear interpolation')
       endif

       c1 = (1._lkind-a1) * (1._lkind-b1)
       c2 =     a1  * (1._lkind-b1)
       c3 =     a1  *     b1
       c4 = (1._lkind-a1) *     b1

! Bilinear interpolation:
       dat2(i,j) = c1*dat1(i1,jc) + c2*dat1(i2,jc) + c3*dat1(i2,jc+1) + c4*dat1(i1,jc+1)

     enddo   !i-loop

5000 continue   ! j-loop
   end subroutine A2A_BILINEAR_

   subroutine READ_RECORD_ (type, Date, Adate, dat)
     character(len=*), intent(in)  :: type
     type (date_type), intent(in)  :: Date
     type (date_type), intent(inout) :: Adate
     real(AMIP_INTERP_KIND_), intent(out) :: dat(mobs,nobs)
     real(AMIP_INTERP_KIND_) :: tmp_dat(360,180)

     integer(I2_KIND) :: idat(mobs,nobs)
     integer :: nrecords, yr, mo, dy, ierr, k
     integer, dimension(:), allocatable :: ryr, rmo, rdy
     character(len=maxc) :: ncfilename, ncfieldname
     type(FmsNetcdfFile_t), pointer :: fileobj
     integer, parameter :: lkind = AMIP_INTERP_KIND_

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
               call A2A_BILINEAR_ (360, 180, tmp_dat, mobs, nobs, dat)
          else
               dat(:,:) = tmp_dat(:,:)
          endif
     else
          call fms2_io_read_data(fileobj, ncfieldname, dat, unlim_dim_level=k)
     endif
     !TODO This assumes that the data is "packed" (has the scale_factor and add_offset attributes)
     ! in fms2_io_read_data the data is unpacked (data_in_file*scale_factor + add_offset)
     ! the line below "packs" the data again. This is needed for reproducibility
     idat =  nint(dat*100._lkind, I2_KIND)

   !---- unpacking of data ----

     if (type(1:3) == 'ice') then
        !---- create fractional [0,1] ice mask
        if (lowercase(trim(data_set)) /= 'amip2' .and. lowercase(trim(data_set)) /= 'hurrell') then
               where ( idat <= ice_crit )
                   dat = 1._lkind
               elsewhere
                   dat = 0._lkind
               endwhere
        else
           dat = dat*0.01_lkind
        endif
     else if (type(1:3) == 'sst') then
        !---- unpack sst ----
        if (lowercase(trim(data_set)) /= 'amip2' .and. lowercase(trim(data_set)) /= 'hurrell') then
               dat = real(idat, AMIP_INTERP_KIND_)*0.01_lkind + real(TFREEZE, AMIP_INTERP_KIND_)
        endif
     endif

     return
   end subroutine READ_RECORD_

   subroutine CLIP_DATA_ (type, dat)
   character(len=*), intent(in)    :: type
   real(AMIP_INTERP_KIND_), intent(inout) :: dat(:,:)
   integer, parameter :: lkind = AMIP_INTERP_KIND_

   if (type(1:3) == 'ice') then
       dat = min(max(dat,0.0_lkind), 1.0_lkind)
   else if (type(1:3) == 'sst') then
       dat = max(real(tice_crit_k, AMIP_INTERP_KIND_),dat)
   endif
   end subroutine CLIP_DATA_

subroutine ZONAL_SST_ (Time, ice, sst)
   type (time_type), intent(in)  :: Time
   real(AMIP_INTERP_KIND_), intent(out) :: ice(mobs,nobs), sst(mobs,nobs)
   real(AMIP_INTERP_KIND_) :: tpi, fdate, eps, ph, sph, sph2, ts
   integer :: j
   integer, parameter :: lkind = AMIP_INTERP_KIND_

! namelist needed
!
!  teq  = sst at equator
!  tdif = equator to pole sst difference
!  tann = amplitude of annual cycle
!  tlag = offset for time of year (for annual cycle)
!

    tpi = 2.0_lkind*real(pi, AMIP_INTERP_KIND_)

    fdate = fraction_of_year (Time)

    eps = sin( tpi*(fdate-real(tlag, AMIP_INTERP_KIND_)) ) * real(tann, AMIP_INTERP_KIND_)

    do j = 1, nobs

        ph  = 0.5_lkind*(LAT_BND_ (j)+LAT_BND_ (j+1))
       sph  = sin(ph)
       sph2 = sph*sph

       ts = real(teq, AMIP_INTERP_KIND_) - real(tdif, AMIP_INTERP_KIND_)*sph2 - eps*sph

       sst(:,j) = ts

    enddo

    where ( sst < real(tice_crit_k, AMIP_INTERP_KIND_) )
       ice = 1.0_lkind
       sst = real(tice_crit_k, AMIP_INTERP_KIND_)
    elsewhere
       ice  = 0.0_lkind
    endwhere
end subroutine ZONAL_SST_

subroutine AMIP_INTERP_TYPE_EQ_ (amip_interp_out, amip_interp_in)
    type(AMIP_INTERP_TYPE_), intent(inout) :: amip_interp_out
    type(AMIP_INTERP_TYPE_), intent(in)    :: amip_interp_in

    if(.not.amip_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'amip_interp_type_eq: amip_interp_type variable on right hand side is unassigned')
    endif

    amip_interp_out%Hintrp     =  amip_interp_in%Hintrp
    amip_interp_out%Hintrp2    =  amip_interp_in%Hintrp2 !< missing assignment statement; added by GPP
    amip_interp_out%data1      =  amip_interp_in%data1
    amip_interp_out%data2      =  amip_interp_in%data2
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%Date1      =  amip_interp_in%Date1
    amip_interp_out%Date2      =  amip_interp_in%Date2
    amip_interp_out%use_climo  =  amip_interp_in%use_climo
    amip_interp_out%use_annual =  amip_interp_in%use_annual
    amip_interp_out%I_am_initialized = .true.
end subroutine AMIP_INTERP_TYPE_EQ_
