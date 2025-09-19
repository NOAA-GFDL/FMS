program test_time_interp_with_interp

  use constants_mod, only: DEG_TO_RAD

  use fms_mod, only: fms_init

  use fms2_io_mod, only: FmsNetcdfFile_t, fms2_io_init, open_file, close_file, register_axis, write_data
  use fms2_io_mod, only: register_field, unlimited, register_variable_attribute

  use horiz_interp_mod, only: horiz_interp_type, horiz_interp_new, horiz_interp_init

  use mpp_mod, only: input_nml_file, mpp_pe, mpp_root_pe, mpp_npes
  use mpp_domains_mod, only: mpp_define_layout, mpp_define_domains, domain2d, mpp_get_compute_domain

  use time_interp_external2_mod, only: init_external_field, time_interp_external_init, time_interp_external

  use time_manager_mod, only: JULIAN, time_type, set_date, set_calendar_type, time_manager_init

  implicit none

  integer, parameter :: rkind = TEST_FMS_KIND_

  integer, parameter :: ntime = 3
  integer, parameter :: nlon_model = 90
  integer, parameter :: nlat_model = 46
  integer, parameter :: nlon_file = nlon_model !/ 2
  integer, parameter :: nlat_file = nlat_model !/ 2

  !time_interp data file
  character(len=128) :: filename='INPUT/aerosol2.climatology.nc'
  character(len=128) :: fieldname='so4_anthro'

  !grid
  real(TEST_FMS_KIND_), allocatable :: lon_model(:,:), lat_model(:,:)
  real(TEST_FMS_KIND_), allocatable :: lon_file(:,:), lat_file(:,:)

  !horiz_interp
  type(horiz_interp_type) :: interp

  !domain
  integer :: isc_m, iec_m, jsc_m, jec_m
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  type(domain2d) :: domain_model
  type(domain2d) :: domain_file
  integer :: layout(2), global_indices(4)

  !time
  type(time_type) :: time
  integer :: times(ntime) = (/0,1,2/)

  !time_interp_external
  integer :: field_id, field2_id
  real(TEST_FMS_KIND_), allocatable :: data_out(:,:)
  real(TEST_FMS_KIND_), allocatable :: data_file(:,:)

  !counter
  integer :: i, j

  !initialize fms
  call fms_init()
  call time_interp_external_init()
  call time_manager_init()
  call horiz_interp_init()

  !set up domain_model
  global_indices = (/1, nlon_model, 1, nlat_model/)
  call mpp_define_layout(global_indices, mpp_npes(), layout)
  call mpp_define_domains(global_indices, layout, domain_model)
  call mpp_get_compute_domain(domain_model, isc, iec, jsc, jec)
  call mpp_get_compute_domain(domain_model, isd, ied, jsd, jed)

  !allocate arrays
  allocate(lon_model(isc:iec+1, jsc:jec+1), lat_model(isc:iec+1, jsc:jec+1))
  allocate(lon_file(nlon_file+1, nlat_file+1), lat_file(nlon_file+1, nlat_file+1))
  allocate(data_out(isc:iec, jsc:jec))
  allocate(data_file(nlon_file, nlat_file))

  !set model coordintes
  do j=jsc, jec+1
    do i=isc, iec+1
      lon_model(i,j) = real(i, rkind) * DEG_TO_RAD
      lat_model(i,j) = real(j, rkind) * DEG_TO_RAD
    end do
  end do

  !set file coordinates
  do j=1, nlat_file+1
    do i=1, nlon_file+1
      lon_file(i,j) = real(i, rkind)
      lat_file(i,j) = real(j, rkind)
    end do
  end do

  !set file data
  do j=1, nlat_file
    do i=1, nlon_file
      data_file(i,j) = i*nlon_file + j
    end do
  end do

  !generate_file
  call generate_file(filename, lon_file, lat_file, times, fieldname, data_file)
  lon_file = lon_file * DEG_TO_RAD
  lat_file = lat_file * DEG_TO_RAD

  !set time
  call set_calendar_type(JULIAN)
  time = set_date(1800, 1, 1, 0, 0, 0)

  call horiz_interp_new(interp, lon_file, lat_file, lon_model, lat_model, &
    interp_method="conservative", is_latlon_in=.false., is_latlon_out=.false.)

  !time_interp
  field_id = init_external_field(filename, fieldname, domain=domain_model, override=.true.)

  call time_interp_external(field_id, time, data_out, horz_interp=interp)

  do i=1, 10
    write(*,*) data_out(i,i) - data_file(i,i)
  end do

contains

  subroutine generate_file(filename, lon, lat, times, fieldname, data)

    character(len=*), intent(in) :: filename
    real(TEST_FMS_KIND_), intent(in) :: lon(:,:), lat(:,:)
    integer, intent(in) :: times(:)
    character(len=*), intent(in) :: fieldname
    real(TEST_FMS_KIND_), intent(inout) :: data(:,:)

    type(FmsNetcdfFile_t) :: fileobj !< fileobj
    integer :: nlon, nlat, i, j

    nlon = size(data, 1)
    nlat = size(data, 2)

    !< Create a file to test with:
    if (mpp_pe() .eq. mpp_root_pe()) then
      if (open_file(fileobj, filename, "overwrite")) then
        call register_axis(fileobj, "lon", nlon)
        call register_axis(fileobj, "lat", nlat)
        call register_axis(fileobj, "time", unlimited)

        call register_field(fileobj, "lon", "double", dimensions=(/"lon"/))
        call register_field(fileobj, "lat", "double", dimensions=(/"lat"/))
        call register_field(fileobj, "time", "double", dimensions=(/"time"/))

        call register_field(fileobj, fieldname, "double", dimensions=(/"lon ", "lat ", "time"/))

        call register_variable_attribute(fileobj, "lon", "cartesian_axis", "X", str_len=1)
        call register_variable_attribute(fileobj, "lat", "cartesian_axis", "Y", str_len=1)

        call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
        call register_variable_attribute(fileobj, "time", "units", "days since 1800-01-01 00:00:00",str_len=30)
        call register_variable_attribute(fileobj, "time", "calendar", "julian", str_len=6)

        !center points
        call write_data(fileobj, "lat", lat(1,1:nlat))
        call write_data(fileobj, "lon", lon(1:nlon,1))
        call write_data(fileobj, "time", times)

        do i=0, size(times) - 1
          call write_data(fileobj, fieldname, data, unlim_dim_level=i+1)
          data = - data
        enddo
        call close_file(fileobj)
      endif
    end if

  end subroutine generate_file

end program test_time_interp_with_interp
