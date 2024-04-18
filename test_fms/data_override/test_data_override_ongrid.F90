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

program test_data_override_ongrid

!> @brief  This programs tests data_override ability to override data for an
!! on grid case and when using bilinear interpolation

use platform_mod,      only: r4_kind, r8_kind
use mpp_domains_mod,   only: mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, &
                             mpp_domains_set_stack_size, mpp_get_compute_domain, domain2d
use mpp_mod,           only: mpp_init, mpp_exit, mpp_pe, mpp_root_pe, mpp_error, FATAL, &
                             input_nml_file, mpp_sync, NOTE
use data_override_mod, only: data_override_init, data_override
use fms2_io_mod
use time_manager_mod,  only: set_calendar_type, time_type, set_date, NOLEAP
use netcdf,            only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
                             nf90_close, nf90_put_att, nf90_clobber, nf90_64bit_offset, nf90_char, &
                             nf90_double, nf90_unlimited
use fms_mod, only: string

implicit none

integer, parameter                         :: lkind = DO_TEST_KIND_
integer, dimension(2)                      :: layout = (/2,3/) !< Domain layout
integer                                    :: nlon             !< Number of points in x axis
integer                                    :: nlat             !< Number of points in y axis
type(domain2d)                             :: Domain           !< Domain with mask table
integer                                    :: is               !< Starting x index
integer                                    :: ie               !< Ending x index
integer                                    :: js               !< Starting y index
integer                                    :: je               !< Ending y index
integer                                    :: nhalox=2, nhaloy=2
integer                                    :: io_status
integer, parameter                         :: ongrid = 1
integer, parameter                         :: bilinear = 2
integer, parameter                         :: scalar = 3
integer                                    :: test_case = ongrid

namelist / test_data_override_ongrid_nml / nhalox, nhaloy, test_case

call mpp_init
call fms2_io_init

read (input_nml_file, test_data_override_ongrid_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>test_data_override_ongrid: Error reading input.nml')



!< Wait for the root PE to catch up
call mpp_sync

!< This is the actual test code:

call set_calendar_type(NOLEAP)

nlon = 360
nlat = 180

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, xhalo=nhalox, yhalo=nhaloy, name='test_data_override_emc')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, is, ie, js, je)

select case (test_case)
case (ongrid)
  call generate_ongrid_input_file ()
case (bilinear)
  call generate_bilinear_input_file ()
case (scalar)
  call generate_scalar_input_file ()
end select

call mpp_sync()
call mpp_error(NOTE, "Finished creating INPUT Files")

!< Initiliaze data_override
call data_override_init(Ocean_domain_in=Domain, mode=lkind)

select case (test_case)
case (ongrid)
  call ongrid_test()
case (bilinear)
  call bilinear_test()
case (scalar)
  call scalar_test()
end select

call mpp_exit

contains

subroutine compare_data(Domain_in, actual_result, expected_result)
type(domain2d), intent(in)              :: Domain_in        !< Domain with mask table
real(lkind), intent(in)                 :: expected_result  !< Expected result from data_override
real(lkind), dimension(:,:), intent(in) :: actual_result    !< Result from data_override
integer                                 :: xsizec, ysizec   !< Size of the compute domain
integer                                 :: xsized, ysized   !< Size of the data domain
integer                                 :: nx, ny           !< Size of acual_result
integer                                 :: nhx, nhy   !< Size of the halos
integer                                 :: i, j             !< Helper indices

!< Data is only expected to be overriden for the compute domain -not at the halos.
call mpp_get_compute_domain(Domain_in, xsize=xsizec, ysize=ysizec)
call mpp_get_data_domain(Domain_in, xsize=xsized, ysize=ysized)

!< Note that actual_result has indices at (1:nx,1:ny) not (is:ie,js:je)
nhx= (xsized-xsizec)/2
nhy = (ysized-ysizec)/2
nx = size(actual_result, 1)
ny = size(actual_result, 2)

do i = 1, nx
   do j = 1, ny
      if (i <= nhx .or. i > (nx-nhx) .or. j <= nhy .or. j > (ny-nhy)) then
         !< This is the result at the halos it should 999.
         if (actual_result(i,j) .ne. 999._lkind) then
            print *, "for i=", i, " and j=", j, " result=", actual_result(i,j)
            call mpp_error(FATAL, "test_data_override_ongrid: Data was overriden in the halos!!")
         endif
      else
         if (actual_result(i,j) .ne. expected_result) then
            print *, "for i=", i, " and j=", j, " result=", actual_result(i,j), " expected=", expected_result
            call mpp_error(FATAL, "test_data_override_ongrid: Result is different from expected answer!")
         endif
      endif
   enddo
enddo

end subroutine

subroutine create_grid_spec_file()
  type(FmsNetcdfFile_t) :: fileobj

  if (open_file(fileobj, 'INPUT/grid_spec.nc', 'overwrite')) then
    call register_axis(fileobj, 'str', 255)
    call register_field(fileobj, 'ocn_mosaic_file', 'char', (/'str'/))
    call write_data(fileobj, 'ocn_mosaic_file', "ocean_mosaic.nc")
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/grid_spec.nc' to write")
  endif
end subroutine create_grid_spec_file

subroutine create_ocean_mosaic_file()
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(2)

  dimnames(1) = 'str'
  dimnames(2) = 'ntiles'
  if (open_file(fileobj, 'INPUT/ocean_mosaic.nc', 'overwrite')) then
    call register_axis(fileobj, dimnames(1) , 255)
    call register_axis(fileobj, dimnames(2), 1)
    call register_field(fileobj, 'gridfiles', 'char', dimnames)
    call write_data(fileobj, 'gridfiles', (/"ocean_hgrid.nc"/))
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/ocean_mosaic.nc' to write")
  endif
end subroutine create_ocean_mosaic_file

subroutine create_ocean_hgrid_file()
  type(FmsNetcdfFile_t) :: fileobj

  real(lkind), allocatable, dimension(:,:) :: xdata, ydata
  integer :: nx, nxp, ny, nyp, i, j
  nx = nlon*2
  nxp = nx+1
  ny = nlat*2
  nyp = ny+1

  allocate(xdata(nxp, nyp))
  xdata(1,:) = 0.0_lkind
  do i = 2, nxp
    xdata(i,:) = xdata(i-1,:) + 0.5
  enddo

  allocate(ydata(nxp, nyp))
  ydata(:,1) = -90.0_lkind
  do i = 2, nyp
    ydata(:,i) = ydata(:, i-1) + 0.5
  enddo

  if (open_file(fileobj, 'INPUT/ocean_hgrid.nc', 'overwrite')) then
    call register_axis(fileobj, "nx", nx)
    call register_axis(fileobj, "ny", ny)
    call register_axis(fileobj, "nxp", nxp)
    call register_axis(fileobj, "nyp", nyp)
    call register_field(fileobj, 'x', 'float', (/'nxp', 'nyp'/))
    call register_field(fileobj, 'y', 'float', (/'nxp', 'nyp'/))
    call register_field(fileobj, 'area', 'float', (/'nx', 'ny'/))
    call write_data(fileobj, "x", xdata)
    call write_data(fileobj, "y", ydata)
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/ocean_hgrid.nc' to write")
  endif
end subroutine create_ocean_hgrid_file

subroutine create_ongrid_data_file()
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(3)
  real(lkind), allocatable, dimension(:,:,:) :: runoff_in
  real(lkind), allocatable, dimension(:)     :: time_data
  integer :: i

  allocate(runoff_in(nlon, nlat, 10))
  allocate(time_data(10))
  do i = 1, 10
    runoff_in(:,:,i) = real(i, lkind)
  enddo
  time_data = (/1., 2., 3., 5., 6., 7., 8., 9., 10., 11./)

  dimnames(1) = 'i'
  dimnames(2) = 'j'
  dimnames(3) = 'time'

  if (open_file(fileobj, 'INPUT/runoff.daitren.clim.1440x1080.v20180328.nc', 'overwrite')) then
    call register_axis(fileobj, "i", nlon)
    call register_axis(fileobj, "j", nlat)
    call register_axis(fileobj, "time", unlimited)

    call register_field(fileobj, "i", "float", (/"i"/))
    call register_variable_attribute(fileobj, "i", "cartesian_axis", "x", str_len=1)

    call register_field(fileobj, "j", "float", (/"j"/))
    call register_variable_attribute(fileobj, "j", "cartesian_axis", "y", str_len=1)

    call register_field(fileobj, "time", "float", (/"time"/))
    call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
    call register_variable_attribute(fileobj, "time", "calendar", "noleap", str_len=6)
    call register_variable_attribute(fileobj, "time", "units", "days since 0001-01-01 00:00:00", str_len=30)

    call register_field(fileobj, "runoff", "float", dimnames)
    call write_data(fileobj, "runoff", runoff_in)
    call write_data(fileobj, "time", time_data)
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/runoff.daitren.clim.1440x1080.v20180328.nc' to write")
  endif
  deallocate(runoff_in)
end subroutine create_ongrid_data_file

subroutine generate_ongrid_input_file()
  !< Create some files needed by data_override!
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_ongrid_data_file()
  endif
  call mpp_sync()
end subroutine

!> @brief Tests ongrid data overrides.
!! In the first case there is no time interpolation
!! In the second case there is time interpolation
subroutine ongrid_test()
  real(lkind)                                :: expected_result  !< Expected result from data_override
  type(time_type)                            :: Time             !< Time
  real(lkind), allocatable, dimension(:,:)   :: runoff           !< Data to be written

  allocate(runoff(is:ie,js:je))

  runoff = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3.
  expected_result = 3._lkind
  call compare_data(Domain, runoff, expected_result)

  !< Run it when time=4
  runoff = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  call compare_data(Domain, runoff, expected_result)

  deallocate(runoff)
end subroutine ongrid_test

!> @brief Creates an input netcdf data file to use for the ongrid data_override test case
!! with either an increasing or decreasing lat, lon grid
subroutine create_bilinear_data_file(increasing_grid)
  logical, intent(in) :: increasing_grid !< .true. if increasing a file with an increasing lat/lon

  type(FmsNetcdfFile_t)           :: fileobj          !< Fms2_io fileobj
  character(len=10)               :: dimnames(3)      !< dimension names for the variable
  real(lkind),      allocatable   :: runoff_in(:,:,:) !< Data to write
  real(lkind),      allocatable   :: time_data(:)     !< Time dimension data
  real(lkind),      allocatable   :: lat_data(:)      !< Lat dimension data
  real(lkind),      allocatable   :: lon_data(:)      !< Lon dimension data
  character(len=:), allocatable   :: filename         !< Name of the file
  integer                         :: factor           !< This is used when creating the grid data
                                                      !! -1 if the grid is decreasing
                                                      !! +1 if the grid is increasing
  integer                         :: i, j, k          !< For looping through variables

  integer :: nlon_data
  integer :: nlat_data

  nlon_data = nlon + 1
  nlat_data = nlat - 1
  allocate(runoff_in(nlon_data, nlat_data, 10))
  allocate(time_data(10))
  allocate(lat_data(nlat_data))
  allocate(lon_data(nlon_data))

  if (.not. increasing_grid) then
    filename = 'INPUT/bilinear_decreasing.nc'
    lon_data(1) = 360.0_lkind
    lat_data(1) = 89.0_lkind
    factor = -1
    do i = 1, nlon_data
      do j = 1, nlat_data
        do k = 1, 10
          runoff_in(i, j, k) = real(362-i, kind=lkind) * 1000._lkind + &
            real(180-j, kind=lkind) + real(k, kind=lkind)/100._lkind
        enddo
      enddo
    enddo
  else
    filename = 'INPUT/bilinear_increasing.nc'
    lon_data(1) = 0.0_lkind
    lat_data(1) = -89.0_lkind
    factor = 1

    do i = 1, nlon_data
      do j = 1, nlat_data
        do k = 1, 10
          runoff_in(i, j, k) = real(i, kind=lkind) * 1000._lkind + real(j, kind=lkind) + &
            real(k, kind=lkind)/100._lkind
        enddo
      enddo
    enddo
  endif

  do i = 2, nlon_data
    lon_data(i) = real(lon_data(i-1) + 1*factor, lkind)
  enddo

  do i = 2, nlat_data
    lat_data(i) =real(lat_data(i-1) + 1*factor, lkind)
  enddo

  time_data = (/1., 2., 3., 5., 6., 7., 8., 9., 10., 11./)

  dimnames(1) = 'i'
  dimnames(2) = 'j'
  dimnames(3) = 'time'

  if (open_file(fileobj, filename, 'overwrite')) then
    call register_axis(fileobj, "i", nlon_data)
    call register_axis(fileobj, "j", nlat_data)
    call register_axis(fileobj, "time", unlimited)

    call register_field(fileobj, "i", "float", (/"i"/))
    call register_variable_attribute(fileobj, "i", "cartesian_axis", "x", str_len=1)

    call register_field(fileobj, "j", "float", (/"j"/))
    call register_variable_attribute(fileobj, "j", "cartesian_axis", "y", str_len=1)

    call register_field(fileobj, "time", "float", (/"time"/))
    call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
    call register_variable_attribute(fileobj, "time", "calendar", "noleap", str_len=6)
    call register_variable_attribute(fileobj, "time", "units", "days since 0001-01-01 00:00:00", str_len=30)

    call register_field(fileobj, "runoff", "float", dimnames)
    call write_data(fileobj, "runoff", runoff_in)
    call write_data(fileobj, "i", lon_data)
    call write_data(fileobj, "j", lat_data)
    call write_data(fileobj, "time", time_data)
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/bilinear_increasing.nc' to write")
  endif
  deallocate(runoff_in)
end subroutine create_bilinear_data_file

!> @brief Generates the input for the bilinear data_override test_case
subroutine generate_bilinear_input_file()
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file ()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_bilinear_data_file(.true.)
    call create_bilinear_data_file(.false.)
  endif
  call mpp_sync()
end subroutine generate_bilinear_input_file

!> @brief Tests bilinear data_override with and increasing and decreasing grid case
!! and comares the output betweeen the cases to ensure it is correct
subroutine bilinear_test()
  type(time_type)                            :: Time              !< Time
  real(lkind), allocatable, dimension(:,:)   :: runoff_decreasing !< Data to be written
  real(lkind), allocatable, dimension(:,:)   :: runoff_increasing !< Data to be written

  integer :: i, j, k
  logical :: success

  allocate(runoff_decreasing(is:ie,js:je))
  allocate(runoff_increasing(is:ie,js:je))

  runoff_decreasing = 999_lkind
  runoff_increasing = 999_lkind
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff_increasing',runoff_increasing, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")
  call data_override('OCN','runoff_decreasing',runoff_decreasing, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")

  do i = is, ie
    do j =  js, je
      if (abs(runoff_decreasing(i,j) - runoff_increasing(i,j)) .gt. 1) then
        call mpp_error(FATAL, "The data is not the same: "// &
        string(i)//","//string(j)//":"// &
        string(runoff_decreasing(i,j))//" vs "//string(runoff_increasing(i,j)))
      endif
    enddo
  enddo
  deallocate(runoff_decreasing, runoff_increasing)
end subroutine bilinear_test

!> @brief Generates the input for the bilinear data_override test_case
subroutine generate_scalar_input_file()
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file ()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_scalar_data_file()
  endif
  call mpp_sync()
end subroutine generate_scalar_input_file

subroutine create_scalar_data_file()
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(1)
  real(lkind), allocatable, dimension(:)     :: co2_in
  real(lkind), allocatable, dimension(:)     :: time_data
  integer :: i

  allocate(co2_in(10))
  allocate(time_data(10))
  do i = 1, 10
    co2_in(i) = real(i, lkind)
  enddo
  time_data = (/1., 2., 3., 5., 6., 7., 8., 9., 10., 11./)

  dimnames(1) = 'time'

  if (open_file(fileobj, 'INPUT/scalar.nc', 'overwrite')) then
    call register_axis(fileobj, "time", unlimited)
    call register_field(fileobj, "time", "float", (/"time"/))
    call register_variable_attribute(fileobj, "time", "cartesian_axis", "T", str_len=1)
    call register_variable_attribute(fileobj, "time", "calendar", "noleap", str_len=6)
    call register_variable_attribute(fileobj, "time", "units", "days since 0001-01-01 00:00:00", str_len=30)

    call register_field(fileobj, "co2", "float", dimnames)
    call write_data(fileobj, "co2", co2_in)
    call write_data(fileobj, "time", time_data)
    call close_file(fileobj)
  else
    call mpp_error(FATAL, "Error opening the file: 'INPUT/scalar.nc' to write")
  endif
  deallocate(co2_in)
end subroutine create_scalar_data_file

subroutine scalar_test()
  real(lkind)                                :: expected_result  !< Expected result from data_override
  type(time_type)                            :: Time             !< Time
  real(lkind)                                :: co2              !< Data to be written

  co2 = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','co2',co2, Time)
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3.
  expected_result = 3._lkind
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overriden to the correct value!")

  !< Run it when time=4
  co2 = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','co2',co2, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overriden to the correct value!")

end subroutine scalar_test

end program test_data_override_ongrid
