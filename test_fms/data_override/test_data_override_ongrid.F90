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
                             input_nml_file, mpp_sync, NOTE, mpp_npes, mpp_get_current_pelist, &
                             mpp_set_current_pelist
use data_override_mod, only: data_override_init, data_override
use fms2_io_mod
use time_manager_mod,  only: set_calendar_type, time_type, set_date, NOLEAP
use netcdf,            only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
                             nf90_close, nf90_put_att, nf90_clobber, nf90_64bit_offset, nf90_char, &
                             nf90_double, nf90_unlimited
use ensemble_manager_mod, only: get_ensemble_size, ensemble_manager_init
use fms_mod, only: string, fms_init, fms_end

implicit none

integer, dimension(2)                      :: layout = (/2,3/) !< Domain layout
integer                                    :: nlon = 360       !< Number of points in x axis
integer                                    :: nlat = 180       !< Number of points in y axis
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
integer, parameter                         :: weight_file = 4
integer, parameter                         :: ensemble_case = 5
integer, parameter                         :: ensemble_same_yaml = 6
integer, parameter                         :: multi_file = 7
integer                                    :: test_case = ongrid
logical                                    :: init_with_mode = .false.
integer                                    :: npes
integer, allocatable                       :: pelist(:)
integer, allocatable                       :: pelist_ens(:)
integer                                    :: ensemble_id
logical                                    :: write_only=.false. !< True if creating the input files only

namelist /test_data_override_ongrid_nml/ nhalox, nhaloy, test_case, init_with_mode, nlon, nlat, layout, &
                                         write_only

call fms_init
call fms2_io_init

read (input_nml_file, test_data_override_ongrid_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>test_data_override_ongrid: Error reading input.nml')

!< Wait for the root PE to catch up
call mpp_sync

!< This is the actual test code:

call set_calendar_type(NOLEAP)

npes = mpp_npes()
allocate(pelist(npes))
call mpp_get_current_pelist(pelist)

select case (test_case)
case (ensemble_case, ensemble_same_yaml)
  call set_up_ensemble_case()
end select

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, xhalo=nhalox, yhalo=nhaloy, name='test_data_override_emc')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, is, ie, js, je)

select case (test_case)
case (ensemble_case, ensemble_same_yaml)
  ! Go back to the full pelist
  call mpp_set_current_pelist(pelist)
end select

if (write_only) then
  select case (test_case)
  case (ongrid)
    call generate_ongrid_input_file ()
  case (bilinear)
    call generate_bilinear_input_file ()
  case (scalar)
    call generate_scalar_input_file ()
  case (weight_file)
    call generate_weight_input_file ()
  case (ensemble_case, ensemble_same_yaml)
    call generate_ensemble_input_file()
  case (multi_file)
    call generate_multi_file_input_file()
  end select

  call mpp_sync()
  call mpp_error(NOTE, "Finished creating INPUT Files")

else
  select case (test_case)
  case (ensemble_case, ensemble_same_yaml)
    !< Go back to the ensemble pelist
    call mpp_set_current_pelist(pelist_ens)
  end select

  if (init_with_mode) then
    ! Initialize data_override_mod with explicit mode arguments
    call data_override_init(Ocean_domain_in=Domain, mode=r4_kind)
    call data_override_init(Ocean_domain_in=Domain, mode=r8_kind)
  else
    ! Initialize data_override_mod with no mode argument
    call data_override_init(Ocean_domain_in=Domain)
  endif

  select case (test_case)
  case (ongrid)
    call ongrid_test_r4
    call ongrid_test_r8
  case (bilinear)
    call bilinear_test_r4
    call bilinear_test_r8
  case (scalar)
    call scalar_test_r4
    call scalar_test_r8
  case (weight_file)
    call weight_file_test_r4
    call weight_file_test_r8
  case (ensemble_case, ensemble_same_yaml)
    call ensemble_test_r4
    call ensemble_test_r8
    call mpp_set_current_pelist(pelist)
  case (multi_file)
    call multi_file_r4
    call multi_file_r8
  end select
endif

call fms_end

contains

subroutine create_grid_spec_file
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

subroutine create_ocean_mosaic_file
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

subroutine create_ocean_hgrid_file
  type(FmsNetcdfFile_t) :: fileobj
  real(r4_kind), allocatable, dimension(:,:) :: xdata, ydata
  integer :: nx, nxp, ny, nyp, i, j

  nx = nlon*2
  nxp = nx+1
  ny = nlat*2
  nyp = ny+1

  allocate(xdata(nxp, nyp))
  xdata(1,:) = 0_r4_kind
  do i = 2, nxp
    xdata(i,:) = xdata(i-1,:) + 0.5_r4_kind
  enddo

  allocate(ydata(nxp, nyp))
  ydata(:,1) = -90.0_r4_kind
  do i = 2, nyp
    ydata(:,i) = ydata(:, i-1) + 0.5_r4_kind
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

subroutine create_ongrid_data_file(is_ensemble)
  logical, intent(in), optional :: is_ensemble
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(3)
  real(r4_kind), allocatable, dimension(:,:,:) :: runoff_in
  real(r4_kind), allocatable, dimension(:)     :: time_data
  integer :: offset
  character(len=256), allocatable :: appendix

  integer :: i

  offset = 0
  appendix = ""

  if (present(is_ensemble)) then
    offset = ensemble_id
    call get_filename_appendix(appendix)
    appendix = "_"//trim(appendix)
  endif

  allocate(runoff_in(nlon, nlat, 10))
  allocate(time_data(10))

  do i = 1, 10
    runoff_in(:,:,i) = real(i+offset, r4_kind)
  enddo

  time_data = (/1_r4_kind, 2_r4_kind, &
                3_r4_kind, 5_r4_kind, &
                6_r4_kind, 7_r4_kind, &
                8_r4_kind, 9_r4_kind, &
                10_r4_kind, 11_r4_kind/)

  dimnames(1) = 'i'
  dimnames(2) = 'j'
  dimnames(3) = 'time'

  if (open_file(fileobj, 'INPUT/runoff.daitren.clim.1440x1080.v20180328'//trim(appendix)//'.nc', 'overwrite')) then
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

subroutine create_multi_file_data_file
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(3)
  real(r4_kind), allocatable, dimension(:,:,:) :: runoff_in
  real(r4_kind), allocatable, dimension(:)     :: time_data
  integer :: i

  integer, parameter :: ntimes = 9
  integer, parameter :: nfiles = 3
  integer :: starting
  integer :: ending

  allocate(runoff_in(nlon, nlat, ntimes))
  allocate(time_data(ntimes))

  time_data = (/1_r4_kind, 2_r4_kind, 3_r4_kind, &
                4_r4_kind, 5_r4_kind, 6_r4_kind, &
                7_r4_kind, 7_r4_kind, 9_r4_kind /)

  do i = 1, ntimes
    runoff_in(:,:,i) = real(i, r4_kind)
  enddo

  dimnames(1) = 'i'
  dimnames(2) = 'j'
  dimnames(3) = 'time'

  do i = 1, nfiles
    if (open_file(fileobj, "INPUT/hadisst_ice.data_yr"//string(i-1)//".nc", "overwrite")) then
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

      starting = 2 * (i - 1) + 1
      ending = starting + 1

      call write_data(fileobj, "runoff", runoff_in(:,:, starting:ending))
      call write_data(fileobj, "time", time_data(starting:ending))
      call close_file(fileobj)
    endif
  enddo

end subroutine create_multi_file_data_file

subroutine generate_ongrid_input_file
  !< Create some files needed by data_override!
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_ongrid_data_file()
  endif
  call mpp_sync()
end subroutine generate_ongrid_input_file

subroutine generate_multi_file_input_file
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_multi_file_data_file()
  endif

  call mpp_sync()
end subroutine generate_multi_file_input_file

!> @brief Creates an input netcdf data file to use for the ongrid data_override test case
!! with either an increasing or decreasing lat, lon grid
subroutine create_bilinear_data_file(increasing_grid)
  logical, intent(in) :: increasing_grid !< .true. if increasing a file with an increasing lat/lon

  type(FmsNetcdfFile_t)         :: fileobj          !< Fms2_io fileobj
  character(len=10)             :: dimnames(3)      !< dimension names for the variable
  real(r4_kind), allocatable    :: runoff_in(:,:,:) !< Data to write
  real(r4_kind), allocatable    :: time_data(:)     !< Time dimension data
  real(r4_kind), allocatable    :: lat_data(:)      !< Lat dimension data
  real(r4_kind), allocatable    :: lon_data(:)      !< Lon dimension data
  character(len=:), allocatable :: filename         !< Name of the file
  integer                       :: factor           !< This is used when creating the grid data
                                                      !! -1 if the grid is decreasing
                                                      !! +1 if the grid is increasing
  integer                       :: i, j, k          !< For looping through variables
  integer                       :: nlon_data, nlat_data

  nlon_data = nlon + 1
  nlat_data = nlat - 1
  allocate(runoff_in(nlon_data, nlat_data, 10))
  allocate(time_data(10))
  allocate(lat_data(nlat_data))
  allocate(lon_data(nlon_data))

  if (.not. increasing_grid) then
    filename = 'INPUT/bilinear_decreasing.nc'
    lon_data(1) = 360.0_r4_kind
    lat_data(1) = 89.0_r4_kind
    factor = -1
    do i = 1, nlon_data
      do j = 1, nlat_data
        do k = 1, 10
          runoff_in(i, j, k) = real(362-i, kind=r4_kind) * 1000._r4_kind + &
            real(180-j, kind=r4_kind) + real(k, kind=r4_kind)/100._r4_kind
        enddo
      enddo
    enddo
  else
    filename = 'INPUT/bilinear_increasing.nc'
    lon_data(1) = 0.0_r4_kind
    lat_data(1) = -89.0_r4_kind
    factor = 1

    do i = 1, nlon_data
      do j = 1, nlat_data
        do k = 1, 10
          runoff_in(i, j, k) = real(i, kind=r4_kind) * 1000._r4_kind + real(j, kind=r4_kind) + &
            real(k, kind=r4_kind)/100._r4_kind
        enddo
      enddo
    enddo
  endif

  do i = 2, nlon_data
    lon_data(i) = real(lon_data(i-1) + 1*factor, r4_kind)
  enddo

  do i = 2, nlat_data
    lat_data(i) =real(lat_data(i-1) + 1*factor, r4_kind)
  enddo

  time_data = (/1_r4_kind, 2_r4_kind, &
                3_r4_kind, 5_r4_kind, &
                6_r4_kind, 7_r4_kind, &
                8_r4_kind, 9_r4_kind, &
                10_r4_kind, 11_r4_kind/)

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
subroutine generate_bilinear_input_file
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file ()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_bilinear_data_file(.true.)
    call create_bilinear_data_file(.false.)
  endif
  call mpp_sync()
end subroutine generate_bilinear_input_file

subroutine generate_weight_input_file
  call create_grid_spec_file ()
  call create_ocean_mosaic_file()
  call create_ocean_hgrid_file()
  call create_bilinear_data_file(.true.)
  call create_weight_file()
end subroutine generate_weight_input_file

subroutine create_weight_file
  type(FmsNetcdfFile_t) :: fileobj
  real(kind=r8_kind), allocatable :: vdata(:,:,:)
  character(len=5) :: dim_names(3)

  dim_names(1) = "nlon"
  dim_names(2) = "nlat"
  if (open_file(fileobj, "INPUT/remap_file.nc", "overwrite")) then
    call register_axis(fileobj, "nlon", nlon)
    call register_axis(fileobj, "nlat", nlat)
    call register_axis(fileobj, "three", 3)
    call register_axis(fileobj, "four", 4)

    dim_names(3) = "three"
    call register_field(fileobj, "index", "int", dim_names)

    dim_names(3) = "four"
    call register_field(fileobj, "weight", "double", dim_names)

    allocate(vdata(nlon,nlat,3))
    vdata(1,:,1) = 1
    vdata(2,:,1) = 2
    vdata(3,:,1) = 3
    vdata(4,:,1) = 4
    vdata(5,:,1) = 5
    vdata(:,1:2,2) = 1
    vdata(:,3,2) = 2
    vdata(:,4,2) = 3
    vdata(:,5,2) = 4
    vdata(:,6,2) = 5
    vdata(:,:,3) = 1
    call write_data(fileobj, "index", vdata)
    deallocate(vdata)

    allocate(vdata(nlon,nlat,4))
    vdata = 0.5_r8_kind
    vdata(:,1,3) = 1_r8_kind
    vdata(:,6,3) = 1_r8_kind
    vdata(:,1,4) = 0_r8_kind
    vdata(:,6,4) = 0_r8_kind

    call write_data(fileobj, "weight", vdata)
    deallocate(vdata)

    call close_file(fileobj)
  endif
end subroutine create_weight_file

!> @brief Generates the input for the bilinear data_override test_case
subroutine generate_scalar_input_file
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file ()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
    call create_scalar_data_file()
  endif
  call mpp_sync()
end subroutine generate_scalar_input_file

subroutine create_scalar_data_file
  type(FmsNetcdfFile_t) :: fileobj
  character(len=10) :: dimnames(1)
  real(r4_kind), allocatable, dimension(:)     :: co2_in
  real(r4_kind), allocatable, dimension(:)     :: time_data
  integer :: i

  allocate(co2_in(10))
  allocate(time_data(10))

  do i = 1, 10
    co2_in(i) = real(i, r4_kind)
  enddo

  time_data = (/1_r4_kind, 2_r4_kind, &
                3_r4_kind, 5_r4_kind, &
                6_r4_kind, 7_r4_kind, &
                8_r4_kind, 9_r4_kind, &
                10_r4_kind, 11_r4_kind/)

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

subroutine set_up_ensemble_case
  integer :: ens_siz(6)
  character(len=10) :: text

  if (npes .ne. 12) &
    call mpp_error(FATAL, "This test requires 12 pes to run")

  if (layout(1)*layout(2) .ne. 6) &
    call mpp_error(FATAL, "The two members of the layout do not equal 6")

  call ensemble_manager_init
  ens_siz = get_ensemble_size()
  if (ens_siz(1) .ne. 2) &
    call mpp_error(FATAL, "This test requires 2 ensembles")

  if (mpp_pe() < 6) then
    !PEs 0-5 are the first ensemble
    ensemble_id = 1
    allocate(pelist_ens(npes/ens_siz(1)))
    pelist_ens = pelist(1:6)
    call mpp_set_current_pelist(pelist_ens)
  else
    !PEs 6-11 are the second ensemble
    ensemble_id = 2
    allocate(pelist_ens(npes/ens_siz(1)))
    pelist_ens = pelist(7:)
    call mpp_set_current_pelist(pelist_ens)
  endif

  write( text,'(a,i2.2)' ) 'ens_', ensemble_id
  call set_filename_appendix(trim(text))

  if (mpp_pe() .eq. mpp_root_pe()) &
  print *, "ensemble_id:", ensemble_id, ":: ", pelist_ens
end subroutine set_up_ensemble_case

subroutine generate_ensemble_input_file
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_grid_spec_file ()
    call create_ocean_mosaic_file()
    call create_ocean_hgrid_file()
  endif

  !< Go back to the ensemble pelist so that each root pe can write its own input file
  call mpp_set_current_pelist(pelist_ens)
  if (mpp_pe() .eq. mpp_root_pe()) then
    call create_ongrid_data_file(is_ensemble=.true.)
  endif
  call mpp_set_current_pelist(pelist)
end subroutine generate_ensemble_input_file

#include "test_data_override_ongrid_r4.fh"
#include "test_data_override_ongrid_r8.fh"

end program test_data_override_ongrid
