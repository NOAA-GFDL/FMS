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

!> @brief  General program to test the different possible reduction methods
program test_reduction_methods
  use fms_mod,           only: fms_init, fms_end
  use testing_utils,     only: allocate_buffer, test_normal, test_openmp, test_halos, no_mask, logical_mask, real_mask
  use platform_mod,      only: r8_kind
  use block_control_mod, only: block_control_type, define_blocks
  use mpp_mod,           only: mpp_sync, FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe, mpp_broadcast, input_nml_file
  use time_manager_mod,  only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
  use diag_manager_mod,  only: diag_manager_init, diag_manager_end, diag_axis_init, register_diag_field, &
                               diag_send_complete, diag_manager_set_time_end, send_data
  use mpp_domains_mod,   only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain, &
                               mpp_get_data_domain

  implicit none

  integer                            :: nx              !< Number of points in the x direction
  integer                            :: ny              !< Number of points in the y direction
  integer                            :: nz              !< Number of points in the z direction
  integer                            :: nw              !< Number of points in the 4th dimension
  integer                            :: layout(2)       !< Layout
  integer                            :: io_layout(2)    !< Io layout
  type(domain2d)                     :: Domain          !< 2D domain
  integer                            :: isc, isd        !< Starting x compute, data domain index
  integer                            :: iec, ied        !< Ending x compute, data domain index
  integer                            :: jsc, jsd        !< Starting y compute, data domaine index
  integer                            :: jec, jed        !< Ending y compute, data domain index
  integer                            :: nhalox          !< Number of halos in x
  integer                            :: nhaloy          !< Number of halos in y
  real(kind=r8_kind), allocatable    :: cdata(:,:,:,:)  !< Data in the compute domain
  real(kind=r8_kind), allocatable    :: ddata(:,:,:,:)  !< Data in the data domain
  real(kind=r8_kind), allocatable    :: crmask(:,:,:,:) !< Mask in the compute domain
  real(kind=r8_kind), allocatable    :: drmask(:,:,:,:) !< Mask in the data domain
  logical,            allocatable    :: clmask(:,:,:,:) !< Logical mask in the compute domain
  logical,            allocatable    :: dlmask(:,:,:,:) !< Logical mask in the data domain
  type(time_type)                    :: Time            !< Time of the simulation
  type(time_type)                    :: Time_step       !< Time of the simulation
  integer                            :: ntimes          !< Number of times
  integer                            :: id_x            !< axis id for the x dimension
  integer                            :: id_y            !< axis id for the y dimension
  integer                            :: id_z            !< axis id for the z dimension
  integer                            :: id_w            !< axis id for the w dimension
  integer                            :: id_var0         !< diag_field id for 0d var
  integer                            :: id_var1         !< diag_field id for 1d var
  integer                            :: id_var2         !< diag_field id for 2d var
  integer                            :: id_var3         !< diag_field id for 3d var
  integer                            :: id_var4         !< diag_field id for 4d var
  integer                            :: io_status       !< Status after reading the namelist
  type(block_control_type)           :: my_block        !< Returns instantiated @ref block_control_type
  logical                            :: message         !< Flag for outputting debug message
  integer                            :: isd1            !< Starting x data domain index (1-based)
  integer                            :: ied1            !< Ending x data domain index (1-based)
  integer                            :: jsd1            !< Starting y data domain index (1-based)
  integer                            :: jed1            !< Ending y data domain index (1-based)
  integer                            :: isw             !< Starting index for each thread in the x direction
  integer                            :: iew             !< Ending index for each thread in the x direction
  integer                            :: jsw             !< Starting index for each thread in the y direction
  integer                            :: jew             !< Ending index for each thread in the y direction
  integer                            :: is1             !< Starting index for each thread in the x direction (1-based)
  integer                            :: ie1             !< Ending index for each thread in the x direction (1-based)
  integer                            :: js1             !< Starting index for each thread in the y direction (1-based)
  integer                            :: je1             !< Ending index for each thread in the y direction (1-based)
  integer                            :: iblock          !< For looping through the blocks
  integer                            :: i               !< For do loops
  logical                            :: used            !< Dummy argument to send_data
  real(kind=r8_kind)                 :: missing_value   !< Missing value to use

  !< Configuration parameters
  integer :: test_case = test_normal !< Indicates which test case to run
  integer :: mask_case = no_mask     !< Indicates which masking option to run

  namelist / test_reduction_methods_nml / test_case, mask_case

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  read (input_nml_file, test_reduction_methods_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_modern_diag: Error reading input.nml')

  Time = set_date(2,1,1,0,0,0)
  Time_step = set_time (3600,0) !< 1 hour
  nx = 96
  ny = 96
  nz = 5
  nw = 2
  layout = (/1, mpp_npes()/)
  io_layout = (/1, 1/)
  nhalox = 2
  nhaloy = 2
  ntimes = 48

  !< Create a lat/lon domain
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain, name='2D domain', xhalo=nhalox, yhalo=nhaloy)
  call mpp_define_io_domain(Domain, io_layout)
  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain, isd, ied, jsd, jed)

  cdata = allocate_buffer(isc, iec, jsc, jec, nz, nw)
  call init_buffer(cdata, isc, iec, jsc, jec, 0)

  select case (test_case)
  case (test_normal)
    if (mpp_pe() .eq. mpp_root_pe()) print *, "Testing the normal send_data calls"
  case (test_halos)
    if (mpp_pe() .eq. mpp_root_pe()) print *, "Testing the send_data calls with halos"
    ddata = allocate_buffer(isd, ied, jsd, jed, nz, nw)
    call init_buffer(ddata, isc, iec, jsc, jec, 2) !< The halos never get set
  case (test_openmp)
    if (mpp_pe() .eq. mpp_root_pe()) print *, "Testing the send_data calls with openmp blocks"
     call define_blocks ('testing_model', my_block, isc, iec, jsc, jec, kpts=0, &
                         nx_block=1, ny_block=4, message=message)
  end select

  select case (mask_case)
  case (logical_mask)
    clmask = allocate_logical_mask(isc, iec, jsc, jec, nz, nw)
    if (mpp_pe() .eq. 0) clmask(isc, jsc, 1, 1) = .False.

    if (test_case .eq. test_halos) then
      dlmask = allocate_logical_mask(isd, ied, jsd, jed, nz, nw)
      if (mpp_pe() .eq. 0) dlmask(1+nhalox, 1+nhaloy, 1, 1) = .False.
    endif
  case (real_mask)
    crmask = allocate_real_mask(isc, iec, jsc, jec, nz, nw)
    if (mpp_pe() .eq. 0) crmask(isc, jsc, 1, 1) = 0_r8_kind

    if (test_case .eq. test_halos) then
      drmask = allocate_real_mask(isd, ied, jsd, jed, nz, nw)
      if (mpp_pe() .eq. 0) drmask(1+nhalox, 1+nhaloy, 1, 1) = 0_r8_kind
    endif
  end select

  !< Register the axis
  id_x  = diag_axis_init('x',  real((/ (i, i = 1,nx) /), kind=r8_kind),  'point_E', 'x', long_name='point_E', &
    Domain2=Domain)
  id_y  = diag_axis_init('y',  real((/ (i, i = 1,ny) /), kind=r8_kind),  'point_N', 'y', long_name='point_N', &
    Domain2=Domain)
  id_z  = diag_axis_init('z',  real((/ (i, i = 1,nz) /), kind=r8_kind),  'point_Z', 'z', long_name='point_Z')
  id_w  = diag_axis_init('w',  real((/ (i, i = 1,nw) /), kind=r8_kind),  'point_W', 'n', long_name='point_W')

  missing_value = -666._r8_kind
  !< Register the fields
  id_var0 = register_diag_field  ('ocn_mod', 'var0', Time, 'Var0d', &
    'mullions', missing_value = missing_value)
  id_var1 = register_diag_field  ('ocn_mod', 'var1', (/id_x/), Time, 'Var1d', &
    'mullions', missing_value = missing_value)
  id_var2 = register_diag_field  ('ocn_mod', 'var2', (/id_x, id_y/), Time, 'Var2d', &
    'mullions', missing_value = missing_value)
  id_var3 = register_diag_field  ('ocn_mod', 'var3', (/id_x, id_y, id_z/), Time, 'Var3d', &
    'mullions', missing_value = missing_value)
  id_var4 = register_diag_field  ('ocn_mod', 'var4', (/id_x, id_y, id_z, id_w/), Time, 'Var4d', &
    'mullions', missing_value = missing_value)

  !< Get the data domain indices (1 based)
  isd1 = isc-isd+1
  jsd1 = jsc-jsd+1
  ied1 = isd1 + iec-isc
  jed1 = jsd1 + jec-jsc

  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))
  do i = 1, ntimes
    Time = Time + Time_step

    call set_buffer(cdata, i)
    used = send_data(id_var0, cdata(1,1,1,1), Time)

    select case(test_case)
    case (test_normal)
      select case (mask_case)
      case (no_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time)
        used = send_data(id_var2, cdata(:,:,1,1), Time)
        used = send_data(id_var3, cdata(:,:,:,1), Time)
      case (real_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time, rmask=crmask(:,1,1,1))
        used = send_data(id_var2, cdata(:,:,1,1), Time, rmask=crmask(:,:,1,1))
        used = send_data(id_var3, cdata(:,:,:,1), Time, rmask=crmask(:,:,:,1))
      case (logical_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time, mask=clmask(:,1,1,1))
        used = send_data(id_var2, cdata(:,:,1,1), Time, mask=clmask(:,:,1,1))
        used = send_data(id_var3, cdata(:,:,:,1), Time, mask=clmask(:,:,:,1))
      end select
    case (test_halos)
      call set_buffer(ddata, i)
      select case (mask_case)
      case (no_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time)
        used = send_data(id_var2, ddata(:,:,1,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1)
        used = send_data(id_var3, ddata(:,:,:,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1)
      case (real_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time, &
          rmask=crmask(:,1,1,1))
        used = send_data(id_var2, ddata(:,:,1,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1, &
          rmask=drmask(:,:,1,1))
        used = send_data(id_var3, ddata(:,:,:,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1, &
          rmask=drmask(:,:,:,1))
      case (logical_mask)
        used = send_data(id_var1, cdata(:,1,1,1), Time, &
          mask=clmask(:,1,1,1))
        used = send_data(id_var2, ddata(:,:,1,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1, &
          mask=dlmask(:,:,1,1))
        used = send_data(id_var3, ddata(:,:,:,1), Time, &
          is_in=isd1, ie_in=ied1, js_in=jsd1, je_in=jed1, &
          mask=dlmask(:,:,:,1))
      end select
    case (test_openmp)
      select case(mask_case)
      case (no_mask)
        used=send_data(id_var1, cdata(:, 1, 1, 1), time)
      case (logical_mask)
        used=send_data(id_var1, cdata(:, 1, 1, 1), time, &
            mask=clmask(:, 1, 1, 1))
      case (real_mask)
        used=send_data(id_var1, cdata(:, 1, 1, 1), time, &
            rmask=crmask(:, 1, 1, 1))
      end select
!$OMP parallel do default(shared) private(iblock, isw, iew, jsw, jew, is1, ie1, js1, je1)
      do iblock=1, 4
        isw = my_block%ibs(iblock)
        jsw = my_block%jbs(iblock)
        iew = my_block%ibe(iblock)
        jew = my_block%jbe(iblock)

      !--- indices for 1-based arrays ---
        is1 = isw-isc+1
        ie1 = iew-isc+1
        js1 = jsw-jsc+1
        je1 = jew-jsc+1

        select case (mask_case)
        case (no_mask)
          used=send_data(id_var2, cdata(is1:ie1, js1:je1, 1, 1), time, is_in=is1, js_in=js1)
          used=send_data(id_var3, cdata(is1:ie1, js1:je1, :, 1), time, is_in=is1, js_in=js1)
        case (real_mask)
          used=send_data(id_var2, cdata(is1:ie1, js1:je1, 1, 1), time, is_in=is1, js_in=js1, &
            rmask=crmask(is1:ie1, js1:je1, 1, 1))
          used=send_data(id_var3, cdata(is1:ie1, js1:je1, :, 1), time, is_in=is1, js_in=js1, &
            rmask=crmask(is1:ie1, js1:je1, :, 1))
        case (logical_mask)
          used=send_data(id_var2, cdata(is1:ie1, js1:je1, 1, 1), time, is_in=is1, js_in=js1, &
            mask=clmask(is1:ie1, js1:je1, 1, 1))
          used=send_data(id_var3, cdata(is1:ie1, js1:je1, :, 1), time, is_in=is1, js_in=js1, &
            mask=clmask(is1:ie1, js1:je1, :, 1))
        end select
      enddo
    end select
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)

  call fms_end

  contains

  !> @brief Allocate the logical mask based on the starting/ending indices
  !! @return logical mask initiliazed to .True.
  function allocate_logical_mask(is, ie, js, je, k, l) &
  result(buffer)
    integer, intent(in) :: is !< Starting x index
    integer, intent(in) :: ie !< Ending x index
    integer, intent(in) :: js !< Starting y index
    integer, intent(in) :: je !< Ending y index
    integer, intent(in) :: k  !< Number of points in the 4th dimension
    integer, intent(in) :: l  !< Number of points in the 5th dimension

    logical, allocatable :: buffer(:,:,:,:)

    allocate(buffer(is:ie, js:je, 1:k, 1:l))
    buffer = .True.
  end function allocate_logical_mask

  !> @brief Allocate the real mask based on the starting/ending indices
  !! @returnreal mask initiliazed to 1_r8_kind
  function allocate_real_mask(is, ie, js, je, k, l) &
  result(buffer)
    integer, intent(in) :: is !< Starting x index
    integer, intent(in) :: ie !< Ending x index
    integer, intent(in) :: js !< Starting y index
    integer, intent(in) :: je !< Ending y index
    integer, intent(in) :: k  !< Number of points in the 4th dimension
    integer, intent(in) :: l  !< Number of points in the 5th dimension
    real(kind=r8_kind), allocatable :: buffer(:,:,:,:)

    allocate(buffer(is:ie, js:je, 1:k, 1:l))
    buffer = 1.0_r8_kind
  end function allocate_real_mask

  !> @brief initiliazed the buffer based on the starting/ending indices
  subroutine init_buffer(buffer, is, ie, js, je, nhalo)
    real(kind=r8_kind), intent(inout) :: buffer(:,:,:,:) !< output buffer
    integer,            intent(in)    :: is              !< Starting x index
    integer,            intent(in)    :: ie              !< Ending x index
    integer,            intent(in)    :: js              !< Starting y index
    integer,            intent(in)    :: je              !< Ending y index
    integer,            intent(in)    :: nhalo           !< Number of halos

    integer :: ii, j, k, l

    do ii = is, ie
      do j = js, je
        do k = 1, size(buffer, 3)
          do l = 1, size(buffer,4)
            buffer(ii-is+1+nhalo, j-js+1+nhalo, k, l) = real(ii, kind=r8_kind)* 1000_r8_kind + &
              real(j, kind=r8_kind)* 10_r8_kind + &
              real(k, kind=r8_kind)
          enddo
        enddo
      enddo
    enddo

  end subroutine init_buffer

  !> @brief Set the buffer based on the time_index
  subroutine set_buffer(buffer, time_index)
    real(kind=r8_kind), intent(inout) :: buffer(:,:,:,:) !< Output buffer
    integer,            intent(in)    :: time_index      !< Time index

    buffer = nint(buffer) + real(time_index, kind=r8_kind)/100_r8_kind

  end subroutine set_buffer

end program test_reduction_methods
