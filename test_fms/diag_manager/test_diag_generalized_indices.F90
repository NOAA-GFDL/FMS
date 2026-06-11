!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief  Test generalized axis permutations for send_data.
!!         Applies predefined permutations to canonical (x,y,z,w) storage
!!         and verifies consistency between data layout and axis metadata.
!!         Assumes default configuration parameters: test_normal + no_mask.
program test_diag_generalized_indices
  use fms_mod,           only: fms_init, fms_end
  use testing_utils,     only: allocate_buffer, permute
  use platform_mod,      only: r8_kind
  use mpp_mod,           only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe
  use time_manager_mod,  only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
  use diag_manager_mod,  only: diag_manager_init, diag_manager_end, diag_axis_init, register_diag_field, &
                               diag_send_complete, diag_manager_set_time_end, send_data
  use mpp_domains_mod,   only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain

  implicit none

  integer                :: nx, ny, nz, nw
  integer                :: layout(2), io_layout(2)
  type(domain2d)         :: Domain
  integer                :: isc, iec, jsc, jec
  integer                :: nhalox, nhaloy
  integer                :: ntimes, i
  type(time_type)        :: Time, Time_step
  real(r8_kind)          :: missing_value
  integer                :: ierr

  ! Axes
  integer                :: id_x, id_y, id_z, id_w
  integer                :: axis(4)

  ! Data
  real(r8_kind), allocatable :: cdata(:,:,:,:)   ! canonical storage: (x,y,z,w)

  ! Permutation test
  integer, parameter :: LAYOUT_XY  = 1
  integer, parameter :: LAYOUT_YX  = 2
  integer, parameter :: LAYOUT_ZX  = 3
  integer, parameter :: LAYOUT_YZX = 4
  integer, parameter :: LAYOUT_ZXY = 5
  integer, parameter :: PERM_TABLE(3,5) = reshape([ &
                           1,2,3, & ! XY (identity)
                           2,1,3, & ! YX
                           3,2,1, & ! ZX
                           2,3,1, & ! YZX
                           3,1,2  & ! ZXY
                        ], [3,5])

  integer :: id_var2_id, id_var2_yx
  integer :: id_var3_id, id_var3_zx, id_var3_yzx, id_var3_zxy

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  Time      = set_date(2,1,1,0,0,0)
  Time_step = set_time(3600,0) ! 1 hour

  nx = 96
  ny = 96
  nz = 5
  nw = 2
  ntimes = 48

  nhalox = 2
  nhaloy = 2
  layout    = (/1, mpp_npes()/)
  io_layout = (/1, 1/)

  ! Domain
  call mpp_define_domains((/1,nx,1,ny/), layout, Domain, name='2D domain', symmetry=.true., &
                          xhalo=nhalox, yhalo=nhaloy)
  call mpp_define_io_domain(Domain, io_layout)
  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)

  ! Data allocation + init (canonical x,y storage)
  cdata = allocate_buffer(isc, iec, jsc, jec, nz, nw)
  call init_buffer(cdata, isc, iec, jsc, jec, 0)

  ! Axes
  id_x = diag_axis_init('x', real((/(i, i=1,nx)/), kind=r8_kind), 'point_E', 'x', long_name='point_E', Domain2=Domain)
  id_y = diag_axis_init('y', real((/(i, i=1,ny)/), kind=r8_kind), 'point_N', 'y', long_name='point_N', Domain2=Domain)
  id_z = diag_axis_init('z', real((/(i, i=1,nz)/), kind=r8_kind), 'point_Z', 'z', long_name='point_Z')
  id_w = diag_axis_init('w', real((/(i, i=1,nw)/), kind=r8_kind), 'point_W', 'n', long_name='point_W')
  axis = [id_x, id_y, id_z, id_w]

  missing_value = -666._r8_kind

  ! Register permuted diagnostic fields ONCE
  id_var2_id  = register_perm_diag_field('var2_id',  'Var2d id',  axis(1:2), LAYOUT_XY)
  id_var3_id  = register_perm_diag_field('var3_id',  'Var3d id',  axis(1:3), LAYOUT_XY)

  id_var2_yx  = register_perm_diag_field('var2_yx',  'Var2d yx',  axis(1:2), LAYOUT_YX)
  id_var3_zx  = register_perm_diag_field('var3_zx',  'Var3d zx',  axis(1:3), LAYOUT_ZX)

  id_var3_yzx = register_perm_diag_field('var3_yzx', 'Var3d yzx', axis(1:3), LAYOUT_YZX)
  id_var3_zxy = register_perm_diag_field('var3_zxy', 'Var3d zxy', axis(1:3), LAYOUT_ZXY)

  if (mpp_pe() == mpp_root_pe()) then
    print *, "Testing generalized indices in default mode (test_normal + no_mask)"
    print *, "  canonical storage is (x,y,z,w)"
    print *, "  sending:"
    print *, "    var2_id with axes (x,y)"
    print *, "    var2_yx with axes (y,x)"
    print *, "    var3_id with axes (x,y,z)"
    print *, "    var3_zx with axes (z,y,x)"
    print *, "    var3_yzx with axes (y,z,x)"
    print *, "    var3_zxy with axes (z,x,y)"
  end if

  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  do i = 1, ntimes
    Time = Time + Time_step
    call set_buffer(cdata, i)

    ! Identity: axes (x,y) / (x,y,z) with canonical storage
    call send_perm_data(id_var2_id, cdata, PERM_TABLE(1:2, LAYOUT_XY), Time)
    call send_perm_data(id_var3_id, cdata, PERM_TABLE(1:3, LAYOUT_XY), Time)

    ! Swap: axes (y,x) / (z,y,x) while canonical storage remains (x,y,...) -> pack to temp and send
    call send_perm_data(id_var2_yx, cdata, PERM_TABLE(1:2, LAYOUT_YX), Time)
    call send_perm_data(id_var3_zx, cdata, PERM_TABLE(1:3, LAYOUT_ZX), Time)

    ! Cyclic: axes (y,z,x) / (z,x,y) while canonical storage remains (x, y, ...) -> pack to temp and send
    call send_perm_data(id_var3_yzx, cdata, PERM_TABLE(1:3, LAYOUT_YZX), Time)
    call send_perm_data(id_var3_zxy, cdata, PERM_TABLE(1:3, LAYOUT_ZXY), Time)

    call diag_send_complete(Time_step)
    call diag_send_complete(Time_step)
  end do

  call diag_manager_end(Time)
  call fms_end

contains

  !> @brief Apply a predefined permutation to an axis array.
  !> Maps canonical axis ordering to a permuted layout using PERM_TABLE.
  !> Supports rank-2 and rank-3 axis subsets.
  subroutine permute_axis(axis_in, perm_id, axis_out)
    integer, intent(in)  :: axis_in(:)
    integer, intent(in)  :: perm_id
    integer, intent(out) :: axis_out(:)

    integer :: order(3)

    order = PERM_TABLE(:, perm_id)

    if (any(order(1:size(axis_out)) > size(axis_in))) then
       call mpp_error(FATAL, "permute_axis: invalid permutation for given rank")
    endif

    axis_out = axis_in(order(1:size(axis_out)))
  end subroutine permute_axis

  !> @brief Register a diagnostic field with permuted axes.
  !> Applies axis permutation before calling register_diag_field.
  function register_perm_diag_field(var_name, long_name, axis, perm_id) result(id_var)
    character(len=*), intent(in)  :: var_name, long_name
    integer,          intent(in)  :: axis(:)
    integer,          intent(in)  :: perm_id

    integer              :: id_var
    integer, allocatable :: axis_perm(:)

    allocate(axis_perm(size(axis)))
    call permute_axis(axis, perm_id, axis_perm)
    id_var = register_diag_field('ocn_mod', var_name, axis_perm, Time, long_name, &
                                 'mullions', missing_value=missing_value)
    deallocate(axis_perm)
  end function register_perm_diag_field

  !> @brief Send data with optional axis permutation.
  !> Applies 2D or 3D permutation to canonical (x,y,z,w) buffers before send_data.
  !> Skips permutation for identity mappings.
  subroutine send_perm_data(id_field, buf, order, Time_in)
    integer, intent(in)         :: id_field
    real(r8_kind), intent(in)   :: buf(:,:,:,:)   ! canonical (x,y,z,w)
    integer, intent(in)         :: order(:)
    type(time_type), intent(in) :: Time_in

    logical :: used_local
    real(r8_kind), allocatable :: tmp2(:,:), tmp3(:,:,:)

    if (size(order) == 2) then
       if (all(order == [1,2])) then
          used_local = send_data(id_field, buf(:,:,1,1), Time_in)
       else
          tmp2 = permute(buf(:,:,1,1), order)
          used_local = send_data(id_field, tmp2, Time_in)
       endif

    else if (size(order) == 3) then
       if (all(order == [1,2,3])) then
          used_local = send_data(id_field, buf(:,:,:,1), Time_in)
       else
          tmp3 = permute(buf(:,:,:,1), order)
          used_local = send_data(id_field, tmp3, Time_in)
       endif

    else
       call mpp_error(FATAL, "send_var_perm: unsupported permutation rank")

    endif
  end subroutine send_perm_data

  !> @brief initialized the buffer based on the starting/ending indices
  subroutine init_buffer(buffer, is, ie, js, je, nhalo)
    real(r8_kind), intent(inout) :: buffer(:,:,:,:)
    integer, intent(in)          :: is, ie, js, je
    integer, intent(in)          :: nhalo

    integer :: ii, j, k, l

    do ii = is, ie
      do j = js, je
        do k = 1, size(buffer, 3)
          do l = 1, size(buffer, 4)
            buffer(ii-is+1+nhalo, j-js+1+nhalo, k, l) = real(ii, kind=r8_kind)*1000._r8_kind + &
                                                        real(j,  kind=r8_kind)*  10._r8_kind + &
                                                        real(k,  kind=r8_kind)
          end do
        end do
      end do
    end do
  end subroutine init_buffer

  !> @brief Set the buffer based on the time_index
  subroutine set_buffer(buffer, time_index)
    real(r8_kind), intent(inout) :: buffer(:,:,:,:)
    integer, intent(in)          :: time_index

    buffer = nint(buffer) + real(time_index, kind=r8_kind)/100._r8_kind
  end subroutine set_buffer

end program test_diag_generalized_indices

