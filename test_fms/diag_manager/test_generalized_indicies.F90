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

!> @brief  Test generalized axis permutations (x/y for now) for send_data
!!         Assumes default configuration parameters: test_normal + no_mask.
program test_generalized_indices
  use fms_mod,           only: fms_init, fms_end
  use testing_utils,     only: allocate_buffer
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

  ! Axes
  integer                :: id_x, id_y, id_z, id_w
  integer                :: axis(4)

  ! Data
  real(r8_kind), allocatable :: cdata(:,:,:,:)   ! canonical storage: (x,y,z,w)

  ! Permutation test
  integer :: p_id(3), p_swap(3)
  integer :: id_var2_id, id_var2_swap
  integer :: id_var3_id, id_var3_swap

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

  ! Define permutations: identity (x,y,z) and swap (y,x,z)
  p_id   = [1,2,3]
  p_swap = [2,1,3]

  ! Register permuted diagnostic fields ONCE
  id_var2_id   = register_diag_field('ocn_mod', 'var2_id',   (/axis(p_id(1)),   axis(p_id(2))/),   Time, 'Var2d id', &
                                     'mullions', missing_value=missing_value)
  id_var2_swap = register_diag_field('ocn_mod', 'var2_swap', (/axis(p_swap(1)), axis(p_swap(2))/), Time, 'Var2d swap', &
                                     'mullions', missing_value=missing_value)

  id_var3_id   = register_diag_field('ocn_mod', 'var3_id',   (/axis(p_id(1)),   axis(p_id(2)),   axis(p_id(3))/),   Time, &
                                     'Var3d id', 'mullions', missing_value=missing_value)
  id_var3_swap = register_diag_field('ocn_mod', 'var3_swap', (/axis(p_swap(1)), axis(p_swap(2)), axis(p_swap(3))/), Time, &
                                     'Var3d swap', 'mullions', missing_value=missing_value)

  if (mpp_pe() == mpp_root_pe()) then
    print *, "Testing generalized indices in default mode (test_normal + no_mask)"
    print *, "  canonical storage is (x,y,z,w)"
    print *, "  sending:"
    print *, "    var2_id   with axes (x,y)"
    print *, "    var2_swap with axes (y,x)"
    print *, "    var3_id   with axes (x,y,z)"
    print *, "    var3_swap with axes (y,x,z)"
  end if

  call diag_manager_set_time_end(set_date(2,1,3,0,0,0))

  do i = 1, ntimes
    Time = Time + Time_step
    call set_buffer(cdata, i)

    ! Identity: axes (x,y) / (x,y,z) with canonical storage
    call send_var2_perm(id_var2_id,   cdata, p_id,   Time)
    call send_var3_perm(id_var3_id,   cdata, p_id,   Time)

    ! Swap: axes (y,x) / (y,x,z) while canonical storage remains (x,y,...) -> pack to temp and send
    call send_var2_perm(id_var2_swap, cdata, p_swap, Time)
    call send_var3_perm(id_var3_swap, cdata, p_swap, Time)

    call diag_send_complete(Time_step)
    call diag_send_complete(Time_step)
  end do

  call diag_manager_end(Time)
  call fms_end

contains

  subroutine send_var2_perm(id_field, buf, p, Time_in)
    integer, intent(in)         :: id_field
    real(r8_kind), intent(in)   :: buf(:,:,:,:)   ! canonical (x,y,z,w)
    integer, intent(in)         :: p(3)
    type(time_type), intent(in) :: Time_in

    logical :: used_local
    real(r8_kind), allocatable :: tmp2(:,:)

    ! Support only identity (1,2,*) and xy-swap (2,1,*) for 2D
    if (p(1)==1 .and. p(2)==2) then
      used_local = send_data(id_field, buf(:,:,1,1), Time_in)
    else if (p(1)==2 .and. p(2)==1) then
      allocate(tmp2(size(buf,2), size(buf,1)))
      tmp2 = transpose(buf(:,:,1,1))
      used_local = send_data(id_field, tmp2, Time_in)
      deallocate(tmp2)
    else
      call mpp_error(FATAL, 'send_var2_perm: only p=(1,2,*) or (2,1,*) implemented')
    end if
  end subroutine send_var2_perm


  subroutine send_var3_perm(id_field, buf, p, Time_in)
    integer, intent(in)         :: id_field
    real(r8_kind), intent(in)   :: buf(:,:,:,:)   ! canonical (x,y,z,w)
    integer, intent(in)         :: p(3)
    type(time_type), intent(in) :: Time_in

    logical :: used_local
    integer :: nxloc, nyloc, nzloc, k
    real(r8_kind), allocatable :: tmp3(:,:,:)

    ! For now, support only keeping z as z
    if (p(3) /= 3) call mpp_error(FATAL, 'send_var3_perm: only permutations with p(3)=3 implemented')

    if (p(1)==1 .and. p(2)==2) then
      used_local = send_data(id_field, buf(:,:,:,1), Time_in)

    else if (p(1)==2 .and. p(2)==1) then
      nxloc = size(buf,1)
      nyloc = size(buf,2)
      nzloc = size(buf,3)

      allocate(tmp3(nyloc, nxloc, nzloc))
      do k = 1, nzloc
        tmp3(:,:,k) = transpose(buf(:,:,k,1))
      end do

      used_local = send_data(id_field, tmp3, Time_in)
      deallocate(tmp3)

    else
      call mpp_error(FATAL, 'send_var3_perm: only p=(1,2,3) or (2,1,3) implemented')
    end if
  end subroutine send_var3_perm


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

end program test_generalized_indices

