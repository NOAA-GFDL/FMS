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
program test_mpp_global_field
  use platform_mod
  use mpp_mod,         only : mpp_init, mpp_error, FATAL, mpp_init_test_requests_allocated
  use mpp_mod,         only : mpp_declare_pelist, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_domains_mod, only : domain2D
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER, XUPDATE, YUPDATE
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field
  use random_numbers_mod, only: randomNumberStream, initializeRandomNumberStream, getRandomNumbers
  use permutable_indices_mod, only: permutable_indices, factorial

  implicit none

  interface arr_compare
    procedure :: arr_compare_2d, arr_compare_3d
  end interface arr_compare

  interface arr_init
    procedure :: arr_init_2d, arr_init_3d
  end interface arr_init

  type test_params_t
    logical :: symmetry
    integer :: position, shift(2)
    character(15) :: name
  end type test_params_t

  type(test_params_t), parameter :: test_params(*) = [ &
      test_params_t(symmetry=.false., position=CENTER, shift=[0,0], name="No symmetry"), &
      test_params_t(symmetry=.true.,  position=CENTER, shift=[0,0], name="Center symmetry"), &
      test_params_t(symmetry=.true.,  position=CORNER, shift=[1,1], name="Corner symmetry"), &
      test_params_t(symmetry=.true.,  position=EAST,   shift=[1,0], name="East symmetry"), &
      test_params_t(symmetry=.true.,  position=NORTH,  shift=[0,1], name="North symmetry")]

  integer, parameter :: nx=20, ny=20, nz=40
  integer, parameter :: whalo=2, ehalo=2, shalo=2, nhalo=2
  integer, parameter :: stackmax=4000000

  FMS_TEST_TYPE_ (FMS_TEST_KIND_), parameter :: zero = 0

  integer :: pe, npes, ierr
  integer :: layout(2)
  integer :: i, p

  !> call mpp_init
  call mpp_init(test_level=mpp_init_test_requests_allocated)

  !> get pe info
  pe = mpp_pe()
  npes = mpp_npes()

  !> initialize mpp domain(s)
  call mpp_domains_init()
  call mpp_domains_set_stack_size(stackmax)

  do i=1, size(test_params)
    ! 2D tests
    do p=1,factorial(2)
      call run_tests_2d(test_params(i), p)
    enddo

    ! 3D tests
    do p=1,factorial(3)
      call run_tests_3d(test_params(i), p)
    enddo
  enddo

  !> exit
  call mpp_domains_exit()
  call MPI_finalize(ierr)

contains

  subroutine arr_init_2d(arr)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(out) :: arr(:,:)
    real(r8_kind) :: unif(size(arr,1), size(arr,2))
    type(randomNumberStream), save :: random_stream

    call getRandomNumbers(random_stream, unif)

    ! Workaround so that when FMS_TEST_TYPE_ is set to `integer`, it resolves to the `int`
    ! intrinsic when called as a function. The `real` keyword matches the name of its typecast
    ! function, so this workaround is only needed for the integral case.
#define integer int
    arr = FMS_TEST_TYPE_ (1e9_r8_kind * (unif - 0.5_r8_kind), FMS_TEST_KIND_)
#undef integer
  end subroutine arr_init_2d

  subroutine arr_init_3d(arr)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(out) :: arr(:,:,:)
    integer :: k

    do k = 1, size(arr, 3)
      call arr_init(arr(:, :, k))
    enddo
  end subroutine arr_init_3d

  subroutine arr_compare_2d(arr0, arr1, msg)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(in), dimension(:,:) :: arr0, arr1
    character(*), intent(in) :: msg

    if (any(arr0.ne.arr1)) then
      call mpp_error(FATAL, "Result from mpp_global_field (2D) does not agree with source data: " // msg)
    endif
  end subroutine arr_compare_2d

  subroutine arr_compare_3d(arr0, arr1, msg)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), intent(in), dimension(:,:,:) :: arr0, arr1
    character(*), intent(in) :: msg

    if (any(arr0.ne.arr1)) then
      call mpp_error(FATAL, "Result from mpp_global_field (3D) does not agree with source data: " // msg)
    endif
  end subroutine arr_compare_3d

  subroutine run_tests_2d(test_params, p)
    type(test_params_t), intent(in) :: test_params
    integer, intent(in) :: p !< Permutation of array indices (ranges from 1 to rank!)

    type(domain2D)  :: domain
    integer         :: i, j
    type(permutable_indices(2)) :: compute, data, global, global_with_halo, global_x, global_y
    integer, allocatable :: pelist(:)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), allocatable  :: global0(:,:), local(:,:), global1(:,:)

    !> set up domain
    call mpp_define_layout([1,nx,1,ny], npes, layout)
    call mpp_define_domains([1,nx,1,ny], layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name=trim(test_params%name), symmetry=test_params%symmetry)

    !> get compute domain
    call mpp_get_compute_domain(domain, compute%lb(1),  compute%ub(1),  compute%lb(2),  compute%ub(2))
    compute%ub = compute%ub + test_params%shift

    !> get data domain
    call mpp_get_data_domain(domain, data%lb(1), data%ub(1), data%lb(2), data%ub(2))
    data%ub = data%ub + test_params%shift

    global%lb = [1, 1]
    global%ub = [nx, ny] + test_params%shift

    global_with_halo%lb = global%lb - [whalo, shalo]
    global_with_halo%ub = global%ub + [ehalo, nhalo]

    global_x%lb = [global%lb(1), compute%lb(2)]
    global_x%ub = [global%ub(1), compute%ub(2)]

    global_y%lb = [compute%lb(1), global%lb(2)]
    global_y%ub = [compute%ub(1), global%ub(2)]

    call compute%permute(p)
    call data%permute(p)
    call global%permute(p)
    call global_with_halo%permute(p)
    call global_x%permute(p)
    call global_y%permute(p)

    !> assign global
    allocate(global0(global_with_halo%lb(1):global_with_halo%ub(1), global_with_halo%lb(2):global_with_halo%ub(2)))
    global0 = zero
    call arr_init(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2)))

    allocate(global1(global%lb(1):global%ub(1), global%lb(2):global%ub(2)))

    !> allocate for global domain
    allocate(local(data%lb(1):data%ub(1), data%lb(2):data%ub(2)))
    local(:,:) = global0(data%lb(1):data%ub(1), data%lb(2):data%ub(2))

    !> test the data on data domain
    global1 = zero
    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2)), global1, trim(test_params%name)//' mpp_global_field on data domain')

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set (only on pe of pelist1)
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate(pelist(0:layout(1)-1))
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=XUPDATE, position=test_params%position)
    call arr_compare(global0(global_x%lb(1):global_x%ub(1),global_x%lb(2):global_x%ub(2)), &
                     global1(global_x%lb(1):global_x%ub(1),global_x%lb(2):global_x%ub(2)), trim(test_params%name)// &
                     ' mpp_global_field xupdate only on data domain')

    !> yupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=YUPDATE, position=test_params%position)
    call arr_compare(global0(global_y%lb(1):global_y%ub(1),global_y%lb(2):global_y%ub(2)), &
                     global1(global_y%lb(1):global_y%ub(1),global_y%lb(2):global_y%ub(2)), trim(test_params%name)// &
                     ' mpp_global_field yupdate only on data domain')

    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2)), global1, trim(test_params%name)//' mpp_global_field on data domain')

    !> test the data on compute domain

    deallocate(local)
    allocate(local(compute%lb(1):compute%ub(1), compute%lb(2):compute%ub(2)))
    local(:,:) = global0(compute%lb(1):compute%ub(1), compute%lb(2):compute%ub(2))

    global1 = zero
    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2)), global1, trim(test_params%name)//' mpp_global_field on compute domain')

    !> xupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=XUPDATE, position=test_params%position)
    call arr_compare(global0(global_x%lb(1):global_x%ub(1),global_x%lb(2):global_x%ub(2)), &
                     global1(global_x%lb(1):global_x%ub(1),global_x%lb(2):global_x%ub(2)), trim(test_params%name)// &
                     ' mpp_global_field xupdate only on compute domain')

    !> yupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=YUPDATE, position=test_params%position)
    call arr_compare(global0(global_y%lb(1):global_y%ub(1),global_y%lb(2):global_y%ub(2)), &
                     global1(global_y%lb(1):global_y%ub(1),global_y%lb(2):global_y%ub(2)), trim(test_params%name)// &
                     ' mpp_global_field yupdate only on compute domain')
  end subroutine run_tests_2d

  subroutine run_tests_3d(test_params, p)
    type(test_params_t), intent(in) :: test_params
    integer, intent(in) :: p !< Permutation of array indices (ranges from 1 to rank!)

    type(domain2D)  :: domain
    integer         :: i, j
    type(permutable_indices(3)) :: compute, data, global, global_with_halo, global_x, global_y
    integer, allocatable :: pelist(:)
    FMS_TEST_TYPE_ (FMS_TEST_KIND_), allocatable  :: global0(:,:,:), local(:,:,:), global1(:,:,:)

    !> set up domain
    call mpp_define_layout([1,nx,1,ny], npes, layout)
    call mpp_define_domains([1,nx,1,ny], layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name=trim(test_params%name), symmetry=test_params%symmetry)

    !> get compute domain
    call mpp_get_compute_domain(domain, compute%lb(1),  compute%ub(1),  compute%lb(2),  compute%ub(2))
    compute%ub(1:2) = compute%ub(1:2) + test_params%shift
    compute%lb(3) = 1
    compute%ub(3) = nz

    !> get data domain
    call mpp_get_data_domain(domain, data%lb(1), data%ub(1), data%lb(2), data%ub(2))
    data%ub(1:2) = data%ub(1:2) + test_params%shift
    data%lb(3) = 1
    data%ub(3) = nz

    global%lb = [1, 1, 1]
    global%ub = [nx, ny, nz]
    global%ub(1:2) = global%ub(1:2) + test_params%shift

    global_with_halo%lb = global%lb - [whalo, shalo, 0]
    global_with_halo%ub = global%ub + [ehalo, nhalo, 0]

    global_x%lb = [global%lb(1), compute%lb(2), global%lb(3)]
    global_x%ub = [global%ub(1), compute%ub(2), global%ub(3)]

    global_y%lb = [compute%lb(1), global%lb(2), global%lb(3)]
    global_y%ub = [compute%ub(1), global%ub(2), global%ub(3)]

    call compute%permute(p)
    call data%permute(p)
    call global%permute(p)
    call global_with_halo%permute(p)
    call global_x%permute(p)
    call global_y%permute(p)

    !> assign global0
    allocate(global0(global_with_halo%lb(1):global_with_halo%ub(1), global_with_halo%lb(2):global_with_halo%ub(2), &
                     global_with_halo%lb(3):global_with_halo%ub(3)))

    global0 = zero
    call arr_init(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2), global%lb(3):global%ub(3)))

    allocate(global1(global%lb(1):global%ub(1), global%lb(2):global%ub(2), global%lb(3):global%ub(3)))

    !> for data domain
    allocate(local(data%lb(1):data%ub(1), data%lb(2):data%ub(2), data%lb(3):data%ub(3)))
    local(:,:,:) = global0(data%lb(1):data%ub(1), data%lb(2):data%ub(2), data%lb(3):data%ub(3))

    !> test the data on data domain
    global1 = zero
    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2), global%lb(3):global%ub(3)), &
                     global1, trim(test_params%name)//' mpp_global_field on data domain')

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set (only on pe of pelist1)
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate(pelist(0:layout(1)-1))
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=XUPDATE, position=test_params%position)
    call arr_compare(global0(global_x%lb(1):global_x%ub(1), global_x%lb(2):global_x%ub(2), &
                     global_x%lb(3):global_x%ub(3)), global1(global_x%lb(1):global_x%ub(1), &
                     global_x%lb(2):global_x%ub(2), global_x%lb(3):global_x%ub(3)),trim(test_params%name)// &
                         & ' mpp_global_field xupdate only on data domain')

    !> yupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=YUPDATE, position=test_params%position)
    call arr_compare(global0(global_y%lb(1):global_y%ub(1), global_y%lb(2):global_y%ub(2), &
                     global_y%lb(3):global_y%ub(3)), global1(global_y%lb(1):global_y%ub(1), &
                     global_y%lb(2):global_y%ub(2), global_y%lb(3):global_y%ub(3)),trim(test_params%name)// &
                     ' mpp_global_field yupdate only on data domain')

    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2), global%lb(3):global%ub(3)), &
                     global1,trim(test_params%name)//' mpp_global_field on data domain')

    !> test the data on compute domain

    deallocate(local)
    allocate(local(compute%lb(1):compute%ub(1), compute%lb(2):compute%ub(2), compute%lb(3):compute%ub(3)))
    local(:,:,:) = global0(compute%lb(1):compute%ub(1), compute%lb(2):compute%ub(2), compute%lb(3):compute%ub(3))

    global1 = zero
    call mpp_global_field(domain, local, global1, position=test_params%position)
    call arr_compare(global0(global%lb(1):global%ub(1), global%lb(2):global%ub(2), global%lb(3):global%ub(3)), &
                     global1, trim(test_params%name)//' mpp_global_field on compute domain')

    !> xupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=XUPDATE, position=test_params%position)
    call arr_compare(global0(global_x%lb(1):global_x%ub(1), global_x%lb(2):global_x%ub(2), &
                     global_x%lb(3):global_x%ub(3)), global1(global_x%lb(1):global_x%ub(1), &
                     global_x%lb(2):global_x%ub(2), global_x%lb(3):global_x%ub(3)), &
                     trim(test_params%name)//' mpp_global_field xupdate only on compute domain')

    !> yupdate
    global1 = zero
    call mpp_global_field(domain, local, global1, flags=YUPDATE, position=test_params%position)
    call arr_compare(global0(global_y%lb(1):global_y%ub(1), global_y%lb(2):global_y%ub(2), &
                     global_y%lb(3):global_y%ub(3)), global1(global_y%lb(1):global_y%ub(1), &
                     global_y%lb(2):global_y%ub(2), global_y%lb(3):global_y%ub(3)), &
                     trim(test_params%name)//' mpp_global_field yupdate only on compute domain')
  end subroutine run_tests_3d
end program test_mpp_global_field
