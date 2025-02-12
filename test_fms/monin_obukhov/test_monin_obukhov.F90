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

! Check monin_obukhov_mod calculations against an array of known answer keys.
! Each answer key should correspond to a particular hardware/compiler/flags combination.

! Express a real array as an integer array via transfer(), and reshape the result
! to match the shape of the original data
#define INT_(arr)  reshape(transfer(arr, [mi]), shape(arr))

! Promote a dimension(n) array to dimension(n,n) by making n copies of the original data
#define ARR_2D_(arr) spread(arr, 2, size(arr))

! Promote a dimension(n) array to dimension(n,n,n) by making n*n copies of the original data
#define ARR_3D_(arr) spread(ARR_2D_(arr), 3, size(arr))

program test_monin_obukhov
  use monin_obukhov_mod
  use mpp_mod, only: mpp_error, FATAL, stdout, mpp_init, mpp_exit, input_nml_file
  use fms_mod, only: check_nml_error
  use platform_mod, only: r4_kind, r8_kind, i4_kind, i8_kind
  use fms_string_utils_mod, only: string

  implicit none

#if MO_TEST_KIND_ == 4
  integer, parameter :: kr = r4_kind
  integer, parameter :: ki = i4_kind
#else
  integer, parameter :: kr = r8_kind
  integer, parameter :: ki = i8_kind
#endif

  integer(ki), parameter :: mi = 0_ki !< Mold for transfer() intrinsic

  !< Shape of 1D arrays passed to monin_obukhov_mod subroutines
  integer, parameter :: n_1d = 5

  integer :: n_answers !< Number of known answer keys
  namelist /metaparams_nml/ n_answers

  type drag_input_t
    real(kr), dimension(n_1d) :: &
      pt     = [268.559120403867_kr, 269.799228886728_kr, 277.443023238556_kr, &
             & 295.79192777341_kr, 293.268717243262_kr], &
      pt0    = [273.42369841804_kr , 272.551410044203_kr, 278.638168565727_kr, &
             & 298.133068766049_kr, 292.898163706587_kr], &
      z      = [29.432779269303_kr, 30.0497139076724_kr, 31.6880000418153_kr, &
             & 34.1873479240475_kr, 33.2184943356517_kr], &
      z0     = [5.86144925739178e-05_kr, 0.0001_kr, 0.000641655193293549_kr, &
             & 3.23383768877187e-05_kr, 0.07_kr], &
      zt     = [3.69403636275411e-05_kr, 0.0001_kr, 1.01735489109205e-05_kr, &
             & 7.63933834969505e-05_kr, 0.00947346982656289_kr], &
      zq     = [5.72575636226887e-05_kr, 0.0001_kr, 5.72575636226887e-05_kr, &
             & 5.72575636226887e-05_kr, 5.72575636226887e-05_kr], &
      speed  = [2.9693638452068_kr, 2.43308757772094_kr, 5.69418282305367_kr, &
             & 9.5608693754561_kr, 4.35302260074334_kr]

    logical, dimension(n_1d) :: avail = [.true., .true., .true., .true., .true.]
  end type

  type stable_mix_input_t
    real(kr), dimension(n_1d) :: rich = [1650.92431853365_kr, 1650.9256285137_kr, &
    & 77.7636819036559_kr, 1.92806556391324_kr, 0.414767442012442_kr]
  end type

  type diff_input_t
    real(kr) :: z      = 19.9982554527751_kr, &
              & u_star = 0.129638955971075_kr, &
              & b_star = 0.000991799765557209_kr
  end type

  type profile_input_t
    real(kr) :: zref   = 10._kr, &
              & zref_t = 2._kr

    real(kr), dimension(n_1d) :: &
      & z      = [29.432779269303_kr, 30.0497139076724_kr, 31.6880000418153_kr, &
               & 34.1873479240475_kr, 33.2184943356517_kr], &
      & z0     = [5.86144925739178e-05_kr, 0.0001_kr, 0.000641655193293549_kr, &
               & 3.23383768877187e-05_kr, 0.07_kr], &
      & zt     = [3.69403636275411e-05_kr, 0.0001_kr, 1.01735489109205e-05_kr, &
               & 7.63933834969505e-05_kr, 0.00947346982656289_kr], &
      & zq     = [5.72575636226887e-05_kr, 0.0001_kr, 5.72575636226887e-05_kr, &
               & 5.72575636226887e-05_kr, 5.72575636226887e-05_kr], &
      & u_star = [0.109462510724615_kr, 0.0932942802513508_kr, 0.223232887323184_kr, &
               & 0.290918439028557_kr, 0.260087579361467_kr], &
      & b_star = [0.00690834676781433_kr, 0.00428178089592372_kr, 0.00121229800895103_kr, &
               & 0.00262353784027441_kr, -0.000570314880866852_kr], &
      & q_star = [0.000110861442197537_kr, 9.44983279664197e-05_kr, 4.17643828631936e-05_kr, &
               & 0.000133135421415819_kr, 9.36317815993945e-06_kr]

    logical :: avail(n_1d) = [.true., .true., .true., .true., .true.]
  end type

  type drag_answers_t
    integer(ki), dimension(n_1d) :: drag_m, drag_t, drag_q, u_star, b_star
  end type

  type stable_mix_answers_t
    integer(ki), dimension(n_1d) :: mix
  end type

  type diff_answers_t
    integer(ki) :: k_m, k_h
  end type

  type profile_answers_t
    integer(ki), dimension(n_1d) :: del_m, del_t, del_q
  end type

  type(drag_input_t), parameter :: drag_input = drag_input_t() !< Input arguments for mo_drag
  type(stable_mix_input_t), parameter :: stable_mix_input = stable_mix_input_t() !< Input arguments for stable_mix
  type(diff_input_t), parameter :: diff_input = diff_input_t() !< Input arguments for mo_diff
  type(profile_input_t), parameter :: profile_input = profile_input_t() !< Input arguments for mo_profile

  ! Entries 1:n of the arrays below contain known answer keys. Entry n+1 contains
  ! the answers that we calculate. Represent answer data using integral arrays,
  ! because Fortran does not guarantee bit-for-bit exactness of real values
  ! stored in namelist files.
  type(drag_answers_t),       allocatable :: drag_answers(:) !< mo_drag answers
  type(stable_mix_answers_t), allocatable :: stable_mix_answers(:) !< stable_mix answers
  type(diff_answers_t),       allocatable :: diff_answers(:) !< mo_diff answers
  type(profile_answers_t),    allocatable :: profile_answers(:) !< mo_profile answers

  namelist /answers_nml/ drag_answers, stable_mix_answers, diff_answers, profile_answers

  call mpp_init

  call monin_obukhov_init
  call read_answers
  call calc_answers

  if (.not.check_answers()) then
    call write_answers
    call mpp_error(FATAL, "monin_obukhov unit tests did not pass with any known answer key")
  endif

  call mpp_exit

  contains

    !< Read answer keys from input.nml
    subroutine read_answers
      integer :: io, ierr

      read (input_nml_file, nml=metaparams_nml, iostat=io)
      ierr = check_nml_error(io, "metaparams_nml")

      allocate(drag_answers(n_answers+1))
      allocate(stable_mix_answers(n_answers+1))
      allocate(diff_answers(n_answers+1))
      allocate(profile_answers(n_answers+1))

      if (n_answers.gt.0) then
        read (input_nml_file, nml=answers_nml, iostat=io)
        ierr = check_nml_error(io, "answers_nml")
      endif
    end subroutine

    !> Store existing answer keys, as well as the answers just calculated, in an
    !> output file.
    subroutine write_answers
      character(:), allocatable :: filename
      integer :: fh

      filename = "OUT.r" // string(MO_TEST_KIND_) // ".nml"
      print "(A)", "Writing newly generated answer key to " // filename

      n_answers = n_answers + 1

      open (newunit=fh, file=filename)
      write (fh, nml=metaparams_nml)
      write (fh, nml=answers_nml)
      close (fh)
    end subroutine

    !> Calculate all answers
    subroutine calc_answers
      call calc_answers_drag
      call calc_answers_stable_mix
      call calc_answers_diff
      call calc_answers_profile
    end subroutine

    !> Calculate 1D answers for mo_drag, and assert that all 2D answers must agree
    !> with the corresponding 1D answers
    subroutine calc_answers_drag
      real(kr), dimension(n_1d) :: drag_m_1d, drag_t_1d, drag_q_1d, u_star_1d, b_star_1d
      real(kr), dimension(n_1d, n_1d) :: drag_m_2d, drag_t_2d, drag_q_2d, u_star_2d, b_star_2d

      drag_m_1d = 0._kr
      drag_t_1d = 0._kr
      drag_q_1d = 0._kr
      u_star_1d = 0._kr
      b_star_1d = 0._kr

      drag_m_2d = 0._kr
      drag_t_2d = 0._kr
      drag_q_2d = 0._kr
      u_star_2d = 0._kr
      b_star_2d = 0._kr

      associate (in => drag_input)
        call mo_drag(in%pt, in%pt0, in%z, in%z0, in%zt, in%zq, in%speed, &
                   & drag_m_1d, drag_t_1d, drag_q_1d, u_star_1d, b_star_1d, in%avail)

        call mo_drag(ARR_2D_(in%pt), ARR_2D_(in%pt0), ARR_2D_(in%z), ARR_2D_(in%z0), &
                   & ARR_2D_(in%zt), ARR_2D_(in%zq), ARR_2D_(in%speed), &
                   & drag_m_2d, drag_t_2d, drag_q_2d, u_star_2d, b_star_2d)
      end associate

      associate(ans => drag_answers(n_answers+1))
        ans%drag_m = INT_(drag_m_1d)
        ans%drag_t = INT_(drag_t_1d)
        ans%drag_q = INT_(drag_q_1d)
        ans%u_star = INT_(u_star_1d)
        ans%b_star = INT_(b_star_1d)

        call answer_validate_2d(ans%drag_m, drag_m_2d)
        call answer_validate_2d(ans%drag_t, drag_t_2d)
        call answer_validate_2d(ans%drag_q, drag_q_2d)
        call answer_validate_2d(ans%u_star, u_star_2d)
        call answer_validate_2d(ans%b_star, b_star_2d)
      end associate
    end subroutine

    !> Calculate 1D answers for stable_mix, and assert that all 2D and 3D answers
    !> must agree with the corresponding 1D answers
    subroutine calc_answers_stable_mix
      real(kr), dimension(n_1d) :: mix_1d
      real(kr), dimension(n_1d, n_1d) :: mix_2d
      real(kr), dimension(n_1d, n_1d, n_1d) :: mix_3d

      mix_1d = 0._kr
      mix_2d = 0._kr
      mix_3d = 0._kr

      associate (in => stable_mix_input)
        call stable_mix(in%rich, mix_1d)
        call stable_mix(ARR_2D_(in%rich), mix_2d)
        call stable_mix(ARR_3D_(in%rich), mix_3d)
      end associate

      associate (ans => stable_mix_answers(n_answers+1))
        ans%mix = INT_(mix_1d)

        call answer_validate_2d(ans%mix, mix_2d)
        call answer_validate_3d(ans%mix, mix_3d)
      end associate
    end subroutine

    !> Calculate answers for mo_diff
    subroutine calc_answers_diff
      real(kr), dimension(1,1,1) :: k_m, k_h

      k_m = 0._kr
      k_h = 0._kr

      associate (in => diff_input)
        ! mo_diff_0d_1
        call mo_diff(in%z, in%u_star, in%b_star, k_m(1,1,1), k_h(1,1,1))

        associate (ans => diff_answers(n_answers+1))
          ans%k_m = transfer(k_m(1,1,1), mi)
          ans%k_h = transfer(k_h(1,1,1), mi)
        end associate

        ! mo_diff_0d_n
        call mo_diff([in%z], in%u_star, in%b_star, k_m(:,1,1), k_h(:,1,1))
        call diff_check(k_m, k_h, "mo_diff_0d_n")

        ! mo_diff_1d_1
        call mo_diff([in%z], [in%u_star], [in%b_star], k_m(:,1,1), k_h(:,1,1))
        call diff_check(k_m, k_h, "mo_diff_1d_1")

        ! mo_diff_1d_n
        call mo_diff(ARR_2D_([in%z]), [in%u_star], [in%b_star], k_m(:,:,1), k_h(:,:,1))
        call diff_check(k_m, k_h, "mo_diff_1d_n")

        ! mo_diff_2d_1
        call mo_diff(ARR_2D_([in%z]), ARR_2D_([in%u_star]), ARR_2D_([in%b_star]), k_m(:,:,1), k_h(:,:,1))
        call diff_check(k_m, k_h, "mo_diff_2d_1")

        ! mo_diff_2d_n
        call mo_diff(ARR_3D_([in%z]), ARR_2D_([in%u_star]), ARR_2D_([in%b_star]), k_m(:,:,:), k_h(:,:,:))
        call diff_check(k_m, k_h, "mo_diff_2d_n")
      end associate
    end subroutine

    subroutine diff_check(k_m, k_h, label)
      real(kr), dimension(1,1,1), intent(in) :: k_m, k_h
      character(*), intent(in) :: label

      associate (ans => diff_answers(n_answers+1))
        if (ans%k_m .ne. transfer(k_m(1,1,1), mi)) then
          call mpp_error(FATAL, label // " test failed: k_m value differs from that of mo_diff_0d_1")
        endif

        if (ans%k_h .ne. transfer(k_h(1,1,1), mi)) then
          call mpp_error(FATAL, label // " test failed: k_h value differs from that of mo_diff_0d_1")
        endif
      end associate
    end subroutine

    !> Calculate 1D answers for mo_profile, and assert that all 2D answers must
    !> agree with the corresponding 1D answers
    subroutine calc_answers_profile
      real(kr), dimension(n_1d)       :: del_m_1d, del_t_1d, del_q_1d
      real(kr), dimension(n_1d, n_1d) :: del_m_2d, del_t_2d, del_q_2d

      del_m_1d = 0._kr
      del_t_1d = 0._kr
      del_q_1d = 0._kr

      del_m_2d = 0._kr
      del_t_2d = 0._kr
      del_q_2d = 0._kr

      associate (in => profile_input)
        call mo_profile(in%zref, in%zref_t, in%z, in%z0, in%zt, in%zq, &
                      & in%u_star, in%b_star, in%q_star, &
                      & del_m_1d, del_t_1d, del_q_1d, in%avail)

        call mo_profile(in%zref, in%zref_t, ARR_2D_(in%z), &
                      & ARR_2D_(in%z0), ARR_2D_(in%zt), ARR_2D_(in%zq), &
                      & ARR_2D_(in%u_star), ARR_2D_(in%b_star), ARR_2D_(in%q_star), &
                      & del_m_2d, del_t_2d, del_q_2d)
      end associate

      associate (ans => profile_answers(n_answers+1))
        ans%del_m = INT_(del_m_1d)
        ans%del_t = INT_(del_t_1d)
        ans%del_q = INT_(del_q_1d)

        call answer_validate_2d(ans%del_m, del_m_2d)
        call answer_validate_2d(ans%del_t, del_t_2d)
        call answer_validate_2d(ans%del_q, del_q_2d)
      end associate
    end subroutine

    !> Check whether the calculated answers agree with a known answer key
    function check_answers() result(res)
      logical :: res
      integer :: i !< Answer key index

      res = .true.

      do i=1, n_answers
        if(check_answer_key(i)) then
          print "(A)", "monin_obukhov tests passed with answer key " // string(i)
          return
        endif
      enddo

      res = .false.
    end function

    !> Check whether the calculated answers agree with answer key i
    function check_answer_key(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      res = check_drag(i) .and. check_stable_mix(i) .and. check_diff(i) .and. check_profile(i)
    end function

    !> Check whether the calculated mo_drag answers agree with answer key i
    function check_drag(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => drag_answers(i), ans1 => drag_answers(n_answers+1))
        res = array_compare_1d(ans0%drag_m, ans1%drag_m) .and. &
              array_compare_1d(ans0%drag_t, ans1%drag_t) .and. &
              array_compare_1d(ans0%drag_q, ans1%drag_q) .and. &
              array_compare_1d(ans0%u_star, ans1%u_star) .and. &
              array_compare_1d(ans0%b_star, ans1%b_star)
      end associate
    end function

    !> Check whether the calculated stable_mix answers agree with answer key i
    function check_stable_mix(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => stable_mix_answers(i), ans1 => stable_mix_answers(n_answers+1))
        res = array_compare_1d(ans0%mix, ans1%mix)
      end associate
    end function

    !> Check whether the calculated mo_diff answers agree with answer key i
    function check_diff(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => diff_answers(i), ans1 => diff_answers(n_answers+1))
        res = (ans0%k_m.eq.ans1%k_m) .and. (ans0%k_h.eq.ans1%k_h)
      end associate
    end function

    !> Check whether the calculated mo_profile answers agree with answer key i
    function check_profile(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => profile_answers(i), ans1 => profile_answers(n_answers+1))
        res = array_compare_1d(ans0%del_m, ans1%del_m) .and. &
            & array_compare_1d(ans0%del_t, ans1%del_t) .and. &
            & array_compare_1d(ans0%del_q, ans1%del_q)
      end associate
    end function

    !< Check whether a pair of integral 1D arrays are equal
    function array_compare_1d(arr1, arr2) result(res)
      integer(ki), intent(in) :: arr1(:), arr2(:)
      logical :: res
      integer :: n, i

      res = .false.

      n = size(arr1, 1)
      if (size(arr2, 1) .ne. n) return

      do i=1, n
        if (arr1(i) .ne. arr2(i)) return
      enddo

      res = .true.
    end function

    !< Check whether a pair of integral 2D arrays are equal
    function array_compare_2d(arr1, arr2) result(res)
      integer(ki), intent(in) :: arr1(:,:), arr2(:,:)
      logical :: res
      integer :: n, i

      res = .false.

      n = size(arr1, 2)
      if (size(arr2, 2) .ne. n) return

      do i=1, n
        if (.not.array_compare_1d(arr1(:, i), arr2(:, i))) return
      enddo

      res = .true.
    end function

    !< Check whether a pair of integral 3D arrays are equal
    function array_compare_3d(arr1, arr2) result(res)
      integer(ki), intent(in) :: arr1(:,:,:), arr2(:,:,:)
      logical :: res
      integer :: n, i

      res = .false.

      n = size(arr1, 3)
      if (size(arr2, 3) .ne. n) return

      do i=1, n
        if (.not.array_compare_2d(arr1(:, :, i), arr2(:, :, i))) return
      enddo

      res = .true.
    end function

    ! Compare an integral 1D reference key array against a real-valued, 2D answer array
    subroutine answer_validate_2d(ref, arr)
      integer(ki), dimension(:), intent(in) :: ref
      real(kr), dimension(:,:), intent(in) :: arr
      integer :: i, n

      n = size(ref)

      if (size(arr, 1).ne.n .or. size(arr, 2).ne.n) then
        call mpp_error(FATAL, "Incorrect array shape")
      endif

      do i = 1,n
        if (.not.array_compare_1d(ref, transfer(arr(:,i), [mi]))) then
          call mpp_error(FATAL, "Array does not match reference value")
        endif
      enddo
    end subroutine

    ! Compare an integral 1D reference key array against a real-valued, 3D answer array
    subroutine answer_validate_3d(ref, arr)
      integer(ki), dimension(:), intent(in) :: ref
      real(kr), dimension(:,:,:), intent(in) :: arr
      integer :: i, j, n

      n = size(ref)

      if (size(arr, 1).ne.n .or. size(arr, 2).ne.n .or. size(arr, 3).ne.n) then
        call mpp_error(FATAL, "Incorrect array shape")
      endif

      do j = 1,n
        do i = 1,n
          if (.not.array_compare_1d(ref, transfer(arr(:,i,j), [mi]))) then
            call mpp_error(FATAL, "Array does not match reference value")
          endif
        enddo
      enddo
    end subroutine
end program test_monin_obukhov
