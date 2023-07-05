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

! Check monin_obukhov_mod calculations against a dictionary of reference answers.

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

  integer(ki), parameter :: mi(1) = [0_ki]

#define INT_(arr)  reshape(transfer(arr, mi), shape(arr))

  integer, parameter :: n_1d = 5, &
                      & diff_ni = 1, &
                      & diff_nj = 1, &
                      & diff_nk = 1

  integer :: n_answers
  namelist /metaparams_nml/ n_answers

  type drag_input_t
    real(ki), dimension(n_1d) :: pt, pt0, z, z0, zt, zq, speed
    logical, dimension(n_1d) :: avail
  end type

  type stable_mix_input_t
    real(ki), dimension(n_1d) :: rich
  end type

  type diff_input_t
    real(ki), dimension(diff_ni, diff_nj, diff_nk) :: z
    real(ki), dimension(diff_ni, diff_nj)          :: u_star, b_star
  end type

  type profile_input_t
    real(ki) :: zref, zref_t
    real(ki), dimension(n_1d) :: z, z0, zt, zq, u_star, b_star, q_star
    logical :: avail(n_1d)
  end type

  type drag_answers_t
    integer(ki), dimension(n_1d) :: drag_m, drag_t, drag_q, u_star, b_star
  end type

  type stable_mix_answers_t
    integer(ki), dimension(n_1d) :: mix
  end type

  type diff_answers_t
    integer(ki), dimension(diff_ni, diff_nj, diff_nk) :: k_m, k_h
  end type

  type profile_answers_t
    integer(ki), dimension(n_1d) :: del_m, del_t, del_q
  end type

  type(drag_input_t)       :: drag_input
  type(stable_mix_input_t) :: stable_mix_input
  type(diff_input_t)       :: diff_input
  type(profile_input_t)    :: profile_input

  type(drag_answers_t),       allocatable :: drag_answers(:)
  type(stable_mix_answers_t), allocatable :: stable_mix_answers(:)
  type(diff_answers_t),       allocatable :: diff_answers(:)
  type(profile_answers_t),    allocatable :: profile_answers(:)

  namelist /answers_nml/ drag_answers, stable_mix_answers, diff_answers, profile_answers

  call mpp_init

  call monin_obukhov_init
  call set_input_params
  call read_answers
  call calc_answers

  if (.not.check_answers()) then
    call write_answers
    call mpp_error(FATAL, "monin_obukhov unit tests did not pass with any known answer key")
  endif

  call mpp_exit

  contains

    subroutine set_input_params
      drag_input%pt     = [268.559120403867_kr, 269.799228886728_kr, 277.443023238556_kr, &
                         & 295.79192777341_kr, 293.268717243262_kr]
      drag_input%pt0    = [273.42369841804_kr , 272.551410044203_kr, 278.638168565727_kr, &
                         & 298.133068766049_kr, 292.898163706587_kr]
      drag_input%z      = [29.432779269303_kr, 30.0497139076724_kr, 31.6880000418153_kr, &
                         & 34.1873479240475_kr, 33.2184943356517_kr]
      drag_input%z0     = [5.86144925739178e-05_kr, 0.0001_kr, 0.000641655193293549_kr, &
                         & 3.23383768877187e-05_kr, 0.07_kr]
      drag_input%zt     = [3.69403636275411e-05_kr, 0.0001_kr, 1.01735489109205e-05_kr, &
                         & 7.63933834969505e-05_kr, 0.00947346982656289_kr]
      drag_input%zq     = [5.72575636226887e-05_kr, 0.0001_kr, 5.72575636226887e-05_kr, &
                         & 5.72575636226887e-05_kr, 5.72575636226887e-05_kr]
      drag_input%speed  = [2.9693638452068_kr, 2.43308757772094_kr, 5.69418282305367_kr, &
                         & 9.5608693754561_kr, 4.35302260074334_kr]
      drag_input%avail  = [.true., .true., .true., .true., .true.]

      stable_mix_input%rich = [1650.92431853365_kr, 1650.9256285137_kr, 77.7636819036559_kr, &
                             & 1.92806556391324_kr, 0.414767442012442_kr]

      profile_input%zref   = 10._kr
      profile_input%zref_t = 2._kr

      profile_input%z      = [29.432779269303_kr, 30.0497139076724_kr, 31.6880000418153_kr, &
                            & 34.1873479240475_kr, 33.2184943356517_kr]
      profile_input%z0     = [5.86144925739178e-05_kr, 0.0001_kr, 0.000641655193293549_kr, &
                            & 3.23383768877187e-05_kr, 0.07_kr]
      profile_input%zt     = [3.69403636275411e-05_kr, 0.0001_kr, 1.01735489109205e-05_kr, &
                            & 7.63933834969505e-05_kr, 0.00947346982656289_kr]
      profile_input%zq     = [5.72575636226887e-05_kr, 0.0001_kr, 5.72575636226887e-05_kr, &
                            & 5.72575636226887e-05_kr, 5.72575636226887e-05_kr]
      profile_input%u_star = [0.109462510724615_kr, 0.0932942802513508_kr, 0.223232887323184_kr, &
                            & 0.290918439028557_kr, 0.260087579361467_kr]
      profile_input%b_star = [0.00690834676781433_kr, 0.00428178089592372_kr, 0.00121229800895103_kr, &
                            & 0.00262353784027441_kr, -0.000570314880866852_kr]
      profile_input%q_star = [0.000110861442197537_kr, 9.44983279664197e-05_kr, 4.17643828631936e-05_kr, &
                            & 0.000133135421415819_kr, 9.36317815993945e-06_kr]
      profile_input%avail = [.true., .true., .true., .true., .true.]

      diff_input%z      = reshape([19.9982554527751_kr], shape(diff_input%z))
      diff_input%u_star = reshape([0.129638955971075_kr], shape(diff_input%u_star))
      diff_input%b_star = reshape([0.000991799765557209_kr], shape(diff_input%b_star))
    end subroutine

    subroutine read_answers
      integer :: io, ierr

      read (input_nml_file, nml=metaparams_nml, iostat=io)
      ierr = check_nml_error(io, "metaparams_nml")

      allocate(drag_answers_t       :: drag_answers(n_answers+1))
      allocate(stable_mix_answers_t :: stable_mix_answers(n_answers+1))
      allocate(diff_answers_t       :: diff_answers(n_answers+1))
      allocate(profile_answers_t    :: profile_answers(n_answers+1))

      if (n_answers.gt.0) then
        read (input_nml_file, nml=answers_nml, iostat=io)
        ierr = check_nml_error(io, "answers_nml")
      endif
    end subroutine

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

    function check_answers() result(res)
      logical :: res
      integer :: i

      res = .true.

      do i=1, n_answers
        if(check_all(i)) then
          print "(A)", "monin_obukhov tests passed with answer key " // string(i)
          return
        endif
      enddo

      res = .false.
    end function

    subroutine calc_answers_drag
      real(kr), dimension(n_1d) :: drag_m, drag_t, drag_q, u_star, b_star

      drag_m = 0._kr
      drag_t = 0._kr
      drag_q = 0._kr
      u_star = 0._kr
      b_star = 0._kr

      associate (in => drag_input)
        call mo_drag(in%pt, in%pt0, in%z, in%z0, in%zt, in%zq, in%speed, &
                   & drag_m, drag_t, drag_q, u_star, b_star, drag_input%avail)
      end associate

      associate(ans => drag_answers(n_answers+1))
        ans%drag_m = INT_(drag_m)
        ans%drag_t = INT_(drag_t)
        ans%drag_q = INT_(drag_q)
        ans%u_star = INT_(u_star)
        ans%b_star = INT_(b_star)
      end associate
    end subroutine

    subroutine calc_answers_stable_mix
      real(kr), dimension(n_1d) :: mix

      mix = 0._kr

      associate (in => stable_mix_input)
        call stable_mix(in%rich, mix)
      end associate

      associate (ans => stable_mix_answers(n_answers+1))
        ans%mix = INT_(mix)
      end associate
    end subroutine

    subroutine calc_answers_diff
      real(kr), dimension(diff_ni, diff_nj, diff_nk) :: k_m, k_h

      k_m = 0._kr
      k_h = 0._kr

      associate (in => diff_input)
        call mo_diff(in%z, in%u_star, in%b_star, k_m, k_h)
      end associate

      associate (ans => diff_answers(n_answers+1))
        ans%k_m = INT_(k_m)
        ans%k_h = INT_(k_h)
      end associate
    end subroutine

    subroutine calc_answers_profile
      real(kr), dimension(n_1d) :: del_m, del_t, del_q

      del_m = 0._kr
      del_t = 0._kr
      del_q = 0._kr

      associate (in => profile_input)
        call mo_profile(in%zref, in%zref_t, in%z, in%z0, in%zt, in%zq, &
                      & in%u_star, in%b_star, in%q_star, &
                      & del_m, del_t, del_q, in%avail)
      end associate

      associate (ans => profile_answers(n_answers+1))
        ans%del_m = INT_(del_m)
        ans%del_t = INT_(del_t)
        ans%del_q = INT_(del_q)
      end associate
    end subroutine

    subroutine calc_answers
      call calc_answers_drag
      call calc_answers_stable_mix
      call calc_answers_diff
      call calc_answers_profile
    end subroutine

    function check_all(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      res = check_drag(i) .and. check_stable_mix(i) .and. check_diff(i) .and. check_profile(i)
    end function

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

    function check_stable_mix(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => stable_mix_answers(i), ans1 => stable_mix_answers(n_answers+1))
        res = array_compare_1d(ans0%mix, ans1%mix)
      end associate
    end function

    function check_diff(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => diff_answers(i), ans1 => diff_answers(n_answers+1))
        res = array_compare_3d(ans0%k_m, ans1%k_m) .and. &
            & array_compare_3d(ans0%k_h, ans1%k_h)
      end associate
    end function

    function check_profile(i) result(res)
      integer, intent(in) :: i !< Answer key to check against
      logical :: res

      associate (ans0 => profile_answers(i), ans1 => profile_answers(n_answers+1))
        res = array_compare_1d(ans0%del_m, ans1%del_m) .and. &
            & array_compare_1d(ans0%del_t, ans1%del_t) .and. &
            & array_compare_1d(ans0%del_q, ans1%del_q)
      end associate
    end function

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
end program test_monin_obukhov
