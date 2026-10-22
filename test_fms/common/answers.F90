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

!> Module for generating and checking reference answers. Example usage:
!! use answers_mod
!! call answers_mod_init(20) ! Maximum of 20 reference answers
!! call fancy_calculation1(..., x1)
!! call register_answer(x1, "Output from fancy_calculation1")
!! call fancy_calculation2(..., x2)
!! call register_answer(x2, "Output from fancy_calculation2")
!! call answers_mod_end
!!
!! If `answers_nml` is present in `input.nml`, the answers encoded therein will be verified,
!! and a fatal error issued if any discrepancy is found.
!! If `answers_nml` is not present, a file named `answers.nml` will be generated containing checksums
!! of all registered answers. In this generation mode, all registered answers are printed to standard
!! output so that they can be manually validated.

module answers_mod
  use platform_mod, only: r8_kind, i8_kind
  use mpp_mod, only: input_nml_file, mpp_error, mpp_pe, mpp_root_pe, mpp_npes, mpp_comm, NOTE, FATAL
  use fms_mod, only: check_nml_error
  use fms_string_utils_mod, only: string
  use mpi_f08, only: mpi_gather, MPI_INTEGER8

  implicit none

  integer(i8_kind), allocatable :: ref(:), answers(:,:)
  logical :: gen_ref_answers
  integer :: n=0, nmax

  namelist /answers_nml/ answers

  contains

  function calc_chksum(arr) result(chksum)
    real(r8_kind), intent(in) :: arr(:,:,:)
    integer(i8_kind) :: chksum
    integer(i8_kind) :: arr_int8(size(arr))
    integer :: i

    arr_int8 = transfer(arr, arr_int8)

    chksum = 0_i8_kind
    do i=1,size(arr)
      chksum = ishftc(chksum, 1)
      chksum = ieor(chksum, arr_int8(i))
    enddo
  end function calc_chksum

  subroutine answers_mod_init(nmax_in)
    integer, intent(in) :: nmax_in !< Maximum number of answers
    integer :: ierr

    nmax = nmax_in

    allocate(answers(nmax,mpp_npes()))
    allocate(ref(nmax))
    ref = 0

    read (input_nml_file, nml=answers_nml, iostat=ierr)
    if (check_nml_error(ierr, "answers_nml").ne.0) then
      gen_ref_answers = .true.
      call mpp_error(NOTE, "Failed to read answers namelist. Generating reference answers instead.")
    else
      gen_ref_answers = .false.
      ref = answers(:, mpp_pe()+1)
    endif
  end subroutine answers_mod_init

  subroutine answers_mod_end
    integer :: answers_out

    if (gen_ref_answers) then
      call mpi_gather(ref, nmax, MPI_INTEGER8, answers, nmax, MPI_INTEGER8, 0, mpp_comm())

      if (mpp_pe().eq.mpp_root_pe()) then
        open(newunit=answers_out, file="answers.nml", status="replace", action="write")
        write(answers_out, nml=answers_nml)
        close(answers_out)
      endif
    endif
  end subroutine answers_mod_end

  subroutine register_answer(arr, msg)
    real(r8_kind), dimension(:,:,:), intent(in) :: arr
    character(*), intent(in) :: msg
    integer :: ierr
    integer(i8_kind) :: chksum

    n = n+1

    if (n.gt.nmax) then
      call mpp_error(FATAL, "register_answer: Too many answers registered")
    endif

    chksum = calc_chksum(arr)

    if (gen_ref_answers) then
      ref(n) = chksum
      call print_matrix(arr(:,:,1), msg)
    else
      if (chksum.ne.ref(n)) then
        call mpp_error(FATAL, "Failed checksum: "//msg)
      endif
    endif
  end subroutine register_answer

  subroutine print_matrix(x, msg)
    real(r8_kind), dimension(:,:), intent(in) :: x
    character(*), intent(in) :: msg
    integer :: i,j
    integer :: iostat

    print "(A)", "========================================================="
    print "(A)", msg
    print "(A)", "================ Begin matrix dump ======================"

    do j=size(x,2),1,-1
      do i=1,size(x,1)
        write(*,'(" ",F15.9)', advance="no", iostat=iostat) x(i,j)
      enddo
      write(*, "(A)") ""
    enddo

    print "(A)", "================ End matrix dump ========================"
    print "(A)", ""
  end subroutine print_matrix
end module answers_mod
