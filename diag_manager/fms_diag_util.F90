module fms_diag_util_mod
!> \author Ganga Purja Pun
!> \email gagna.purjapun@noaa.gov
!! \brief Contains routines for the modern diag_manager
!!
!! \description
#ifdef use_yaml
  use platform_mod
  use fms_mod, only: fms_error_handler
  implicit none
  private

  public :: compare_two_sets_of_bounds, real_copy_set, check_indices_order, init_mask_3d

  contains

  !> @brief Compares the corresponding bounding indices of the first set with the second set.
  !> @return .TRUE. if any comparison returns true; i.e. the box bounded by the indices of the first set
  !! is out side the box bounded by the indices of the second set.
  LOGICAL FUNCTION compare_two_sets_of_bounds(bounds_a, bounds_b, error_str)
    integer, intent(in) :: bounds_a(:) !< First array with order: (/imin, imax, jmin, jmax, kmin, kmax/)
    integer, intent(in) :: bounds_b(:) !< Second array with the same order as the first
    character(*), intent(out) :: error_str

    compare_two_sets_of_bounds = .FALSE.

    if (size(bounds_a) .ne. size(bounds_b)) then
      compare_two_sets_of_bounds = .TRUE.
      error_str = 'diag_util_mod::compare_two_sets_of_bounds Error: sizes of sets do not match'
    else
      if ((size(bounds_a) .ne. 6) .and. (size(bounds_b) .ne. 6)) then
        compare_two_sets_of_bounds = .TRUE.
        error_str = 'diag_util_mod::compare_two_sets_of_bounds Error: sizes of sets must be 6'
      end if
    end if

    IF (bounds_a(1) .lt. bounds_b(1) .OR. bounds_a(2) .gt. bounds_b(2) .OR. &
      bounds_a(3) .lt. bounds_b(3) .OR. bounds_a(4) .gt. bounds_b(4) .OR. &
      bounds_a(5) .lt. bounds_b(5) .OR. bounds_a(6) .gt. bounds_b(6)) THEN
      compare_two_sets_of_bounds = .TRUE.
      error_str ='First set of bounds=   :   ,   :   ,   :     Second set of bounds=   :   ,   :   ,   :    '
      WRITE(error_str(21:23),'(i3)') bounds_a(1)
      WRITE(error_str(25:27),'(i3)') bounds_a(2)
      WRITE(error_str(29:31),'(i3)') bounds_a(3)
      WRITE(error_str(33:35),'(i3)') bounds_a(4)
      WRITE(error_str(37:39),'(i3)') bounds_a(5)
      WRITE(error_str(41:43),'(i3)') bounds_a(6)
      WRITE(error_str(68:70),'(i3)') bounds_b(1)
      WRITE(error_str(72:74),'(i3)') bounds_b(2)
      WRITE(error_str(76:78),'(i3)') bounds_b(3)
      WRITE(error_str(80:82),'(i3)') bounds_b(4)
      WRITE(error_str(84:86),'(i3)') bounds_b(5)
      WRITE(error_str(88:90),'(i3)') bounds_b(6)
    ELSE
      compare_two_sets_of_bounds = .FALSE.
      error_str = ''
    END IF
  END FUNCTION compare_two_sets_of_bounds

  !> @brief Checks improper combinations of is, ie, js, and je.
  !> @return Returns .false. if there is no error else .true.
  !> @note send_data works in either one or another of two modes.
  ! 1. Input field is a window (e.g. FMS physics)
  ! 2. Input field includes halo data
  ! It cannot handle a window of data that has halos.
  ! (A field with no windows or halos can be thought of as a special case of either mode.)
  ! The logic for indexing is quite different for these two modes, but is not clearly separated.
  ! If both the beggining and ending indices are present, then field is assumed to have halos.
  ! If only beggining indices are present, then field is assumed to be a window.
  !> @par
  ! There are a number of ways a user could mess up this logic, depending on the combination
  ! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
  function check_indices_order(is_in, ie_in, js_in, je_in, error_msg) result(rslt)
    integer, intent(in), optional :: is_in, ie_in, js_in, je_in !< Indices passed to fms_diag_accept_data()
    character(len=*), intent(inout), optional :: error_msg !< An error message used only for testing purpose!!!

    character(len=128) :: err_module_name !< Stores the module name to be used in error calls
    logical :: rslt !< Return value

    rslt = .false. !< If no error occurs.

    err_module_name = 'diag_util_mod:check_indices_order'

    IF ( PRESENT(ie_in) ) THEN
      IF ( .NOT.PRESENT(is_in) ) THEN
        rslt = fms_error_handler(trim(err_module_name), 'ie_in present without is_in', error_msg)
        IF (rslt) return
      END IF
      IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
        rslt = fms_error_handler(trim(err_module_name),&
          & 'is_in and ie_in present, but js_in present without je_in', error_msg)
        IF (rslt) return
      END IF
    END IF

    IF ( PRESENT(je_in) ) THEN
      IF ( .NOT.PRESENT(js_in) ) THEN
        rslt = fms_error_handler(trim(err_module_name), 'je_in present without js_in', error_msg)
        IF (rslt) return
      END IF
      IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
        rslt = fms_error_handler(trim(err_module_name),&
          & 'js_in and je_in present, but is_in present without ie_in', error_msg)
        IF (rslt) return
      END IF
    END IF
  end function check_indices_order

  !> @brief Copies input data to output data with proper type if the input data is present
  !! else sets the output data to a given value val if it is present.
  !! If the value val and the input data are not present, the output data is untouched.
  subroutine real_copy_set(out_data, in_data, val, err_msg)
    real, intent(out) :: out_data !< Proper type copy of in_data
    class(*), intent(in), optional :: in_data !< Data to copy to out_data
    real, intent(in), optional :: val !< Default value to assign to out_data if in_data is absent
    character(len=*), intent(out), optional :: err_msg !< Error message to pass back to caller

    IF ( PRESENT(err_msg) ) err_msg = ''

    IF ( PRESENT(in_data) ) THEN
      SELECT TYPE (in_data)
      TYPE IS (real(kind=r4_kind))
        out_data = in_data
      TYPE IS (real(kind=r8_kind))
        out_data = real(in_data)
      CLASS DEFAULT
        if (fms_error_handler('diag_util_mod:real_copy_set',&
          & 'The in_data is not one of the supported types of real(kind=4) or real(kind=8)', err_msg)) THEN
          return
        end if
      END SELECT
    ELSE
      if (present(val)) out_data = val
    END IF
  end subroutine real_copy_set

  !> @brief Allocates outmask(second argument) with sizes of the first three dimensions of field(first argument).
  !! Initializes the outmask depending on presence/absence of inmask and rmask.
  !! Uses and sets rmask_threshold.
  subroutine init_mask_3d(field, outmask, rmask_threshold, inmask, rmask, err_msg)
    class(*), intent(in) :: field(:,:,:,:)  !< Dummy variable whose sizes only in the first three
                                            !! dimensions are important
    logical, allocatable, intent(inout) :: outmask(:,:,:) !< Output logical mask
    real, intent(inout) :: rmask_threshold !< Holds the values 0.5_r4_kind or 0.5_r8_kind, or related threhold values
                                           !! needed to be passed to the math/buffer update functions.
    logical, intent(in), optional :: inmask(:,:,:) !< Input logical mask
    class(*), intent(in), optional :: rmask(:,:,:) !< Floating point input mask value
    character(len=*), intent(out), optional :: err_msg !< Error message to relay back to caller

    character(len=256) :: err_msg_local !< Stores locally generated error message
    integer :: status !< Stores status of memory allocation call

    ! Initialize character strings
    err_msg_local = ''
    if (present(err_msg)) err_msg = ''

    ! Check if outmask is allocated
    if (allocated(outmask)) deallocate(outmask)
    ALLOCATE(outmask(SIZE(field, 1), SIZE(field, 2), SIZE(field, 3)), STAT=status)
    IF ( status .NE. 0 ) THEN
      WRITE (err_msg_local, FMT='("Unable to allocate outmask(",I5,",",I5,",",I5,"). (STAT: ",I5,")")')&
            & SIZE(field, 1), SIZE(field, 2), SIZE(field, 3), status
      if (fms_error_handler('diag_util_mod:init_mask_3d', trim(err_msg_local), err_msg)) then
        return
      end if
    END IF

    IF ( PRESENT(inmask) ) THEN
      outmask = inmask
    ELSE
      outmask = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) THEN
      SELECT TYPE (rmask)
        TYPE IS (real(kind=r4_kind))
            WHERE (rmask < real(rmask_threshold, kind=r4_kind)) outmask = .FALSE.
            rmask_threshold = real(rmask_threshold, kind=r4_kind)
         TYPE IS (real(kind=r8_kind))
            WHERE ( rmask < real(rmask_threshold, kind=r8_kind) ) outmask = .FALSE.
            rmask_threshold = real(rmask_threshold, kind=r8_kind)
        CLASS DEFAULT
          if (fms_error_handler('diag_util_mod:init_mask_3d',&
            & 'The rmask is not one of the supported types of real(kind=4) or real(kind=8)', err_msg)) then
          end if
      END SELECT
    END IF
  end subroutine init_mask_3d
#endif
end module fms_diag_util_mod
