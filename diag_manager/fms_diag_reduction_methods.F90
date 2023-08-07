!> \author Ganga Purja Pun
!> \email GFDL.Climate.Model.Info@noaa.gov
!! \brief Contains routines for the modern diag manager
!! These routines are meant to be used for checks and in reduction methods.

module fms_diag_reduction_methods_mod
  implicit none
  private

  public :: check_indices_order

  contains

  !> @brief Checks improper combinations of is, ie, js, and je.
  !> @return The error message, empty string if no errors were found
  !> @note accept_data works in either one or another of two modes.
  !! 1. Input field is a window (e.g. FMS physics)
  !! 2. Input field includes halo data
  !! It cannot handle a window of data that has halos.
  !! (A field with no windows or halos can be thought of as a special case of either mode.)
  !! The logic for indexing is quite different for these two modes, but is not clearly separated.
  !! If both the beggining and ending indices are present, then field is assumed to have halos.
  !! If only beggining indices are present, then field is assumed to be a window.
  !> @par
  !! There are a number of ways a user could mess up this logic, depending on the combination
  !! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
  function check_indices_order(is_in, ie_in, js_in, je_in) result(error_msg)
    integer, intent(in), optional :: is_in, ie_in, js_in, je_in !< Indices passed to fms_diag_accept_data()
    character(len=128) :: error_msg !< An error message used only for testing purpose!!!

    error_msg = ""
    IF ( PRESENT(ie_in) ) THEN
      IF ( .NOT.PRESENT(is_in) ) THEN
        error_msg = 'ie_in present without is_in'
        return
      END IF
      IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
        error_msg = 'is_in and ie_in present, but js_in present without je_in'
        return
      END IF
    END IF

    IF ( PRESENT(je_in) ) THEN
      IF ( .NOT.PRESENT(js_in) ) THEN
        error_msg = 'je_in present without js_in'
        return
      END IF
      IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
        error_msg = 'js_in and je_in present, but is_in present without ie_in'
        return
      END IF
    END IF
  end function check_indices_order

end module fms_diag_reduction_methods_mod