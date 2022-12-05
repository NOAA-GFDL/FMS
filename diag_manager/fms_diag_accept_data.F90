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
!> @defgroup diag_data_mod diag_data_mod
!> @ingroup diag_manager
!! @brief Type descriptions and global variables for the diag_manager modules.
!! @author Tom Robinson
module fms_diag_accept_data_mod

implicit none
private
contains
!> Used to accept data from the send_data functions.  If this is in an openmp region with more than 
!! one thread, the data is buffered in the field object and processed later.  If onle a single thread
!! is being used, then the processing can be done and stored in the buffer object.  The hope is that
!! the increase in memory footprint related to buffering can be handled by the shared memory of the
!! multithreaded case.
!! \note If some of the diag manager is offloaded in the future, then is should be treated similarly
!! to the multi-threaded option for processing later
logical function fms_diag_accept_data (diag_field_id, field, time, is_in, js_in, ks_in, &
                  mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
  INTEGER, INTENT(in) :: diag_field_id !< The ID of the input diagnostic field
  CLASS(*), DIMENSION(:,:,:), INTENT(in) :: field !< The data for the input diagnostic
  REAL, INTENT(in), OPTIONAL :: weight !< The weight used for averaging
  TYPE (time_type), INTENT(in), OPTIONAL :: time !< The current time
  INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in !< Indicies for the variable
  LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask !< The location of the mask
  CLASS(*), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: rmask !< The masking values
  CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg !< An error message returned

  integer :: omp_num_threads !< Number of openmp threads
  integer :: omp_level !< The openmp active level

!> initialize the number of threads and level to be 0
  omp_num_threads = 0
  omp_level = 0
#if defined(_OPENMP)
  omp_num_threads = omp_get_num_threads()
  omp_level = omp_get_level()
#endif
  main_if: if (omp_num_threads > 1 .AND. omp_level > 0) then !If this is true, only buffer data
    ! Make sure buffer is clear
    ! allocate buffer
    ! add to buffer
    fms_diag_accept_data = .TRUE.
    return
  else
    ! Loop through fields and do averages
    fms_diag_accept_data = .TRUE.
    return 
  end if main_if
!> Return false if nothing is done
  fms_diag_accept_data = .FALSE.
  return 
end function fms_diag_accept_data
end module fms_diag_accept_data_mod
