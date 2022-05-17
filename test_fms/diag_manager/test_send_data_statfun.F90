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

!> @brief  This programs tests
program test_send_data_statfun
   use mpp_mod, only: mpp_init, mpp_set_stack_size, mpp_init_test_requests_allocated
   use mpp_io_mod, only: mpp_io_init
   use fms_mod, ONLY: error_mesg, FATAL,NOTE

   USE time_manager_mod, ONLY: time_type

   use diag_data_mod, only:  input_fields, output_fields
   use fms_diag_weight_procs_mod !! only:  FmsWeightProcCfg_t
   use fms_send_data_statfun_mod

   implicit  none
   integer :: sum                                                  !< Temp sum of vaalues of id sets
   !!
   integer :: ierr                                                 !< An error flag
   !!
   logical :: test_passed                                          !< Flag indicating if the test_passed
   !! These fields below used to initialize diag object data. TBD
   integer, dimension(2) :: axes
   !!type (diag_fields_type)  :: diag_field
   logical :: temp_result
   REAL, DIMENSION(10,10,10) :: field
   TYPE(statfun_idx_cfg_t) :: idx_cfg

   INTEGER:: diag_field_id, out_num
   INTEGER:: sample !!diurnal_index
   REAL :: weight
   TYPE (time_type):: time
   INTEGER:: is_in, js_in, ks_in,ie_in,je_in, ke_in
   LOGICAL, DIMENSION(10,10,10) :: mask
   REAL, DIMENSION(10,10,10) :: rmask
   CHARACTER(len=128) :: err_msg, err_msg_local
   integer, dimension(3) ::l_start, l_end
   REAL :: missvalue
   LOGICAL :: missvalue_present

   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask

   call mpp_init(mpp_init_test_requests_allocated)
   call mpp_io_init()
   call mpp_set_stack_size(145746)

   call error_mesg('test_send_data_statfun', 'Test has started',NOTE)
   diag_field_id = 1
   out_num = 1
   sample = 1

   allocate(input_fields(1))
   allocate(output_fields(1))

   call init_input_output_fields_cfg(diag_field_id, out_num)

   call init_field_values (field, 10,10,10)

   test_passed = .true.  !! will be set to false if there are any issues.

   idx_cfg%is = 3
   idx_cfg%js = 3
   idx_cfg%ks = 3
   idx_cfg%ie = 8
   idx_cfg%je = 8
   idx_cfg%ke = 8
   idx_cfg%hi = 1
   idx_cfg%hj = 2
   idx_cfg%f1 = 1
   idx_cfg%f2 = 1
   idx_cfg%f3 = 1
   idx_cfg%f4 = 1
   l_start(1) = 3
   l_start(2) = 2
   l_start(3) = 3
   l_end(1)  = 8
   l_end(2) = 8
   l_end(3) = 8

   missvalue = 1.0e-5
   missvalue_present = .true.

   temp_result = average_the_field(diag_field_id, field, out_num, mask, weight, &
   & sample, missvalue, missvalue_present, l_start, l_end, idx_cfg, err_msg, err_msg_local )
   IF (temp_result .eqv. .FALSE.) THEN
      DEALLOCATE(oor_mask)
      RETURN
   endif

   call error_mesg('test_send_data_statfun', 'Test has finished',NOTE)

   call MPI_finalize(ierr)

CONTAINS

   subroutine init_input_output_fields_cfg( field_id, out_num )
      integer, intent(in):: field_id
      integer, intent(in) :: out_num
      integer :: sample
      sample = 1

      !!values of input_fields used by the tests
      input_fields(diag_field_id)%field_name = 'X'
      input_fields(diag_field_id)%mask_variant = .false.
      input_fields(diag_field_id)%module_name = 'test'
      input_fields(diag_field_id)%issued_mask_ignore_warning = .true.

      !!values of input_fields used by the tests
      output_fields(out_num)%pow_value = 2

      output_fields(out_num)%output_name = 'X'
      output_fields(out_num)%phys_window = .false.
      output_fields(out_num)%need_compute = .false.
      output_fields(out_num)%reduced_k_range = .false.
      allocate(output_fields(out_num)%buffer(10,10,10,1))
      allocate(output_fields(out_num)%counter(10,10,10,1))
      allocate(output_fields(out_num)%num_elements(1))
      allocate(output_fields(out_num)%count_0d(1))
   end subroutine init_input_output_fields_cfg

   subroutine init_field_values (field, nx, ny, nz)
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: field
      INTEGER, INTENT(in) :: NX,NY,NZ

      INTEGER i,j,k

      DO i = 1, NX
         DO j = 1, NY
            DO k = 1, NZ
               field(i,j,k) = 1.1
            END DO
         END DO
      END DO
   end subroutine init_field_values

   subroutine check_field_value(field, buffer, counter, NX, NY, NZ)
      REAL, DIMENSION(:,:,:), INTENT(IN) :: field
      INTEGER, INTENT(in) :: NX,NY,NZ
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: buffer
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: counter
   end subroutine check_field_value




end program test_send_data_statfun


