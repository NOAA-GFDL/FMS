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

   use diag_data_mod, only:  input_fields, output_fields
   use fms_send_data_statfun_mod, ONLY: fms_diag_field_procs_t

   implicit  none
   integer :: sum                                                  !< Temp sum of vaalues of id sets
   !!
   integer :: ierr                                                 !< An error flag
   !!
   logical :: test_passed                                          !< Flag indicating if the test_passed
   !! These fields below used to initialize diag object data. TBD
   logical :: temp_result
   REAL, DIMENSION(10,10,10) :: field

   INTEGER:: diag_field_id, out_num
   INTEGER:: sample !!diurnal_index
   REAL :: weight
   INTEGER:: hi, hj  !!for halo sizes
   LOGICAL, DIMENSION(10,10,10) :: mask
   REAL, DIMENSION(10,10,10) :: rmask
   CHARACTER(len=128) :: err_msg, err_msg_local
   integer, dimension(3) ::l_start, l_end
   REAL :: missvalue
   LOGICAL :: missvalue_present

   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
   TYPE(fms_diag_field_procs_t), allocatable :: sprocs_obj


   call mpp_init(mpp_init_test_requests_allocated)
   call mpp_io_init()
   call mpp_set_stack_size(145746)

   call error_mesg('test_send_data_statfun', 'Test has started',NOTE)
   diag_field_id = 1
   out_num = 1
   sample = 1
   weight = 1.0

   allocate(input_fields(1))
   allocate(output_fields(1))

   call init_input_output_fields_cfg(diag_field_id, out_num,1)

   call init_field_values (field, 10,10,10)

   call init_output_field_values(out_num)

   test_passed = .true.  !! will be set to false if there are any issues.


   hi = 0 !!halo size i
   hj = 0 !!halo size j
   l_start(1) = 1 !!local (to PE) start inddex
   l_start(2) = 1
   l_start(3) = 1
   l_end(1)  = 10
   l_end(2) = 10
   l_end(3) = 10


   allocate(sprocs_obj)

   call sprocs_obj%initialize(1+hi, 1+hj, 1, 10 - hi, 10 - hj, 10,&
   &  hi, hj, 1 + hi, 10 - hi, 1 + hj, 10 - hj, &
   &  output_fields(out_num)%pow_value, output_fields(out_num)%phys_window, &
   &  output_fields(out_num)%need_compute,output_fields(out_num)%reduced_k_range, &
   &  output_fields(out_num)%time_rms,  output_fields(out_num)%time_max, &
   &  output_fields(out_num)%time_min,  output_fields(out_num)%time_sum)


   missvalue = 1.0e-5

   !! Case: mask_var=false & missval not present & mask not present & not_reduced_k_range
   missvalue_present = .false.
   temp_result = sprocs_obj%average_the_field(diag_field_id, field, out_num, &
      & output_fields(out_num)%buffer, output_fields(out_num)%counter, &
      & output_fields(out_num)%ntval,  output_fields(out_num)%output_name, &
      & mask, weight, sample, missvalue, missvalue_present, &
      l_start, l_end, err_msg, err_msg_local )
   call check_results_1(output_fields(out_num)%buffer(:,:,:,sample))

   missvalue_present = .true.

   ! temp_result = average_the_field(diag_field_id, field, out_num, mask, weight, &
   !& sample, missvalue, missvalue_present, l_start, l_end, idx_cfg, err_msg, err_msg_local )
   !IF (temp_result .eqv. .FALSE.) THEN
   !   DEALLOCATE(oor_mask)
   !endif

   call print_output_field_values(1)

   call error_mesg('test_send_data_statfun', 'Test has finished',NOTE)

   call MPI_finalize(ierr)

CONTAINS

   subroutine init_input_output_fields_cfg( field_id, out_num, pow_val )
      integer, intent(in):: field_id
      integer, intent(in) :: out_num
      integer, intent(in) ::pow_val
      integer :: sample
      sample = 1

      !!values of input_fields used by the tests
      input_fields(diag_field_id)%field_name = 'X'
      input_fields(diag_field_id)%mask_variant = .false.
      input_fields(diag_field_id)%module_name = 'test'
      input_fields(diag_field_id)%issued_mask_ignore_warning = .true.

      !!values of input_fields used by the tests
      output_fields(out_num)%pow_value = pow_val

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
               field(i,j,k) = one_dim_from_three(i,j,k,NX,NY,NZ)
            END DO
         END DO
      END DO
   end subroutine init_field_values

   subroutine init_output_field_values (onum)
      INTEGER, INTENT(IN) :: onum
      output_fields(onum)%buffer = 0
      output_fields(onum)%counter = 0
      output_fields(onum)%count_0d = 0
   end subroutine init_output_field_values

   subroutine check_field_value(field, buffer, counter, NX, NY, NZ)
      REAL, DIMENSION(:,:,:), INTENT(IN) :: field
      INTEGER, INTENT(in) :: NX,NY,NZ
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: buffer
      REAL, DIMENSION(:,:,:,:), INTENT(IN) :: counter
   end subroutine check_field_value

   subroutine print_output_field_values (onum)
      INTEGER, INTENT(IN) :: onum
      INTEGER i,j,k
      DO i = 1,10
         DO j=i,10
            print "(10f10.1)",output_fields(onum)%buffer(i,j,:,1)
         end do
      end do
   end subroutine print_output_field_values

   subroutine check_results_1(buff)
      REAL, DIMENSION(:,:,:), INTENT(IN) :: buff
      INTEGER :: NX,NY,NZ
      INTEGER :: i,j,k
      LOGICAL :: pass

      pass = .true.
      NX = size(buff,1)
      NY= size(buff,2)
      NZ= size(buff,3)
      DO i = 1, NX
         DO j = 1, NY
            DO k = 1, NZ
               if ( one_dim_from_three(i,j,k,NX,NY,NZ)  /= buff(i,j,k) ) then
                  pass = .false.
               end if
            END DO
         END DO
      END DO
      if ( pass .eqv. .false.) then
         call error_mesg('check_results_1', 'Test has failed',FATAL)
      end if
   end subroutine check_results_1

   pure integer function one_dim_from_three(i,j,k,NX,NY,NZ)
      INTEGER, INTENT(IN) :: i, j, k
      INTEGER, INTENT(IN) :: NX, NY, NZ
      one_dim_from_three =  (k-1) * NX * NY + (j-1) * NX + i
   end function one_dim_from_three

end program test_send_data_statfun


