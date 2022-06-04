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

!> @defgroup fms_send_data_statfun_mod fms_send_data_statfun_mod
!> @ingroup diag_manager
!> @brief fms_send_data_statfun_mod contains functions used by the
!! diag_managers send_data routines tto prepare the filed/data for output.
!!  The functions include agregation, min and max, etc.
!!
!! <TT>fms_send_data_stattfun_mod</TT>  defines functions used by the
!! diag_managers send_data routines tto prepare the filed/data for output.
!!  The functions include agregation, min and max, etc.
!!
!> @file
!> @brief File for @ref fms_send_data_statfun_mod
!> @addtogroup fms_diag_dlinked_list_mod
!> @{
MODULE fms_send_data_statfun_mod

   use platform_mod
   USE mpp_mod, ONLY: mpp_pe, mpp_root_pe

   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,&
   & fms_error_handler
   USE diag_data_mod, ONLY:  input_fields, output_fields, debug_diag_manager
   use diag_util_mod, ONLY: check_out_of_bounds, update_bounds

   IMPLICIT NONE

   TYPE, public :: fms_diag_field_procs_t
      private
      INTEGER :: f1,f2,f3,f4
      INTEGER :: is, js, ks, ie, je, ke
      INTEGER :: hi !< halo size in x direction
      INTEGER :: hj !< halo size in y direction
      INTEGER :: pow_value   ! Power value for rms or pow(x) calculations
      LOGICAL :: phys_window   !!TODO: name change PHYS_WINDOWS -->> ? OMP subsetted data, See output_fields
      LOGICAL :: need_compute
      LOGICAL :: reduced_k_range
      LOGICAL :: time_rms, time_max, time_min, time_sum
   contains
      procedure :: initialize => initialize_imp
      procedure :: average_the_field => average_the_field_imp
      procedure :: sample_the_field => sample_the_field_imp
   end type fms_diag_field_procs_t

CONTAINS

   SUBROUTINE initialize_imp (this, is, js , ks, ie, je, ke, hi, hj, f1, f2, f3, f4, &
       & pow_value,  phys_window, need_compute, reduced_k_range, &
       & time_rms, time_max, time_min, time_sum )
    CLASS(fms_diag_field_procs_t), INTENT(inout)  :: this
    INTEGER :: is, js, ks, ie, je, ke
    INTEGER :: hi, hj
    INTEGER :: f1, f2, f3, f4
    INTEGER :: pow_value
    LOGICAL :: phys_window , need_compute , reduced_k_range
    LOGICAL :: time_rms, time_max, time_min, time_sum

    this%is = is
    this%js = js
    this%ks = ks
    this%ie = ie
    this%je = je
    this%ke = ke

    this%hi = hi
    this%hj = hj

    this%f1 = f1
    this%f2 = f2
    this%f3 = f3
    this%f4 = f4

    this%pow_value = pow_value
    this%phys_window = phys_window
    this%need_compute = need_compute
    this%reduced_k_range =reduced_k_range

    !Is this output field the rms? If so, then average is also .TRUE.
    this%time_rms = time_rms
    this%time_max = time_max
    this%time_min = time_min
    ! Sum output over time interval
    this%time_sum = time_sum

    END SUBROUTINE initialize_imp


   FUNCTION AVERAGE_THE_FIELD_IMP(this, diag_field_id, field, out_num, ofb, ofc, &
   & mask, weight1, sample, missvalue, missvalue_present, &
   & l_start, l_end, err_msg,  err_msg_local ) result( succeded )
      CLASS(fms_diag_field_procs_t) , INTENT(inout) :: this
      INTEGER, INTENT(in) :: diag_field_id
      REAL, DIMENSION(:,:,:), INTENT(in) :: field
      INTEGER, INTENT(in) :: out_num
      REAL, allocatable, DIMENSION(:,:,:,:), INTENT(inout) :: ofb
      !class(*),  pointer, INTENT(inout) :: ofb_in  !!TODO:
      REAL, allocatable, DIMENSION(:,:,:,:), INTENT(inout) :: ofc
      LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
      REAL, INTENT(in) :: weight1
      INTEGER, INTENT(in) :: sample
      REAL, INTENT(in), OPTIONAL :: missvalue
      LOGICAL, INTENT(in) :: missvalue_present
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
      CHARACTER(len=*), INTENT(inout), OPTIONAL :: err_msg
      CHARACTER(len=*), INTENT(inout) :: err_msg_local

      LOGICAL :: succeded

      !!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
      CHARACTER(len=128):: error_string

      !! TODO: if possible
      !!TYPE(FmsWeightProcCfg_t), allocatable :: weight_procs

       ! Power value for rms or pow(x) calculations
      INTEGER :: pow_value, is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
      LOGICAL :: phys_window , need_compute , reduced_k_range

      INTEGER :: ksr, ker
      INTEGER :: i, j, k,  i1, j1, k1

      INTEGER :: numthreads
      INTEGER :: active_omp_level
#if defined(_OPENMP)
      INTEGER :: omp_get_num_threads !< OMP function
      INTEGER :: omp_get_level !< OMP function
#endif


      !REAL, DIMENSION(:,:,:,:), pointer :: ofb
      !select type (ofb_in)
      ! type is ( buff_r_4d_t)
      !   ofb => ofb_in%buffer
      ! class default
      !   stop 'Error in type selection'
      !end select


      ksr= l_start(3)
      ker= l_end(3)
      is = this%is
      js = this%js
      ks = this%ks
      ie = this%ie
      je = this%je
      ke = this%ke
      hi = this%hi
      hj = this%hj
      f1 = this%f1
      f2 = this%f2
      f3 = this%f3
      f4 = this%f4
      pow_value = this%pow_value
      phys_window = this%phys_window
      reduced_k_range = this%reduced_k_range
      need_compute = this%need_compute

!$OMP CRITICAL
      input_fields(diag_field_id)%numthreads = 1
      active_omp_level=0
#if defined(_OPENMP)
      input_fields(diag_field_id)%numthreads = omp_get_num_threads()
      input_fields(diag_field_id)%active_omp_level = omp_get_level()
#endif
      numthreads = input_fields(diag_field_id)%numthreads
      active_omp_level = input_fields(diag_field_id)%active_omp_level
!$OMP END CRITICAL

      MASK_VAR_IF: IF ( input_fields(diag_field_id)%mask_variant ) THEN
         IF ( need_compute ) THEN
            WRITE (error_string,'(a,"/",a)')  &
            & TRIM(input_fields(diag_field_id)%module_name), &
            & TRIM(output_fields(out_num)%output_name)
            IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
            & ', regional output NOT supported with mask_variant', err_msg)) THEN
               !!DEALLOCATE(oor_mask)
               succeded = .FALSE.
               RETURN
            END IF
         END IF

         ! Should reduced_k_range data be supported with the mask_variant option   ?????
         ! If not, error message should be produced and the reduced_k_range loop below eliminated
         MASK_PR_1_IF: IF ( PRESENT(mask) ) THEN
            MISSVAL_PR_1_IF: IF ( missvalue_present ) THEN !!(section: mask_varian .eq. true + mask present)
               IF ( debug_diag_manager ) THEN
                  CALL update_bounds(out_num, is-hi, ie-hi, this%js-hj, this%je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               !!
               IF( numthreads>1 .AND. phys_window ) then
                  REDU_KR1_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                 ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) +&
                                 & (field(i-is+1+hi, j-js+1+hj, k) * weight1 ) ** pow_value
                                 ofc(i-hi,j-hj,k1,sample) = ofc(i-hi,j-hj,k1,sample) + weight1
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR1_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) + &
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                                 ofc(i-hi,j-hj,k,sample) = ofc(i-hi,j-hj,k,sample) + weight1
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR1_IF
               ELSE
!$OMP CRITICAL
                  REDU_KR2_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                 ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) + &
                                 & (field(i-is+1+hi, j-js+1+hj, k) * weight1 ) ** pow_value
                                 ofc(i-hi,j-hj,k1,sample) = ofc(i-hi,j-hj,k1,sample) + weight1
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR2_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN !!USE WHERE
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) +  &
                                 & ( field(i-is+1+hi,j-js+1+hj,k) * weight1 ) ** pow_value
                                 ofc(i-hi,j-hj,k,sample) = ofc(i-hi,j-hj,k,sample) + weight1
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR2_IF
!$OMP END CRITICAL
               END IF
            ELSE MISSVAL_PR_1_IF
               WRITE (error_string,'(a,"/",a)')&
               & TRIM(input_fields(diag_field_id)%module_name), &
               & TRIM(output_fields(out_num)%output_name)
               IF(fms_error_handler('diag_manager_mod::send_data_3d', &
               & 'module/output_field '//TRIM(error_string)//', variable mask but no missing value defined', &
               & err_msg)) THEN
                  succeded = .FALSE.
                  RETURN
               END IF
            END IF  MISSVAL_PR_1_IF
         ELSE MASK_PR_1_IF ! no mask present
            WRITE (error_string,'(a,"/",a)')&
            & TRIM(input_fields(diag_field_id)%module_name), &
            & TRIM(output_fields(out_num)%output_name)
            IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
            & ', variable mask but no mask given', err_msg)) THEN
               succeded = .FALSE.
               RETURN
            END IF
         END IF MASK_PR_1_IF
      ELSE MASK_VAR_IF
         MASK_PR_2_IF: IF ( PRESENT(mask) ) THEN
            MISSVAL_PR_2_IF: IF ( missvalue_present ) THEN !!section:(mask_var false +mask present +missval prsnt)
               NDCMP_RKR_1_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                    ofb(i1,j1,k1,sample) = ofb(i1,j1,k1,sample) +&
                                    & ( field(i-is+1+hi,j-js+1+hj,k) * weight1 ) ** pow_value
                                 ELSE
                                    ofb(i1,j1,k1,sample) = missvalue
                                 END IF
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                    ofb(i1,j1,k1,sample) = ofb(i1,j1,k1,sample) + &
                                    & ( field(i-is+1+hi,j-js+1+hj,k) * weight1 ) ** pow_value
                                 ELSE
                                    ofb(i1,j1,k1,sample) = missvalue
                                 END IF
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  ENDIF
!$OMP CRITICAL
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                           output_fields(out_num)%num_elements(sample) = &
                           & output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                        END IF
                     END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_1_IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                 ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) + &
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k1,sample)= missvalue
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                 ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) + &
                                 & ( field(i-is+1+hi,j-js+1+hj,k) * weight1 ) **  pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k1,sample)= missvalue
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_1_IF
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) + &
                                 & ( field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k,sample)= missvalue
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) + &
                                 & ( field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k,sample)= missvalue
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_1_IF
!$OMP CRITICAL
               IF ( need_compute .AND. .NOT.phys_window ) THEN
                  IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3))) ) &
                  & output_fields(out_num)%count_0d(sample) =&
                  & output_fields(out_num)%count_0d(sample) + weight1
               ELSE
                  IF ( ANY(mask(f1:f2,f3:f4,ks:ke)) ) output_fields(out_num)%count_0d(sample) =&
                  & output_fields(out_num)%count_0d(sample)+weight1
               END IF
!$OMP END CRITICAL
            ELSE MISSVAL_PR_2_IF !! (section: mask_varian .eq. false + mask present + miss value not present)
               IF (   (.NOT.ALL(mask(f1:f2,f3:f4,ks:ke)) .AND. mpp_pe() .EQ. mpp_root_pe()).AND.&
               &  .NOT.input_fields(diag_field_id)%issued_mask_ignore_warning ) THEN
                  ! <ERROR STATUS="WARNING">
                  !   Mask will be ignored since missing values were not specified for field <field_name>
                  !   in module <module_name>
                  ! </ERROR>
                  CALL error_mesg('diag_manager_mod::send_data_3d',&
                  & 'Mask will be ignored since missing values were not specified for field '//&
                  & trim(input_fields(diag_field_id)%field_name)//' in module '//&
                  & trim(input_fields(diag_field_id)%module_name), WARNING)
                  input_fields(diag_field_id)%issued_mask_ignore_warning = .TRUE.
               END IF
               NDCMP_RKR_2_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,sample)=  ofb(i1,j1,:,sample)+ &
                              & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)) * weight1 ) ** pow_value
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,sample) = ofb(i1,j1,:,sample) + &
                              & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)) * weight1 ) ** pow_value
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                           output_fields(out_num)%num_elements(sample)=&
                           & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                        END IF
                     END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_2_IF
                  IF (numthreads>1 .AND. phys_window) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) = ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                     & (field(f1:f2,f3:f4,ksr:ker) * weight1) ** pow_value
                  ELSE
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) = ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                     & (field(f1:f2,f3:f4,ksr:ker) * weight1) ** pow_value
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_2_IF
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '') THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                     & (field(f1:f2,f3:f4,ks:ke) * weight1) ** pow_value
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                     & (field(f1:f2,f3:f4,ks:ke) * weight1) ** pow_value
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_2_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
               & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_2_IF
         ELSE MASK_PR_2_IF !!(section: mask_variant .eq. false + mask not present + missvalue)
            MISSVAL_PR_3_IF: IF (missvalue_present ) THEN
               NDCMP_RKR_3_IF: IF ( need_compute ) THEN
                  NTAPW_IF: If( numthreads>1 .AND. phys_window ) then
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                    ofb(i1,j1,k1,sample) = ofb(i1,j1,k1,sample) + &
                                    & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                                 ELSE
                                    ofb(i1,j1,k1,sample) = missvalue
                                 END IF
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE NTAPW_IF
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                    ofb(i1,j1,k1,sample) = ofb(i1,j1,k1,sample) + &
                                    & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                                 ELSE
                                    ofb(i1,j1,k1,sample) = missvalue
                                 END IF
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF NTAPW_IF
!$OMP CRITICAL
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj) THEN
                           output_fields(out_num)%num_elements(sample) =&
                           & output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                        END IF
                     END DO
                  END DO
                  IF ( .NOT.phys_window ) THEN
                     DO k = l_start(3), l_end(3)
                        DO j=l_start(2)+hj, l_end(2)+hj
                           DO i=l_start(1)+hi, l_end(1)+hi
                              IF ( field(i,j,k) /= missvalue ) THEN
                                 output_fields(out_num)%count_0d(sample) = &
                                 & output_fields(out_num)%count_0d(sample) + weight1
                                 EXIT
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_3_IF
                  if( numthreads>1 .AND. phys_window ) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                 ofb(i-hi,j-hj,k1,sample) =  ofb(i-hi,j-hj,k1,sample) + &
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k1,sample) = missvalue
                              END IF
                           END DO
                        END DO
                     END DO
                  else
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                 ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) +&
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k1,sample) = missvalue
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO k = ksr, ker
                     k1=k-ksr+1
                     DO j=f3, f4
                        DO i=f1, f2
                           !! TODO: verify this below
                           IF ( field(i,j,k) /= missvalue ) THEN
                              output_fields(out_num)%count_0d(sample) = &
                              & output_fields(out_num)%count_0d(sample) + weight1
                              EXIT
                           END IF
                        END DO
                     END DO
                  END DO
!$OMP END CRITICAL
               ELSE NDCMP_RKR_3_IF
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) +&
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k,sample) = missvalue
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                 ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) +&
                                 & (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                              ELSE
                                 ofb(i-hi,j-hj,k,sample) = missvalue
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO k=ks, ke
                     DO j=f3, f4
                        DO i=f1, f2
                           IF ( field(i,j,k) /= missvalue ) THEN
                              output_fields(out_num)%count_0d(sample) = &
                              & output_fields(out_num)%count_0d(sample)  + weight1
                              EXIT
                           END IF
                        END DO
                     END DO
                  END DO
!$OMP END CRITICAL
               END IF NDCMP_RKR_3_IF
            ELSE MISSVAL_PR_3_IF !!(section: mask_variant .eq. false + mask not present + missvalue not present)
               NDCMP_RKR_4_IF: IF ( need_compute ) THEN
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,sample) = ofb(i1,j1,:,sample) + &
                              & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)) * weight1 ) ** pow_value
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,sample)= ofb(i1,j1,:,sample) +&
                              & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)) * weight1 ) ** pow_value
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                           output_fields(out_num)%num_elements(sample) =&
                           & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                        END IF
                     END DO
                  END DO
!$OMP END CRITICAL
                  ! Accumulate time average
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_4_IF
                  ksr= l_start(3)
                  ker= l_end(3)
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                     & (field(f1:f2,f3:f4,ksr:ker) * weight1) ** pow_value
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                     & (field(f1:f2,f3:f4,ksr:ker) * weight1) ** pow_value
!$OMP END CRITICAL
                  END IF

               ELSE NDCMP_RKR_4_IF
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF (fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                     & (field(f1:f2,f3:f4,ks:ke) * weight1) ** pow_value
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                     & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                     & (field(f1:f2,f3:f4,ks:ke) * weight1) ** pow_value
                     !!
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_4_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
               & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_3_IF
         END IF MASK_PR_2_IF ! if mask present
      END IF MASK_VAR_IF

!$OMP CRITICAL
      IF ( .NOT.need_compute .AND. .NOT.reduced_k_range )&
      & output_fields(out_num)%num_elements(sample) =&
      & output_fields(out_num)%num_elements(sample) + (ie-is+1)*(je-js+1)*(ke-ks+1)
      IF ( reduced_k_range ) &
      & output_fields(out_num)%num_elements(sample) = output_fields(out_num)%num_elements(sample) +&
      & (ie-is+1)*(je-js+1)*(ker-ksr+1)
!$OMP END CRITICAL

      succeded = .TRUE.
      RETURN

   END FUNCTION AVERAGE_THE_FIELD_IMP

   FUNCTION SAMPLE_THE_FIELD_IMP (this, diag_field_id, field, out_num, &
   & mask, sample, missvalue, missvalue_present, &
   & l_start, l_end, err_msg,  err_msg_local) result( succeded )
      CLASS(fms_diag_field_procs_t), INTENT(inout)  :: this
      INTEGER, INTENT(in) :: diag_field_id
      REAL, DIMENSION(:,:,:), INTENT(in) :: field
      INTEGER, INTENT(in) :: out_num
      LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
      INTEGER, INTENT(in) :: sample
      REAL, INTENT(in), OPTIONAL :: missvalue   !!TODO: optional?
      LOGICAL, INTENT(in) :: missvalue_present
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
      CHARACTER(len=*), INTENT(inout), OPTIONAL :: err_msg
      CHARACTER(len=256), INTENT(inout) :: err_msg_local
      LOGICAL :: succeded

      !!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
      CHARACTER(len=128):: error_string

      ! Power value for rms or pow(x) calculations
      INTEGER :: pow_value, is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
      LOGICAL :: phys_window , need_compute , reduced_k_range


      INTEGER :: ksr, ker
      INTEGER :: i, j, k, i1, j1, k1

      LOGICAL :: time_rms, time_max, time_min, time_sum
#if defined(_OPENMP)
      INTEGER :: omp_get_num_threads !< OMP function
      INTEGER :: omp_get_level !< OMP function
#endif

      ksr= l_start(3)
      ker= l_end(3)

      is = this%is
      js = this%js
      ks = this%ks
      ie = this%ie
      je = this%je
      ke = this%ke
      hi = this%hi
      hj = this%hj
      f1 = this%f1
      f2 = this%f2
      f3 = this%f3
      f4 = this%f4

      time_rms = this%time_rms
      time_max = this%time_max
      time_min = this%time_min
      time_sum = this%time_sum

      reduced_k_range = this%reduced_k_range
      need_compute = this%need_compute

      ASSOCIATE( OFB => output_fields(out_num)%buffer)

         ! Add processing for Max and Min
         TIME_IF: IF ( time_max ) THEN
            MASK_PRSNT_1_IF: IF ( PRESENT(mask) ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                              & field(i-is+1+hi,j-js+1+hj,k)>OFB(i1,j1,k1,sample)) THEN
                                 OFB(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND. &
                  & field(f1:f2,f3:f4,ksr:ker) > OFB(is-hi:ie-hi,js-hj:je-hj,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                  & field(f1:f2,f3:f4,ks:ke)>OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
               END IF
            ELSE MASK_PRSNT_1_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF(l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              IF ( field(i-is+1+hi,j-js+1+hj,k) > OFB(i1,j1,k1,sample) ) THEN
                                 OFB(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker) > OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke) > OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
               END IF
            END IF MASK_PRSNT_1_IF
            output_fields(out_num)%count_0d(sample) = 1
            !END TIME MAX
         ELSE IF ( time_min ) THEN TiME_IF
            MASK_PRSNT_2_IF: IF ( PRESENT(mask) ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                              & field(i-is+1+hi,j-js+1+hj,k) < OFB(i1,j1,k1,sample) ) THEN
                                 OFB(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND.&
                  & field(f1:f2,f3:f4,ksr:ker) < OFB(is-hi:ie-hi,js-hj:je-hj,:,sample)) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                  & field(f1:f2,f3:f4,ks:ke) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
               END IF
            ELSE MASK_PRSNT_2_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              IF ( field(i-is+1+hi,j-js+1+hj,k) < OFB(i1,j1,k1,sample) ) THEN
                                 OFB(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker) < OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
               END IF
            END IF MASK_PRSNT_2_IF
            output_fields(out_num)%count_0d(sample) = 1

            !! END_TIME_MIN
         ELSE IF ( time_sum ) THEN TIME_IF
            MASK_PRSNT_3_IF: IF ( PRESENT(mask) ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                 OFB(i1,j1,k1,sample) = &
                                    OFB(i1,j1,k1,sample) + &
                                    field(i-is+1+hi,j-js+1+hj,k)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = &
                  &   OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                  &   field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke) ) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = &
                  &  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) + &
                  &  field(f1:f2,f3:f4,ks:ke)
               END IF
            ELSE MASK_PRSNT_3_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              OFB(i1,j1,k1,sample) = &
                              &    OFB(i1,j1,k1,sample) + &
                              &    field(i-is+1+hi,j-js+1+hj,k)
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = &
                  &  OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                  &  field(f1:f2,f3:f4,ksr:ker)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = &
                  &    OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) + &
                  &    field(f1:f2,f3:f4,ks:ke)
               END IF
            END IF MASK_PRSNT_3_IF
            output_fields(out_num)%count_0d(sample) = 1
            !END time_sum
         ELSE TIME_IF !! ( not average, not min, not max, not sum )
            output_fields(out_num)%count_0d(sample) = 1
            IF ( need_compute ) THEN
               DO j = js, je
                  DO i = is, ie
                     IF (l_start(1)+hi<= i .AND. i<= l_end(1)+hi .AND. l_start(2)+hj<= j .AND. j<= l_end(2)+hj) THEN
                        i1 = i-l_start(1)-hi+1
                        j1 = j-l_start(2)-hj+1
                        OFB(i1,j1,:,sample) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))
                     END IF
                  END DO
               END DO
               ! instantaneous output
            ELSE IF ( reduced_k_range ) THEN
               ksr = l_start(3)
               ker = l_end(3)
               OFB(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
            ELSE
               IF ( debug_diag_manager ) THEN
                  !!TODO update_bounds and chck_out_of_bounds may need mods for new diag
                  CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
            END IF

            IF ( PRESENT(mask) .AND. missvalue_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                              & OFB(i1,j1,k1,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  DO k=ksr, ker
                     k1= k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                           IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) ) &
                           & OFB(i-hi,j-hj,k1,sample)= missvalue
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                           & OFB(i-hi,j-hj,k,sample)= missvalue
                        END DO
                     END DO
                  END DO
               END IF
            END IF
         END IF TIME_IF
      END ASSOCIATE

      succeded = .TRUE.
      RETURN

   END FUNCTION SAMPLE_THE_FIELD_IMP
END MODULE fms_send_data_statfun_mod
!> @}
! close documentation grouping
