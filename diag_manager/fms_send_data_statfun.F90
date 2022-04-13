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
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,&
   & fms_error_handler
   USE diag_data_mod, ONLY:  input_fields, output_fields, debug_diag_manager
   use diag_util_mod, ONLY: check_out_of_bounds

   TYPE STATFUN_IDX_CFG_T
      INTEGER :: f1,f2,f3,f4
      INTEGER ::  is, js, ks, ie, je, ke
      INTEGER :: hi !< halo size in x direction
      INTEGER :: hj !< halo size in y direction
   END TYPE STATFUN_IDX_CFG_T

   ABSTRACT INTERFACE
      FUNCTION weigh_field ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL :: weigh_field
      END FUNCTION weigh_field
   END INTERFACE

   ABSTRACT INTERFACE
      FUNCTION weigh_field_3d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:,:,:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d
      END FUNCTION weigh_field_3d
   END INTERFACE


   ABSTRACT INTERFACE
      FUNCTION weigh_field_1d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d
      END FUNCTION weigh_field_1d
   END INTERFACE

CONTAINS

   PURE FUNCTION weigh_field_0d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weigh_field_0d_p1 = field_val * weight
   END FUNCTION weigh_field_0d_p1

   PURE FUNCTION weigh_field_0d_p2 ( field_val, weight,  pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      REAL :: fTw
      fTw =  field_val * weight
      weigh_field_0d_p2 = fTw * fTw
   END FUNCTION weigh_field_0d_p2

   PURE FUNCTION weigh_field_0d_pp ( field_val, weight, pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weigh_field_0d_pp = (field_val * weight) ** pow_value
   END FUNCTION weigh_field_0d_pp

   PURE FUNCTION weigh_field_1d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_p1
      ALLOCATE(weigh_field_1d_p1(size(field_val)))
      DO i = 1, size(field_val)
         weigh_field_1d_p1(i) = field_val(i) * weight
      END DO
   END FUNCTION weigh_field_1d_p1

   PURE FUNCTION weigh_field_1d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_p2
      !!TODO :verify the allocate
      ALLOCATE(weigh_field_1d_p2, mold=field_val)
      DO i = 1, size(field_val)
         weigh_field_1d_p2(i) = field_val(i) * field_val(i) * weight * weight
      END DO
   END FUNCTION weigh_field_1d_p2

   PURE FUNCTION weigh_field_1d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_pp
      ALLOCATE(weigh_field_1d_pp(size(field_val)))
      DO i = 1, size(field_val)
         weigh_field_1d_pp(i) = (field_val(i) * weight) ** pow_value
      END DO
   END FUNCTION weigh_field_1d_pp

   PURE FUNCTION weigh_field_3d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p1
      ALLOCATE(weigh_field_3d_p1, mold=field_val)
      weigh_field_3d_p1 = field_val * weight
   END FUNCTION weigh_field_3d_p1

   PURE FUNCTION weigh_field_3d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p2
      ALLOCATE(weigh_field_3d_p2, mold=field_val)
      weigh_field_3d_p2 = field_val * field_val * weight * weight
   END FUNCTION weigh_field_3d_p2

   PURE FUNCTION weigh_field_3d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_pp
      ALLOCATE(weigh_field_3d_pp, mold=field_val)
      weigh_field_3d_pp = (field_val * weight) ** pow_value
   END FUNCTION weigh_field_3d_pp


   FUNCTION AVERAGE_THE_FIELD (diag_field_id, field, out_num, &
   & mask, weight1, sample, missvalue, missvalue_present, &
   & l_start, l_end, idx_cfg, err_msg,  err_msg_local) result( succeded )
      INTEGER, INTENT(in) :: diag_field_id
      REAL, DIMENSION(:,:,:), INTENT(in) :: field
      INTEGER, INTENT(in) :: out_num
      LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
      REAL, INTENT(in) :: weight1
      INTEGER, INTENT(in) :: sample
      REAL, INTENT(in), OPTIONAL :: missvalue   !!TODO: optional?
      LOGICAL, INTENT(in) :: missvalue_present
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
      CHARACTER(len=*), INTENT(inout), OPTIONAL :: err_msg
      TYPE(STATFUN_IDX_CFG_T), INTENT(in) :: idx_cfg
      CHARACTER(len=256), INTENT(inout) :: err_msg_local
      LOGICAL :: succeded

      !!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
      CHARACTER(len=128):: error_string

      ! Power value for rms or pow(x) calculations
      INTEGER :: pow_value, ksr, ker, is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
      LOGICAL :: phys_window , need_compute , reduced_k_range

      INTEGER :: i, j, k,  i1, j1, k1

      INTEGER :: numthreads
      INTEGER :: active_omp_level
#if defined(_OPENMP)
      INTEGER :: omp_get_num_threads !< OMP function
      INTEGER :: omp_get_level !< OMP function
#endif

      !!A pointer to the field weighing function.
      procedure (weigh_field), pointer :: fwf_0d_ptr => null ()
      !! A pointer to the 3D filed weighn function.
      procedure (weigh_field_1d), pointer :: fwf_1d_ptr => null ()
!! A pointer to the 3D filed weighn function.
      procedure (weigh_field_3d), pointer :: fwf_3d_ptr => null ()

      pow_value = output_fields(out_num)%pow_value

      if ( pow_value == 1 ) then
         fwf_0d_ptr => weigh_field_0d_p1
         fwf_1D_ptr => weigh_field_1d_p1
         fwf_3D_ptr => weigh_field_3d_p1
      else if ( pow_value == 2 ) then
         fwf_0d_ptr => weigh_field_0d_p2
         fwf_1d_ptr => weigh_field_1d_p2
         fwf_3D_ptr => weigh_field_3d_p2
      else
         fwf_0d_ptr => weigh_field_0d_pp
         fwf_1D_ptr => weigh_field_1d_pp
         fwf_3D_ptr => weigh_field_3d_pp
      end if

      phys_window = output_fields(out_num)%phys_window
      need_compute = output_fields(out_num)%need_compute
      reduced_k_range = output_fields(out_num)%reduced_k_range

      ksr= l_start(3)
      ker= l_end(3)
      is = idx_cfg%is
      js = idx_cfg%js
      ks = idx_cfg%ks
      ie = idx_cfg%ie
      je = idx_cfg%je
      ke = idx_cfg%ke
      hi = idx_cfg%hi
      hj = idx_cfg%hj
      f1 = idx_cfg%f1
      f2 = idx_cfg%f2
      f3 = idx_cfg%f3
      f4 = idx_cfg%f4

      ASSOCIATE( ofb => output_fields(out_num)%buffer , &
      & ofc => output_fields(out_num)%counter)

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
            IF ( PRESENT(mask) ) THEN
               IF ( missvalue_present ) THEN !!TODO: (section: section( mask_varian .eq. true + mask present) )
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
                  IF( numthreads>1 .AND. phys_window ) then
                     REDU_KR1_IF: IF ( reduced_k_range ) THEN
                        DO k= ksr, ker
                           k1= k - ksr + 1
                           DO j=js, je
                              DO i=is, ie
                                 IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                    ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) +&
                                    & fwf_0d_ptr (field(i-is+1+hi, j-js+1+hj, k), weight1, pow_value)
                                    ofc(i-hi,j-hj,k1,sample) = ofc(i-hi,j-hj,k1,sample) + weight1
                                 END IF
                              END DO
                           END DO
                        END DO
                     ELSE
                        DO k=ks, ke
                           DO j=js, je
                              DO i=is, ie
                                 IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                    ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) + &
                                    & fwf_0d_ptr (field(i-is+1+hi,j-js+1+hj,k), weight1,  pow_value)
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
                                    & fwf_0d_ptr (field(i-is+1+hi, j-js+1+hj, k),  weight1, pow_value)
                                    ofc(i-hi,j-hj,k1,sample) = ofc(i-hi,j-hj,k1,sample) + weight1
                                 END IF
                              END DO
                           END DO
                        END DO
                     ELSE
                        DO k=ks, ke
                           DO j=js, je
                              DO i=is, ie
                                 IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                    ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) +  &
                                    & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                                    ofc(i-hi,j-hj,k,sample) = ofc(i-hi,j-hj,k,sample) + weight1
                                 END IF
                              END DO
                           END DO
                        END DO
                     END IF REDU_KR2_IF
!$OMP END CRITICAL
                  END IF
               ELSE
                  WRITE (error_string,'(a,"/",a)')&
                  & TRIM(input_fields(diag_field_id)%module_name), &
                  & TRIM(output_fields(out_num)%output_name)
                  IF(fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                  & ', variable mask but no missing value defined', err_msg)) THEN
                     succeded = .FALSE.
                     RETURN
                  END IF
               END IF
            ELSE  ! no mask present
               WRITE (error_string,'(a,"/",a)')&
               & TRIM(input_fields(diag_field_id)%module_name), &
               & TRIM(output_fields(out_num)%output_name)
               IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
               & ', variable mask but no mask given', err_msg)) THEN
                  succeded = .FALSE.
                  RETURN
               END IF
            END IF
         ELSE !! MASK_VAR_IF
            MASK_PRESENT_IF: IF ( PRESENT(mask) ) THEN
               MVAL_PRESENT_IF: IF ( missvalue_present ) THEN !!section:(mask_var false +mask present +missval prsnt)
                  IF ( need_compute ) THEN
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
                                       & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k),  weight1, pow_value)
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
                                       & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                                 output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  ELSE IF ( reduced_k_range ) THEN
                     IF (numthreads>1 .AND. phys_window) then
                        DO k=ksr, ker
                           k1 = k - ksr + 1
                           DO j=js, je
                              DO i=is, ie
                                 IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                    ofb(i-hi,j-hj,k1,sample) = ofb(i-hi,j-hj,k1,sample) + &
                                    & fwf_0d_ptr (field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                                    & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k) ,weight1, pow_value)
                                 ELSE
                                    ofb(i-hi,j-hj,k1,sample)= missvalue
                                 END IF
                              END DO
                           END DO
                        END DO
!$OMP END CRITICAL
                     END IF
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
                     IF (numthreads>1 .AND. phys_window) then
                        DO k=ks, ke
                           DO j=js, je
                              DO i=is, ie
                                 IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                    ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) + &
                                    & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                                    & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                                 ELSE
                                    ofb(i-hi,j-hj,k,sample)= missvalue
                                 END IF
                              END DO
                           END DO
                        END DO
!$OMP END CRITICAL
                     END IF
                  END IF
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
               ELSE !!MVAL_PRESENT_IF (section: mask_varian .eq. false + mask present + miss value not present)
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
                  IF ( need_compute ) THEN
                     IF (numthreads>1 .AND. phys_window) then
                        DO j = js, je !!TODO: What is diff in next two secs?
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1 =  j-l_start(2)-hj+1
                                 ofb(i1,j1,:,sample)=  ofb(i1,j1,:,sample)+ &
                                 & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
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
                                 & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
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
                  ELSE IF ( reduced_k_range ) THEN
                     IF (numthreads>1 .AND. phys_window) then
                        ksr= l_start(3)  !!TODO. What is deff in these two secs below.
                        ker= l_end(3)
                        ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) = ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                        & fwf_3d_ptr (field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
                     ELSE
!$OMP CRITICAL
                        ksr= l_start(3)
                        ker= l_end(3)
                        ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) = ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                        & fwf_3D_ptr (field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
!$OMP END CRITICAL
                     END IF
                  ELSE
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
                     !!TODO: What is difference in two below?
                     IF (numthreads>1 .AND. phys_window) then
                        ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                        & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                        & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
                     ELSE
!$OMP CRITICAL
                        ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                        & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                           fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1 , pow_value)
!$OMP END CRITICAL
                     END IF
                  END IF
!$OMP CRITICAL
                  IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                  & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
               END IF MVAL_PRESENT_IF
            ELSE !!MASK_PRESENT_IF (section: mask_variant .eq. false + mask not present + missvalue)
               MVAL_PRESENT2_IF: IF (missvalue_present ) THEN
                  IF ( need_compute ) THEN
                     if( numthreads>1 .AND. phys_window ) then
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
                                       & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                           k1 = k - l_start(3) + 1
                           DO j = js, je
                              DO i = is, ie
                                 IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                                 & j <= l_end(2)+hj) THEN
                                    i1 = i-l_start(1)-hi+1
                                    j1=  j-l_start(2)-hj+1
                                    IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                       ofb(i1,j1,k1,sample) = ofb(i1,j1,k1,sample) + &
                                       & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                                    ELSE
                                       ofb(i1,j1,k1,sample) = missvalue
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END DO
!$OMP END CRITICAL
                     END IF
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
                        outer0: DO k = l_start(3), l_end(3)
                           DO j=l_start(2)+hj, l_end(2)+hj
                              DO i=l_start(1)+hi, l_end(1)+hi
                                 IF ( field(i,j,k) /= missvalue ) THEN
                                    output_fields(out_num)%count_0d(sample) = &
                                    output_fields(out_num)%count_0d(sample) + weight1
                                    EXIT outer0
                                 END IF
                              END DO
                           END DO
                        END DO outer0
                     END IF
!$OMP END CRITICAL
                  ELSE IF ( reduced_k_range ) THEN
                     if( numthreads>1 .AND. phys_window ) then
                        ksr= l_start(3)
                        ker= l_end(3)
                        DO k = ksr, ker
                           k1 = k - ksr + 1
                           DO j=js, je
                              DO i=is, ie
                                 IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                    ofb(i-hi,j-hj,k1,sample) =  ofb(i-hi,j-hj,k1,sample) + &
                                    & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                                    & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                                 ELSE
                                    ofb(i-hi,j-hj,k1,sample) = missvalue
                                 END IF
                              END DO
                           END DO
                        END DO
!$OMP END CRITICAL
                     END IF
!$OMP CRITICAL
                     outer3: DO k = ksr, ker
                        k1=k-ksr+1
                        DO j=f3, f4
                           DO i=f1, f2
                              !! TODO: verify this below
                              IF ( field(i,j,k) /= missvalue ) THEN
                                 output_fields(out_num)%count_0d(sample) = &
                                 output_fields(out_num)%count_0d(sample) + weight1
                                 EXIT outer3
                              END IF
                           END DO
                        END DO
                     END DO outer3
!$OMP END CRITICAL
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
                     IF( numthreads > 1 .AND. phys_window ) then
                        DO k=ks, ke
                           DO j=js, je
                              DO i=is, ie
                                 IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                    ofb(i-hi,j-hj,k,sample) = ofb(i-hi,j-hj,k,sample) +&
                                    & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
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
                                    & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                                 ELSE
                                    ofb(i-hi,j-hj,k,sample) = missvalue
                                 END IF
                              END DO
                           END DO
                        END DO
!$OMP END CRITICAL
                     END IF
!$OMP CRITICAL
                     outer1: DO k=ks, ke
                        DO j=f3, f4
                           DO i=f1, f2
                              IF ( field(i,j,k) /= missvalue ) THEN
                                 output_fields(out_num)%count_0d(sample) = &
                                 output_fields(out_num)%count_0d(sample)  + weight1
                                 EXIT outer1
                              END IF
                           END DO
                        END DO
                     END DO outer1
!$OMP END CRITICAL
                  END IF
               ELSE !MVAL_PRESENT2_IF (section:  mask_variant .eq. false + mask not present + missvalue not present)
                  NEED_COMP_IF: IF ( need_compute ) THEN
                     IF( numthreads > 1 .AND. phys_window ) then
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 ofb(i1,j1,:,sample) = ofb(i1,j1,:,sample) + &
                                 & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
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
                                 & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
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
                  ELSE IF ( reduced_k_range ) THEN !!NEED_COMP_IF
                     ksr= l_start(3)
                     ker= l_end(3)
                     IF( numthreads > 1 .AND. phys_window ) then
                        ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                        & ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                        & fwf_3d_ptr(field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
                     ELSE
!$OMP CRITICAL
                        ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                        & ofb(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                        & fwf_3d_ptr(field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
!$OMP END CRITICAL
                     END IF

                  ELSE
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
                        & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
                     ELSE
!$OMP CRITICAL
                        ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                        & ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                        & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
                        !!
!$OMP END CRITICAL
                     END IF
                  END IF NEED_COMP_IF
!$OMP CRITICAL
                  IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                  & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
               END IF MVAL_PRESENT2_IF
            END IF MASK_PRESENT_IF ! if mask present
         END IF MASK_VAR_IF
      END ASSOCIATE

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

   END FUNCTION AVERAGE_THE_FIELD

   FUNCTION SAMPLE_THE_FIELD (diag_field_id, field, out_num, &
    & mask, sample, missvalue, missvalue_present, &
    & l_start, l_end, idx_cfg, err_msg,  err_msg_local) result( succeded )
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
       TYPE(STATFUN_IDX_CFG_T), INTENT(in) :: idx_cfg
       CHARACTER(len=256), INTENT(inout) :: err_msg_local
       LOGICAL :: succeded

       !!LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
       CHARACTER(len=128):: error_string

       ! Power value for rms or pow(x) calculations
       INTEGER :: pow_value, ksr, ker, is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
       LOGICAL :: phys_window , need_compute , reduced_k_range

       INTEGER :: i, j, k, i1, j1, k1

       LOGICAL :: time_rms, time_max, time_min, time_sum
#if defined(_OPENMP)
      INTEGER :: omp_get_num_threads !< OMP function
      INTEGER :: omp_get_level !< OMP function
#endif

!!A pointer to the field weighing function.
      procedure (weigh_field), pointer :: fwf_0d_ptr => null ()
!! A pointer to the 3D filed weighn function.
      procedure (weigh_field_1d), pointer :: fwf_1d_ptr => null ()
!! A pointer to the 3D filed weighn function.
      procedure (weigh_field_3d), pointer :: fwf_3d_ptr => null ()

      pow_value = output_fields(out_num)%pow_value

      if ( pow_value == 1 ) then
         fwf_0d_ptr => weigh_field_0d_p1
         fwf_1D_ptr => weigh_field_1d_p1
         fwf_3D_ptr => weigh_field_3d_p1
      else if ( pow_value == 2 ) then
         fwf_0d_ptr => weigh_field_0d_p2
         fwf_1d_ptr => weigh_field_1d_p2
         fwf_3D_ptr => weigh_field_3d_p2
      else
         fwf_0d_ptr => weigh_field_0d_pp
         fwf_1D_ptr => weigh_field_1d_pp
         fwf_3D_ptr => weigh_field_3d_pp
      end if

      phys_window = output_fields(out_num)%phys_window
      need_compute = output_fields(out_num)%need_compute
      reduced_k_range = output_fields(out_num)%reduced_k_range

      ! Is this output field the rms?
      ! If so, then average is also .TRUE.
      time_rms = output_fields(out_num)%time_rms
      ! Looking for max and min value of this field over the sampling interval?
      time_max = output_fields(out_num)%time_max
      time_min = output_fields(out_num)%time_min
      ! Sum output over time interval
      time_sum = output_fields(out_num)%time_sum

      ksr= l_start(3)
      ker= l_end(3)
      is = idx_cfg%is
      js = idx_cfg%js
      ks = idx_cfg%ks
      ie = idx_cfg%ie
      je = idx_cfg%je
      ke = idx_cfg%ke
      hi = idx_cfg%hi
      hj = idx_cfg%hj
      f1 = idx_cfg%f1
      f2 = idx_cfg%f2
      f3 = idx_cfg%f3
      f4 = idx_cfg%f4

      ASSOCIATE( OFB => output_fields(out_num)%buffer)

         ! Add processing for Max and Min
         TIME_IF: IF ( time_max ) THEN
            IF ( PRESENT(mask) ) THEN
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
            ELSE
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
            END IF
            output_fields(out_num)%count_0d(sample) = 1
            !END TIME MAX
         ELSE IF ( time_min ) THEN   !!TiME_IF
            IF ( PRESENT(mask) ) THEN
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
            ELSE
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
            END IF
            output_fields(out_num)%count_0d(sample) = 1

            !! END_TIME_MIN
         ELSE IF ( time_sum ) THEN  !!TIME_IF
            IF ( PRESENT(mask) ) THEN
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
            ELSE
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
            END IF
            output_fields(out_num)%count_0d(sample) = 1
            !END time_sum
         ELSE  ! ( not average, not min, not max, not sum )  !!TIME_IF
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

   END FUNCTION SAMPLE_THE_FIELD
END MODULE fms_send_data_statfun_mod
