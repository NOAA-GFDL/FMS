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

  TYPE STATFUN_CFG_T
    REAL :: weight1
    REAL :: missvalue
    INTEGER :: f1,f2,f3,f4

    INTEGER ::  is, js, ks, ie, je, ke
    INTEGER :: hi !< halo size in x direction
    INTEGER :: hj !< halo size in y direction
    INTEGER :: twohi !< halo size in x direction
    INTEGER :: twohj !< halo size in y direction
    INTEGER :: sample !< index along the diurnal time axis
    LOGICAL :: missvalue_present

  END TYPE STATFUN_CFG_T

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
  weigh_field_0d_pp = (field_val * weight) ** pow_val
  END FUNCTION weigh_field_0d_pp

PURE FUNCTION weigh_field_1d_p1 ( field_val, weight, pow_value )
  REAL, INTENT(in) :: field_val(:)
  REAL, INTENT(in) :: weight
  INTEGER, INTENT(in) :: pow_value
  INTEGER :: i
  REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_p1
  ALLOCATE(weigh_field_1d_p1(size(field_val)))
  DO i = i, size(field_val)
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
  DO i = i, size(field_val)
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
  DO i = i, size(field_val)
    weigh_field_1d_pp(i) = (field_val(i) * weight) ** pow_value
  END DO
  END FUNCTION weigh_field_1d_pp

  PURE FUNCTION weigh_field_3d_p1 ( field_val, weight, pow_value )
  REAL, INTENT(in) :: field_val(:,:,:)
  REAL, INTENT(in) :: weight
  INTEGER, INTENT(in) :: pow_value
  INTEGER :: i
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p1
  ALLOCATE(weigh_field_3d_p1, mold=field_val)
    weigh_field_3d_p1 = field_val * weight
  END FUNCTION weigh_field_3d_p1

  PURE FUNCTION weigh_field_3d_p2 ( field_val, weight, pow_value )
  REAL, INTENT(in) :: field_val(:,:,:)
  REAL, INTENT(in) :: weight
  INTEGER, INTENT(in) :: pow_value
  INTEGER :: i
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p2
  ALLOCATE(weigh_field_3d_p2, mold=field_val)
  weigh_field_3d_p2 = field_val * fiedl_val * weight * weight
  END FUNCTION weigh_field_3d_p2

  PURE FUNCTION weigh_field_3d_pp ( field_val, weight, pow_value )
  REAL, INTENT(in) :: field_val(:,:,:)
  REAL, INTENT(in) :: weight
  INTEGER, INTENT(in) :: pow_value
  INTEGER :: i
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_pp
  ALLOCATE(weigh_field_3d_pp, mold=field_val)
    weigh_field_3d_pp = (field_val * weight) ** pow_val
  END FUNCTION weigh_field_3d_pp


  FUNCTION AVERAGE_THE_FIELD (diag_field_id, field, out_num, mask, weight, &
      & l_start, l_end, cfg, err_msg,  err_msg_local) result( succeded )
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: field
    INTEGER, INTENT(in) :: out_num
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
    INTEGER, DIMENSION(3) :: l_start !< local start indices on 3 axes for regional output
    INTEGER, DIMENSION(3) :: l_end !< local end indices on 3 axes for regional output
    CHARACTER(len=*), INTENT(inout), OPTIONAL :: err_msg
    TYPE(STATFUN_CFG_T), INTENT(in) :: cfg
    CHARACTER(len=256), INTENT(inout) :: err_msg_local
    LOGICAL :: succeded

    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
    CHARACTER(len=128):: error_string


    ! Power value for rms or pow(x) calculations
    INTEGER :: pow_value , ksr, ker, is, js, ks, ie, je, ke, hi, hj, sample, f1, f2, f3, f4
    LOGICAL :: phys_window , need_compute , reduced_k_range
    REAL :: weight1, missvalue

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

  if ( pow_val == 1 ) then
    fwf_0d_ptr => weigh_field_0d_p1
    fwf_1D_ptr => weigh_field_1d_p1
    fwf_3D_ptr => weigh_field_3d_p1
  else if ( pow_val == 2 ) then
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
    is = cfg%is
    js = cfg%js
    ks = cfg%ks
    ie = cfg%ie
    je = cfg%je
    ke = cfg%ke
    hi = cfg%hi
    hj = cfg%hj
    sample = cfg%sample
    f1 = cfg%f1
    f2 = cfg%f2
    f3 = cfg%f3
    f4 = cfg%f4
    weight1 = cfg%weight1
    missvalue = cfg%missvalue



    IF ( input_fields(diag_field_id)%mask_variant ) THEN
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
        IF ( cfg%missvalue_present ) THEN
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
            IF ( reduced_k_range ) THEN
              NORMAL1: DO k= ksr, ker
                k1= k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                      IF ( pow_value /= 1 ) THEN
                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                          & (field(i-is+1+hi, j-js+1+hj, k) * weight1)**(pow_value)
                      ELSE
                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                          & field(i-is+1+hi, j-js+1+hj, k) * weight1
                      END IF
                      output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                        & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                    END IF
                  END DO
                END DO
              END DO NORMAL1
            ELSE
              NORMAL2: DO k=ks, ke
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN  !!Includes removal of ( pow_value /= 1 ) ifhten
                        ASSOCIATE( ofb => output_fields(out_num)%buffer (i-hi,j-hj,k,sample) , &
                          & ofc => output_fields(out_num)%counter(i-hi,j-hj,k,sample))
                          ofb = ofb + (field(i-is+1+hi,j-js+1+hj,k) * weight1) ** pow_value
                          ofc = ofc + weight1
                        END ASSOCIATE
                    END IF
                  END DO
                END DO
              END DO NORMAL2
            END IF
          ELSE  !$OMP CRITICAL
            IF ( reduced_k_range ) THEN
              DO k= ksr, ker
                k1= k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) + &
                          & fwf_0d_ptr (field(i-is+1+hi, j-js+1+hj, k),  weight1, pow_value)
                      output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                        & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                    END IF
                  END DO
                END DO
              END DO
            ELSE
              DO k=ks, ke
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = &
                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +  &
                          & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                      output_fields(out_num)%counter(i-hi,j-hj,k,sample) =&
                        &output_fields(out_num)%counter(i-hi,j-hj,k,sample) + weight1
                    END IF
                  END DO
                END DO
              END DO
            END IF  !$OMP END CRITICAL
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
    ELSE ! mask_variant=false
      IF ( PRESENT(mask) ) THEN
        IF ( cfg%missvalue_present ) THEN
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
                          output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                            & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                            & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k),  weight1, pow_value)
                      ELSE
                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                      END IF
                    END IF
                  END DO
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO k = l_start(3), l_end(3)
                k1 = k-l_start(3)+1
                DO j = js, je
                  DO i = is, ie
                    IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                      & j <= l_end(2)+hj ) THEN
                      i1 = i-l_start(1)-hi+1
                      j1=  j-l_start(2)-hj+1
                      IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                          output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                            & output_fields(out_num)%buffer(i1,j1,k1,sample) + &
                            & fwf_0d_ptr( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                      ELSE
                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                      END IF
                    END IF
                  END DO
                END DO
              END DO  !$OMP END CRITICAL
            ENDIF  !$OMP CRITICAL
            DO j = js, je
              DO i = is, ie
                IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                  & j <= l_end(2)+hj ) THEN
                  output_fields(out_num)%num_elements(sample) = &
                    output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                END IF
              END DO
            END DO  !$OMP END CRITICAL
          ELSE IF ( reduced_k_range ) THEN
            IF (numthreads>1 .AND. phys_window) then
              DO k=ksr, ker
                k1 = k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)+ &
                          & fwf_0d_ptr (field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                    END IF
                  END DO
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO k=ksr, ker
                k1 = k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                        & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)+ &
                        & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k) ,weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                    END IF
                  END DO
                END DO
              END DO  !$OMP END CRITICAL
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
                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = &
                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)+ &
                          & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                    END IF
                  END DO
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO k=ks, ke
                DO j=js, je
                  DO i=is, ie
                    IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                        & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)+ &
                        & fwf_0d_ptr ( field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                    END IF
                  END DO
                END DO
              END DO  !$OMP END CRITICAL
            END IF
          END IF  !$OMP CRITICAL
          IF ( need_compute .AND. .NOT.phys_window ) THEN
            IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3))) ) &
              & output_fields(out_num)%count_0d(sample) =&
              & output_fields(out_num)%count_0d(sample) + weight1
          ELSE
            IF ( ANY(mask(f1:f2,f3:f4,ks:ke)) ) output_fields(out_num)%count_0d(sample) =&
              & output_fields(out_num)%count_0d(sample)+weight1
          END IF!$OMP END CRITICAL

        ELSE ! missing value NOT present
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
              DO j = js, je
                DO i = is, ie
                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                    & j <= l_end(2)+hj ) THEN
                    i1 = i-l_start(1)-hi+1
                    j1 =  j-l_start(2)-hj+1
                      output_fields(out_num)%buffer(i1,j1,:,sample)= &
                        & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                        & weigh_field_pp_1d(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1,pow_value)
                        !& (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)

                  END IF
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO j = js, je
                DO i = is, ie
                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                    & j <= l_end(2)+hj ) THEN
                    i1 = i-l_start(1)-hi+1
                    j1 =  j-l_start(2)-hj+1
                    output_fields(out_num)%buffer(i1,j1,:,sample)= &
                        & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                        & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)),weight1, pow_value)
                        !!& (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)

                  END IF
                END DO
              END DO  !$OMP END CRITICAL
            END IF  !$OMP CRITICAL
            DO j = js, je
              DO i = is, ie
                IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                  & j <= l_end(2)+hj ) THEN
                  output_fields(out_num)%num_elements(sample)=&
                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1

                END IF
              END DO
            END DO !$OMP END CRITICAL
          ELSE IF ( reduced_k_range ) THEN
            IF (numthreads>1 .AND. phys_window) then
              ksr= l_start(3)
              ker= l_end(3)
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                  & fwf_3d_ptr (field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
                  !!& (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)

            ELSE  !$OMP CRITICAL
              ksr= l_start(3)
              ker= l_end(3)
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                  & fwf_3D_ptr( field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
              END IF  !$OMP END CRITICAL
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
            !!TODO: DWhat is difference in two below?
            IF (numthreads>1 .AND. phys_window) then
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                  & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
            ELSE  !$OMP CRITICAL
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                  fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1 , pow_value)
                  !!& (field(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
            END IF
          END IF  !$OMP CRITICAL
          IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
            & output_fields(out_num)%count_0d(sample) + weight1  !$OMP END CRITICAL
        END IF
      ELSE ! mask NOT present
        IF ( cfg%missvalue_present ) THEN
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
                          output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                            & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                            & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                      ELSE
                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                      END IF
                    END IF
                  END DO
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO k = l_start(3), l_end(3)
                k1 = k - l_start(3) + 1
                DO j = js, je
                  DO i = is, ie
                    IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                      & j <= l_end(2)+hj) THEN
                      i1 = i-l_start(1)-hi+1
                      j1=  j-l_start(2)-hj+1
                      IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                        output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                          & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                          & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                          !& (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                      ELSE
                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                      END IF
                    END IF
                  END DO
                END DO
              END DO  !$OMP END CRITICAL
            END IF  !$OMP CRITICAL
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
                      output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample)&
                        & + weight1
                      EXIT outer0
                    END IF
                  END DO
                END DO
              END DO outer0
            END IF !$OMP END CRITICAL
          ELSE IF ( reduced_k_range ) THEN
            if( numthreads>1 .AND. phys_window ) then
              ksr= l_start(3)
              ker= l_end(3)
              DO k = ksr, ker
                k1 = k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                          & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                          !& (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)

                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                    END IF
                  END DO
                END DO
              END DO
            else !$OMP CRITICAL
              ksr= l_start(3)
              ker= l_end(3)
              DO k = ksr, ker
                k1 = k - ksr + 1
                DO j=js, je
                  DO i=is, ie
                    IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                        & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                        & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                    END IF
                  END DO
                END DO
              END DO !$OMP END CRITICAL
            END IF !$OMP CRITICAL
            outer3: DO k = ksr, ker
              k1=k-ksr+1
              DO j=f3, f4
                DO i=f1, f2
                  IF ( field(i,j,k) /= missvalue ) THEN
                    output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) &
                      & + weight1
                    EXIT outer3
                  END IF
                END DO
              END DO
            END DO outer3 !$OMP END CRITICAL
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
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                          & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                          !& (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                    END IF
                  END DO
                END DO
              END DO
            ELSE  !$OMP CRITICAL
              DO k=ks, ke
                DO j=js, je
                  DO i=is, ie
                    IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                          & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                          & fwf_0d_ptr(field(i-is+1+hi,j-js+1+hj,k), weight1, pow_value)
                    ELSE
                      output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                    END IF
                  END DO
                END DO
              END DO !$OMP END CRITICAL
            END IF !$OMP CRITICAL
            outer1: DO k=ks, ke
              DO j=f3, f4
                DO i=f1, f2
                  IF ( field(i,j,k) /= missvalue ) THEN
                    output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) &
                      & + weight1
                    EXIT outer1
                  END IF
                END DO
              END DO
            END DO outer1 !$OMP END CRITICAL
          END IF
        ELSE ! no missing value defined, No mask
          IF ( need_compute ) THEN
            IF( numthreads > 1 .AND. phys_window ) then
              DO j = js, je
                DO i = is, ie
                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                    & j <= l_end(2)+hj ) THEN
                    i1 = i-l_start(1)-hi+1
                    j1=  j-l_start(2)-hj+1
                    output_fields(out_num)%buffer(i1,j1,:,sample)= &
                        & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                        & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
                        !!& (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)

                  END IF
                END DO
              END DO
            ELSE !$OMP CRITICAL
              DO j = js, je
                DO i = is, ie
                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                    & j <= l_end(2)+hj ) THEN
                    i1 = i-l_start(1)-hi+1
                    j1=  j-l_start(2)-hj+1
                    output_fields(out_num)%buffer(i1,j1,:,sample)= &
                        & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                        & fwf_1d_ptr(field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3)), weight1, pow_value)
                  END IF
                END DO
              END DO  !$OMP END CRITICAL
            END IF  !$OMP CRITICAL
            DO j = js, je
              DO i = is, ie
                IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                  & j <= l_end(2)+hj ) THEN
                  output_fields(out_num)%num_elements(sample) =&
                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                END IF
              END DO
            END DO  !$OMP END CRITICAL
            ! Accumulate time average
          ELSE IF ( reduced_k_range ) THEN
            ksr= l_start(3)
            ker= l_end(3)
            IF( numthreads > 1 .AND. phys_window ) then
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                  & fwf_3d_ptr(field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
            ELSE  !$OMP CRITICAL
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                  & fwf_3d_ptr(field(f1:f2,f3:f4,ksr:ker), weight1, pow_value)
                  !& (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
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
            IF( numthreads > 1 .AND. phys_window ) then
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                  & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
            ELSE  !$OMP CRITICAL
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                  & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                  & fwf_3d_ptr(field(f1:f2,f3:f4,ks:ke), weight1, pow_value)
              !!  !$OMP END CRITICAL
            END IF
          END IF  !$OMP CRITICAL
          IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
            & output_fields(out_num)%count_0d(sample) + weight1  !$OMP END CRITICAL
        END IF
      END IF ! if mask present
    END IF  !if mask_variant  !$OMP CRITICAL
    IF ( .NOT.need_compute .AND. .NOT.reduced_k_range )&
      & output_fields(out_num)%num_elements(sample) =&
      & output_fields(out_num)%num_elements(sample) + (ie-is+1)*(je-js+1)*(ke-ks+1)
    IF ( reduced_k_range ) &
      & output_fields(out_num)%num_elements(sample) = output_fields(out_num)%num_elements(sample) +&
      & (ie-is+1)*(je-js+1)*(ker-ksr+1)  !$OMP END CRITICAL

      succeded = .TRUE.
      RETURN

  END FUNCTION AVERAGE_THE_FIELD
END MODULE fms_send_data_statfun_mod
