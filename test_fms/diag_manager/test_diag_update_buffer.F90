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

!> @brief  This programs tests the update of field data buffers with
!! the  "math" functions in module fms_diag_fieldbuff_update_mod. It mimics
!! the daig_manager::send_4d operation of calling those functions.
program test_diag_update_buffer
   use platform_mod
   use mpp_mod, only: mpp_init, mpp_set_stack_size, mpp_init_test_requests_allocated
   use fms_mod, ONLY: fms_init, fms_end, error_mesg, FATAL,NOTE
   use diag_data_mod, ONLY: VERY_LARGE_AXIS_LENGTH
   USE fms_diag_outfield_mod, ONLY: fmsDiagOutfield_type, fmsDiagOutfieldIndex_type
   USE fms_diag_fieldbuff_update_mod, ONLY: fieldbuff_update, fieldbuff_copy_missvals, &
   & fieldbuff_copy_fieldvals
   USE fms_diag_time_reduction_mod, ONLY: fmsDiagTimeReduction_type, time_average, time_rms
   USE fms_diag_bbox_mod, ONLY: fmsDiagIbounds_type

   implicit  none

   !! Class diag_buffer_type is here only for temporary use for modern diag_manager
   !! development until the real buffer class is sufficiently ready and merged.
   TYPE diagTestBuffer_type
      CLASS(*), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: buffer
      CLASS(*), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: counter
      CLASS(*), ALLOCATABLE, DIMENSION(:)  :: count_0d
      INTEGER, ALLOCATABLE, dimension(:) :: num_elements
   END TYPE diagTestBuffer_type

   integer,parameter :: SZ=10    !< Field data this size in all spatiall dims.
   integer,parameter :: SL=2     !< Field data this size in 4th dim
   integer,parameter :: NDI=1    !< Number of diurnal elemes
   CLASS(*), ALLOCATABLE :: r4_datapoint, i8_datapoint !< to be allocated of rype data (e.g. r4. i8)
                              !! to be used thought.

   TYPE(fmsDiagIbounds_type) ::  buff_bounds

   !!Diag_manager::send_data uses CLASS(*) in function signature, SO
   !! we mimic the resulting operations. The set of ClASS(*) data needs to be allocated of same
   !! type in order to be able to call the math/buffer update functions.
   CLASS(*), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: field_data
   CLASS(*), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: rmask
   CLASS(*), ALLOCATABLE, TARGET :: missvalue
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: mask
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: oor_mask
   TYPE(diagTestBuffer_type), ALLOCATABLE, TARGET :: buff_obj

   !! In principle, the field_data can be r4,r8,i4,i8,but we will only rest r4,i8
   !!These belwo will be pointers to the data
   REAL (kind=r4_kind),dimension (:,:,:,:),pointer::field_r4_ptr => null() !< Ptr to r4 field data array
   REAL (kind=r4_kind),dimension (:,:,:,:),pointer::rmask_r4_ptr => null() !< Ptr to r4 field data rmask array
   REAL (kind=r4_kind),pointer::missval_r4_ptr => null() !< Ptr to r4 missing value data.
   INTEGER (kind=i8_kind),dimension (:,:,:,:),pointer::field_i8_ptr => null() !< Ptr to i8 field data array
   INTEGER (kind=i8_kind),dimension (:,:,:,:),pointer::rmask_i8_ptr => null() !< Ptr to i8 field data rmask array
   INTEGER (kind=i8_kind),pointer::missval_i8_ptr => null()  !< Ptr to i8 missing value data.

   !! Typed pointers to buffer class(*) data will be needed
   REAL (kind=r4_kind),dimension (:,:,:,:,:),pointer::ofb_r4_ptr => null() !<Ptr to r4 buffer member of buffer obj.
   REAL (kind=r4_kind),dimension (:,:,:,:,:),pointer::ofc_r4_ptr => null()  !<Ptr to r4 counter  member of buffer obj.
   REAL (kind=r4_kind),dimension (:),pointer::ofb0d_r4_ptr => null() !< Ptr to r4 count0d member of buffer obj.
   !! Typed pointers to buffer class(*) data will be needed
   INTEGER (kind=i8_kind),dimension (:,:,:,:,:),pointer::ofb_i8_ptr => null() !<Ptr to i8 buffer member of buffer obj.
   INTEGER (kind=i8_kind),dimension (:,:,:,:,:),pointer::ofc_i8_ptr => null() !<Ptr to i8 counter member of buffer obj.
   INTEGER (kind=i8_kind),dimension (:),pointer::ofb0d_i8_ptr => null() !<Ptr to i8 count0d  member of buffer obj.

   integer :: ierr               !< An error flag (required by some mpp routines)
   logical :: test_passed        !< Flag indicating if the test_passed
   logical :: temp_result !< Set to result of one of the update functions.

   CHARACTER(LEN=*), PARAMETER :: module_name1 = "modX" !< Some dummy valuel
   CHARACTER(LEN=*), PARAMETER:: field_name1 = "fieldX" !< Some dummy valuel
   CHARACTER(LEN=*), PARAMETER :: output_name1 = "ofieldX" !< Some dummy valuel
   INTEGER :: time_reduction_type1 = time_average !<The first time reduction type.
   INTEGER :: output_freq1 = 1 !< The firs output frequencey type.
   INTEGER :: pow_value     ! Power value for rms or pow(x) calculations
   LOGICAL :: phys_window   ! OMP subsetted data? See output_fields
   LOGICAL :: need_compute  ! if local_output, does the current PE take part in send_data?
   LOGICAL :: mask_variant
   INTEGER :: num_elems
   LOGICAL :: reduced_k_range
   TYPE(fmsDiagTimeReduction_type), allocatable :: time_reduction !!Replaces LOGICAL::time_rms,time_max,time_min...

   INTEGER:: diag_field_id
   INTEGER:: sample !!diurnal_index
   REAL(kind=r4_kind) ::   weight
   INTEGER:: hi, hj  !!for halo sizes
   integer num_threads
   integer active_omp_level
   logical issued_mask_ignore_warning

   CHARACTER(len=256) :: err_msg, err_msg_local
   integer, dimension(3) :: l_start, l_end

   LOGICAL :: missvalue_present = .false.

   TYPE(fmsDiagOutfield_type), ALLOCATABLE :: ofield_cfg
   TYPE(fmsDiagOutfieldIndex_type), ALLOCATABLE :: ofield_index_cfg

   call fms_init

   call error_mesg('test_update_buffers_with_field', 'Test has started',NOTE)

   !! Allocate the field data and associated data
   allocate(real(kind=r4_kind) :: r4_datapoint)
   allocate(integer(kind=i8_kind) :: i8_datapoint)

   call allocate_input_data_and_ptrs(r4_datapoint, field_data, rmask, missvalue, mask, SZ,SZ,SZ,SL)

   call allocate_buffer_obj(r4_datapoint, buff_obj, SZ, SZ, SZ, SL, NDI)
   call buff_bounds%reset(VERY_LARGE_AXIS_LENGTH, 0)

   call init_field_values (field_data)

   !!TODO: (MDM) Can switch to final diang_manager buffer_object type in modern diag effort.

   !!In this version, we will meerely set type specific pointers to data. Some will be
   !! null, but at the end either the r4 pointers are non-null or the i8 pointers are not null
   SELECT TYPE ( field_data )
    TYPE IS (real(kind=r4_kind))
      SELECT TYPE ( rmask )
       TYPE IS (real(kind=r4_kind))
         SELECT TYPE ( missvalue )
          TYPE IS (real(kind=r4_kind))
            field_r4_ptr => field_data
            rmask_r4_ptr => rmask
            missval_r4_ptr => missvalue
         END SELECT
      END SELECT
    TYPE IS (integer(kind=i8_kind))
      SELECT TYPE ( rmask )
       TYPE IS (INTEGER(kind=i8_kind))
         SELECT TYPE ( missvalue )
          TYPE IS (INTEGER(kind=i8_kind))
            field_i8_ptr => field_data
            rmask_i8_ptr => rmask
            missval_i8_ptr => missvalue
         END SELECT
      END SELECT
    CLASS DEFAULT
      CALL error_mesg ('test_update_buffers_with_field','ptr assignemnt unsupported type', FATAL)
   END SELECT

   SELECT TYPE (  ofb => buff_obj%buffer )
    TYPE IS (real(kind=r4_kind))
      SELECT TYPE ( ofc => buff_obj%counter )
       TYPE IS (real(kind=r4_kind))
         SELECT TYPE ( ofb0d => buff_obj%count_0d )
          TYPE IS (real(kind=r4_kind))
            ofb_r4_ptr =>  ofb
            ofc_r4_ptr =>  ofc
            ofb0d_r4_ptr => ofb0d
         END SELECT
      END SELECT
    TYPE IS (integer(kind=i8_kind))
      SELECT TYPE ( ofc => buff_obj%counter )
       TYPE IS (INTEGER(kind=i8_kind))
         SELECT TYPE ( ofb0d => buff_obj%count_0d )
          TYPE IS (INTEGER(kind=i8_kind))
            ofb_i8_ptr =>  ofb
            ofc_i8_ptr =>  ofc
            ofb0d_i8_ptr => ofb0d
         END SELECT
      END SELECT
    CLASS DEFAULT
      CALL error_mesg ('diag_manager_mod::send_data_4d', 'ptr assigenment error', FATAL)
   END SELECT


   diag_field_id = 1
   sample = 1
   weight = 1.0
   missvalue = 1.0e-5_r4_kind
   pow_value = 1
   phys_window = .false.
   need_compute = .false.
   mask_variant = .false.
   reduced_k_range = .false.
   num_elems = 0
   num_threads = 1
   active_omp_level = 0
   issued_mask_ignore_warning = .false.
   mask = .true.


   call init_buff_values_1 (buff_obj%buffer, buff_obj%counter, buff_obj%count_0d, buff_obj%num_elements)

   hi = 0 !!halo size i
   hj = 0 !!halo size j
   l_start(1) = 1 !!local (to PE) start inddex
   l_start(2) = 1
   l_start(3) = 1
   l_end(1)  = SZ
   l_end(2) = SZ
   l_end(3) = SZ


   ALLOCATE( ofield_cfg )
   call ofield_cfg%initialize_for_ut(module_name1, field_name1, output_name1, pow_value, &
      & phys_window, need_compute, mask_variant, reduced_k_range , &
      & num_elems, time_reduction_type1, output_freq1 )
   ALLOCATE( ofield_index_cfg )
   CALL init_ofield_index_cfg(ofield_index_cfg, 1+hi, 1+hj, 1, SZ - hi, SZ - hj, SZ,&
      &  hi, hj, 1 + hi, SZ - hi, 1 + hj, SZ - hj)

   !!First make sure buffer vals are all zero
   call check_results_2(ofb_r4_ptr, 1, 0)

   !! Update the buffer values with the fieldbuff_update function.
   !! Case: mask_var=false & missval not present & mask not present & not_reduced_k_range
   test_passed = .true.  !! will be set to false if there are any issues.

   temp_result = fieldbuff_update(ofield_cfg, ofield_index_cfg, field_r4_ptr, sample, &
      & ofb_r4_ptr, ofc_r4_ptr, buff_bounds, &
      & ofb0d_r4_ptr (sample), buff_obj%num_elements(sample), &
      & mask, weight, missval_r4_ptr, &
      & num_threads, active_omp_level, &
      & issued_mask_ignore_warning, &
      & l_start, l_end, err_msg, err_msg_local )

   call check_results_1(ofb_r4_ptr, 1, "Buffer_update_test01")
   call print_output_field_values( buff_obj%buffer, 1 )

   !! ************ 2ND TEST: **********************
   !!First make sure buffer vals are all zero
   ofb_r4_ptr = 0
   call check_results_2(ofb_r4_ptr, 1, 0)

   !! Update the buffer values with the copy_fieldvals function.
   ! missvalue_present = .true. TBD
   !!call print_output_field_values( buff_obj%buffer, 1 )
   temp_result = fieldbuff_copy_fieldvals(ofield_cfg, ofield_index_cfg, field_r4_ptr, sample, &
      & ofb_r4_ptr, buff_bounds, &
      & ofb0d_r4_ptr(sample),  mask, missval_r4_ptr, &
      & l_start, l_end, err_msg, err_msg_local )

   !!call print_output_field_values(  buff_obj%buffer, 1 )

   call check_results_1(ofb_r4_ptr, 1, "Buffer_update_test02")

  call error_mesg('test_diag_update_buffer', 'Test has finished',NOTE)

   call fms_end

CONTAINS

   !> @brief Initialized an fms_diag_outfield_index_type by calling member function of
   !! fms_diag_outfield_index_type input object.
   SUBROUTINE init_ofield_index_cfg(idx_cfg, is, js , ks, ie, je, ke, hi, hj, f1, f2, f3, f4)
      type(fmsDiagOutfieldIndex_type), INTENT(inout)  :: idx_cfg !< The object to initialize.
      INTEGER, INTENT(in) :: is, js, ks  !< Var with same name in fms_diag_outfield_index_type
      INTEGER, INTENT(in) ::  ie, je, ke !< Var with same name in fms_diag_outfield_index_type
      INTEGER, INTENT(in) :: hi, hj !< Var with same name in fms_diag_outfield_index_type
      INTEGER, INTENT(in)  :: f1, f2, f3, f4 !< Var with same name in fms_diag_outfield_index_type
      call idx_cfg%initialize ( is, js , ks, ie, je, ke, hi, hj, f1, f2, f3, f4)
   end subroutine init_ofield_index_cfg

   SUBROUTINE init_field_values (field)
      CLASS(*), DIMENSION(:,:,:,:), INTENT(INOUT) :: field
      INTEGER :: NX,NY,NZ, NL
      INTEGER :: i,j,k,l
      INTEGER :: itemp
      NX = size(field,1)
      NY= size(field,2)
      NZ= size(field,3)
      NL= size(field,4)
      DO l = 1, NL
         DO k = 1, NZ
            DO j = 1, NY
               DO i = 1, NX
                  itemp = get_array_index_from_4D(i,j,k,l,NX,NY,NZ)
                  SELECT TYPE ( field)
                   TYPE IS (real(kind=r4_kind))
                     field(i,j,k,l) = real(itemp, kind=r4_kind)
                   TYPE IS (integer(kind=i8_kind))
                     field(i,j,k,l) = int(itemp, kind=i8_kind)
                  END SELECT
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE init_field_values

  !> @brief Init to zero the buffer, counter , an
   SUBROUTINE init_buff_values_1 (buffer, counter, count_0d, num_elems)
      CLASS(*), DIMENSION(:,:,:,:,:), INTENT(INOUT) :: buffer !< The actual buffer array of the buffer class.
      CLASS(*), DIMENSION(:,:,:,:,:), INTENT(INOUT) :: counter  !< The actual buffer array of the buffer class.
      CLASS(*), DIMENSION(:), INTENT(INOUT) :: count_0d  !< A counter used in time averaging.
      INTEGER, DIMENSION(:), INTENT(INOUT) :: num_elems  !< A counter used in time averaging.
      INTEGER, PARAMETER :: sample = 1  !< The diurnal sample.

      SELECT TYPE ( buffer)
       TYPE IS (real(kind=r4_kind))
         buffer  = 0
       TYPE IS (integer(kind=i8_kind))
         buffer = 0
      END SELECT

      SELECT TYPE ( counter)
       TYPE IS (real(kind=r4_kind))
         counter = 0
       TYPE IS (integer(kind=i8_kind))
         counter = 0
      END SELECT

      SELECT TYPE ( count_0d)
       TYPE IS (real(kind=r4_kind))
         count_0d  = 0
       TYPE IS (integer(kind=i8_kind))
         count_0d = 0
      end select

      num_elems = 0
   END SUBROUTINE init_buff_values_1


   SUBROUTINE print_output_field_values (buffer, onum)
      CLASS(*), ALLOCATABLE, DIMENSION(:,:,:,:,: ) :: buffer
      INTEGER, INTENT(IN) :: onum
      INTEGER :: i,j,k
      INTEGER :: ti
      REAL :: tr
      print *, "Start of print_output_field_values"
      k = 1
      DO j  =1 ,10
         DO i = 1,10
            SELECT TYPE ( buffer)
             TYPE IS (real(kind=r4_kind))
               !print "(10f10.1)", buffer(:,j,k,1,1)
               tr =  buffer(i,j,k,1,1)
               print "(f10.1)", tr
             TYPE IS (integer(kind=i8_kind))
               !print "(10I10)", buffer(:,j,k,1,1)
               !print "(I8))", buffer(i,j,k,1,1)
               print "(I8)", ti
            END SELECT
         end do
         print *, "************************"
      end do
      print *, "End of print_output_field_values"
   END SUBROUTINE print_output_field_values

!> @brief Verify that the buffer data is equal to the expected index value
   SUBROUTINE check_results_1(buff, sample, test_name)
      CLASS(*), DIMENSION(:,:,:,:,:), INTENT(IN) :: buff !< The 5D buffer
      INTEGER, INTENT(in) :: sample !< The diurnal sample
      CHARACTER(*), INTENT(in) :: test_name !< The test name
      INTEGER :: NX,NY,NZ, NL
      INTEGER :: i,j,k,l
      LOGICAL :: pass
      integer :: idx
      real :: bv
      pass = .true.
      NX = size(buff,1)
      NY= size(buff,2)
      NZ= size(buff,3)
      NL= size(buff,4)

      DO l = 1, NL
         DO k = 1, NZ
            DO j = 1, NY
               DO i = 1, NX
                  SELECT TYPE ( buff)
                   TYPE IS (real(kind=r4_kind))
                      idx =  get_array_index_from_4D(i,j,k,l,NX,NY,NZ)
                      bv = buff(i,j,k,l,sample)
                     if ( get_array_index_from_4D(i,j,k,l,NX,NY,NZ)  /= buff(i,j,k,l,sample) ) then
                        pass = .false.
                     endif
                   TYPE IS (integer(kind=i8_kind))
                     if ( get_array_index_from_4D(i,j,k,l,NX,NY,NZ)  /= buff(i,j,k,l,sample) ) then
                        pass = .false.
                     endif
                  END SELECT
               END DO
            END DO
         END DO
      END DO
      if ( pass .eqv. .false.) then
         call error_mesg('check_results_1', test_name//" has failed.",FATAL)
      end if
   end subroutine check_results_1

   SUBROUTINE check_results_2(buff, sample, val)
      CLASS(*), DIMENSION(:,:,:,:,:), INTENT(IN) :: buff
      INTEGER, INTENT(in) :: sample
      INTEGER, INTENT(in) :: val
      INTEGER :: NX,NY,NZ, NL
      INTEGER :: i,j,k,l
      LOGICAL :: pass

      pass = .true.
      NX = size(buff,1)
      NY= size(buff,2)
      NZ= size(buff,3)
      NL= size(buff,4)
      DO l = 1, NL
         DO k = 1, NZ
            DO j = 1, NY
               DO i = 1, NX
                  SELECT TYPE ( buff)
                   TYPE IS (real(kind=r4_kind))
                     if (  buff(i,j,k,l,sample) /= val ) then
                        pass = .false.
                     endif
                   TYPE IS (integer(kind=i8_kind))
                     if ( buff(i,j,k,l,sample)  /= val  ) then
                        pass = .false.
                     endif
                  END SELECT
               END DO
            END DO
         END DO
      END DO
      if ( pass .eqv. .false.) then
         call error_mesg('check_results_2', 'Test has failed',FATAL)
      end if
   end subroutine check_results_2

   !> @brief Calculate the unique index into a 4D array given the first four indecies
   !! i,j,k,l and the with in the fist three dimensions.
   pure integer function get_array_index_from_4D(i,j,k, l, NX,NY,NZ)
      INTEGER, INTENT(IN) :: i, j, k, l !< The three spatial dimentsions plus another
      INTEGER, INTENT(IN) :: NX, NY, NZ !< The size of the spatial dimentions.
      get_array_index_from_4D =  (l-1)* (NX * NY * NZ) + (k-1) * NX * NY + (j-1) * NX + i
   end function get_array_index_from_4D

   subroutine allocate_input_data_and_ptrs(datapoint, field_data, rmask, missvalue, mask, NX,NY,NZ, NL)
      CLASS(*), INTENT(in) :: datapoint !!The type of data we want
      CLASS(*), ALLOCATABLE, INTENT(inout), DIMENSION(:,:,:,:) :: field_data
      CLASS(*), ALLOCATABLE, INTENT(inout), DIMENSION(:,:,:,:) :: rmask
      CLASS(*), ALLOCATABLE,  INTENT(inout) :: missvalue
      LOGICAL, ALLOCATABLE, INTENT(inout), DIMENSION(:,:,:,:) :: mask
      INTEGER , INTENT(in) :: NX,NY,NZ, NL
      select type (datapoint)
       type is (integer(kind=i8_kind))
         allocate(integer(kind=i8_kind) :: field_data(NX,NY,NZ,NL))
         allocate(integer(kind=i8_kind) :: rmask(NX,NY,NZ,NL))
         allocate(integer(kind=i8_kind) :: missvalue)
       type is (real(kind=r4_kind))
         allocate(real(kind=r4_kind) :: field_data(NX,NY,NZ,NL))
         allocate(real(kind=r4_kind) :: rmask(NX,NY,NZ,NL))
         allocate(real(kind=r4_kind) :: missvalue)
       class default
         call error_mesg("allocate input data", "The input data type is not a  r4 or i8", FATAL)
      end select

      allocate(mask(NX,NY,NZ,NL))
   END  subroutine allocate_input_data_and_ptrs


   subroutine allocate_buffer_obj( data_point, bo, NX,NY,NZ, NL, NDI)
      TYPE(diagTestBuffer_type), INTENT(inout), allocatable :: bo
      CLASS(*), INTENT(in) :: data_point !< Sample point allocated to the type being tested.
      INTEGER, INTENT(IN) :: NX, NY, NZ !< The three spatial dimensions.
      INTEGER, INTENT(IN) :: NL !< Size of the 4th dimentions
      INTEGER, INTENT(IN) :: NDI !< Diurnal axis length,
      allocate (bo)
      select type (data_point)
       type is (integer(kind=i8_kind))
         allocate(integer(kind=i8_kind) :: bo%buffer(NX,NY,NZ,NL, NDI))
         allocate(integer(kind=i8_kind) :: bo%counter(NX,NY,NZ,NL, NDI))
         allocate(integer(kind=i8_kind) :: bo%count_0d(NDI))
       type is (real(kind=r4_kind))
         allocate(real(kind=r4_kind) :: bo%buffer(NX,NY,NZ,NL,NDI))
         allocate(real(kind=r4_kind) :: bo%counter(NX,NY,NZ,NL,NDI))
         allocate(real(kind=r4_kind) :: bo%count_0d(NDI))
       class default
         call error_mesg("allocate buffer obj", "The input data type is not a  r4 or i8", FATAL)
      end select

      allocate( bo%num_elements(NDI))

   END subroutine allocate_buffer_obj
end program test_diag_update_buffer


