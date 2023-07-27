!> \author Ganga Purja Pun
!> \email GFDL.Climate.Model.Info@noaa.gov
!! \brief Contains routines for the modern diag manager
!! These routines are meant to be used for checks and in reduction methods.

module fms_diag_reduction_methods_mod
  use platform_mod
  use mpp_mod, only: mpp_error
  use fms_mod, only: fms_error_handler
  use fms_diag_bbox_mod
  use fms_diag_output_buffer_mod
  use diag_data_mod, only: debug_diag_manager, time_max, time_min, time_average, time_sum
  use fms_diag_field_object_mod, only: fmsDiagField_type

  implicit none
  private

#ifdef use_yaml
  public :: compare_two_sets_of_bounds, real_copy_set, check_indices_order, init_mask_3d
  public :: fms_diag_update_extremum, update_scalar_extremum, update_array_extremum
  public :: fms_diag_time_average, fms_diag_mask_variant_do_sum, sum_scalar_field_data
  public :: fms_diag_no_mask_variant_do_sum, update_buffer_obj_num_elements, set_buffer_obj_count_0d
  public :: update_buffer_obj_count_0d
  contains

  !> @brief Compares the corresponding bounding indices of the first set with the second set.
  !> @return .TRUE. if any comparison returns true; i.e. the box bounded by the indices of the first set
  !! is out side the box bounded by the indices of the second set.
  LOGICAL FUNCTION compare_two_sets_of_bounds(bounds_a, bounds_b, error_str)
    integer, intent(in) :: bounds_a(:) !< First array with order: (/imin, imax, jmin, jmax, kmin, kmax/)
    integer, intent(in) :: bounds_b(:) !< Second array with the same order as the first
    character(*), intent(out) :: error_str !< Error message to report back

    compare_two_sets_of_bounds = .FALSE.

    if (size(bounds_a) .ne. size(bounds_b)) then
      compare_two_sets_of_bounds = .TRUE.
      error_str = 'fms_diag_reduction_methods_mod::compare_two_sets_of_bounds Error: sizes of sets do not match'
    else
      if ((size(bounds_a) .ne. 6) .and. (size(bounds_b) .ne. 6)) then
        compare_two_sets_of_bounds = .TRUE.
        error_str = 'fms_diag_reduction_methods_mod::compare_two_sets_of_bounds Error: sizes of sets must be 6'
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
  function check_indices_order(is_in, ie_in, js_in, je_in, error_msg) result(rslt)
    integer, intent(in), optional :: is_in, ie_in, js_in, je_in !< Indices passed to fms_diag_accept_data()
    character(len=*), intent(inout), optional :: error_msg !< An error message used only for testing purpose!!!

    character(len=52) :: err_module_name !< Stores the module name to be used in error calls
    logical :: rslt !< Return value

    rslt = .false. !< If no error occurs.

    err_module_name = 'fms_diag_reduction_methods_mod::check_indices_order'

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

  !> @brief Copies input data to output data with specific type and precision
  !! if the input data is present else sets the output data to a given value val if it is present.
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
        if (fms_error_handler('fms_diag_reduction_methods_mod::real_copy_set',&
          & 'The in_data is not one of the supported types of real(kind=4) or real(kind=8)', err_msg)) THEN
          return
        end if
      END SELECT
    ELSE
      if (present(val)) then
        out_data = val
      else
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::real_copy_set both in_data and val can be absent')
      end if
    END IF
  end subroutine real_copy_set

  !> @brief Allocates `outmask'(second argument) with sizes of the first three dimensions of
  !! the field(first argument).
  !! Initializes the `outmask' depending on presence/absence of `inmask' and `rmask'.
  !! Uses `rmask_threshold' to set the `outmask'.
  subroutine init_mask_3d(field, outmask, rmask_threshold, inmask, rmask, err_msg)
    class(*), intent(in) :: field(:,:,:,:)  !< Dummy variable whose sizes only in the first three
                                            !! dimensions are important
    logical, allocatable, intent(inout) :: outmask(:,:,:) !< Output logical mask
    class(*), intent(in) :: rmask_threshold !< Holds the values 0.5_r4_kind or 0.5_r8_kind, or related threhold values
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
      if (fms_error_handler('fms_diag_reduction_methods_mod::init_mask_3d', trim(err_msg_local), err_msg)) then
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
          select type (rmask_threshold)
          type is (real(kind=r4_kind))
            WHERE (rmask < rmask_threshold) outmask = .FALSE.
          class default
            call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::init_mask_3d'//&
              ' types of rmask and rmask_threshold do not match')
          end select
        TYPE IS (real(kind=r8_kind))
          select type (rmask_threshold)
          type is (real(kind=r8_kind))
            WHERE (rmask < rmask_threshold) outmask = .FALSE.
          class default
            call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::init_mask_3d'//&
              ' types of rmask and rmask_threshold do not match')
          end select
        CLASS DEFAULT
          call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::init_mask_3d'//&
            & ' The rmask is not one of the supported types of real(kind=4) or real(kind=8)')
      END SELECT
    END IF
  end subroutine init_mask_3d

  !> @brief Updates the buffer with the field data based on the value of the flag passed:
  !! time_min for minimum; time_max for maximum.
  subroutine fms_diag_update_extremum(flag, buffer_obj, field_data, recon_bounds, l_start, &
    l_end, is_regional, reduced_k_range, sample, mask, fieldName, hasDiurnalAxis, err_msg)
    integer, intent(in) :: flag !< Flag to indicate what to update: time_min for minimum; time_max for maximum
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    type(fmsDiagBoundsHalos_type), intent(inout) :: recon_bounds !< Indices of bounds in the first three dimension
                                                                 !! of the field data
    integer, intent(in) :: l_start(:) !< Local starting indices for the first three dimensions
    integer, intent(in) :: l_end(:)   !< Local ending indices for the first three dimensions
    logical, intent(in) :: is_regional !< Flag indicating if the current PE takes part in send_data
    logical, intent(in) :: reduced_k_range !< Flag indicating if the field has zbounds
    integer, intent(in) :: sample !< Index along the diurnal time axis
    logical, intent(in) :: mask(:,:,:,:) !< Must be out of range mask
    character(len=*), intent(in) :: fieldName !< Field name for error reporting
    logical, intent(in) :: hasDiurnalAxis !< Flag to indicate if the buffer has a diurnal axis
    character(len=*), intent(inout), optional :: err_msg !< Error mesage to report back

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    integer :: ksr, ker !< Reduced indices in the K dimension
    integer :: i, j, k !< For loops
    integer :: i1, j1, k1 !< Intermediate computed indices
    character(len=128) :: err_msg_local !< Stores local error message
    class(*), pointer :: ptr_buffer(:,:,:,:,:) !< Pointer to 5D buffer for remapping
    type(fmsDiagIbounds_type) :: IJKBounds !< Bounding object for the I, J, and K indices

    !> Get the `bounds3D` member of the `recon_bounds`
    IJKBounds = recon_bounds%get_bounds3D() !< Assignment of data structure with intrinsic type members may work!!!

    !> Unpack recon_bounds
    is = IJKBounds%get_imin()
    js = IJKBounds%get_jmin()
    ks = IJKBounds%get_kmin()
    ie = IJKBounds%get_imax()
    je = IJKBounds%get_jmax()
    ke = IJKBounds%get_kmax()
    hi = recon_bounds%get_hi()
    f1 = recon_bounds%get_fis()
    f2 = recon_bounds%get_fie()
    hj = recon_bounds%get_hj()
    f3 = recon_bounds%get_fjs()
    f4 = recon_bounds%get_fje()

    if (flag .ne. 3 .and. flag .ne. 4) then
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::fms_diag_update_extremum: flag must be either 3 or 4.")
    end if

    !! TODO: remap buffer before passing to subroutines update_scalar_extremum and update_array_extremum
    ptr_buffer => buffer_obj%remap_buffer(fieldName, hasDiurnalAxis)

    ! Update buffer
    regional_if: IF (is_regional) THEN
      DO k = l_start(3), l_end(3)
        k1 = k - l_start(3) + 1
        DO j = js, je
          DO i = is, ie
            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
              & j <= l_end(2)+hj ) THEN
              i1 = i-l_start(1)-hi+1
              j1=  j-l_start(2)-hj+1
              select type (buffer_obj)
              type is (outputBuffer0d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              type is (outputBuffer1d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              type is (outputBuffer2d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              type is (outputBuffer3d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              type is (outputBuffer4d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              type is (outputBuffer5d_type)
                call update_scalar_extremum(flag, field_data, ptr_buffer, mask, sample, &
                  recon_bounds, (/i,j,k/), (/i1,j1,k1/))
              class default
                call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
                  ' regional buffer_obj is not one of the support buffer types: outputBuffer0d_type'//&
                  ' outputBuffer1d_type outputBuffer2d_type outputBuffer3d_type'//&
                  ' outputBuffer4d_type outputBuffer5d_type')
              end select
            end if
          END DO
        END DO
      END DO
    ELSE !< if not regional
      reduced_k_range_if: IF (reduced_k_range) THEN
        call IJKBounds%set_kbounds(l_start(3), l_end(3))
        select type (buffer_obj)
        type is (outputBuffer0d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer1d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer2d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer3d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer4d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer5d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        class default
          call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum in reduced_k_range_if'//&
            ' regional buffer_obj is not one of the support buffer types: outputBuffer0d_type'//&
            ' outputBuffer1d_type outputBuffer2d_type outputBuffer3d_type'//&
            ' outputBuffer4d_type outputBuffer5d_type')
        end select
      ELSE !< does not have reduced_k_range
        debug_diag_if: IF ( debug_diag_manager ) THEN
          ! Compare bounds {is-hi, ie-hi, js-hj, je-hj, ks, ke} with the bounds of first three dimensions of the buffer
          if (compare_two_sets_of_bounds((/is-hi, ie-hi, js-hj, je-hj, ks, ke/), &
            (/LBOUND(ptr_buffer,1), UBOUND(ptr_buffer,1), LBOUND(ptr_buffer,2), UBOUND(ptr_buffer,2), &
            LBOUND(ptr_buffer,3), UBOUND(ptr_buffer,3)/), err_msg_local)) THEN
            IF ( fms_error_handler('fms_diag_object_mod::fms_diag_update_extremum', err_msg_local, err_msg) ) THEN
              RETURN
            END IF
          END IF
        END IF debug_diag_if

        !> If no error above, do update the buffer
        select type (buffer_obj)
        type is (outputBuffer0d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer1d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer2d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer3d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer4d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        type is (outputBuffer5d_type)
          call update_array_extremum(flag, field_data, ptr_buffer, mask, sample, recon_bounds, reduced_k_range)
        class default
          call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
            ' regional buffer_obj is not one of the support buffer types: outputBuffer0d_type'//&
            ' outputBuffer1d_type outputBuffer2d_type outputBuffer3d_type'//&
            ' outputBuffer4d_type outputBuffer5d_type')
        end select
      END IF reduced_k_range_if
    end if regional_if

    ! Reset counter count_0d of the buffer object
    call set_buffer_obj_count_0d(buffer_obj, 1.0)
  end subroutine fms_diag_update_extremum

  !> @brief Sets the counter `count_0d` of a buffer object to a given value
  subroutine set_buffer_obj_count_0d(buffer_obj, val)
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    class(*) intent(in) :: val !< Reset value

    select type (buffer_obj)
    type is (outputBuffer0d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer1d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer2d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer3d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer4d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer5d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    class default
      call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported buffer type')
    end select
  end subroutine set_buffer_obj_count_0d

  !> @brief Increments the counter `count_0d` of a buffer object by a given value
  subroutine update_buffer_obj_count_0d(buffer_obj, val)
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    class(*) intent(in) :: val !< Increment value

    select type (buffer_obj)
    type is (outputBuffer0d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer1d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer2d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer3d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer4d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    type is (outputBuffer5d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r4_kind)
      type is (real(kind=r8_kind))
        real_counter(sample) = real_counter(sample) + real(val, kind=r8_kind)
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum'//&
          ' Unsupported type of buffer_obj%count_0d')
      end select
    class default
      call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported buffer type')
    end select
  end subroutine update_buffer_obj_count_0d
 
  !> @brief Updates individual element of the buffer associated with indices in running_indx1 and running_indx2
  subroutine update_scalar_extremum(flag, field_data, buffer, mask, sample, recon_bounds, &
    running_indx1, running_indx2)
    integer, intent(in) :: flag !< Flag indicating maximum(time_max) or minimum(time_min)
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    class(*), intent(inout) :: buffer(:,:,:,:,:) !< Remapped output buffer
    logical, intent(in) :: mask(:,:,:,:) !< Update mask
    integer, intent(in) :: sample !< diurnal sample index
    type(fmsDiagBoundsHalos_type), intent(in) :: recon_bounds !< Holds starting and ending indices in the
                                                              !! I, J, and K dimensions and
                                                              !! halo sizes in the I, and J dimensions
    integer, intent(in) :: running_indx1(3) !< Holds indices i, j, and k
    integer, intent(in) :: running_indx2(3) !< Holds indices i1, j1, and k1

    type(fmsDiagIbounds_type) :: IJKBounds !< Bounding object for the I, J, and K indices
    integer :: i, j, k !< Unpack running_indx1 to
    integer :: i1, j1, k1 !< Unpack running_indx2 to
    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensiions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions

    !> Check flag for unsupported operation
    if (flag .ne. time_max .and. flag .ne. time_min) then
      call mpp_error(FATAL, "fms_diag_reduction_methods_mod::fms_diag_scalar_extremum &
        unsupported reduction method")
    endif

    ! Initialize i, j, and k
    i = running_indx1(1)
    j = running_indx1(2)
    k = running_indx1(3)

    ! Initialize i1, j1, and k1
    i1 = running_indx2(1)
    j1 = running_indx2(2)
    k1 = running_indx2(3)

    !> Get the `bounds3D` member of the `recon_bounds`
    IJKBounds = recon_bounds%get_bounds3D() !< Assignment of data structure with intrinsic type members may work!!!

    !> Unpack index bounds
      is = IJKBounds%get_imin()
      js = IJKBounds%get_jmin()
      ks = IJKBounds%get_kmin()
      ie = IJKBounds%get_imax()
      je = IJKBounds%get_jmax()
      ke = IJKBounds%get_kmax()
      hi = recon_bounds%get_hi()
      hj = recon_bounds%get_hj()

    ! Select proper type and update the buffer
    select type (field_data)
    type is (real(kind=r4_kind))
      select type (buffer)
      type is (real(kind=r4_kind))
        minimum_if: if (flag .eq. time_min) then
        !> Update the buffer with the current minimum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) <&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        else !< if not minimum, check for maximum
        !> Update the buffer with the current maximum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) >&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (real(kind=r8_kind))
      select type (buffer)
      type is (real(kind=r8_kind))
        minimum_if: if (flag .eq. time_min) then
        ! Update the buffer with the current minimum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) <&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        else
        ! Update the buffer with the current maximum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) >&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        endif minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (integer(kind=i4_kind))
      select type (buffer)
      type is (integer(kind=i4_kind))
        minimum_if: if (flag .eq. time_min) then
        ! Update the buffer with the current minimum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) <&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        else
        ! Update the buffer with the current maximum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) >&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        endif minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (integer(kind=i8_kind))
      select type (buffer)
      type is (integer(kind=i8_kind))
        minimum_if: if (flag .eq. time_min) then
        ! Update the buffer with the current minimum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) <&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        else
        ! Update the buffer with the current maximum
          where (mask(i-is+1+hi,j-js+1+hj,k,:) .AND. field_data(i-is+1+hi,j-js+1+hj,k,:) >&
            buffer(i1,j1,k1,:,sample))
            buffer(i1,j1,k1,:,sample) = field_data(i-is+1+hi,j-js+1+hj,k,:)
          end where
        end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    class default
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum unsupported field data type")
    end select
  end subroutine update_scalar_extremum

  !> @brief Updates a chunk of the buffer defined by the bounds in recon_bounds
  subroutine update_array_extremum(flag, field_data, buffer, mask, sample, recon_bounds, reduced_k_range)
    integer :: flag !< Flag indicating maximum(time_max) or minimum(time_min)
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    class(*), intent(inout) :: buffer(:,:,:,:,:) !< Remapped output buffer
    logical, intent(in) :: mask(:,:,:,:) !< Updated mask
    integer, intent(in) :: sample !< diurnal sample index
    type(fmsDiagBoundsHalos_type), intent(in) :: recon_bounds !< Object to hold starting and ending indices
                                                              !! in the I, J, and K dimensions; also holds
                                                              !! halo sizes in the I, and J dimensions
    logical, intent(in) :: reduced_k_range !< Flag indicating if the range in the K dimension is present

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    type(fmsDiagIbounds_type) :: IJKBounds !< Bounding object for the I, J, and K indices

    !> Check flag for unsupported operation
    if (flag .ne. time_max .and. flag .ne. time_min) then
      call mpp_error(FATAL, "fms_diag_reduction_methods_mod::fms_diag_scalar_extremum &
        unsupported reduction method")
    endif

    !> Get the `bounds3D` member of the `recon_bounds`
    IJKBounds = recon_bounds%get_bounds3D() !< Assignment of data structure with intrinsic type members may work!!!

    !> Unpack bounds (/is, js, ks, ie, je, ke, hi, f1, f2, hj, f3, f4/)
    is = IJKBounds%get_imin()
    js = IJKBounds%get_jmin()
    ks = IJKBounds%get_kmin()
    ie = IJKBounds%get_imax()
    je = IJKBounds%get_jmax()
    ke = IJKBounds%get_kmax()
    hi = recon_bounds%get_hi()
    f1 = recon_bounds%get_fis()
    f2 = recon_bounds%get_fie()
    hj = recon_bounds%get_hj()
    f3 = recon_bounds%get_fjs()
    f4 = recon_bounds%get_fje()

    ! Select proper type and update the buffer
    select type (field_data)
    type is (real(kind=r4_kind))
      select type (buffer)
      type is (real(kind=r4_kind))
        minimum_if: if (flag .eq. time_min) then
        !> Update the buffer with the current minimum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        else
        !> Update the buffer with the current maximum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:)>&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (real(kind=r8_kind))
      select type (buffer)
      type is (real(kind=r8_kind))
          minimum_if: if (flag .eq. time_min) then
          !> Update the buffer with the current minimum
            if (reduced_k_range) then
              ! recon_bounds must have ks = ksr and ke = ker
              WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
                buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
                buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
            else
              WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
                buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
                buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
            end if
          else
          !> Update the buffer with the current maximum
            if (reduced_k_range) then
              ! recon_bounds must have ks = ksr and ke = ker
              WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
                buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
                buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
            else
              WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:)>&
                buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
                buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
            end if
          end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (integer(kind=i4_kind))
      select type (buffer)
      type is (integer(kind=i4_kind))
        minimum_if: if (flag .eq. time_min) then
        !> Update the buffer with the current minimum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        else
        !> Update the buffer with the current maximum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:)>&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    type is (integer(kind=i8_kind))
      select type (buffer)
      type is (integer(kind=i8_kind))
        minimum_if: if (flag .eq. time_min) then
        !> Update the buffer with the current minimum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        else
        !> Update the buffer with the current maximum
          if (reduced_k_range) then
            ! recon_bounds must have ks = ksr and ke = ker
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:) <&
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          else
            WHERE (mask(f1:f2,f3:f4,ks:ke,:) .AND. field_data(f1:f2,f3:f4,ks:ke,:)>&
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample)) &
              buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field_data(f1:f2,f3:f4,ks:ke,:)
          end if
        end if minimum_if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum"//&
          " buffer type does not match with field_data type.")
      end select
    class default
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum unsupported field data type")
    end select
  end subroutine update_array_extremum

  !> @brief Performs actual computation for sum or average for mask variant field depending
  !! on the flag value passed: time_sum for sum; time_average for average
  subroutine fms_diag_mask_variant_do_sum(flag, field_data, buffer_obj, recon_bounds, kr_start, kr_end, &
    reduced_k_range, pow_val, mask, weight, hasDiurnalAxis, field_name, sample)
    integer, intent(in) :: flag !< Flag indicating a reduction method: sum(time_sum) or average(time_average)
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    type(fmsDiagBoundsHalos_type), intent(inout) :: recon_bounds !< Indices of bounds in the first three dimensions
                                                                 !! of the field data
    integer, intent(in) :: kr_start(:) !< Local starting indices for the first three dimensions
    integer, intent(in) :: kr_end(:)   !< Local ending indices for the first three dimensions
    logical, intent(in) :: reduced_k_range !< Flag indicating if the field has zbounds
    integer, intent(in) :: pow_val !< Power value as an exponent
    logical, intent(in) :: mask(:,:,:,:) !< Must be out of range mask
    real, intent(in)    :: weight !< Must be a updated weight
    logical, intent(in) :: hasDiurnalAxis !< Flag to indicate if the buffer has a diurnal axis
    character(len=*), intent(in) :: field_name !< Field name for error reporting
    integer, intent(in) :: sample !< Index along the diurnal time axis

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    integer :: ksr, ker !< Reduced indices in the K dimension
    integer :: i, j, k !< Running indices for looping
    integer :: k1 !< Secondary running index
    class(*), pointer :: ptr_buffer(:,:,:,:,:) !< Pointer to 5D buffer for remapping
    character(len=52) :: code_module_name !< String to hold module and routine names
    type(fmsDiagIbounds_type) :: IJKBounds !< Index bounds in the I, J, and K dimensions

    !> Check flag for time_sum and time_average
    if (flag /= time_sum .and. flag /= time_average) then
      call mpp_error( FATAL, TRIM(str_modname)//" flag must be a parameter either time_sum or time_average.")
    endif

    !> Initialize code_module_name
    code_module_name = "fms_diag_reduction_methods_mod::fms_diag_mask_do_sum"

    !> Remap buffer
    ptr_buffer => buffer_obj%remap_buffer(fieldName, hasDiurnalAxis)

    !> Unpack bounds to individual indices
    IJKBounds = recon_bounds%get_bounds3D()
    is = IJKBounds%get_imin()
    js = IJKBounds%get_jmin()
    ks = IJKBounds%get_kmin()
    ie = IJKBounds%get_imax()
    je = IJKBounds%get_jmax()
    ke = IJKBounds%get_kmax()
    hi = recon_bounds%get_hi()
    hj = recon_bounds%get_hj()

    if (reduced_k_range) then
      ksr = kr_start(3)
      ker = kr_end(3)
    else
      ksr = ks
      ker = ke
    end if

    DO k= ksr, ker
      if (reduced_k_range) then
        k1= k - ksr + 1
      else
        k1 = k
      end if
      DO j=js, je
         DO i=is, ie
            IF (mask(i-is+1+hi, j-js+1+hj, k, :)) THEN
              select type (buffer_obj)
                type is (outputBuffer0d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                type is (outputBuffer1d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                type is (outputBuffer2d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                type is (outputBuffer3d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                type is (outputBuffer4d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                type is (outputBuffer5d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i-hi, j-hj, k1, 1/), weight, pow_val, sample)
                  if (flag == time_average) then
                    select type (counter_type => buffer_obj%counter)
                      type is (real(kind=r4_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      type is (real(kind=r8_kind))
                        counter_type(i-hi, j-hj, k1, :, sample) = counter_type(i-hi, j-hj, k1, :, sample) + weight
                      class default
                        call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported type and must be either &
                          r4 or r8 type')
                    end select
                  endif
                class default
                  call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported dimensional buffer type')
              end select
            END IF
         END DO
      END DO
    END DO
  end subroutine fms_diag_mask_variant_do_sum

  !> @brief Updates the counter `num_elements` of buffer object
  !! which is polymorphic in dimensionality
  subroutine update_buffer_obj_num_elements(buffer_obj, incr, sample)
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    integer, intent(in) :: incr !< Value to increment by
    integer, intent(in) :: sample !< Index along the diurnal time axis

    select type (buffer_obj)
      type is (outputBuffer0d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      type is (outputBuffer1d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      type is (outputBuffer2d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      type is (outputBuffer3d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      type is (outputBuffer4d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      type is (outputBuffer5d_type)
        buffer_obj%num_elements(sample) = buffer_obj%num_elements(sample) + incr
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::update_buffer_num_elemets &
          Unsupported dimensional buffer type')
    end select
  end subroutine update_buffer_obj_num_elements

  !> @brief Computes sum for the non mask variant case in time averaging.
  subroutine fms_diag_no_mask_variant_do_sum(field_data, buffer_obj, recon_bounds, l_start, l_end, &
    pow_val, mask, weight, field_name, hasDiurnalAxis, sample, missing_value)
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    type(fmsDiagBoundsHalos_type), intent(inout) :: recon_bounds !< Indices of bounds in the first three dimensions
                                                                 !! of the field data
    integer, intent(in) :: l_start(:) !< Local starting indices for the first three dimensions
    integer, intent(in) :: l_end(:)   !< Local ending indices for the first three dimensions
    integer, intent(in) :: pow_val !< Power value as an exponent
    logical, intent(in) :: mask(:,:,:,:) !< Must be out of range mask
    real, intent(in)    :: weight !< Must be a updated weight
    character(len=*), intent(in) :: field_name !< Field name for error reporting
    logical, intent(in) :: hasDiurnalAxis !< Flag to indicate if the buffer has a diurnal axis
    integer, intent(in) :: sample !< Index along the diurnal time axis
    class(*), intent(in) :: missing_value !< Missing value of the field data

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    integer :: ksr, ker !< Reduced indices in the K dimension
    integer :: i, j, k !< Running indices for looping
    integer :: i1, j1, k1 !< Secondary running indices
    class(*), pointer :: ptr_buffer(:,:,:,:,:) !< Pointer to 5D buffer for remapping
    type(fmsDiagIbounds_type) :: IJKBounds !< Index bounds in the I, J, and K dimensions
    character(len=64) :: code_module_name !< Holds module and routine names for error reporting

    code_module_name = "fms_diag_reduction_methods_mod::fms_diag_no_mask_variant_do_sum"

    !> Remap buffer
    ptr_buffer => buffer_obj%remap_buffer(field_name, hasDiurnalAxis)

    !> Unpack bounds to individual indices
    IJKBounds = recon_bounds%get_bounds3D()
    is = IJKBounds%get_imin()
    js = IJKBounds%get_jmin()
    ks = IJKBounds%get_kmin()
    ie = IJKBounds%get_imax()
    je = IJKBounds%get_jmax()
    ke = IJKBounds%get_kmax()
    hi = recon_bounds%get_hi()
    hj = recon_bounds%get_hj()

    DO k = l_start(3), l_end(3)
      k1 = k-l_start(3)+1
      DO j = js, je
        DO i = is, ie
          IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
            & j <= l_end(2)+hj ) THEN
            i1 = i-l_start(1)-hi+1
            j1=  j-l_start(2)-hj+1
            IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
              select type (buffer_obj)
                type is (outputBuffer0d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                type is (outputBuffer1d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                type is (outputBuffer2d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                type is (outputBuffer3d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                type is (outputBuffer4d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                type is (outputBuffer5d_type)
                  call sum_scalar_field_data(field_data, (/i-is+1+hi, j-js+1+hj, k, 1/), &
                    ptr_buffer, (/i1, j1, k1, 1/), weight, pow_val, sample)
                class default
                  call mpp_error(FATAL, TRIM(code_module_name)//' Unsupported dimensional buffer type')
              end select
            ELSE
              select type (buffer_obj) !< Select dimensional buffer type
                type is (outputBuffer0d_type) !< Scalar buffer
                  select type (ptr_buffer) !< Select type of scalar buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of scalar buffer type
                type is (outputBuffer1d_type) !< 1D buffer
                  select type (ptr_buffer) !< Select type of 1D buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of 1D buffer type
                type is (outputBuffer2d_type) !< 2D buffer
                  select type (ptr_buffer) !< Select type of 2D buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of 2D buffer type
                type is (outputBuffer3d_type) !< 3D buffer
                  select type (ptr_buffer) !< Select type of 3D buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of 3D buffer type
                type is (outputBuffer4d_type) !< 4D buffer
                  select type (ptr_buffer) !< Select type of 4D buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of 4D buffer type
                type is (outputBuffer5d_type) !< 5D buffer
                  select type (ptr_buffer) !< Select type of 5D buffer
                    type is (real(kind=r4_kind))
                      select type (missing_value)
                        type is (real(kind=r4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r4')
                      end select
                    type is (real(kind=r8_kind))
                      select type (missing_value)
                        type is (real(kind=r8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be r8')
                      end select
                    type is (integer(kind=i4_kind))
                      select type (missing_value)
                        type is (integer(kind=i4_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i4')
                      end select
                    type is (integer(kind=i8_kind))
                      select type (missing_value)
                        type is (integer(kind=i8_kind))
                          ptr_buffer(i1, j1, k1, 1, sample) = missing_value
                        class default
                          call mpp_error(FATAL, TRIM(code_module_name)//&
                            ' Unsupported type in missing value; must be i8')
                      end select
                    class default
                      call mpp_error(FATAL, TRIM(code_module_name)//&
                        ' Unsupported buffer type; must be either i4, i8, r4 or r8')
                  end select !< End of selection of 5D buffer type
                class default
                  call mpp_error(FATAL, TRIM(code_module_name)//&
                    ' Unsupported dimensional buffer type')
              end select !< End of selection of dimensional buffer
            END IF
          END IF
        END DO
      END DO
    END DO
  end subroutine

  !> @brief Adds single field data to output buffer element.
  !! Order of indices in field_indices and buf_indices must be in (/i, j, k, l/)
  !! The argument `weight` is not involved in any calculations related to buffer type i4 or i8
  subroutine sum_scalar_field_data(field_data, field_indices, buffer, buf_indices, weight, power_val, sample)
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    integer, intent(in) :: field_indices(4) !< Indices of field data
    class(*), intent(inout) :: buffer(:,:,:,:,:) !< Remapped output buffer
    integer, intent(in) :: buf_indices(4) !< Indices of buffer data
    real, intent(in)    :: weight !< Must be a updated weight
    integer, intent(in) :: power_val !< Power value as an exponent
    integer, intent(in) :: sample !< Index along the diurnal time axis

    character(len=48) :: mod_info !< Holds code module information

    mod_info = 'fms_diag_reduction_methods_mod::sum_scalar_field_data'
    
    !> Add field data to buffer; if `weight` is not equal to 1.0, it is a weighted sum/average.
    select type (field_data)
    type is (real(kind=r4_kind))
      select type (buffer)
      type is (real(kind=r4_kind))
        select case (power_val)
        case (1)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight
        case (2)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight *&
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight
        case default
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            (field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight)**(power_val)
        end select
      class default
         call mpp_error(FATAL, TRIM(mod_info)//' Buffer type not supported and must be either r4, r8, i4 or i8')
      end select
    type is (real(kind=r8_kind))
      select type (buffer)
      type is (real(kind=r8_kind))
        select case (power_val)
        case (1)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight
        case (2)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight *&
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight
        case default
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            (field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) * weight)**(power_val)
        end select
      class default
        call mpp_error(FATAL, TRIM(mod_info)//' Buffer type not supported and must be either r4, r8, i4 or i8')
      end select
    type is (integer(kind=i4_kind))
      select type (buffer)
      type is (integer(kind=i4_kind))
        select case (power_val)
        case (1)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4))
        case (2)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) *&
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4))
        case default
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            (field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)))**(power_val)
        end select
      class default
        call mpp_error(FATAL, TRIM(mod_info)//' Buffer type not supported and must be either r4, r8, i4 or i8')
      end select
    type is (integer(kind=i8_kind))
      select type (buffer)
      type is (integer(kind=i8_kind))
        select case (power_val)
        case (1)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4))
        case (2)
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)) *&
            field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4))
        case default
          buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) = &
            buffer(buf_indices(1), buf_indices(2), buf_indices(3), buf_indices(4), sample) + &
            (field_data(field_indices(1), field_indices(2), field_indices(3), field_indices(4)))**(power_val)
        end select
      class default
        call mpp_error(FATAL, TRIM(mod_info)//' Buffer type not supported and must be either r4, r8, i4 or i8')
      end select
    class default
      call mpp_error(FATAL, TRIM(mod_info)//' Unsupported field data type and must be either r4, r8, i4 or i8')
    end select
  end subroutine sum_scalar_field_data

  !> @brief A wrapper to perform time average or sum of the input field data
  !! depending on the flag value passed: time_sum for sum; time_average for average
  subroutine fms_diag_time_average(flag, field, buffer_obj, field_data, recon_bounds, l_start, &
    l_end, is_regional, reduced_k_range, sample, mask, fieldName, has_diurnal_axis, phys_win, &
    weight, pow_val, err_msg)
    integer, intent(in) :: flag !< Flag indicating reduction method: sum(time_sum) or average(time_average)
    type(fmsDiagField_type), intent(in) :: field !< Field object
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Buffer object
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    type(fmsDiagBoundsHalos_type), intent(inout) :: recon_bounds !< Holds indices of bounds in the I, J, and K
                                                                 !! dimensions and halo sizes in the I and
                                                                 !! J dimensions
    integer, intent(in) :: l_start(:) !< Local starting indices for the first three dimensions
    integer, intent(in) :: l_end(:)   !< Local ending indices for the first three dimensions
    logical, intent(in) :: is_regional !< Flag indicating if the current PE takes part in send_data
    logical, intent(in) :: reduced_k_range !< Flag indicating if the field has zbounds
    integer, intent(in) :: sample !< Index along the diurnal time axis
    logical, intent(in) :: mask(:,:,:,:) !< Must be out of range mask
    character(len=*), intent(in) :: fieldName !< Field name for error reporting
    logical, intent(in) :: has_diurnal_axis !< Flag to indicate if the buffer has a diurnal axis
    logical, intent(in) :: phys_win !< Flag indicating if the field is a physics window
    real, intent(in)    :: weight !< Must be a updated weight
    integer, intent(in) :: pow_val !< Power value as an exponent
    character(len=*), intent(inout), optional :: err_msg !< Error message

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    integer :: ksr, ker !< Reduced indices in the K dimension
    integer :: i, j, k, i1, j1, k1 !< For loops
    class(*), pointer :: ptr_buffer(:,:,:,:,:) !< Pointer to 5D buffer for remapping
    character(len=128) :: err_msg_local !< Stores local error message
    character(len=128) :: error_string !< Holds partial error text
    character(len=52) :: str_modname !< Holds names of this subroutine and the module it is in
    type(fmsDiagIbounds_type) :: IJKBounds !< Index bounds in the I, J, and K dimensions

    !> Initialize the local error reporting strings
    error_string = ''
    str_modname = "fms_diag_reduction_methods_mod::fms_diag_update_sum"

    !> Unpack recon_bounds to respective indices
    IJKBounds = recon_bounds%get_bounds3D()
    is = IJKBounds%get_imin()
    js = IJKBounds%get_jmin()
    ks = IJKBounds%get_kmin()
    ie = IJKBounds%get_imax()
    je = IJKBounds%get_jmax()
    ke = IJKBounds%get_kmax()
    hi = recon_bounds%get_hi()
    f1 = recon_bounds%get_fis()
    f2 = recon_bounds%get_fie()
    hj = recon_bounds%get_hj()
    f3 = recon_bounds%get_fjs()
    f4 = recon_bounds%get_fje()

    if (flag /= time_average .and. flag /= time_sum) then
      call mpp_error( FATAL, TRIM(str_modname)//" flag must be a parameter either time_sum or time_average.")
    end if

    IF (field%is_mask_variant()) THEN
      IF (is_regional) THEN
         WRITE (error_string,'(a,"/",a)')  &
              & TRIM(field%get_modname()), &
              & TRIM(fieldName)
         IF (fms_error_handler(TRIM(str_modname), 'module/field_name '//TRIM(error_string)//&
              & ', regional output NOT supported with mask_variant', err_msg)) THEN
            RETURN
         END IF
      END IF

      ! Should reduced_k_range data be supported with the mask_variant option   ?????
      ! If not, error message should be produced and the reduced_k_range loop below eliminated
      IF (PRESENT(mask)) THEN
         IF (field%has_missing_value()) THEN
            IF (debug_diag_manager) THEN
              ! Compare bounds {is-hi, ie-hi, js-hj, je-hj, ks, ke} with the bounds of first three dimensions of the buffer
              if (compare_two_sets_of_bounds((/is-hi, ie-hi, js-hj, je-hj, ks, ke/), &
                (/LBOUND(ptr_buffer,1), UBOUND(ptr_buffer,1), LBOUND(ptr_buffer,2), UBOUND(ptr_buffer,2), &
                LBOUND(ptr_buffer,3), UBOUND(ptr_buffer,3)/), err_msg_local)) THEN
                IF (fms_error_handler(TRIM(str_modname), err_msg_local, err_msg)) THEN
                  RETURN
                END IF
              END IF
            END IF
            IF(phys_win) then
              call fms_diag_mask_variant_do_sum(flag, field_data, buffer_obj, recon_bounds, l_start, l_end, &
                reduced_k_range, pow_val, mask, weight, has_diurnal_axis, fieldName, sample)
            ELSE
              call fms_diag_mask_variant_do_sum(flag, field_data, buffer_obj, recon_bounds, l_start, l_end, &
                reduced_k_range, pow_val, mask, weight, has_diurnal_axis, fieldName, sample)
            END IF
         ELSE
            WRITE (error_string,'(a,"/",a)')&
                 & TRIM(field%get_modname()), &
                 & TRIM(fieldName)
            IF(fms_error_handler(TRIM(str_modname), 'module/field_name '//TRIM(error_string)//&
                 & ', variable mask but no missing value defined', err_msg)) THEN
               RETURN
            END IF
         END IF
      ELSE  ! no mask present
         WRITE (error_string,'(a,"/",a)')&
              & TRIM(field%get_modname()), &
              & TRIM(fieldName)
         IF(fms_error_handler(TRIM(str_modname),'module/field_name'//TRIM(error_string)//&
              & ', variable mask but no mask given', err_msg)) THEN
            RETURN
         END IF
      END IF
   ELSE ! mask_variant=false
      IF ( PRESENT(mask) ) THEN
         IF (field%has_missing_value()) THEN
            IF (is_regional) THEN
               IF (phys_window) then
                 call fms_diag_no_mask_variant_do_sum(field_data, buffer_obj, recon_bounds, &
                   l_start, l_end, pow_val, mask, weight, fieldName, has_diurnal_axis, &
                   sample, field%get_missing_value())
               ELSE
                 call fms_diag_no_mask_variant_do_sum(field_data, buffer_obj, recon_bounds, &
                   l_start, l_end, pow_val, mask, weight, fieldName, has_diurnal_axis, &
                   sample, field%get_missing_value())
               ENDIF
               DO j = js, je
                  DO i = is, ie
                     IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                        call update_buffer_obj_num_elements(buffer_obj, (l_end(3)-l_start(3)+1), sample)
                     END IF
                  END DO
               END DO
!======================================================================================================
               !!TODO: everything below has not been refactored/updated
            ELSE IF ( reduced_k_range ) THEN
               IF (numthreads>1 .AND. phys_window) then
                  DO k=ksr, ker
                     k1 = k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                           IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ksr, ker
                     k1 = k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                           IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        DEALLOCATE(field_out)
                        DEALLOCATE(oor_mask)
                        RETURN
                     END IF
                  END IF
               END IF
               IF (numthreads>1 .AND. phys_window) then
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            END IF
            IF ( need_compute .AND. .NOT.phys_window ) THEN
               IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3))) ) &
                    & output_fields(out_num)%count_0d(sample) =&
                    & output_fields(out_num)%count_0d(sample) + weight1
            ELSE
               IF ( ANY(mask(f1:f2,f3:f4,ks:ke)) ) output_fields(out_num)%count_0d(sample) =&
                    & output_fields(out_num)%count_0d(sample)+weight1
            END IF
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
                           IF ( pow_value /= 1 ) THEN
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                   & (field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                           ELSE
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                   & field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           END IF
                        END IF
                     END DO
                  END DO
               ELSE
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                           i1 = i-l_start(1)-hi+1
                           j1 =  j-l_start(2)-hj+1
                           IF ( pow_value /= 1 ) THEN
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                   & (field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                           ELSE
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                   & field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           END IF
                        END IF
                     END DO
                  END DO
               END IF
               DO j = js, je
                  DO i = is, ie
                     IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                        output_fields(out_num)%num_elements(sample)=&
                             & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1

                     END IF
                  END DO
               END DO
            ELSE IF ( reduced_k_range ) THEN
               IF (numthreads>1 .AND. phys_window) then
                  ksr= l_start(3)
                  ker= l_end(3)
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                          & (field_out(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                          & field_out(f1:f2,f3:f4,ksr:ker)*weight1
                  END IF
               ELSE
                  ksr= l_start(3)
                  ker= l_end(3)
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                          & (field_out(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                          & field_out(f1:f2,f3:f4,ksr:ker)*weight1
                  END IF
               END IF
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '') THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        DEALLOCATE(field_out)
                        DEALLOCATE(oor_mask)
                        RETURN
                     END IF
                  END IF
               END IF
               IF (numthreads>1 .AND. phys_window) then
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & (field_out(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & field_out(f1:f2,f3:f4,ks:ke)*weight1
                  END IF
               ELSE
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & (field_out(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & field_out(f1:f2,f3:f4,ks:ke)*weight1
                  END IF
               END IF
            END IF
            IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                 & output_fields(out_num)%count_0d(sample) + weight1
         END IF
      ELSE ! mask NOT present
         IF ( missvalue_present ) THEN
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
                              IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                 IF ( pow_value /= 1 ) THEN
                                    output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                         & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                         & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                 ELSE
                                    output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                         & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                         & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                                 END IF
                              ELSE
                                 output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                 IF ( pow_value /= 1 ) THEN
                                    output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                         & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                         & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                 ELSE
                                    output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                         & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                         & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                                 END IF
                              ELSE
                                 output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
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
                           IF ( field_out(i,j,k) /= missvalue ) THEN
                              output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample)&
                                                                      & + weight1
                              EXIT outer0
                           END IF
                        END DO
                     END DO
                  END DO outer0
               END IF
            ELSE IF ( reduced_k_range ) THEN
               if( numthreads>1 .AND. phys_window ) then
                  ksr= l_start(3)
                  ker= l_end(3)
                  DO k = ksr, ker
                     k1 = k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                           IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               else
                  ksr= l_start(3)
                  ker= l_end(3)
                  DO k = ksr, ker
                     k1 = k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                           IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
               outer3: DO k = ksr, ker
                  k1=k-ksr+1
                  DO j=f3, f4
                     DO i=f1, f2
                        IF ( field_out(i,j,k) /= missvalue ) THEN
                           output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) &
                                                                   & + weight1
                           EXIT outer3
                        END IF
                     END DO
                  END DO
               END DO outer3
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        DEALLOCATE(field_out)
                        DEALLOCATE(oor_mask)
                        RETURN
                     END IF
                  END IF
               END IF
               IF( numthreads > 1 .AND. phys_window ) then
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           IF ( field_out(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                              IF ( pow_value /= 1 ) THEN
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & (field_out(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                              ELSE
                                 output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                      & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                      & field_out(i-is+1+hi,j-js+1+hj,k) * weight1
                              END IF
                           ELSE
                              output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
               outer1: DO k=ks, ke
                  DO j=f3, f4
                     DO i=f1, f2
                        IF ( field_out(i,j,k) /= missvalue ) THEN
                           output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) &
                                                                   & + weight1
                           EXIT outer1
                        END IF
                     END DO
                  END DO
               END DO outer1
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
                           IF ( pow_value /= 1 ) THEN
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                   & (field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                           ELSE
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                   & field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           END IF
                        END IF
                     END DO
                  END DO
               ELSE
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                           i1 = i-l_start(1)-hi+1
                           j1=  j-l_start(2)-hj+1
                           IF ( pow_value /= 1 ) THEN
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                   & (field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                           ELSE
                              output_fields(out_num)%buffer(i1,j1,:,sample)= &
                                   & output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                   & field_out(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           END IF
                        END IF
                     END DO
                  END DO
               END IF

               DO j = js, je
                  DO i = is, ie
                     IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                        & j <= l_end(2)+hj ) THEN
                        output_fields(out_num)%num_elements(sample) =&
                             & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                     END IF
                  END DO
               END DO
               ! Accumulate time average
            ELSE IF ( reduced_k_range ) THEN
               ksr= l_start(3)
               ker= l_end(3)
               IF( numthreads > 1 .AND. phys_window ) then
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                          & (field_out(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                          & field_out(f1:f2,f3:f4,ksr:ker)*weight1
                  END IF
               ELSE
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                          & (field_out(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                          & field_out(f1:f2,f3:f4,ksr:ker)*weight1
                  END IF
               END IF
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        DEALLOCATE(field_out)
                        DEALLOCATE(oor_mask)
                        RETURN
                     END IF
                  END IF
               END IF
               IF( numthreads > 1 .AND. phys_window ) then
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & (field_out(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & field_out(f1:f2,f3:f4,ks:ke)*weight1
                  END IF
               ELSE
                  IF ( pow_value /= 1 ) THEN
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & (field_out(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                  ELSE
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                          & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                          & field_out(f1:f2,f3:f4,ks:ke)*weight1
                  END IF
               END IF
            END IF
            IF ( .NOT.phys_win ) call update_buffer_obj_count_0d(buffer_obj, weight)
         END IF
      END IF ! if mask present
   END IF  !if mask_variant
   IF ( .NOT.is_regional .AND. .NOT.reduced_k_range )&
     call update_buffer_obj_num_elements(buffer_obj, (ie-is+1)*(je-js+1)*(ke-ks+1), sample)
   IF ( reduced_k_range ) &
     call update_buffer_obj_num_elements(buffer_obj, (ie-is+1)*(je-js+1)*(ker-ksr+1), sample)
        

  end subroutine fms_diag_time_average
#endif
end module fms_diag_reduction_methods_mod
