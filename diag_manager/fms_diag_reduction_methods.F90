!> \author Ganga Purja Pun
!> \email gagna.purjapun@noaa.gov
!! \brief Contains routines for the modern diag_manager
!!
!! \description

module fms_diag_reduction_methods_mod
  use platform_mod
  use mpp_mod, only: mpp_error
  use fms_mod, only: fms_error_handler
  use fms_diag_bbox_mod
  use fms_diag_output_buffer_mod
  use diag_data_mod, only: debug_diag_manager

  implicit none
  private

#ifdef use_yaml
  public :: compare_two_sets_of_bounds, real_copy_set, check_indices_order, init_mask_3d
  public :: fms_diag_update_extremum, update_scalar_extremum, update_array_extremum
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
        if (fms_error_handler('fms_diag_reduction_methods_mod::real_copy_set',&
          & 'The in_data is not one of the supported types of real(kind=4) or real(kind=8)', err_msg)) THEN
          return
        end if
      END SELECT
    ELSE
      if (present(val)) out_data = val
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
            if (fms_error_handler('fms_diag_reduction_methods_mod::init_mask_3d', 'type mismatch', err_msg)) then
              return
            end if
          end select
        TYPE IS (real(kind=r8_kind))
          select type (rmask_threshold)
          type is (real(kind=r8_kind))
            WHERE (rmask < rmask_threshold) outmask = .FALSE.
          class default
            if (fms_error_handler('fms_diag_reduction_methods_mod::init_mask_3d', 'type mismatch', err_msg)) then
              return
            end if
          end select
        CLASS DEFAULT
          if (fms_error_handler('fms_diag_reduction_methods_mod::init_mask_3d',&
            & 'The rmask is not one of the supported types of real(kind=4) or real(kind=8)', err_msg)) then
            return
          end if
      END SELECT
    END IF
  end subroutine init_mask_3d

  !> @brief Updates the buffer with the field data based on the value of the flag passed:
  !! 0 for minimum; 1 for maximum.
  subroutine fms_diag_update_extremum(flag, buffer_obj, field_data, recon_bounds, l_start, &
    l_end, is_regional, reduced_k_range, sample, mask, fieldName, hasDiurnalAxis, err_msg)
    integer, intent(in) :: flag !< Flag to indicate what to update: 0 for minimum; 1 for maximum
    class(fmsDiagOutputBuffer_class), intent(inout) :: buffer_obj !< Remapped buffer to update
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    type(fmsDiagBoundsHalos_type), intent(inout) :: recon_bounds !< Indices of bounds in the first three dimension
                                                                 !! of the field data
    integer, intent(in) :: l_start(:) !< Local starting indices for the first three dimensions
    integer, intent(in) :: l_end(:)   !< Local ending indices for the first three dimensions
    logical, intent(in) :: is_regional
    logical, intent(in) :: reduced_k_range
    integer :: sample !< Index along the diurnal time axis
    logical, intent(in) :: mask(:,:,:,:) !< Must be out of range mask
    character(len=*), intent(in) :: fieldName !< Field name for error reporting
    logical :: hasDiurnalAxis !< Flag to indicate if the buffer has a diurnal axis
    character(len=*), intent(inout) :: err_msg

    integer :: is, js, ks !< Starting indices in the I, J, and K dimensions
    integer :: ie, je, ke !< Ending indices in the I, J, and K dimensions
    integer :: hi, hj !< Halo sizes in the I, and J dimensions
    integer :: f1, f2 !< Updated starting and ending indices in the I dimension
    integer :: f3, f4 !< Updated starting and ending indices in the J dimension
    integer :: ksr, ker !< Reduced indices in the K dimension
    integer :: i, j, k, i1, j1, k1 !< For loops
    character(len=128) :: err_msg_local !< Stores local error message
    class(*), pointer :: ptr_buffer(:,:,:,:,:) !< Pointer to 5D buffer for remapping

    !> Unpack recon_bounds
    is = recon_bounds%bounds3D%get_imin()
    js = recon_bounds%bounds3D%get_jmin()
    ks = recon_bounds%bounds3D%get_kmin()
    ie = recon_bounds%bounds3D%get_imax()
    je = recon_bounds%bounds3D%get_jmax()
    ke = recon_bounds%bounds3D%get_kmax()
    hi = recon_bounds%get_hi()
    f1 = recon_bounds%get_fis()
    f2 = recon_bounds%get_fie()
    hj = recon_bounds%get_hj()
    f3 = recon_bounds%get_fjs()
    f4 = recon_bounds%get_fje()

    if (flag .ne. 0 .and. flag .ne. 1) then
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::fms_diag_update_extremum: flag must be either 0 or 1.")
    end if

    !! TODO: remap buffer before passing to subroutines update_scalar_extremum and update_array_extremum
    ptr_buffer => buffer_obj%remap_buffer(fieldName, hasDiurnalAxis)

    ! Update buffer
    IF (is_regional) THEN
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
                call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum &
                  unsupported buffer type')
              end select
            end if
          END DO
        END DO
      END DO
    ELSE
      IF (reduced_k_range) THEN
        call recon_bounds%bounds3D%set_kbounds(l_start(3), l_end(3))
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
          call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported buffer type')
        end select
      ELSE
        IF ( debug_diag_manager ) THEN
          ! Compare bounds {is-hi, ie-hi, js-hj, je-hj, ks, ke} with the bounds of first three dimensions of the buffer
          if (compare_two_sets_of_bounds((/is-hi, ie-hi, js-hj, je-hj, ks, ke/), &
            (/LBOUND(ptr_buffer,1), UBOUND(ptr_buffer,1), LBOUND(ptr_buffer,2), UBOUND(ptr_buffer,2), &
            LBOUND(ptr_buffer,3), UBOUND(ptr_buffer,3)/), err_msg_local)) THEN
            IF ( fms_error_handler('fms_diag_object_mod::fms_diag_update_extremum', err_msg_local, err_msg) ) THEN
              !if (associated(field_data)) deallocate(field_data)
              !if (allocated(mask)) deallocate(mask)
              RETURN
            END IF
          END IF
        END IF
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
          call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported buffer type')
        end select
      END IF
    end if
    ! Reset counter count_0d of the buffer object
    select type (buffer_obj)
    type is (outputBuffer0d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    type is (outputBuffer1d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    type is (outputBuffer2d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    type is (outputBuffer3d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    type is (outputBuffer4d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    type is (outputBuffer5d_type)
      select type (real_counter => buffer_obj%count_0d)
      type is (real(kind=r4_kind))
        real_counter(sample) = 1.0_r4_kind
      type is (real(kind=r8_kind))
        real_counter(sample) = 1.0_r8_kind
      class default
        call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported intrinsic type')
      end select
    class default
      call mpp_error(FATAL, 'fms_diag_reduction_methods_mod::fms_diag_update_extremum unsupported buffer type')
    end select
  end subroutine fms_diag_update_extremum

  !> @brief Updates individual element of buffer
  subroutine update_scalar_extremum(flag, field_data, buffer, mask, sample, recon_bounds, &
    running_indx1, running_indx2)
    integer, intent(in) :: flag !< 0 for minimum; 1 for maximum
    class(*), intent(in) :: field_data(:,:,:,:) !< Field data
    class(*), intent(inout) :: buffer(:,:,:,:,:) !< Remapped output buffer
    logical, intent(in) :: mask(:,:,:,:) !< Update mask
    integer, intent(in) :: sample !< diurnal sample index
    type(fmsDiagBoundsHalos_type), intent(in) :: recon_bounds !< Holds starting and ending indices in the
                                                              !! I, J, and K dimensions and
                                                              !! halo sizes in the I, and J dimensions
    integer, intent(in) :: running_indx1(3) !< Holds indices i, j, and k
    integer, intent(in) :: running_indx2(3) !< Holds indices i1, j1, and k1

    integer :: i, j, k
    integer :: i1, j1, k1
    integer :: is, js, ks
    integer :: ie, je, ke
    integer :: hi, hj

    ! Initialize i, j, and k
    i = running_indx1(1)
    j = running_indx1(2)
    k = running_indx1(3)

    ! Initialize i1, j1, and k1
    i1 = running_indx2(1)
    j1 = running_indx2(2)
    k1 = running_indx2(3)

    !> Unpack bounds (/is, js, ks, ie, je, ke, hi, f1, f2, hj, f3, f4/)
      is = recon_bounds%bounds3D%get_imin()
      js = recon_bounds%bounds3D%get_jmin()
      ks = recon_bounds%bounds3D%get_kmin()
      ie = recon_bounds%bounds3D%get_imax()
      je = recon_bounds%bounds3D%get_jmax()
      ke = recon_bounds%bounds3D%get_kmax()
      hi = recon_bounds%get_hi()
      hj = recon_bounds%get_hj()

    ! Select proper type and update the buffer
    select type (field_data)
    type is (real(kind=r4_kind))
      select type (buffer)
      type is (real(kind=r4_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum type mismatch")
      end select
    type is (real(kind=r8_kind))
      select type (buffer)
      type is (real(kind=r8_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum type mismatch")
      end select
    type is (integer(kind=i4_kind))
      select type (buffer)
      type is (integer(kind=i4_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum type mismatch")
      end select
    type is (integer(kind=i8_kind))
      select type (buffer)
      type is (integer(kind=i8_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum type mismatch")
      end select
    class default
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_scalar_extremum unsupported field data type")
    end select
  end subroutine update_scalar_extremum

  !> @brief Updates a chunk of buffer
  subroutine update_array_extremum(flag, field_data, buffer, mask, sample, recon_bounds, reduced_k_range)
    integer :: flag !< 0 for minimum; 1 for extremum
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

    !> Unpack bounds (/is, js, ks, ie, je, ke, hi, f1, f2, hj, f3, f4/)
    is = recon_bounds%bounds3D%get_imin()
    js = recon_bounds%bounds3D%get_jmin()
    ks = recon_bounds%bounds3D%get_kmin()
    ie = recon_bounds%bounds3D%get_imax()
    je = recon_bounds%bounds3D%get_jmax()
    ke = recon_bounds%bounds3D%get_kmax()
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
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum type mismatch")
      end select
    type is (real(kind=r8_kind))
      select type (buffer)
      type is (real(kind=r8_kind))
          if (flag .eq. 0) then
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
          end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum type mismatch")
      end select
    type is (integer(kind=i4_kind))
      select type (buffer)
      type is (integer(kind=i4_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum type mismatch")
      end select
    type is (integer(kind=i8_kind))
      select type (buffer)
      type is (integer(kind=i8_kind))
        if (flag .eq. 0) then
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
        end if
      class default
        call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum type mismatch")
      end select
    class default
      call mpp_error( FATAL, "fms_diag_reduction_methods_mod::update_array_extremum unsupported field data type")
    end select
  end subroutine update_array_extremum
#endif
end module fms_diag_reduction_methods_mod