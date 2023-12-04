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
!> @defgroup coupler_types_mod coupler_types_mod
!> @ingroup coupler
!> @brief This module contains type declarations for the coupler.
!> @author Richard Slater, John Dunne

!> @addtogroup coupler_types_mod
!> @{

!> @addtogroup coupler_types_mod
!> @{

  !> @brief Rescales the fields in the fields in the elements of a coupler_2d_bc_type
  !!
  !! Rescales the fields in the elements of a coupler_2d_bc_type by multiplying by a factor scale.
  !! If scale is 0, this is a direct assignment to 0, so that NaNs will not persist.
  subroutine CT_RESCALE_DATA_2D_(var, scale, halo_size, bc_index, field_index,&
      & exclude_flux_type, only_flux_type, pass_through_ice)
    type(coupler_2d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being rescaled
    real(FMS_CP_KIND_),                       intent(in)    :: scale   !< A scaling factor to multiply fields by
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default or
                                                           !! the full arrays if scale is 0.
    integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types
                                                           !! of fluxes to exclude from this copy.
    character(len=*), optional, intent(in)    :: only_flux_type !< A string describing which types
                                                           !! of fluxes to include from this copy.
    logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                           !! value of pass_through ice matches this

    logical :: do_bc
    integer :: i, j, m, n, n1, n2, halo
    integer, parameter :: kindl = FMS_CP_KIND_

    if (present(bc_index)) then
      if (bc_index > var%num_bcs)&
          & call mpp_error(FATAL, "CT_rescale_data_2d: bc_index is present and exceeds var%num_bcs.")
      if (present(field_index)) then ; if (field_index > var%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
          & call mpp_error(FATAL, "CT_rescale_data_2d: field_index is present and exceeds num_fields for" //&
          & trim(var%FMS_CP_BC_TYPE_(bc_index)%field_name) )
      endif
    elseif (present(field_index)) then
      call mpp_error(FATAL, "CT_rescale_data_2d: bc_index must be present if field_index is present.")
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size

    n1 = 1
    n2 = var%num_bcs
    if (present(bc_index)) then
      n1 = bc_index
      n2 = bc_index
    endif

    if (n2 >= n1) then
      ! A more consciencious implementation would include a more descriptive error messages.
      if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
          & call mpp_error(FATAL, "CT_rescale_data_2d: Excessive i-direction halo size.")
      if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
          & call mpp_error(FATAL, "CT_rescale_data_2d: Excessive j-direction halo size.")
    endif

    do n = n1, n2
      do_bc = .true.
      if (do_bc .and. present(exclude_flux_type))&
          & do_bc = .not.(trim(var%FMS_CP_BC_TYPE_(n)%flux_type) == trim(exclude_flux_type))
      if (do_bc .and. present(only_flux_type))&
          & do_bc = (trim(var%FMS_CP_BC_TYPE_(n)%flux_type) == trim(only_flux_type))
      if (do_bc .and. present(pass_through_ice))&
          & do_bc = (pass_through_ice .eqv. var%FMS_CP_BC_TYPE_(n)%pass_through_ice)
      if (.not.do_bc) cycle

      do m = 1, var%FMS_CP_BC_TYPE_(n)%num_fields
        if (present(field_index)) then
          if (m /= field_index) cycle
        endif
        if ( associated(var%FMS_CP_BC_TYPE_(n)%field(m)%values) ) then
          if (scale == 0.0_kindl) then
            if (present(halo_size)) then
              do j=var%jsc-halo,var%jec+halo
                do i=var%isc-halo,var%iec+halo
                  var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j) = 0.0_kindl
                enddo
              enddo
            else
              var%FMS_CP_BC_TYPE_(n)%field(m)%values(:,:) = 0.0_kindl
            endif
          else
            do j=var%jsc-halo,var%jec+halo
              do i=var%isc-halo,var%iec+halo
                var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j) = scale * var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j)
              enddo
            enddo
          endif
        endif
      enddo
    enddo
  end subroutine CT_RESCALE_DATA_2D_

  !> @brief Rescales the fields in the elements of a coupler_3d_bc_type
  !!
  !! This subroutine rescales the fields in the elements of a coupler_3d_bc_type by multiplying by a
  !! factor scale.  If scale is 0, this is a direct assignment to 0, so that NaNs will not persist.
  subroutine CT_RESCALE_DATA_3D_(var, scale, halo_size, bc_index, field_index,&
      & exclude_flux_type, only_flux_type, pass_through_ice)
    type(coupler_3d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being rescaled
    real(FMS_CP_KIND_),                       intent(in)    :: scale   !< A scaling factor to multiply fields by
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default or
                                                           !! the full arrays if scale is 0.
    integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
    integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
    character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types
                                                         !! of fluxes to exclude from this copy.
    character(len=*), optional, intent(in)    :: only_flux_type !< A string describing which types of
                                                         !! fluxes to include from this copy.
    logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                         !! value of pass_through ice matches this

    logical :: do_bc
    integer :: i, j, k, m, n, n1, n2, halo
    integer, parameter :: kindl = FMS_CP_KIND_

    if (present(bc_index)) then
      if (bc_index > var%num_bcs)&
          & call mpp_error(FATAL, "CT_rescale_data_2d: bc_index is present and exceeds var%num_bcs.")
      if (present(field_index)) then ; if (field_index > var%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
          & call mpp_error(FATAL, "CT_rescale_data_2d: field_index is present and exceeds num_fields for" //&
          & trim(var%FMS_CP_BC_TYPE_(bc_index)%field_name) )
      endif
    elseif (present(field_index)) then
      call mpp_error(FATAL, "CT_rescale_data_2d: bc_index must be present if field_index is present.")
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size

    n1 = 1
    n2 = var%num_bcs
    if (present(bc_index)) then
      n1 = bc_index
      n2 = bc_index
    endif

    if (n2 >= n1) then
      ! A more consciencious implementation would include a more descriptive error messages.
      if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
          & call mpp_error(FATAL, "CT_rescale_data_3d: Excessive i-direction halo size.")
      if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
          & call mpp_error(FATAL, "CT_rescale_data_3d: Excessive j-direction halo size.")
    endif

    do n = n1, n2
      do_bc = .true.
      if (do_bc .and. present(exclude_flux_type))&
          & do_bc = .not.(trim(var%FMS_CP_BC_TYPE_(n)%flux_type) == trim(exclude_flux_type))
      if (do_bc .and. present(only_flux_type))&
          & do_bc = (trim(var%FMS_CP_BC_TYPE_(n)%flux_type) == trim(only_flux_type))
      if (do_bc .and. present(pass_through_ice))&
          & do_bc = (pass_through_ice .eqv. var%FMS_CP_BC_TYPE_(n)%pass_through_ice)
      if (.not.do_bc) cycle

      do m = 1, var%FMS_CP_BC_TYPE_(n)%num_fields
        if (present(field_index)) then
          if (m /= field_index) cycle
        endif
        if ( associated(var%FMS_CP_BC_TYPE_(n)%field(m)%values) ) then
          if (scale == 0.0_kindl) then
            if (present(halo_size)) then
              do k=var%ks,var%ke
                do j=var%jsc-halo,var%jec+halo
                  do i=var%isc-halo,var%iec+halo
                    var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j,k) = 0.0_kindl
                  enddo
                enddo
              enddo
            else
              var%FMS_CP_BC_TYPE_(n)%field(m)%values(:,:,:) = 0.0_kindl
            endif
          else
            do k=var%ks,var%ke
              do j=var%jsc-halo,var%jec+halo
                do i=var%isc-halo,var%iec+halo
                  var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j,k) = scale * var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j,k)
                enddo
              enddo
            enddo
          endif
        endif
      enddo
    enddo
  end subroutine CT_RESCALE_DATA_3D_

  !> @brief Increment data in the elements of a coupler_2d_bc_type with weighted averages of elements of a
  !! coupler_3d_bc_type
  !!
  !! Increments the data in the elements of a coupler_2d_bc_type with the weighed average of the
  !! elements of a coupler_3d_bc_type. Both must have the same horizontal array sizes and the
  !! normalized weight array must match the array sizes of the coupler_3d_bc_type.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index is present and exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "bc_index must be present if field_index is present."
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "There is an k-direction computational domain size mismatch."
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "weights array must be the i-size of a computational or data domain."
  !! @throw FATAL, "weights array must be the j-size of a computational or data domain."
  subroutine CT_INCREMENT_DATA_2D_3D_(var_in, weights, var, halo_size, bc_index, field_index,&
      & scale_factor, scale_prev, exclude_flux_type, only_flux_type, pass_through_ice)
    type(coupler_3d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to add to the other type
    real(FMS_CP_KIND_), dimension(:,:,:),     intent(in)    :: weights !< An array of normalized weights for the 3d-data
                                                        !! to increment the 2d-data.  There is no renormalization,
                                                        !! so if the weights do not sum to 1 in the 3rd dimension
                                                        !! there may be adverse consequences!
    type(coupler_2d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being incremented
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
    integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
    real(FMS_CP_KIND_),             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that
                                                                             !! is being added
    real(FMS_CP_KIND_),             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that
                                                                      !! is already here and is being added
    character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types
                                                         !! of fluxes to exclude from this increment.
    character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types
                                                         !! of fluxes to include from this increment.
    logical,          optional, intent(in)    :: pass_through_ice !< If true, only increment BCs whose
                                                         !! value of pass_through ice matches this

    real(FMS_CP_KIND_) :: scale, sc_prev
    logical :: increment_bc
    integer :: i, j, k, m, n, n1, n2, halo
    integer :: io1, jo1, iow, jow, kow  ! Offsets to account for different index conventions.
    integer, parameter :: kindl = FMS_CP_KIND_

    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor
    sc_prev = 1.0_kindl
    if (present(scale_prev)) sc_prev = scale_prev

    if (present(bc_index)) then
      if (bc_index > var_in%num_bcs)&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: bc_index is present and exceeds var_in%num_bcs.")
      if (present(field_index)) then ; if (field_index > var_in%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: field_index is present and exceeds num_fields for" //&
          & trim(var_in%FMS_CP_BC_TYPE_(bc_index)%field_name) )
      endif
    elseif (present(field_index)) then
      call mpp_error(FATAL, "CT_increment_data_2d_3d: bc_index must be present if field_index is present.")
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size

    n1 = 1
    n2 = var_in%num_bcs
    if (present(bc_index)) then
      n1 = bc_index
      n2 = bc_index
    endif

    if (n2 >= n1) then
      ! A more consciencious implementation would include a more descriptive error messages.
      if ((var_in%iec-var_in%isc) /= (var%iec-var%isc))&
          & call mpp_error(FATAL, &
                           &  "CT_increment_data_2d_3d: There is an i-direction computational domain size mismatch.")
      if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc))&
          & call mpp_error(FATAL, &
                           &  "CT_increment_data_2d_3d: There is a j-direction computational domain size mismatch.")
      if ((1+var_in%ke-var_in%ks) /= size(weights,3))&
          & call mpp_error(FATAL, &
                           &  "CT_increment_data_2d_3d: There is a k-direction size mismatch with the weights array.")
      if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo))&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive i-direction halo size for the input structure.")
      if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo))&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive j-direction halo size for the input structure.")
      if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive i-direction halo size for the output structure.")
      if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
          & call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive j-direction halo size for the output structure.")

      if ((1+var%iec-var%isc) == size(weights,1)) then
        iow = 1 - var%isc
      elseif ((1+var%ied-var%isd) == size(weights,1)) then
        iow = 1 - var%isd
      elseif ((1+var_in%ied-var_in%isd) == size(weights,1)) then
        iow = 1 + (var_in%isc - var_in%isd) - var%isc
      else
        call mpp_error(FATAL, &
                   &  "CT_increment_data_2d_3d: weights array must be the i-size of a computational or data domain.")
      endif
      if ((1+var%jec-var%jsc) == size(weights,2)) then
        jow = 1 - var%jsc
      elseif ((1+var%jed-var%jsd) == size(weights,2)) then
        jow = 1 - var%jsd
      elseif ((1+var_in%jed-var_in%jsd) == size(weights,2)) then
        jow = 1 + (var_in%jsc - var_in%jsd) - var%jsc
      else
        call mpp_error(FATAL, &
                   &  "CT_increment_data_2d_3d: weights array must be the j-size of a computational or data domain.")
      endif

      io1 = var_in%isc - var%isc
      jo1 = var_in%jsc - var%jsc
      kow = 1 - var_in%ks
    endif

    do n = n1, n2
      increment_bc = .true.
      if (increment_bc .and. present(exclude_flux_type))&
          & increment_bc = .not.(trim(var_in%FMS_CP_BC_TYPE_(n)%flux_type) == trim(exclude_flux_type))
      if (increment_bc .and. present(only_flux_type))&
          & increment_bc = (trim(var_in%FMS_CP_BC_TYPE_(n)%flux_type) == trim(only_flux_type))
      if (increment_bc .and. present(pass_through_ice))&
          & increment_bc = (pass_through_ice .eqv. var_in%FMS_CP_BC_TYPE_(n)%pass_through_ice)
      if (.not.increment_bc) cycle

      do m = 1, var_in%FMS_CP_BC_TYPE_(n)%num_fields
        if (present(field_index)) then
          if (m /= field_index) cycle
        endif
        if ( associated(var%FMS_CP_BC_TYPE_(n)%field(m)%values) ) then
          do k=var_in%ks,var_in%ke
            do j=var%jsc-halo,var%jec+halo
              do i=var%isc-halo,var%iec+halo
                var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j) = sc_prev * var%FMS_CP_BC_TYPE_(n)%field(m)%values(i,j) +&
                    & (scale * weights(i+iow,j+jow,k+kow)) * var_in%FMS_CP_BC_TYPE_(n)%field(m)%values(i+io1,j+io1,k)
              enddo
            enddo
          enddo
        endif
      enddo
    enddo
  end subroutine CT_INCREMENT_DATA_2D_3D_

  !> @brief Extract a 2d field from a coupler_2d_bc_type
  !!
  !! Extract a single 2-d field from a coupler_2d_bc_type into a two-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  subroutine CT_EXTRACT_DATA_2D_(var_in, bc_index, field_index, array_out,&
      & scale_factor, halo_size, idim, jdim)
    type(coupler_2d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    real(FMS_CP_KIND_), dimension(1:,1:),     intent(out)   :: array_out !< The recipient array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    real(FMS_CP_KIND_), optional, intent(in) :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list

    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_extract_data_2d):'
    character(len=400)      :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, halo, i_off, j_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) then
      array_out(:,:) = 0.0_kindl
      return
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var_in%num_bcs)&
        & call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
    if (field_index > var_in%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var_in%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var_in%isc)
    else
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_out,1), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var_in%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
    else
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_out,2), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var_in%jsc-halo)
    endif

    do j=var_in%jsc-halo,var_in%jec+halo
      do i=var_in%isc-halo,var_in%iec+halo
        array_out(i+i_off,j+j_off) = scale * var_in%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j)
      enddo
    enddo
  end subroutine CT_EXTRACT_DATA_2D_

  !> @brief Extract a single k-level of a 3d field from a coupler_3d_bc_type
  !!
  !! Extract a single k-level of a 3-d field from a coupler_3d_bc_type into a two-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The extracted k-index of 'k' is outside of the valid range of 'ks' to 'ke'"
  subroutine CT_EXTRACT_DATA_3D_2D_(var_in, bc_index, field_index, k_in, array_out,&
      & scale_factor, halo_size, idim, jdim)
    type(coupler_3d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    integer,                    intent(in)    :: k_in      !< The k-index to extract
    real(FMS_CP_KIND_), dimension(1:,1:),     intent(out)   :: array_out !< The recipient array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    real(FMS_CP_KIND_), optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list
    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_extract_data_3d_2d):'
    character(len=400)      :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, halo, i_off, j_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) then
      array_out(:,:) = 0.0_kindl
      return
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var_in%num_bcs)&
        & call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
    if (field_index > var_in%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var_in%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var_in%isc)
    else
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_out,1), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var_in%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
    else
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_out,2), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var_in%jsc-halo)
    endif

    if ((k_in > var_in%ke) .or. (k_in < var_in%ks)) then
      write (error_msg, *) trim(error_header), ' The extracted k-index of ', k_in,&
          & ' is outside of the valid range of ', var_in%ks, ' to ', var_in%ke
      call mpp_error(FATAL, trim(error_msg))
    endif

    do j=var_in%jsc-halo,var_in%jec+halo
      do i=var_in%isc-halo,var_in%iec+halo
        array_out(i+i_off,j+j_off) = scale * var_in%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j,k_in)
      enddo
    enddo
  end subroutine CT_EXTRACT_DATA_3D_2D_

  !> @brief Extract single 3d field from a coupler_3d_bc_type
  !!
  !! Extract a single 3-d field from a coupler_3d_bc_type into a three-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with k-dimension size 'n' does not match the data of size 'd'"
  subroutine CT_EXTRACT_DATA_3D_(var_in, bc_index, field_index, array_out,&
      & scale_factor, halo_size, idim, jdim)
    type(coupler_3d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    real(FMS_CP_KIND_), dimension(1:,1:,1:),  intent(out)   :: array_out !< The recipient array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    real(FMS_CP_KIND_), optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list

    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_extract_data_3d):'
    character(len=400) :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, k, halo, i_off, j_off, k_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) then
      array_out(:,:,:) = 0.0_kindl
      return
    endif

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var_in%num_bcs)&
        & call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
    if (field_index > var_in%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var_in%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var_in%isc)
    else
      if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_out,1), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%iec - var_in%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var_in%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
    else
      if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_out,2), ' does not match the data of size ',&
            & (2*halo + 1 + var_in%jec - var_in%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var_in%jsc-halo)
    endif

    if (size(array_out,3) /= 1 + var_in%ke - var_in%ks) then
      write (error_msg, *) trim(error_header), ' The target array with k-dimension size ',&
          & size(array_out,3), ' does not match the data of size ',&
          & (1 + var_in%ke - var_in%ks)
      call mpp_error(FATAL, trim(error_msg))
    endif
    k_off = 1 - var_in%ks

    do k=var_in%ks,var_in%ke
      do j=var_in%jsc-halo,var_in%jec+halo
        do i=var_in%isc-halo,var_in%iec+halo
          array_out(i+i_off,j+j_off,k+k_off) = scale * var_in%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j,k)
        enddo
      enddo
    enddo
  end subroutine CT_EXTRACT_DATA_3D_

  !> @brief Set single 2d field in coupler_3d_bc_type
  !!
  !! Set a single 2-d field in a coupler_3d_bc_type from a two-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  subroutine CT_SET_DATA_2D_(array_in, bc_index, field_index, var,&
      & scale_factor, halo_size, idim, jdim)
    real(FMS_CP_KIND_), dimension(1:,1:),     intent(in)   :: array_in   !< The source array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    type(coupler_2d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to set
    real(FMS_CP_KIND_), optional, intent(in)  :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list
    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_set_data_2d):'
    character(len=400) :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, halo, i_off, j_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) return

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var%num_bcs) &
        call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
    if (field_index > var%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%iec-var%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var%isc)
    else
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_in,1), ' does not match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%jec-var%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
    else
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_in,2), ' does not match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var%jsc-halo)
    endif

    do j=var%jsc-halo,var%jec+halo
      do i=var%isc-halo,var%iec+halo
        var%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j) = scale * array_in(i+i_off,j+j_off)
      enddo
    enddo
  end subroutine CT_SET_DATA_2D_

  !> @brief Set one k-level of a single 3d field in a coupler_3d_bc_type
  !!
  !! This subroutine sets a one k-level of a single 3-d field in a coupler_3d_bc_type from a
  !! two-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The k-index of 'k' is outside of the valid range of 'ks' to 'ke'"
  subroutine CT_SET_DATA_2D_3D_(array_in, bc_index, field_index, k_out, var,&
      & scale_factor, halo_size, idim, jdim)
    real(FMS_CP_KIND_), dimension(1:,1:),     intent(in)    :: array_in  !< The source array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    integer,                    intent(in)    :: k_out     !< The k-index to set
    type(coupler_3d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to be set
    real(FMS_CP_KIND_), optional, intent(in)  :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list

    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_set_data_3d_2d):'
    character(len=400)      :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, halo, i_off, j_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) return

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var%num_bcs)&
        & call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
    if (field_index > var%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%iec-var%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var%isc)
    else
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_in,1), ' does not match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%jec-var%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
    else
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_in,2), ' does not match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var%jsc-halo)
    endif

    if ((k_out > var%ke) .or. (k_out < var%ks)) then
      write (error_msg, *) trim(error_header), ' The k-index of ', k_out,&
          & ' is outside of the valid range of ', var%ks, ' to ', var%ke
      call mpp_error(FATAL, trim(error_msg))
    endif

    do j=var%jsc-halo,var%jec+halo
      do i=var%isc-halo,var%iec+halo
        var%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j,k_out) = scale * array_in(i+i_off,j+j_off)
      enddo
    enddo
  end subroutine CT_SET_DATA_2D_3D_

  !> @brief Set a single 3d field in a coupler_3d_bc_type
  !!
  !! This subroutine sets a single 3-d field in a coupler_3d_bc_type from a three-dimensional array.
  !!
  !! @throw FATAL, "bc_index is present and exceeds var_in%num_bcs."
  !! @throw FATAL, "field_index exceeds num_fields for var_in%FMS_CP_BC_TYPE_(bc_incdx)%field_name"
  !! @throw FATAL, "Excessive i-direction halo size for the input structure."
  !! @throw FATAL, "Excessive j-direction halo size for the input structure."
  !! @throw FATAL, "Disordered i-dimension index bound list"
  !! @throw FATAL, "Disordered j-dimension index bound list"
  !! @throw FATAL, "The declared i-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "The declared j-dimension size of 'n' does not match the actual size of 'a'"
  !! @throw FATAL, "There is an i-direction computational domain size mismatch."
  !! @throw FATAL, "There is an j-direction computational domain size mismatch."
  !! @throw FATAL, "The target array with i-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with j-dimension size 'n' is too small to match the data of size 'd'"
  !! @throw FATAL, "The target array with K-dimension size 'n' is too small to match the data of size 'd'"
  subroutine CT_SET_DATA_3D_(array_in, bc_index, field_index, var,&
      & scale_factor, halo_size, idim, jdim)
    real(FMS_CP_KIND_), dimension(1:,1:,1:),  intent(in)    :: array_in  !< The source array for the field; its size
                                                           !! must match the size of the data being copied
                                                           !! unless idim and jdim are supplied.
    integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                           !! that is being copied
    integer,                    intent(in)    :: field_index !< The index of the field in the
                                                           !! boundary condition that is being copied
    type(coupler_3d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to be set
    real(FMS_CP_KIND_), optional, intent(in)  :: scale_factor !< A scaling factor for the data that is being added
    integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                           !! the first dimension of the output array
                                                           !! in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                           !! the second dimension of the output array
                                                           !! in a non-decreasing list

    character(len=*), parameter :: error_header =&
        & '==>Error from coupler_types_mod (CT_set_data_3d):'
    character(len=400) :: error_msg

    real(FMS_CP_KIND_) :: scale
    integer :: i, j, k, halo, i_off, j_off, k_off
    integer, parameter :: kindl = FMS_CP_KIND_

    if (bc_index <= 0) return

    halo = 0
    if (present(halo_size)) halo = halo_size
    scale = 1.0_kindl
    if (present(scale_factor)) scale = scale_factor

    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo))&
        & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

    if (bc_index > var%num_bcs)&
        & call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
    if (field_index > var%FMS_CP_BC_TYPE_(bc_index)%num_fields)&
        & call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
        & trim(var%FMS_CP_BC_TYPE_(bc_index)%field_name) )

    ! Do error checking on the i-dimension and determine the array offsets.
    if (present(idim)) then
      if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,1) /= (1+idim(4)-idim(1))) then
        write (error_msg, *) trim(error_header), ' The declared i-dimension size of ',&
            & (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%iec-var%isc) /= (idim(3)-idim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an i-direction computational domain size mismatch.")
      if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & (1+idim(4)-idim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      i_off = (1-idim(1)) + (idim(2)-var%isc)
    else
      if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
        write (error_msg, *) trim(error_header), ' The target array with i-dimension size ',&
            & size(array_in,1), ' does not match the data of size ',&
            & (2*halo + 1 + var%iec - var%isc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      i_off = 1 - (var%isc-halo)
    endif

    ! Do error checking on the j-dimension and determine the array offsets.
    if (present(jdim)) then
      if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
        write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
        call mpp_error(FATAL, trim(error_msg))
      endif
      if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
        write (error_msg, *) trim(error_header), ' The declared j-dimension size of ',&
            & (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
        call mpp_error(FATAL, trim(error_msg))
      endif
      if ((var%jec-var%jsc) /= (jdim(3)-jdim(2)))&
          & call mpp_error(FATAL, trim(error_header)//" There is an j-direction computational domain size mismatch.")
      if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo))&
          & call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & (1+jdim(4)-jdim(1)), ' is too small to match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif

      j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
    else
      if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
        write (error_msg, *) trim(error_header), ' The target array with j-dimension size ',&
            & size(array_in,2), ' does not match the data of size ',&
            & (2*halo + 1 + var%jec - var%jsc)
        call mpp_error(FATAL, trim(error_msg))
      endif
      j_off = 1 - (var%jsc-halo)
    endif

    if (size(array_in,3) /= 1 + var%ke - var%ks) then
      write (error_msg, *) trim(error_header), ' The target array with k-dimension size ',&
          & size(array_in,3), ' does not match the data of size ',&
          & (1 + var%ke - var%ks)
      call mpp_error(FATAL, trim(error_msg))
    endif
    k_off = 1 - var%ks

    do k=var%ks,var%ke
      do j=var%jsc-halo,var%jec+halo
        do i=var%isc-halo,var%iec+halo
          var%FMS_CP_BC_TYPE_(bc_index)%field(field_index)%values(i,j,k) = scale * array_in(i+i_off,j+j_off,k+k_off)
        enddo
      enddo
    enddo
  end subroutine CT_SET_DATA_3D_
