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

!> @defgroup fms_diag_bbox_mod fms_diag_bbox_mod
!> @ingroup diag_manager
!> @brief fms_diag_bbox_mod defines classes encapsulating bounding boxes
!!   and interval bounds.
!!
!> @author Miguel Zuniga
!!
!> @file
!> @brief File for @ref fms_diag_bbox_mod
!> @addtogroup fms_diag_bbox_mod
!> @{
MODULE fms_diag_bbox_mod

   USE fms_mod, ONLY: error_mesg, FATAL, fms_error_handler, string

   implicit none

!> @brief Data structure holding a 3D bounding box. It is commonlyused to
!! represent the interval bounds or limits of a 3D sub-array such as the
!! array index bounds of the spatial component a diag_manager field output
!! buffer array.
   TYPE, public :: fmsDiagIbounds_type
      INTEGER :: imin !< Lower i bound.
      INTEGER :: imax !< Upper i bound.
      INTEGER :: jmin !< Lower j bound.
      INTEGER :: jmax !< Upper j bound.
      INTEGER :: kmin !< Lower k bound.
      INTEGER :: kmax !< Upper k bound.
      logical :: has_halos !< .True. if the buffer has halos
      integer :: nhalo_I !< Number of halos in i
      integer :: nhalo_J !< Number of halos in j
   contains
      procedure :: reset => reset_bounds
      procedure :: reset_bounds_from_array_4D
      procedure :: reset_bounds_from_array_5D
      procedure :: update_bounds
      procedure :: set_bounds
      procedure :: rebase_input
      procedure :: rebase_output
      procedure :: get_imin
      procedure :: get_imax
      procedure :: get_jmin
      procedure :: get_jmax
      procedure :: get_kmin
      procedure :: get_kmax
      procedure :: update_index
   END TYPE fmsDiagIbounds_type

   !> @brief Data structure holding starting and ending indices in the I, J, and
   !! K dimensions. It also has extra members related to halo sizes and updated indices
   !! in I and J dimensions.
   type, public :: fmsDiagBoundsHalos_type
      private
      type(fmsDiagIbounds_type) :: bounds3D !< Holds starting and ending indices of
                                            !! the I, J, and K dimensions
      integer :: hi !< Halo size in the I dimension
      integer :: hj !< Halo size in the J dimension
      integer :: fis !< Updated starting index in the I dimension
      integer :: fie !< Updated ending index in the I dimension
      integer :: fjs !< Updated starting index in the J dimension
      integer :: fje !< Updated ending index in the J dimension
      contains
      procedure :: get_hi
      procedure :: get_hj
      procedure :: get_fis
      procedure :: get_fie
      procedure :: get_fjs
      procedure :: get_fje
   end type fmsDiagBoundsHalos_type

   public :: recondition_indices, determine_if_block_is_in_region

   integer, parameter :: xdimension = 1 !< Parameter defining the x dimension
   integer, parameter :: ydimension = 2 !< Parameter defining the y dimension
   integer, parameter :: zdimension = 3 !< Parameter defininf the z dimension

CONTAINS

!> @brief The PEs grid points are divided further into "blocks". This function determines if a block
! has data for a given subregion and dimension
!! @return .true. if the a subergion is inside a block
logical pure function determine_if_block_is_in_region(subregion_start, subregion_end, bounds, dim)
  integer,                   intent(in) :: subregion_start !< Begining of the subregion
  integer,                   intent(in) :: subregion_end   !< Ending of the subregion
  type(fmsDiagIbounds_type), intent(in) :: bounds          !< Starting and ending of the subregion
  integer,                   intent(in) :: dim             !< Dimension to check

  integer :: block_start !< Begining index of the block
  integer :: block_end   !< Ending index of the block

  determine_if_block_is_in_region = .true.
  select case (dim)
  case (xdimension)
    block_start = bounds%imin
    block_end = bounds%imax
  case (ydimension)
    block_start = bounds%jmin
    block_end = bounds%jmax
  case (zdimension)
    block_start = bounds%kmin
    block_end = bounds%kmax
  end select

  if (block_start < subregion_start .and. block_end < subregion_start) then
    determine_if_block_is_in_region = .false.
    return
  endif

  if (block_start > subregion_end   .and. block_end > subregion_end) then
    determine_if_block_is_in_region = .false.
    return
  endif

  determine_if_block_is_in_region = .true.
end function determine_if_block_is_in_region

   !> @brief Gets imin of fmsDiagIbounds_type
   !! @return copy of integer member imin
   pure integer function get_imin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%imin
   end function get_imin

   !> @brief Gets imax of fmsDiagIbounds_type
   !! @return copy of integer member imax
   pure integer function get_imax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%imax
   end function get_imax

   !> @brief Gets jmin of fmsDiagIbounds_type
   !! @return copy of integer member jmin
   pure integer function get_jmin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%jmin
   end function get_jmin

   !> @brief Gets jmax of fmsDiagIbounds_type
   !! @return copy of integer member jmax
   pure integer function get_jmax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%jmax
   end function get_jmax


   !> @brief Gets kmin of fmsDiagIbounds_type
   !! @return copy of integer member kmin
   pure integer function get_kmin (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%kmin
   end function get_kmin

   !> @brief Gets kmax of fmsDiagIbounds_type
   !! @return copy of integer member kmax
   pure integer function get_kmax (this) result(rslt)
      class (fmsDiagIbounds_type), intent(in) :: this !< The !< ibounds instance
      rslt = this%kmax
   end function get_kmax

   !> @brief Updates the starting and ending index of a given dimension
   subroutine update_index(this, starting_index, ending_index, dim, ignore_halos)
     class (fmsDiagIbounds_type), intent(inout) :: this           !< The bounding box to update
     integer,                     intent(in)    :: starting_index !< Starting index to update to
     integer,                     intent(in)    :: ending_index   !< Ending index to update to
     integer,                     intent(in)    :: dim            !< Dimension to update
     logical,                     intent(in)    :: ignore_halos   !< If .true. halos will be ignored
                                                                  !! i.e output buffers can ignore halos as
                                                                  !! they do not get updates. The indices of the
                                                                  !! Input buffers need to add the number of halos
                                                                  !! so math is done only on the compute domain

     integer :: nhalox !< Number of halos in x
     integer :: nhaloy !< Number of halos in y

     if (ignore_halos) then
      nhalox = 0
      nhaloy = 0
     else
      nhalox= this%nhalo_I
      nhaloy= this%nhalo_J
     endif
     select case(dim)
     case (xdimension)
      this%imin = starting_index + nhalox
      this%imax = ending_index + nhalox
     case (ydimension)
      this%jmin = starting_index + nhaloy
      this%jmax = ending_index + nhaloy
     case (zdimension)
      this%kmin = starting_index
      this%kmax = ending_index
     end select
   end subroutine

   !> @brief Gets the halo size of fmsDiagBoundsHalos_type in the I dimension
   !! @return copy of integer member hi
   pure integer function get_hi (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%hi
   end function get_hi

   !> @brief Gets the halo size of fmsDiagBoundsHalos_type in the J dimension
   !! @return copy of integer member hj
   pure integer function get_hj (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%hj
   end function get_hj

   !> @brief Gets the updated index `fis' of fmsDiagBoundsHalos_type in the I dimension
   !! @return copy of integer member `fis'
   pure integer function get_fis (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%fis
   end function get_fis

   !> @brief Gets the updated index `fie' of fmsDiagBoundsHalos_type in the I dimension
   !! @return copy of integer member `fie'
   pure integer function get_fie (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%fie
   end function get_fie

   !> @brief Gets the updated index `fjs' of fmsDiagBoundsHalos_type in the I dimension
   !! @return copy of integer member `fjs'
   pure integer function get_fjs (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%fjs
   end function get_fjs

   !> @brief Gets the updated index `fje' of fmsDiagBoundsHalos_type in the I dimension
   !! @return copy of integer member `fje'
   pure integer function get_fje (this) result(rslt)
      class (fmsDiagBoundsHalos_type), intent(in) :: this !< Calling object
      rslt = this%fje
   end function get_fje

   !> @brief Reset the instance bounding lower and upper bounds to lower_val and upper_val, respectively.
   SUBROUTINE reset_bounds (this, lower_val, upper_val)
      class (fmsDiagIbounds_type), target, intent(inout) :: this   !< ibounds instance
      integer, intent(in) :: lower_val  !< value for the lower bounds in each dimension
      integer, intent(in) :: upper_val  !< value for the upper bounds in each dimension
      this%imin = lower_val
      this%jmin = lower_val
      this%kmin = lower_val
      this%imax = upper_val
      this%jmax = upper_val
      this%kmax = upper_val
   END SUBROUTINE reset_bounds

   !> @brief Update the the first three (normally  x, y, and z)  min and
   !! max boundaries (array indices) of the instance bounding box
   !! the six specified bounds values.
   SUBROUTINE update_bounds(this, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k)
      CLASS  (fmsDiagIbounds_type), intent(inout) :: this !<The bounding box of the output field buffer inindex space.
      INTEGER, INTENT(in) :: lower_i !< Lower i bound.
      INTEGER, INTENT(in) :: upper_i !< Upper i bound.
      INTEGER, INTENT(in) :: lower_j !< Lower j bound.
      INTEGER, INTENT(in) :: upper_j !< Upper j bound.
      INTEGER, INTENT(in) :: lower_k !< Lower k bound.
      INTEGER, INTENT(in) :: upper_k !< Upper k bound.
      this%imin = MIN(this%imin, lower_i)
      this%imax = MAX(this%imax, upper_i)
      this%jmin = MIN(this%jmin, lower_j)
      this%jmax = MAX(this%jmax, upper_j)
      this%kmin = MIN(this%kmin, lower_k)
      this%kmax = MAX(this%kmax, upper_k)
   END SUBROUTINE update_bounds

   !> @brief Sets the bounds of a bounding region
   !! @return empty string if sucessful or error message if unsucessful
   function set_bounds(this, field_data, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k, has_halos) &
      result(error_msg)
      CLASS  (fmsDiagIbounds_type), intent(inout) :: this                !< The bounding box of the field
      class(*),                     intent(in)    :: field_data(:,:,:,:) !< Field data
      INTEGER,                      INTENT(in)    :: lower_i             !< Lower i bound.
      INTEGER,                      INTENT(in)    :: upper_i             !< Upper i bound.
      INTEGER,                      INTENT(in)    :: lower_j             !< Lower j bound.
      INTEGER,                      INTENT(in)    :: upper_j             !< Upper j bound.
      INTEGER,                      INTENT(in)    :: lower_k             !< Lower k bound.
      INTEGER,                      INTENT(in)    :: upper_k             !< Upper k bound.
      LOGICAL,                      INTENT(in)    :: has_halos           !< .true. if the field has halos

      character(len=150) :: error_msg !< Error message to output

      integer :: nhalos_2 !< 2 times the number of halo points
      integer :: nhalox   !< Number of halos in x
      integer :: nhaloy   !< Number of halos in y

      error_msg = ""
      this%kmin = lower_k
      this%kmax = upper_k
      this%has_halos = has_halos
      this%nhalo_I = 0
      this%nhalo_J = 0
      if (has_halos) then
         !upper_i-lower_i+1 is the size of the compute domain
         !ubound(field_data,1) is the size of the data domain
         nhalos_2 = ubound(field_data,1)-(upper_i-lower_i+1)
         if (mod(nhalos_2, 2) .ne. 0) then
           error_msg = "There are non-symmetric halos in the first dimension"
           return
         endif
         nhalox = nhalos_2/2
         this%nhalo_I = nhalox

         nhalos_2 = ubound(field_data,2)-(upper_j-lower_j + 1)
         if (mod(nhalos_2, 2) .ne. 0) then
           error_msg = "There are non-symmetric halos in the second dimension"
           return
         endif
         nhaloy = nhalos_2/2
         this%nhalo_J = nhaloy

         this%imin = 1 + nhalox
         this%imax = ubound(field_data,1) - nhalox
         this%jmin = 1 + nhaloy
         this%jmax = ubound(field_data,2) - nhaloy
      else
         this%imin = lower_i
         this%imax = upper_i
         this%jmin = lower_j
         this%jmax = upper_j
      endif

   end function set_bounds
   !> @brief Reset the instance bounding box with the bounds determined from the
   !! first three dimensions of the 5D "array" argument
   SUBROUTINE reset_bounds_from_array_4D(this, array)
      CLASS (fmsDiagIbounds_type), INTENT(inout) :: this !< The instance of the bounding box.
      class(*), INTENT( in), DIMENSION(:,:,:,:) :: array !< The 4D input array.
      this%imin = LBOUND(array,1)
      this%imax = UBOUND(array,1)
      this%jmin = LBOUND(array,2)
      this%jmax = UBOUND(array,2)
      this%kmin = LBOUND(array,3)
      this%kmax = UBOUND(array,3)

      this%has_halos = .false.
      this%nhalo_I = 0
      this%nhalo_J = 0
   END SUBROUTINE  reset_bounds_from_array_4D

   !> @brief Reset the instance bounding box with the bounds determined from the
   !! first three dimensions of the 5D "array" argument
   SUBROUTINE reset_bounds_from_array_5D(this, array)
      CLASS (fmsDiagIbounds_type), INTENT(inout) :: this !< The instance of the bounding box.
      CLASS(*), INTENT( in), DIMENSION(:,:,:,:,:) :: array !< The 5D input array.
      this%imin = LBOUND(array,1)
      this%imax = UBOUND(array,1)
      this%jmin = LBOUND(array,2)
      this%jmax = UBOUND(array,2)
      this%kmin = LBOUND(array,3)
      this%kmax = UBOUND(array,3)
   END SUBROUTINE  reset_bounds_from_array_5D

  !> @brief Updates indices based on presence/absence of input indices is, js, ks, ie, je, and ke.
  ! Computes halo sizes in the I and J dimensions.
  ! This routine is intended to be used in diag manager.
  !> @return .false. if there is no error else .true.
  function recondition_indices(indices, field, is_in, js_in, ks_in, &
   ie_in, je_in, ke_in, err_msg) result(ierr)
   type(fmsDiagBoundsHalos_type), intent(inout) :: indices !< Stores indices in order:
                                             !! (/is, js, ks, ie, je, ke, hi, fis, fie, hj, fjs, fje/)
   class(*), intent(in) :: field(:,:,:,:) !< Dummy variable; only the sizes of the first 3 dimensions are used
   integer, intent(in), optional :: is_in, js_in, ks_in, ie_in, je_in, ke_in !< User input indices
   character(len=*), intent(out), optional :: err_msg !< Error message to pass back to caller
   logical :: ierr !< Error flag

   integer :: is, js, ks, ie, je, ke !< Local indices to update
   integer   :: hi !< halo size in the I dimension
   integer   :: hj !< halo size in the J dimension
   integer   :: twohi, twohj !< Temporary storages
   integer   :: fis, fie, fjs, fje !< ! Updated starting and ending indices in the I and J dimensions
   integer :: n1, n2, n3 !< Sizes of the first 3 dimenstions indicies of the data

   ierr = .false.
   if (present(err_msg)) err_msg = ''

   ! If is, js, or ks not present default them to 1
   is = 1
   js = 1
   ks = 1

   IF ( PRESENT(is_in) ) is = is_in
   IF ( PRESENT(js_in) ) js = js_in
   IF ( PRESENT(ks_in) ) ks = ks_in

   n1 = SIZE(field, 1)
   n2 = SIZE(field, 2)
   n3 = SIZE(field, 3)

   ie = is + n1 - 1
   je = js + n2 - 1
   ke = ks + n3 - 1

   IF ( PRESENT(ie_in) ) ie = ie_in
   IF ( PRESENT(je_in) ) je = je_in
   IF ( PRESENT(ke_in) ) ke = ke_in

   twohi = n1 - (ie - is + 1)
   IF ( MOD(twohi, 2) /= 0 ) THEN
     ierr = fms_error_handler('diag_util_mod:recondition_indices', &
       'non-symmetric halos in first dimension', err_msg)
     IF (ierr) RETURN
   END IF

   twohj = n2 - (je - js + 1)
   IF ( MOD(twohj, 2) /= 0 ) THEN
     ierr = fms_error_handler('diag_util_mod:recondition_indices', &
       'non-symmetric halos in second dimension', err_msg)
     IF (ierr) RETURN
   END IF

   hi = twohi/2
   hj = twohj/2

   ! The next line is necessary to ensure that is, ie, js, ie are relative to field(1:,1:)
   ! But this works only when there is no windowing.
   IF ( PRESENT(ie_in) .AND. PRESENT(je_in) ) THEN
     is = 1 + hi
     ie = n1 - hi
     js = 1 + hj
     je = n2 - hj
   END IF

   ! Used for field, mask and rmask bounds
   fis = 1 + hi
   fie = n1 - hi
   fjs = 1 + hj
   fje = n2 - hj

   ! Update indices
   indices%bounds3D%imin = is
   indices%bounds3D%imax = ie
   indices%bounds3D%jmin = js
   indices%bounds3D%jmax = je
   indices%bounds3D%kmin = ks
   indices%bounds3D%kmax = ke
   indices%hi = hi
   indices%hj = hj
   indices%fis = fis
   indices%fie = fie
   indices%fjs = fjs
   indices%fje = fje
 end function recondition_indices

 !> @brief Rebase the ouput bounds for a given dimension based on the starting and ending indices of
 !! a subregion. This is for when blocking is used.
 subroutine rebase_output(bounds_out, starting, ending, dim)
   CLASS (fmsDiagIbounds_type), INTENT(inout) :: bounds_out !< Bounds to rebase
   integer,                     intent(in)    :: starting   !< Starting index of the dimension
   integer,                     intent(in)    :: ending     !< Ending index of the dimension
   integer,                     intent(in)    :: dim        !< Dimension to update

   !> The starting index is going to be either "starting" if only a section of the
   !! block is in the subregion or bounds_out%[]min if the whole section of the block is in the
   !! subregion. The -starting+1 s needed so that indices start as 1 since the output buffer has
   !! indices 1:size of a subregion

   !> The ending index is going to be either bounds_out%[]max if the whole section of the block
   !! is in the subregion or bounds_out%[]min + size of the subregion if only a section of the
   !! block is in the susbregion
   select case (dim)
   case (xdimension)
      bounds_out%imin = max(starting, bounds_out%imin)-starting+1
      bounds_out%imax = min(bounds_out%imax, bounds_out%imin + ending-starting)
   case (ydimension)
      bounds_out%jmin =  max(starting, bounds_out%jmin)-starting+1
      bounds_out%jmax = min(bounds_out%jmax, bounds_out%jmin + ending-starting)
   case (zdimension)
      bounds_out%kmin =max(starting, bounds_out%kmin)-starting+1
      bounds_out%kmax = min(bounds_out%kmax, bounds_out%kmin + ending-starting)
   end select
 end subroutine

 !> @brief Rebase the input bounds for a given dimension based on the starting and ending indices
 !! of a subregion. This is for when blocking is used
 subroutine rebase_input(bounds_in, bounds, starting, ending, dim)
   CLASS (fmsDiagIbounds_type), INTENT(inout) :: bounds_in  !< Bounds to rebase
   CLASS (fmsDiagIbounds_type), INTENT(in)    :: bounds     !< Original indices (i.e is_in, ie_in,
                                                            !! passed into diag_manager)
   integer,                     intent(in)    :: starting   !< Starting index of the dimension
   integer,                     intent(in)    :: ending     !< Ending index of the dimension
   integer,                     intent(in)    :: dim        !< Dimension to update

   !> The starting index is going to be either "starting" if only a section of the
   !! block is in the subregion or starting-bounds%imin+1 if the whole section of the block is in the
   !! subregion.

   !> The ending index is going to be either bounds_out%[]max if the whole section of the block
   !! is in the subregion or bounds%[]min + size of the subregion if only a section of the
   !! block is in the susbregion
   select case (dim)
   case (xdimension)
      bounds_in%imin = min(abs(starting-bounds%imin+1), starting)
      bounds_in%imax = min(bounds_in%imax, (bounds_in%imin + ending-starting))
   case (ydimension)
      bounds_in%jmin = min(abs(starting-bounds%jmin+1), starting)
      bounds_in%jmax = min(bounds_in%jmax, (bounds_in%jmin + ending-starting))
   case (zdimension)
      bounds_in%kmin = min(abs(starting-bounds%kmin+1), starting)
      bounds_in%kmax = min(bounds_in%kmax, (bounds_in%kmin + ending-starting))
   end select
 end subroutine

  END MODULE fms_diag_bbox_mod
  !> @}
  ! close documentation grouping
