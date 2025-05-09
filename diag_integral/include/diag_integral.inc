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
!> @defgroup diag_integral_mod diag_integral_mod
!> @ingroup diag_integral
!!
!! @author Fei Liu <Fei.Liu@noaa.gov>
!!
!! @brief This module computes and outputs global and / or hemispheric physics
!!        integrals.
!###############################################################################
!> @brief Perform a 2 dimensional summation of named field
!!
!! @implements sum_diag_integral_field
!!
!! <b> Template: </b>
!!
!! @code{.f90}
!! call sum_field_2d (name, field_data2d, is, js)
!! @endcode
!!
!! <b> Parameters: </b>
!!
!! @code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: field_data2d(:,:)
!! integer, optional, intent(in) :: is, js
!! @endcode
!!
!! @param [in] <name> Name of the field to be integrated
!! @param [in] <field_data2d> field of integrands to be summed over
!! @param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine SUM_FIELD_2D_ (name, field_data2d, is, js)

character(len=*),  intent(in) :: name !< Name of the field to be integrated
real(FMS_DI_KIND_), intent(in) :: field_data2d(:,:) !< field of integrands to be summed over
integer, optional, intent(in) :: is !< starting i indices over which summation is to occur
integer, optional, intent(in) :: js !< starting j indices over which summation is to occur

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer :: field           !< index of desired integral
      integer :: i1  !< location indices of current data in
                                 !! processor-global coordinates
      integer :: j1  !< location indices of current data in
                                 !! processor-global coordinates
      integer :: i2  !< location indices of current data in
                                 !! processor-global coordinates
      integer :: j2  !< location indices of current data in
                                 !! processor-global coordinates


!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)  then
        call error_mesg ('diag_integral_mod', &
                                    'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
     i1 = 1;  if (present(is)) i1 = is
     j1 = 1;  if (present(js)) j1 = js
     i2 = i1 + size(field_data2d,1) - 1
     j2 = j1 + size(field_data2d,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral and add the
!    values at this set of grid points to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count(field) +   &
                            size(field_data2d,1)*size(field_data2d,2)
      field_sum   (field) = field_sum   (field) +  &
                            sum (real(field_data2d,r8_kind) * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine SUM_FIELD_2D_



!###############################################################################
!> @brief Perform a 3 dimensional summation of named field
!!
!! @implements sum_diag_integral_field
!!
!! <b> Template: </b>
!!
!! @code{.f90}
!! call sum_field_3d (name, field_data3d, is, js)
!! @endcode
!!
!! <b> Parameters: </b>
!!
!! @code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: field_data3d(:,:,:)
!! integer, optional, intent(in) :: is, js
!! @endcode
!!
!! @param [in] <name> Name of the field to be integrated
!! @param [in] <field_data3d> field of integrands to be summed over
!! @param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine SUM_FIELD_3D_ (name, field_data3d, is, js)

character(len=*),  intent(in) :: name !< Name of the field to be integrated
real(FMS_DI_KIND_),  intent(in) :: field_data3d(:,:,:) !< field of integrands to be summed over
integer, optional, intent(in) :: is !< starting i,j indices over which summation is to occur
integer, optional, intent(in) :: js !< starting i,j indices over which summation is to occur

!-------------------------------------------------------------------------------
! local variables:
!     data2
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
      real(r8_kind), dimension (size(field_data3d,1),  &
                       size(field_data3d,2)) :: data2

      integer :: field !< index of desired integral
      integer :: i1 !< location indices of current data in
                                !! processor-global coordinates
      integer :: j1 !< location indices of current data in
                                !! processor-global coordinates
      integer :: i2 !< location indices of current data in
                                !! processor-global coordinates
      integer :: j2 !< location indices of current data in
                                !! processor-global coordinates


!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)   then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
      i1 = 1;  if (present(is)) i1 = is
      j1 = 1;  if (present(js)) j1 = js
      i2 = i1 + size(field_data3d,1) - 1
      j2 = j1 + size(field_data3d,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. sum first
!    in the vertical and then add the values at this set of grid points
!    to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) +   &
                            size(field_data3d,1)*size(field_data3d,2)
      data2 = sum(real(field_data3d,r8_kind),3)
      field_sum   (field) = field_sum   (field) +  &
                            sum (data2 * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine SUM_FIELD_3D_



!###############################################################################
!> @brief Perform a 3 dimensional weighted summation of named field
!!
!! @implements sum_diag_integral_field
!!
!! <b> Template: </b>
!!
!! @code{.f90}
!! call sum_field_wght_3d (name, field_data3d, wt, is, js)
!! @endcode
!!
!! <b> Parameters: </b>
!!
!! @code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: field_data3d(:,:,:), wt(:,:,:)
!! integer, optional, intent(in) :: is, js
!! @endcode
!!
!! @param [in] <name> Name of the field to be integrated
!! @param [in] <field_data3d> field of integrands to be summed over
!! @param [in] <wt> the weight function to be evaluated at summation
!! @param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine SUM_FIELD_WGHT_3D_ (name, field_data3d, wt, is, js)

character(len=*),   intent(in) :: name !< Name of the field to be integrated
real(FMS_DI_KIND_), intent(in) :: field_data3d(:,:,:) !< field of integrands to be summed over
real(FMS_DI_KIND_), intent(in) :: wt(:,:,:) !< the weight function to be evaluated at summation
integer, optional,  intent(in) :: is !< starting i indices over which summation is to occur
integer, optional,  intent(in) :: js !< starting j indices over which summation is to occur

!-------------------------------------------------------------------------------
! local variables:
!     data2
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
      real(r8_kind), dimension (size(field_data3d,1),size(field_data3d,2)) :: data2
      integer :: field !< index of desired integral
      integer :: i1 !< location indices of current data in
                                       !! processor-global coordinates
      integer :: j1 !< location indices of current data in
                                       !! processor-global coordinates
      integer :: i2 !< location indices of current data in
                                       !! processor-global coordinates
      integer :: j2 !< location indices of current data in
                                       !! processor-global coordinates

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)   then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
      i1 = 1;  if (present(is)) i1 = is
      j1 = 1;  if (present(js)) j1 = js
      i2 = i1 + size(field_data3d,1) - 1
      j2 = j1 + size(field_data3d,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. sum first
!    in the vertical (including a vertical weighting factor) and then
!    add the values at this set of grid points to the accumulation
!    array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) +   &
                            size(field_data3d,1)*size(field_data3d,2)
      data2 = vert_diag_integral (real(field_data3d,r8_kind), real(wt,r8_kind))
      field_sum(field) = field_sum   (field) +  &
                         sum (data2 * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine SUM_FIELD_WGHT_3D_



!###############################################################################
!> @brief Perform a 2 dimensional hemispherical summation of named field
!!
!! @implements sum_diag_integral_field
!!
!! <b> Template: </b>
!!
!! @code{.f90}
!! call sum_field_2d_hemi (name, field_data2d, is, ie, js, je)
!! @endcode
!!
!! <b> Parameters: </b>
!!
!! @code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: field_data2d(:,:)
!! integer,           intent(in) :: is, js, ie, je
!! @endcode
!!
!! @param [in] <name> Name of the field to be integrated
!! @param [in] <field_data2d> field of integrands to be summed over
!! @param [in] <is, js, ie, je> starting/ending i,j indices over which summation
!!        is to occur
!!
subroutine SUM_FIELD_2D_HEMI_ (name, field_data2d, is, ie, js, je)

character(len=*),  intent(in) :: name !< Name of the field to be integrated
real(FMS_DI_KIND_),intent(in) :: field_data2d(:,:) !< field of integrands to be summed over
integer,           intent(in) :: is !< starting/ending i,j indices over which summation
                                                !! is to occur
integer,           intent(in) :: js !< starting/ending i,j indices over which summation
                                                !! is to occur
integer,           intent(in) :: ie !< starting/ending i,j indices over which summation
                                                !! is to occur
integer,           intent(in) :: je !< starting/ending i,j indices over which summation
                                                !! is to occur

!-------------------------------------------------------------------------------
! local variables:
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
   integer :: field !< index of desired integral
   integer :: i1 !< location indices of current data in
                             !! processor-global coordinates
   integer :: j1 !< location indices of current data in
                             !! processor-global coordinates
   integer :: i2 !< location indices of current data in
                             !! processor-global coordinates
   integer :: j2 !< location indices of current data in
                             !! processor-global coordinates

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)    then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!    define the processor-global indices of the current data. this form
!    is needed to handle case of 2d domain decomposition with physics
!    window smaller than processor domain size.
!-------------------------------------------------------------------------------
      i1 = mod ( (is-1), size(field_data2d,1) ) + 1
      i2 = i1 + size(field_data2d,1) - 1

!-------------------------------------------------------------------------------
!    for a hemispheric sum, sum one jrow at a time in case a processor
!    has data from both hemispheres.
!-------------------------------------------------------------------------------
      j1 = mod ( (js-1) ,size(field_data2d,2) ) + 1
      j2 = j1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. include hemi-
!    spheric factor of 2 in field_count. add the data values at this
!    set of grid points to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) + 2* (i2-i1+1)*(j2-j1+1)
      field_sum   (field) = field_sum   (field) +  &
                            sum ( real(field_data2d(i1:i2,j1:j2),r8_kind) * area(is:ie,js:je))

!$OMP END CRITICAL

end subroutine SUM_FIELD_2D_HEMI_


!> @}
! close documentation grouping
