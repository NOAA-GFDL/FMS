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
!> @defgroup drifters_io_mod drifters_io_mod
!> @ingroup drifters
!> @brief Saves drifter data for postprocessing and restarts

!> @file
!> @brief File for @ref drifters_io_mod

!> @addtogroup drifters_io_mod
!> @{
module drifters_io_mod

  use netcdf
  use netcdf_nf_data
  use netcdf_nf_interfaces
  use netcdf4_nf_interfaces

  implicit none
  private

  public :: drifters_io_type, drifters_io_new, drifters_io_del, drifters_io_set_time_units
  public :: drifters_io_set_position_names, drifters_io_set_position_units, drifters_io_set_field_names
  public :: drifters_io_set_field_units, drifters_io_write

  ! Globals
  integer, parameter, private   :: MAX_STR_LEN = 128
  ! Include variable "version" to be written to log file.
#include<file_version.h>

  real :: drfts_eps_t = 10.*epsilon(1.)

!> @}
  !> @brief IO data for drifters.
  !> @ingroup drifters_input_mod
  type drifters_io_type
     real                 :: time
     integer              :: it !< time index
     integer              :: it_id !< infinite axis index
     integer              :: ncid
     integer              :: nc_positions, nc_fields, nc_ids, nc_time, nc_index_time
     logical              :: enddef
  end type drifters_io_type
!> @addtogroup drifters_io_mod
!> @{
contains

!###############################################################################
  subroutine drifters_io_new(self, filename, nd, nf, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: nd  !< number of dims
    integer, intent(in)           :: nf  !< number of fields
    character(len=*), intent(out) :: ermesg

    integer ier, nc_it_id, nc_nd, nc_nf
    integer :: size1(1), size2(2)

    ermesg=''
    self%enddef = .FALSE.

    ier = nf_create(filename, NF_CLOBBER, self%ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_create ('//filename//') '//nf_strerror(ier)

    ! global attributes
    ier = nf_put_att_text(self%ncid, NF_GLOBAL, 'version', len_trim(version), trim(version))


    ! dimensions
    ier = nf_def_dim(self%ncid, 'np', NF_UNLIMITED, nc_it_id)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (it_id) '//nf_strerror(ier)

    ier = nf_def_dim(self%ncid, 'nf', nf, nc_nf)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (nf) '//nf_strerror(ier)

    ier = nf_def_dim(self%ncid, 'nd', nd, nc_nd)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_dim (nd) '//nf_strerror(ier)

    ! variables
    size1 = (/nc_it_id/)
    ier = nf_def_var(self%ncid, 'index_time', NF_INT, 1, size1, self%nc_index_time)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (index_time)'//nf_strerror(ier)

    ier = nf_def_var(self%ncid, 'time', NF_DOUBLE, 1, size1, self%nc_time)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (time)'//nf_strerror(ier)

    ier = nf_def_var(self%ncid, 'ids', NF_INT, 1, size1, self%nc_ids)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (ids)'//nf_strerror(ier)

    size2 = (/nc_nd, nc_it_id/)
    ier = nf_def_var(self%ncid, 'positions', NF_DOUBLE, 2, size2, self%nc_positions)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (positions)'//nf_strerror(ier)

    size2 = (/nc_nf, nc_it_id/)
    ier = nf_def_var(self%ncid, 'fields', NF_DOUBLE, 2, size2, self%nc_fields)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_new::nf_def_var (fields)'//nf_strerror(ier)

    self%time  = -huge(1.)
    self%it    = -1
    self%it_id = 1

  end subroutine drifters_io_new

!###############################################################################
  subroutine drifters_io_del(self, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(out) :: ermesg

    integer ier

    ermesg = ''

    ier = nf_close(self%ncid)
    if(ier/=NF_NOERR) ermesg = 'drifters_io_del::nf_close '//nf_strerror(ier)

  end subroutine drifters_io_del

!###############################################################################
  subroutine drifters_io_set_time_units(self, name, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: ermesg

    integer ier

    ermesg = ''
    ier = nf_put_att_text(self%ncid, NF_GLOBAL, &
         & 'time_units', len_trim(name), trim(name))
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_set_time_units::failed to add time_units attribute ' &
         & //nf_strerror(ier)

  end subroutine drifters_io_set_time_units

!###############################################################################
  subroutine drifters_io_set_position_names(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'name_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_positions, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_position_names::failed to add name attribute to positions '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_position_names

!###############################################################################
  subroutine drifters_io_set_position_units(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'unit_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_positions, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_position_names::failed to add unit attribute to positions '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_position_units

!###############################################################################
  subroutine drifters_io_set_field_names(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'name_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_fields, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_field_names::failed to add name attribute to fields '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_field_names

!###############################################################################
  subroutine drifters_io_set_field_units(self, names, ermesg)
    type(drifters_io_type)        :: self
    character(len=*), intent(in)  :: names(:)
    character(len=*), intent(out) :: ermesg

    integer n, ier, i
    character(len=128) :: attname

    n = size(names)
    ermesg = ''

    do i = 1, n
       write(attname, '(i6)' ) i
       attname = 'unit_'//adjustl(attname)
       ier = nf_put_att_text(self%ncid, self%nc_fields, &
            & trim(attname), len_trim(names(i)), trim(names(i)))
       if(ier/=NF_NOERR) &
            & ermesg = 'drifters_io_set_field_units::failed to add unit attribute to fields '//nf_strerror(ier)
    enddo

  end subroutine drifters_io_set_field_units
!###############################################################################

  subroutine drifters_io_write(self, time, np, nd, nf, ids, positions, fields, ermesg)
    type(drifters_io_type)        :: self
    real, intent(in)              :: time
    integer, intent(in)           :: np    !< number of dirfters
    integer, intent(in)           :: nd    !< number of dimensions
    integer, intent(in)           :: nf    !< number of fields
    integer, intent(in)           :: ids(np)          !< of size np
    real, intent(in)              :: positions(nd,np) !< nd times np
    real, intent(in)              :: fields(nf,np)    !< nf times np
    character(len=*), intent(out) :: ermesg

    integer ier, i
    integer :: start1(1), len1(1), start2(2), len2(2)
    integer :: it_indices(np)
    real    :: time_array(np)

    ermesg = ''

    if(.not. self%enddef) then
       ier = nf_enddef(self%ncid)
       if(ier/=NF_NOERR) then
            ermesg = 'drifters_io_write::nf_enddef failure. No data will be written. '//nf_strerror(ier)
            return
       endif
       self%enddef = .TRUE.
    endif

    if(abs(time - self%time) > drfts_eps_t) then
       self%it = self%it + 1
       self%time = time
    endif

    start1(1) = self%it_id
    len1(1)   = np

    it_indices = (/(self%it,i=1,np)/)
    ier = nf_put_vara_int( self%ncid, self%nc_index_time, start1, len1, it_indices )
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write index_time: ' //nf_strerror(ier)

    time_array = (/(time,i=1,np)/)
    ier = nf90_put_var( self%ncid, self%nc_time, time_array, start1, len1 )
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write time: ' //nf_strerror(ier)

    ier = nf_put_vara_int(self%ncid, self%nc_ids, start1, len1, ids)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write ids: '//nf_strerror(ier)

    start2(1) = 1
    start2(2) = self%it_id

    len2(1)   = nd
    len2(2)   = np

    ier = nf90_put_var(self%ncid, self%nc_positions, positions, start2, len2)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write positions: '//nf_strerror(ier)

    len2(1)   = nf
    len2(2)   = np

    ier = nf90_put_var(self%ncid, self%nc_fields, fields, start2, len2)
    if(ier/=NF_NOERR) &
         & ermesg = 'drifters_io_write::failed to write fields: '//nf_strerror(ier)

    self%it_id = self%it_id + np

  end subroutine drifters_io_write

end module drifters_io_mod
!> @}
! close documentation grouping
