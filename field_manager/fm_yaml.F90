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

!> @defgroup fm_yaml_mod fm_yaml_mod
!> @ingroup fm_yaml
!> @brief Reads entries from a field table yaml into a
!! nested object for use in the field manager.
!!
!> @author Eric Stofferahn
!!

!> @file
!> @brief File for @ref fm_yaml_mod

!> @addtogroup fm_yaml_mod
!> @{
module fm_yaml_mod
#ifdef use_yaml

use yaml_parser_mod
use mpp_mod, only: mpp_error, fatal
implicit none
private

!> @}

public :: build_fmTable

!> @brief This type represents a subparameter block for a given variable parameter.
!> This type contains the name of the associated parameter and the subparameter key/value pairs
!> @ingroup fm_yaml_mod
type, public :: fmAttr_t
  integer                                     :: id                    !< block id of this var
  character(len=:), allocatable               :: paramname             !< name of associated parameter
  character(len=:), dimension(:), allocatable :: keys                  !< name of the attribute
  character(len=:), dimension(:), allocatable :: values                !< value of the attribute
end type fmAttr_t

!> @brief This type represents the entries for a given variable, e.g. dust.
!> This type contains the name of the variable, the block id, the key/value pairs for the
!> variable's parameters, and any applicable subparameters
!> @ingroup fm_yaml_mod
type, public :: fmVar_t
  integer                                     :: id                    !< block id of this var
  character(len=:), allocatable               :: name                  !< name of the variable
  character(len=:), dimension(:), allocatable :: keys                  !< names of params
  character(len=:), dimension(:), allocatable :: values                !< values of params
  type (fmAttr_t), allocatable                :: attributes(:)         !< attributes in this var
end type fmVar_t

!> @brief This type represents the entries for a given model, e.g. land, ocean, atmosphere.
!> This type contains the name of the model, the block id, and the variables within this model
!> @ingroup fm_yaml_mod
type, public :: fmModel_t
  integer                       :: id                  !< block id of this model
  character(len=:), allocatable :: name                !< name of the model
  type (fmVar_t), allocatable   :: variables(:)        !< variables in this model
end type fmModel_t

!> @brief This type represents the entries for a specific field type, e.g. a tracer.
!> This type contains the name of the field type, the block id, and the models within this field type
!> @ingroup fm_yaml_mod
type, public :: fmType_t
  integer                       :: id                  !< block id of this type
  character(len=:), allocatable :: name                !< name of the type
  type (fmModel_t), allocatable :: models(:)           !< models in this type
end type fmType_t

!> @brief This type contains the field types within a field table.
!> @ingroup fm_yaml_mod
type, public :: fmTable_t
  type (fmType_t), allocatable :: types(:) !< field types in this table
end type fmTable_t

contains

!> @addtogroup fm_yaml_mod
!> @{

!> @brief Subroutine to populate an fmTable by reading a yaml file, given an optional filename.
subroutine build_fmTable(fmTable, filename)
  type(fmTable_t), intent(out)           :: fmTable  !< the field table
  character(len=*), intent(in), optional :: filename !< the name of the yaml file
  integer                                :: yfid     !< file id of the yaml file
  integer                                :: ntypes   !< number of field types attached to this table
  integer                                :: i        !< Loop counter

  if (.not. present(filename)) then
    yfid = open_and_parse_file("field_table.yaml")
  else
    yfid = open_and_parse_file(trim(filename))
  endif

  ntypes = get_num_blocks(yfid, "field_table", 0)
  allocate(fmTable%types(ntypes))

  ! Gets the block ids for the associated types of fmTable.
  call get_block_ids(yfid, "field_table", fmTable%types(:)%id)

  do i=1,ntypes
    call build_fmType(fmTable%types(i), yfid)
  enddo
end subroutine build_fmTable

!> @brief Populates an fmType, which is assumed to already have its `id` parameter set.
subroutine build_fmType(fmType, yfid)
  type(fmType_t), intent(inout) :: fmType    !< type object
  integer, intent(in)           :: yfid      !< file id of the yaml file
  integer, dimension(1)         :: key_ids   !< array of key ids
  character(len=256)            :: key_name  !< the name of a key
  character(len=256)            :: key_value !< the value of a key
  integer                       :: nmodels   !< number of models attached to this type
  integer                       :: i         !< Loop counter

  nmodels = get_num_blocks(yfid, "modlist", fmType%id)
  allocate(fmType%models(nmodels))

  ! Gets the block ids for the associated models of fmType.
  call get_block_ids(yfid, "modlist", fmType%models(:)%id, fmType%id)

  if (get_nkeys(yfid, fmType%id).ne.1) then
    call mpp_error(FATAL, "fm_yaml_mod: A single `field_type` key is expected")
  endif

  call get_key_ids(yfid, fmType%id, key_ids)
  call get_key_name(yfid, key_ids(1), key_name)
  call get_key_value(yfid, key_ids(1), key_value)

  if (trim(key_name).ne."field_type") then
    call mpp_error(FATAL, "fm_yaml_mod: A single `field_type` key is expected")
  endif

  fmType%name = trim(key_value)

  do i=1,nmodels
    call build_fmModel(fmType%models(i), yfid)
  enddo
end subroutine build_fmType

!> @brief Populates an fmModel, which is assumed to already have its `id` parameter set.
subroutine build_fmModel(fmModel, yfid)
  type(fmModel_t), intent(inout) :: fmModel   !< model object
  integer, intent(in)            :: yfid      !< file id of the yaml file
  integer, dimension(1)          :: key_ids   !< array of key ids
  character(len=256)             :: key_name  !< the name of a key
  character(len=256)             :: key_value !< the value of a key
  integer                        :: nvars     !< number of variables attached to this model
  integer                        :: i         !< Loop counter

  nvars = get_num_blocks(yfid, "varlist", fmModel%id)
  allocate(fmModel%variables(nvars))

  ! gets the block ids for the associated variables of fmModel.
  call get_block_ids(yfid, "varlist", fmModel%variables(:)%id, fmModel%id)

  if (get_nkeys(yfid, fmModel%id).ne.1) then
    call mpp_error(FATAL, "fm_yaml_mod: A single `model_type` key is expected")
  endif

  call get_key_ids(yfid, fmModel%id, key_ids)
  call get_key_name(yfid, key_ids(1), key_name)
  call get_key_value(yfid, key_ids(1), key_value)

  if (trim(key_name).ne."model_type") then
    call mpp_error(FATAL, "fm_yaml_mod: A single `model_type` key is expected")
  endif

  fmModel%name = trim(key_value)

  do i=1,nvars
    call build_fmVar(fmModel%variables(i), yfid)
  enddo
end subroutine build_fmModel

!> @brief Populates an fmVar and creates any associated fmAttrs
subroutine build_fmVar(fmVar, yfid)
  type(fmVar_t), intent(inout) :: fmVar                 !< variable object
  integer, intent(in)          :: yfid                  !< file id of the yaml file
  integer                      :: nkeys                 !< number of keys defined for this var
  integer, allocatable         :: key_ids(:)            !< array of key ids
  character(len=256)           :: key_name              !< the name of a key
  character(len=256)           :: key_value             !< the value of a key
  integer                      :: nattrs                !< number of attribute blocks attached to this var
  integer                      :: nmethods              !< total number of methods attached to this var
  integer                      :: maxln                 !< max string length of method names
  integer                      :: maxlv                 !< max string length of method values
  character(:), allocatable    :: attr_method_keys(:)   !< Keys of methods defined in attribute blocks
  character(:), allocatable    :: attr_method_values(:) !< Values of methods defined in attribute blocks
  integer                      :: i_name                !< Index of the key containing the variable's name
  integer                      :: i, j                  !< Loop indices

  ! Read attribute blocks attached to this variable
  call fmVar_read_attrs(fmVar, yfid, attr_method_keys, attr_method_values)
  nattrs = size(attr_method_keys)

  nkeys = get_nkeys(yfid, fmVar%id)
  allocate(key_ids(nkeys))
  call get_key_ids(yfid, fmVar%id, key_ids)

  maxln = len(attr_method_keys)
  maxlv = len(attr_method_values)
  i_name = -1

  do i=1,nkeys
    call get_key_name(yfid, key_ids(i), key_name)
    call get_key_value(yfid, key_ids(i), key_value)

    if (trim(key_name) .eq. "variable") then
      if (i_name .ne. -1) then
        call mpp_error(FATAL, "fm_yaml_mod: A variable can have only one `variable` key")
      endif

      fmVar%name = trim(key_value)
      i_name = i
    else
      maxln = max(maxln, len_trim(key_name))
      maxlv = max(maxlv, len_trim(key_value))
    endif
  enddo

  if (i_name .eq. -1) then
    call mpp_error(FATAL, "fm_yaml_mod: Every variable must have a `variable` key")
  endif

  ! Number of methods is the number of keys (excluding `variable`), plus one for each attribute block.
  nmethods = nkeys - 1 + nattrs

  allocate(character(len=maxln)::fmVar%keys(nmethods))
  allocate(character(len=maxlv)::fmVar%values(nmethods))

  j = 1
  do i=1,nkeys
    if (i.eq.i_name) cycle ! Exclude `variable` key

    call get_key_name(yfid, key_ids(i), key_name)
    call get_key_value(yfid, key_ids(i), key_value)
    fmVar%keys(j) = trim(key_name)
    fmVar%values(j) = trim(key_value)

    j = j + 1
  enddo

  ! Add methods defined within attribute blocks.
  fmVar%keys(j:) = attr_method_keys
  fmVar%values(j:) = attr_method_values
end subroutine build_fmVar

!> @brief Reads the attribute blocks attached to a variable and populates the associated fmAttr structures.
!! Returns two arrays containing key/value pairs of all methods defined via attribute blocks.
subroutine fmVar_read_attrs(fmVar, yfid, method_keys, method_values)
  type(fmVar_t), intent(inout)           :: fmVar            !< variable object
  integer, intent(in)                    :: yfid             !< file id of the yaml file
  character(:), allocatable, intent(out) :: method_keys(:)   !< Method keys (names of attribute blocks)
  character(:), allocatable, intent(out) :: method_values(:) !< Method values from attribute blocks
  integer                                :: nattrs           !< number of attribute blocks
  integer                                :: nkeys            !< number of keys in an attribute block
  integer, allocatable                   :: key_ids(:)       !< array of key ids
  character(len=256)                     :: key_name         !< the name of a key
  character(len=256)                     :: key_value        !< the value of a key
  integer                                :: maxln_m          !< max string length of method names
  integer                                :: maxlv_m          !< max string length of method values
  integer                                :: maxln_a          !< max string length of subparameter names
  integer                                :: maxlv_a          !< max string length of subparameter values
  integer,allocatable                    :: name_key_id(:)   !< Indices of attribute `value` keys
  integer                                :: i, j, k          !< Loop counters

  nattrs = get_num_unique_blocks(yfid, fmVar%id)
  allocate(fmVar%attributes(nattrs))
  allocate(name_key_id(nattrs))

  ! gets the block ids for the associated attributes of fmVar.
  call get_unique_block_ids(yfid, fmVar%attributes(:)%id, fmVar%id)

  maxln_m = 0
  maxlv_m = 0
  name_key_id = -1

  do i=1,nattrs
    associate (fmAttr => fmVar%attributes(i))
      call get_block_name(yfid, fmAttr%id, key_value)
      fmAttr%paramname = trim(key_value)

      nkeys = get_nkeys(yfid, fmAttr%id)
      allocate(key_ids(nkeys))
      call get_key_ids(yfid, fmAttr%id, key_ids)

      maxln_a = 0
      maxlv_a = 0

      do j=1,nkeys
        call get_key_name(yfid, key_ids(j), key_name)
        call get_key_value(yfid, key_ids(j), key_value)

        if (trim(key_name) .eq. "value") then
          if (name_key_id(i) .ne. -1) then
            call mpp_error(FATAL, "fm_yaml_mod: A variable attribute block can only have one `value` key")
          endif

          maxln_m = max(maxln_m, len(fmAttr%paramname))
          maxlv_m = max(maxlv_m, len_trim(key_value))

          name_key_id(i) = key_ids(j)
        else
          maxln_a = max(maxln_a, len_trim(key_name))
          maxlv_a = max(maxlv_a, len_trim(key_value))
        endif
      enddo

      if (name_key_id(i) .eq. -1) then
        call mpp_error(FATAL, "fm_yaml_mod: Every variable attribute must have a `value` key")
      endif

      allocate(character(len=maxln_a)::fmAttr%keys(nkeys - 1))
      allocate(character(len=maxlv_a)::fmAttr%values(nkeys - 1))

      k = 1
      do j=1,nkeys
        if (key_ids(j).eq.name_key_id(i)) cycle

        call get_key_name(yfid, key_ids(j), key_name)
        call get_key_value(yfid, key_ids(j), key_value)
        fmAttr%keys(k) = trim(key_name)
        fmAttr%values(k) = trim(key_value)

        k = k + 1
      enddo

      deallocate(key_ids)
    end associate
  enddo

  allocate(character(len=maxln_m)::method_keys(nattrs))
  allocate(character(len=maxlv_m)::method_values(nattrs))

  do i=1,nattrs
    method_keys(i) = fmVar%attributes(i)%paramname
    call get_key_value(yfid, name_key_id(i), method_values(i))
  enddo
end subroutine fmVar_read_attrs

#endif
end module fm_yaml_mod

!> @}
! close documentation grouping
