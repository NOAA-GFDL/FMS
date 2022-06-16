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
implicit none
private

integer :: i, table_i, type_i, model_i, var_i, var_j, attr_j !< counters

!> @}
! close documentation grouping

type, public :: fmTable_t
  integer                       :: yfid                    !< file id of a yaml file
  character(len=11)             :: blockname="field_table" !< name of the root block
  integer                       :: nchildren               !< number of field types
  integer, allocatable          :: child_ids(:)            !< array of type ids
  type (fmType_t), allocatable  :: children(:)             !< models in this type
  contains
    procedure :: get_blocks => get_blocks_fmTable_t
    procedure :: destruct => destruct_fmTable_t
    procedure :: create_children => create_children_fmTable_t
end type fmTable_t

type, public :: fmType_t
  integer                       :: yfid                !< file id of a yaml file
  integer                       :: id                  !< block id of this type
  character(len=:), allocatable :: name                !< name of the type
  character(len=7)              :: blockname="modlist" !< name of the root block
  integer                       :: nchildren           !< number of model types
  integer, allocatable          :: child_ids(:)        !< array of model ids
  type (fmModel_t), allocatable :: children(:)         !< models in this type
  contains
    procedure :: get_blocks => get_blocks_fmType_t
    procedure :: destruct => destruct_fmType_t
    procedure :: get_name => get_name_fmType_t
    procedure :: create_children => create_children_fmType_t
end type fmType_t

type, public :: fmModel_t
  integer                       :: yfid                !< file id of a yaml file
  integer                       :: id                  !< block id of this model
  character(len=:), allocatable :: name                !< name of the model
  character(len=7)              :: blockname="varlist" !< name of the root block
  integer                       :: nchildren           !< number of var types
  integer, allocatable          :: child_ids(:)        !< array of var ids
  type (fmVar_t), allocatable   :: children(:)         !< variables in this model
  contains
    procedure :: get_blocks => get_blocks_fmModel_t
    procedure :: destruct => destruct_fmModel_t
    procedure :: get_name => get_name_fmModel_t
    procedure :: create_children => create_children_fmModel_t
end type fmModel_t

type, public :: fmVar_t
  integer                                     :: yfid                  !< file id of a yaml file
  integer                                     :: id                    !< block id of this var
  character(len=:), allocatable               :: name                  !< name of the variable
  integer, dimension(:), allocatable          :: key_ids               !< key ids for params
  character(len=:), dimension(:), allocatable :: keys                  !< names of params
  character(len=:), dimension(:), allocatable :: values                !< values of params
  character(len=9)                            :: blockname="subparams" !< name of the root block
  integer                                     :: nchildren             !< number of attributes
  integer, allocatable                        :: child_ids(:)          !< array of attribute ids
  type (fmAttr_t), allocatable                :: children(:)           !< attributes in this var
  contains
    procedure :: get_blocks => get_blocks_fmVar_t
    procedure :: destruct => destruct_fmVar_t
    procedure :: get_names_and_props => get_name_fmVar_t
    procedure :: create_children => create_children_fmVar_t
end type fmVar_t

type, public :: fmAttr_t
  integer                                     :: yfid                  !< file id of a yaml file
  integer                                     :: id                    !< block id of this var
  character(len=:), allocatable               :: paramname             !< name of associated parameter  
  character(len=:), dimension(:), allocatable :: keys                  !< name of the variable
  character(len=:), dimension(:), allocatable :: values                !< name of the variable
  contains
    procedure :: destruct => destruct_fmAttr_t
    procedure :: get_names_and_props => get_name_fmAttr_t
end type fmAttr_t

interface fmTable_t
  module procedure construct_fmTable_t
end interface fmTable_t

interface fmType_t
  module procedure construct_fmType_t
end interface fmType_t

interface fmModel_t
  module procedure construct_fmModel_t
end interface fmModel_t

interface fmVar_t
  module procedure construct_fmVar_t
end interface fmVar_t

interface fmAttr_t
  module procedure construct_fmAttr_t
end interface fmAttr_t

contains

function construct_fmTable_t(filename) result(this)
  type (fmTable_t)                       :: this     !< the field table
  character(len=*), intent(in), optional :: filename !< the name of the yaml file

  if (.not. present(filename)) then
    this%yfid = open_and_parse_file("field_table.yaml")
  else
    this%yfid = open_and_parse_file(trim(filename))
  endif
  this%nchildren = get_num_blocks(this%yfid, this%blockname)
  allocate(this%child_ids(this%nchildren))
end function construct_fmTable_t

function construct_fmType_t(in_yfid, in_id) result(this)
  type (fmType_t)     :: this    !< the type object
  integer, intent(in) :: in_yfid !< yaml file id
  integer, intent(in) :: in_id   !< block_id of type from parent

  this%yfid = in_yfid
  this%id = in_id
  this%nchildren = get_num_blocks(this%yfid, this%blockname, this%id)
  allocate(this%child_ids(this%nchildren))
end function construct_fmType_t

function construct_fmModel_t(in_yfid, in_id) result(this)
  type (fmModel_t)    :: this    !< the model object
  integer, intent(in) :: in_yfid !< yaml file id
  integer, intent(in) :: in_id   !< block_id of model from parent

  this%yfid = in_yfid
  this%id = in_id
  this%nchildren = get_num_blocks(this%yfid, this%blockname, this%id)
  allocate(this%child_ids(this%nchildren))
end function construct_fmModel_t

function construct_fmVar_t(in_yfid, in_id) result(this)
  type (fmVar_t)      :: this    !< the var object
  integer, intent(in) :: in_yfid !< yaml file id
  integer, intent(in) :: in_id   !< block_id of var from parent

  this%yfid = in_yfid
  this%id = in_id
  this%nchildren = get_num_blocks(this%yfid, this%blockname, this%id)
  allocate(this%child_ids(this%nchildren))
end function construct_fmVar_t

function construct_fmAttr_t(in_yfid, in_id) result(this)
  type (fmAttr_t)      :: this    !< the var object
  integer, intent(in) :: in_yfid !< yaml file id
  integer, intent(in) :: in_id   !< block_id of var from parent

  this%yfid = in_yfid
  this%id = in_id
end function construct_fmAttr_t

subroutine destruct_fmTable_t(this)
  class (fmTable_t) :: this       !< the field table

  if (allocated(this%child_ids)) deallocate(this%child_ids)
  if (allocated(this%children)) then
    do table_i=1,this%nchildren
      call destruct_fmType_t(this%children(table_i))
    end do
  end if
  if (allocated(this%children)) deallocate(this%children)
end subroutine destruct_fmTable_t

subroutine destruct_fmType_t(this)
  class (fmType_t) :: this !< type object

  if (allocated(this%name)) deallocate(this%name)
  if (allocated(this%child_ids)) deallocate(this%child_ids)
  if (allocated(this%children)) then
    do type_i=1,this%nchildren
      call destruct_fmModel_t(this%children(type_i))
    end do
  end if
  if (allocated(this%children)) deallocate(this%children)
end subroutine destruct_fmType_t

subroutine destruct_fmModel_t(this)
  class (fmModel_t) :: this !< model object

  if (allocated(this%name)) deallocate(this%name)
  if (allocated(this%child_ids)) deallocate(this%child_ids)
  if (allocated(this%children)) then
    do model_i=1,this%nchildren
      call destruct_fmVar_t(this%children(model_i))
    end do
  end if
  if (allocated(this%children)) deallocate(this%children)
end subroutine destruct_fmModel_t

subroutine destruct_fmVar_t(this)
  class (fmVar_t) :: this !< variable object

  if (allocated(this%name)) deallocate(this%name)
  if (allocated(this%key_ids)) deallocate(this%key_ids)
  if (allocated(this%keys)) deallocate(this%keys)
  if (allocated(this%values)) deallocate(this%values)
  if (allocated(this%child_ids)) deallocate(this%child_ids)
  if (allocated(this%children)) then
    do var_i=1,this%nchildren
      call destruct_fmAttr_t(this%children(var_i))
    end do
  end if
  if (allocated(this%children)) deallocate(this%children)
end subroutine destruct_fmVar_t

subroutine destruct_fmAttr_t(this)
  class (fmAttr_t) :: this !< variable object

  if (allocated(this%paramname)) deallocate(this%paramname)
  if (allocated(this%keys)) deallocate(this%keys)
  if (allocated(this%values)) deallocate(this%values)
end subroutine destruct_fmAttr_t

subroutine get_blocks_fmTable_t(this)
  class (fmTable_t) :: this !< field table object

  call get_block_ids(this%yfid, this%blockname, this%child_ids)
end subroutine get_blocks_fmTable_t

subroutine get_blocks_fmType_t(this)
  class (fmType_t) :: this !< type object

  call get_block_ids(this%yfid, this%blockname, this%child_ids, this%id)
end subroutine get_blocks_fmType_t

subroutine get_name_fmType_t(this)
  class (fmType_t)     :: this !< type object
  integer              :: nkeys !< numkeys
  integer, allocatable :: key_ids(:) !< array of key ids
  character(len=256)   :: key_value  !< the value of a key

  nkeys = get_nkeys(this%yfid, this%id)
  allocate(key_ids(nkeys))
  call get_key_ids(this%yfid, this%id, key_ids)
  call get_key_value(this%yfid, key_ids(1), key_value)
  this%name = trim(key_value)
end subroutine get_name_fmType_t

subroutine get_blocks_fmModel_t(this)
  class (fmModel_t) :: this !< model object

  call get_block_ids(this%yfid, this%blockname, this%child_ids, this%id)
end subroutine get_blocks_fmModel_t

subroutine get_name_fmModel_t(this)
  class (fmModel_t)    :: this !< model object
  integer              :: nkeys !< numkeys
  integer, allocatable :: key_ids(:) !< array of key ids
  character(len=256)   :: key_value  !< the value of a key

  nkeys = get_nkeys(this%yfid, this%id)
  allocate(key_ids(nkeys))
  call get_key_ids(this%yfid, this%id, key_ids)
  call get_key_value(this%yfid, key_ids(1), key_value)
  this%name = trim(key_value)
end subroutine get_name_fmModel_t

subroutine get_blocks_fmVar_t(this)
  class (fmVar_t) :: this !< variable object

  call get_block_ids(this%yfid, this%blockname, this%child_ids, this%id)
end subroutine get_blocks_fmVar_t

subroutine get_name_fmVar_t(this)
  class (fmVar_t)      :: this       !< variable object
  integer              :: nkeys      !< numkeys
  integer              :: maxln      !< max string length names
  integer              :: maxlv      !< max string length values
  integer, allocatable :: key_ids(:) !< array of key ids
  character(len=256)   :: key_name   !< the name of a key
  character(len=256)   :: key_value  !< the value of a key

  nkeys = get_nkeys(this%yfid, this%id)
  allocate(key_ids(nkeys))
  call get_key_ids(this%yfid, this%id, key_ids)
  call get_key_value(this%yfid, key_ids(1), key_value)
  this%name = trim(key_value)
  if (nkeys .gt. 1) then
    maxln = 0
    maxlv = 0
    do var_j=2,nkeys
      call get_key_name(this%yfid, key_ids(var_j), key_name)
      call get_key_value(this%yfid, key_ids(var_j), key_value)
      maxln = max(maxln, len_trim(key_name))
      maxlv = max(maxlv, len_trim(key_value))
    end do
    allocate(this%key_ids(nkeys-1))
    allocate(character(len=maxln)::this%keys(nkeys-1))
    allocate(character(len=maxlv)::this%values(nkeys-1))
    do var_j=2,nkeys
      this%key_ids(var_j-1) = key_ids(var_j)
      call get_key_name(this%yfid, key_ids(var_j), key_name)
      call get_key_value(this%yfid, key_ids(var_j), key_value)
      this%keys(var_j-1) = trim(key_name)
      this%values(var_j-1) = trim(key_value)
    end do
  end if
end subroutine get_name_fmVar_t

subroutine get_name_fmAttr_t(this)
  class (fmAttr_t)     :: this       !< variable object
  integer              :: nkeys      !< numkeys
  integer              :: maxln      !< max string length names
  integer              :: maxlv      !< max string length values
  integer, allocatable :: key_ids(:) !< array of key ids
  character(len=256)   :: key_name   !< the name of a key
  character(len=256)   :: key_value  !< the value of a key
  character(len=256)   :: paramname  !< the value of a key

  call get_key_name(this%yfid, this%id-1, paramname)
  allocate(character(len=len_trim(paramname))::this%paramname)
  this%paramname = trim(paramname)
  nkeys = get_nkeys(this%yfid, this%id)
  allocate(key_ids(nkeys))
  call get_key_ids(this%yfid, this%id, key_ids)
  maxln = 0
  maxlv = 0
  do attr_j=1,nkeys
    call get_key_name(this%yfid, key_ids(attr_j), key_name)
    call get_key_value(this%yfid, key_ids(attr_j), key_value)
    maxln = max(maxln, len_trim(key_name))
    maxlv = max(maxlv, len_trim(key_value))
  end do
  allocate(character(len=maxln)::this%keys(nkeys))
  allocate(character(len=maxlv)::this%values(nkeys))
  do attr_j=1,nkeys
    call get_key_name(this%yfid, key_ids(attr_j), key_name)
    call get_key_value(this%yfid, key_ids(attr_j), key_value)
    this%keys(attr_j) = trim(key_name)
    this%values(attr_j) = trim(key_value)
  end do
end subroutine get_name_fmAttr_t

subroutine create_children_fmTable_t(this)
  class (fmTable_t) :: this !< the field table

  allocate(this%children(this%nchildren))
  do table_i=1,this%nchildren
    this%children(table_i) = fmType_t(this%yfid, this%child_ids(table_i))
    call this%children(table_i)%get_blocks
    call this%children(table_i)%get_name
    call this%children(table_i)%create_children
  end do
end subroutine create_children_fmTable_t

subroutine create_children_fmType_t(this)
  class (fmType_t) :: this !< type object

  allocate(this%children(this%nchildren))
  do type_i=1,this%nchildren
    this%children(type_i) = fmModel_t(this%yfid, this%child_ids(type_i))
    call this%children(type_i)%get_blocks
    call this%children(type_i)%get_name
    call this%children(type_i)%create_children
  end do
end subroutine create_children_fmType_t

subroutine create_children_fmModel_t(this)
  class (fmModel_t) :: this !< model object

  allocate(this%children(this%nchildren))
  do model_i=1,this%nchildren
    this%children(model_i) = fmVar_t(this%yfid, this%child_ids(model_i))
    call this%children(model_i)%get_blocks
    call this%children(model_i)%get_names_and_props
    call this%children(model_i)%create_children
  end do
end subroutine create_children_fmModel_t

subroutine create_children_fmVar_t(this)
  class (fmVar_t) :: this !< var object

  if (this%nchildren .gt. 0) then
    allocate(this%children(this%nchildren))
    do var_i=1,this%nchildren
      this%children(var_i) = fmAttr_t(this%yfid, this%child_ids(var_i))
      call this%children(var_i)%get_names_and_props
  end do
  end if
end subroutine create_children_fmVar_t
#endif
end module fm_yaml_mod
