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

!********************* Sample field table required: *********************
! "TRACER", "ocean_mod", "biotic1"
!           "diff_horiz", "linear", "slope=ok"
!           "longname", "biotic one" /
! "TRACER", "ocean_mod", "age_ctl" /
! "TRACER", "atmos_mod","radon"
!           "longname","radon-222"
!           "units","VMR*1E21"
!           "profile_type","fixed","surface_value=0.0E+00"
!           "convection","all"/
! "TRACER", "land_mod", "sphum"
!           "longname",     "specific humidity"
!            "units",        "kg/kg" /
!***********************************************************************

program test_field_manager

use field_manager_mod
use mpp_mod, only : mpp_exit, mpp_pe, mpp_root_pe, mpp_error, NOTE, FATAL

implicit none


integer :: i, j, nfields, num_methods, model
character(len=fm_string_len) :: field_type, field_name, str, name_field_type, path
character(len=512) :: method_name, method_control
real :: param
integer :: flag, index
logical :: success
type(method_type), dimension(20) :: methods

call field_manager_init(nfields)

! Dump the list of fields produced from reading the field_table

! Here are the lists that propagate off the root "/"
! By calling fm_dump_list with a single argument you only get the
! lists branching off this argument in the list.
write(*,*) "Here's a baseline listing"
success = fm_dump_list("/")
if (.not. success ) call mpp_error(FATAL, "ERROR: unable to get baseline listing")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! By adding the optional .true. argument you get a recursive listing of the fields.
write(*,*) "Here's a recursive listing"
success = fm_dump_list("/", .true.)
if (.not. success ) call mpp_error(FATAL, "ERROR: unable to get recursive listing")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Using fm_dump_list with a blank first argument returns the last field accessed by field manager.
write(*,*) 'Dumping last field changed to by field_manager using fm_change_list'
success = fm_dump_list("", .true.)
if (.not. success ) call mpp_error(FATAL, "ERROR: unable to get the last field acessed")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Change list to look at the land model fields
write(*,*) 'Changing list to land_mod'
success = fm_change_list("/land_mod")
if (.not. success ) call mpp_error(FATAL, "ERROR: unable to look at the land model fields")
write(*,*) 'Dumping last list changed to by field_manager using fm_change_list i.e list of land model fields'
success = fm_dump_list("", .true.)
if (.not. success ) call mpp_error(FATAL, "ERROR: unable to look at the last field changed in land model fields")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Now let's modify some of the field entries.
!
!In this example we add a field ( convection = 'off' ) to the radon list
write(*,*) "ADDING convection = off TO RADON LIST"
!if ( fm_change_list('/atmos_mod/tracer/radon')) then
if ( fm_exists('/atmos_mod/tracer/radon')) then
   write(*,*) "'/atmos_mod/tracer/radon' exists "
   success = fm_change_list('/atmos_mod/tracer/radon')
   if (.not. success ) call mpp_error(FATAL, "ERROR: unable to change to the radon list")
! The next line creates a new field branching off radon.
   index = fm_new_value('convection','off')
endif

success = fm_query_method('radon',method_name,method_control)
if (success ) then
call mpp_error(NOTE, "Method names for radon is/are "//trim(method_name))
call mpp_error(NOTE, "Method controls for radon is/are "//trim(method_control))
else
call mpp_error(NOTE, "There is no atmos model radon field defined in the field_table")
endif
! Dump the listing of the modified tracer
success = fm_dump_list("/atmos_mod/tracer/radon", .true.)
if (.not. success ) call mpp_error(FATAL, "There is no atmos model radon field defined in the field_table")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


! Find out what the current path is. Should be '/atmos_mod/tracer/radon' as set in fm_change_list above.
path = fm_get_current_list()
write(*,*) 'Current path is ',trim(path)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Now let's modify the value of the field we just added.
write(*,*) "MODIFYING RADON FIELD CONVECTION ATTRIBUTE TO convection = RAS_off "
index = fm_new_value('convection','RAS_off')

! Dump the listing of the modified tracer
success = fm_dump_list("/atmos_mod/tracer/radon", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) "ORIGINAL OCEAN MODEL TRACER FIELDS"

! Dump the listing of the original ocean model tracers
success = fm_dump_list("/ocean_mod/tracer", .true.)
if (.not. success ) call mpp_error(FATAL, "Unable to dump the listing of the original ocean model tracers")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


index = fm_get_length("/ocean_mod/tracer")
write(*,*) "The length of the current list '/ocean_mod/tracer' is ",index," i.e."
success = fm_dump_list("/ocean_mod/tracer")
if (.not. success ) call mpp_error(FATAL, "Unable to get the length of the current list")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Find out what type of field this is. Possibilities are real, integer, string, logical, and list
name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'The type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type

success = fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',str)
if (.not. success ) call mpp_error(FATAL, "Unable to get the value of this field")
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (character) ',str


write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) "MODIFYING BIOTIC1 FIELD slope ATTRIBUTE TO slope = 0.95 "
if ( fm_change_list('/ocean_mod/tracer/biotic1/diff_horiz/linear')) &
   index = fm_new_value('slope',0.95, index = 1)

! Dump the listing of the modified ocean model tracer attribute
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
if (.not. success) call mpp_error(FATAL, "Unable to dump the listing of the modified ocean model tracer attribute")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'Now the type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type
success =  fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',param)
if (.not. success) call mpp_error(FATAL, "Unable to get the value of biotic1 slope")
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (real) ',param
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) 'Changing the name of biotic1 to biotic_control'
success = fm_modify_name('/ocean_mod/tracer/biotic1', 'biotic_control')
if (.not. success) call mpp_error(FATAL, "Unable to change the name of biotic1 to biotic_control")

! Dump the listing of the modified tracer
success = fm_dump_list("/ocean_mod/tracer/biotic_control", .true.)
if (.not. success) call mpp_error(FATAL, "Unable to dump the listing of the modified tracers")

! Double check to show that the tracer has been renamed and the original doesn't exist anymore.
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
if (success ) call mpp_error(NOTE, "Ocean model tracer biotic1 still exists.")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


if ( fm_change_list("/ocean_mod/tracer/age_ctl") ) then
success = fm_dump_list("", .true.)
write(*,*) "Now we'll add a new list to this list"
index = fm_new_list("units",create = .true.)

success = fm_dump_list("", .true.)
if (.not. success) call mpp_error(FATAL, "Unable to dump the listing to the recently changed field")

write(*,*) "Now we'll give it a value"
if (success) then
    index = fm_new_value('units','days')
else
    call mpp_error(FATAL, "Unable to give field new value")
endif

success = fm_dump_list("", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

else
call mpp_error(FATAL, "Unable to modify BIOTIC1 field")

endif
!errorcode = 121
!CALL MPI_ERROR_STRING(errorcode, string, resultlen, ierror)
!write(*,*) string
call field_manager_end

call mpp_exit

end program test_field_manager
