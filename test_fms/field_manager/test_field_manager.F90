program test_field_manager

use field_manager_mod
use mpp_mod, only : mpp_exit, mpp_pe, mpp_root_pe, mpp_error, NOTE

implicit none
!#include "mpif.h"


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
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! By adding the optional .true. argument you get a recursive listing of the fields.
write(*,*) "Here's a recursive listing"
success = fm_dump_list("/", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Using fm_dump_list with a blank first argument returns the last field accessed by field manager.
write(*,*) 'Dumping last field changed to by field_manager using fm_change_list'
success = fm_dump_list("", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Change list to look at the land model fields
write(*,*) 'Changing list to land_mod'
success = fm_change_list("/land_mod")
write(*,*) 'Dumping last list changed to by field_manager using fm_change_list i.e list of land model fields'
success = fm_dump_list("", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Now let's modify some of the field entries.
! 
!In this example we add a field ( convection = 'off' ) to the radon list
write(*,*) "ADDING convection = off TO RADON LIST"
!if ( fm_change_list('/atmos_mod/tracer/radon')) then
if ( fm_exists('/atmos_mod/tracer/radon')) then
   write(*,*) "'/atmos_mod/tracer/radon' exists "
   success = fm_change_list('/atmos_mod/tracer/radon')
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
if (.not. success ) call mpp_error(NOTE, "There is no atmos model radon field defined in the field_table")
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
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


index = fm_get_length("/ocean_mod/tracer") 
write(*,*) "The length of the current list '/ocean_mod/tracer' is ",index," i.e."
success = fm_dump_list("/ocean_mod/tracer")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

! Find out what type of field this is. Possibilities are real, integer, string, logical, and list
name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'The type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type

success = fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',str)
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (character) ',str


write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) "MODIFYING BIOTIC1 FIELD slope ATTRIBUTE TO slope = 0.95 "
if ( fm_change_list('/ocean_mod/tracer/biotic1/diff_horiz/linear')) &
   index = fm_new_value('slope',0.95, index = 1)

! Dump the listing of the modified ocean model tracer attribute
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

name_field_type = fm_get_type('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope')
write(*,*) 'Now the type for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is ',name_field_type
success =  fm_get_value('/ocean_mod/tracer/biotic1/diff_horiz/linear/slope',param)
write(*,*) 'The value for /ocean_mod/tracer/biotic1/diff_horiz/linear/slope is (real) ',param
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

write(*,*) 'Changing the name of biotic1 to biotic_control'
success = fm_modify_name('/ocean_mod/tracer/biotic1', 'biotic_control')

! Dump the listing of the modified tracer
success = fm_dump_list("/ocean_mod/tracer/biotic_control", .true.)

! Double check to show that the tracer has been renamed and the original doesn't exist anymore. 
success = fm_dump_list("/ocean_mod/tracer/biotic1", .true.)
if (.not. success ) call mpp_error(NOTE, "Ocean model tracer biotic1 does not exist anymore.")
write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'


if ( fm_change_list("/ocean_mod/tracer/age_ctl") ) then
success = fm_dump_list("", .true.)
write(*,*) "Now we'll add a new list to this list"
index = fm_new_list("units",create = .true.)

success = fm_dump_list("", .true.)

write(*,*) "Now we'll give it a value"
if (success) index = fm_new_value('units','days')

success = fm_dump_list("", .true.)


write(*,*) '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
endif
!errorcode = 121
!CALL MPI_ERROR_STRING(errorcode, string, resultlen, ierror)
!write(*,*) string
call field_manager_end

call mpp_exit

end program test_field_manager
