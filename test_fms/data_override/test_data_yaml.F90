program test_data_yaml

use data_override_mod
use fms_mod, only: fms_init, fms_end

implicit none

#ifndef use_yaml
type(data_type), dimension(100) :: data_table
#else
type(data_type), dimension(:), allocatable :: data_table2
#endif

integer :: ntable
integer :: i

call fms_init()
#ifndef use_yaml
call read_table(data_table, ntable)
do i = 1, ntable
   print *, "Entry number:", i
   print *, "gridname:", trim(data_table(i)%gridname)
   print *, "fieldname_code:", trim(data_table(i)%fieldname_code)
   print *, "fieldname_file:", trim(data_table(i)%fieldname_file)
   print *, "file_name:", trim(data_table(i)%file_name)
   print *, "interpol_method:", trim(data_table(i)%interpol_method)
   print *, "factor:", data_table(i)%factor
   print *, "grid:", data_table(i)%lon_start, data_table(i)%lon_end, data_table(i)%lat_start, data_table(i)%lat_end
   print *, "region_type:", data_table(i)%region_type
   print *, ""
enddo
#else
call read_table_yaml(data_table2, ntable)
do i = 1, ntable
   print *, "Entry number:", i
   print *, "gridname:", trim(data_table2(i)%gridname)
   print *, "fieldname_code:", trim(data_table2(i)%fieldname_code)
   print *, "fieldname_file:", trim(data_table2(i)%fieldname_file)
   print *, "file_name:", trim(data_table2(i)%file_name)
   print *, "interpol_method:", trim(data_table2(i)%interpol_method)
   print *, "factor:", data_table2(i)%factor
   print *, "grid:", data_table2(i)%lon_start, data_table2(i)%lon_end, data_table2(i)%lat_start, data_table2(i)%lat_end
   print *, "region_type:", data_table2(i)%region_type
   print *, ""
enddo
#endif

call fms_end()

end program test_data_yaml
