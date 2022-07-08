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
!> @defgroup fms_yaml_output_mod fms_yaml_output_mod
!> @ingroup parser
!> @author Tom Robinson
!> @description Writes a 3 tiered yaml where the first and second tier can have 1 key that has an 
!! array of values.  This is usefule for writing a diag_output.yaml Here is an example:
!! \verbatim
!! ---
!! basedate: 0 0 0 0 0 0 0
!! experiment: c192L65_exp
!! files:
!! - name: atmos_daily
!!   freq: days
!!   vars:
!!   - name: temp
!!     units: K
!!   - name: P
!!     longname: pressure
!!   - name: wind
!!     units: m/s
!!     longname: wind velocity magnitude
!! - name: atmos_monthly
!!   freq: months
!!   vars:
!!   - name: temp
!!     units: K
!! ...
!! \endverbatim
!! In this example, basedate, experiment, and files are the level 1 (top level) keys.  files is the 
!! level2key in the keys struct.  The array of structs for the files with name, freq has a size of 2, and
!! vars is the level2key.  The key3 and var3 arrays should have a size of 4, and a3each=(\3,1\)
!! corresponding to the number of elements in each array within the yaml.

!> @file
!> @brief File for @ref fms_yaml_output_mod

!> @addtogroup fms_yaml_output_mod
!> @{
module fms_yaml_output_mod
#ifdef use_yaml

use iso_c_binding
implicit none

integer, parameter :: string_len_parameter = 255 !< Number of characters in the keys and strings.
                                                 !! Must match whats in yaml_output_functions.c 
!> Keys for the output yaml on a given level corresponding to the struct in yaml_output_functions.c
type, bind(c) :: fmsYamlOutKeys_type
        character (c_char) :: key1 (string_len_parameter)
        character (c_char) :: key2 (string_len_parameter)
        character (c_char) :: key3 (string_len_parameter)
        character (c_char) :: key4 (string_len_parameter)
        character (c_char) :: key5 (string_len_parameter)
        character (c_char) :: key6 (string_len_parameter)
        character (c_char) :: key7 (string_len_parameter)
        character (c_char) :: key8 (string_len_parameter)
        character (c_char) :: key9 (string_len_parameter)
        character (c_char) :: key10 (string_len_parameter)
        character (c_char) :: key11 (string_len_parameter)
        character (c_char) :: key12 (string_len_parameter)
        character (c_char) :: key13 (string_len_parameter)
        character (c_char) :: key14 (string_len_parameter)
        character (c_char) :: key15 (string_len_parameter)
        character (c_char) :: level2key (string_len_parameter)
end type fmsYamlOutKeys_type
!> Values for the output yaml on a given level corresponding to the struct in yaml_output_functions.c
type, bind(c) :: fmsYamlOutValues_type
        character (c_char) :: val1 (string_len_parameter)
        character (c_char) :: val2 (string_len_parameter)
        character (c_char) :: val3 (string_len_parameter)
        character (c_char) :: val4 (string_len_parameter)
        character (c_char) :: val5 (string_len_parameter)
        character (c_char) :: val6 (string_len_parameter)
        character (c_char) :: val7 (string_len_parameter)
        character (c_char) :: val8 (string_len_parameter)
        character (c_char) :: val9 (string_len_parameter)
        character (c_char) :: val10 (string_len_parameter)
        character (c_char) :: val11 (string_len_parameter)
        character (c_char) :: val12 (string_len_parameter)
        character (c_char) :: val13 (string_len_parameter)
        character (c_char) :: val14 (string_len_parameter)
        character (c_char) :: val15 (string_len_parameter)
end type fmsYamlOutValues_type


interface
subroutine write_yaml_from_struct_3 (yamlname, a1size, keys, vals, a2size, key2, val2, a3size, a3each,&
                       key3, val3) bind(C, name="write_yaml_from_struct_3")
use iso_c_binding
import fmsYamlOutKeys_type, fmsYamlOutValues_type
character (c_char) :: yamlname !< The output yaml file name
integer (c_int), value :: a1size !< The size of the first yaml array (only supports 1)
type (fmsYamlOutKeys_type) :: keys(a1size) !< Top level yaml keys
type (fmsYamlOutValues_type) :: vals(a1size) !< Values corresponding to keys
integer (c_int), value :: a2size !< The size of the second yaml array
type (fmsYamlOutKeys_type) :: key2(a2size) !< Second level keys
type (fmsYamlOutValues_type) :: val2(a2size) !< Values corresponding to key2
integer (c_int), value :: a3size !< The size of the third yaml array
integer (c_int) :: a3each (a2size) !< Array that has the number of elements for each level 2 array's 
                                   !! third level elements
type (fmsYamlOutKeys_type) :: key3(a3size) !< Third level keys
type (fmsYamlOutValues_type) :: val3(a3size) !< Values corresponding to keys2

end subroutine write_yaml_from_struct_3

end interface

#endif
end module fms_yaml_output_mod
!> @}
! close documentation grouping
