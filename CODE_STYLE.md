# Coding Style

## General

* Trim all trailing whitespace from every line (some editors can do this
  automatically).
* No <Tab> characters.
* Supply a header for each file with a description of the file and the author(s)
  name or GitHub ID.
* A copy of the [Gnu Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0.en.html)
  must be included at the top of each file.
* Documentation should be written so that it can be parsed by [Doxygen](http://www.doxygen.nl/).
* All variables should be defined, and include units. Unit-less variables should be marked `unitless`
* Provide detailed descriptions of modules, interfaces, functions, and subroutines
* Define all function/subroutine arguments, and function results (see below)
* Follow coding style of the current file, as much as possible.

## Fortran

### General

* Use Fortran 95 standard or newer
* Two space indentation
* Use `KIND` parameters from platform_mod
* Never use implicit variables (i.e., always specify `IMPLICIT NONE`)
* Lines must be <= 120 characters long (including comments)
* logical, compound logical, and relational if statements may be one line,
  using “&” for line continuation if necessary:
  ```Fortran
  if(file_exists(fileName)) call open_file(fileObj,fileName, is_restart=.false)
  ```
* Avoid the use of `GOTO` statements
* Avoid the use of Fortran reserved words as variables (e.g. `DATA`, `NAME`)
* Avoid the use of `COMMON` blocks

### Derived types

* Type names must be in CapitalWord format.
* Variables names must be in underscore_word format.
* All member variables must be private.
* Doxygen description on the line before the type definition.
* Inline doxygen descriptions for all member variables.

## Functions
* If a function has a result variable, it should be declared on its own line,
  and the variable should not be declared with a specific intent.
* Inline doxygen descriptions for all arguments, except the result variable.
* Doxygen description on the line(s) before the function definition.  This must
  specify what the function is returning using the `@return` doxygen keyword.

## Blocks
* terminate `do` loops with `enddo`
* terminate block `if`, `then` statements with `endif`

## OpenMP

* Directives should start at the beginning of the line, and be in lowercase.
* All openMP directives should specify default(none), and then explicitly list
  all shared and private variables.
* All critical sections must have a unique name.

## Precision
* Precision of all real arguments should be explicitly defined as r4_kind, r8_kind,
  or as any other precision parameters defined in platform_mod.
* The precision of real numerical values should be consistent with the precision
  of the associated variable.  For example, if the variable `a` has been declared
  as r4_kind, then `a=1.4_r4_kind`, not a=1.4.

## Macros
* All letters in the macro name are capitalized
* All macro names end with an underscore "_"
* All precision related macro names start with the letters "FMS"
* Macro names should be unique to each module.  For example,
  `FMS_AU_KIND_` is used in the axis_utils_mod.  `FMS_HI_KIND_` is used
  in the horiz_interp_mod

## .fh files
* The .fh header files contain macro definitions.
* .fh files containing precision related macro definitions should be named
  with `_r4.fh` and `_r8.fh` extensions in the include subdirectory found
  in the module directory.

## Fortran Example

```Fortran example.F90 file

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

!> @file
!! @brief Example code
!! @author <developer>
!! @email gfdl.climate.model.info@noaa.gov

module example_mod
  use platform_mod, only r4_kind, r8_kind, i4_kind, i8_kind
  use util_mod, only: util_func1
  implicit none
  private

  public :: sub1
  public :: func1

  !> @brief Doxygen description of type.
  type,public :: CustomType
    private
    integer(kind=i4_kind) :: a_var !< Inline doxygen description.
    real(kind=r8_kind),dimension(:),allocatable :: b_arr !< long description
                                                         !! continued on
                                                         !! multiple lines.
  endtype CustomType

  contains

  !> @brief Doxygen description.
  subroutine sub1(arg1, arg2, &
    & arg3)
    real(kind=r4_kind),intent(in) :: arg1 !< Inline doxygen description.
    integer(kind=i8_kind),intent(inout) :: arg2 !< Inline doxygen description.
    character(len=*),intent(inout) :: arg3 !< Long inline doxygen
                                           !! description.

    arg1=2.456_r4_kind
  end subroutine sub1

  !> @brief Doxygen description
  !! @return Function return value.
  function func1(arg1, arg2) result(res)
    integer(kind=i4_kind),intent(in) :: arg1 !< Inline doxygen description
    integer(kind=i4_kind),intent(in) :: arg2 !< Inline doxygen description
    integer(kind=r8_kind) :: res

    res=real(arg1,r8_kind) * 3.14_r8_kind
  end function func1

#include "example_r4.fh"
#include "example_r8.fh"

end module example_mod
```
```Fortran example_r4.fh file
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

!> @file
!! @brief Example _r4.fh file containing macro definitions
!! @author <developer>
!! @email gfdl.climate.model.info@noaa.gov

#undef   FMS_EX_KIND_
#define  FMS_EX_KIND_ r4_kind

#undef  FMS_EX_SUBROUTINE_
#define FMS_EX_SUBROUTINE_ fms_ex_subroutine_r4

#include "example.inc"
```
```Fortran example_r8.fh file
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

!> @file
!! @brief Example file _r8.fh file containing macro definitions
!! @author <developer>
!! @email gfdl.climate.model.info@noaa.gov

#undef   FMS_EX_KIND_
#define  FMS_EX_KIND_ r8_kind

#undef  FMS_EX_SUBROUTINE_
#define FMS_EX_SUBROUTINE_ fms_ex_subroutine_r8

#include "example.inc"
```
``` Fortran example.inc file
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
!> @file
!! @brief Example .inc file containing subroutine definitions/declarations
!! @author <developer>
!! @email gfdl.climate.model.info@noaa.gov

subroutine FMS_EX_SUBROUTINE_(arg1, arg2, arg3)
  real(FMS_EX_KIND_) :: arg1, arg2
  integer(i4_kind) :: arg3
  integer, parameter :: lkind=FMS_EX_KIND_

  arg1 = arg2 / 4.0_lkind

end subroutine FMS_EX_SUBROUTINE_
```

## C/C++

### General
* C code is written in GNU style.  Each new level in a program block is indented
  by 2 spaces. Braces start on a new line, and are also indented by 2 spaces.
* See the [Gnome C coding style guide](https://developer.gnome.org/programming-guidelines/stable/c-coding-style.html.en)
  for more information
