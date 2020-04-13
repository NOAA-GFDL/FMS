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
* Use `KIND` parameters from intrinsic fortran modules such as iso_fortran_env
  or iso_c_binding to ensure portability
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
* Functions should include a result variable on its own line, that does not have
  a specific intent.
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

## Fortran Example

```Fortran

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
  use, intrinsic :: iso_fortran_env, only: INT32, REAL32
  use util_mod, only: util_func1
  implicit none
  private

  public :: sub1
  public :: func1

  !> @brief Doxygen description of type.
  type,public :: CustomType
    private
    integer(kind=INT32) :: a_var !< Inline doxygen description.
    real(kind=REAL32),dimension(:),allocatable :: b_arr !< long description
                                                        !! continued on
                                                        !! multiple lines.
  endtype CustomType

  contains

  !> @brief Doxygen description.
  subroutine sub1(arg1, &
    & arg2, &
    & arg3)
    real(kind=REAL32),intent(in) :: arg1 !< Inline doxygen description.
    integer(kind=INT32),intent(inout) :: arg2 !< Inline doxygen description.
    character(len=*),intent(inout) :: arg3 !< Long inline doxygen
                                           !! description.
  end subroutine sub1

  !> @brief Doxygen description
  !! @return Function return value.
  function func1(arg1, &
    & arg2) &
    & result(res)
    integer(kind=INT32),intent(in) :: arg1 !< Inline doxygen description
    integer(kind=INT32),intent(in) :: arg2 !< Inline doxygen description
    integer(kind=INT32) :: res
  end function func1

end module example_mod
```

## C/C++

### General
* C code is written in GNU style.  Each new level in a program block is indented
  by 2 spaces. Braces start on a new line, and are also indented by 2 spaces.
* See the [Gnome C coding style guide](https://developer.gnome.org/programming-guidelines/stable/c-coding-style.html.en)
  for more information
