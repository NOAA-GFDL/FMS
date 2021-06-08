# Documentation Style Guide

Best practices for documenting FMS code with Doxygen.

### Basics

For .F90 files:
- `!>` Starts or continues a multi-line doxygen comment
- `!!` Continues a comment
- `!<` Starts a comment(usually used after variables/parameters)
- `@commandname` to use a given doxygen command

For .c files, javadoc-style comments are used (`///` or multi-line with `/**` and `*/`), and the command specifier is instead `\commandname`.

Certain doxygen commands must be terminated with an empty comment line, mainly `@brief` and `@author`
, otherwise will include any additional documentation below it.

Some html commands, such as `<TT>` and `<br>` are allowed. Full list can be found
[here](https://www.doxygen.nl/manual/htmlcmds.html). Coding examples can be formatted by
either tabbing, or within the `@code` commands(see full subroutine documentation
below for examples of both). Links to other pages within the documentation can be created with
`@ref name`.

The following previously used commands are invalid:
- `@email`, email can be included in `@author` line
- `@description`, `@details` can be used, but is also implied by `!>` if a brief description is already given
- `@example`, not actually invalid, but should not be used in subroutine/function descriptions as it is for presenting examples of source code separately. Some examples of what can be used for programming examples in descriptions are below.
### Documenting subroutines and functions:

The first `!>` comment above a routine implicitly starts the brief description, and the second `!>`
will start the full description (as long as the brief description is terminated).

Simple subroutine/function documentation:

<!--- @code{.F90} -->
```
!> description
!! continued description
function foo()
  integer :: bar !< variable description
```
<!--- @endcode -->

Full subroutine/function documentation:
<!--- @code{.F90} -->
```
!> @brief short description
!! continued short description
!!
!> Long description
!! continued long description
!! <br>Example 1:
!!
!! 		foo()
!!
!! Example 2:
!! @code{.F90}
!! foo()
!! @endcode
subroutine foo()
  integer :: bar !< variable description
  ...
```
<!--- @endcode -->

### Documenting Interfaces
Doxygen supports documenting interfaces, just not Fortran's version of an interface. As a workaround,
the `@page interface_name Interface name to display` command can be used to create sections for each
interface. The descriptions work the same as routines/functions, but any parameters/return values
must be defined with commands as the interface definition is not actually getting parsed.

Example:
<!--- @code{.F90} -->
```
!> @page foo foo Interface
!> @brief short description
!!
!> Longer description
!! @param[inout] arg1 description
!! @returns retval description
```
<!--- @endcode -->

### Documenting Type Definitions
Type definitions are only parsed through doxygen if they are defined in a certain way. They must
include the visibility(public or private) and `::` otherwise will likely be missed by the parser.
Other than that they are documented the same as subroutines.

Example:
<!--- @code{.F90} -->
```
!> description
type, public :: typename
    integer :: num !< description
end type typename
```
<!--- @endcode -->

### Module Grouping
To get the modules parsed correctly, each module page is created through a grouping defined at the
top of the file. Each FMS module is then also in a group corresponding to its subdirectory, with
these groups defined in `docs/grouping.h`. The groups can also contain any other relevant
files, such as additional documentation files, header, or include files.  Files or additional
documentation can be added to either directory or module groups by adding `@ingroup group_name`
to a documentation block.
