# FMS Doxygen Documentation Guide

### Basics
For .F90 files:
- `!>` Starts or continues a multi-line doxygen comment
- `!!` Continues a comment
- `!<` Starts a comment
- `@commandname` to use a given doxygen command

For .c files, multi-line javadoc-style blocks are used, and the command specifier is instead `\commandname`.

Certain doxygen commands must be terminated with an empty comment line, such as `@brief` or `@author`.

Simple html commands, namely `<TT>` and `<br>` are allowed. Coding examples can be formatted by
either tabbing as seen above, or within the `@code` commands(see full subroutine documentation 
below for an example).

The following previously used commands are invalid:
- `@email`, email can be included in @author line
- `@description`, @details can be used, but is also implied by `!>` if a brief description is already given
- `@example`, not actually invalid, but should not be used in subroutine/function descriptions as it is for presenting examples of source code. @code{.F90} can be used to enclose programming examples for routines

### Documenting subroutines and functions:

The first `!>` comment above a routine implicitly starts the brief description, and the second `!>`
will start the full description (as long as the brief description is terminated).

Simple subroutine/function documentation:

```
!> brief description
!! continued description
function foo()
  integer :: bar !< variable description
```

Full subroutine/function documentation:
```
!> @brief short description
!! continued short description
!!
!> @author author_name @link author@email.com @endlink
!!
!! Long description
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
### Documenting Interfaces
Fortran supports documenting interfaces, just not fortran interfaces. As a workaround, the 
`@page interface_name "Interface name to display"` command is used to create sections for each 
interface. The descriptions work the same as routines/functions, but any parameters/return values
must be defined with commands as nothing is actually getting parsed.

Example:
```
!> @page foo "foo Interface"
!> @brief short description
!!
!> Longer description
!! @param [inout] arg1 description
!! @returns some_return_value
``` 

### Module Grouping
Each FMS module is in a group corresponding to subdirectory, with these groups are defined in
`docs/grouping.h`. The groups can also contain any other relevant files, such as additional
documentation files, header, or include files. Files be added by adding `@ingroup dir_name` to 
a documentation block.
