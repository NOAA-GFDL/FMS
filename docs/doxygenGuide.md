# Documentation Style Guide

Best practices for documenting FMS code with Doxygen.

### Basics

For .F90 files:
- `!>` Starts or continues a multi-line doxygen comment
- `!!` Continues a comment
- `!<` Starts a same-line comment (used for documenting variable/parameter declarations)
- `@commandname` to use a given doxygen command

For .c files, javadoc-style comments are used (`///` or multi-line with `/**` and `*/`), and the command specifier is instead `\commandname`.

Certain doxygen commands must be terminated with an empty comment line, mainly `@brief` and `@author`
, otherwise will include any additional documentation below it.

The following previously used commands are invalid:
- `@email`, email can be included in `@author` line
- `@description`, `@details` can be used, but is also implied by `!>` if a brief description is already given
- `@example`, not actually invalid, but should not be used in subroutine/function descriptions as it is for presenting examples of source code separately. Some examples of what can be used for programming examples in descriptions are below.

Some html commands, such as `<TT>` and `<br>` are allowed. Full list can be found
[**here**](https://www.doxygen.nl/manual/htmlcmds.html). Coding examples can be formatted by
either tabbing, or within the `@code` commands(see full subroutine documentation
below for examples of both). Links to other pages within the documentation can be created with
`@ref name`.

### Documenting Subroutines and Functions

The first `!>` comment above a routine implicitly starts the brief description, and the second `!>`
will start the full description (as long as the brief description is terminated).

Simple subroutine/function documentation:

@code{.F90}
!> description
!! continued description
function foo()
  integer :: bar !< variable description
  ...
@endcode

Full subroutine/function documentation:

@code{.F90}
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
!! endcode (`@` ommitted to not close this example)
subroutine foo()
  integer :: bar !< variable description
  ...
@endcode

### Module Grouping
To get the modules parsed correctly, each module page is created through a grouping defined at the
top of the file. Each FMS module is then also in a group corresponding to its subdirectory, with
these groups defined in `docs/grouping.h`. The groups can also contain any other relevant
files, such as additional documentation files, header, or include files.  Files or additional
documentation can be added to either directory or module groups by adding `@ingroup group_name`
to a specific objects documentation block, or to enclose any section of code with the following:

@code{.F90}
!> @addtogroup foo_mod
!> @{
< any documented code within here will be included in foo_mod >
!> @}
@endcode

### Documenting Interfaces and Type Definitions

Interfaces and type definitions can be documentated the in a similar way as subroutines/functions.
In the documentation they will both be listed data types and will show any included subroutines
or variables with links if available.

For interfaces and type definitions, duplicates can appear for subroutines/functions and variables
contained in interfaces and typedefs if they are included in the module group as well.
This can be avoided by ending the module group before the
interface is defined with `!> @}` and adding `@ingroup module_name_mod` to the interface
documentation blocks. This will add the interface/type to the module page without any redundant
items. After any interface or type definitions are defined, `!> @addtogroup mod_name_mod` and a
following `!> @{` must be added after or else the rest of the module will be excluded.

Example:
@code{.F90}
!> @}
! closes module group started at beginning of file

!> @brief short description
!!
!> Longer description
!! @param[inout] arg1 description
!! @returns retval description
!> @ingroup foo_mod
interface foo
  ...
end interface foo

!> description
!> @ingroup foo_mod
type, public :: typename
  integer :: num !< description
end type typename

< any other interface/type definitions and documentation >

! add rest of module
!> @addtogroup foo_mod
!> @{
...
@endcode
