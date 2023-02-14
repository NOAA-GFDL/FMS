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

!> @defgroup fms_diag_time_reduction_mod fms_diag_time_reduction_mod
!> @ingroup diag_manager
!> @brief fms_diag_time_reduction_mod defines classes encapsulating the diag_manager
!! time redution types.
!!
!> @author Miguel Zuniga
!!
!> @file
!> @brief File for @ref fms_diag_time_reduction_mod
!> @addtogroup fms_diag_time_reduction_mod
!> @{
MODULE fms_diag_time_reduction_mod

  USE diag_data_mod, only: EVERY_TIME
  USE fms_mod, ONLY: error_mesg, FATAL

  implicit none

  !!TODO: (Future effort) Note that time_diurnal processing is a little different
  !! and more complex than the other reduction methods, and therefore refactoring its
  !! processing may simplify the overall related codebase. The refactoring,
  !! if possible, may be done elsewhere in the diag_manager.

  !!These parametes are the possible kinds of time reduction operations.
  INTEGER, PARAMETER :: time_none    = 0 !< There is no reduction method
  INTEGER, PARAMETER :: time_average = 1 !< The reduction method is average
  INTEGER, PARAMETER :: time_rms     = 2 !< The reduction method is rms
  INTEGER, PARAMETER :: time_max     = 3 !< The reduction method is max
  INTEGER, PARAMETER :: time_min     = 4 !< The reduction method is min
  INTEGER, PARAMETER :: time_sum     = 5 !< The reudction method is sum
  INTEGER, PARAMETER :: time_diurnal = 6 !< The reduction method is diurnal
  INTEGER, PARAMETER :: time_power   = 7 !< The reduction method is power

  !> @brief Class fmsDiagTimeReduction_type has an encapsulation of the "Fortran enum" time
  !! reduction integer parameters, plus an encapsulation of the groupings of
  !! the time reduction types. It is intended to provide some of the functionality
  !! that was coded in the legacy function diag_data.F90:init_output_fields.
  !! The functionality in the end is used by send_data in (EFFICIENT) do loops calling
  !! the weighting or math functions to update buffers.
  !! The integer parameters above are the legal time reduction types,
  !! but they are not necessarily mutually exclusive in some contexts.
  !!
  !> @addtogroup fms_diag_time_reduction_mod
  TYPE fmsDiagTimeReduction_type
     integer , private :: the_time_reduction !< The time reduction type, as an  integer defined above.
     logical , private :: time_averaging !< Set true iff time_average, time_rms, time_power or time_diurnal is true
     logical , private :: time_ops !< Set true iff time_min, time_max, time_rms or time_average is true.
   CONTAINS
     procedure, public :: do_time_averaging => do_time_averaging_imp
     procedure, public :: has_time_ops => has_time_ops_imp
     procedure, public :: is_time_none => is_time_none_imp
     procedure, public :: is_time_average => is_time_average_imp
     procedure, public :: is_time_rms => is_time_rms_imp
     procedure, public :: is_time_max => is_time_max_imp
     procedure, public :: is_time_min => is_time_min_imp
     procedure, public :: is_time_sum => is_time_sum_imp
     procedure, public :: is_time_diurnal => is_time_diurnal_imp
     procedure, public :: is_time_power => is_time_power_imp
     procedure, public :: initialize
     procedure, public :: copy
  END TYPE fmsDiagTimeReduction_type

  !> @brief This interface is for the class constructor.
  !> @addtogroup fms_diag_time_reduction_mod
  interface  fmsDiagTimeReduction_type
     procedure  :: fmsDiagTimeReduction_type_constructor
  end interface fmsDiagTimeReduction_type

CONTAINS

  !> @brief The class contructors. Just allocates the class and calls an initializer
  !! @return An allocated instance of fmsDiagTimeReduction_type, which is nitialized using
  !! provided values for arguments dt and out_freqeuncy.
  function fmsDiagTimeReduction_type_constructor(dt, out_frequency) result(time_redux)
    integer, intent(in) :: dt  !< The redution type (time_rms, time_power, etc)
    integer, intent(in) :: out_frequency  !< The output frequency.
    class (fmsDiagTimeReduction_type), allocatable :: time_redux !< The instance of the fmsDiagTimeReduction_type
    !!class allocated and returned by this constructor.
    allocate(time_redux)
    call time_redux%initialize(dt, out_frequency)
  end function fmsDiagTimeReduction_type_constructor

  !> @brief Initialize the object. As an alternative to the constructor, one can
  !! allocate an fmsDiagTimeReduction_type instance, then call its initialize function.
  subroutine initialize(this, dt, out_frequency)
    class (fmsDiagTimeReduction_type), intent(inout) :: this !< The fmsDiagTimeReduction_type object
    integer, intent(in) :: dt !< The redution type (time_rms, time_porer, etc)
    integer, intent(in) :: out_frequency !< The output frequency.

    this%the_time_reduction = dt

    !! Set the time_averaging flag
    !! See legacy init_ouput_fields function, lines 1470ff
    IF(( dt .EQ. time_average) .OR. (dt .EQ. time_rms) .OR. (dt .EQ. time_power) .OR. &
         & (dt .EQ. time_diurnal)) THEN
       this%time_averaging = .true.
    ELSE
       this%time_averaging= .false.
       IF((dt .NE. time_max) .AND. (dt .ne. time_min) .AND. (dt .NE. time_sum) &
            & .AND. (dt .NE. time_none)) THEN
          CALL error_mesg('fmsDiagTimeReduction_type: initialize', &
               & 'time_averaging=.false. but reduction type not compatible', FATAL)
       ENDIF
    END IF

    !!TODO: (MDM) Add other checks? E.g. If time_averaging == .false., then
    !!  out_frequency == EVERY_TIME

    IF((dt .EQ. time_min) .OR. (dt .EQ. time_max) .OR. &
         & ( dt .EQ. time_average) .OR. (dt .EQ. time_sum)  ) THEN
       this%time_ops = .true.
    ELSE
       this%time_ops = .false.
    END IF
  end subroutine initialize

  !> @brief Copy the source time reduction object into the this object.
  subroutine copy(this, source)
    class (fmsDiagTimeReduction_type),intent(inout):: this !< The fmsDiagTimeReduction_type object
    class (fmsDiagTimeReduction_type),intent(in):: source !< The fmsDiagTimeReduction_type object
    this%the_time_reduction = source%the_time_reduction
    this%time_averaging = source%time_averaging
    this%time_ops = source%time_ops
  end subroutine copy

  !> \brief Returns true if any of time_min, time_max, time_rms or time_average is true.
  !! @return true if any of time_min, time_max, time_rms or time_average is true.
  pure function has_time_ops_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: has_time_ops_imp
    has_time_ops_imp = this%time_ops
  end function has_time_ops_imp

  !> \brief Returns true iff time_averaging is true.
  !! @return true iff time_averaging is true.
  pure function do_time_averaging_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: do_time_averaging_imp
    do_time_averaging_imp = this%time_averaging
  end function do_time_averaging_imp

  !> \brief Returns true iff the_time_reduction is time_average
  !! @return true iff the_time_reduction is time_average
  pure function is_time_average_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_average_imp
    is_time_average_imp = this%the_time_reduction .EQ. time_average
  end function is_time_average_imp

  !> \brief Returns true iff the_time_reduction is time_none
  !! @return true iff the_time_reduction is time_none
  pure function is_time_none_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_none_imp
    is_time_none_imp = (this%the_time_reduction .EQ. time_none)
  end function is_time_none_imp

  !> \brief Returns true iff the_time_reduction is time_rms
  !! @return true iff the_time_reduction is time_rms
  pure function is_time_rms_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_rms_imp
    is_time_rms_imp = this%the_time_reduction .EQ. time_rms
  end function is_time_rms_imp

  !> \brief Returns true iff the_time_reduction is time_max
  !! @return true iff the_time_reduction is time_max
  pure function is_time_max_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_max_imp
    is_time_max_imp = this%the_time_reduction .EQ. time_max
  end function is_time_max_imp

  !> \brief Returns true iff the_time_reduction is time_min
  !! @return true iff the_time_reduction is time_min
  pure function is_time_min_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_min_imp
    is_time_min_imp = this%the_time_reduction .EQ. time_min
  end function is_time_min_imp

  !> \brief Returns true iff the_time_reduction is time_sum
  !! @return true iff the_time_reduction is time_sum
  pure function is_time_sum_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_sum_imp
    is_time_sum_imp = this%the_time_reduction .EQ. time_sum
  end function is_time_sum_imp

  !> \brief Returns true iff the_time_reduction is time_diurnal
  !! @return true iff the_time_reduction is time_diurnal
  pure function is_time_diurnal_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_diurnal_imp
    is_time_diurnal_imp = this%the_time_reduction .EQ. time_diurnal
  end function is_time_diurnal_imp

  !> \brief Returns true iff the_time_reduction is time_power
  !! @return true iff the_time_reduction is time_power
  pure function is_time_power_imp (this)
    class (fmsDiagTimeReduction_type), intent(in) :: this  !<The object this function is bound to.
    logical :: is_time_power_imp
    is_time_power_imp = this%the_time_reduction .EQ. time_power
  end function is_time_power_imp

end module fms_diag_time_reduction_mod
!> @}
! close documentation grouping
