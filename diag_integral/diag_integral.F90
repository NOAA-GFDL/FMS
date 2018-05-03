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

!> \file
!! \brief Contains the \ref diag_integral_mod module

                     module diag_integral_mod
#include <fms_platform.h>



!###############################################################################
!> \defgroup diag_integral_mod diag_integral_mod
!!
!! \author Fei Liu <Fei.Liu@noaa.gov>
!!
!! \brief This module computes and outputs global and / or hemispheric physics
!!        integrals.
!!
!! <b> Modules Included: </b>
!!
!! <table>
!!   <tr>
!!     <th> Module Name </th>
!!     <th> Included Values </th>
!!   </tr>
!!   <tr>
!!     <td> time_manager_mod </td>
!!     <td> time_type, get_time, set_time, time_manager_init, operator(+),
!!          operator(-), operator(==), operator(>=), operator(/=) </td>
!!   </tr>
!!   <tr>
!!     <td> mpp_mod </td>
!!     <td> input_nml_file <\td>
!!   </tr>
!!   <tr>
!!     <td> fms_mod </td>
!!     <td> open_file, file_exist, error_mesg, open_namelist_file,
!!          check_nml_error, fms_init, mpp_pe, mpp_root_pe, FATAL,
!!          write_version_number, stdlog, close_file </td>
!!   </tr>
!!   <tr>
!!     <td> constants_mod </td>
!!     <td> radius, constants_init </td>
!!   </tr>
!!   <tr>
!!     <td> mpp_mod </td>
!!     <td> mpp_sum, mpp_init </td>
!!   </tr>
!! </table>
!!
!! <b> Public Interfaces: </b>
!!
!! - sum_diag_integral_field
!!
!! <b> Public Subroutines: </b>
!!
!! - diag_integral_init
!! - diag_integral_field_init
!! - diag_integral_output
!! - diag_integral_end
!! - sum_field_2d
!! - sum_field_3d
!! - sum_field_wght_3d
!! - sum_field_2d_hemi
!!
!! <b> Private Functions: </b>
!!
!! - set_axis_time
!! - get_field_index
!! - get_axis_time
!! - diag_integral_alarm
!! - vert_diag_integral
!!
!! <b> Private Subroutines: </b>
!!
!! - write_field_averages
!! - format_text_init
!! - format_data_init
!!

use time_manager_mod, only:  time_type, get_time, set_time,  &
                             time_manager_init, &
                             operator(+),  operator(-),      &
                             operator(==), operator(>=),     &
                             operator(/=)
use mpp_mod,          only:  input_nml_file
use fms_mod,          only:  open_file, file_exist, error_mesg, &
                             open_namelist_file, check_nml_error, &
                             fms_init, &
                             mpp_pe, mpp_root_pe,&
                             FATAL, write_version_number, &
                             stdlog, close_file
use constants_mod,    only:  radius, constants_init
use mpp_mod,          only:  mpp_sum, mpp_init

!-------------------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!-------------------------------------------------------------------------------
!------ interfaces ------

public      &
          diag_integral_init, diag_integral_field_init, &
          sum_diag_integral_field, diag_integral_output,  &
          diag_integral_end



!###############################################################################
!> \defgroup sum_diag_integral_field sum_diag_integral_field
!! \ingroup diag_integral_mod
!!
!! This interface can be called in any one of three ways:
!!
!! \code{.f90}
!! call sum_diag_integral_field (name, data, is, js)
!! call sum_diag_integral_field (name, data, wt, is, js)
!! call sum_diag_integral_field (name, data, is, ie, js, je)
!! \endcode
!!
!! in the first option above, `data` may be either
!! \code{.f90}
!! real, intent(in) :: data(:,:)   ![ sum_field_2d ]
!! real, intent(in) :: data(:,:,:) ![ sum_field_3d ]
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: wt(:,:,:)
!! integer, optional, intent(in) :: is, ie, js, je
!! \endcode
!!
!! \param [in] <name> name associated with integral
!! \param [in] <data> field of integrands to be summed over
!! \param [in] <wt> vertical weighting factor to be applied to integrands
!!        when summing
!! \param [in] <is, ie, js, je> starting/ending i,j indices over which summation
!!        is to occur
!!
interface sum_diag_integral_field
   module procedure sum_field_2d,   &
                    sum_field_2d_hemi, &
                    sum_field_3d,   &
                    sum_field_wght_3d
end interface



private         &

!   from diag_integral_init:
          set_axis_time,  &

!   from diag_integral_field_init and sum_diag_integral_field:
          get_field_index, &

!   from diag_integral_output and diag_integral_end:
          write_field_averages,  &

!   from write_field_averages:
          format_text_init, format_data_init, &
          get_axis_time,     &

!   from diag_integral_output:
          diag_integral_alarm, &

!   from sum_diag_integral_field:
          vert_diag_integral

!-------------------------------------------------------------------------------
!------ namelist -------

integer, parameter  ::    &
                      mxch = 64    ! maximum number of characters in
                                   ! the optional output file name
real                ::    &
         output_interval = -1.0    ! time interval at which integrals
                                   ! are to be output
character(len=8)    ::    &
            time_units = 'hours'   ! time units associated with
                                   ! output_interval
character(len=mxch) ::    &
                 file_name = ' '   ! optional integrals output file name
logical             ::    &
           print_header = .true.   ! print a header for the integrals
                                   ! file ?
integer             ::    &
       fields_per_print_line = 4   ! number of fields to write per line
                                   ! of output


namelist / diag_integral_nml /      &
                                output_interval, time_units,  &
                                file_name, print_header, &
                                fields_per_print_line

!-------------------------------------------------------------------------------
!------- public data ------


!-------------------------------------------------------------------------------
!------- private data ------

!-------------------------------------------------------------------------------
!    variables associated with the determination of when integrals
!    are to be written.
!         Next_alarm_time  next time at which integrals are to be
!                          written
!         Alarm_interval   time interval between writing integrals
!         Zero_time        time_type variable set to (0,0); used as
!                          flag to indicate integrals are not being
!                          output
!         Time_init_save   initial time associated with experiment;
!                          used as a base for defining time
!-------------------------------------------------------------------------------
type (time_type) :: Next_alarm_time, Alarm_interval, Zero_time
type (time_type) :: Time_init_save

!-------------------------------------------------------------------------------
!    variables used in determining weights associated with each
!    contribution to the integrand.
!        area         area of each grid box
!        idim         x dimension of grid on local processor
!        jdim         y dimension of grid on local processor
!        field_size   number of columns on global domain
!        sum_area     surface area of globe
!-------------------------------------------------------------------------------
real, allocatable, dimension(:,:) :: area
integer                           :: idim, jdim, field_size
real                              :: sum_area

!-------------------------------------------------------------------------------
!    variables used to define the integral fields:
!      max_len_name     maximum length of name associated with integral
!      max_num_field    maximum number of integrals allowed
!      num_field        number of integrals that have been activated
!      field_name(i)    name associated with integral i
!      field_format(i)  output format for integral i
!      field_sum(i)     integrand for integral i
!      field_count(i)   number of values in integrand i
!-------------------------------------------------------------------------------
integer, parameter          :: max_len_name   = 12
integer, parameter          :: max_num_field = 32
integer                     :: num_field = 0
character(len=max_len_name) :: field_name   (max_num_field)
character(len=16)           :: field_format (max_num_field)
real                        :: field_sum    (max_num_field)
integer                     :: field_count  (max_num_field)

!-------------------------------------------------------------------------------
!    variables defining output formats.
!       format_text       format statement for header
!       format_data       format statement for data output
!       do_format_data    a data format needs to be generated ?
!       nd                number of characters in data format statement
!       nt                number of characters in text format statement
!-------------------------------------------------------------------------------
character(len=160) :: format_text, format_data
logical            :: do_format_data = .true.
integer            :: nd, nt

!-------------------------------------------------------------------------------
!    miscellaneous variables.
!-------------------------------------------------------------------------------
integer :: diag_unit = 0             ! unit number for output file
logical :: module_is_initialized = .false.
                                     ! module is initialized ?



                           contains



!###############################################################################
!
!                     PUBLIC SUBROUTINES
!
!###############################################################################



!###############################################################################
!> \fn diag_integral_init
!!
!! \brief diag_integral_init is the constructor for diag_integral_mod.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call diag_integral_init (Time_init, Time, blon, blat)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type (time_type),  intent(in), optional :: Time_init, Time
!! real,dimension(:,:), intent(in), optional :: blon, blat, area_in
!! \endcode
!!
!! \param [in] <Time_init> Initial time to start the integral
!! \param [in] <Time> current time
!! \param [in] <latb> array of model latitudes at cell boundaries [radians]
!! \param [in] <lonb> array of model longitudes at cell boundaries [radians]
!!
subroutine diag_integral_init (Time_init, Time, blon, blat, area_in)

type (time_type),  intent(in), optional :: Time_init, Time
real,dimension(:,:), intent(in), optional :: blon, blat, area_in

!-------------------------------------------------------------------------------
! local variables:
!       r2
!       rsize
!       unit
!       io
!       ierr
!       seconds
!       nc
!       i,j
!-------------------------------------------------------------------------------
      real    :: rsize
      integer :: unit, io, ierr, nc, logunit
      integer :: field_size_local
      real    :: sum_area_local

!-------------------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------------------
      call fms_init
      call mpp_init
      call constants_init
      call time_manager_init

!-------------------------------------------------------------------------------
!    if this is the initialization call, proceed. if this was simply
!    a verification of previous initialization, return.
!-------------------------------------------------------------------------------
      if (present(Time_init) .and. present(Time) .and. &
          present(blon) .and. present(blat) ) then

!-------------------------------------------------------------------------------
!    read namelist.
!-------------------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=diag_integral_nml, iostat=io)
        ierr = check_nml_error(io,'diag_integral_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=diag_integral_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'diag_integral_nml')
        end do
10      call close_file (unit)
#endif
      endif

!-------------------------------------------------------------------------------
!    write version number and namelist to logfile.
!-------------------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                       write (logunit, nml=diag_integral_nml)

!-------------------------------------------------------------------------------
!    save the initial time to time-stamp the integrals which will be
!    calculated.
!-------------------------------------------------------------------------------
      Time_init_save = Time_init

!-------------------------------------------------------------------------------
!    define the model grid sizes and the total number of columns on
!    the processor. sum over all processors and store the global
!    number of columns in field_size.
!-------------------------------------------------------------------------------
      idim = size(blon,1) - 1
      jdim = size(blon,2) - 1
      field_size_local = idim*jdim
      rsize = real(field_size_local)
      call mpp_sum (rsize)
      field_size = nint(rsize)

!-------------------------------------------------------------------------------
!    define an array to hold the surface area of each grid column
!    so that the integrals may be weighted properly. sum over the
!    processor, and then over all processors, storing the total
!    global surface area in sum_area.
!-------------------------------------------------------------------------------
      allocate (area(idim,jdim))

      area = area_in

      sum_area_local = sum(area)
      sum_area = sum_area_local
      call mpp_sum (sum_area)

!-------------------------------------------------------------------------------
!    if integral output is  to go to a file, open the file on unit
!    diag_unit.
!-------------------------------------------------------------------------------
      if (file_name(1:1) /= ' ' ) then
        nc = len_trim(file_name)
        diag_unit = open_file (file_name(1:nc), action='write')
      endif

!-------------------------------------------------------------------------------
!    define the variables needed to control the time interval of
!    output. Zero time is a flag indicating that the alarm is not set,
!    i.e., integrals are not desired.  otherwise set the next time to
!    output integrals to be at the value of nml variable
!    output_interval from now.
!-------------------------------------------------------------------------------
      Zero_time = set_time (0,0)
      if (output_interval >= -0.01) then
        Alarm_interval = set_axis_time (output_interval, time_units)
        Next_alarm_time = Time + Alarm_interval
      else
        Alarm_interval = Zero_time
      endif
      Next_alarm_time = Time + Alarm_interval

!-------------------------------------------------------------------------------
!    deallocate the local array and mark the module as initialized.
!-------------------------------------------------------------------------------
      module_is_initialized = .true.
   endif  ! (present optional arguments)

end subroutine diag_integral_init



!###############################################################################
!> \fn diag_integral_field_init
!!
!! \brief diag_integral_field_init registers and intializes an integral field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call diag_integral_field_init (name, format)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*), intent(in) :: name, format
!! \endcode
!!
!! \param [in] <name> Name of the field to be integrated
!! \param [in] <format> Output format of the field to be integrated
!!
subroutine diag_integral_field_init (name, format)

character(len=*), intent(in) :: name, format

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer :: field   ! index assigned to the current integral

!-------------------------------------------------------------------------------
!    note: no initialization is required for this interface. all needed
!    variables are initialized in the source.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!    make sure the integral name is not too long.
!-------------------------------------------------------------------------------
      if (len(name) > max_len_name )  then
        call error_mesg ('diag_integral_mod',  &
                ' integral name too long', FATAL)
      endif

!-------------------------------------------------------------------------------
!    check to be sure the integral name has not already been
!    initialized.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field /= 0)   then
        call error_mesg ('diag_integral_mod', &
                             'integral name already exists', FATAL)
      endif

!-------------------------------------------------------------------------------
!    prepare to register the integral. make sure that there are not
!    more integrals registered than space was provided for; if so, exit.
!-------------------------------------------------------------------------------
      num_field = num_field + 1
      if (num_field > max_num_field)  then
        call error_mesg ('diag_integral_mod', &
                              'too many fields initialized', FATAL)
      endif

!-------------------------------------------------------------------------------
!    register the name and output format desired for the given integral.
!    initialize its value and the number of grid points that have been
!    counted to zero.
!-------------------------------------------------------------------------------
      field_name   (num_field) = name
      field_format (num_field) = format
      field_sum    (num_field) = 0.0
      field_count  (num_field) = 0

end subroutine diag_integral_field_init



!###############################################################################
!> \fn sum_field_2d
!! \implements sum_diag_integral_field
!!
!! \brief Perform a 2 dimensional summation of named field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call sum_field_2d (name, data, is, js)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: data(:,:)
!! integer, optional, intent(in) :: is, js
!! \endcode
!!
!! \param [in] <name> Name of the field to be integrated
!! \param [in] <data> field of integrands to be summed over
!! \param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine sum_field_2d (name, data, is, js)

character(len=*),  intent(in) :: name
real,              intent(in) :: data(:,:)
integer, optional, intent(in) :: is, js

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer :: field           ! index of desired integral
      integer :: i1, j1, i2, j2  ! location indices of current data in
                                 ! processor-global coordinates

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)  then
        call error_mesg ('diag_integral_mod', &
                                    'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
     i1 = 1;  if (present(is)) i1 = is
     j1 = 1;  if (present(js)) j1 = js
     i2 = i1 + size(data,1) - 1
     j2 = j1 + size(data,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral and add the
!    values at this set of grid points to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count(field) +   &
                            size(data,1)*size(data,2)
      field_sum   (field) = field_sum   (field) +  &
                            sum (data * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine sum_field_2d



!###############################################################################
!> \fn sum_field_3d
!! \implements sum_diag_integral_field
!!
!! \brief Perform a 3 dimensional summation of named field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call sum_field_3d (name, data, is, js)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: data(:,:,:)
!! integer, optional, intent(in) :: is, js
!! \endcode
!!
!! \param [in] <name> Name of the field to be integrated
!! \param [in] <data> field of integrands to be summed over
!! \param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine sum_field_3d (name, data, is, js)

character(len=*),  intent(in) :: name
real,              intent(in) :: data(:,:,:)
integer, optional, intent(in) :: is, js

!-------------------------------------------------------------------------------
! local variables:
!     data2
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
      real, dimension (size(data,1),  &
                       size(data,2)) :: data2

      integer :: field
      integer :: i1, j1, i2, j2


!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)   then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
      i1 = 1;  if (present(is)) i1 = is
      j1 = 1;  if (present(js)) j1 = js
      i2 = i1 + size(data,1) - 1
      j2 = j1 + size(data,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. sum first
!    in the vertical and then add the values at this set of grid points
!    to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) +   &
                            size(data,1)*size(data,2)
      data2 = sum(data,3)
      field_sum   (field) = field_sum   (field) +  &
                            sum (data2 * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine sum_field_3d



!###############################################################################
!> \fn sum_field_wght_3d
!! \implements sum_diag_integral_field
!!
!! \brief Perform a 3 dimensional weighted summation of named field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call sum_field_wght_3d (name, data, wt, is, js)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: data(:,:,:), wt(:,:,:)
!! integer, optional, intent(in) :: is, js
!! \endcode
!!
!! \param [in] <name> Name of the field to be integrated
!! \param [in] <data> field of integrands to be summed over
!! \param [in] <wt> the weight function to be evaluated at summation
!! \param [in] <is, js> starting i,j indices over which summation is to occur
!!
subroutine sum_field_wght_3d (name, data, wt, is, js)

character(len=*),  intent(in) :: name
real,              intent(in) :: data(:,:,:), wt(:,:,:)
integer, optional, intent(in) :: is, js

!-------------------------------------------------------------------------------
! local variables:
!     data2
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
      real, dimension (size(data,1),size(data,2)) :: data2
      integer :: field, i1, j1, i2, j2

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)   then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!   define the processor-global indices of the current data. use the
!   value 1 for the initial grid points, if is and js are not input.
!-------------------------------------------------------------------------------
      i1 = 1;  if (present(is)) i1 = is
      j1 = 1;  if (present(js)) j1 = js
      i2 = i1 + size(data,1) - 1
      j2 = j1 + size(data,2) - 1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. sum first
!    in the vertical (including a vertical weighting factor) and then
!    add the values at this set of grid points to the accumulation
!    array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) +   &
                            size(data,1)*size(data,2)
      data2 = vert_diag_integral (data, wt)
      field_sum(field) = field_sum   (field) +  &
                         sum (data2 * area(i1:i2,j1:j2))

!$OMP END CRITICAL

end subroutine sum_field_wght_3d



!###############################################################################
!> \fn sum_field_2d_hemi
!! \implements sum_diag_integral_field
!!
!! \brief Perform a 2 dimensional hemispherical summation of named field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call sum_field_2d_hemi (name, data, is, ie, js, je)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! real,              intent(in) :: data(:,:)
!! integer,           intent(in) :: is, js, ie, je
!! \endcode
!!
!! \param [in] <name> Name of the field to be integrated
!! \param [in] <data> field of integrands to be summed over
!! \param [in] <is, js, ie, je> starting/ending i,j indices over which summation
!!        is to occur
!!
subroutine sum_field_2d_hemi (name, data, is, ie, js, je)

character(len=*),  intent(in) :: name
real,              intent(in) :: data(:,:)
integer,           intent(in) :: is, js, ie, je

!-------------------------------------------------------------------------------
! local variables:
!     field           ! index of desired integral
!     i1, j1, i2, j2  ! location indices of current data in
!                       processor-global coordinates
!-------------------------------------------------------------------------------
   integer :: field, i1, j1, i2, j2

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    obtain the index of the current integral. make certain it is valid.
!-------------------------------------------------------------------------------
      field = get_field_index (name)
      if (field == 0)    then
        call error_mesg ('diag_integral_mod', &
                               'field does not exist', FATAL)
      endif

!-------------------------------------------------------------------------------
!    define the processor-global indices of the current data. this form
!    is needed to handle case of 2d domain decomposition with physics
!    window smaller than processor domain size.
!-------------------------------------------------------------------------------
      i1 = mod ( (is-1), size(data,1) ) + 1
      i2 = i1 + size(data,1) - 1

!-------------------------------------------------------------------------------
!    for a hemispheric sum, sum one jrow at a time in case a processor
!    has data from both hemispheres.
!-------------------------------------------------------------------------------
      j1 = mod ( (js-1) ,size(data,2) ) + 1
      j2 = j1

!-------------------------------------------------------------------------------
!    increment the count of points toward this integral. include hemi-
!    spheric factor of 2 in field_count. add the data values at this
!    set of grid points to the accumulation array.
!-------------------------------------------------------------------------------
!$OMP CRITICAL
      field_count (field) = field_count (field) + 2* (i2-i1+1)*(j2-j1+1)
      field_sum   (field) = field_sum   (field) +  &
                            sum (data(i1:i2,j1:j2)*area(is:ie,js:je))

!$OMP END CRITICAL

end subroutine sum_field_2d_hemi



!###############################################################################
!> \fn diag_integral_output
!!
!! \brief diag_integral_output determines if this is a timestep on which
!!        integrals are to be written. if not, it returns; if so, it calls
!!        write_field_averages.
!!
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call diag_integral_output (Time)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type (time_type), intent(in) :: Time
!! \endcode
!!
!! \param [in] <Time> integral time stamp at the current time
!!
subroutine diag_integral_output (Time)

type (time_type), intent(in) :: Time

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    see if integral output is desired at this time.
!-------------------------------------------------------------------------------
      if ( diag_integral_alarm(Time) ) then

!-------------------------------------------------------------------------------
!    write the integrals by calling write_field_averages. upon return
!    reset the alarm to the next diagnostics time.
!-------------------------------------------------------------------------------
        call write_field_averages (Time)
        Next_alarm_time = Next_alarm_time + Alarm_interval
      endif

end subroutine diag_integral_output



!###############################################################################
!> \fn diag_integral_end
!!
!! \brief diag_integral_end is the destructor for diag_integral_mod.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call diag_integral_end (Time)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type (time_type), intent(in) :: Time
!! \endcode
!!
!! \param [in] <Time> integral time stamp at the current time
!!
subroutine diag_integral_end (Time)

type (time_type), intent(in) :: Time

!-------------------------------------------------------------------------------
!    be sure module has been initialized.
!-------------------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('diag_integral_mod',   &
              'module has not been initialized', FATAL )
      endif

!-------------------------------------------------------------------------------
!    if the alarm interval was set to Zero_time (meaning no integral
!    output during the model run) call write_field_averages to output
!    the integrals valid over the entire period of integration.
!-------------------------------------------------------------------------------
      if (Alarm_interval == Zero_time ) then
!       if (Alarm_interval /= Zero_time ) then
!       else
        call write_field_averages (Time)
      endif

!-------------------------------------------------------------------------------
!    deallocate module variables.
!-------------------------------------------------------------------------------
      deallocate (area)

!-------------------------------------------------------------------------------
!    mark the module as uninitialized.
!-------------------------------------------------------------------------------
      module_is_initialized = .false.

end subroutine diag_integral_end



!###############################################################################
!
!                    PRIVATE SUBROUTINES
!
!###############################################################################



!###############################################################################
!> \fn set_axis_time
!!
!! \brief Function to convert input time to a time_type
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! time = set_axis_time (atime, units)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! real,             intent(in) :: atime
!! character(len=*), intent(in) :: units
!! type(time_type)  :: Time
!! \endcode
!!
!! \param [in] <atime> integral time stamp at the current time
!! \param [in] <units> input units, not used
!! \param [out] <Time>
!!
function set_axis_time (atime, units) result (Time)

real,             intent(in) :: atime
character(len=*), intent(in) :: units
type(time_type)  :: Time

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer          :: sec     ! seconds corresponding to the input
                                  ! variable atime
      integer          :: day = 0 ! day component of time_type variable

!-------------------------------------------------------------------------------
!    convert the input time to seconds, regardless of input units.
!-------------------------------------------------------------------------------
      if (units(1:3) == 'sec') then
         sec = int(atime + 0.5)
      else if (units(1:3) == 'min') then
         sec = int(atime*60. + 0.5)
      else if (units(1:3) == 'hou') then
         sec = int(atime*3600. + 0.5)
      else if (units(1:3) == 'day') then
         sec = int(atime*86400. + 0.5)
      else
         call error_mesg('diag_integral_mod', &
                         'Invalid units sent to set_axis_time', FATAL)
      endif

!-------------------------------------------------------------------------------
!    convert the time in seconds to a time_type variable.
!-------------------------------------------------------------------------------
      Time = set_time (sec, day)

end function set_axis_time



!###############################################################################
!> \fn get_field_index
!!
!! \brief get_field_index returns returns the index associated with an
!!   integral name.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! index = get_field_index (name)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),  intent(in) :: name
!! integer                       :: index
!! \endcode
!!
!! \param [in] <name> Name associated with an integral
!! \param [out] <index>
!!
function get_field_index (name) result (index)

character(len=*),  intent(in) :: name
integer                       :: index

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer :: nc
      integer :: i

      nc = len_trim (name)
      if (nc > max_len_name)  then
        call error_mesg ('diag_integral_mod',  &
                                        'name too long', FATAL)
      endif

!-------------------------------------------------------------------------------
!    search each field name for the current string. when found exit
!    with the index. if not found index will be 0 upon return, which
!    initiates error condition.
!-------------------------------------------------------------------------------
      index = 0
      do i = 1, num_field
        if (name(1:nc) ==     &
                       field_name(i) (1:len_trim(field_name(i))) ) then
          index = i
          exit
        endif
      end do

end function get_field_index



!###############################################################################
!> \fn write_field_averages
!!
!! \brief Subroutine to sum multiple fields, average them and then write the
!!        result to an output file.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call  write_field_averages (Time)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type (time_type), intent(in) :: Time
!! \endcode
!!
!! \param [in] <Time> integral time stamp at the current time
!!
subroutine write_field_averages (Time)

type (time_type), intent(in) :: Time

!-------------------------------------------------------------------------------
! local variables:
!      field_avg
!      xtime
!      rcount
!      nn
!      ninc
!      nst
!      nend
!      fields_to_print
!      i
!      kount
!-------------------------------------------------------------------------------
      real    :: field_avg(max_num_field)
      real    :: xtime, rcount
      integer :: nn, ninc, nst, nend, fields_to_print
      integer :: i, kount
      integer(LONG_KIND) :: icount

!-------------------------------------------------------------------------------
!    each header and data format may be different and must be generated
!    as needed.
!-------------------------------------------------------------------------------
      fields_to_print = 0
      do i = 1, num_field

!-------------------------------------------------------------------------------
!    increment the fields_to_print counter.  sum the integrand and the
!    number of data points contributing to it over all processors.
!-------------------------------------------------------------------------------
        fields_to_print = fields_to_print + 1
        rcount = real(field_count(i))
        call mpp_sum (rcount)
        call mpp_sum (field_sum(i))
        icount = rcount

!-------------------------------------------------------------------------------
!    verify that all the data expected for an integral has been
!    obtained.
!-------------------------------------------------------------------------------
        if (icount == 0 ) call error_mesg &
                     ('diag_integral_mod',  &
                      'field_count equals zero for field_name ' //  &
                       field_name(i)(1:len_trim(field_name(i))), FATAL )
        kount = icount/field_size
        if ((field_size)*kount /= icount) then
           print*,"name,pe,kount,field_size,icount,rcount=",trim(field_name(i)),mpp_pe(),kount,field_size,icount,rcount
           call error_mesg &
                 ('diag_integral_mod',  &
                  'field_count not a multiple of field_size', FATAL )
        endif
!-------------------------------------------------------------------------------
!    define the global integral for field i. reinitialize the point
!    and data accumulators.
!-------------------------------------------------------------------------------
        field_avg(fields_to_print) = field_sum(i)/  &
                                     (sum_area*float(kount))
        field_sum  (i) = 0.0
        field_count(i) = 0
      end do

!-------------------------------------------------------------------------------
!    only the root pe will write out data.
!-------------------------------------------------------------------------------
      if ( mpp_pe() /= mpp_root_pe() ) return

!-------------------------------------------------------------------------------
!    define the time associated with the integrals just calculated.
!-------------------------------------------------------------------------------
      xtime = get_axis_time (Time-Time_init_save, time_units)

!-------------------------------------------------------------------------------
!    generate the new header and data formats.
!-------------------------------------------------------------------------------
      nst = 1
      nend = fields_per_print_line
      ninc = (num_field-1)/fields_per_print_line + 1
      do nn=1, ninc
        nst = 1 + (nn-1)*fields_per_print_line
        nend = MIN (nn*fields_per_print_line, num_field)
        if (print_header)  call format_text_init (nst, nend)
        call format_data_init (nst, nend)
        if (diag_unit /= 0) then
          write (diag_unit,format_data(1:nd)) &
                 xtime, (field_avg(i),i=nst,nend)
        else
          write (*, format_data(1:nd)) &
                 xtime, (field_avg(i),i=nst,nend)
        endif
      end do

end subroutine write_field_averages



!###############################################################################
!> \fn format_text_init
!!
!! \brief format_text_init generates the header records to be output in the
!!    integrals file.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call  format_text_init (nst_in, nend_in)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! integer, intent(in), optional :: nst_in, nend_in
!! \endcode
!!
!! \param [in] <nst_in, nend_in> starting/ending integral index which will be
!!        included in this format statement
!!
subroutine format_text_init (nst_in, nend_in)

integer, intent(in), optional :: nst_in, nend_in

!-------------------------------------------------------------------------------
! local variables:
!        i
!        nc
!        nst
!        nend
!-------------------------------------------------------------------------------
      integer :: i, nc, nst, nend

!-------------------------------------------------------------------------------
!    only the root pe need execute this routine, since only it will
!    be outputting integrals.
!-------------------------------------------------------------------------------
      if (mpp_pe() /= mpp_root_pe()) return

!-------------------------------------------------------------------------------
!    define the starting and ending integral indices that will be
!    included in this format statement.
!-------------------------------------------------------------------------------
      if (present (nst_in) ) then
        nst = nst_in
        nend = nend_in
      else
        nst = 1
        nend = num_field
      endif

!-------------------------------------------------------------------------------
!    define the first 11 characters in the format statement.
!-------------------------------------------------------------------------------
      nt = 11
      format_text(1:nt) = "('#    time"

!-------------------------------------------------------------------------------
!    generate the rest of the format statement, which will cover
!    integral indices nst to nend. if satndard printout is desired,
!    cycle through the loop.
!-------------------------------------------------------------------------------
      do i=nst,nend
        nc = len_trim(field_name(i))
        format_text(nt+1:nt+nc+5) =  '     ' // field_name(i)(1:nc)
        nt = nt+nc+5
      end do

!-------------------------------------------------------------------------------
!    include the end of the format statement.
!-------------------------------------------------------------------------------
      format_text(nt+1:nt+2) = "')"
      nt = nt+2

!-------------------------------------------------------------------------------
!    write the format statement to either an output file or to stdout.
!-------------------------------------------------------------------------------
      if (diag_unit /= 0) then
        write (diag_unit, format_text(1:nt))
      else
        write (*, format_text(1:nt))
      endif

end subroutine format_text_init



!###############################################################################
!> \fn format_data_init
!!
!! \brief format_text_init generates the format to be output in the
!!    integrals file.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! call  format_data_init (nst_in, nend_in)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! integer, intent(in), optional :: nst_in, nend_in
!! \endcode
!!
!! \param [in] <nst_in, nend_in> starting/ending integral index which will be
!!        included in this format statement
!!
subroutine format_data_init (nst_in, nend_in)

integer, intent(in), optional :: nst_in, nend_in

!-------------------------------------------------------------------------------
! local variables:
!        i
!        nc
!        nst
!        nend
!-------------------------------------------------------------------------------
      integer :: i, nc, nst, nend

!-------------------------------------------------------------------------------
!    define the start of the format, which covers the time stamp of the
!    integrals. this section is 9 characters long.
!-------------------------------------------------------------------------------
      nd = 9
      format_data(1:nd) = '(1x,f10.2'

!-------------------------------------------------------------------------------
!    define the indices of the integrals that are to be written by this
!    format statement.
!-------------------------------------------------------------------------------
      if ( present (nst_in) ) then
        nst = nst_in
        nend = nend_in
      else
        nst = 1
        nend = num_field
      endif

!-------------------------------------------------------------------------------
!    complete the data format. use the format defined for the
!    particular integral in setting up the format statement.
!-------------------------------------------------------------------------------
      do i=nst,nend
         nc = len_trim(field_format(i))
         format_data(nd+1:nd+nc+5) =  ',1x,' // field_format(i)(1:nc)
         nd = nd+nc+5
      end do

!-------------------------------------------------------------------------------
!    close the format statement.
!-------------------------------------------------------------------------------
      format_data(nd+1:nd+1) = ')'
      nd = nd + 1

end subroutine format_data_init



!###############################################################################
!> \fn get_axis_time
!!
!! \brief Function to convert the time_type input variable into units of
!!    units and returns it in atime.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! atime = get_axis_time (Time, units)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type(time_type),  intent(in) :: Time
!! character(len=*), intent(in) :: units
!! real                         :: atime
!! \endcode
!!
!! \param [in] <Time> integral time stamp
!! \param [in] <units> input units of time_type
!! \param [out] <atime>
!!
function get_axis_time (Time, units) result (atime)

type(time_type),  intent(in) :: Time
character(len=*), intent(in) :: units
real                         :: atime

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
      integer      :: sec, day  ! components of time_type variable

      call get_time (Time, sec, day)
      if (units(1:3) == 'sec') then
         atime = float(sec) + 86400.*float(day)
      else if (units(1:3) == 'min') then
         atime = float(sec)/60. + 1440.*float(day)
      else if (units(1:3) == 'hou') then
         atime = float(sec)/3600. + 24.*float(day)
      else if (units(1:3) == 'day') then
         atime = float(sec)/86400. + float(day)
      endif

end function get_axis_time



!###############################################################################
!> \fn diag_integral_alarm
!!
!! \brief Function to check if it is time to write integrals.
!!   if not writing integrals, return.
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! result = diag_integral_alarm (Time)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! type (time_type), intent(in) :: Time
!! logical                      :: answer
!! \endcode
!!
!! \param [in] <Time> current time
!! \param [out] <answer>
!!
function diag_integral_alarm (Time) result (answer)

type (time_type), intent(in) :: Time
logical                      :: answer

      answer = .false.
      if (Alarm_interval == Zero_time) return
      if (Time >= Next_alarm_time) answer = .true.

end function diag_integral_alarm



!###############################################################################
!> \fn vert_diag_integral
!!
!! \brief Function to perform a weighted integral in the vertical
!!        direction of a 3d data field
!!
!! <b> Template: </b>
!!
!! \code{.f90}
!! data2 = vert_diag_integral (data, wt)
!! \endcode
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! real, dimension (:,:,:),         intent(in) :: data, wt
!! real, dimension (size(data,1),size(data,2)) :: data2
!! \endcode
!!
!! \param [in] <data> integral field data arrays
!! \param [in] <wt> integral field weighting functions
!! \param [out] <data2>
function vert_diag_integral (data, wt) result (data2)

real, dimension (:,:,:),         intent(in) :: data, wt
real, dimension (size(data,1),size(data,2)) :: data2

!-------------------------------------------------------------------------------
! local variables:
!       wt2
!-------------------------------------------------------------------------------
      real, dimension(size(data,1),size(data,2)) :: wt2

!-------------------------------------------------------------------------------
      wt2 = sum(wt,3)
      if (count(wt2 == 0.) > 0)  then
        call error_mesg ('diag_integral_mod',  &
                             'vert sum of weights equals zero', FATAL)
      endif
      data2 = sum(data*wt,3) / wt2

end function vert_diag_integral





                    end module diag_integral_mod
