module diag_util_mod

use diag_data_mod, only  : output_fields, input_fields, files, do_diag_field_log, diag_log_unit, &
                           VERY_LARGE_AXIS_LENGTH, time_zero, VERY_LARGE_FILE_FREQ, END_OF_RUN, &
                           EVERY_TIME, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, &
                           DIAG_YEARS, base_time, time_unit_list, max_files, base_year, base_month, &
                           base_day, base_hour, base_minute, base_second, num_files, max_files, &
                           max_fields_per_file, max_out_per_in_field, max_input_fields,num_input_fields, &
                           max_output_fields, num_output_fields, coord_type, mix_snapshot_average_fields, &
                           global_descriptor
use diag_axis_mod, only  : get_diag_axis_data, get_axis_global_length, get_diag_axis_cart, &
                           get_domain1d, get_domain2d, diag_subaxes_init, diag_axis_init, get_diag_axis, &
                           get_axis_aux, get_axes_shift, get_diag_axis_name
use diag_output_mod, only: diag_flush, diag_field_out, diag_output_init, write_axis_meta_data, &
                           write_field_meta_data, done_meta_data
use fms_mod, only        : error_mesg, FATAL, WARNING, mpp_pe, mpp_root_pe, lowercase, fms_error_handler
use fms_io_mod, only     : get_tile_string, return_domain
use mpp_domains_mod,only : domain1d, domain2d, mpp_get_compute_domain, null_domain1d,&
                           null_domain2d, operator(/=), operator(==), mpp_modify_domain, mpp_get_domain_components, &
                           mpp_get_ntile_count, mpp_get_current_ntile, mpp_get_tile_id, mpp_mosaic_defined, &
                           mpp_get_tile_npes
use time_manager_mod,only: time_type, operator(==), operator(>), NO_CALENDAR, increment_date, &
                           increment_time, get_calendar_type, get_date, get_time, leap_year, &
                           operator(-),  operator(<), operator(>=)
use     mpp_io_mod, only : mpp_close
use     mpp_mod,    only : mpp_npes
use  constants_mod, only : SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE
implicit none
private

public get_subfield_size, log_diag_field_info, update_bounds, check_out_of_bounds, &
       check_bounds_are_exact_dynamic, check_bounds_are_exact_static, init_file, diag_time_inc, &
       find_input_field, init_input_field, init_output_field, diag_data_out, write_static, &
       check_duplicate_output_fields, get_date_dif

character(len=128),private  :: version = '$Id: diag_util.F90,v 16.0.2.2.2.1.2.2 2008/09/19 21:32:40 z1l Exp $'
character(len=128),private  :: tagname = '$Name: perth_2008_10 $'

contains

subroutine get_subfield_size(axes, outnum)
! Get size, start and end indices for output_fields(outnum), fill in
! output_fields(outnum)%output_grid%(start_indx, end_indx)

integer, intent(in) :: axes(:) ! axes of the input_field
integer, intent(in) :: outnum  ! position in array output_fields
real, allocatable   :: global_lat(:), global_lon(:), global_depth(:)
integer             :: global_axis_size
integer             :: i,xbegin,xend,ybegin,yend,xbegin_l,xend_l,ybegin_l,yend_l 
character(len=1)    :: cart
type(domain2d)      :: Domain2, Domain2_new
type(domain1d)      :: Domain1,Domain1x,Domain1y,Domain1x_new,Domain1y_new
real                :: start(3), end(3) ! start and end coordinates in 3 axes
integer             :: gstart_indx(3), gend_indx(3) ! global start and end indices of output domain in 3 axes 
real, allocatable   :: subaxis_x(:), subaxis_y(:), subaxis_z(:) !containing local coordinates in x,y,z axes
character(len=128)  :: msg
integer             :: ishift, jshift

!initilization for local output
start = -1.e10; end=-1.e10 ! initially out of (lat/lon/depth) range
gstart_indx = -1; gend_indx=-1

! get axis data (lat, lon, depth) and indices
   start= output_fields(outnum)%output_grid%start
   end = output_fields(outnum)%output_grid%end

do i = 1,size(axes(:))   
   global_axis_size = get_axis_global_length(axes(i))
   output_fields(outnum)%output_grid%subaxes(i) = -1
   call get_diag_axis_cart(axes(i), cart)
   select case(cart)
   case ('X')
      if(i.ne.1) &
           call error_mesg ('diag_util, get subfield size', 'wrong order of axes, X should come first',FATAL)
      allocate(global_lon(global_axis_size))
      call get_diag_axis_data(axes(i),global_lon)
      gstart_indx(i) = get_index(start(i),global_lon)
      gend_indx(i) = get_index(end(i),global_lon)
      allocate(subaxis_x(gstart_indx(i):gend_indx(i)))
      subaxis_x=global_lon(gstart_indx(i):gend_indx(i))   
   case ('Y')
      if(i.ne.2) &
           call error_mesg ('diag_util, get subfield size', 'wrong order of axes, Y should come second',FATAL)
      allocate(global_lat(global_axis_size))
      call get_diag_axis_data(axes(i),global_lat)
      gstart_indx(i) = get_index(start(i),global_lat)
      gend_indx(i) = get_index(end(i),global_lat)
      allocate(subaxis_y(gstart_indx(i):gend_indx(i)))
      subaxis_y=global_lat(gstart_indx(i):gend_indx(i))
   case ('Z')
      if(start(i)*end(i)<0) &
           call error_mesg ('diag_util, get subfield size','wrong values in vertical axis of region',FATAL)
      if(start(i)>0 .and. end(i)>0) then 
         allocate(global_depth(global_axis_size))
         call get_diag_axis_data(axes(i),global_depth)
         gstart_indx(i) = get_index(start(i),global_depth)
         gend_indx(i) = get_index(end(i),global_depth)
         allocate(subaxis_z(gstart_indx(i):gend_indx(i)))
         subaxis_z=global_depth(gstart_indx(i):gend_indx(i))
         output_fields(outnum)%output_grid%subaxes(i) = &
              diag_subaxes_init(axes(i),subaxis_z, gstart_indx(i),gend_indx(i))
         deallocate(subaxis_z,global_depth)
      else ! regional vertical axis is the same as global vertical axis
         gstart_indx(i) = 1
         gend_indx(i) = global_axis_size
         output_fields(outnum)%output_grid%subaxes(i) = axes(i)
         if(i /= 3) &
              call error_mesg ('diag_util, get subfield size','i should equal 3 for z axis', FATAL)
      endif      
   case default
       call error_mesg ('diag_util, get_subfield_size', 'Wrong axis_cart', FATAL)
   end select
enddo
do i = 1,size(axes(:))
   if(gstart_indx(i)== -1 .or. gend_indx(i)== -1) then
      write(msg,'(a,I2)') ' check region bounds for axis ', i
      call error_mesg ('diag_util, get_subfield_size', 'can not find gstart_indx/gend_indx for ' &
           //trim(output_fields(outnum)%output_name)//','//trim(msg), FATAL)
   endif
enddo

! get domain and compute_domain(xbegin,xend,ybegin,yend)
xbegin=-1; xend=-1
ybegin=-1; yend=-1

Domain2 = get_domain2d(axes)
if(Domain2 /= NULL_DOMAIN2D) then
   call mpp_get_compute_domain(Domain2,xbegin,xend,ybegin,yend)
   call mpp_get_domain_components(Domain2, Domain1x, Domain1y)
else
   do i = 1, MIN(size(axes(:)),2)    
      Domain1 = get_domain1d(axes(i))
      if(Domain1 /= NULL_DOMAIN1D) then
         call get_diag_axis_cart(axes(i),cart)
         select case(cart)
         case ('X')
            Domain1x = get_domain1d(axes(i))
            call mpp_get_compute_domain(Domain1x,xbegin,xend)   
         case ('Y')
            Domain1y = get_domain1d(axes(i))
            call mpp_get_compute_domain(Domain1y,ybegin,yend)
         case default ! do nothing here
         end select
      else
         call error_mesg ('diag_util, get_subfield_size', 'NO domain available', FATAL)
      endif
   enddo
endif      

call get_axes_shift(axes, ishift, jshift)
xend = xend+ishift
yend = yend+jshift

if(xbegin== -1 .or. xend==-1 .or. ybegin==-1 .or. yend==-1) &
   call error_mesg ('diag_util, get_subfield_size', 'wrong compute domain indices',FATAL)  

! get the area containing BOTH compute domain AND local output area
if(gstart_indx(1)> xend .or. xbegin > gend_indx(1)) then
   output_fields(outnum)%output_grid%l_start_indx(1) = -1
   output_fields(outnum)%output_grid%l_end_indx(1) = -1
   output_fields(outnum)%need_compute = .false. ! not involved
elseif (gstart_indx(2)> yend .or. ybegin > gend_indx(2)) then
   output_fields(outnum)%output_grid%l_start_indx(2) = -1
   output_fields(outnum)%output_grid%l_end_indx(2) = -1
   output_fields(outnum)%need_compute = .false. ! not involved
else
   output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1))
   output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2))
   output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1))
   output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2))
   output_fields(outnum)%need_compute = .true.  ! involved in local output
endif

if(output_fields(outnum)%need_compute) then
! need to modify domain1d and domain2d for subaxes
   xbegin_l = output_fields(outnum)%output_grid%l_start_indx(1)
   xend_l = output_fields(outnum)%output_grid%l_end_indx(1)
   ybegin_l = output_fields(outnum)%output_grid%l_start_indx(2)
   yend_l = output_fields(outnum)%output_grid%l_end_indx(2)
   call mpp_modify_domain(Domain2, Domain2_new, xbegin_l,xend_l, ybegin_l,yend_l, &
                          gstart_indx(1),gend_indx(1), gstart_indx(2),gend_indx(2))
   call mpp_get_domain_components(Domain2_new, Domain1x_new, Domain1y_new)

   output_fields(outnum)%output_grid%subaxes(1) = &
        diag_subaxes_init(axes(1),subaxis_x, gstart_indx(1),gend_indx(1),Domain1x_new,Domain2_new)
   output_fields(outnum)%output_grid%subaxes(2) = &
        diag_subaxes_init(axes(2),subaxis_y, gstart_indx(2),gend_indx(2),Domain1y_new,Domain2_new)
   do i = 1, size(axes(:))
      if(output_fields(outnum)%output_grid%subaxes(i) == -1) then  
         write(msg,'(a,"/",I4)') 'at i = ',i
         call error_mesg ('diag_util, get_subfield_size '//trim(output_fields(outnum)%output_name),&
              'error '//trim(msg), FATAL)   
      endif
   enddo
! local start index should start from 1
   output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1)) - xbegin + 1   
   output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2)) - ybegin + 1
   output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1)) - xbegin + 1 
   output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2)) - ybegin + 1
   if(size(axes(:))>2) then
      output_fields(outnum)%output_grid%l_start_indx(3) = gstart_indx(3)
      output_fields(outnum)%output_grid%l_end_indx(3) = gend_indx(3)
   else
      output_fields(outnum)%output_grid%l_start_indx(3) = 1
      output_fields(outnum)%output_grid%l_end_indx(3) = 1
   endif
endif
deallocate(subaxis_x, global_lon)
deallocate(subaxis_y, global_lat)

end subroutine get_subfield_size 
!---------------------------------------------------------------------------------------------------
function get_index(number, array)
  real, intent(in)               :: number
  real, intent(in), dimension(:) :: array
  integer                        :: get_index, i, n
  logical                        :: found
! Find index i of array such that array(i) is closest to number
! array must be  monotonouslly ordered

  n = size(array(:))
! check if array is monotonous
  do i = 2, n-1
     if((array(i-1)<array(i) .and. array(i)>array(i+1)).or.(array(i-1)>array(i) .and. array(i)<array(i+1))) &
          call error_mesg('diag_util', 'get_index, array NOT monotonously ordered',FATAL) 
  enddo
  get_index = -1
  found = .false.
! search in increasing array 
  do i = 1, n-1                
     if ((array(i)<=number) .and. (array(i+1)>= number)) then
        if(number - array(i) <= array(i+1) - number) then
           get_index = i
           found=.true.
        else
           get_index = i+1
           found=.true.
        endif
        exit
     endif     
  enddo 
! if not found, search in decreasing array
  if(.not.found) then
      do i = 1, n-1
         if ((array(i)>=number) .and. (array(i+1)<= number)) then
            if(array(i)-number <= number-array(i+1)) then
               get_index = i              
            else
               get_index = i+1               
            endif
            exit
         endif
      enddo
   endif
  end function get_index

! </FUNCTION>
!---------------------------------------------------------------------------------------------------
! <SUBROUTINE NAME="log_diag_field_info">
!   <OVERVIEW>
!     Writes brief diagnostic field info to the log file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If <TT>init_verbose</TT> namelist parameter is true, adds a line briefly 
!     describing diagnostic field to the log file. Normally users should not call
!     this subroutine directly, since it is called by register_static_field and register_diag_field
!     if do_not_log is not set to true. It is used, however, in LM3 to avoid excessive
!     log due to a number of fields registered for each of the tile types. LM3 code uses
!     do_not_log parameter in the registration calls, and calls this subroutine to
!     log field information under generic name.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call log_diag_field_info ( module_name, field_name, axes, long_name, units,
!     missing_value, range, dynamic )
!   </TEMPLATE>
subroutine log_diag_field_info ( module_name, field_name, axes, long_name, units, &
                                 missing_value, range, dynamic)
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
character(len=*), optional, intent(in) :: long_name, units
real   , optional, intent(in)          :: missing_value, range(2)
logical, optional, intent(in)          :: dynamic

! ---- local vars
character(len=256) :: lmodule, lfield, lname, lunits
character(len=64)  :: lmissval, lmin, lmax
character(len=8)   :: numaxis, timeaxis
character(len=1)   :: sep = '|'
character(len=256) :: axis_name, axes_list
integer :: i

if (.not.do_diag_field_log)    return
if (mpp_pe().ne.mpp_root_pe()) return

lmodule = trim(module_name)
lfield = trim(field_name)
lname  = ''; if(present(long_name)) lname  = trim(long_name)
lunits = ''; if(present(units))     lunits = trim(units)

write (numaxis,'(i1)') size(axes)

if (present(missing_value)) then
   write (lmissval,*) missing_value
else
   lmissval = ''
endif

if (present(range)) then
   write (lmin,*) range(1)
   write (lmax,*) range(2)
else
   lmin = ''
   lmax = ''
endif

if (present(dynamic)) then
   if (dynamic) then
      timeaxis = 'T'
   else
      timeaxis = 'F'
   endif
else
   timeaxis = ''
endif

axes_list=''
do i = 1,size(axes)
   call get_diag_axis_name(axes(i),axis_name)
   if(trim(axes_list)/='')axes_list=trim(axes_list)//','
   axes_list=trim(axes_list)//trim(axis_name)
enddo

!write (diag_log_unit,'(8(a,a),a)') &
write (diag_log_unit,'(777a)') &
             trim(lmodule),  sep, trim(lfield),  sep, trim(lname),    sep, &
             trim(lunits),   sep, trim(numaxis), sep, trim(timeaxis), sep, &
             trim(lmissval), sep, trim(lmin),    sep, trim(lmax),     sep, &
             trim(axes_list)

end subroutine log_diag_field_info
! </SUBROUTINE>
!===================================================================================================
 subroutine update_bounds(out_num, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k)
 integer, intent(in) :: out_num, lower_i, upper_i, lower_j, upper_j, lower_k, upper_k

 output_fields(out_num)%imin = min(output_fields(out_num)%imin, lower_i)
 output_fields(out_num)%imax = max(output_fields(out_num)%imax, upper_i)
 output_fields(out_num)%jmin = min(output_fields(out_num)%jmin, lower_j)
 output_fields(out_num)%jmax = max(output_fields(out_num)%jmax, upper_j)
 output_fields(out_num)%kmin = min(output_fields(out_num)%kmin, lower_k)
 output_fields(out_num)%kmax = max(output_fields(out_num)%kmax, upper_k)

end subroutine update_bounds
!===================================================================================================
 subroutine check_out_of_bounds(out_num, diag_field_id, err_msg)
 integer, intent(in) :: out_num, diag_field_id
 character(len=*), intent(out) :: err_msg
 character(len=128)  :: error_string1, error_string2

 if(output_fields(out_num)%imin < lbound(output_fields(out_num)%buffer,1) .or. &
    output_fields(out_num)%imax > ubound(output_fields(out_num)%buffer,1) .or. &
    output_fields(out_num)%jmin < lbound(output_fields(out_num)%buffer,2) .or. &
    output_fields(out_num)%jmax > ubound(output_fields(out_num)%buffer,2) .or. &
    output_fields(out_num)%kmin < lbound(output_fields(out_num)%buffer,3) .or. &
    output_fields(out_num)%kmax > ubound(output_fields(out_num)%buffer,3) ) then
    write(error_string1,'(a,"/",a)') trim(input_fields(diag_field_id)%module_name),trim(output_fields(out_num)%output_name)
    error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
    write(error_string2(15:17),'(i3)') lbound(output_fields(out_num)%buffer,1)
    write(error_string2(19:21),'(i3)') ubound(output_fields(out_num)%buffer,1)
    write(error_string2(23:25),'(i3)') lbound(output_fields(out_num)%buffer,2)
    write(error_string2(27:29),'(i3)') ubound(output_fields(out_num)%buffer,2)
    write(error_string2(31:33),'(i3)') lbound(output_fields(out_num)%buffer,3)
    write(error_string2(35:37),'(i3)') ubound(output_fields(out_num)%buffer,3)
    write(error_string2(54:56),'(i3)') output_fields(out_num)%imin
    write(error_string2(58:60),'(i3)') output_fields(out_num)%imax
    write(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
    write(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
    write(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
    write(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
    err_msg = 'module/output_field='//trim(error_string1)//'  Bounds of buffer exceeded.  '//trim(error_string2)
!   imax, imin, etc need to be reset in case the program is not terminated.
    output_fields(out_num)%imax = 0
    output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%jmax = 0
    output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
    output_fields(out_num)%kmax = 0
    output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
 else
    err_msg = ''
 endif

 end subroutine check_out_of_bounds
!===================================================================================================
 subroutine check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg)
 integer,          intent(in)  :: out_num, diag_field_id
 type(time_type),  intent(in)  :: Time
 character(len=*), intent(out) :: err_msg
 character(len=128)            :: error_string1, error_string2
 logical                       :: do_check

 err_msg = ''

! Check bounds only when the value of Time changes. When windows are used, a change in Time indicates
! that a new loop through the windows has begun, so a check of the previous loop can be done.

 if(Time == output_fields(out_num)%Time_of_prev_field_data) then
    do_check = .false.
 else
    if(output_fields(out_num)%Time_of_prev_field_data == Time_zero) then
       do_check = .false. ! It may or may not be OK to check, I don't know how to tell. Check will be done on subsequent calls anyway.
    else
       do_check = .true.
    endif
    output_fields(out_num)%Time_of_prev_field_data = Time
 endif

 if(do_check) then
   if(output_fields(out_num)%imin /= lbound(output_fields(out_num)%buffer,1) .or. &
      output_fields(out_num)%imax /= ubound(output_fields(out_num)%buffer,1) .or. &
      output_fields(out_num)%jmin /= lbound(output_fields(out_num)%buffer,2) .or. &
      output_fields(out_num)%jmax /= ubound(output_fields(out_num)%buffer,2) .or. &
      output_fields(out_num)%kmin /= lbound(output_fields(out_num)%buffer,3) .or. &
      output_fields(out_num)%kmax /= ubound(output_fields(out_num)%buffer,3) ) then
        write(error_string1,'(a,"/",a)') trim(input_fields(diag_field_id)%module_name),trim(output_fields(out_num)%output_name)
        error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
        write(error_string2(15:17),'(i3)') lbound(output_fields(out_num)%buffer,1)
        write(error_string2(19:21),'(i3)') ubound(output_fields(out_num)%buffer,1)
        write(error_string2(23:25),'(i3)') lbound(output_fields(out_num)%buffer,2)
        write(error_string2(27:29),'(i3)') ubound(output_fields(out_num)%buffer,2)
        write(error_string2(31:33),'(i3)') lbound(output_fields(out_num)%buffer,3)
        write(error_string2(35:37),'(i3)') ubound(output_fields(out_num)%buffer,3)
        write(error_string2(54:56),'(i3)') output_fields(out_num)%imin
        write(error_string2(58:60),'(i3)') output_fields(out_num)%imax
        write(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
        write(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
        write(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
        write(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
        err_msg = trim(error_string1)//' Bounds of data do not match those of buffer. '//trim(error_string2)
   endif
   output_fields(out_num)%imax = 0
   output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
   output_fields(out_num)%jmax = 0
   output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
   output_fields(out_num)%kmax = 0
   output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
 endif

 end subroutine check_bounds_are_exact_dynamic
!===================================================================================================
 subroutine check_bounds_are_exact_static(out_num, diag_field_id, err_msg)
 integer, intent(in) :: out_num, diag_field_id
 character(len=*), intent(out) :: err_msg
 character(len=128)  :: error_string1, error_string2

 err_msg = ''

 if(output_fields(out_num)%imin /= lbound(output_fields(out_num)%buffer,1) .or. &
    output_fields(out_num)%imax /= ubound(output_fields(out_num)%buffer,1) .or. &
    output_fields(out_num)%jmin /= lbound(output_fields(out_num)%buffer,2) .or. &
    output_fields(out_num)%jmax /= ubound(output_fields(out_num)%buffer,2) .or. &
    output_fields(out_num)%kmin /= lbound(output_fields(out_num)%buffer,3) .or. &
    output_fields(out_num)%kmax /= ubound(output_fields(out_num)%buffer,3) ) then
      write(error_string1,'(a,"/",a)') trim(input_fields(diag_field_id)%module_name),trim(output_fields(out_num)%output_name)
      error_string2 ='Buffer bounds=   :   ,   :   ,   :     Actual bounds=   :   ,   :   ,   :   '
      write(error_string2(15:17),'(i3)') lbound(output_fields(out_num)%buffer,1)
      write(error_string2(19:21),'(i3)') ubound(output_fields(out_num)%buffer,1)
      write(error_string2(23:25),'(i3)') lbound(output_fields(out_num)%buffer,2)
      write(error_string2(27:29),'(i3)') ubound(output_fields(out_num)%buffer,2)
      write(error_string2(31:33),'(i3)') lbound(output_fields(out_num)%buffer,3)
      write(error_string2(35:37),'(i3)') ubound(output_fields(out_num)%buffer,3)
      write(error_string2(54:56),'(i3)') output_fields(out_num)%imin
      write(error_string2(58:60),'(i3)') output_fields(out_num)%imax
      write(error_string2(62:64),'(i3)') output_fields(out_num)%jmin
      write(error_string2(66:68),'(i3)') output_fields(out_num)%jmax
      write(error_string2(70:72),'(i3)') output_fields(out_num)%kmin
      write(error_string2(74:76),'(i3)') output_fields(out_num)%kmax
      err_msg = trim(error_string1)//' Bounds of data do not match those of buffer. '//trim(error_string2)
   endif
   output_fields(out_num)%imax = 0
   output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
   output_fields(out_num)%jmax = 0
   output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
   output_fields(out_num)%kmax = 0
   output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH

 end subroutine check_bounds_are_exact_static
!===================================================================================================
subroutine init_file(name, output_freq, output_units, format, time_units, long_name, &
     new_file_freq, new_file_freq_units, start_time,file_duration,file_duration_units)

character(len=*), intent(in)          :: name, long_name
integer, intent(in)                   :: output_freq, output_units, format, time_units
integer, intent(in), optional         :: new_file_freq, new_file_freq_units
integer, intent(in), optional         :: file_duration, file_duration_units
type(time_type), intent(in), optional :: start_time
integer                               :: new_file_freq1, new_file_freq_units1
integer                               :: file_duration1, file_duration_units1
real, dimension(1)                    :: tdata
character(len=128)                    :: time_units_str
! Get a number for this file
num_files = num_files + 1
if(num_files >= max_files) then
   call error_mesg('diag_util, init_file', ' max_files exceeded, incease max_files', FATAL)
endif

new_file_freq1 = VERY_LARGE_FILE_FREQ
if(present(new_file_freq)) new_file_freq1 = new_file_freq
if (get_calendar_type() == NO_CALENDAR) then
  new_file_freq_units1 = DIAG_DAYS
else 
   new_file_freq_units1 = DIAG_YEARS
endif
if(present(new_file_freq_units)) new_file_freq_units1 = new_file_freq_units 
file_duration1 = new_file_freq1
if(present(file_duration)) file_duration1 = file_duration 
file_duration_units1 = new_file_freq_units1
if(present(file_duration_units)) file_duration_units1 = file_duration_units

files(num_files)%name = trim(name)
files(num_files)%output_freq = output_freq
files(num_files)%output_units = output_units
files(num_files)%format = format
files(num_files)%time_units = time_units
files(num_files)%long_name = trim(long_name)
files(num_files)%num_fields = 0
files(num_files)%local = .false.
files(num_files)%last_flush = base_time
files(num_files)%file_unit = -1
files(num_files)%new_file_freq = new_file_freq1
files(num_files)%new_file_freq_units = new_file_freq_units1
files(num_files)%duration = file_duration1
files(num_files)%duration_units = file_duration_units1
if(present(start_time)) then 
   files(num_files)%start_time = start_time
else
   files(num_files)%start_time = base_time
endif
files(num_files)%next_open=diag_time_inc(files(num_files)%start_time,new_file_freq1,new_file_freq_units1)
files(num_files)%close_time = diag_time_inc(files(num_files)%start_time,file_duration1, &
     file_duration_units1)
if(files(num_files)%close_time>files(num_files)%next_open) call error_mesg('init_file', 'close time '//&
     'GREATER than'//' next_open time, check file duration, file frequency in '// &
     files(num_files)%name, FATAL)
! add time_axis_id and time_bounds_id here
write(time_units_str, 11) trim(time_unit_list(files(num_files)%time_units)), base_year, &
     base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
files(num_files)%time_axis_id = diag_axis_init (trim(long_name), tdata, time_units_str, 'T', &
     trim(long_name) , set_name=trim(name) )
!---- register axis for storing time boundaries
files(num_files)%time_bounds_id = diag_axis_init( 'nv',(/1.,2./),'none','N','vertex number', &
     set_name=trim(name))
end subroutine init_file

!---------------------------------------------------------------------------------------------------
function diag_time_inc(time, output_freq, output_units, err_msg)

type (time_type)                        :: diag_time_inc
type (time_type), intent(in)            :: time
integer, intent(in)                     :: output_freq, output_units
character(len=*), intent(out), optional :: err_msg
character(len=128)                      :: error_message_local

if(present(err_msg)) err_msg = ''
error_message_local = ''

! special values for output frequency are -1 for output at end of run
! and 0 for every timestep.  Need to check for these here?
! Return zero time increment, hopefully this value is never used

if (output_freq == END_OF_RUN .or. output_freq == EVERY_TIME) then
    diag_time_inc = time
    return
endif

! Make sure calendar was not set after initialization

if(output_units == DIAG_SECONDS) then
   if (get_calendar_type() == NO_CALENDAR) then
       diag_time_inc = increment_time(time, output_freq, 0, err_msg=error_message_local)
   else
       diag_time_inc = increment_date(time, 0, 0, 0, 0, 0, output_freq, err_msg=error_message_local)
   endif
else if(output_units == DIAG_MINUTES) then
   if (get_calendar_type() == NO_CALENDAR) then
       diag_time_inc = increment_time(time, nint(output_freq*SECONDS_PER_MINUTE), 0, err_msg=error_message_local)
   else
       diag_time_inc = increment_date(time, 0, 0, 0, 0, output_freq, 0, err_msg=error_message_local)
   endif
else if(output_units == DIAG_HOURS) then
   if (get_calendar_type() == NO_CALENDAR) then
       diag_time_inc = increment_time(time, nint(output_freq*SECONDS_PER_HOUR), 0, err_msg=error_message_local)
   else
       diag_time_inc = increment_date(time, 0, 0, 0, output_freq, 0, 0, err_msg=error_message_local)
   endif
else if(output_units == DIAG_DAYS) then
   if (get_calendar_type() == NO_CALENDAR) then
       diag_time_inc = increment_time(time, 0, output_freq, err_msg=error_message_local)
   else
       diag_time_inc = increment_date(time, 0, 0, output_freq, 0, 0, 0, err_msg=error_message_local)
   endif
else if(output_units == DIAG_MONTHS) then
   if (get_calendar_type() == NO_CALENDAR) then
       error_message_local = 'output units of months NOT allowed with no calendar'
   else
       diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0, err_msg=error_message_local)
   endif
else if(output_units == DIAG_YEARS) then
   if (get_calendar_type() == NO_CALENDAR) then
       error_message_local = 'output units of years NOT allowed with no calendar'
   else
       diag_time_inc = increment_date(time, output_freq, 0, 0, 0, 0, 0, err_msg=error_message_local)
   endif
else 
    error_message_local = 'illegal output units'
endif

if(error_message_local /= '') then
  if(fms_error_handler('diag_time_inc',error_message_local,err_msg)) return
endif

end function diag_time_inc

!---------------------------------------------------------------------------------------------------
function find_file(name)
integer                      :: find_file
character(len=*), intent(in) :: name
integer                      :: i

find_file = -1
do i = 1, num_files
   if(trim(files(i)%name) == trim(name)) then
      find_file = i
      return
   end if
end do
end function find_file
!---------------------------------------------------------------------------------------------------

function find_input_field(module_name, field_name)
integer                      :: find_input_field
character(len=*), intent(in) :: module_name, field_name
integer                      :: i

find_input_field = -1
do i = 1, num_input_fields
   if(trim(input_fields(i)%module_name) == trim(module_name) .and. &
      lowercase(trim(input_fields(i)%field_name)) == &
      lowercase(trim(field_name))) then 
      find_input_field = i
      return
   endif
end do
end function find_input_field
!---------------------------------------------------------------------------------------------------

subroutine init_input_field(module_name, field_name)

character(len=*), intent(in) :: module_name, field_name

! Get a number for this input_field if not already set up
if(find_input_field(module_name, field_name) < 0) then
   num_input_fields = num_input_fields + 1
   if(num_input_fields > max_input_fields) then
      call error_mesg('diag_util,init_input_field', 'max_input_fields exceeded, increase it via diag_manager_nml',&
           FATAL)
   end if
else
! If this is already initialized don't need to do anything
   return
end if

input_fields(num_input_fields)%module_name = trim(module_name)
input_fields(num_input_fields)%field_name = trim(field_name)
input_fields(num_input_fields)%num_output_fields = 0
! Set flag that this field has not been registered
input_fields(num_input_fields)%register = .false.
input_fields(num_input_fields)%local = .false.
input_fields(num_input_fields)%standard_name = 'none'
end subroutine init_input_field

!---------------------------------------------------------------------------------------------------

subroutine init_output_field(module_name, field_name, output_name, output_file,&
   time_method, pack, local_coord)

character(len=*), intent(in)           :: module_name, field_name, output_name, output_file
character(len=*), intent(in)           :: time_method
integer, intent(in)                    :: pack
type(coord_type), intent(in), optional :: local_coord
integer                                :: out_num, in_num, file_num, &
                                          num_fields, i, method_selected, l1
character(len=128)                     :: error_msg
character(len=8)                       :: t_method

! Get a number for this output field
num_output_fields = num_output_fields + 1
if(num_output_fields > max_output_fields) then
   call error_mesg('diag_util', 'max_output_fields exceeded, increase it via diag_manager_nml', FATAL)
endif
out_num = num_output_fields

! First, find the index to the associated input field
in_num = find_input_field(module_name, field_name)
if(in_num < 0) then
   write (error_msg,'(a,"/",a)') trim(module_name),trim(field_name)
   call error_mesg('diag_util,init_output_field', &
      'module_name/field_name '//trim(error_msg)//' NOT registered', FATAL)
endif

! Add this output field into the list for this input field
input_fields(in_num)%num_output_fields = &
   input_fields(in_num)%num_output_fields + 1
if(input_fields(in_num)%num_output_fields > max_out_per_in_field) then
   call error_mesg('diag_util,init_output_field', 'max_out_per_in_field exceeded, increase max_out_per_in_field', FATAL)
endif
input_fields(in_num)%output_fields(input_fields(in_num)%num_output_fields) &
   = out_num

! Also put pointer to input field in this output field
output_fields(out_num)%input_field = in_num

! Next, find the number for the corresponding file
if(trim(output_file).eq.'null') then
   file_num = max_files
else
   file_num = find_file(output_file)
   if(file_num < 0) then
      call error_mesg('diag_util,init_output_field', 'file '//trim(output_file)//' is NOT found in diag_table', FATAL)
   end if
endif

! Insert this field into list for this file
files(file_num)%num_fields = files(file_num)%num_fields + 1
if(files(file_num)%num_fields > max_fields_per_file) then
   call error_mesg('diag_util,init_output_field', 'max_fields_per_file exceeded, increase max_fields_per_file ', FATAL)
endif
num_fields = files(file_num)%num_fields
files(file_num)%fields(num_fields) = out_num
output_fields(out_num)%count_0d = 0

! Set the file for this output field
output_fields(out_num)%output_file = file_num

! Enter the other data for this output field
output_fields(out_num)%output_name = trim(output_name)
output_fields(out_num)%pack = pack
output_fields(out_num)%num_elements = 0
output_fields(out_num)%total_elements = 0
output_fields(out_num)%region_elements = 0
output_fields(out_num)%imax = 0
output_fields(out_num)%jmax = 0
output_fields(out_num)%kmax = 0
output_fields(out_num)%imin = VERY_LARGE_AXIS_LENGTH
output_fields(out_num)%jmin = VERY_LARGE_AXIS_LENGTH
output_fields(out_num)%kmin = VERY_LARGE_AXIS_LENGTH
method_selected = 0
! Initialize all time method to false
output_fields(out_num)%time_average = .false.
output_fields(out_num)%time_min = .false.
output_fields(out_num)%time_max = .false. 
output_fields(out_num)%time_ops = .false.
output_fields(out_num)%written_once = .false.
! cannot time average fields output every time
if (files(file_num)%output_freq == EVERY_TIME) then
   output_fields(out_num)%time_average = .false.
   method_selected = method_selected+1
   t_method = 'point'
else
   t_method = lowercase(time_method)
   select case (trim(t_method))
   case('.true.')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('mean')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('average')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('avg')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('.false.')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'
   case('none')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'
   case('point')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'   
   case ('maximum')
      output_fields(out_num)%time_max = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'max') &
           output_fields(out_num)%output_name = trim(output_name)//'_max'
      method_selected = method_selected+1
      t_method = 'max'        
   case ('max')
      output_fields(out_num)%time_max = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'max') &
           output_fields(out_num)%output_name = trim(output_name)//'_max'
      method_selected = method_selected+1
      t_method = 'max'
   case ('minimum')
      output_fields(out_num)%time_min = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'min') &
           output_fields(out_num)%output_name = trim(output_name)//'_min'
      method_selected = method_selected+1
      t_method = 'min'        
   case ('min')
      output_fields(out_num)%time_min = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'min') &
           output_fields(out_num)%output_name = trim(output_name)//'_min'
      method_selected = method_selected+1
      t_method = 'min'     
   end select
endif
output_fields(out_num)%time_ops = output_fields(out_num)%time_min.or.output_fields(out_num)%time_max &
     .or.output_fields(out_num)%time_average

output_fields(out_num)%phys_window = .false.
! need to initialize grid_type = -1(start, end, l_start_indx,l_end_indx etc...)
if(present(local_coord)) then
   input_fields(in_num)%local = .true.
   output_fields(out_num)%output_grid%start(1) = local_coord%xbegin
   output_fields(out_num)%output_grid%start(2) = local_coord%ybegin
   output_fields(out_num)%output_grid%start(3) = local_coord%zbegin
   output_fields(out_num)%output_grid%end(1) = local_coord%xend
   output_fields(out_num)%output_grid%end(2) = local_coord%yend
   output_fields(out_num)%output_grid%end(3) = local_coord%zend
   do i = 1,3
      output_fields(out_num)%output_grid%l_start_indx(i) = -1
      output_fields(out_num)%output_grid%l_end_indx(i) = -1
      output_fields(out_num)%output_grid%subaxes(i) = -1
   enddo
   output_fields(out_num)%local_output = .true.
   output_fields(out_num)%need_compute = .false.
else
   output_fields(out_num)%local_output = .false.
   output_fields(out_num)%need_compute = .false.
endif

if (method_selected /= 1) call error_mesg('init_output_field','improper &
     &time method in diag_table for output field:'//trim(output_name),FATAL)
output_fields(out_num)%time_method = trim(t_method)

end subroutine init_output_field

!---------------------------------------------------------------------------------------------------

subroutine opening_file(file, time)

! WARNING: Assumes that all data structures are fully initialized
  integer, intent(in)           :: file
  type(time_type), intent(in)   :: time  
  integer                       :: j, field_num, axes(5), input_field_num, num_axes, k
  character(len=128)            :: time_units, timeb_units, avg, error_string, filename, aux_name,fieldname
  logical                       :: time_ops, aux_present, match_aux_name
  integer                       :: time_axis_id(1), field_num1, time_bounds_id(1)
  character(len=7)              :: prefix
  character (len = 7)           :: avg_name = 'average'
  character(len=128)            :: suffix, base_name
  integer                       :: position
  character(len=32)             :: time_name, timeb_name,time_longname, timeb_longname, cart_name
  type(domain1d)                :: domain
  integer                       :: dir, edges
  real, dimension(2)            :: data
  character(len=256)            :: fname
  integer                       :: ntileMe, nfiles_in_set
  integer, allocatable          :: tile_id(:)
  type(domain2d)                :: domain2
  logical                       :: all_scalar_or_1d


  aux_present = .false.; match_aux_name = .false.
! it's unlikely that a file starts with word "rregion", need to check anyway.
  if (len(files(file)%name) >=7 .and. .not. files(file)%local) then
     prefix = files(file)%name(1:7)
     if(lowercase(prefix) == 'rregion') &
          call error_mesg ('diag_util opening_file', 'file name should not start with' &
          //' word "rregion"', WARNING)
  endif
  
! Here is where time_units string must be set up; time since base date
  write(time_units, 11) trim(time_unit_list(files(file)%time_units)), base_year, &
       base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
  base_name = files(file)%name
  if(files(file)%new_file_freq < VERY_LARGE_FILE_FREQ) then
     position = INDEX(files(file)%name, '%')
     if(position>0) then
        base_name = base_name(1:position-1)
     else
        call error_mesg ('diag_util opening_file', 'file name '//TRIM(files(file)%name)// &
             ' does not contain % for time stamp string', FATAL) 
     endif
     suffix = get_time_string(files(file)%name, time)
  else
     suffix = ' '
  endif 
! Add CVS tag as prefix of filename  (currently not implemented)
!  i1 = INDEX(tagname,':') + 2
!  i2 = len_trim(tagname) - 2
!  if(i2 <=i1)  call error_mesg ('diag_util opening_file','error in CVS tagname index',FATAL)
!  prefix2 = tagname(i1:i2)//'_'
  if(files(file)%local) then      
! prepend "rregion" to all local files for post processing, the prefix will be removed in postprocessing
     filename = 'rregion'//trim(base_name)//trim(suffix)
  else
!     filename = trim(prefix2)//trim(base_name)//trim(suffix)
     filename = trim(base_name)//trim(suffix)
  endif

! Loop through all fields with this file to output axes
! JWD: This is a klooge; need something more robust
  nfiles_in_set = mpp_npes()

  domain2 = NULL_DOMAIN2D
  all_scalar_or_1d = .true.
  do j = 1, files(file)%num_fields
     field_num = files(file)%fields(j)
     num_axes = output_fields(field_num)%num_axes
     if(num_axes>1)then
        all_scalar_or_1d = .false.
        domain2 = get_domain2d ( output_fields(field_num)%axes(1:num_axes) )
        if(domain2 /= NULL_DOMAIN2D) exit
     endif
  enddo
  if(.not.all_scalar_or_1d) then
     if(domain2 == NULL_DOMAIN2D) call return_domain(domain2)
     if(domain2 == NULL_DOMAIN2D)then
        call error_mesg ('diag_util opening_file', &
             'Domain not defined through set_domain interface; cannot retrieve tile info', FATAL)
     endif
     nfiles_in_set = mpp_get_tile_npes(domain2)
     if(mpp_get_ntile_count(domain2) > 1)then
        ntileMe = mpp_get_current_ntile(domain2)
        allocate(tile_id(ntileMe))
        tile_id = mpp_get_tile_id(domain2)
        fname = trim(filename)
        call get_tile_string(filename, trim(fname)//'.tile' , tile_id(1))
        deallocate(tile_id)
     endif
  endif

  call diag_output_init(filename, files(file)%format, global_descriptor, &
       files(file)%long_name, time_units, files(file)%file_unit, nfiles_in_set, all_scalar_or_1d) 
  files(file)%bytes_written = 0 
! Does this file contain time_average fields?
  time_ops = .false.
  do j = 1, files(file)%num_fields
     field_num = files(file)%fields(j)
     if(output_fields(field_num)%time_ops) then
        time_ops = .true.
        exit
     endif
  enddo
! Loop through all fields with this file to output axes
  do j = 1, files(file)%num_fields
     field_num = files(file)%fields(j)
     input_field_num = output_fields(field_num)%input_field
     if (.not.input_fields(input_field_num)%register) then
        write (error_string,'(a,"/",a)')  &
             trim(input_fields(input_field_num)%module_name), &
             trim(input_fields(input_field_num)%field_name)
        if(mpp_pe() .eq. mpp_root_pe()) &
             call error_mesg ('diag_util opening_file', &
             'module/field_name '//trim(error_string)//' NOT registered', WARNING)  
        cycle
     endif
! Put the time axis in the axis field
     num_axes = output_fields(field_num)%num_axes
     axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
! make sure that axis_id are not -1
     do k = 1,num_axes
        if(axes(k)<0) then
           write(error_string,'(a)') output_fields(field_num)%output_name
           call error_mesg ('diag_util opening_file','output_name '//trim(error_string)// &
                ' has axis_id = -1', FATAL)
        endif        
     enddo
! check if aux is present in any axes
     if (.not. aux_present) then
        do k = 1,num_axes
           aux_name = get_axis_aux(axes(k))
           if (trim(aux_name) /= 'none') then
              aux_present = .true.
              exit
           endif
        end do
     endif

     axes(num_axes + 1) = files(file)%time_axis_id
     call write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 1), time_ops)
     if(time_ops) then
        axes(num_axes + 2) = files(file)%time_bounds_id
        call write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 2))     
     endif
  end do

! Looking for the first NON-static field in a file
  field_num1 = files(file)%fields(1)
  do j = 1, files(file)%num_fields
     field_num = files(file)%fields(j)
     if(output_fields(field_num)%time_ops) then
        field_num1 = field_num
        exit
     endif
  enddo
  do j = 1, files(file)%num_fields
     field_num = files(file)%fields(j)
     input_field_num = output_fields(field_num)%input_field
     if (.not.input_fields(input_field_num)%register) cycle
! Make sure that 1 file contains either time_average or instantaneous fields
! cannot have both time_average and instantaneous in 1 file
     if(.not. mix_snapshot_average_fields) then
        if((output_fields(field_num)%time_ops.neqv.output_fields(field_num1)%time_ops) .and. &
             .not.output_fields(field_num1)%static .and. .not.output_fields(field_num)%static) then
           if(mpp_pe() == mpp_root_pe()) call error_mesg ('diag_util opening_file','file '// &
                trim(files(file)%name)//' can NOT have BOTH time average AND instantaneous fields.'// &
                ' Create a new file or set mix_snapshot_average_fields=.true.' , FATAL)
        endif
     endif    
! check if any field has the same name as aux_name
     if(aux_present .and. .not. match_aux_name) then
        fieldname = output_fields(field_num)%output_name
        if(index(aux_name, trim(fieldname))>0) match_aux_name = .true.   
     endif

! Put the time axis in the axis field
     num_axes = output_fields(field_num)%num_axes
     axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
     if (.not. output_fields(field_num)%static) then
        num_axes=num_axes+1
        axes(num_axes) = files(file)%time_axis_id
     endif
     if(output_fields(field_num)%time_average) then
        avg = avg_name
     else if(output_fields(field_num)%time_max) then
        avg = avg_name
     else if(output_fields(field_num)%time_min) then
        avg = avg_name
     else
        avg = " "
     end if
     if(input_fields(input_field_num)%missing_value_present) then
        if(len_trim(input_fields(input_field_num)%interp_method)>0) then
           output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit, &
                output_fields(field_num)%output_name, axes(1:num_axes), &
                input_fields(input_field_num)%units, &
                input_fields(input_field_num)%long_name, &
                input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                input_fields(input_field_num)%missing_value, avg_name = avg,&
                time_method=output_fields(field_num)%time_method,&
                standard_name = input_fields(input_field_num)%standard_name, &
                interp_method = input_fields(input_field_num)%interp_method)
        else
           output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit, &
                output_fields(field_num)%output_name, axes(1:num_axes), &
                input_fields(input_field_num)%units, &
                input_fields(input_field_num)%long_name, &
                input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                input_fields(input_field_num)%missing_value, avg_name = avg,&
                time_method=output_fields(field_num)%time_method,&
                standard_name = input_fields(input_field_num)%standard_name)
        endif        
! NEED TO TAKE CARE OF TIME AVERAGING INFO TOO BOTH CASES
     else
        if(len_trim(input_fields(input_field_num)%interp_method) >0) then
           output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit, &
                output_fields(field_num)%output_name, axes(1:num_axes), &
                input_fields(input_field_num)%units, &
                input_fields(input_field_num)%long_name, &
                input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                avg_name = avg,&
                time_method=output_fields(field_num)%time_method, &
                standard_name = input_fields(input_field_num)%standard_name, &
                interp_method = input_fields(input_field_num)%interp_method)
        else
           output_fields(field_num)%f_type = write_field_meta_data(files(file)%file_unit, &
                output_fields(field_num)%output_name, axes(1:num_axes), &
                input_fields(input_field_num)%units, &
                input_fields(input_field_num)%long_name, &
                input_fields(input_field_num)%range, output_fields(field_num)%pack,&
                avg_name = avg,&
                time_method=output_fields(field_num)%time_method, &
                standard_name = input_fields(input_field_num)%standard_name)
        endif
     endif
  end do

! If any of the fields in the file are time averaged, need to output the axes
! Use double precision since time axis is double precision
   if(time_ops) then
      time_axis_id(1) = files(file)%time_axis_id
      files(file)%f_avg_start = write_field_meta_data(files(file)%file_unit, &
           avg_name // '_T1', time_axis_id, time_units, &
           "Start time for average period", pack=1)
      files(file)%f_avg_end = write_field_meta_data(files(file)%file_unit, &
           avg_name // '_T2', time_axis_id, time_units, &
           "End time for average period", pack=1)
      files(file)%f_avg_nitems = write_field_meta_data(files(file)%file_unit, &
           avg_name // '_DT', time_axis_id,     &
           trim(time_unit_list(files(file)%time_units)), &
           "Length of average period", pack=1)
  endif

  if (time_ops) then
      time_axis_id(1) = files(file)%time_axis_id
      time_bounds_id(1) = files(file)%time_bounds_id
      call get_diag_axis( time_axis_id(1), time_name, time_units, time_longname, &
           &cart_name, dir, edges, Domain, data)
      call get_diag_axis( time_bounds_id(1), timeb_name, timeb_units, timeb_longname, &
           &cart_name, dir, edges, Domain, data)     
      files(file)%f_bounds =  write_field_meta_data(files(file)%file_unit, &
           trim(time_name)//'_bounds', (/time_bounds_id,time_axis_id/),     &
           trim(time_unit_list(files(file)%time_units)), &
           trim(time_name)//' axis boundaries', pack=1)      
   end if
! Let lower levels know that all meta data has been sent
   call done_meta_data(files(file)%file_unit)
   if(aux_present .and. .not. match_aux_name) &
        call error_mesg ('diag_util opening_file','one axis has auxiliary but the corresponding field'// &
        ' is NOT found in file '//files(file)%name , WARNING)
 end subroutine opening_file
!- -------------------------------------------------------------------------------------------------

function get_time_string(filename, current_time)
! This function determines a string based on current time.
! This string is used as suffix in output file name

type(time_type), intent(in)    :: current_time
character(len=128), intent(in) :: filename
character(len=128)             :: get_time_string
integer                        :: yr1, mo1, dy1, hr1, mi1, sc1  ! get from current time
integer                        :: yr2, dy2, hr2, mi2            ! for computing next_level time unit
integer                        :: yr1_s, mo1_s, dy1_s, hr1_s, mi1_s, sc1_s ! actual values to write string
character(len=20)              :: yr, mo, dy, hr, mi, sc        ! string of current time (output)
integer                        :: abs_sec, abs_day              ! component of current_time
integer                        :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer                        :: julian_day, i, position, len, first_percent
character(len=10)              :: format
character(len=1)               :: width  ! width of the field in format write
character(len=128)             :: filetail

format = '("_",i*.*)'
call get_date(current_time, yr1, mo1, dy1, hr1, mi1, sc1)
len = len_trim(filename)
first_percent = index(filename, '%')
filetail = filename(first_percent:len)
! compute year string 
position = INDEX(filetail, 'yr')
if(position>0) then
   width = filetail(position-1:position-1)
   yr1_s = yr1
   format(7:9) = width//'.'//width
   write(yr, format) yr1_s   
   yr2 = 0
else  
   yr = ' '
   yr2 = yr1 - 1
endif
! compute month string 
position = INDEX(filetail, 'mo')
if(position>0) then   
   width = filetail(position-1:position-1)
   mo1_s = yr2*12 + mo1  
   format(7:9) = width//'.'//width
   write(mo, format) mo1_s
else
   mo = ' '
endif
! compute day string        
if(LEN_TRIM(mo) > 0) then       !  month present
   dy1_s = dy1 
   dy2 = dy1_s - 1
elseif(LEN_TRIM(yr) >0 )  then  ! no month, year present
! compute julian day
   if(mo1 == 1) then
      dy1_s = dy1
   else
      julian_day = 0
      do i = 1, mo1-1
         julian_day = julian_day + days_per_month(i)
      enddo
      if(leap_year(current_time) .and. mo1>2) julian_day = julian_day + 1
      julian_day = julian_day + dy1
      dy1_s = julian_day
   endif
   dy2 = dy1_s - 1
else                            ! no month, no year
   call get_time(current_time, abs_sec, abs_day)
   dy1_s = abs_day  
   dy2 = dy1_s 
endif
position = INDEX(filetail, 'dy')
if(position>0) then 
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(dy, format) dy1_s
else
   dy = ' '
endif
! compute hour string
if(LEN_TRIM(dy) > 0) then
   hr1_s = hr1
else
   hr1_s = dy2*24 + hr1
endif
hr2 = hr1_s
position = INDEX(filetail, 'hr')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(hr, format) hr1_s
else
   hr = ' '
endif
! compute minute string
if(LEN_TRIM(hr) > 0) then
   mi1_s = mi1
else
   mi1_s = hr2*60 + mi1
endif
mi2 = mi1_s
position = INDEX(filetail, 'mi')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(mi, format) mi1_s
else
   mi = ' '
endif
! compute second string
if(LEN_TRIM(mi) > 0) then
   sc1_s = sc1
else
   sc1_s = nint(mi2*SECONDS_PER_MINUTE) + sc1
endif
position = INDEX(filetail, 'sc')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(sc, format) sc1_s
else
   sc = ' '
endif
get_time_string = trim(yr)//trim(mo)//trim(dy)//trim(hr)//trim(mi)//trim(sc)
end function get_time_string
      
!--------------------------------------------------------------------------

function get_date_dif(t2, t1, units)

real                        :: get_date_dif
type(time_type), intent(in) :: t2, t1
integer, intent(in)         :: units
integer                     :: dif_seconds, dif_days
type(time_type)             :: dif_time

! Compute time axis label value
if(t2 < t1)   call error_mesg('get_date_dif', &
                't2 is less than t1', FATAL)

dif_time = t2 - t1

call get_time(dif_time, dif_seconds, dif_days)

if(units == DIAG_SECONDS) then
   get_date_dif = dif_seconds + SECONDS_PER_DAY * dif_days
else if(units == DIAG_MINUTES) then
   get_date_dif = 1440 * dif_days + dif_seconds / SECONDS_PER_MINUTE
else if(units == DIAG_HOURS) then
   get_date_dif = 24 * dif_days + dif_seconds / SECONDS_PER_HOUR
else if(units == DIAG_DAYS) then
   get_date_dif = dif_days + dif_seconds / SECONDS_PER_DAY
else if(units == DIAG_MONTHS) then
   call error_mesg('diag_data_out', 'months not supported as output units', FATAL)
else if(units == DIAG_YEARS) then
call error_mesg('diag_data_out', 'years not supported as output units', FATAL)
else
   call error_mesg('diag_data_out', 'illegal time units', FATAL)
end if

end function get_date_dif
!------------------------------------
subroutine diag_data_out(file, field, dat, time, final_call_in, static_write_in)

integer, intent(in)          :: file, field
real, intent(inout)          :: dat(:, :, :)
type(time_type), intent(in)  :: time
logical, optional, intent(in):: final_call_in, static_write_in
logical                      :: final_call, do_write, static_write
integer                      :: i, num
real                         :: dif, time_data(2, 1, 1), dt_time(1, 1, 1), start_dif, end_dif

do_write = .true.
final_call = .false.
if(present(final_call_in)) final_call = final_call_in
static_write = .false.
if(present(static_write_in)) static_write = static_write_in
dif = get_date_dif(time, base_time, files(file)%time_units)
! get file_unit, open new file and close curent file if necessary
If(.not. static_write .or. files(file)%file_unit<0) call check_and_open(file, time, do_write)
if(.not. do_write) return  ! no need to write data
call diag_field_out(files(file)%file_unit,output_fields(field)%f_type, dat, dif)
! record number of bytes written to this file
files(file)%bytes_written = files(file)%bytes_written + (size(dat,1)*size(dat,2)*size(dat,3))*(8/output_fields(field)%pack)
if(.not. output_fields(field)%written_once) output_fields(field)%written_once = .true.
! *** inserted this line because start_dif < 0 for static fields ***
if (.not. output_fields(field)%static) then 
   start_dif = get_date_dif(output_fields(field)%last_output, base_time,files(file)%time_units)
   if(.not. mix_snapshot_average_fields) then
      end_dif = get_date_dif(output_fields(field)%next_output, base_time, files(file)%time_units)
   else
      end_dif = dif
   endif
endif

! Need to write average axes out;
do i = 1, files(file)%num_fields
   num = files(file)%fields(i)
   if(output_fields(num)%time_ops .and. &
      input_fields(output_fields(num)%input_field)%register) then
      if(num == field) then
! Output the axes if this is first time-averaged field
         time_data(1, 1, 1) = start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_start, &
            time_data(1:1,:,:), dif)
         time_data(2, 1, 1) = end_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_end, &
            time_data(2:2,:,:), dif)
! Compute the length of the average
         dt_time(1, 1, 1) = end_dif - start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_nitems, &
            dt_time(1:1,:,:), dif)
!
! Include boundary variable for CF compliance
!
         call diag_field_out(files(file)%file_unit, files(file)%f_bounds, &
            time_data(1:2,:,:), dif)         
         exit
      endif
   end if
end do

! If write time is greater (equal for the last call) than last_flush for this file, flush it
if(final_call) then
   if(time >= files(file)%last_flush) then
      call diag_flush(files(file)%file_unit)
      files(file)%last_flush = time
   endif
else
   if(time > files(file)%last_flush) then
      call diag_flush(files(file)%file_unit)
      files(file)%last_flush = time
   endif
endif

end subroutine diag_data_out

!---------------------------------------------------------------------------------------------------
subroutine check_and_open(file, time, do_write )
! this routine checks if it is time to open a new file. If yes, it first closes the
! current file, open a new file and returns file_unit
! previous diag_manager_end is replaced by closing_file and output_setup by opening_file.

integer, intent(in)         :: file
type(time_type), intent(in) :: time
logical, intent(out)        :: do_write

if(time>=files(file)%start_time) then 
   if(files(file)%file_unit < 0)then ! need to open a new file
      call opening_file(file, time)
      do_write =.true.
   else
      do_write =.true.
      if(time > files(file)%close_time .and. time<files(file)%next_open ) then
         do_write = .false. ! file still open but receives NO MORE data
      elseif(time>files(file)%next_open) then ! need to close current file and open a new one 
         call write_static(file)  ! write all static fields and close this file
         call opening_file(file, time)        
         files(file)%start_time = files(file)%next_open
         files(file)%close_time = diag_time_inc(files(file)%start_time,files(file)%duration, &
            files(file)%duration_units)  
         files(file)%next_open = diag_time_inc(files(file)%next_open, files(file)%new_file_freq, &
              files(file)%new_file_freq_units)
         if (files(file)%close_time>files(file)%next_open) call error_mesg('check_and_open',&
              files(file)%name// &
              ' has close time GREATER than next_open time, check file duration and frequency',FATAL)
      endif ! no need to open new file, simply return file_unit
   endif  
else
   do_write = .false.
endif
end subroutine check_and_open
!---------------------------------------------------------------------------------------------------
subroutine write_static(file)
! Output all static fields in this file
integer, intent(in) :: file
integer             :: j, i, input_num

do j = 1, files(file)%num_fields
   i = files(file)%fields(j)
   input_num = output_fields(i)%input_field
! skip fields that were not registered
   if (.not.input_fields(input_num)%register) cycle
! only output static fields here
   if (.not.output_fields(i)%static) cycle
   call diag_data_out(file, i, output_fields(i)%buffer, files(file)%last_flush, .true., .true.)
end do
! Close up this file   
call mpp_close(files(file)%file_unit)
files(file)%file_unit = -1
end subroutine write_static
!---------------------------------------------------------------------------------------------------
subroutine check_duplicate_output_fields(err_msg)
! pair(output_name and output_file) should be unique in output_fields
character(len=*), intent(out), optional :: err_msg
integer            :: i, j, tmp_file
character(len=128) :: tmp_name
character(len=256) :: err_msg_local
! Do the checking when more than 1 output_fileds present

if(present(err_msg)) err_msg=''
if(num_output_fields <= 1) return 
err_msg_local = ''

i_loop: do i = 1, num_output_fields-1
  tmp_name = trim(output_fields(i)%output_name)
  tmp_file =  output_fields(i)%output_file
  do j = i+1, num_output_fields
    if((tmp_name == trim(output_fields(j)%output_name)).and.(tmp_file == output_fields(j)%output_file)) then
      err_msg_local = ' output_field "'//trim(tmp_name)//'" duplicated in file "'//trim(files(tmp_file)%name)//'"'
      exit i_loop
    endif
  enddo
enddo i_loop
if(err_msg_local /= '') then
  if(fms_error_handler(' ERROR in diag_table',err_msg_local,err_msg)) return
endif
end subroutine check_duplicate_output_fields

!---------------------------------------------------------------------------------------------------
end module diag_util_mod
