#include <fms_platform.h>

module fms_io_mod

!
!
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
! G.T. Nong
! </CONTACT>

! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </CONTACT>
!
! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </REVIEWER>

! <REVIEWER EMAIL="Bruce.Wyman@noaa.gov">
! B. Wyman
! </REVIEWER>

!<DESCRIPTION>
! This module is for writing and reading restart data in NetCDF format.
! fms_io_init must be called before the first write_data/read_data call
! For writing, fms_io_exit must be called after ALL write calls have
! been made. Typically, fms_io_init and fms_io_exit are placed in the
! main (driver) program while read_data and write_data can be called where needed.
! Presently, two combinations of threading and fileset are supported, users can choose
! one line of the following by setting namelist:
!
! With the introduction of netCDF restart files, there is a need for a global
! switch to turn on/off netCDF restart options in all of the modules that deal with
! restart files. Here two more namelist variables (logical type) are introduced to fms_io
!
! fms_netcdf_override
! fms_netcdf_restart
!
! because default values of both flags are .true., the default behavior of the entire model is 
! to use netCDF IO mode. To turn off netCDF restart, simply set fms_netcdf_restart to .false.
!
! Fei.Liu@noaa.gov
! 05222006
! Read distributed files in NetCDF is available. Details can be found in read_data_3d_new
! <PRE>
!threading_read='multi', threading_write='multi', fileset_write='multi' (default)
!threading_read='multi', threading_write='single', fileset_write='single'
! </PRE>
!</DESCRIPTION>
! <NAMELIST NAME="fms_io_nml">  
! <DATA NAME="threading_read" TYPE="character">
! threading_read can be 'single' or 'multi'
! </DATA>
! <DATA NAME="threading_write" TYPE="character">
! threading_write can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fileset_write" TYPE="character">
! fileset_write can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fms_netcdf_override" TYPE="logical">
!   .true. : fms_netcdf_restart overrides individual do_netcdf_restart value (default behavior)
!   .false.: individual module settings has a precedence over the global setting, therefore fms_netcdf_restart is ignored
! </DATA>
! <DATA NAME="fms_netcdf_restart" TYPE="logical">
!   .true. : all modules deal with restart files will operate under netCDF mode (default behavior)
!   .false.: all modules deal with restart files will operate under binary mode
!   This flag is effective only when fms_netcdf_override is .true. When fms_netcdf_override is .false., individual
!   module setting takes over.
! </DATA>
! <DATA NAME="time_stamped_restart" TYPE="logical">
!   .true. : time_stamp will be added to the restart file name as a prefix when 
!            optional argument time_stamp is passed into routine save_restart. 
!   .false.: time_stmp will not be added to the restart file name even though
!            time_stamp is passed into save_restart.
!    default is true.
! </DATA>
!</NAMELIST>
  
use mpp_io_mod,      only: mpp_open, mpp_close, mpp_io_init, mpp_io_exit, mpp_read, mpp_write
use mpp_io_mod,      only: mpp_write_meta, mpp_get_info, mpp_get_atts, mpp_get_fields
use mpp_io_mod,      only: mpp_get_axes, mpp_get_axis_data, mpp_get_att_char, mpp_get_att_name
use mpp_io_mod,      only: mpp_get_att_real_scalar
use mpp_io_mod,      only: fieldtype, axistype, atttype, default_field, default_axis, default_att
use mpp_io_mod,      only: MPP_NETCDF, MPP_ASCII, MPP_MULTI, MPP_SINGLE, MPP_OVERWR, MPP_RDONLY
use mpp_io_mod,      only: MPP_IEEE32, MPP_NATIVE, MPP_DELETE, MPP_APPEND, MPP_SEQUENTIAL, MPP_DIRECT
use mpp_io_mod,      only: MAX_FILE_SIZE
use mpp_domains_mod, only: domain2d, domain1d, NULL_DOMAIN1D, NULL_DOMAIN2D, operator( == ), CENTER
use mpp_domains_mod, only: mpp_get_domain_components, mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_domain_shift, mpp_get_global_domain, mpp_global_field, mpp_domain_is_root_pe
use mpp_domains_mod, only: mpp_get_ntile_count, mpp_get_current_ntile, mpp_get_tile_id, mpp_mosaic_defined
use mpp_mod,         only: mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_npes, stdlog, stdout
use mpp_mod,         only: mpp_broadcast, ALL_PES, mpp_chksum, mpp_get_current_pelist, mpp_npes, lowercase

use platform_mod, only: r8_kind

implicit none
private


integer, parameter, private :: max_split_file = 50
integer, parameter, private :: max_fields=400
integer, parameter, private :: max_axes=40
integer, parameter, private :: max_atts=20
integer, parameter, private :: max_domains = 10
integer, parameter, private :: MAX_TIME_LEVEL_REGISTER = 2
integer, parameter, private :: MAX_TIME_LEVEL_WRITE = 20
integer, parameter          :: max_axis_size=10000

type var_type
   private
   character(len=128)                     :: name
   character(len=128)                     :: longname
   character(len=128)                     :: units
   real, dimension(:,:,:,:), _ALLOCATABLE :: buffer _NULL 
   logical                                :: domain_present
   integer                                :: domain_idx
   logical                                :: is_dimvar
   type(fieldtype)                        :: field
   type(axistype)                         :: axis
   integer                                :: position
   integer                                :: ndim
   integer                                :: siz(4)      ! X/Y/Z/T extent of fields (data domain 
                                                         ! size for distributed writes;global size for reads)
   integer                                :: gsiz(4)     ! global X/Y/Z/T extent of fields 
   integer                                :: csiz(4)     ! actual data size in the file
   integer                                :: id_axes(3)  ! store index for x/y/z axistype.
   logical                                :: initialized ! indicate if the field is read or not in routine save_state.
   logical                                :: mandatory   ! indicate if the field is mandatory to be when restart.
end type var_type

type Ptr0Dr
   real,                   pointer :: p => NULL()
end type Ptr0Dr

type Ptr1Dr
   real, dimension(:),     pointer :: p => NULL()
end type Ptr1Dr

type Ptr2Dr
   real, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Dr

type Ptr3Dr
   real, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Dr

type Ptr0Di
   integer,                   pointer :: p => NULL()
end type Ptr0Di

type Ptr1Di
   integer, dimension(:),     pointer :: p => NULL()
end type Ptr1Di

type Ptr2Di
   integer, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Di

type Ptr3Di
   integer, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Di

type restart_file_type
   private
   integer                                  :: unit ! mpp_io unit for netcdf file
   character(len=128)                       :: name
   integer                                  :: nvar, natt, max_ntime
   logical                                  :: is_root_pe
   integer                                  :: tile_count
   type(var_type), dimension(:),   pointer  :: var  => NULL()
   type(Ptr0Dr),   dimension(:,:), pointer  :: p0dr => NULL()
   type(Ptr1Dr),   dimension(:,:), pointer  :: p1dr => NULL()
   type(Ptr2Dr),   dimension(:,:), pointer  :: p2dr => NULL()
   type(Ptr3Dr),   dimension(:,:), pointer  :: p3dr => NULL()
   type(Ptr0Di),   dimension(:,:), pointer  :: p0di => NULL()
   type(Ptr1Di),   dimension(:,:), pointer  :: p1di => NULL()
   type(Ptr2Di),   dimension(:,:), pointer  :: p2di => NULL()
   type(Ptr3Di),   dimension(:,:), pointer  :: p3di => NULL()
end type restart_file_type

interface read_data
   module procedure read_data_3d_new
   module procedure read_data_2d_new
   module procedure read_data_1d_new
   module procedure read_data_scalar_new
   module procedure read_data_i3d_new
   module procedure read_data_i2d_new
   module procedure read_data_i1d_new
   module procedure read_data_iscalar_new
   module procedure read_data_2d, read_ldata_2d, read_idata_2d
   module procedure read_data_3d, read_data_4d
#ifdef OVERLOAD_C8
   module procedure read_cdata_2d,read_cdata_3d,read_cdata_4d
#endif
   module procedure read_data_text
end interface

interface write_data
   module procedure write_data_3d_new
   module procedure write_data_2d_new
   module procedure write_data_1d_new
   module procedure write_data_scalar_new
   module procedure write_data_i3d_new
   module procedure write_data_i2d_new
   module procedure write_data_i1d_new
   module procedure write_data_iscalar_new
   module procedure write_data_2d, write_ldata_2d, write_idata_2d
   module procedure write_data_3d, write_data_4d
#ifdef OVERLOAD_C8
   module procedure write_cdata_2d,write_cdata_3d,write_cdata_4d
#endif
end interface

interface register_restart_field
   module procedure register_restart_field_r0d
   module procedure register_restart_field_r1d
   module procedure register_restart_field_r2d
   module procedure register_restart_field_r3d
   module procedure register_restart_field_i0d
   module procedure register_restart_field_i1d
   module procedure register_restart_field_i2d
   module procedure register_restart_field_i3d
   module procedure register_restart_field_r0d_2level
   module procedure register_restart_field_r1d_2level
   module procedure register_restart_field_r2d_2level
   module procedure register_restart_field_r3d_2level
   module procedure register_restart_field_i0d_2level
   module procedure register_restart_field_i1d_2level
   module procedure register_restart_field_i2d_2level
   module procedure register_restart_field_i3d_2level
end interface

interface reset_field_pointer
   module procedure reset_field_pointer_r0d
   module procedure reset_field_pointer_r1d
   module procedure reset_field_pointer_r2d
   module procedure reset_field_pointer_r3d
   module procedure reset_field_pointer_i0d
   module procedure reset_field_pointer_i1d
   module procedure reset_field_pointer_i2d
   module procedure reset_field_pointer_i3d
   module procedure reset_field_pointer_r0d_2level
   module procedure reset_field_pointer_r1d_2level
   module procedure reset_field_pointer_r2d_2level
   module procedure reset_field_pointer_r3d_2level
   module procedure reset_field_pointer_i0d_2level
   module procedure reset_field_pointer_i1d_2level
   module procedure reset_field_pointer_i2d_2level
   module procedure reset_field_pointer_i3d_2level
end interface

interface restore_state
   module procedure restore_state_all
   module procedure restore_state_one_field
end interface

interface query_initialized
   module procedure query_initialized_name
   module procedure query_initialized_r2d
end interface

interface get_global_att_value
  module procedure get_global_att_value_text
  module procedure get_global_att_value_real
end interface

integer :: num_files_r = 0 ! number of currently opened files for reading
integer :: num_files_w = 0 ! number of currently opened files for writing
integer :: num_domains = 0 ! number of domains in array_domain
integer :: num_registered_files ! mumber of files registered by calling register_restart_file


integer :: thread_r, thread_w, fset_w, form
logical :: module_is_initialized = .FALSE.

character(len=32) :: pelist_name
character(len=5)  :: pe_name
character(len=128):: error_msg  

  
!------ private data, pointer to current 2d domain ------
! entrained from fms_mod.  This will be deprecated in the future.
type(domain2D), pointer, private :: Current_domain =>NULL()

integer, private :: is,ie,js,je      ! compute domain
integer, private :: isd,ied,jsd,jed  ! data domain
integer, private :: isg,ieg,jsg,jeg  ! global domain
character(len=128),      dimension(:), allocatable         :: registered_file ! file names registered through register_restart_file 
type(restart_file_type), dimension(:), allocatable         :: files_read  ! store files that are read through read_data
type(restart_file_type), dimension(:), allocatable, target :: files_write ! store files that are written through write_data
type(domain2d), dimension(max_domains), save       :: array_domain
type(domain1d), dimension(max_domains), save       :: domain_x, domain_y
public  :: read_data, write_data, fms_io_init, fms_io_exit, field_size
public  :: open_namelist_file, open_restart_file, open_ieee32_file, close_file 
public  :: set_domain, nullify_domain, get_domain_decomp, return_domain
public  :: open_file, open_direct_file
public  :: get_restart_io_mode, get_tile_string, string
public  :: get_mosaic_tile_grid, get_mosaic_tile_file
public  :: get_global_att_value
public  :: file_exist, field_exist
public  :: register_restart_field, save_restart, restore_state
public  :: restart_file_type, query_initialized
public  :: reset_field_name, reset_field_pointer
private :: lookup_field_r, lookup_axis, unique_axes

!--- public interface ---
interface string
   module procedure string_from_integer
   module procedure string_from_real
end interface

!--- namelist interface
logical           :: fms_netcdf_override = .true.
logical           :: fms_netcdf_restart  = .true.
character(len=32) :: threading_read      = 'multi'
character(len=32) :: threading_write     = 'multi'
character(len=32) :: fileset_write       = 'multi'
character(len=32) :: format              = 'netcdf' 
logical           :: read_all_pe         = .TRUE.
character(len=64) :: iospec_ieee32       = '-N ieee_32'
integer           :: max_files_w         = 40
integer           :: max_files_r         = 40
logical           :: read_data_bug       = .false.
logical           :: time_stamp_restart  = .true.
  namelist /fms_io_nml/ fms_netcdf_override, fms_netcdf_restart, &
       threading_read, threading_write, &
       fileset_write, format, read_all_pe, iospec_ieee32,max_files_w,max_files_r, &
       read_data_bug, time_stamp_restart


character(len=128) :: version = '$Id: fms_io.F90,v 16.0 2008/07/30 22:45:32 fms Exp $'
character(len=128) :: tagname = '$Name: perth $'

contains

! <SUBROUTINE NAME="get_restart_io_mode">
! <DESCRIPTION>
! With the introduction of netCDF restart files, there is a need for a global
! switch to turn on/off netCDF restart options in all of the modules that deal with
! restart files. Here two more namelist variables (logical type) are introduced to fms_io
!
! fms_netcdf_override
! fms_netcdf_restart
!
! because default values of both flags are .true., the default behavior of the entire model is 
! to use netCDF IO mode. To turn off netCDF restart, simply set fms_netcdf_restart to .false.
! 
! </DESCRIPTION>
! <TEMPLATE>
!  call get_fms_io_mode(do_netcdf_restart)
! </TEMPLATE>
! <INOUT NAME="do_netcdf_restart" TYPE="logical">
!  This the input argument that contains the individual module setting of restart IO mode.
!  Upon return from this subroutine, this output argument contains the actual setting of restart IO mode
!  the calling module will be using
! </INOUT>
! </SUBROUTINE>
subroutine get_restart_io_mode(do_netcdf_restart)

  logical, intent(inout)  :: do_netcdf_restart

  if(fms_netcdf_override) do_netcdf_restart = fms_netcdf_restart
  
end subroutine get_restart_io_mode
!.....................................................................
! <SUBROUTINE NAME="fms_io_init">
!   <DESCRIPTION>
! Initialize fms_io module
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_init()
!   </TEMPLATE>
subroutine fms_io_init()
    
  integer                            :: i, unit, io_status, logunit
  integer, allocatable, dimension(:) :: pelist


  if (module_is_initialized) return
  call mpp_io_init()

  call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
  read(unit,fms_io_nml,iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>fms_io_init: Error reading input.nml')
  endif
  call mpp_close (unit)
  if (mpp_pe() == mpp_root_pe()) then
    logunit = stdlog() ; write(logunit, fms_io_nml)
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
  end if
! take namelist options if present

  select case (threading_read) 
  case ('multi')
     thread_r = MPP_MULTI
  case ('single')
     thread_r = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: threading_read should be multi/single but you chose'//trim(threading_read))
  end select
! take namelist options if present

  select case (fileset_write) 
  case ('multi')
     fset_w = MPP_MULTI
  case ('single')
     fset_w = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: fileset_write should be multi/single but you chose'//trim(fileset_write))
  end select

  select case (threading_write) 
  case ('multi')
     thread_w = MPP_MULTI
  case ('single')
     thread_w = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: threading_write should be multi/single but you chose'//trim(threading_write))
  end select

  select case(format)
  case ('netcdf')
     form=MPP_NETCDF
  case default
     call mpp_error(FATAL,'fms_io_init: only NetCDF format currently supported in fms_io')
  end select

! Initially allocate  files_write and files_read
  allocate(files_write(max_files_w),files_read(max_files_r))
  allocate(registered_file(max_files_w))

  allocate(pelist(mpp_npes()))        
  call mpp_get_current_pelist(pelist,pelist_name)
  write(pe_name,'(a,i4.4)' )'.', mpp_pe()    
  deallocate(pelist)

  do i = 1, max_domains
     array_domain(i) = NULL_DOMAIN2D
  enddo
  !---- initialize module domain2d pointer ----
  nullify (Current_domain)
  module_is_initialized = .TRUE.
  
end subroutine fms_io_init

! </SUBROUTINE>
! <SUBROUTINE NAME="fms_io_exit">
!   <DESCRIPTION>
! This routine is called after ALL fields have been written to temporary files
! The result NETCDF files are created here.
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_exit
!   </TEMPLATE>

subroutine fms_io_exit()
    integer                             :: num_x_axes, num_y_axes, num_z_axes
    integer                             :: unit
    real, dimension(max_axis_size)      :: axisdata
    real(r8_kind)                       :: tlev  
    real,           allocatable         :: global_data(:,:,:)
    integer,        dimension(max_axes) :: id_x_axes, siz_x_axes
    integer,        dimension(max_axes) :: id_y_axes, siz_y_axes
    integer,        dimension(max_axes) :: id_z_axes, siz_z_axes
    type(axistype), dimension(max_axes) :: x_axes, y_axes, z_axes
    type(axistype)                      :: t_axes
    type(var_type), pointer, save       :: cur_var=>NULL()
    integer                             :: i, j, k, kk
    character(len=256)                  :: filename
    character(len=10)                   :: axisname
    logical                             :: domain_present

    if( .NOT.module_is_initialized )return !make sure it's only called once per PE

    do i=1,max_axis_size
       axisdata(i) = i
    enddo

    ! each field has an associated domain type (may be undefined).
    ! each file only needs to write unique axes (i.e. if 2 fields share an identical axis, then only write the axis once)
    ! unique axes are defined by the global size and domain decomposition (i.e. can support identical axis sizes with
    ! different domain decomposition)

    do i = 1, num_files_w
       filename = files_write(i)%name

       !--- check if any field in this file present domain.
       domain_present = .false.
       do j = 1, files_write(i)%nvar
          if (files_write(i)%var(j)%domain_present) then
              domain_present = .true.
              exit
          end if
       end do

       !--- get the unique axes for all the fields.
       num_x_axes = unique_axes(files_write(i), 1, id_x_axes, siz_x_axes, domain_x)
       num_y_axes = unique_axes(files_write(i), 2, id_y_axes, siz_y_axes, domain_y)
       num_z_axes = unique_axes(files_write(i), 3, id_z_axes, siz_z_axes          )

       if( domain_present ) then
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=thread_w,&
               fileset=fset_w, is_root_pe=files_write(i)%is_root_pe)
       else  ! global data
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=MPP_SINGLE,&
               fileset=MPP_SINGLE, is_root_pe=files_write(i)%is_root_pe)
       end if

       do j = 1, num_x_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'xaxis_',j
          else
             write(axisname,'(a,i2)') 'xaxis_',j
          endif          
          if(id_x_axes(j) > 0) then
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),domain=domain_x(id_x_axes(j)),cartesian='X',    &
                  is_root_pe=files_write(i)%is_root_pe )
          else
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),cartesian='X', is_root_pe=files_write(i)%is_root_pe )
          endif             
       end do

       do j = 1, num_y_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'yaxis_',j
          else
             write(axisname,'(a,i2)') 'yaxis_',j
          endif          
          if(id_y_axes(j) > 0) then
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),domain=domain_y(id_y_axes(j)),cartesian='Y',    &
                  is_root_pe=files_write(i)%is_root_pe )
          else
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),cartesian='Y', is_root_pe=files_write(i)%is_root_pe )
          endif             
       end do

       do j = 1, num_z_axes
          if (j < 10) then
             write(axisname,'(a,i1)') 'zaxis_',j
          else
             write(axisname,'(a,i2)') 'zaxis_',j
          endif          
          call mpp_write_meta(unit,z_axes(j),axisname,'none',axisname, &
               data=axisdata(1:siz_z_axes(j)),cartesian='Z', is_root_pe=files_write(i)%is_root_pe )
       end do


       ! write time axis  (comment out if no time axis)
       call mpp_write_meta(unit,t_axes,&
            'Time','time level','Time',cartesian='T', is_root_pe=files_write(i)%is_root_pe)

       ! write metadata for fields
       do j = 1, files_write(i)%nvar
          cur_var => files_write(i)%var(j)
          call mpp_write_meta(unit,cur_var%field, (/x_axes(cur_var%id_axes(1)), &
               y_axes(cur_var%id_axes(2)), z_axes(cur_var%id_axes(3)), t_axes/), cur_var%name, &
               'none',cur_var%name,pack=1, is_root_pe=files_write(i)%is_root_pe)
       enddo

       ! write values for ndim of spatial axes
       do j = 1, num_x_axes
          call mpp_write(unit,x_axes(j), is_root_pe=files_write(i)%is_root_pe)
       enddo
       do j = 1, num_y_axes
          call mpp_write(unit,y_axes(j), is_root_pe=files_write(i)%is_root_pe)
       enddo
       do j = 1, num_z_axes
          call mpp_write(unit,z_axes(j), is_root_pe=files_write(i)%is_root_pe)
       enddo

       ! write data of each field
       do k = 1, files_write(i)%max_ntime
          do j = 1, files_write(i)%nvar
             cur_var => files_write(i)%var(j)
             tlev=k
             ! If some fields only have one time level, we do not need to write the second level, just keep
             ! the data missing.
             ! If some fields only have one time level, we just write out 0 to the other level
             if(k > cur_var%siz(4)) then
                cur_var%buffer(:,:,:,1) = 0.0
                kk = 1
             else
                kk = k
             end if
             if(cur_var%domain_present) then
                if(thread_w == MPP_MULTI) then
                   call mpp_write(unit, cur_var%field,array_domain(cur_var%domain_idx), cur_var%buffer(:,:,:,kk), tlev)
                else
                   allocate(global_data(cur_var%gsiz(1), cur_var%gsiz(2), cur_var%gsiz(3)) )
                   call mpp_global_field(array_domain(cur_var%domain_idx), cur_var%buffer(:,:,:,kk), global_data, &
                        position=cur_var%position,tile_count=files_write(i)%tile_count)

                   call mpp_write(unit,  cur_var%field, global_data, tlev, is_root_pe=files_write(i)%is_root_pe)
                   deallocate(global_data)
                end if
             else if (thread_w == MPP_MULTI .or. (files_write(i)%is_root_pe.and.thread_w == MPP_SINGLE)) then
                call mpp_write(unit, cur_var%field, cur_var%buffer(:,:,:,kk), tlev, is_root_pe=files_write(i)%is_root_pe)
             end if
          enddo ! end j loop
       enddo ! end k loop
       call mpp_close(unit)
    enddo ! end i loop

    !--- release the memory 

    do i = 1,  num_files_w
       do j = 1, files_write(i)%nvar
          deallocate(files_write(i)%var(j)%buffer)
       end do
    end do

  cur_var=>NULL()
  module_is_initialized = .false.
  num_files_w = 0
  num_files_r = 0    

end subroutine fms_io_exit
!.....................................................................
! </SUBROUTINE>

! <SUBROUTINE NAME="write_data">
    !<DESCRIPTION>
    ! This subroutine performs writing "fieldname" to file "filename". All values of "fieldname" 
    ! will be written to a temporary file. The final NETCDF file will be created only at a later step
    ! when the user calls fms_io_exit. Therefore, make sure that fms_io_exit is called after all
    ! fields have been written by this subroutine.
    !</DESCRIPTION>
!   <TEMPLATE>
! call write_data(filename, fieldname, data, domain)
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <IN NAME="data"  TYPE="real">
!   array containing data of fieldname
!   </IN>
!   <IN NAME="domain"  TYPE="domain, optional">
!   domain of fieldname
!   </IN>
!=================================================================================
subroutine write_data_i3d_new(filename, fieldname, data, domain,append_pelist_name, &
                              no_domain, position, tile_count, data_default, mosaicfile)

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: append_pelist_name, no_domain
  integer, intent(in), optional :: position, tile_count, data_default
  character(len=*), intent(in), optional :: mosaicfile
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)

  call write_data_3d_new(filename, fieldname, real(data), domain,append_pelist_name, &
                         no_domain, position, tile_count, data_default=default_data, mosaicfile=mosaicfile)
end subroutine write_data_i3d_new
!.....................................................................
subroutine write_data_i2d_new(filename, fieldname, data, domain,append_pelist_name, &
                              no_domain, position, tile_count, data_default, mosaicfile)

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: append_pelist_name, no_domain
  integer, intent(in), optional :: position, tile_count, data_default
  character(len=*), intent(in), optional :: mosaicfile
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_2d_new(filename, fieldname, real(data), domain,append_pelist_name, &
                         no_domain, position, tile_count, data_default=default_data, mosaicfile=mosaicfile)

end subroutine write_data_i2d_new
!.....................................................................
subroutine write_data_i1d_new(filename, fieldname, data, domain, append_pelist_name, &
                              no_domain, tile_count, data_default, mosaicfile)
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:), intent(in) :: data
  logical, intent(in), optional :: append_pelist_name, no_domain
  integer, intent(in), optional :: tile_count, data_default
  character(len=*), intent(in), optional :: mosaicfile
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_1d_new(filename, fieldname, real(data), domain, append_pelist_name, &
                         no_domain, tile_count, data_default=default_data, mosaicfile=mosaicfile)
end subroutine write_data_i1d_new
!.....................................................................
subroutine write_data_iscalar_new(filename, fieldname, data, domain, append_pelist_name, &
                                  no_domain, tile_count, data_default, mosaicfile)
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, intent(in) :: data
  logical, intent(in), optional :: append_pelist_name, no_domain
  integer, intent(in), optional :: tile_count, data_default
  character(len=*), intent(in), optional :: mosaicfile
  real :: default_data

  default_data = 0
  if(present(data_default)) default_data = real(data_default)
  call write_data_scalar_new(filename, fieldname, real(data), domain, append_pelist_name, &
                             no_domain, tile_count, data_default=default_data, mosaicfile=mosaicfile)

end subroutine write_data_iscalar_new
!.....................................................................
subroutine write_data_3d_new(filename, fieldname, data, domain, append_pelist_name, no_domain, &
                             position, tile_count, data_default, mosaicfile)

  character(len=*),         intent(in)         :: filename, fieldname 
  real, dimension(:,:,:),   intent(in)         :: data
  type(domain2d), optional, intent(in), target :: domain
  real,           optional, intent(in)         :: data_default
  logical,        optional, intent(in)         :: append_pelist_name, no_domain   
  integer,        optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)       :: mosaicfile

  !--- local variables
  real,               allocatable :: tmp_buffer(:,:,:,:)
  integer                         :: index_field ! position of the fieldname in the list of fields
  integer                         :: index_file  ! position of the filename in the list of files_write
  logical                         :: append_pelist, is_no_domain
  character(len=256)              :: fname, filename2
  real                            :: default_data
  integer                         :: length, i, domain_idx
  integer                         :: ishift, jshift
  integer                         :: gxsize, gysize
  integer                         :: cxsize, cysize
  integer                         :: dxsize, dysize
  type(domain2d), pointer, save   :: d_ptr   =>NULL()
  type(var_type), pointer, save   :: cur_var =>NULL()
  type(restart_file_type), pointer, save :: cur_file =>NULL()

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_3d_new): need to call fms_io_init')  

  append_pelist = .false.
  if(present(append_pelist_name)) append_pelist = append_pelist_name

  if(PRESENT(data_default))then
     default_data=data_default
  else
     default_data=0.
  endif

  if(present(tile_count) .AND. .not. present(domain)) call mpp_error(FATAL, &
         'fms_io write_data: when tile_count is present, domain must be present')

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     is_no_domain = no_domain
  end if

  if(is_no_domain) then
!     comments the following to allow passing domain when writing 1-D or scalar variables. 
!     if(PRESENT(domain)) &
!       call mpp_error(FATAL, 'fms_io(write_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')  
  else if(PRESENT(domain))then
     d_ptr => domain
  else if (ASSOCIATED(Current_domain) .and. .NOT. is_no_domain ) then
     d_ptr => Current_domain
  endif

  !--- remove .nc from file name
  length = len_trim(filename)
  if(filename(length-2:length) == '.nc') then
     filename2 = filename(1:length-3)
  else
     filename2 = filename(1:length)
  end if

  if(present(append_pelist_name)) then
     if(append_pelist_name) filename2 = trim(filename2)//'.'//trim(pelist_name)
  end if

  !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count,
  !JWD:  I just don't see how the physics can remain "tile neutral"
  !z1l:  one solution is add one more public interface called set_tile_count
  call get_mosaic_tile_file(filename2, fname, is_no_domain, domain, tile_count, mosaicfile)

  ! Check if filename has been open  or not
  index_file = -1
  do i=1,num_files_w
     if (trim(files_write(i)%name) == trim(fname)) then
        index_file = i
        cur_file => files_write(index_file)
        exit
     endif
  enddo

  if (index_file < 0) then 
     if(num_files_w == max_files_w) &  ! need to have bigger max_files_w
          call mpp_error(FATAL,'fms_io(write_data_3d_new): max_files_w exceeded, increase it via fms_io_nml')    
     ! record the file name in array files_write
     num_files_w=num_files_w + 1
     index_file = num_files_w
     cur_file => files_write(index_file)
     cur_file%name = trim(fname)         
     cur_file%tile_count=1
     if(present(tile_count)) cur_file%tile_count = tile_count
     if(ASSOCIATED(d_ptr))then
        cur_file%is_root_pe = mpp_domain_is_root_pe(d_ptr)
     else
        cur_file%is_root_pe = mpp_pe() == mpp_root_pe()
     endif
     cur_file%max_ntime = 1
     !-- allocate memory
     allocate(cur_file%var(max_fields) )
     do i = 1, max_fields
        cur_file%var(i)%name           = 'none'
        cur_file%var(i)%domain_present = .false.
        cur_file%var(i)%domain_idx     = -1
        cur_file%var(i)%is_dimvar      = .false.
        cur_file%var(i)%position       = CENTER
        cur_file%var(i)%siz(:)         = 0
        cur_file%var(i)%gsiz(:)        = 0
        cur_file%var(i)%id_axes(:)     = -1
     end do
  endif

  ! check if the field is new or not and get position and dimension of the field
  index_field = -1
  do i = 1, cur_file%nvar
     if(trim(cur_file%var(i)%name) == trim(fieldname)) then
        index_field = i
        exit
     end if 
  end do

  if(index_field > 0) then
     cur_var   => cur_file%var(index_field)
     cur_var%siz(4) =  cur_var%siz(4) + 1
     if(cur_file%max_ntime < cur_var%siz(4) ) cur_file%max_ntime = cur_var%siz(4)
     ! the time level should be no larger than MAX_TIME_LEVEL_WRITE ( =20) for write_data.
     if( cur_var%siz(4) > MAX_TIME_LEVEL_WRITE ) call mpp_error(FATAL, 'fms_io(write_data_3d_new): ' // &
          'the time level of field '//trim(cur_var%name)//' in file '//trim(cur_file%name)// &
          ' is greater than MAX_TIME_LEVEL_WRITE(=20), increase MAX_TIME_LEVEL_WRITE or check your code')
  else 
     cur_file%nvar = cur_file%nvar +1
     if(cur_file%nvar>max_fields) then
        write(error_msg,'(I3,"/",I3)') cur_file%nvar, max_fields 
        call  mpp_error(FATAL,'fms_io(write_data_3d_new): max_fields exceeded, needs increasing, nvar/max_fields=' &
             //trim(error_msg))
     endif
     index_field =  cur_file%nvar
     cur_var   => cur_file%var(index_field)
     cur_var%siz(1)  = size(data,1)
     cur_var%siz(2)  = size(data,2)
     cur_var%siz(3)  = size(data,3)
     cur_var%siz(4)  = 1
     cur_var%gsiz(3) = cur_var%siz(3)
     cur_var%name = fieldname
     cur_var%ndim = 3
     if(present(position)) cur_var%position = position
     if(ASSOCIATED(d_ptr)) then
        cur_var%domain_present = .true.
        domain_idx = lookup_domain(d_ptr)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io(write_data_3d_new), 1: max_domains exceeded,' &
                //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = d_ptr
           call mpp_get_domain_components(array_domain(domain_idx), domain_x(domain_idx), domain_y(domain_idx), &
                tile_count=tile_count)
        endif
        cur_var%domain_idx = domain_idx
        call mpp_get_domain_shift ( array_domain(domain_idx), ishift, jshift, position)
        call mpp_get_global_domain(array_domain(domain_idx), xsize=gxsize,ysize=gysize,tile_count=tile_count)
        call mpp_get_compute_domain(array_domain(domain_idx), xsize = cxsize, ysize = cysize, tile_count=tile_count)
        call mpp_get_data_domain   (array_domain(domain_idx), xsize = dxsize, ysize = dysize, tile_count=tile_count)
        if (ishift .NE. 0) then
           cxsize = cxsize+ishift; dxsize = dxsize+ishift; gxsize = gxsize + ishift
        end if
        if (jshift .NE. 0) then
           cysize = cysize+jshift; dysize = dysize+jshift; gysize = gysize + jshift
        endif
        if( (cur_var%siz(1) .NE. cxsize .AND. cur_var%siz(1) .NE. dxsize ) .OR. &
            (cur_var%siz(2) .NE. cysize .AND. cur_var%siz(2) .NE. dysize ) ) then
            call mpp_error(FATAL, 'fms_io(write_data_3d_new): data should be on either computer domain '//&
              'or data domain when domain is present for field '//trim(fieldname)//' of file '//trim(filename) )
        end if
        cur_var%gsiz(1)   = gxsize
        cur_var%gsiz(2)   = gysize
     else
        cur_var%gsiz(1) = size(data,1)
        cur_var%gsiz(2) = size(data,2)
        cur_var%domain_present=.false.
     endif
  end if

  ! copy the data to the buffer 
  ! if the time level is greater than the size(cur_var%buffer,4), 
  ! need to increase the buffer size 

  if(cur_var%siz(4) == 1) then
     allocate(cur_var%buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), cur_var%siz(4)) )
  else
     allocate(tmp_buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), size(cur_var%buffer,4)) )
     tmp_buffer = cur_var%buffer
     deallocate(cur_var%buffer)
     allocate(cur_var%buffer(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3), cur_var%siz(4)) )
     cur_var%buffer(:,:,:,1:size(tmp_buffer,4)) = tmp_buffer 
     deallocate(tmp_buffer)
  endif

  cur_var%buffer(:,:,:,cur_var%siz(4)) = data ! copy current data to buffer for future write out 

  d_ptr =>NULL()
  cur_var =>NULL()
  cur_file =>NULL()

end subroutine write_data_3d_new
! </SUBROUTINE>  

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r0d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                                      position, tile_count, data_default, longname, units)
  type(restart_file_type),    intent(inout)      :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,                       intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: mandatory
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r0d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r0d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 1/), index_field, domain,      &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 0
  register_restart_field_r0d = index_field  

  return    

end function register_restart_field_r0d

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r1d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real, dimension(:),         intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r1d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r1d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), 1, 1, 1/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units )

  fileObj%p1dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 1
  register_restart_field_r1d = index_field  

  return    
  
end function register_restart_field_r1d

!-------------------------------------------------------------------------------
!
!   The routine will register a 2-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r2d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:),   intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r2d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r2d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), 1, 1/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 2
  register_restart_field_r2d = index_field    

  return    

end function register_restart_field_r2d


!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r3d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:,:), intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r3d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_r3d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), size(data,3), 1/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3dr(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 3
  register_restart_field_r3d = index_field   

  return    

end function register_restart_field_r3d

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar integer restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i0d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,                    intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i0d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i0d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 1/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 0
  register_restart_field_i0d = index_field
  
  return
   
end function register_restart_field_i0d

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D integer restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i1d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer, dimension(:),      intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i1d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i1d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), 1, 1, 1/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 1
  register_restart_field_i1d = index_field
  
  return
  
end function register_restart_field_i1d


!-------------------------------------------------------------------------------
!
!   The routine will register a 2-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i2d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:),   intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i2d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i2d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), 1, 1/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 2
  register_restart_field_i2d = index_field
  
  return 

end function register_restart_field_i2d

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with one time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i3d(fileObj, filename, fieldname, data, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:,:), intent(in), target :: data 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i3d

  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(register_restart_field_i3d): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data,1), size(data,2), size(data,3), 1/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3di(fileObj%var(index_field)%siz(4), index_field)%p => data
  fileObj%var(index_field)%ndim = 3
  register_restart_field_i3d = index_field
  
  return  

end function register_restart_field_i3d

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r0d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,                       intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r0d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r0d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 2/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0dr(1, index_field)%p => data1
  fileObj%p0dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 0
  register_restart_field_r0d_2level = index_field
  
  return  

end function register_restart_field_r0d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r1d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:),     intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r1d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r1d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), 1, 1, 2/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1dr(1, index_field)%p => data1
  fileObj%p1dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 1
  register_restart_field_r1d_2level = index_field
  
  return    

end function register_restart_field_r1d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r2d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:),   intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r2d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r2d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), 1, 2/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2dr(1, index_field)%p => data1
  fileObj%p2dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 2
  register_restart_field_r2d_2level = index_field
  
  return    

end function register_restart_field_r2d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D real restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_r3d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  real,     dimension(:,:,:), intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_r3d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_r3d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), size(data1,3), 2/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3dr(1, index_field)%p => data1
  fileObj%p3dr(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 3
  register_restart_field_r3d_2level = index_field
  
  return    

end function register_restart_field_r3d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a scalar integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i0d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,                    intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i0d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i0d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/1, 1, 1, 2/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p0di(1, index_field)%p => data1
  fileObj%p0di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 0
  register_restart_field_i0d_2level = index_field
  
  return   
 
end function register_restart_field_i0d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 1-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i1d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:),     intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i1d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i1d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), 1, 1, 2/), index_field, domain, &
                       instance_name, mandatory, no_domain=.true., position=position, tile_count=tile_count, &
                       data_default=data_default, longname=longname, units=units)
  fileObj%p1di(1, index_field)%p => data1
  fileObj%p1di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 1
  register_restart_field_i1d_2level = index_field
  
  return  
  
end function register_restart_field_i1d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i2d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:),   intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i2d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i2d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), 1, 2/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p2di(1, index_field)%p => data1
  fileObj%p2di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 2
  register_restart_field_i2d_2level = index_field
  
  return    

end function register_restart_field_i2d_2level

!-------------------------------------------------------------------------------
!
!   The routine will register a 3-D integer restart file field with two time level
!
!-------------------------------------------------------------------------------
function register_restart_field_i3d_2level(fileObj, filename, fieldname, data1, data2, domain, instance_name, mandatory, &
                             no_domain, position, tile_count, data_default, longname, units)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),           intent(in)         :: filename, fieldname
  integer,  dimension(:,:,:), intent(in), target :: data1, data2 
  type(domain2d),   optional, intent(in), target :: domain
  real,             optional, intent(in)         :: data_default
  logical,          optional, intent(in)         :: no_domain   
  integer,          optional, intent(in)         :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  integer                                        :: index_field
  integer                                        :: register_restart_field_i3d_2level

  if(.not.module_is_initialized) call mpp_error(FATAL, &
      'fms_io(register_restart_field_i3d_2level): need to call fms_io_init')  
  call setup_one_field(fileObj, filename, fieldname, (/size(data1,1), size(data1,2), size(data1,3), 2/), &
                       index_field, domain, instance_name, mandatory, no_domain, &
                       position, tile_count, data_default, longname, units)
  fileObj%p3di(1, index_field)%p => data1
  fileObj%p3di(2, index_field)%p => data2
  fileObj%var(index_field)%ndim = 3
  register_restart_field_i3d_2level = index_field
  
  return    

end function register_restart_field_i3d_2level

!-------------------------------------------------------------------------------
!
!  saves all registered variables to restart files. Those variables are set 
!  through register_restart_field
!
!-------------------------------------------------------------------------------

subroutine save_restart(fileObj, time_stamp, directory )
  type(restart_file_type), intent(inout) :: fileObj
  character(len=*), intent(in), optional :: directory
  character(len=*), intent(in), optional :: time_stamp
  ! Arguments: 
  !  (in)      directory  - The directory where the restart file goes.
  !  (in)      time_stamp - character format of the time of this restart file.
  character(len=256) :: dir
  character(len=256) :: restartpath          ! The restart file path (dir/file).
  character(len=80)  :: restartname          ! The restart file name (no dir).
  character(len=8)   :: suffix               ! A suffix (like _2) that is appended to the name of files after the first.
  integer            :: var_sz, size_in_file ! The size in bytes of each variable and of the variables already in a file.
  integer            :: start_var, next_var  ! The starting variables of the current and next files.
  integer            :: unit                 ! The mpp unit of the open file.
  real, dimension(max_axis_size)      :: axisdata
  integer,        dimension(max_axes) :: id_x_axes, siz_x_axes
  integer,        dimension(max_axes) :: id_y_axes, siz_y_axes
  integer,        dimension(max_axes) :: id_z_axes, siz_z_axes
  type(axistype), dimension(max_axes) :: x_axes, y_axes, z_axes
  type(axistype)                      :: t_axes            
  integer                             :: num_var_axes
  type(axistype), dimension(4)        :: var_axes
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: num_x_axes, num_y_axes, num_z_axes
  integer                             :: naxes_x, naxes_y, naxes_z
  integer                             :: nfiles, i, j, k, l, siz
  logical                             :: domain_present
  real(r8_kind)                       :: tlev  
  character(len=10)                   :: axisname  

  real, allocatable, dimension(:,:,:) :: r3d, global_r3d
  real, allocatable, dimension(:,:)   :: r2d, global_r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(save_restart): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  do i=1,max_axis_size
     axisdata(i) = i
  enddo

  dir = "RESTART"
  if(present(directory)) dir = directory

  restartname = fileObj%name
  nfiles = 0
  if(time_stamp_restart) then
     if (PRESENT(time_stamp)) then
        restartname = trim(time_stamp)//"."//trim(restartname)
     endif
  end if
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(restartname)
  else
     restartpath = trim(restartname)
  end if
  !--- check if any field in this file present domain.
  domain_present = .false.
  do j = 1, fileObj%nvar
     if (fileObj%var(j)%domain_present) then
        domain_present = .true.
        exit
     end if
  end do
  num_x_axes = unique_axes(fileObj, 1, id_x_axes, siz_x_axes, domain_x)
  num_y_axes = unique_axes(fileObj, 2, id_y_axes, siz_y_axes, domain_y)
  num_z_axes = unique_axes(fileObj, 3, id_z_axes, siz_z_axes          )
  next_var = 1
  size_in_file = 0
  do j = 1, num_x_axes
     size_in_file = size_in_file + siz_x_axes(j)
  end do
  do j = 1, num_y_axes
     size_in_file = size_in_file + siz_y_axes(j)
  end do
  do j = 1, num_z_axes
     size_in_file = size_in_file + siz_z_axes(j)
  end do
  size_in_file = 8*(size_in_file*2+1000)

  do while (next_var <= fileObj%nvar )
     start_var = next_var

     do j=start_var,fileObj%nvar
        cur_var => fileObj%var(j)
        var_sz = 8*cur_var%csiz(1)*cur_var%csiz(2)*cur_var%csiz(3)
        if ((j==start_var) .OR. (size_in_file < MAX_FILE_SIZE-var_sz)) then
           size_in_file = size_in_file + var_sz
        else 
           exit
        endif
     enddo
     next_var = j
     ! For distribute write, normally will not over the limit. 
     if( nfiles > 0 ) then
        if(fset_w == MPP_MULTI .AND. domain_present) call mpp_error(FATAL, "fms_io_mod(save_restart): "// &
             "For distribute write(fileset_write='multi'), the file size should not be very large and need to be split")
        if (nfiles < 10) then
           write(suffix,'("_",I1)') nfiles
        else if(nfiles < 100) then
           write(suffix,'("_",I2)') nfiles
        else
           call mpp_error(FATAL, "fms_io(save_restart): num_files should be less than 100")
        endif
        !--- remove .nc from restartpath and attach suffix.
        siz = len_trim(restartpath)
        if(restartpath(siz-2:siz) == ".nc") then
           restartpath = restartpath(1:siz-3)//trim(suffix)
        else      
           restartpath = trim(restartpath) // trim(suffix)
        end if
     end if
     if( domain_present ) then
        call mpp_open(unit,trim(restartpath),action=MPP_OVERWR,form=form,threading=thread_w,&
             fileset=fset_w, is_root_pe=fileObj%is_root_pe)
     else  ! global data
        call mpp_open(unit,trim(restartpath),action=MPP_OVERWR,form=form,threading=MPP_SINGLE,&
             fileset=MPP_SINGLE, is_root_pe=fileObj%is_root_pe)
     end if

     ! write_out x_axes
     naxes_x = 0
     do j = 1, num_x_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(1) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_x = naxes_x + 1
        if (naxes_x < 10) then
           write(axisname,'(a,i1)') 'xaxis_',naxes_x
        else
           write(axisname,'(a,i2)') 'xaxis_',naxes_x
        endif
        if(id_x_axes(j) > 0) then
           call mpp_write_meta(unit,x_axes(naxes_x),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),domain=domain_x(id_x_axes(j)),cartesian='X',    &
                is_root_pe=fileObj%is_root_pe )
        else
           call mpp_write_meta(unit,x_axes(naxes_x),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),cartesian='X', is_root_pe=fileObj%is_root_pe )
        endif
     end do

     ! write out y_axes
     naxes_y = 0
     do j = 1, num_y_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(2) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_y = naxes_y + 1
        if (naxes_y < 10) then
           write(axisname,'(a,i1)') 'yaxis_',naxes_y
        else
           write(axisname,'(a,i2)') 'yaxis_',naxes_y
        endif
        if(id_y_axes(j) > 0) then
           call mpp_write_meta(unit,y_axes(naxes_y),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),domain=domain_y(id_y_axes(j)),cartesian='Y',    &
                is_root_pe=fileObj%is_root_pe )
        else
           call mpp_write_meta(unit,y_axes(naxes_y),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),cartesian='Y', is_root_pe=fileObj%is_root_pe )
        endif
     end do

     ! write out z_axes
     naxes_z = 0
     do j = 1, num_z_axes
        ! make sure this axis is used by some variable 
        do l=start_var,next_var-1
           if( fileObj%var(l)%id_axes(3) == j ) exit
        end do
        if(l == next_var) cycle  
        naxes_z = naxes_z + 1
        if (naxes_z < 10) then
           write(axisname,'(a,i1)') 'zaxis_',naxes_z
        else
           write(axisname,'(a,i2)') 'zaxis_',naxes_z
        endif
        call mpp_write_meta(unit,z_axes(naxes_z),axisname,'none',axisname, &
             data=axisdata(1:siz_z_axes(j)),cartesian='Z', is_root_pe=fileObj%is_root_pe )
     end do

     ! write out time axis  
     call mpp_write_meta(unit,t_axes,&
          'Time','time level','Time',cartesian='T', is_root_pe=fileObj%is_root_pe)
     ! write metadata for fields
     do j = start_var,next_var-1
        cur_var => fileObj%var(j)
        if(cur_var%siz(4) > 1 .AND. cur_var%siz(4) .NE. fileObj%max_ntime ) call mpp_error(FATAL, &
         "fms_io(save_restart): "//trim(cur_var%name)//" in file "//trim(fileObj%name)// &
         " has more than one time level, but number of time level is not equal to max_ntime")

        if(cur_var%ndim == 0) then
           num_var_axes = 1
           var_axes(1) = t_axes
        else if(cur_var%ndim == 1) then
           num_var_axes = 1
           var_axes(1) = x_axes(cur_var%id_axes(1))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 2
              var_axes(2) = t_axes
           end if
        else if(cur_var%ndim == 2) then
           num_var_axes = 2
           var_axes(1) = x_axes(cur_var%id_axes(1))
           var_axes(2) = y_axes(cur_var%id_axes(2))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 3
              var_axes(3) = t_axes
           end if
        else if(cur_var%ndim == 3) then
           num_var_axes = 3
           var_axes(1) = x_axes(cur_var%id_axes(1))
           var_axes(2) = y_axes(cur_var%id_axes(2))
           var_axes(3) = z_axes(cur_var%id_axes(3))
           if(cur_var%siz(4) == fileObj%max_ntime) then
              num_var_axes = 4
              var_axes(4) = t_axes
           end if
        end if
        call mpp_write_meta(unit,cur_var%field, var_axes(1:num_var_axes), cur_var%name, &
                 cur_var%units,cur_var%longname,pack=1, is_root_pe=fileObj%is_root_pe)
     enddo

     ! write values for ndim of spatial axes
     do j = 1, naxes_x
        call mpp_write(unit,x_axes(j), is_root_pe=fileObj%is_root_pe)
     enddo
     do j = 1, naxes_y
        call mpp_write(unit,y_axes(j), is_root_pe=fileObj%is_root_pe)
     enddo
     do j = 1, naxes_z
        call mpp_write(unit,z_axes(j), is_root_pe=fileObj%is_root_pe)
     enddo

     ! write data of each field
     do k = 1, fileObj%max_ntime
        do j=start_var,next_var-1
           cur_var => fileObj%var(j)
           tlev=k
           ! If some fields only have one time level, we do not need to write the second level, just keep
           ! the data missing.
           if(k <= cur_var%siz(4)) then
              if(cur_var%domain_present) then  ! one 2-D or 3-D case possible present domain
                 if(thread_w == MPP_MULTI) then
                    if( Associated(fileObj%p2dr(k,j)%p) ) then
                       call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), fileObj%p2dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), fileObj%p3dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = fileObj%p2di(k,j)%p
                       call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), r2d, tlev)
                       deallocate(r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = fileObj%p3di(k,j)%p
                       call mpp_write(unit, cur_var%field, array_domain(cur_var%domain_idx), r3d, tlev)
                       deallocate(r3d)
                    else
                       call mpp_error(FATAL, "fms_io(save_restart): domain is present and thread_w  "// &
                            "is MPP_MULTI, field "//trim(cur_var%name)//" of file "//trim(fileObj%name)// &
                            ", but none of p2dr, p3dr, p2di and p3di is associated") 
                    end if
                 else 
                    if( Associated(fileObj%p2dr(k,j)%p) ) then
                       allocate(global_r2d(cur_var%gsiz(1), cur_var%gsiz(2)) )
                       call mpp_global_field(array_domain(cur_var%domain_idx), fileObj%p2dr(k,j)%p, global_r2d, &
                            position=cur_var%position,tile_count=fileObj%tile_count)
                       call mpp_write(unit,  cur_var%field, global_r2d, tlev, is_root_pe=fileObj%is_root_pe)
                       deallocate(global_r2d)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       allocate(global_r3d(cur_var%gsiz(1), cur_var%gsiz(2), cur_var%gsiz(3)) )
                       call mpp_global_field(array_domain(cur_var%domain_idx), fileObj%p3dr(k,j)%p, global_r3d, &
                            position=cur_var%position,tile_count=fileObj%tile_count)
                       call mpp_write(unit,  cur_var%field, global_r3d, tlev, is_root_pe=fileObj%is_root_pe)
                       deallocate(global_r3d)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(global_r2d(cur_var%gsiz(1), cur_var%gsiz(2)) )
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = fileObj%p2di(k,j)%p
                       call mpp_global_field(array_domain(cur_var%domain_idx), r2d, global_r2d, &
                            position=cur_var%position,tile_count=fileObj%tile_count)
                       call mpp_write(unit,  cur_var%field, global_r2d, tlev, is_root_pe=fileObj%is_root_pe)
                       deallocate(global_r2d, r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(global_r3d(cur_var%gsiz(1), cur_var%gsiz(2), cur_var%gsiz(3)) )
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = fileObj%p3di(k,j)%p
                       call mpp_global_field(array_domain(cur_var%domain_idx), r3d, global_r3d, &
                            position=cur_var%position,tile_count=fileObj%tile_count)
                       call mpp_write(unit,  cur_var%field, global_r3d, tlev, is_root_pe=fileObj%is_root_pe)
                       deallocate(global_r3d, r3d)
                    else
                       call mpp_error(FATAL, "fms_io(save_restart): domain is present and thread_w  "// &
                            "is MPP_SINGLE, field "//trim(cur_var%name)//" of file "//trim(fileObj%name)// &
                            ", but none of p2dr, p3dr, p2di and p3di is associated") 
                    end if
                 end if
              else if (thread_w == MPP_MULTI .or. (fileObj%is_root_pe.and.thread_w == MPP_SINGLE)) then     
                 if ( Associated(fileObj%p0dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p0dr(k,j)%p, tlev, is_root_pe=fileObj%is_root_pe)
                 else if ( Associated(fileObj%p1dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p1dr(k,j)%p, tlev, is_root_pe=fileObj%is_root_pe)
                 else if ( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p2dr(k,j)%p, tlev, is_root_pe=fileObj%is_root_pe)
                 else if ( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_write(unit, cur_var%field, fileObj%p3dr(k,j)%p, tlev, is_root_pe=fileObj%is_root_pe)
                 else if ( Associated(fileObj%p0di(k,j)%p) ) then
                    r0d =  fileObj%p0di(k,j)%p
                    call mpp_write(unit, cur_var%field, r0d,                  tlev, is_root_pe=fileObj%is_root_pe)
                 else if ( Associated(fileObj%p1di(k,j)%p) ) then
                    allocate(r1d(cur_var%siz(1)) )
                    r1d = fileObj%p1di(k,j)%p
                    call mpp_write(unit, cur_var%field, r1d,                  tlev, is_root_pe=fileObj%is_root_pe)
                    deallocate(r1d)
                 else if ( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = fileObj%p2di(k,j)%p
                    call mpp_write(unit, cur_var%field, r2d,                  tlev, is_root_pe=fileObj%is_root_pe)
                    deallocate(global_r2d, r2d)
                 else if ( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = fileObj%p3di(k,j)%p
                    call mpp_write(unit, cur_var%field, r3d,                  tlev, is_root_pe=fileObj%is_root_pe)
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(save_restart): There is no pointer associated with the data of  field "// &
                         trim(cur_var%name)//" of file "//trim(fileObj%name) )
                 end if
              end if
           end if
        enddo ! end j loop
     enddo ! end k loop
     call mpp_close(unit)
     nfiles = nfiles+1
  enddo

  cur_var =>NULL()

end subroutine save_restart

!-------------------------------------------------------------------------------
!
!    This subroutine reads the model state from previously
!    generated files.  All restart variables are read from the first
!    file in the input filename list in which they are found.

subroutine restore_state_all(fileObj, directory)
  type(restart_file_type), intent(inout)       :: fileObj
  character(len=*),      intent(in), optional  :: directory

! Arguments: 
!  (in)      directory - The directory where the restart or save
!                        files should be found. The default is 'INPUT'

  character(len=128) :: dir    
  character(len=256) :: restartpath ! The restart file path (dir/file). 
  character(len=200) :: filepath    ! The path (dir/file) to the file being opened.
  character(len=8)   :: suffix      ! A suffix (like "_2") that is added to any
                                    ! additional restart files.
  character(len=80)  :: varname     ! A variable's name.
  integer            :: num_restart ! The number of restart files that have already
                                    ! been opened.
  integer            :: nfile       ! The number of files (restart files and others
                                    ! explicitly in filename) that are open.
  integer   :: unit(max_split_file) ! The mpp unit of all open files.
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: ndim, nvar, natt, ntime, tlev, siz
  type(fieldtype), allocatable        :: fields(:)
  logical                             :: fexist, domain_present
  integer                             :: j, n, l, k, missing_fields, domain_idx
  real, allocatable, dimension(:,:,:) :: r3d
  real, allocatable, dimension(:,:)   :: r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(restore_state_all): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  dir = 'INPUT'
  if(present(directory)) dir = directory

  num_restart = 0
  nfile = 0
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(fileObj%name)
  else
     restartpath = trim(fileObj%name)
  end if
  !--- first open all the restart files
  !--- NOTE: For distributed restart file, we are assuming there is only one file exist.
  inquire (file=trim(restartpath)//trim(pe_name), exist=fexist)     
  if(fexist) then
     nfile = 1
     call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
          fileset=MPP_MULTI)
  else
     do while(.true.)
        if (num_restart < 10) then
           write(suffix,'("_",I1)') num_restart
        else
           write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) then
           siz = len_trim(restartpath)
           if(restartpath(siz-2:siz) == ".nc") then
              filepath = restartpath(1:siz-3)//trim(suffix)
           else      
              filepath = trim(restartpath) // trim(suffix)
           end if
        else
           filepath = trim(restartpath)
        end if
        inquire (file=trim(filepath), exist=fexist) 
        if(.not. fexist) inquire(file=trim(filepath)//".nc", exist=fexist)
        if(fexist) then
           nfile = nfile + 1
           if(nfile > max_split_file) call mpp_error(FATAL, &
                "fms_io(restore_state_all): nfile is larger than max_split_file, increase max_split_file")
           call mpp_open(unit(nfile), trim(filepath), form=form,action=MPP_RDONLY,threading=thread_r, &
                fileset=MPP_SINGLE)
        else
           exit
        end if
        num_restart = num_restart + 1
     end do
  end if
  if(nfile == 0) call mpp_error(FATAL, "fms_io(restore_state_all): unable to find any restart files "// &
       "specified by "//trim(restartpath))


  ! Read each variable from the first file in which it is found.
  do n=1,nfile
     call mpp_get_info(unit(n), ndim, nvar, natt, ntime)

     allocate(fields(nvar))
     call mpp_get_fields(unit(n),fields(1:nvar))

     missing_fields = 0

     do j=1,fileObj%nvar
        cur_var => fileObj%var(j)
        domain_present = cur_var%domain_present
        domain_idx = cur_var%domain_idx
        do l=1, nvar
           call mpp_get_atts(fields(l),name=varname)
           if (lowercase(trim(varname)) == lowercase(trim(cur_var%name))) then
              cur_var%initialized = .true.
              do k = 1, cur_var%siz(4)
                 tlev = k
                 if(domain_present) then        
                    if( Associated(fileObj%p2dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p2dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p3dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = 0
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), r2d, tlev)
                       fileObj%p2di(k,j)%p = r2d
                       deallocate(r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = 0
                       call mpp_read(unit(n), fields(l), array_domain(domain_idx), r3d, tlev)
                       fileObj%p3di(k,j)%p = r3d
                       deallocate(r3d)
                    else
                       call mpp_error(FATAL, "fms_io(restore_state_all): domain is present for the field "//trim(varname)// &
                            " of file "//trim(fileObj%name)//", but none of p2dr, p3dr, p2di and p3di is associated") 
                    end if
                 else
                    if( Associated(fileObj%p0dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p0dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p1dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p1dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p2dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p2dr(k,j)%p, tlev)
                    else if( Associated(fileObj%p3dr(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), fileObj%p3dr(k,j)%p, tlev) 
                    else if( Associated(fileObj%p0di(k,j)%p) ) then
                       call mpp_read(unit(n), fields(l), r0d, tlev)
                       fileObj%p0di(k,j)%p = r0d
                    else if( Associated(fileObj%p1di(k,j)%p) ) then
                       allocate(r1d(cur_var%siz(1)) )
                       call mpp_read(unit(n), fields(l), r1d, tlev)                
                       fileObj%p1di(k,j)%p = r1d
                       deallocate(r1d)
                    else if( Associated(fileObj%p2di(k,j)%p) ) then
                       allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                       r2d = 0
                       call mpp_read(unit(n), fields(l), r2d, tlev)                
                       fileObj%p2di(k,j)%p = r2d
                       deallocate(r2d)
                    else if( Associated(fileObj%p3di(k,j)%p) ) then
                       allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                       r3d = 0
                       call mpp_read(unit(n), fields(l), r3d, tlev)                
                       fileObj%p3di(k,j)%p = r3d
                       deallocate(r3d)
                    else
                       call mpp_error(FATAL, "fms_io(restore_state_all): There is no pointer "//&
                            "associated with the data of  field "// trim(varname)//" of file "//trim(fileObj%name) )
                    end if
                 end if
              end do
              exit ! Start search for next restart variable.
           endif
        enddo
        if (l>nvar) missing_fields = missing_fields+1
     enddo

     deallocate(fields)
     if (missing_fields == 0) exit
  enddo

  do n=1,nfile
     call close_file(unit(n))
  enddo

  ! check whether all fields have been found
  do j = 1, fileObj%nvar
     if( .NOT. fileObj%var(j)%initialized ) then
        if( fileObj%var(j)%mandatory ) then
           call mpp_error(FATAL, "fms_io(restore_state_all): unable to find mandatory variable "// &
                trim(fileObj%var(j)%name)//" in restart file "//trim(fileObj%name) )
        end if
     end if
  end do
  cur_var =>NULL()

end subroutine restore_state_all

!-------------------------------------------------------------------------------
!
!    This subroutine reads the model state from previously
!    generated files.  All restart variables are read from the first
!    file in the input filename list in which they are found.

subroutine restore_state_one_field(fileObj, id_field, directory)
  type(restart_file_type), intent(inout)       :: fileObj
  integer,                 intent(in)          :: id_field
  character(len=*),      intent(in), optional  :: directory

! Arguments: 
!  (in)      directory - The directory where the restart or save
!                        files should be found. The default is 'INPUT'

  character(len=128) :: dir    
  character(len=256) :: restartpath ! The restart file path (dir/file). 
  character(len=200) :: filepath    ! The path (dir/file) to the file being opened.
  character(len=8)   :: suffix      ! A suffix (like "_2") that is added to any
                                    ! additional restart files.
  character(len=80)  :: varname     ! A variable's name.
  integer            :: num_restart ! The number of restart files that have already
                                    ! been opened.
  integer            :: nfile       ! The number of files (restart files and others
                                    ! explicitly in filename) that are open.
  integer   :: unit(max_split_file) ! The mpp unit of all open files.
  type(var_type), pointer, save       :: cur_var=>NULL()
  integer                             :: ndim, nvar, natt, ntime, tlev, siz
  type(fieldtype), allocatable        :: fields(:)
  logical                             :: fexist, domain_present
  integer                             :: j, n, l, k, missing_fields, domain_idx
  real, allocatable, dimension(:,:,:) :: r3d
  real, allocatable, dimension(:,:)   :: r2d
  real, allocatable, dimension(:)     :: r1d  
  real                                :: r0d

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(restore_state_one_field): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  dir = 'INPUT'
  if(present(directory)) dir = directory

  num_restart = 0
  nfile = 0
  if(len_trim(dir) > 0) then
     restartpath = trim(dir)//"/"// trim(fileObj%name)
  else
     restartpath = trim(fileObj%name)
  end if
  !--- first open all the restart files
  !--- NOTE: For distributed restart file, we are assuming there is only one file exist.
  inquire (file=trim(restartpath)//trim(pe_name), exist=fexist)     
  if(fexist) then
     nfile = 1
     call mpp_open(unit(nfile), trim(restartpath), form=form,action=MPP_RDONLY,threading=thread_r, &
          fileset=MPP_MULTI)
  else
     do while(.true.)
        if (num_restart < 10) then
           write(suffix,'("_",I1)') num_restart
        else
           write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) then
           siz = len_trim(restartpath)
           if(restartpath(siz-2:siz) == ".nc") then
              filepath = restartpath(1:siz-3)//trim(suffix)
           else      
              filepath = trim(restartpath) // trim(suffix)
           end if
        else
           filepath = trim(restartpath)
        end if
        inquire (file=trim(filepath), exist=fexist) 
        if(.not. fexist) inquire(file=trim(filepath)//".nc", exist=fexist)
        if(fexist) then
           nfile = nfile + 1
           if(nfile > max_split_file) call mpp_error(FATAL, &
                "fms_io(restore_state_one_field): nfile is larger than max_split_file, increase max_split_file")
           call mpp_open(unit(nfile), trim(filepath), form=form,action=MPP_RDONLY,threading=thread_r, &
                fileset=MPP_SINGLE)
        else
           exit
        end if
        num_restart = num_restart + 1
     end do
  end if
  if(nfile == 0) call mpp_error(FATAL, "fms_io(restore_state_one_field): unable to find any restart files "// &
       "specified by "//trim(restartpath))


  ! Read each variable from the first file in which it is found.
  do n=1,nfile
     call mpp_get_info(unit(n), ndim, nvar, natt, ntime)

     allocate(fields(nvar))
     call mpp_get_fields(unit(n),fields(1:nvar))

     missing_fields = 0
     j = id_field
     cur_var => fileObj%var(j)
     domain_present = cur_var%domain_present
     domain_idx = cur_var%domain_idx
     do l=1, nvar
        call mpp_get_atts(fields(l),name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(cur_var%name))) then
           cur_var%initialized = .true.
           do k = 1, cur_var%siz(4)
              tlev = k
              if(domain_present) then        
                 if( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p2dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), fileObj%p3dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = 0
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), r2d, tlev)
                    fileObj%p2di(k,j)%p = r2d
                    deallocate(r2d)
                 else if( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = 0
                    call mpp_read(unit(n), fields(l), array_domain(domain_idx), r3d, tlev)
                    fileObj%p3di(k,j)%p = r3d
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(restore_state_one_field): domain is present for the field "//trim(varname)// &
                         " of file "//trim(fileObj%name)//", but none of p2dr, p3dr, p2di and p3di is associated") 
                 end if
              else
                 if( Associated(fileObj%p0dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p0dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p1dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p1dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p2dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p2dr(k,j)%p, tlev)
                 else if( Associated(fileObj%p3dr(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), fileObj%p3dr(k,j)%p, tlev) 
                 else if( Associated(fileObj%p0di(k,j)%p) ) then
                    call mpp_read(unit(n), fields(l), r0d, tlev)
                    fileObj%p0di(k,j)%p = r0d
                 else if( Associated(fileObj%p1di(k,j)%p) ) then
                    allocate(r1d(cur_var%siz(1)) )
                    call mpp_read(unit(n), fields(l), r1d, tlev)                
                    fileObj%p1di(k,j)%p = r1d
                    deallocate(r1d)
                 else if( Associated(fileObj%p2di(k,j)%p) ) then
                    allocate(r2d(cur_var%siz(1), cur_var%siz(2)) )
                    r2d = 0
                    call mpp_read(unit(n), fields(l), r2d, tlev)                
                    fileObj%p2di(k,j)%p = r2d
                    deallocate(r2d)
                 else if( Associated(fileObj%p3di(k,j)%p) ) then
                    allocate(r3d(cur_var%siz(1), cur_var%siz(2), cur_var%siz(3)) )
                    r3d = 0
                    call mpp_read(unit(n), fields(l), r3d, tlev)                
                    fileObj%p3di(k,j)%p = r3d
                    deallocate(r3d)
                 else
                    call mpp_error(FATAL, "fms_io(restore_state_one_field): There is no pointer "// &
                         "associated with the data of  field "//trim(varname)//" of file "//trim(fileObj%name) )
                 end if
              end if
           end do
           exit ! Start search for next restart variable.
        endif
     enddo
     if (l>nvar) missing_fields = missing_fields+1
     deallocate(fields)
     if (missing_fields == 0) exit
  enddo

  do n=1,nfile
     call close_file(unit(n))
  enddo

  ! check whether the field have been found
  if( .NOT. fileObj%var(id_field)%initialized ) then
     if( fileObj%var(id_field)%mandatory ) then
        call mpp_error(FATAL, "fms_io(restore_state_one_field): unable to find mandatory variable "// &
             trim(fileObj%var(id_field)%name)//" in restart file "//trim(fileObj%name) )
     end if
  end if
  cur_var =>NULL()

end subroutine restore_state_one_field

!-------------------------------------------------------------------------------
!
!     This routine will setup one entry to be written out 
!
!-------------------------------------------------------------------------------
subroutine setup_one_field(fileObj, filename, fieldname, field_siz, index_field,  domain, instance_name, mandatory, &
                           no_domain, position, tile_count, data_default, longname, units, mosaicfile)
  type(restart_file_type), intent(inout)         :: fileObj
  character(len=*),         intent(in)           :: filename, fieldname 
  integer, dimension(:),    intent(in)           :: field_siz
  integer,                  intent(out)          :: index_field
  type(domain2d), optional, intent(in), target   :: domain
  real,           optional, intent(in)           :: data_default
  logical,        optional, intent(in)           :: no_domain   
  integer,        optional, intent(in)           :: position, tile_count
  character(len=*), optional, intent(in)         :: instance_name
  logical,          optional, intent(in)         :: mandatory
  character(len=*), optional, intent(in)         :: longname, units
  character(len=*), optional, intent(in)         :: mosaicfile

  !--- local variables
  integer                         :: i, domain_idx
  integer                         :: ishift, jshift
  integer                         :: gxsize, gysize
  integer                         :: cxsize, cysize
  integer                         :: dxsize, dysize
  real                            :: default_data
  logical                         :: is_no_domain = .false.
  character(len=256)              :: fname, filename2
  type(domain2d), pointer, save   :: d_ptr   =>NULL()
  type(var_type), pointer, save   :: cur_var =>NULL()
  integer                         :: length

  if(ANY(field_siz < 1)) then
     call mpp_error(FATAL, "fms_io(setup_one_field): each entry of field_size should be a positive integer")
  end if

  if(PRESENT(data_default))then
     default_data=data_default
  else
     default_data=0.
  endif

  if(present(tile_count) .AND. .not. present(domain)) call mpp_error(FATAL, &
         'fms_io(setup_one_field): when tile_count is present, domain must be present')

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     is_no_domain = no_domain
  end if
  if(is_no_domain) then
!     comments the following to allow passing domain when writing 1-D or scalar variables. 
!     if(PRESENT(domain)) &
!       call mpp_error(FATAL, 'fms_io(setup_one_field): no_domain cannot be .true. when optional argument domain is present.')
  else if(PRESENT(domain))then
     d_ptr => domain
  else if (ASSOCIATED(Current_domain)) then
     d_ptr => Current_domain
  endif

  !--- remove .nc from file name
  length = len_trim(filename)
  if(filename(length-2:length) == '.nc') then
     filename2 = filename(1:length-3)
  else
     filename2 = filename(1:length)
  end if

  if(present(instance_name)) then
     if(len_trim(instance_name) > 0) filename2 = trim(filename2)//'.'//trim(instance_name)
  end if

  !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count,
  !JWD:  I just don't see how the physics can remain "tile neutral"
  !z1l:  one solution is add one more public interface called set_tile_count
  call get_mosaic_tile_file(filename2, fname, is_no_domain, domain, tile_count, mosaicfile)

  if(Associated(fileObj%var) ) then
     ! make sure the consistency of file name
     if(trim(fileObj%name) .NE. trim(fname)) call mpp_error(FATAL, 'fms_io(setup_one_field): filename = '// &
         trim(fname)//' is not consistent with the filename of the restart object = '//trim(fileObj%name) )
  else
     allocate(fileObj%var(max_fields) )
     allocate(fileObj%p0dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p1dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p2dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p3dr(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p0di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p1di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p2di(MAX_TIME_LEVEL_REGISTER, max_fields))
     allocate(fileObj%p3di(MAX_TIME_LEVEL_REGISTER, max_fields))
     !--- make sure fname is not used in other restart_file_type object.
     do i = 1, num_registered_files
        if(trim(fname) == trim(registered_file(i)) ) call mpp_error(FATAL, &
          'fms_io(setup_one_field): '//trim(fname)//' is already registered with other restart_file_type data')
     end do
     num_registered_files = num_registered_files + 1
     registered_file(num_registered_files) = trim(fname)
     fileObj%name = trim(fname)         
     fileObj%tile_count=1
     if(present(tile_count)) fileObj%tile_count = tile_count
     if(ASSOCIATED(d_ptr))then
        fileObj%is_root_pe = mpp_domain_is_root_pe(d_ptr)
     else
        fileObj%is_root_pe = mpp_pe() == mpp_root_pe()
     endif
     fileObj%max_ntime = field_siz(4)
     fileObj%nvar      = 0
     !-- allocate memory
     do i = 1, max_fields
        fileObj%var(i)%name           = 'none'
        fileObj%var(i)%domain_present = .false.
        fileObj%var(i)%domain_idx     = -1
        fileObj%var(i)%is_dimvar      = .false.
        fileObj%var(i)%position       = CENTER
        fileObj%var(i)%siz(:)         = 0
        fileObj%var(i)%gsiz(:)        = 0
        fileObj%var(i)%id_axes(:)     = -1
        fileObj%var(i)%longname       = "";
        fileObj%var(i)%units          = "none";
        fileObj%var(i)%mandatory      = .true.
        fileObj%var(i)%initialized    = .false.
     end do
  endif

  ! check if the field is new or not and get position and dimension of the field
  index_field = -1
  do i = 1, fileObj%nvar
     if(trim(fileObj%var(i)%name) == trim(fieldname)) then
        index_field = i
        exit
     end if 
  end do

  if(index_field > 0) then
     cur_var   => fileObj%var(index_field)
     if(cur_var%siz(1) .NE. field_siz(1) .OR. cur_var%siz(2) .NE. field_siz(2) .OR. cur_var%siz(3) .NE. field_siz(3) ) &
        call mpp_error(FATAL, 'fms_io(setup_one_field): field size mismatch for field '// &
                       trim(fieldname)//' of file '//trim(filename) )

     cur_var%siz(4) =  cur_var%siz(4) + field_siz(4)
     if(fileObj%max_ntime < cur_var%siz(4) ) fileObj%max_ntime = cur_var%siz(4)
     ! the time level should be no larger than MAX_TIME_LEVEL_REGISTER ( = 2) 
     if( cur_var%siz(4) > MAX_TIME_LEVEL_REGISTER ) call mpp_error(FATAL, 'fms_io(setup_one_field): ' // &
          'the time level of field '//trim(cur_var%name)//' in file '//trim(fileObj%name)// &
          ' is greater than MAX_TIME_LEVEL_REGISTER(=2), increase MAX_TIME_LEVEL_REGISTER or check your code')
  else 
     fileObj%nvar = fileObj%nvar +1
     if(fileObj%nvar>max_fields) then
        write(error_msg,'(I3,"/",I3)') fileObj%nvar, max_fields 
        call  mpp_error(FATAL,'fms_io(setup_one_field): max_fields exceeded, needs increasing, nvar/max_fields=' &
             //trim(error_msg))
     endif
     index_field =  fileObj%nvar
     cur_var   => fileObj%var(index_field)
     cur_var%siz(:)  = field_siz(:)
     cur_var%gsiz(3) = field_siz(3)
     cur_var%csiz(3) = field_siz(3)
     cur_var%name = fieldname
     if(present(mandatory)) cur_var%mandatory = mandatory
     if(present(longname)) then
        cur_var%longname = longname
     else
        cur_var%longname = fieldname
     end if
     if(present(units))    cur_var%units    = units
     if(present(position)) cur_var%position = position
     if(ASSOCIATED(d_ptr)) then
        cur_var%domain_present = .true.
        domain_idx = lookup_domain(d_ptr)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io(setup_one_field), 1: max_domains exceeded,' &
                //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = d_ptr
           call mpp_get_domain_components(array_domain(domain_idx), domain_x(domain_idx), domain_y(domain_idx), &
                tile_count=tile_count)
        endif
        cur_var%domain_idx = domain_idx
        call mpp_get_domain_shift ( array_domain(domain_idx), ishift, jshift, position)
        call mpp_get_global_domain(array_domain(domain_idx), xsize=gxsize,ysize=gysize,tile_count=tile_count)
        call mpp_get_compute_domain(array_domain(domain_idx), xsize = cxsize, ysize = cysize, tile_count=tile_count)
        call mpp_get_data_domain   (array_domain(domain_idx), xsize = dxsize, ysize = dysize, tile_count=tile_count)
        if (ishift .NE. 0) then
           cxsize = cxsize+ishift; dxsize = dxsize+ishift; gxsize = gxsize + ishift
        end if
        if (jshift .NE. 0) then
           cysize = cysize+jshift; dysize = dysize+jshift; gysize = gysize + jshift
        endif
        if( (cur_var%siz(1) .NE. cxsize .AND. cur_var%siz(1) .NE. dxsize ) .OR. &
            (cur_var%siz(2) .NE. cysize .AND. cur_var%siz(2) .NE. dysize ) ) then
            call mpp_error(FATAL, 'fms_io(setup_one_field): data should be on either computer domain '//&
              'or data domain when domain is present for field '//trim(fieldname)//' of file '//trim(filename) )
        end if

        cur_var%gsiz(1)   = gxsize
        cur_var%gsiz(2)   = gysize
        if(thread_w == MPP_MULTI) then
           call mpp_get_compute_domain(array_domain(domain_idx), xsize=cxsize,ysize=cysize,tile_count=tile_count)
           cur_var%csiz(1)   = cxsize
           cur_var%csiz(2)   = cysize
        else
           cur_var%csiz(1)   = cur_var%gsiz(1)
           cur_var%csiz(2)   = cur_var%gsiz(2)
        end if
     else
        cur_var%gsiz(1:2) = field_siz(1:2)
        cur_var%csiz(1:2) = field_siz(1:2)
        cur_var%domain_present=.false.
     endif
  end if

  d_ptr =>NULL()
  cur_var =>NULL()

end subroutine setup_one_field


!.....................................................................
subroutine write_data_2d_new(filename, fieldname, data, domain,append_pelist_name, &
                             no_domain, position,tile_count, data_default, mosaicfile)

  character(len=*), intent(in)                 :: filename, fieldname 
  real, dimension(:,:), intent(in)             :: data
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  real, intent(in), optional                   :: data_default
  type(domain2d), intent(in), optional         :: domain
  logical, intent(in), optional                :: append_pelist_name, no_domain
  integer, intent(in), optional                :: position, tile_count
  character(len=*), intent(in), optional       :: mosaicfile
 
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_2d_new):need to call fms_io_init first')
  data_3d(:,:,1) = data(:,:)

  call write_data_3d_new(filename, fieldname, data_3d, domain, append_pelist_name, &
                         no_domain, position, tile_count, data_default, mosaicfile)

end subroutine write_data_2d_new

! ........................................................
subroutine write_data_1d_new(filename, fieldname, data,domain,append_pelist_name, &
                             no_domain, tile_count, data_default, mosaicfile)
  
  type(domain2d), intent(in), optional   :: domain
  character(len=*), intent(in)           :: filename, fieldname 
  real, dimension(:), intent(in)         :: data
  real, dimension(size(data(:)),1,1)     :: data_3d
  real, intent(in), optional             :: data_default
  logical, intent(in), optional          :: append_pelist_name, no_domain
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile

  
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_1d_new): module not initialized')  
  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(write_data_1d_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if
  data_3d(:,1,1) = data(:)  
  call write_data_3d_new(filename, fieldname, data_3d,domain,append_pelist_name, &
                         no_domain=.true., tile_count=tile_count, data_default=data_default, mosaicfile=mosaicfile)  
end subroutine write_data_1d_new

! ..........................................................
subroutine write_data_scalar_new(filename, fieldname, data, domain, append_pelist_name, &
                                 no_domain, tile_count, data_default, mosaicfile)

  type(domain2d), intent(in), optional   :: domain
  character(len=*), intent(in)           :: filename, fieldname 
  real, intent(in)                       :: data
  real, dimension(1,1,1)                 :: data_3d
  real, intent(in), optional             :: data_default
  logical, intent(in), optional          :: append_pelist_name, no_domain
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile
    
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_scalar_new):  module not initialized: '//fieldname)  
  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(write_data_scalar_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  data_3d(1,1,1) = data
  call write_data_3d_new(filename, fieldname, data_3d,domain,append_pelist_name, &
                         no_domain=.true., tile_count=tile_count, data_default=data_default, mosaicfile=mosaicfile)
end subroutine write_data_scalar_new

! ..........................................................

function lookup_field_r(nfile,fieldname)
! Given fieldname, this function returns the field position in the model's fields list

  integer, intent(in)          :: nfile
  character(len=*), intent(in) :: fieldname
  integer                      :: lookup_field_r
  integer                      :: j

  lookup_field_r=-1
  do j = 1, files_read(nfile)%nvar
     if (trim(files_read(nfile)%var(j)%name) == trim(fieldname)) then
        lookup_field_r = j
        exit 
     endif
  enddo
  return
end function lookup_field_r


!..........................................................

function lookup_domain(domain)
! given domain, this function returns the position of domain in array_domain or -1 if not found

  type(domain2d), intent(in) :: domain
  integer                    :: i, lookup_domain
  lookup_domain = -1
  do i =1, num_domains
     if(domain == array_domain(i)) then
        lookup_domain = i
        exit
     endif
  enddo
end function lookup_domain
!.........................................................
function lookup_axis(axis_sizes,siz,domains,dom)

! Given axis size (global), this function returns the axis id  

  integer, intent(in)      :: axis_sizes(:), siz
  type(domain1d), optional :: domains(:)
  type(domain1d), optional :: dom
  integer :: lookup_axis
  integer :: j


  lookup_axis=-1
  do j=1,size(axis_sizes(:))
     if (siz == axis_sizes(j)) then
        if (PRESENT(domains)) then
           if (dom == domains(j)) then 
              lookup_axis = j
              exit
           endif
        else
           lookup_axis = j
           exit
        endif
     endif
  enddo
  if (lookup_axis == -1) call mpp_error(FATAL,'fms_io(lookup_axis): could not find axis in set of axes')
end function lookup_axis
!.....................................................................
! <SUBROUTINE NAME="field_size">
!<DESCRIPTION>
! Given filename and fieldname, this subroutine returns the size of field
!</DESCRIPTION>
!   <TEMPLATE>
! call field_size(filename, fieldname, siz) 
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <OUT NAME="siz" TYPE="integer" DIM="(*)">
!    siz must be a dimension(4) array to retrieve the size of the field
!   </OUT>
!   <OUT NAME="field_found" TYPE="logical, optional">
!    if this flag is present, field_size will not abort if
!    called for a non-existent field.
!    Instead it will return T or F depending on
!    whether or not the field was found.
!   </OUT>
subroutine field_size(filename, fieldname, siz, append_pelist_name, field_found, domain )

  character(len=*), intent(in)         :: filename, fieldname
  integer,       intent(inout)         :: siz(:)
  logical,        intent(in), optional :: append_pelist_name
  logical,       intent(out), optional :: field_found  
  type(domain2d), intent(in), optional :: domain
  integer                              :: i, nfile, unit
  logical                              :: file_opened, found, is_exist
  character(len=256)                   :: fname
  if (size(siz(:)) < 4) call mpp_error(FATAL,'fms_io(field_size): size array must be >=4 to receive field size of ' &
       //trim(fieldname)//' in file '// trim(filename))

! Need to check if filename has been opened or not

  fname = trim(filename)
  if (PRESENT(append_pelist_name)) then
     if (append_pelist_name) then
        fname = trim(filename)//trim(pelist_name)
     endif
  endif        
  nfile = 0
  file_opened=.false.
  do i=1,num_files_r
     if (trim(files_read(i)%name) == trim(fname))  then
        nfile = i
        file_opened = .true.
        exit ! file is already opened
     endif
  enddo
!Need to open the file now, Only works for single NetCDF files for now ...
  found= .false.
  siz=-1
  if (.not. file_opened) then
     inquire (file=trim(fname)//trim(pe_name), exist=is_exist)
     if(.not. is_exist) inquire (file=trim(fname)//'.nc'//trim(pe_name), exist=is_exist)
     if(is_exist) then
        call mpp_open(unit,trim(fname),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_MULTI, &
             fileset=MPP_MULTI)
     else 
        inquire (file=trim(fname), exist=is_exist)
        if(.not. is_exist) inquire (file=trim(fname)//'.nc', exist=is_exist)
        if(is_exist) then     
           call mpp_open(unit,trim(fname),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_SINGLE, &
                fileset=MPP_SINGLE)
        end if
     end if
     if(is_exist) then
        call get_size(unit,fieldname,siz,found)
        call mpp_close(unit)
     end if
  else
     do i=1, files_read(nfile)%nvar
        if (trim(fieldname) == trim(files_read(nfile)%var(i)%name)) then
           found = .true.
           siz = files_read(nfile)%var(i)%siz(:)
           exit
        endif
     enddo
     if (.not. found) then
        call get_size(files_read(nfile)%unit,fieldname,siz,found)
     endif
 endif

 if (.not. found) then 
! Perhaps the variable is in a "tile" file.
     call get_mosaic_tile_file(filename, fname, is_no_domain= .false., domain=domain)

     nfile = 0
     file_opened=.false.
     do i=1,num_files_r
        if (trim(files_read(i)%name) == trim(fname))  then
           nfile = i
           file_opened = .true.
           exit ! file is already opened
        endif
     enddo
     found= .false.
     siz=-1
     if (.not. file_opened) then
        ! File is not open yet. 
        inquire (file=trim(fname)//trim(pe_name), exist=is_exist)
        if(is_exist) then
           call mpp_open(unit,trim(fname),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_MULTI, &
                fileset=MPP_MULTI)
        else
           inquire (file=trim(fname), exist=is_exist)
           if(.not. is_exist) call mpp_error(FATAL, 'fms_io(field_size): file '//trim(filename)// &
                ' does not exist with the consideration tile number and distribute file')
           call mpp_open(unit,trim(fname),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_SINGLE, &
                fileset=MPP_SINGLE)
        end if
        call get_size(unit, fieldname, siz, found)
        call mpp_close(unit)
     else
        do i=1, files_read(nfile)%nvar
           if (trim(fieldname) == trim(files_read(nfile)%var(i)%name)) then
              found = .true.
              siz = files_read(nfile)%var(i)%siz(:)
              exit
           endif
        enddo
! Scan the already open file in order to get the size of the field.
        if (.not. found) then
           call get_size(files_read(nfile)%unit, fieldname, siz, found)
        endif
     endif
  endif

  if( PRESENT(field_found) )then
     field_found = found
  else if (.not. found .and. mpp_pe() == mpp_root_pe() )then
     call mpp_error(FATAL, 'fms_io(field_size): field '//trim(fieldname)//' NOT found in file '//trim(filename))
  end if

  return
end subroutine field_size
! </SUBROUTINE>

subroutine get_size(unit, fieldname, siz, found)
integer,          intent(in)    :: unit
character(len=*), intent(in)    :: fieldname
integer,          intent(inout) :: siz(:)
logical,          intent(out)   :: found

  character(len=128)             :: name
  character(len=1)               :: cart
  integer                        :: i, ndim, nvar, natt, ntime, siz_in(4), j, len
  type(fieldtype)                :: fields(max_fields)
  type(axistype)                 :: axes(max_fields)
     found = .false.
     call mpp_get_info(unit,ndim,nvar,natt,ntime)
     if (nvar > max_fields) then
        write(error_msg,'(I3,"/",I3)') nvar,max_fields 
        call  mpp_error(FATAL,'fms_io(field_size): max_fields too small, needs increasing, nvar/max_fields=' &
             //trim(error_msg))!//' in file '//trim(filename))
     endif
     call mpp_get_fields(unit,fields(1:nvar))
     do i=1, nvar
        call mpp_get_atts(fields(i),name=name)
        if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
           call mpp_get_atts(fields(i),ndim=ndim)
           call mpp_get_atts(fields(i),axes=axes(1:ndim))
           call mpp_get_atts(fields(i),siz=siz_in)
           siz = siz_in
           siz(4) = ntime
           if(ndim == 1) then
              call mpp_get_atts(axes(1), len=siz(1))
           end if
           do j = 1, ndim
              call mpp_get_atts(axes(j),len=len)
              call get_axis_cart(axes(j),cart)
              select case (cart)
              case ('X')
                 siz(1) = len
              case('Y')
                 siz(2) = len
              case('Z')
                 siz(3) = len
              case('T')
                 siz(4) = len
              end select
           enddo
           found = .true.
           exit
        endif
     enddo

     if(.not. found) then
        call mpp_get_axes(unit,axes(1:ndim))
        do i=1, ndim
           call mpp_get_atts(axes(i),name=name, len= siz_in(1))
           if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
              siz(1)= siz_in(1)
              found = .true.
              exit
           endif
        enddo
     endif
end subroutine get_size

! <SUBROUTINE NAME="read_data">
!<DESCRIPTION>
! This routine performs reading "fieldname" stored in "filename". The data values of fieldname
! will be stored in "data" at the end of this routine. For fieldname with multiple timelevel
! just repeat the routine with explicit timelevel in each call.
!</DESCRIPTION>
!   <TEMPLATE>
! call read_data(filename,fieldname,data,domain,timelevel)
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character" DIM="(*)">
!    File name
!   </IN>
!   <IN NAME="fieldname" TYPE="character" DIM="(*)">
!    Field  name
!   </IN>
!   <IN NAME="domain"  TYPE="domain, optional">
!   domain of fieldname
!   </IN>
!   <IN NAME="timelevel" TYPE="integer, optional">
!     time level of fieldname
!   </IN>
!   <OUT NAME="data"  TYPE="real">
!   array containing data of fieldname
!   </OUT>
!=====================================================================================
subroutine read_data_i3d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, &
                             no_domain,position, tile_count, mosaicfile)
  character(len=*),           intent(in)   :: filename, fieldname
  integer, dimension(:,:,:), intent(inout) :: data ! 3 dimensional data    
  type(domain2d), intent(in),   optional   :: domain
  integer, intent(in),          optional   :: timelevel
  logical, intent(in),          optional   :: append_pelist_name, no_domain 
  integer, intent(in) ,         optional   :: position, tile_count
  character(len=*), intent(in), optional   :: mosaicfile

  real, dimension(size(data,1),size(data,2),size(data,3)) :: r_data
  r_data = 0
  call read_data_3d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, &
                        no_domain,position, tile_count, mosaicfile)
  data = CEILING(r_data)
end subroutine read_data_i3d_new

subroutine read_data_i2d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, &
                             no_domain,position, tile_count, mosaicfile)
  character(len=*),         intent(in)   :: filename, fieldname
  integer, dimension(:,:), intent(inout) :: data ! 2 dimensional data    
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in),        optional   :: timelevel
  logical, intent(in),        optional   :: append_pelist_name , no_domain
  integer, intent(in) ,       optional   :: position, tile_count
  character(len=*), intent(in), optional :: mosaicfile
  real, dimension(size(data,1),size(data,2)) :: r_data

  r_data = 0
  call read_data_2d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, &
                        no_domain,position, tile_count, mosaicfile)
  data = CEILING(r_data)
end subroutine read_data_i2d_new
!.....................................................................
subroutine read_data_i1d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, &
                             no_domain, tile_count, mosaicfile)
  character(len=*), intent(in)           :: filename, fieldname
  integer, dimension(:), intent(inout)   :: data ! 1 dimensional data    
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: append_pelist_name, no_domain 
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile

  real, dimension(size(data,1))        :: r_data

  call read_data_1d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, &
                        no_domain, tile_count, mosaicfile)
  data = CEILING(r_data)
end subroutine read_data_i1d_new
!.....................................................................
subroutine read_data_iscalar_new(filename,fieldname,data,domain,timelevel,append_pelist_name, &
                                 no_domain, tile_count, mosaicfile)
  character(len=*), intent(in)           :: filename, fieldname
  integer, intent(inout)                 :: data     
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: append_pelist_name, no_domain 
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile

  real                                 :: r_data
  call read_data_scalar_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, &
                            no_domain, tile_count, mosaicfile)
  data = CEILING(r_data)
end subroutine read_data_iscalar_new
!=====================================================================================
subroutine read_data_3d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, &
      no_domain, position, tile_count, mosaicfile)
  character(len=*),                  intent(in) :: filename, fieldname
  real, dimension(:,:,:),         intent(inout) :: data ! 3 dimensional data    
  type(domain2d), target, optional,  intent(in) :: domain
  integer,                optional,  intent(in) :: timelevel
  logical,                optional,  intent(in) :: append_pelist_name, no_domain 
  integer,                optional,  intent(in) :: position, tile_count
  character(len=*),       optional,  intent(in) :: mosaicfile

  character(len=256)            :: fname
  integer                       :: unit, siz_in(4)
  integer                       :: file_index  ! index of the opened file in array files
  integer                       :: tlev=1
  integer                       :: index_field ! position of the fieldname in the list of variables
  integer                       :: cxsize, cysize
  integer                       :: dxsize, dysize
  integer                       :: gxsize, gysize
  integer                       :: ishift, jshift
  logical                       :: is_no_domain = .false.
  logical                       :: read_dist
  type(domain2d), pointer, save :: d_ptr =>NULL()

! read disttributed files is used when reading restart files that are NOT mppnccombined. In this
! case PE 0 will read file_res.nc.0000, PE 1 will read file_res.nc.0001 and so forth.
! 
! namelist to be used with read_dist_files: threading_read=multi,
! threading_write=multi, fileset_write=multi.

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_3d_new):  module not initialized')  
  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
!z1l:  comment out the following error check to allow reading 1-D or scalar variable with passing domain.
!!$     if(PRESENT(domain) .AND. no_domain) &
!!$       call mpp_error(FATAL, 'fms_io(read_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')
     is_no_domain = no_domain
  endif  
 
  if(PRESENT(domain))then
     d_ptr => domain
  elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
     d_ptr => Current_domain
  endif

  if(.not. PRESENT(domain) .and. .not. ASSOCIATED(Current_domain) ) is_no_domain = .true.

  call get_file_name(filename, fname, read_dist, is_no_domain, domain, append_pelist_name, tile_count, mosaicfile)
  call get_file_unit(fname, unit, file_index, read_dist)

  siz_in(3) = size(data,3)
  if( read_dist .or. is_no_domain) then
     if(associated(d_ptr) .AND. (.NOT. is_no_domain) ) then  !-- read_dist will be true, global size will be compute domain size.
        call mpp_get_compute_domain(d_ptr, xsize = gxsize, ysize = gysize, tile_count=tile_count)
        call mpp_get_domain_shift  (d_ptr, ishift, jshift, position)
        if (ishift .NE. 0)  gxsize = gxsize+ishift
        if (jshift .NE. 0)  gysize = gysize+jshift
     else
        gxsize = size(data,1)
        gysize = size(data,2)
     end if
  else  ! d_ptr must not be null
     call mpp_get_compute_domain(d_ptr, xsize = cxsize, ysize = cysize, tile_count=tile_count)
     call mpp_get_data_domain   (d_ptr, xsize = dxsize, ysize = dysize, tile_count=tile_count)
     call mpp_get_global_domain (d_ptr, xsize = gxsize, ysize = gysize, tile_count=tile_count)
     call mpp_get_domain_shift  (d_ptr, ishift, jshift, position)
     if (ishift .NE. 0) then
        cxsize = cxsize+ishift; dxsize = dxsize+ishift; gxsize = gxsize+ishift
     endif
     if (jshift .NE. 0) then
        cysize = cysize+jshift; dysize = dysize+jshift; gysize = gysize+jshift
     endif
     if( (size(data,1) .NE. cxsize .AND. size(data,1) .NE. dxsize) .OR. &
         (size(data,2) .NE. cysize .AND. size(data,2) .NE. dysize) )then
       call mpp_error(FATAL,'fms_io(read_data_3d_new): data should be on either computer domain '//&
                            'or data domain when domain is present. '//&
                            'shape(data)=',shape(data),'  cxsize,cysize,dxsize,dysize=',(/cxsize,cysize,dxsize,dysize/))
     end if 
  endif

  if (PRESENT(timelevel)) then
     tlev = timelevel
  else
     tlev = 1
  endif  

  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call get_field_id(unit, file_index, fieldname, index_field, is_no_domain, .false. )
     siz_in = files_read(file_index)%var(index_field)%siz
     if(files_read(file_index)%var(index_field)%is_dimvar ) then
        if (.not. read_dist) then
           if (siz_in(1) /= gxsize) &
                call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                ' in file '//trim(filename)//' field size mismatch 2')
        endif
     else
        if (siz_in(1) /= gxsize .or. siz_in(2) /= gysize .or. siz_in(3) /= size(data,3)) then
           PRINT *, gxsize, gysize, size(data, 3), siz_in(1), siz_in(2), siz_in(3)
           call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                ' in file '//trim(filename)//': field size mismatch 1')
        endif
     end if
     if ( tlev < 1 .or. files_read(file_index)%max_ntime < tlev)  then
        write(error_msg,'(I5,"/",I5)') tlev, files_read(file_index)%max_ntime
        call mpp_error(FATAL,'fms_io(read_data_3d_new): time level out of range, time level/max_time_level=' &
             //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
     endif

     
     if(is_no_domain) then
        if (files_read(file_index)%var(index_field)%is_dimvar) then
           call mpp_get_axis_data(files_read(file_index)%var(index_field)%axis,data(:,1,1))
        else
           call mpp_read(unit,files_read(file_index)%var(index_field)%field,data(:,:,:),tlev)
        endif
     else 
        call mpp_read(unit,files_read(file_index)%var(index_field)%field,d_ptr,data,tlev,tile_count)
     endif
  endif  

  d_ptr =>NULL()

  return
end subroutine read_data_3d_new

!=====================================================================================
!--- we assume any text data are at most 2-dimensional and level is for first dimension
subroutine read_data_text(filename,fieldname,data,level)
  character(len=*), intent(in)   :: filename, fieldname
  character(len=*), intent(out)  :: data
  integer, intent(in) , optional :: level
  logical                        :: file_opened, fexist, read_dist      
  integer                        :: lev, unit, index_field
  integer                        :: file_index
  character(len=256)             :: fname

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_text):  module not initialized')  

  file_opened=.false.
  fexist = .false.
  if (PRESENT(level)) then
     lev = level
  else
     lev = 1
  endif  

  call get_file_name(filename, fname, read_dist, .true. )
  call get_file_unit(fname, unit, file_index, read_dist)

! Get info of this file and field   
  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call get_field_id(unit, file_index, fieldname, index_field, .true., .true. )     

     if ( lev < 1 .or. lev > files_read(file_index)%var(index_field)%siz(1) )  then
        write(error_msg,'(I5,"/",I5)') lev, files_read(file_index)%var(index_field)%siz(1)
        call mpp_error(FATAL,'fms_io(read_data_text): text level out of range, level/max_level=' &
             //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
     endif

     call mpp_read(unit,files_read(file_index)%var(index_field)%field,data, level=level)
  endif  
  return
end subroutine read_data_text
!.............................................................. 
! </SUBROUTINE>
subroutine read_data_2d_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain,position,tile_count, mosaicfile)
  character(len=*), intent(in)                 :: filename, fieldname
  real, dimension(:,:), intent(inout)          :: data     !2 dimensional data 
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  type(domain2d), intent(in), optional         :: domain
  integer, intent(in) , optional               :: timelevel
  logical, intent(in), optional                :: append_pelist_name, no_domain
  integer, intent(in) , optional               :: position, tile_count
  character(len=*), intent(in), optional       :: mosaicfile
  integer                                      :: isc,iec,jsc,jec,isd,ied,jsd,jed
  integer :: isg,ieg,jsg,jeg
  integer                                      :: xsize_c,ysize_c,xsize_d,ysize_d
  integer                                      :: xsize_g,ysize_g, ishift, jshift

!#ifdef use_CRI_pointers
!  pointer( p, data_3d )
!  p = LOC(data)
!#endif
  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain,position,tile_count, mosaicfile)

  if(PRESENT(domain)) then
     call mpp_get_global_domain( domain,isg,ieg,jsg,jeg,xsize=xsize_g,ysize=ysize_g, tile_count=tile_count)
     call mpp_get_compute_domain( domain,isc,iec,jsc,jec,xsize=xsize_c,ysize=ysize_c, tile_count=tile_count)
     call mpp_get_data_domain( domain,isd,ied,jsd,jed,xsize=xsize_d,ysize=ysize_d, tile_count=tile_count)
     call mpp_get_domain_shift  (domain, ishift, jshift, position)
     if((size(data,1)==xsize_c+ishift) .and. (size(data,2)==ysize_c+jshift)) then !on_comp_domain
        data(:,:) = data_3d(:,:,1)
     else if((size(data,1)==xsize_d+ishift) .and. (size(data,2)==ysize_d+jshift)) then !on_data_domain
        data(isc-isd+1:iec-isd+1+ishift,jsc-jsd+1:jec-jsd+1+jshift) = &
               data_3d(isc-isd+1:iec-isd+1+ishift,jsc-jsd+1:jec-jsd+1+jshift,1)
     else if((size(data,1)==xsize_g+ishift) .and. (size(data,2)==ysize_g+jshift)) then !on_global_domain
        data(:,:) = data_3d(:,:,1)
     else 
        call mpp_error(FATAL,'error in read_data_2d_new, field '//trim(fieldname)// &
                      ' in file '//trim(filename)//' data must be in compute or data domain')
     endif
  else     
     data(:,:) = data_3d(:,:,1)
  endif

end subroutine read_data_2d_new
!.....................................................................
subroutine read_data_1d_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain, tile_count, mosaicfile)
  character(len=*), intent(in)           :: filename, fieldname
  real, dimension(:), intent(inout)      :: data     !1 dimensional data 
  real, dimension(size(data,1),1,1)      :: data_3d
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: append_pelist_name, no_domain
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile  
#ifdef use_CRI_pointers
  pointer( p, data_3d )
  p = LOC(data)
#endif

  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(read_data_1d_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain=.true., tile_count=tile_count, mosaicfile=mosaicfile)

end subroutine read_data_1d_new
!.....................................................................

subroutine read_data_scalar_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain, tile_count, mosaicfile)

! this subroutine is for reading a single number
  character(len=*), intent(in)           :: filename, fieldname
  real, intent(inout)                    :: data     !zero dimension data 
  real, dimension(1,1,1)                 :: data_3d
  type(domain2d), intent(in), optional   :: domain
  integer, intent(in) , optional         :: timelevel
  logical, intent(in), optional          :: append_pelist_name, no_domain
  integer, intent(in), optional          :: tile_count
  character(len=*), intent(in), optional :: mosaicfile  

  if(present(no_domain)) then
     if(.NOT. no_domain) call mpp_error(FATAL, 'fms_io(read_data_scalar_new): no_domain should be true for field ' &
                                 //trim(fieldname)//' of file '//trim(filename) )
  end if

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain=.true., tile_count=tile_count, mosaicfile=mosaicfile)

  data = data_3d(1,1,1)

end subroutine read_data_scalar_new
!.....................................................................

function unique_axes(file, index, id_axes, siz_axes, dom)
  type(restart_file_type),   intent(inout)           :: file
  integer,                      intent(in)           :: index
  integer, dimension(:),       intent(out)           :: id_axes
  integer, dimension(:),       intent(out)           :: siz_axes
  type(domain1d), dimension(:), intent(in), optional :: dom
  integer                                            :: unique_axes
  type(var_type), pointer, save :: cur_var => NULL()
  integer :: i,j
  logical :: found
  
  unique_axes=0

  if(index <0 .OR. index > 3) call mpp_error(FATAL,"unique_axes(fms_io_mod): index should be 1, 2 or 3")

  do i = 1, file%nvar
     cur_var => file%var(i)
     if(cur_var%ndim < index) cycle
     found = .false.
     do j = 1, unique_axes
        if(siz_axes(j) == cur_var%gsiz(index) ) then
           if(PRESENT(dom)) then
              if(cur_var%domain_idx == id_axes(j) ) then
                 found = .true.
                 exit
              else if(cur_var%domain_idx >0 .AND. id_axes(j) >0) then
                 if(dom(cur_var%domain_idx) == dom(id_axes(j)) ) then
                    found = .true.
                    exit
                 end if
              end if
           else
              found = .true.
              exit
           end if
        end if  
     end do
     if(found) then
        cur_var%id_axes(index) = j
     else
        unique_axes = unique_axes+1
        if(unique_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') unique_axes, max_axes
           if(index == 1 ) then
              call mpp_error(FATAL,'# x axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else if(index == 2 ) then
              call mpp_error(FATAL,'# y axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else
              call mpp_error(FATAL,'# z axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           end if
        endif
        id_axes(unique_axes)   = cur_var%domain_idx
        siz_axes(unique_axes) = cur_var%gsiz(index)        
        cur_var%id_axes(index) = unique_axes
     end if
  end do       

  cur_var => NULL()
  
  return

end function unique_axes

  !#######################################################################
  !#######################################################################
  !   --------- routines for reading distributed data ---------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !
  ! reading can be done either by all PEs (default) or by only the root PE
  ! this is controlled by namelist variable "read_all_pe".
  
  ! By default, array data is expected to be declared in data domain and no_halo
  !is NOT needed, however IF data is decalared in COMPUTE domain then optional NO_HALO should be .true.

  !#######################################################################

subroutine read_data_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  real,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  real, dimension(isg:ieg,jsg:jeg)           :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"  
end subroutine read_data_2d

!#######################################################################

subroutine read_ldata_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  logical, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  logical, dimension(isg:ieg,jsg:jeg)        :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"
end subroutine read_ldata_2d
!#######################################################################

subroutine read_idata_2d ( unit, data, end)

  integer, intent(in)                        :: unit
  integer, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end
  integer, dimension(isg:ieg,jsg:jeg)        :: gdata
  integer                                    :: len
  logical                                    :: no_halo

  include "read_data_2d.inc"
end subroutine read_idata_2d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_2d ( unit, data, end)
  
  integer, intent(in)                           :: unit
  complex,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional                :: end
  complex, dimension(isg:ieg,jsg:jeg)           :: gdata
  integer                                       :: len
  logical                                       :: no_halo

  include "read_data_2d.inc"
end subroutine read_cdata_2d
#endif

!#######################################################################

subroutine read_data_3d ( unit, data, end)

  integer, intent(in)                           :: unit
  real,    intent(out), dimension(isd:,jsd:,:)  :: data
  logical, intent(out), optional                :: end  
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer                                       :: len
  logical                                       :: no_halo

  include "read_data_3d.inc"
end subroutine read_data_3d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_3d ( unit, data, end)

  integer, intent(in)                              :: unit
  complex, intent(out), dimension(isd:,jsd:,:)     :: data
  logical, intent(out), optional                   :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer                                          :: len
  logical                                          :: no_halo

  include "read_data_3d.inc"  
end subroutine read_cdata_3d
#endif

!#######################################################################

subroutine read_data_4d ( unit, data, end)

  integer, intent(in)                                        :: unit
  real,    intent(out), dimension(isd:,jsd:,:,:)             :: data
  logical, intent(out), optional                             :: end
  real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer                                                    :: len
  logical                                                    :: no_halo
! WARNING: memory usage with this routine could be costly
   
  include "read_data_4d.inc"  
end subroutine read_data_4d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine read_cdata_4d ( unit, data, end)

  integer, intent(in)                                           :: unit
  complex, intent(out), dimension(isd:,jsd:,:,:)                :: data
  logical, intent(out), optional                                :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer                                                       :: len
  logical                                                       :: no_halo
! WARNING: memory usage with this routine could be costly
 
  include "read_data_4d.inc"
end subroutine read_cdata_4d
#endif

!#######################################################################
!     -------- routines for writing distributed data --------
! before calling these routines the domain decompostion must be set
! by calling "set_domain" with the appropriate domain2d data type
!#######################################################################
subroutine write_data_2d ( unit, data )
  integer, intent(in)                       :: unit
  real,    intent(in), dimension(isd:,jsd:) :: data  
  real, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_data_2d

!#######################################################################

subroutine write_ldata_2d ( unit, data )

  integer, intent(in)                       :: unit
  logical, intent(in), dimension(isd:,jsd:) :: data
  logical, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_ldata_2d

!#######################################################################
subroutine write_idata_2d ( unit, data )

  integer, intent(in)                       :: unit
  integer, intent(in), dimension(isd:,jsd:) :: data
  integer, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_idata_2d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_2d ( unit, data )

  integer, intent(in)                       :: unit
  complex, intent(in), dimension(isd:,jsd:) :: data
  complex, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_cdata_2d
#endif

!#######################################################################

subroutine write_data_3d ( unit, data )

  integer, intent(in) :: unit
  real,    intent(in), dimension(isd:,jsd:,:) :: data
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
    
  include "write_data.inc"
end subroutine write_data_3d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_3d ( unit, data )

  integer, intent(in) :: unit
  complex, intent(in), dimension(isd:,jsd:,:) :: data
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

  include "write_data.inc"
end subroutine write_cdata_3d
#endif

!#######################################################################
subroutine write_data_4d ( unit, data )

  integer, intent(in) :: unit
  real,    intent(in), dimension(isd:,jsd:,:,:) :: data
  real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: n

  if (.not.associated(Current_domain))  &
       call mpp_error(FATAL,'fms_io(write_data_4d): need to call set_domain ')

! get the global data and write only on root pe
! do this one field at a time to save memory
  do n = 1, size(data,4)
     call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
  enddo
  if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata
end subroutine write_data_4d

!#######################################################################

#ifdef OVERLOAD_C8
subroutine write_cdata_4d ( unit, data )

  integer, intent(in) :: unit
  complex,    intent(in), dimension(isd:,jsd:,:,:) :: data
  complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: n

  if (.not.associated(Current_domain)) call mpp_error(FATAL, 'fms_io(write_cdata_4d): need to call set_domain')

! get the global data and write only on root pe
! do this one field at a time to save memory
  do n = 1, size(data,4)
     call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
  enddo
  if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata
end subroutine write_cdata_4d
#endif

!#######################################################################
! private routines (read_eof,do_read)
! this routine is called when an EOF is found while
! reading a distributed data file using read_data

subroutine read_eof (end_found)
  logical, intent(out), optional :: end_found
  
  if (present(end_found))then
     end_found = .true.
  else
     call mpp_error(FATAL,'fms_io(read_eof): unexpected EOF')
  endif
end subroutine read_eof

!#######################################################################
! determines if current pe should read data
! checks namelist variable read_all_pe

function do_read ( )
  logical :: do_read
  do_read = mpp_pe() == mpp_root_pe() .or. read_all_pe
end function do_read

!!#######################################################################

subroutine reset_field_name(fileObj, id_field, name)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  character(len=*),        intent(in)         :: name

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_name): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_name): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )

  fileObj%var(id_field)%name = trim(name)

end subroutine reset_field_name 

!#######################################################################

subroutine reset_field_pointer_r0d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real,                    intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r0d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p0dr(1, id_field)%p => data

end subroutine reset_field_pointer_r0d

!#######################################################################

subroutine reset_field_pointer_r1d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:),      intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r1d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p1dr(1, id_field)%p => data

end subroutine reset_field_pointer_r1d


!#######################################################################
subroutine reset_field_pointer_r2d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:),    intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p2dr(1, id_field)%p => data

end subroutine reset_field_pointer_r2d

!#######################################################################

subroutine reset_field_pointer_r3d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:,:),  intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r3d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p3dr(1, id_field)%p => data

end subroutine reset_field_pointer_r3d

!#######################################################################

subroutine reset_field_pointer_i0d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer,                 intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i0d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p0di(1, id_field)%p => data

end subroutine reset_field_pointer_i0d

!#######################################################################

subroutine reset_field_pointer_i1d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:),   intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i1d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p1di(1, id_field)%p => data

end subroutine reset_field_pointer_i1d


!#######################################################################
subroutine reset_field_pointer_i2d(fileObj, id_field, data)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:,:), intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p2di(1, id_field)%p => data

end subroutine reset_field_pointer_i2d

!#######################################################################

subroutine reset_field_pointer_i3d(fileObj, id_field, data)
  type(restart_file_type),   intent(inout)      :: fileObj
  integer,                   intent(in)         :: id_field
  integer, dimension(:,:,:), intent(in), target :: data   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i3d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 1) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d): one-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not one level" )

  fileObj%p3di(1, id_field)%p => data

end subroutine reset_field_pointer_i3d

!#######################################################################

subroutine reset_field_pointer_r0d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real,                    intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r0d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r0d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p0dr(1, id_field)%p => data1
  fileObj%p0dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r0d_2level

!#######################################################################

subroutine reset_field_pointer_r1d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:),      intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r1d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r1d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p1dr(1, id_field)%p => data1
  fileObj%p1dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r1d_2level

!#######################################################################

subroutine reset_field_pointer_r2d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:),    intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r2d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r2d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p2dr(1, id_field)%p => data1
  fileObj%p2dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r2d_2level

!#######################################################################

subroutine reset_field_pointer_r3d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  real, dimension(:,:,:),  intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_r3d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_r3d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p3dr(1, id_field)%p => data1
  fileObj%p3dr(2, id_field)%p => data2

end subroutine reset_field_pointer_r3d_2level

!#######################################################################

subroutine reset_field_pointer_i0d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer,                 intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i0d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i0d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p0di(1, id_field)%p => data1
  fileObj%p0di(2, id_field)%p => data2

end subroutine reset_field_pointer_i0d_2level

!#######################################################################

subroutine reset_field_pointer_i1d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:),   intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i1d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i1d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p1di(1, id_field)%p => data1
  fileObj%p1di(2, id_field)%p => data2

end subroutine reset_field_pointer_i1d_2level

!#######################################################################

subroutine reset_field_pointer_i2d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type), intent(inout)      :: fileObj
  integer,                 intent(in)         :: id_field
  integer, dimension(:,:), intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i2d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i2d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p2di(1, id_field)%p => data1
  fileObj%p2di(2, id_field)%p => data2

end subroutine reset_field_pointer_i2d_2level

!#######################################################################

subroutine reset_field_pointer_i3d_2level(fileObj, id_field, data1, data2)
  type(restart_file_type),   intent(inout)      :: fileObj
  integer,                   intent(in)         :: id_field
  integer, dimension(:,:,:), intent(in), target :: data1, data2   

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(reset_field_pointer_i3d_2level): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id_field < 0 .OR. id_field > fileObj%nvar) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d_2level): id_field should be positive integer and "// &
         "no larger than number of fields in the file "//trim(fileObj%name) )
  if(fileObj%var(id_field)%siz(4) .NE. 2) call mpp_error(FATAL, &
         "fms_io(reset_field_pointer_i3d_2level): two-level reset_field_pointer is called, but "//&
         "field "//trim(fileObj%var(id_field)%name)//" of file "//trim(fileObj%name)//" is not two level" )

  fileObj%p3di(1, id_field)%p => data1
  fileObj%p3di(2, id_field)%p => data2

end subroutine reset_field_pointer_i3d_2level

!#########################################################################
!   This function returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise. 
!
! Arguments: name - A pointer to the field that is being queried.
!  (in)  fileObj - The control structure returned by a previous call to
!                  register_restart_field
function query_initialized_name(fileObj, name)
  type(restart_file_type)      :: fileObj
  character(len=*), intent(in) :: name

  logical :: query_initialized_name

  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_name): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_name = .false.
  do m=1,fileObj%nvar
    if (trim(name) == fileObj%var(m)%name) then
      if (fileObj%var(m)%initialized) query_initialized_name = .true.
      exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (m<=fileObj%nvar) then
     fileObj%var(m)%initialized = .true.
  else if(mpp_pe() == mpp_root_pe()) then
    call mpp_error(NOTE,"fms_io(query_initialized_name): Unknown restart variable "//name// &
                        " queried for initialization.")
  end if

end function query_initialized_name


!   This function returns 1 if the field pointed to by f_ptr has
! initialized from a restart file, and 0 otherwise.  If f_ptr is
! NULL, it tests whether the entire restart file has been success-
! fully read.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
function query_initialized_r2d(fileObj, f_ptr, name)
  type(restart_file_type),      intent(in) :: fileObj 
  real, dimension(:,:), target, intent(in) :: f_ptr
  character(len=*),             intent(in) :: name

  logical :: query_initialized_r2d
  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_r2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_r2d = .false.
  do m=1, fileObj%nvar
     if (ASSOCIATED(fileObj%p2dr(1,m)%p,f_ptr)) then
        if (fileObj%var(m)%initialized) query_initialized_r2d = .true.
        exit
     endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if
  ! queried again.
  if (m<=fileObj%nvar) then
     fileObj%var(m)%initialized = .true.
  else
     query_initialized_r2d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() ) call mpp_error(NOTE, "fms_io(query_initialized_r2d): Unable to find "// &
          trim(name)//" queried by pointer, "//"probably because of the suspect comparison of pointers by ASSOCIATED.")
     query_initialized_r2d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() .AND. query_initialized_r2d) call mpp_error(NOTE, &
          "fms_io(query_initialized_r2d): "//trim(name)// " initialization confirmed by name.")
  endif

  return

end function query_initialized_r2d


!#######################################################################
!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action 
! open_namelist_file  MPP_ASCII   MPP_RDONLY  
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################


! <FUNCTION NAME="open_namelist_file">
!   <DESCRIPTION>
! Opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml".
!   </DESCRIPTION>
! <IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_namelist_file (file) result (unit)
  character(len=*), intent(in), optional :: file
  integer :: unit

  if (.not.module_is_initialized) call fms_io_init ( )    
  if (present(file)) then
     call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  else
     call mpp_open ( unit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  endif
end function open_namelist_file
! </FUNCTION>

! <FUNCTION NAME="open_restart_file">
!   <DESCRIPTION> 
! Opens single restart file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_restart_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit  
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---

  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error(FATAL,'fms_io(open_restart_file): action should be either read or write in file'//trim(file))
  end select
  
  call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
       access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

end function open_restart_file
! </FUNCTION>


! <FUNCTION NAME="open_direct_file">
!   <DESCRIPTION> 
! Opens single direct access file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>

  function open_direct_file (file, action, recl) result (unit)
    character(len=*), intent(in) :: file, action
    integer,          intent(in) :: recl
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    select case (lowercase(trim(action)))
    case ('read')
       mpp_action = MPP_RDONLY
    case ('write')
       mpp_action = MPP_OVERWR
    case default
       call mpp_error(FATAL,'invalid option for argument action')
    end select

    call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
         access=MPP_DIRECT, threading=MPP_SINGLE, nohdrs=.true., recl=recl )

  end function open_direct_file
! </FUNCTION>

! <FUNCTION NAME=" open_ieee32_file">
!   <DESCRIPTION>  
! Opens single 32-bit ieee file for reading by all PEs or 
! writing by root PE only (writing is not recommended)
! the file has no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_ieee32_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit  
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---
  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error (FATAL,'fms_io(open_ieee32_file): action should be either read or write in file'//trim(file))
  end select
  
  if (iospec_ieee32(1:1) == ' ') then
     call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
          nohdrs=.true. )
  else
     call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
          nohdrs=.true., iospec=iospec_ieee32 )
  endif
end function open_ieee32_file
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME=" close_file">
!   <DESCRIPTION>
!  Closes files that are opened by: open_namelist_file, open restart_file,
! and open_ieee32_file. Users should use mpp_close for other cases.
!   </DESCRIPTION>
!<IN NAME="unit" TYPE="integer">
! unit number of the file to be closed
! </IN>
!<IN NAME="status" TYPE="character, optional">
! action to be performed: can be 'delete'
! </IN>

subroutine close_file (unit, status)
  integer,          intent(in)           :: unit
  character(len=*), intent(in), optional :: status
  
  if (.not.module_is_initialized) call fms_io_init ( )
  if (unit == stdlog()) return    
  if (present(status)) then
     if (lowercase(trim(status)) == 'delete') then
        call mpp_close (unit, action=MPP_DELETE)
     else
        call mpp_error(FATAL,'fms_io(close_file): status should be DELETE')
     endif
  else
     call mpp_close (unit)
  endif
end subroutine close_file
! </FUNCTION>

!#######################################################################
  

! <SUBROUTINE NAME="set_domain">
!   <DESCRIPTION>
! set_domain is called to save the domain2d data type prior to
! calling the distributed data I/O routines, read_data and write_data.
!   </DESCRIPTION>
! <IN NAME="Domain2" TYPE="domain2D">
! domain to be passed to routines in fms_io_mod, Current_domain will point to
! this Domain2
! </IN>
subroutine set_domain (Domain2)
    
  type(domain2D), intent(in), target :: Domain2

  if (.NOT.module_is_initialized) call fms_io_init ( )

!  --- set_domain must be called before a read_data or write_data ---
  if (associated(Current_domain)) nullify (Current_domain)
  Current_domain => Domain2
  
  !  --- module indexing to shorten read/write routines ---
  
  call mpp_get_compute_domain (Current_domain,is ,ie ,js ,je )
  call mpp_get_data_domain    (Current_domain,isd,ied,jsd,jed)
  call mpp_get_global_domain  (Current_domain,isg,ieg,jsg,jeg)
end subroutine set_domain
!#######################################################################
! </SUBROUTINE>

! <SUBROUTINE NAME="nullify_domain">
subroutine nullify_domain ()
!   <DESCRIPTION>
! Use to nulify domain that has been assigned by set_domain.
!   </DESCRIPTION>
  if (.NOT.module_is_initialized) call fms_io_init ( )

!  --- set_domain must be called before a read_data or write_data ---

  if (associated(Current_domain)) nullify (Current_domain)
  is=0;ie=0;js=0;je=0
  isd=0;ied=0;jsd=0;jed=0
  isg=0;ieg=0;jsg=0;jeg=0
end subroutine nullify_domain
! </SUBROUTINE>

! <SUBROUTINE NAME="return_domain">
!   <DESCRIPTION>
! This routine is the reverse of set_domain above. This routine is called when 
! users want to retrieve the domain2d that is used in fms_io_mod
!   </DESCRIPTION>
! <OUT NAME="domain2" TYPE="domain2D">
! domain returned from  fms_io_mod.
! </OUT>
subroutine return_domain(domain2)
  type(domain2D), intent(inout) :: domain2

  if (associated(Current_domain)) then
     domain2 = Current_domain
  else
     domain2 = NULL_DOMAIN2D       
  endif
end subroutine return_domain
! </SUBROUTINE>

!#######################################################################
! this will be a private routine with the next release
! users should get the domain decomposition from the domain2d data type

!#######################################################################
! <SUBROUTINE NAME="get_domain_decomp">
!   <DESCRIPTION>
! This will be a private routine with the next release.
! Users should get the domain decomposition from the domain2d data type.
!   </DESCRIPTION>
! <OUT NAME="x" TYPE="integer">
! array containing beginning and ending indices of global and compute domain in x direction
! </OUT>
! <OUT NAME="y" TYPE="integer">
! array containing beginning and ending indices of global and compute domain in y direction
! </OUT>
subroutine get_domain_decomp ( x, y )

  integer, intent(out), dimension(4) :: x, y
  
  if (mpp_pe() == mpp_root_pe())  call mpp_error(NOTE, &
       'subroutine get_domain_decomp will be removed with the next release')
  x = (/ isg, ieg, is, ie /)
  y = (/ jsg, jeg, js, je /)

end subroutine get_domain_decomp
! </SUBROUTINE>

subroutine get_axis_cart(axis, cart)      

  type(axistype), intent(in) :: axis
  character(len=1), intent(out) :: cart
  character(len=1) :: axis_cart
  character(len=16), dimension(2) :: lon_names, lat_names
  character(len=16), dimension(3) :: z_names
  character(len=16), dimension(2) :: t_names
  character(len=16), dimension(2) :: lon_units, lat_units
  character(len=8) , dimension(4) :: z_units
  character(len=3) , dimension(4) :: t_units
  character(len=32) :: name
  integer :: i

  lon_names = (/'lon','x  '/)
  lat_names = (/'lat','y  '/)
  z_names = (/'depth ','height','z     '/)
  t_names = (/'time','t   '/)
  lon_units = (/'degrees_e   ', 'degrees_east'/)
  lat_units = (/'degrees_n    ', 'degrees_north'/)
  z_units = (/'cm ','m  ','pa ','hpa'/)
  t_units = (/'sec', 'min','hou','day'/)  
  call mpp_get_atts(axis,cartesian=axis_cart)
  cart = 'N'  
  if (axis_cart == 'x' ) cart = 'X'
  if (axis_cart == 'y' ) cart = 'Y'
  if (axis_cart == 'z' ) cart = 'Z'
  if (axis_cart == 't' ) cart = 'T'
  if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
     call mpp_get_atts(axis,name=name)
     name = lowercase(name)
     do i=1,size(lon_names(:))
        if (lowercase(name(1:3)) == trim(lon_names(i))) cart = 'X'
     enddo
     do i=1,size(lat_names(:))
        if (name(1:3) == trim(lat_names(i))) cart = 'Y'
     enddo
     do i=1,size(z_names(:))
        if (name == trim(z_names(i))) cart = 'Z'
     enddo
     do i=1,size(t_names(:))
        if (name(1:3) == t_names(i)) cart = 'T'
     enddo
  end if

  if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
     call mpp_get_atts(axis,units=name)
     name = lowercase(name)
     do i=1,size(lon_units(:))
        if (trim(name) == trim(lon_units(i))) cart = 'X'
     enddo
     do i=1,size(lat_units(:))
        if (trim(name) == trim(lat_units(i))) cart = 'Y'
     enddo
     do i=1,size(z_units(:))
        if (trim(name) == trim(z_units(i))) cart = 'Z'
     enddo
     do i=1,size(t_units(:))
        if (name(1:3) == trim(t_units(i))) cart = 'T'
     enddo
  end if
  
  return
end subroutine get_axis_cart


! The following function is here as a last resort.
! This is copied from what was utilities_mod in order that redundant code 
! could be deleted.

 function open_file ( file, form, action, access, threading, recl ) &
             result ( unit )

 character(len=*), intent(in) :: file 
 character(len=*), intent(in), optional :: form, action, access, threading
 integer         , intent(in), optional :: recl 
 integer  :: unit 

 character(len=32) :: form_local, action_local, access_local, thread_local
 character(len=32) :: action_ieee32
 logical :: open, no_headers, do_ieee32
 integer :: mpp_format, mpp_action, mpp_access, mpp_thread
!-----------------------------------------------------------------------

   if ( .not. module_is_initialized ) then
        call fms_io_init ( )
!        do_init = .false.
   endif

!   ---- return stdlog if this is the logfile ----

    if (trim(file) == 'logfile.out') then
       unit = stdlog()
       return
    endif

!   ---- is this file open and connected to a unit ?? ---- 

   inquire (file=trim(file), opened=open, number=unit)

!  cannot open a file that is already open
!  except for the log file

   if ( open .and. unit >= 0 ) then
      call mpp_error (FATAL, 'open_file in fms_mod : '// &
                       'file '//trim(file)//' is already open')
   endif   

!  --- defaults ---

   form_local   = 'formatted';  if (present(form))      form_local   = form
   access_local = 'sequential'; if (present(access))    access_local = access
   thread_local = 'single';     if (present(threading)) thread_local = threading
   no_headers   = .true.
   do_ieee32    = .false.

   if (present(action)) then    ! must be present
      action_local = action
   else
      call mpp_error (FATAL, 'open_file in fms_mod : argument action not present')
   endif


!   --- file format ---

    select case (lowercase(trim(form_local)))
       case ('formatted')
           mpp_format = MPP_ASCII
       case ('ascii')
           mpp_format = MPP_ASCII
       case ('unformatted')
           mpp_format = MPP_NATIVE
       case ('native')
           mpp_format = MPP_NATIVE
       case ('ieee32')
           do_ieee32 = .true.
       case ('netcdf')
           mpp_format = MPP_NETCDF
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument form')
    end select

!   --- action (read,write,append) ---

    select case (lowercase(trim(action_local)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case ('append')
           mpp_action = MPP_APPEND
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument action')
    end select

!   --- file access (sequential,direct) ---

    select case (lowercase(trim(access_local)))
       case ('sequential')
           mpp_access = MPP_SEQUENTIAL
       case ('direct')
           mpp_access = MPP_DIRECT
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument access')
    end select

!   --- threading (single,multi) ---

    select case (lowercase(trim(thread_local)))
       case ('single')
           mpp_thread = MPP_SINGLE
       case ('multi')
           mpp_thread = MPP_MULTI
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument thread')
           if (trim(file) /= '_read_error.nml') no_headers = .false.
    end select

!   ---------------- open file -----------------------

    if ( .not.do_ieee32 ) then
       call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                       access=mpp_access, threading=mpp_thread,        &
                       nohdrs=no_headers, recl=recl )
    else
     ! special open for ieee32 file
     ! fms_mod has iospec value
     ! pass local action flag to open changing append to write
       action_ieee32 = action_local
       if (lowercase(trim(action_ieee32)) == 'append') action_ieee32 = 'write'
       unit = open_ieee32_file ( file, action_ieee32 )
    endif

!-----------------------------------------------------------------------

 end function open_file

  !#######################################################################

  function string_from_integer(n)
    integer, intent(in) :: n
    character(len=16) :: string_from_integer

    if(n<0) then
       call mpp_error(FATAL, 'fms_io_mod: n should be non-negative integer, contact developer')
    else if( n<10 ) then
       write(string_from_integer,'(i1)') n
    else if( n<100 ) then
       write(string_from_integer,'(i2)') n
    else if( n<1000 ) then
       write(string_from_integer,'(i3)') n
    else if( n<10000 ) then
       write(string_from_integer,'(i4)') n
    else if( n<100000 ) then
       write(string_from_integer,'(i5)') n
    else if( n<1000000 ) then
       write(string_from_integer,'(i6)') n
    else if( n<10000000 ) then
       write(string_from_integer,'(i7)') n
    else if( n<100000000 ) then
       write(string_from_integer,'(i8)') n
    else
       call mpp_error(FATAL, 'fms_io_mod: n is too big, contact developer')
    end if

    return

  end function string_from_integer

  !#######################################################################
  function string_from_real(a)
    real, intent(in) :: a
    character(len=32) :: string_from_real

    write(string_from_real,*) a

    return

  end function string_from_real

  !#######################################################################

 subroutine get_tile_string(str_out, str_in, tile, str2_in)
    character(len=*), intent(inout)        :: str_out
    character(len=*), intent(in)           :: str_in
    integer,          intent(in)           :: tile
    character(len=*), intent(in), optional :: str2_in

    if(tile > 0 .AND. tile < 9) then
       write(str_out,'(a,i1)') trim(str_in), tile
    else if(tile >= 10 .AND. tile < 99) then
       write(str_out,'(a,i2)') trim(str_in), tile
    else
       call mpp_error(FATAL, "FMS_IO: get_tile_string: tile must be a positive number less than 100")
    end if

    if(present(str2_in)) str_out=trim(str_out)//trim(str2_in)

 end subroutine get_tile_string


  !#####################################################################
  subroutine get_mosaic_tile_file(file_in, file_out, is_no_domain, domain, tile_count, mosaicfile)
    character(len=*), intent(in)                   :: file_in
    character(len=*), intent(out)                  :: file_out
    logical,          intent(in)                   :: is_no_domain
    type(domain2D),   intent(in), optional, target :: domain
    character(len=*), intent(in), optional         :: mosaicfile
    integer,          intent(in), optional         :: tile_count 
    character(len=256)                             :: basefile, tilename
    integer                                        :: lens, ntiles, ntileMe, tile
    integer, dimension(:), allocatable             :: tile_id
    type(domain2d), pointer, save                  :: d_ptr =>NULL()


    if(index(file_in, '.nc', back=.true.)==0) then
       basefile = trim(file_in)
    else
       lens = len_trim(file_in)
       if(file_in(lens-2:lens) .NE. '.nc') call mpp_error(FATAL, &
            'fms_io_mod: .nc should be at the end of file '//trim(file_in))
       basefile = file_in(1:lens-3)
    end if

    if(mpp_mosaic_defined())then
       !--- get the tile name
       ntiles = 1
       if(PRESENT(domain))then
          ntiles = mpp_get_ntile_count(domain)
          d_ptr => domain
       elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
          ntiles = mpp_get_ntile_count(Current_domain)
          d_ptr => Current_domain
       endif
       if(ntiles > 1 )then
          ntileMe = mpp_get_current_ntile(d_ptr)
          allocate(tile_id(ntileMe))
          tile_id = mpp_get_tile_id(d_ptr)
          tile = 1
          if(present(tile_count)) tile = tile_count
          if(present(mosaicfile)) then
             !--- read tilename from mosaic file
             call read_data(mosaicfile, "gridtiles", tilename, level=tile)
          else
             tilename = 'tile'//string(tile_id(tile))
          end if
          deallocate(tile_id)
          if(index(basefile,'.'//trim(tilename),back=.true.) == 0)then
             basefile = trim(basefile)//'.'//trim(tilename);
          end if
       end if
    endif

    file_out = trim(basefile)//'.nc'

    d_ptr =>NULL()

  end subroutine get_mosaic_tile_file

  !#############################################################################
  subroutine get_mosaic_tile_grid(grid_file, mosaic_file, domain, tile_count) 
    character(len=*), intent(out)          :: grid_file
    character(len=*), intent(in)           :: mosaic_file
    type(domain2D),   intent(in)           :: domain
    integer,          intent(in), optional :: tile_count 
    integer                                :: tile, ntileMe
    integer, dimension(:), allocatable     :: tile_id

    tile = 1
    if(present(tile_count)) tile = tile_count
    ntileMe = mpp_get_current_ntile(domain)
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(domain)     
    call read_data(mosaic_file, "gridfiles", grid_file, level=tile_id(tile) )
    grid_file = 'INPUT/'//trim(grid_file)
    deallocate(tile_id)

  end subroutine get_mosaic_tile_grid

  !#############################################################################
  ! return false if the attribute is not find in the file.
  function get_global_att_value_text(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    character(len=*), intent(inout) :: attvalue
    logical                         :: get_global_att_value_text
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_text = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = trim(mpp_get_att_char(global_atts(i)))
          get_global_att_value_text = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_text

  !#############################################################################
  ! return false if the attribute is not find in the file.
  function get_global_att_value_real(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    real,             intent(inout) :: attvalue
    logical                         :: get_global_att_value_real
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_real = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = mpp_get_att_real_scalar(global_atts(i))
          get_global_att_value_real = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_real

  !#############################################################################
  ! This routine will get the actual file name, as well as if read_dist is true or false.
  subroutine get_file_name(orig_file, actual_file, read_dist, is_no_domain, domain, &
                          append_pelist_name, tile_count, mosaicfile  )
    character(len=*),                 intent(in) :: orig_file
    character(len=*),                intent(out) :: actual_file
    logical,                         intent(out) :: read_dist
    logical,                          intent(in) :: is_no_domain
    type(domain2D), target, optional, intent(in) :: domain
    logical,                optional, intent(in) :: append_pelist_name  
    integer,                optional, intent(in) :: tile_count  
    character(len=*),       optional, intent(in) :: mosaicfile

    logical             :: fexist


    fexist = .false.
    read_dist=.false.

    !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count, 
    !JWD:  I just don't see how the physics can remain "tile neutral"
    !z1l:  one solution is add one more public interface called set_tile_count
    call get_mosaic_tile_file(orig_file, actual_file, is_no_domain, domain, tile_count, mosaicfile)

    if (PRESENT(append_pelist_name)) then
       if (append_pelist_name) then
          if(is_no_domain) call mpp_error(FATAL, &
               'fms_io(get_file_name): when append_pelist_name is set to true, no_domain should not be set to true') 
          actual_file = trim(actual_file)//trim(pelist_name)
       endif
    endif

    inquire (file=trim(actual_file)//trim(pe_name), exist=fexist)
    if(.not. fexist) inquire (file=trim(actual_file)//'.nc'//trim(pe_name), exist=fexist)
    if(fexist) then
       read_dist = .true.  
    else
       read_dist = .false.
       inquire (file=trim(actual_file), exist=fexist)
       if(.not. fexist) inquire (file=trim(actual_file)//'.nc', exist=fexist)
    endif
    if( .not. fexist) then 
       call mpp_error(FATAL,'fms_io(get_file_name): file '//trim(orig_file)//' and distributed file not found')
    end if

  end subroutine get_file_name
  

  !#############################################################################
  subroutine get_file_unit(filename, unit, index_file, read_dist )
    character(len=*),                  intent(in) :: filename
    integer,                          intent(out) :: unit, index_file 
    logical,                           intent(in) :: read_dist

    logical  :: file_opened
    integer  :: i

    ! Need to check if filename has been opened or not
    file_opened=.false.
    do i=1,num_files_r
       if (files_read(i)%name == trim(filename))  then
          index_file = i
          unit = files_read(index_file)%unit
          return 
       endif
    enddo

    ! need to open the file now
    ! Increase num_files_r and set file_type 
    if(num_files_r == max_files_r) &  ! need to have bigger max_files_r
         call mpp_error(FATAL,'fms_io(get_file_unit): max_files_r exceeded, increase it via fms_io_nml')
    num_files_r=num_files_r + 1 
    if (read_dist .and. thread_r == MPP_SINGLE) then
       call mpp_error(FATAL,'fms_io(get_file_unit): single-threaded read from distributed fileset not allowed' &
            //'change threading_read to MULTI')
    endif
    if(read_dist) then
       call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=thread_r, &
            fileset=MPP_MULTI)
    else
       call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=thread_r, &
            fileset=MPP_SINGLE)
    end if
    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var (max_fields) )
    index_file = num_files_r
    files_read(index_file)%unit = unit   

  end subroutine get_file_unit

  !#############################################################################
  subroutine get_field_id(unit, index_file, fieldname, index_field, is_no_domain, is_not_dim)
    integer,          intent(in) :: unit
    integer,          intent(in) :: index_file
    character(len=*), intent(in) :: fieldname
    integer,         intent(out) :: index_field
    logical,          intent(in) :: is_no_domain
    logical,          intent(in) :: is_not_dim

    character(len=128)                     :: name
    type(axistype),  dimension(max_axes)   :: axes
    type(fieldtype), dimension(max_fields) :: fields
    integer                                :: i, j, ndim, nvar, natt, var_dim
    integer                                :: siz_in(4)

    index_field = -1
    do j = 1, files_read(index_file)%nvar
       if (trim(files_read(index_file)%var(j)%name) == trim(fieldname)) then
          index_field = j
          return
       endif
    enddo

    !--- fieldname is not read, so need to get fieldname from file 
    files_read(index_file)%nvar = files_read(index_file)%nvar + 1
    if(files_read(index_file)%nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar, max_fields 
       call  mpp_error(FATAL,'fms_io(get_field_id): max_fields exceeded, needs increasing, nvar/max_fields=' &
            //trim(error_msg))
    endif
    call mpp_get_info(unit, ndim, nvar, natt, files_read(index_file)%max_ntime)
    if(files_read(index_file)%max_ntime < 1)  files_read(index_file)%max_ntime = 1
    if(nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar,max_fields
       call mpp_error(FATAL,'fms_io(get_field_id): max_fields too small needs increasing,nvar/max_fields=' &
            //trim(error_msg)//'in file'//trim(files_read(index_file)%name))
    endif
    call mpp_get_fields(unit, fields(1:nvar))     
    siz_in = 1
    index_field = files_read(index_file)%nvar
    files_read(index_file)%var(index_field)%is_dimvar = .false.

    do i=1, nvar
       call mpp_get_atts(fields(i),name=name,ndim=var_dim,siz=siz_in)
       if (lowercase(trim(name)) == lowercase(trim(fieldname))) then ! found the variable
          if(var_dim .lt.3) then
             do j=var_dim+1,3
                siz_in(j)=1
             enddo
          endif
          files_read(index_file)%var(index_field)%name    = fieldname
          files_read(index_file)%var(index_field)%field   = fields(i)
          files_read(index_file)%var(index_field)%siz(:)  = siz_in
          files_read(index_file)%var(index_field)%gsiz(:) = siz_in
          return
       endif
    enddo

    !--- the fieldname may be a dimension variable.
    if( .not. is_not_dim) then
       if (ndim > max_axes) then
          write(error_msg,'(I3,"/",I3)') ndim, max_axes 
          call  mpp_error(FATAL,'fms_io(get_field_id): max_axes exceeded, needs increasing, ndim/max_fields=' &
               //trim(error_msg)//' in file '//trim(files_read(index_file)%name))             
       endif
       call mpp_get_axes(unit, axes(1:ndim))
       do i=1,ndim
          call mpp_get_atts(axes(i), name=name, len = siz_in(1)) 
          if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
             if(.not. is_no_domain) call mpp_error(FATAL, &
                  'fms_io(get_field_id): the field is a dimension variable, no_domain should be true.')
             files_read(index_file)%var(index_field)%is_dimvar = .true.
             files_read(index_file)%var(index_field)%name      = fieldname
             files_read(index_file)%var(index_field)%axis      = axes(i)
             files_read(index_file)%var(index_field)%siz(:)    = siz_in
             files_read(index_file)%var(index_field)%gsiz(:)   = siz_in
             return
          endif
       enddo
    end if
    !--- the field is not in the file when reaching here.
    call mpp_error(FATAL, 'fms_io(get_field_id): field '//trim(fieldname)// &
                   ' NOT found in file '//trim(files_read(index_file)%name))

  end subroutine get_field_id

!#######################################################################
! check the existence of the given file name
! if the file_name string has zero length or the
! first character is blank return a false result
! <FUNCTION NAME="file_exist">

!   <OVERVIEW>
!     Checks the existence of a given file name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Checks the existence of the given file name.
!     If the file_name string has zero length or the
!     first character is blank return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     file_exist ( file_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If file_name exists the result 
!     is true, otherwise false is returned.
!     If the length of character string "file_name" is zero or the first
!     character is blank, then the returned value will be false.
!     When reading a file, this function is often used in conjunction with
!     routine open_file.
!   </OUT>
!   <ERROR MSG="set_domain not called" STATUS="FATAL">
!     Before calling write_data you must first call set_domain with domain2d data 
!     type associated with the distributed data you are writing.
!   </ERROR>

 function file_exist (file_name, domain)
  character(len=*), intent(in) :: file_name
  type(domain2d), intent(in), optional :: domain
  logical  file_exist
  character(len=5)                      :: pe_name
  character(len=256)                    :: fname

   file_exist = .false.
   if (len_trim(file_name) == 0) return
   if (file_name(1:1) == ' ')    return

   if(present(domain)) then
      call get_mosaic_tile_file(file_name, fname, .false.,  domain)
   else
      fname = file_name
   end if

   inquire (file=trim(fname), exist=file_exist)
   !--- also check to see if it is distributed data
   if(.not. file_exist) then
      write(pe_name,'(a,i4.4)' )'.', mpp_pe()    
     inquire (file=trim(fname)//trim(pe_name), exist=file_exist)
   end if


 end function file_exist
! </FUNCTION>


!#######################################################################
! <FUNCTION NAME="field_exist">

!   <OVERVIEW>
!     check if a given field name exists in a given file name. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     check if a given field name exists in a given file name. 
!     If the field_name string has zero length or the
!     first character is blank return a false result.
!     if the file file_name don't exist, return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     field_exist ( file_name, field_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <IN NAME="field_name"  TYPE="character" >
!     A field name that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If field exists in the 
!     file file_name, the result is true, otherwise false is returned.
!     If the length of character string "field_name" is zero or the first
!     character is blank, then the returned value will be false.
!     if the file file_name don't exist, return a false result.
!   </OUT>

 function field_exist (file_name, field_name)
  character(len=*), intent(in) :: file_name
  character(len=*), intent(in) :: field_name
  logical                      :: field_exist
  integer                      :: unit, ndim, nvar, natt, ntime, i
  character(len=64)            :: name
  type(fieldtype), allocatable :: fields(:)
  character(len=5)             :: pe_name
  logical                      :: is_exist

   field_exist = .false.
   if (len_trim(field_name) == 0) return
   if (field_name(1:1) == ' ')    return


   write(pe_name,'(a,i4.4)' )'.', mpp_pe()   
   inquire (file=trim(file_name)//trim(pe_name), exist=is_exist)
   if(.not. is_exist) inquire (file=trim(file_name)//'.nc'//trim(pe_name), exist=is_exist)
   if(is_exist) then
      call mpp_open(unit, trim(file_name), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, &
           fileset = MPP_MULTI)
   else
      inquire (file=trim(file_name), exist=is_exist)
      if(.not. is_exist) inquire (file=trim(file_name)//'.nc', exist=is_exist)
      if(is_exist) then
         !--- open the file file_name
         call mpp_open(unit, trim(file_name), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, &
              fileset = MPP_SINGLE)
      end if
   end if
   if(.not. is_exist) return

    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit,fields)

    do i=1, nvar
       call mpp_get_atts(fields(i),name=name)
       if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
    enddo

    deallocate(fields)
    call mpp_close(unit)

    return

 end function field_exist
! </FUNCTION>

end module fms_io_mod


#ifdef test_fms_io

 program fms_io_test

 use mpp_mod,         only: mpp_pe, mpp_npes, mpp_root_pe, mpp_init, mpp_exit
 use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE, mpp_chksum
 use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_mosaic
 use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
 use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
 use mpp_domains_mod, only: mpp_domains_set_stack_size
 use mpp_io_mod,      only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
 use fms_io_mod,      only: read_data, write_data, fms_io_init, fms_io_exit
 use fms_io_mod,      only: file_exist, register_restart_field, save_restart, restore_state
 use fms_io_mod,      only: restart_file_type
 use mpp_io_mod,      only: MAX_FILE_SIZE

 implicit none

 integer :: sizex_latlon_grid = 144
 integer :: sizey_latlon_grid = 90
 integer :: size_cubic_grid = 48
 integer :: nz = 10, nt = 2, halo = 1
 integer :: stackmax =4000000
 integer :: num_step = 4 ! number of time steps to run, this is used for intermediate run.
                         ! set num_step = 0 for no intermediate run.
 logical :: do_write=.true. ! set this to false for high resolution and single file,
                            ! split file capability is not implemented for write_data

 namelist /test_fms_io_nml/ sizex_latlon_grid, sizey_latlon_grid, size_cubic_grid, &
                            nz, nt, halo, num_step, stackmax, do_write

 integer           :: unit, io_status, step
 character(len=20) :: time_stamp

 type data_storage_type
    real,    allocatable, dimension(:,:,:,:) :: data1_r3d, data2_r3d, data1_r3d_read, data2_r3d_read
    real,    allocatable, dimension(:,:,:)   :: data1_r2d, data2_r2d, data1_r2d_read, data2_r2d_read
    real,    allocatable, dimension(:,:)     :: data1_r1d, data2_r1d, data1_r1d_read, data2_r1d_read
    real,    allocatable, dimension(:)       :: data1_r0d, data2_r0d, data1_r0d_read, data2_r0d_read
    integer, allocatable, dimension(:,:,:,:) :: data1_i3d, data2_i3d, data1_i3d_read, data2_i3d_read
    integer, allocatable, dimension(:,:,:)   :: data1_i2d, data2_i2d, data1_i2d_read, data2_i2d_read
    integer, allocatable, dimension(:,:)     :: data1_i1d, data2_i1d, data1_i1d_read, data2_i1d_read
    integer, allocatable, dimension(:)       :: data1_i0d, data2_i0d, data1_i0d_read, data2_i0d_read
 end type data_storage_type
 
 type(data_storage_type), save :: latlon_data
 type(data_storage_type), save :: cubic_data
 type(domain2d),          save :: domain_latlon
 type(domain2d),          save :: domain_cubic
 type(restart_file_type), save :: restart_latlon
 type(restart_file_type), save :: restart_cubic
 integer                       :: ntile_latlon = 1
 integer                       :: ntile_cubic = 6
 integer                       :: npes

 character(len=128) :: file_latlon, file_cubic

 call mpp_init
 npes = mpp_npes()

 call mpp_domains_init  

 call fms_io_init

 if (file_exist('input.nml') )then
    call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
    read(unit,test_fms_io_nml,iostat=io_status)

    if (io_status > 0) then
     call mpp_error(FATAL,'=>test_fms_io: Error reading test_fms_io_nml')
  endif
    call mpp_close (unit)
 end if



 write(stdout(), test_fms_io_nml )
  call mpp_domains_set_stack_size(stackmax)
 !--- currently we assume at most two time level will be written to restart file.
 if(nt > 2) call mpp_error(FATAL, "test_fms_io: test_fms_io_nml variable nt should be no larger than 2")

 file_latlon   = "test.res.latlon_grid.save_restore.nc"
 file_cubic    = "test.res.cubic_grid.save_restore.nc"

 call setup_test_restart(restart_latlon, "latlon_grid", ntile_latlon, latlon_data, file_latlon, domain_latlon)
 call setup_test_restart(restart_cubic,  "cubic_grid", ntile_cubic, cubic_data, file_cubic, domain_cubic )

 if(file_exist('INPUT/'//trim(file_latlon), domain_latlon)) then
    call restore_state(restart_latlon)
    call compare_restart("latlon_grid save_restore", latlon_data)
 end if
 if(file_exist('INPUT/'//trim(file_cubic), domain_cubic) ) then
    call restore_state(restart_cubic)
    call compare_restart("cubic_grid save_restore", cubic_data)
 end if
 
 !---copy data
 if(mod(npes,ntile_latlon) == 0) call copy_restart_data(latlon_data)
 if(mod(npes,ntile_cubic) == 0 ) call copy_restart_data(cubic_data)

 do step = 1, num_step
    write(time_stamp, '(a,I4.4)') "step", step
    call save_restart(restart_latlon, time_stamp)
    call save_restart(restart_cubic, time_stamp)
 end do
 call save_restart(restart_latlon)
 call save_restart(restart_cubic)

 if(mod(npes,ntile_latlon) == 0) call release_storage_memory(latlon_data)
 if(mod(npes,ntile_cubic) == 0 ) call release_storage_memory(cubic_data)

 if(mod(npes,ntile_cubic) == 0 ) call mpp_error(NOTE, "test_fms_io: restart test is done for latlon_grid")
 if(mod(npes,ntile_cubic) == 0 ) call mpp_error(NOTE, "test_fms_io: restart test is done for cubic_grid")

 call fms_io_exit
 call mpp_domains_exit
 call mpp_exit

contains

  !******************************************************************************
  subroutine setup_test_restart(restart_data, type, ntiles, storage, file, domain)
    type(restart_file_type),   intent(inout) :: restart_data
    character(len=*), intent(in)             :: type
    integer,          intent(in)             :: ntiles
    type(data_storage_type), intent(inout)   :: storage
    character(len=*), intent(in)             :: file  
    type(domain2d),   intent(inout)          :: domain
    character(len=128)                       :: file_r
    character(len=128)                       :: file_w
    integer                                  :: pe, npes_per_tile, tile
    integer                                  :: num_contact
    integer                                  :: n, layout(2)
    integer, allocatable, dimension(:,:)     :: global_indices, layout2D
    integer, allocatable, dimension(:)       :: pe_start, pe_end
    integer, dimension(1)                    :: tile1, tile2
    integer, dimension(1)                    :: istart1, iend1, jstart1, jend1
    integer, dimension(1)                    :: istart2, iend2, jstart2, jend2
    integer                                  :: i, j, k, nx, ny
    integer                                  :: isc, iec, jsc, jec
    integer                                  :: isd, ied, jsd, jed


    file_r = "INPUT/test.res."//trim(type)//".read_write.nc"
    file_w = "RESTART/test.res."//trim(type)//".read_write.nc"

    select case(type)
    case("latlon_grid")
       nx = sizex_latlon_grid
       ny = sizey_latlon_grid
    case("cubic_grid")
       nx = size_cubic_grid
       ny = size_cubic_grid
    case default
       call mpp_error(FATAL, "test_fms_io: "//type//" is not a valid option")
    end select

    pe   = mpp_pe()
    if(mod(npes,ntiles) .NE. 0) then
       call mpp_error(NOTE, "test_fms_io: npes can not be divided by ntiles, no test will be done for "//trim(type))
       return
    end if
    npes_per_tile = npes/ntiles
    tile = pe/npes_per_tile + 1

    call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
    allocate(global_indices(4,ntiles), layout2D(2,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)       = layout
       pe_start(n)         = (n-1)*npes_per_tile
       pe_end(n)           = n*npes_per_tile-1
    end do
    num_contact = 0
    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                           istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                           pe_start, pe_end, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, name = type  )
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    allocate(storage%data1_r3d(isd:ied, jsd:jed, nz, nt), storage%data1_r3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data2_r3d(isd:ied, jsd:jed, nz, nt), storage%data2_r3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data1_i3d(isd:ied, jsd:jed, nz, nt), storage%data1_i3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data2_i3d(isd:ied, jsd:jed, nz, nt), storage%data2_i3d_read(isd:ied, jsd:jed, nz, nt) )
    allocate(storage%data1_r2d(isd:ied, jsd:jed,     nt), storage%data1_r2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data2_r2d(isd:ied, jsd:jed,     nt), storage%data2_r2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data1_i2d(isd:ied, jsd:jed,     nt), storage%data1_i2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data2_i2d(isd:ied, jsd:jed,     nt), storage%data2_i2d_read(isd:ied, jsd:jed,     nt) )
    allocate(storage%data1_r1d(                  nz, nt), storage%data1_r1d_read(                  nz, nt) )
    allocate(storage%data2_r1d(                  nz, nt), storage%data2_r1d_read(                  nz, nt) )
    allocate(storage%data1_i1d(                  nz, nt), storage%data1_i1d_read(                  nz, nt) )
    allocate(storage%data2_i1d(                  nz, nt), storage%data2_i1d_read(                  nz, nt) )
    allocate(storage%data1_r0d(                      nt), storage%data1_r0d_read(                      nt) )
    allocate(storage%data2_r0d(                      nt), storage%data2_r0d_read(                      nt) )
    allocate(storage%data1_i0d(                      nt), storage%data1_i0d_read(                      nt) )
    allocate(storage%data2_i0d(                      nt), storage%data2_i0d_read(                      nt) )

    storage%data1_r3d = 0; storage%data1_r3d_read = 0; storage%data2_r3d = 0; storage%data2_r3d_read = 0
    storage%data1_i3d = 0; storage%data1_i3d_read = 0; storage%data2_i3d = 0; storage%data2_i3d_read = 0
    storage%data1_r2d = 0; storage%data1_r2d_read = 0; storage%data2_r2d = 0; storage%data2_r2d_read = 0
    storage%data1_i2d = 0; storage%data1_i2d_read = 0; storage%data2_i2d = 0; storage%data2_i2d_read = 0
    storage%data1_r1d = 0; storage%data1_r1d_read = 0; storage%data2_r1d = 0; storage%data2_r1d_read = 0
    storage%data1_i1d = 0; storage%data1_i1d_read = 0; storage%data2_i1d = 0; storage%data2_i1d_read = 0
    storage%data1_r0d = 0; storage%data1_r0d_read = 0; storage%data2_r0d = 0; storage%data2_r0d_read = 0
    storage%data1_i0d = 0; storage%data1_i0d_read = 0; storage%data2_i0d = 0; storage%data2_i0d_read = 0
    do n = 1, nt
       storage%data1_r0d(n) =  tile + n*1e-3
       storage%data2_r0d(n) = -tile - n*1e-3
       storage%data1_i0d(n) =  tile*1e3 + n
       storage%data2_i0d(n) = -tile*1e3 - n
       do k = 1, nz
          storage%data1_r1d(k,n) =   tile*1e3 + n + k*1e-3
          storage%data2_r1d(k,n) =  -tile*1e3 - n - k*1e-3
          storage%data1_i1d(k,n) =   tile*1e6 + n*1e3 + k
          storage%data2_i1d(k,n) =  -tile*1e6 - n*1e3 - k
          do j = jsc, jec
             do i = isc, iec
                storage%data1_r3d(i,j,k,n) =  tile*1e6 + n*1e3 + k + i*1e-3 + j*1e-6; 
                storage%data2_r3d(i,j,k,n) = -tile*1e6 - n*1e3 - k - i*1e-3 - j*1e-6; 
                storage%data1_i3d(i,j,k,n) =  tile*1e9 + n*1e8 + k*1e6 + i*1e3 + j; 
                storage%data2_i3d(i,j,k,n) = -tile*1e9 - n*1e8 - k*1e6 - i*1e3 - j; 
             end do
          end do
       end do

       do j = jsc, jec
          do i = isc, iec
             storage%data1_r2d(i,j,n) =  tile*1e1 + n + i*1e-3 + j*1e-6; 
             storage%data2_r2d(i,j,n) = -tile*1e1 - n - i*1e-3 - j*1e-6; 
             storage%data1_i2d(i,j,n) =  tile*1e7 + n*1e6 + i*1e3 + j; 
             storage%data2_i2d(i,j,n) = -tile*1e7 - n*1e6 - i*1e3 - j; 
          end do
       end do
    end do
    if(file_exist(file_r, domain)) then
       do n = 1, nt
          call read_data(file_r, "data1_r3d", storage%data1_r3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data2_r3d", storage%data2_r3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data1_i3d", storage%data1_i3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data2_i3d", storage%data2_i3d_read(:,:,:,n), domain, timelevel = n )
          call read_data(file_r, "data1_r2d", storage%data1_r2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data2_r2d", storage%data2_r2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data1_i2d", storage%data1_i2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data2_i2d", storage%data2_i2d_read(:,:,  n), domain, timelevel = n )
          call read_data(file_r, "data1_r1d", storage%data1_r1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data2_r1d", storage%data2_r1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data1_i1d", storage%data1_i1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data2_i1d", storage%data2_i1d_read(:,    n), domain, timelevel = n )
          call read_data(file_r, "data1_r0d", storage%data1_r0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data2_r0d", storage%data2_r0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data1_i0d", storage%data1_i0d_read(      n), domain, timelevel = n )
          call read_data(file_r, "data2_i0d", storage%data2_i0d_read(      n), domain, timelevel = n )
       end do
       call compare_restart(type//" read_write", storage)
    end if


    !--- high resolution restart is not implemented for write data
    if(do_write ) then 
       do n = 1, nt
          call write_data(file_w, "data1_r3d", storage%data1_r3d(:,:,:,n), domain )
          call write_data(file_w, "data2_r3d", storage%data2_r3d(:,:,:,n), domain )
          call write_data(file_w, "data1_i3d", storage%data1_i3d(:,:,:,n), domain )
          call write_data(file_w, "data2_i3d", storage%data2_i3d(:,:,:,n), domain )
          call write_data(file_w, "data1_r2d", storage%data1_r2d(:,:,  n), domain )
          call write_data(file_w, "data2_r2d", storage%data2_r2d(:,:,  n), domain )
          call write_data(file_w, "data1_i2d", storage%data1_i2d(:,:,  n), domain )
          call write_data(file_w, "data2_i2d", storage%data2_i2d(:,:,  n), domain )
          call write_data(file_w, "data1_r1d", storage%data1_r1d(:,    n), domain )
          call write_data(file_w, "data2_r1d", storage%data2_r1d(:,    n), domain )
          call write_data(file_w, "data1_i1d", storage%data1_i1d(:,    n), domain )
          call write_data(file_w, "data2_i1d", storage%data2_i1d(:,    n), domain )
          call write_data(file_w, "data1_r0d", storage%data1_r0d(      n), domain )
          call write_data(file_w, "data2_r0d", storage%data2_r0d(      n), domain )
          call write_data(file_w, "data1_i0d", storage%data1_i0d(      n), domain )
          call write_data(file_w, "data2_i0d", storage%data2_i0d(      n), domain )
       end do
    end if

    !--- test register_restart_field, save_restart, restore_state

    call register_restart_field(restart_data, file, "data1_r3d", storage%data1_r3d_read(:,:,:,1), &
                                domain, longname="first data_r3d",units="none")
    call register_restart_field(restart_data, file, "data1_r3d", storage%data1_r3d_read(:,:,:,2), &
                                domain, longname="first data_r3d",units="none")
    call register_restart_field(restart_data, file, "data2_r3d", storage%data2_r3d_read(:,:,:,1), storage%data2_r3d_read(:,:,:,2), &
                                domain, longname="second data_i3d", units="none")

    call register_restart_field(restart_data, file, "data1_i3d", storage%data1_i3d_read(:,:,:,1), &
                                domain, longname="first data_i3d",units="none")
    call register_restart_field(restart_data, file, "data1_i3d", storage%data1_i3d_read(:,:,:,2), &
                                domain, longname="first data_i3d",units="none")
    call register_restart_field(restart_data, file, "data2_i3d", storage%data2_i3d_read(:,:,:,1), storage%data2_i3d_read(:,:,:,2), &
                                domain, longname="second data_i3d", units="none")

    call register_restart_field(restart_data, file, "data1_r2d", storage%data1_r2d_read(:,:,  1), &
                                domain, longname="first data_r2d",units="none")
    call register_restart_field(restart_data, file, "data1_r2d", storage%data1_r2d_read(:,:,  2), &
                                domain, longname="first data_r2d",units="none")
    call register_restart_field(restart_data, file, "data2_r2d", storage%data2_r2d_read(:,:,  1), storage%data2_r2d_read(:,:,2), &
                                domain, longname="second data_i2d", units="none")

    call register_restart_field(restart_data, file, "data1_i2d", storage%data1_i2d_read(:,:,  1), &
                                domain, longname="first data_i2d",units="none")
    call register_restart_field(restart_data, file, "data1_i2d", storage%data1_i2d_read(:,:,  2), &
                                domain, longname="first data_i2d",units="none")
    call register_restart_field(restart_data, file, "data2_i2d", storage%data2_i2d_read(:,:,  1), storage%data2_i2d_read(:,:,2), &
                                domain, longname="second data_i2d", units="none")

    call register_restart_field(restart_data, file, "data1_r1d", storage%data1_r1d_read(:,    1), &
                                domain, longname="first data_r1d",units="none")
    call register_restart_field(restart_data, file, "data1_r1d", storage%data1_r1d_read(:,    2), &
                                domain, longname="first data_r1d",units="none")
    call register_restart_field(restart_data, file, "data2_r1d", storage%data2_r1d_read(:,    1), storage%data2_r1d_read(:,  2), &
                                domain, longname="second data_i1d", units="none")

    call register_restart_field(restart_data, file, "data1_i1d", storage%data1_i1d_read(:,    1), &
                                domain, longname="first data_i1d",units="none")
    call register_restart_field(restart_data, file, "data1_i1d", storage%data1_i1d_read(:,    2), &
                                domain, longname="first data_i1d",units="none")
    call register_restart_field(restart_data, file, "data2_i1d", storage%data2_i1d_read(:,    1), storage%data2_i1d_read(:,  2), &
                                domain, longname="second data_i1d", units="none")


    call register_restart_field(restart_data, file, "data1_r0d", storage%data1_r0d_read(      1), &
                                domain, longname="first data_r0d",units="none")
    call register_restart_field(restart_data, file, "data1_r0d", storage%data1_r0d_read(      2), &
                                domain, longname="first data_r0d",units="none")
    call register_restart_field(restart_data, file, "data2_r0d", storage%data2_r0d_read(      1), storage%data2_r0d_read(    2), &
                                domain, longname="second data_i0d", units="none")

    call register_restart_field(restart_data, file, "data1_i0d", storage%data1_i0d_read(      1), &
                                domain, longname="first data_i0d",units="none")
    call register_restart_field(restart_data, file, "data1_i0d", storage%data1_i0d_read(      2), &
                                domain, longname="first data_i0d",units="none")
    call register_restart_field(restart_data, file, "data2_i0d", storage%data2_i0d_read(      1), storage%data2_i0d_read(    2), &
                                domain, longname="second data_i0d", units="none")

  end subroutine setup_test_restart

  subroutine compare_restart(type, storage)
    character(len=*), intent(in)             :: type
    type(data_storage_type), intent(inout)   :: storage

       call compare_data_r4d(storage%data1_r3d, storage%data1_r3d_read, type//" data1_r3d")
       call compare_data_r4d(storage%data2_r3d, storage%data2_r3d_read, type//" data2_r3d")
       call compare_data_i4d(storage%data1_i3d, storage%data1_i3d_read, type//" data1_i3d")
       call compare_data_i4d(storage%data2_i3d, storage%data2_i3d_read, type//" data2_i3d")
       call compare_data_r3d(storage%data1_r2d, storage%data1_r2d_read, type//" data1_r2d")
       call compare_data_r3d(storage%data2_r2d, storage%data2_r2d_read, type//" data2_r2d")
       call compare_data_i3d(storage%data1_i2d, storage%data1_i2d_read, type//" data1_i2d")
       call compare_data_i3d(storage%data2_i2d, storage%data2_i2d_read, type//" data2_i2d")
       call compare_data_r2d(storage%data1_r1d, storage%data1_r1d_read, type//" data1_r1d")
       call compare_data_r2d(storage%data2_r1d, storage%data2_r1d_read, type//" data2_r1d")
       call compare_data_i2d(storage%data1_i1d, storage%data1_i1d_read, type//" data1_i1d")
       call compare_data_i2d(storage%data2_i1d, storage%data2_i1d_read, type//" data2_i1d")
       call compare_data_r1d(storage%data1_r0d, storage%data1_r0d_read, type//" data1_r0d")
       call compare_data_r1d(storage%data2_r0d, storage%data2_r0d_read, type//" data2_r0d")
       call compare_data_i1d(storage%data1_i0d, storage%data1_i0d_read, type//" data1_i0d")
       call compare_data_i1d(storage%data2_i0d, storage%data2_i0d_read, type//" data2_i0d")

  end subroutine compare_restart

  subroutine release_storage_memory(storage)
    type(data_storage_type), intent(inout)   :: storage

    deallocate(storage%data1_r3d, storage%data2_r3d, storage%data1_r3d_read, storage%data2_r3d_read)
    deallocate(storage%data1_i3d, storage%data2_i3d, storage%data1_i3d_read, storage%data2_i3d_read)
    deallocate(storage%data1_r2d, storage%data2_r2d, storage%data1_r2d_read, storage%data2_r2d_read)
    deallocate(storage%data1_i2d, storage%data2_i2d, storage%data1_i2d_read, storage%data2_i2d_read)
    deallocate(storage%data1_r1d, storage%data2_r1d, storage%data1_r1d_read, storage%data2_r1d_read)
    deallocate(storage%data1_i1d, storage%data2_i1d, storage%data1_i1d_read, storage%data2_i1d_read)
    deallocate(storage%data1_r0d, storage%data2_r0d, storage%data1_r0d_read, storage%data2_r0d_read)
    deallocate(storage%data1_i0d, storage%data2_i0d, storage%data1_i0d_read, storage%data2_i0d_read)

  end subroutine release_storage_memory

  subroutine copy_restart_data(storage)
    type(data_storage_type), intent(inout)   :: storage

    storage%data1_r3d_read = storage%data1_r3d; storage%data2_r3d_read = storage%data2_r3d
    storage%data1_i3d_read = storage%data1_i3d; storage%data2_i3d_read = storage%data2_i3d
    storage%data1_r2d_read = storage%data1_r2d; storage%data2_r2d_read = storage%data2_r2d
    storage%data1_i2d_read = storage%data1_i2d; storage%data2_i2d_read = storage%data2_i2d
    storage%data1_r1d_read = storage%data1_r1d; storage%data2_r1d_read = storage%data2_r1d
    storage%data1_i1d_read = storage%data1_i1d; storage%data2_i1d_read = storage%data2_i1d
    storage%data1_r0d_read = storage%data1_r0d; storage%data2_r0d_read = storage%data2_r0d
    storage%data1_i0d_read = storage%data1_i0d; storage%data2_i0d_read = storage%data2_i0d

    return

  end subroutine copy_restart_data

  subroutine compare_data_r4d( a, b, string )
    real, intent(in), dimension(:,:,:,:) :: a, b
    character(len=*), intent(in)         :: string
    integer(LONG_KIND)                   :: sum1, sum2
    integer                              :: i, j, k, l
    integer, parameter                   :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) .or. size(a,4) .ne. size(b,4) ) &
         call mpp_error(FATAL,'compare_data_r4d: size of a and b does not match')

    do l = 1, size(a,4)
       do k = 1, size(a,3)
          do j = 1, size(a,2)
             do i = 1, size(a,1)
                if(a(i,j,k,l) .ne. b(i,j,k,l)) then
                   write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,f18.9,a,f18.9)')" at pe ", mpp_pe(), &
                        ", at point (",i,", ", j, ", ", k, ", ", l, "), a = ", a(i,j,k,l), ", b = ", b(i,j,k,l)
                   call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
                endif
             enddo
          enddo
       enddo
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r4d

  subroutine compare_data_i4d( a, b, string )
    integer, intent(in), dimension(:,:,:,:) :: a, b
    character(len=*), intent(in)            :: string
    real                                    :: real_a(size(a,1),size(a,2),size(a,3),size(a,4))
    real                                    :: real_b(size(b,1),size(b,2),size(b,3),size(b,4))

    real_a = a 
    real_b = b
    call compare_data_r4d(real_a, real_b, string)

  end subroutine compare_data_i4d


  subroutine compare_data_r3d( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in)       :: string
    integer(LONG_KIND)                 :: sum1, sum2
    integer                            :: i, j, l
    integer, parameter                 :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_data_r3d: size of a and b does not match')

    do l = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,l) .ne. b(i,j,l)) then
                write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", l, "), a = ", a(i,j,l), ", b = ", b(i,j,l)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r3d

  subroutine compare_data_i3d( a, b, string )
    integer, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in)          :: string
    real                                  :: real_a(size(a,1),size(a,2),size(a,3))
    real                                  :: real_b(size(b,1),size(b,2),size(b,3))

    real_a = a 
    real_b = b
    call compare_data_r3d(real_a, real_b, string)

  end subroutine compare_data_i3d


  subroutine compare_data_r2d( a, b, string )
    real, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in)     :: string
    integer(LONG_KIND)               :: sum1, sum2
    integer                          :: i, l
    integer, parameter               :: stdunit = 6

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
         call mpp_error(FATAL,'compare_data_r2d: size of a and b does not match')

    do l = 1, size(a,2)
       do i = 1, size(a,1)
          if(a(i,l) .ne. b(i,l)) then
             write(stdunit,'(a,i3,a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
                  ", at point (",i, ", ", l, "), a = ", a(i,l), ", b = ", b(i,l)
             call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
          endif
       enddo
    end do
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r2d

  subroutine compare_data_i2d( a, b, string )
    integer, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in)        :: string
    real                                :: real_a(size(a,1),size(a,2))
    real                                :: real_b(size(b,1),size(b,2))

    real_a = a 
    real_b = b
    call compare_data_r2d(real_a, real_b, string)

  end subroutine compare_data_i2d

  subroutine compare_data_r1d( a, b, string )
    real, intent(in), dimension(:) :: a, b
    character(len=*), intent(in)   :: string
    integer(LONG_KIND)             :: sum1, sum2
    integer                        :: l
    integer, parameter             :: stdunit = 6

    if(size(a,1) .ne. size(b,1) ) &
         call mpp_error(FATAL,'compare_data_r1d: size of a and b does not match')

    do l = 1, size(a(:))
       if(a(l) .ne. b(l)) then
          write(stdunit,'(a,i3,a,i3,a,f16.9,a,f16.9)')" at pe ", mpp_pe(), &
               ", at point (",l, "), a = ", a(l), ", b = ", b(l)
          call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
       endif
    enddo
    sum1 = mpp_chksum( a, (/mpp_pe()/) )
    sum2 = mpp_chksum( b, (/mpp_pe()/) )

    if( sum1.EQ.sum2 )then
       if( mpp_pe() .EQ. mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
       !--- in some case, even though checksum agree, the two arrays 
       !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
       !--- hence we need to check the value point by point.
    else
       call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_data_r1d

  subroutine compare_data_i1d( a, b, string )
    integer, intent(in), dimension(:) :: a, b
    character(len=*), intent(in)      :: string
    real                              :: real_a(size(a(:)))
    real                              :: real_b(size(b(:)))

    real_a = a 
    real_b = b
    call compare_data_r1d(real_a, real_b, string)

  end subroutine compare_data_i1d

end program fms_io_test
#endif
