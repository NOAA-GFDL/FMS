#include <fms_platform.h>

module fms_io_mod

!
!
! <CONTACT EMAIL="Giang.Nong@noaa.gov">
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
! <PRE>
!threading_read='multi', fileset_read='single', threading_write='multi', fileset_write='multi' (default)
!threading_read='multi', fileset_read='single', threading_write='single', fileset_write='single'
! </PRE>
!</DESCRIPTION>
! <NAMELIST NAME="fms_io_nml">  
! <DATA NAME="threading_read" TYPE="character">
! threading_read can be 'single' or 'multi'
! </DATA>
! <DATA NAME="threading_write" TYPE="character">
! threading_write can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fileset_read" TYPE="character">
! fileset_read can be 'single' or 'multi'
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
!</NAMELIST>
  
use mpp_io_mod, only : mpp_open, mpp_close, mpp_io_init, mpp_io_exit, &
     MPP_NETCDF, MPP_ASCII, MPP_MULTI, MPP_SINGLE, &
     mpp_read, mpp_write, mpp_write_meta, &
     mpp_get_info, mpp_get_atts, MPP_IEEE32, &
     MPP_OVERWR, fieldtype, axistype, atttype, &
     MPP_RDONLY, MPP_NATIVE, MPP_DELETE, MPP_APPEND, &
     default_field, default_axis, default_att, &
     mpp_get_fields, MPP_SEQUENTIAL, MPP_DIRECT, mpp_get_axes, &
     mpp_get_axis_data
use mpp_domains_mod, only : domain2d, domain1d, mpp_get_domain_components, &
     mpp_get_compute_domain, mpp_get_data_domain, &
     mpp_get_global_domain, NULL_DOMAIN1D, &
     NULL_DOMAIN2D, mpp_global_field, operator( == )
  
use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_npes, &
     stdlog, stdout, mpp_broadcast, ALL_PES, &
     mpp_chksum, mpp_sync, mpp_get_current_pelist, mpp_npes, lowercase
implicit none
private

integer                     :: max_files_w, max_files_r    !  max_files will change as needed
integer, parameter, private :: max_fields=150
integer, parameter, private :: max_axes=40
integer, parameter, private :: max_atts=20
integer, parameter, private :: max_domains = 10
type buff_type
   real, dimension(:,:,:), _ALLOCATABLE :: buffer _NULL
end type buff_type

type file_type
   integer                                :: unit ! mpp_io unit for netcdf file
   integer                                :: ndim, nvar, natt, max_ntime
   logical                                :: domain_present(max_fields)
   character(len=128)                     :: filename
   integer                                :: siz(max_fields,4)  ! X/Y/Z/T extent of fields (data domain 
!size for distributed writes;global size for reads)
   integer                                :: gsiz(max_fields,4) ! global X/Y/Z/T extent of fields
   integer                                :: unit_tmpfile(max_fields)
   character(len=128)                     :: fieldname(max_fields)
   type(buff_type), dimension(:), _ALLOCATABLE :: field_buffer _NULL
   type(fieldtype), dimension(max_fields) :: fields
   type(axistype),  dimension(max_axes)   :: axes   ! spatial axes
   type(atttype),  dimension(max_atts)    :: atts
   integer                                :: domain_idx(max_fields) 
   logical                                :: is_dimvar(max_fields)
end type file_type

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
   module procedure read_cdata_2d,read_cdata_3d,read_cdata_4d
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
   module procedure write_cdata_2d,write_cdata_3d,write_cdata_4d
end interface

integer, private :: num_files_r=0 ! number of currently opened files for reading
integer, private :: num_files_w=0 ! number of currently opened files for writing
integer, private :: num_domains = 0 ! number of domains in array_domain
logical, private :: fms_netcdf_override = .true., fms_netcdf_restart = .true.
character(len=32), private :: threading_read, fileset_read, threading_write, &
     fileset_write, format ! global i/o settings
integer, private :: thread_r, thread_w, fset_r, fset_w, form
logical, private :: module_is_initialized = .FALSE.
logical, private :: read_all_pe = .TRUE.
integer, allocatable, dimension(:) :: pelist
character(len=32) :: pelist_name
character(len=128):: error_msg  
character(len=64) :: iospec_ieee32 = '-N ieee_32'
  
!------ private data, pointer to current 2d domain ------
! entrained from fms_mod.  This will be deprecated in the future.
type(domain2D), pointer, private :: Current_domain =>NULL()

integer, private :: is,ie,js,je      ! compute domain
integer, private :: isd,ied,jsd,jed  ! data domain
integer, private :: isg,ieg,jsg,jeg  ! global domain

type(file_type), dimension(:), allocatable            :: files_read
type(file_type), dimension(:), allocatable            :: files_write
type(file_type), save                                 :: default_file
type(domain2d), dimension(max_domains), private, save :: array_domain
public  :: read_data, write_data, fms_io_init, fms_io_exit, field_size
public  :: open_namelist_file, open_restart_file, open_ieee32_file, close_file 
public  :: set_domain, nullify_domain, get_domain_decomp, return_domain
public  :: open_file, open_direct_file
public  :: get_restart_io_mode
private :: lookup_field_w, lookup_axis, unique_axes

character(len=128) :: version = '$Id: fms_io.F90,v 12.0 2005/04/14 17:56:29 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'

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

! <SUBROUTINE NAME="fms_io_init">
!   <DESCRIPTION>
! Initialize fms_io module
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_init()
!   </TEMPLATE>
subroutine fms_io_init()
    
  IMPLICIT NONE
    
  integer  :: i,j, unit, io_status
  logical :: file_exist
  namelist /fms_io_nml/ fms_netcdf_override, fms_netcdf_restart, &
       threading_read, fileset_read, threading_write, &
       fileset_write, format, read_all_pe, iospec_ieee32,max_files_w,max_files_r

  call mpp_io_init()
  if (module_is_initialized) return
! Initialize values of max_files_w and max_files_r  
  max_files_w = 40; max_files_r = 40
  threading_read='multi';fileset_read='single';format='netcdf'
  threading_write='multi';fileset_write='multi'

  call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
  read(unit,fms_io_nml,iostat=io_status)
  write(stdlog(), fms_io_nml)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>fms_io_init: Error reading input.nml')
  endif
  call mpp_close (unit)

! take namelist options if present

  select case (fileset_read) 
  case ('multi')
     fset_r = MPP_MULTI
  case ('single')
     fset_r = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: fileset_read should be multi/single but you chose'//trim(fileset_read))
  end select

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
  default_file%unit = -1
  default_file%ndim = 0  
  default_file%nvar = 0
  default_file%max_ntime = -1
  default_file%natt = 0
  default_file%siz(:,:) = 0
  default_file%gsiz(:,:) = 0
  default_file%unit_tmpfile(:) = -1
  default_file%filename = 'none'
  default_file%fieldname = 'none' 
  default_file%domain_present(:)=.false.
  default_file%domain_idx(:)= -1
  default_file%fields(:) = default_field
  default_file%axes(:) = default_axis
  default_file%atts(:) = default_att      
  default_file%is_dimvar(:) = .false.
! Initially allocate  files_write and files_read
  allocate(files_write(max_files_w),files_read(max_files_r))
  files_write(:) = default_file
  files_read(:)  = default_file

  do i = 1, max_domains
     array_domain(i) = NULL_DOMAIN2D
  enddo
  !---- initialize module domain2d pointer ----
  nullify (Current_domain)
  module_is_initialized = .TRUE.
  write (stdlog(),'(/,80("="),/(a))') trim(version), trim(tagname)
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
  use platform_mod, only: r8_kind
  IMPLICIT NONE
  integer, parameter :: max_axis_size=10000
  integer :: i,j,k,unit,unit2,index_field
  integer :: num_axes, t_axis_id, x_axis_id, y_axis_id, z_axis_id
  integer, dimension(4) :: size_field, global_size ! x/y/z/t
  integer :: siz_x_axes(max_axes), siz_y_axes(max_axes), siz_z_axes(max_axes), max_t_size
  integer :: x_axes(max_axes), y_axes(max_axes), z_axes(max_axes)
  integer :: num_x_axes, num_y_axes, num_z_axes
  type(domain1d) :: domain_x(max_fields), domain_y(max_fields), x_domains(max_axes), y_domains(max_axes)
  real, dimension(max_axis_size) :: axisdata
  real(r8_kind) :: tlev  
  character (len=128) :: axisname,filename, fieldname,temp_name
  type(domain2D) :: domain
  integer :: domain_idx
  if( .NOT.module_is_initialized )return !make sure it's only called once per PE
  do i=1,max_axis_size
     axisdata(i) = i
  enddo

! each field has an associated domain type (may be undefined).
! each file only needs to write unique axes (i.e. if 2 fields share an identical axis, then only write the axis once)
! unique axes are defined by the global size and domain decomposition (i.e. can support identical axis sizes with
! different domain decomposition)
  
  do i = 1, num_files_w
! determine maximum number of time levels for this file
     files_write(i)%max_ntime = 0
     do j=1,files_write(i)%nvar
        files_write(i)%max_ntime = max(files_write(i)%max_ntime,files_write(i)%siz(j,4))
     enddo
     allocate(files_write(i)%field_buffer(files_write(i)%nvar))
     filename = files_write(i)%filename
     call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=thread_w,&
          fileset=fset_w)
     do j = 1, max_fields
        if (files_write(i)%domain_present(j)) then
           domain_idx = files_write(i)%domain_idx(j)
           call mpp_get_domain_components(array_domain(domain_idx), domain_x(j), domain_y(j))
        else
           domain_x(j) = NULL_DOMAIN1D
           domain_y(j) = NULL_DOMAIN1D
        endif
     enddo
     siz_x_axes = -1;siz_y_axes=-1;siz_z_axes=-1
     x_axes = -1; y_axes = -1; z_axes=-1
     x_domains(:) = NULL_DOMAIN1D; y_domains(:) = NULL_DOMAIN1D
     x_axes = unique_axes(files_write(i)%gsiz(:,1),domain_x(:))
     do j=1,max_axes
        if (x_axes(j) > 0) then
           siz_x_axes(j) = files_write(i)%gsiz(x_axes(j),1) ! global array sizes
           x_domains(j) = domain_x(x_axes(j))
        endif
     end do
     y_axes = unique_axes(files_write(i)%gsiz(:,2),domain_y(:))
     do j=1,max_axes
        if (y_axes(j) > 0) then
           siz_y_axes(j) = files_write(i)%gsiz(y_axes(j),2) ! global array sizes
           y_domains(j) = domain_y(y_axes(j))
        endif
     end do
     z_axes = unique_axes(files_write(i)%gsiz(:,3))
     do j=1,max_axes
        if (z_axes(j) > 0) then
           siz_z_axes(j) = files_write(i)%gsiz(z_axes(j),3) ! global array sizes
        endif
     end do     
     num_axes=0
     num_x_axes=0
     j=1
     do while (x_axes(j) > 0) 
        if (j < 10) then
           write(axisname,'(a,i1)') 'xaxis_',j
        else
           write(axisname,'(a,i2)') 'xaxis_',j
        endif
        num_axes=num_axes+1
        if (files_write(i)%domain_present(x_axes(j))) then
           call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),domain=domain_x(x_axes(j)),cartesian='X')
        else
!    if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
           call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                data=axisdata(1:siz_x_axes(j)),cartesian='X')
        endif
        j=j+1
        num_x_axes=num_x_axes+1
        if (num_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') num_axes, max_axes
           call mpp_error(FATAL,'# x axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
        endif
     enddo
     num_y_axes=0
     j=1
     do while (y_axes(j) > 0) 
        if (j < 10) then
           write(axisname,'(a,i1)') 'yaxis_',j
        else
           write(axisname,'(a,i2)') 'yaxis_',j
        endif
        num_axes=num_axes+1
        if (files_write(i)%domain_present(y_axes(j))) then
           call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),domain=domain_y(y_axes(j)),cartesian='Y')
        else
!       if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
           call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                data=axisdata(1:siz_y_axes(j)),cartesian='Y')
        endif
        j=j+1
        num_y_axes=num_y_axes+1
        if (num_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') num_axes, max_axes
           call mpp_error(FATAL,'# y axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
        endif
     enddo     
     num_z_axes=0
     j=1
     do while (z_axes(j) > 0) 
        if (j < 10) then
           write(axisname,'(a,i1)') 'zaxis_',j
        else
           write(axisname,'(a,i2)') 'zaxis_',j
        endif
        num_axes=num_axes+1
!       if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
        call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
             data=axisdata(1:siz_z_axes(j)),cartesian='Z')
        j=j+1
        num_z_axes=num_z_axes+1
        if (num_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') num_axes, max_axes
           call mpp_error(FATAL,'# z axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
        endif
     enddo

! write time axis  (comment out if no time axis)
!       if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
     call mpp_write_meta(unit,files_write(i)%axes(num_axes+1),&
          'Time','time level','Time',cartesian='T')
     t_axis_id = num_axes+1

! write metadata for fields
     do j = 1, files_write(i)%nvar
        size_field = files_write(i)%gsiz(j,:) 
        x_axis_id = lookup_axis(siz_x_axes,size_field(1),x_domains,domain_x(j))            
        y_axis_id = lookup_axis(siz_y_axes,size_field(2),y_domains,domain_y(j))+num_x_axes            
        z_axis_id = lookup_axis(siz_z_axes,size_field(3))+num_x_axes+num_y_axes            
!          x_axis_id = x_axes(j)
!          y_axis_id = y_axes(j)+num_x_axes
!          z_axis_id = z_axes(j)+num_x_axes+num_y_axes
!          if((mpp_pe()==mpp_root_pe().and.thread_w == MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
        call mpp_write_meta(unit,files_write(i)%fields(j), &
             (/files_write(i)%axes(x_axis_id),&
             files_write(i)%axes(y_axis_id),&
             files_write(i)%axes(z_axis_id),&
             files_write(i)%axes(t_axis_id)/),files_write(i)%fieldname(j),&
             'none',files_write(i)%fieldname(j),pack=1)
     enddo

! write values for ndim of spatial axes
     do j = 1, num_axes
!          if((mpp_pe()==mpp_root_pe().and.thread_w == MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
        call mpp_write(unit,files_write(i)%axes(j))
     enddo

! retrieve and write data of each field
     do k = 1, files_write(i)%max_ntime
        do j = 1, files_write(i)%nvar
           size_field = files_write(i)%siz(j,:) ! local size here
           global_size = files_write(i)%gsiz(j,:) ! global size here
           unit2 = files_write(i)%unit_tmpfile(j)
           if (k == 1) then
              temp_name = trim(files_write(i)%filename)//'_'//&
                   trim(files_write(i)%fieldname(j))//'_tmp'
              call mpp_close(unit2)
              call mpp_open(unit2,temp_name,form=MPP_NATIVE,nohdrs=.true.,threading=MPP_MULTI, &
                   fileset=MPP_MULTI, action=MPP_RDONLY) 
              if(thread_w.eq.MPP_SINGLE .and. mpp_pe() == mpp_root_pe()) then
                 allocate(files_write(i)%field_buffer(j)%buffer(global_size(1),global_size(2),&
                      global_size(3)))
              else
                 allocate(files_write(i)%field_buffer(j)%buffer(size_field(1),size_field(2),&
                      size_field(3)))
              endif
              files_write(i)%unit_tmpfile(j) = unit2
           endif
           if ( k <= size_field(4)) then
              if((mpp_pe() == mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
                   read(unit2) files_write(i)%field_buffer(j)%buffer
           else
              files_write(i)%field_buffer(j)%buffer = 0.0
           endif           
           tlev=k
           if(files_write(i)%domain_present(j)) then
!                if (num_x_axes > 1 .or. num_y_axes > 1) call mpp_error(FATAL,&
!                     'restart data need to be on same grid when domain flag present')
              domain=array_domain(files_write(i)%domain_idx(j))
              if (mpp_pe() == mpp_root_pe().and. thread_w==MPP_SINGLE) then
                 call mpp_write(unit,files_write(i)%fields(j),&
                      files_write(i)%field_buffer(j)%buffer,tlev)
              else if (thread_w == MPP_MULTI) then
                 call mpp_write(unit,files_write(i)%fields(j),domain,&
                      files_write(i)%field_buffer(j)%buffer,tlev)
              endif
           else
              if (thread_w == MPP_MULTI .or. ((mpp_pe() == mpp_root_pe()).and.thread_w == MPP_SINGLE)) then
                 call mpp_write(unit,files_write(i)%fields(j),&
                      files_write(i)%field_buffer(j)%buffer,tlev)
              endif
           endif
        enddo ! end j loop
     enddo ! end k loop
     call mpp_close(unit)
     do j = 1, files_write(i)%nvar
        call mpp_close(files_write(i)%unit_tmpfile(j), action = MPP_DELETE)
     enddo
     deallocate(files_write(i)%field_buffer)
  enddo ! end i loop
!  call mpp_io_exit() !don't call it here... this may not be the last I/O
  module_is_initialized = .false.
  num_files_w = 0
  num_files_r = 0    
end subroutine fms_io_exit
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
subroutine write_data_i3d_new(filename, fieldname, data, domain,append_pelist_name, no_domain)
  IMPLICIT NONE

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: append_pelist_name, no_domain

  call write_data_3d_new(filename, fieldname, real(data), domain,append_pelist_name, no_domain)

end subroutine write_data_i3d_new
subroutine write_data_i2d_new(filename, fieldname, data, domain,append_pelist_name, no_domain)
  IMPLICIT NONE

  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: append_pelist_name, no_domain

  call write_data_2d_new(filename, fieldname, real(data), domain,append_pelist_name, no_domain)

end subroutine write_data_i2d_new
subroutine write_data_i1d_new(filename, fieldname, data, domain, append_pelist_name, no_domain)
  IMPLICIT NONE
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, dimension(:), intent(in) :: data
  logical, intent(in), optional :: append_pelist_name, no_domain

  call write_data_1d_new(filename, fieldname, real(data), domain, append_pelist_name, no_domain)

end subroutine write_data_i1d_new
subroutine write_data_iscalar_new(filename, fieldname, data, domain, append_pelist_name, no_domain)
  IMPLICIT NONE
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  integer, intent(in) :: data
  logical, intent(in), optional :: append_pelist_name, no_domain

  call write_data_scalar_new(filename, fieldname, real(data), domain, append_pelist_name, no_domain)

end subroutine write_data_iscalar_new
!=================================================================================
subroutine write_data_3d_new(filename, fieldname, data, domain,append_pelist_name, no_domain)
  IMPLICIT NONE

  character(len=*), intent(in) :: filename, fieldname 
  real, dimension(:,:,:), intent(in) :: data
  type(domain2d), intent(in), optional :: domain
  real, dimension(:,:,:), pointer ::global_data =>NULL()
  logical, intent(in), optional :: append_pelist_name, no_domain   
  character(len=128) :: temp_name     ! temp_name: name of the temporary file
  integer :: i, domain_idx
  integer :: nfile  ! index of the currently open file in array files
  integer :: index_field ! position of the fieldname in the list of fields
  integer :: unit2 ! unit of temporary file
  integer :: gxsize, gysize
  logical :: file_open = .false., is_no_domain = .false.
  character(len=256) :: fname  
  type(file_type), dimension(:), pointer  :: new_files_write

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_3d_new): need to call fms_io_init')  
  fname = trim(filename)

  if (PRESENT(append_pelist_name)) then
     if (append_pelist_name) then
        allocate(pelist(mpp_npes()))        
        call mpp_get_current_pelist(pelist,pelist_name)
        fname = trim(filename)//trim(pelist_name)
        deallocate(pelist)
     endif
  endif

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     if(PRESENT(domain) .AND. no_domain) &
       call mpp_error(FATAL, 'fms_io(write_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')
     is_no_domain = no_domain
  endif

! Check if filename has been open  or not
  if(num_files_w == 0) then
     file_open=.false.
  else
     file_open=.false.
     do i=1,num_files_w
        nfile = i
        if (trim(files_write(i)%filename) == trim(fname)) then
           file_open = .true. !file already open
           exit
        endif
     enddo
  endif

  if (.not.file_open) then 
     if(num_files_w == max_files_w) &  ! need to have bigger max_files_w
          call mpp_error(FATAL,'fms_io write_data: max_files_w exceeded, increase it via fms_io_nml')    
! record the file name in array files_write
     num_files_w=num_files_w + 1
     nfile = num_files_w           
     files_write(nfile)%filename = trim(fname)         
  endif

! check if the field is new or not and get position and dimension of the field
  index_field =  lookup_field_w(nfile,fieldname)
  if(index_field < 0) then  ! this field is new in file filename
! open temporary file for writing data only, each field is written to a separate file
     temp_name = trim(fname)//'_'//trim(fieldname)//'_tmp'
     call mpp_open(unit2,temp_name,form=MPP_NATIVE,nohdrs=.true.,threading=MPP_MULTI,&
          fileset=MPP_MULTI, action=MPP_OVERWR)               
     files_write(nfile)%nvar = files_write(nfile)%nvar +1
     if(files_write(nfile)%nvar>max_fields) then
        write(error_msg,'(I3,"/",I3)') files_write(nfile)%nvar, max_fields 
        call  mpp_error(FATAL,'fms_io,write_data_3d_new: max_fields exceeded, needs increasing, nvar/max_fields=' &
             //trim(error_msg))
     endif
     index_field = files_write(nfile)%nvar
     files_write(nfile)%fieldname(index_field) = fieldname
     files_write(nfile)%siz(index_field,1) = size(data,1)
     files_write(nfile)%siz(index_field,2) = size(data,2)
     files_write(nfile)%siz(index_field,3) = size(data,3)
     files_write(nfile)%siz(index_field,4) = 1
     files_write(nfile)%unit_tmpfile(index_field) = unit2

     if(PRESENT(domain)) then
        files_write(nfile)%domain_present(index_field)=.true.
        domain_idx = lookup_domain(domain)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io,write_data_3d_new, 1: max_domains exceeded,' &
           //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = domain
        endif
        files_write(nfile)%domain_idx(index_field) = domain_idx
     else if (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
! if domain flag not present, alternatively use domain defined by previous call to set_domain. 
! This is needed in the case of atmospheric physics modules which don't have access to the 
! Domain2d information
        files_write(nfile)%domain_present(index_field)=.true.
        domain_idx = lookup_domain(Current_domain)
        if(domain_idx == -1) then
           num_domains = num_domains + 1
           if(num_domains > max_domains) call  mpp_error(FATAL,'fms_io,write_data_3d_new, 2: max_domains exceeded,' &
           //' needs increasing')
           domain_idx = num_domains
           array_domain(domain_idx) = Current_domain
        endif
        files_write(nfile)%domain_idx(index_field) = domain_idx
     else
        files_write(nfile)%domain_present(index_field)=.false.
     endif

     if (files_write(nfile)%domain_present(index_field)) then
        call mpp_get_global_domain(array_domain(files_write(nfile)%domain_idx(index_field)),xsize=gxsize,ysize=gysize)
        files_write(nfile)%gsiz(index_field,1) = gxsize
        files_write(nfile)%gsiz(index_field,2) = gysize
        files_write(nfile)%gsiz(index_field,3) = size(data,3)
     else
        files_write(nfile)%gsiz(index_field,1) = size(data,1)
        files_write(nfile)%gsiz(index_field,2) = size(data,2)
        files_write(nfile)%gsiz(index_field,3) = size(data,3)
     endif
  else     ! this field  already exists
! get previously stored unit2
     unit2=files_write(nfile)%unit_tmpfile(index_field)
! Increase time level of the fieldname
     files_write(nfile)%siz(index_field,4) =  files_write(nfile)%siz(index_field,4) +1
  endif
  if(thread_w.eq.MPP_SINGLE .and. files_write(nfile)%domain_present(index_field) ) then
     gxsize = files_write(nfile)%gsiz(index_field,1)
     gysize = files_write(nfile)%gsiz(index_field,2)
     allocate (global_data(gxsize,gysize,size(data,3)))
     call mpp_global_field(array_domain(files_write(nfile)%domain_idx(index_field)),data,global_data)
     if(mpp_pe() == mpp_root_pe()) write(unit2) global_data
     deallocate(global_data)
! write data to temporary storage without halos
  else 
     write(unit2) data ! write data to a temporary file
  endif
end subroutine write_data_3d_new
! </SUBROUTINE>  

subroutine write_data_2d_new(filename, fieldname, data, domain,append_pelist_name, no_domain)

  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname 
  real, dimension(:,:), intent(in) :: data
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  type(domain2d), intent(in), optional :: domain
  logical, intent(in), optional :: append_pelist_name, no_domain
  
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_2d_new):need to call fms_io_init first')
  data_3d(:,:,1) = data(:,:)
  call write_data_3d_new(filename, fieldname, data_3d, domain, append_pelist_name, no_domain)
end subroutine write_data_2d_new

! ........................................................
subroutine write_data_1d_new(filename, fieldname, data,domain,append_pelist_name, no_domain)
  
  IMPLICIT NONE
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  real, dimension(:), intent(in) :: data
  real, dimension(size(data(:)),1,1) :: data_3d
  logical, intent(in), optional :: append_pelist_name, no_domain
  
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_1d_new): module not initialized')  
  data_3d(:,1,1) = data(:)  
  call write_data_3d_new(filename, fieldname, data_3d,domain,append_pelist_name, no_domain)  
end subroutine write_data_1d_new

! ..........................................................
subroutine write_data_scalar_new(filename, fieldname, data, domain, append_pelist_name, no_domain)

  IMPLICIT NONE
  type(domain2d), intent(in), optional :: domain
  character(len=*), intent(in) :: filename, fieldname 
  real, intent(in) :: data
  real, dimension(1,1,1) :: data_3d
  logical, intent(in), optional :: append_pelist_name, no_domain
    
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(write_data_scalar_new):  module not initialized: '//fieldname)  
  data_3d(1,1,1) = data
  call write_data_3d_new(filename, fieldname, data_3d,domain,append_pelist_name, no_domain)
end subroutine write_data_scalar_new
! ..........................................................

function lookup_field_w(nfile,fieldname)
  IMPLICIT NONE
! Given fieldname, this function returns the field position in the model's fields list

  integer, intent(in) :: nfile
  character(len=*), intent(in) :: fieldname
  integer :: lookup_field_w
  integer :: j
  character(len=128) :: name

  lookup_field_w=-1
  do j = 1, files_write(nfile)%nvar
     name = files_write(nfile)%fieldname(j)
     if (trim(name) == trim(fieldname)) then
        lookup_field_w = j
        exit 
     endif
  enddo
  return
end function lookup_field_w
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

  IMPLICIT NONE
  integer, intent(in) :: axis_sizes(:), siz
  type(domain1d), optional :: domains(:)
  type(domain1d), optional :: dom
  integer :: lookup_axis
  integer :: i,j
  character(len=128) :: name

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
subroutine field_size(filename, fieldname, siz, append_pelist_name, field_found )

  character(len=*), intent(in) :: filename, fieldname
  integer, intent(inout) :: siz(:)
  logical, intent(in), optional :: append_pelist_name
  logical, intent(out), optional :: field_found
  
  character(len=128) :: name
  character(len=1) :: cart
  integer :: i, nfile, unit, ndim, nvar, natt, ntime, siz_in(4), j, len
  logical :: file_opened, found
  character(len=256) :: fname

  type(fieldtype) :: fields(max_fields)
  type(axistype) :: axes(max_fields)
  if (size(siz(:)) < 4) call mpp_error(FATAL,'fms_io(field_size): size array must be >=4 to receive field size of ' &
       //trim(fieldname)//' in file '// trim(filename))

! Need to check if filename has been opened or not

  fname = trim(filename)
  if (PRESENT(append_pelist_name)) then
     if (append_pelist_name) then
        allocate(pelist(mpp_npes()))
        call mpp_get_current_pelist(pelist,pelist_name)
        fname = trim(filename)//trim(pelist_name)
        deallocate(pelist)
     endif
  endif        
  nfile = 0
  file_opened=.false.
  do i=1,num_files_r
     if (trim(files_read(i)%filename) == trim(fname))  then
        nfile = i
        file_opened = .true.
        exit ! file is already opened
     endif
  enddo
!Need to open the file now, Only works for single NetCDF files for now ...
  found= .false.
  siz=-1
  if (.not. file_opened) then
     call mpp_open(unit,trim(fname),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_SINGLE, &
          fileset=MPP_SINGLE)
     call mpp_get_info(unit,ndim,nvar,natt,ntime)
     if (nvar > max_fields) then
        write(error_msg,'(I3,"/",I3)') nvar,max_fields 
        call  mpp_error(FATAL,'fms_io(field_size): max_fields too small, needs increasing, nvar/max_fields=' &
             //trim(error_msg)//' in file '//trim(filename))
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
     call mpp_close(unit)
  else
     do i=1, files_read(nfile)%nvar
        if (trim(fieldname) == trim(files_read(nfile)%fieldname(i))) then
           found = .true.
           siz = files_read(nfile)%siz(i,:)
           exit
        endif
     enddo
     if (.not. found) then
        call mpp_get_info(files_read(nfile)%unit,ndim,nvar,natt,ntime)
        call mpp_get_fields(files_read(nfile)%unit,fields(1:nvar))
        do i=1, nvar
           call mpp_get_atts(fields(i),name=name)
           if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
              call mpp_get_atts(fields(i),ndim=ndim)
              call mpp_get_atts(fields(i),axes=axes(1:ndim))
              call mpp_get_atts(fields(i),siz=siz_in)
              siz = siz_in
              siz(4) = ntime
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
              found=.true.
              exit
           endif
        enddo
     endif
     if(.not. found) then         
        do i=1, files_read(nfile)%ndim
           call mpp_get_atts(files_read(nfile)%axes(i),name=name, len= siz_in(1))
           if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
              siz(1)= siz_in(1)
              found = .true.
              exit
           endif
        enddo
     endif
 endif
 if( PRESENT(field_found) )then
     field_found = found
 else if (.not. found .and. mpp_pe() == mpp_root_pe() )then
     call mpp_error(FATAL, 'fms_io(field_size): field '//trim(fieldname)// ' NOT found in file '//trim(filename))
 end if
  call mpp_sync() !is this needed?
  return
end subroutine field_size
! </SUBROUTINE>


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
subroutine read_data_i3d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  integer, dimension(:,:,:), intent(out) :: data ! 3 dimensional data    
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain 

  real, dimension(size(data,1),size(data,2),size(data,3)) :: r_data
  call read_data_3d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, no_domain)
  data = CEILING(r_data)
end subroutine read_data_i3d_new
subroutine read_data_i2d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  integer, dimension(:,:), intent(out) :: data ! 2 dimensional data    
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name , no_domain

  real, dimension(size(data,1),size(data,2)) :: r_data
  call read_data_2d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, no_domain)
  data = CEILING(r_data)
end subroutine read_data_i2d_new
subroutine read_data_i1d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  integer, dimension(:), intent(out) :: data ! 1 dimensional data    
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain 

  real, dimension(size(data,1)) :: r_data
  call read_data_1d_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, no_domain)
  data = CEILING(r_data)
end subroutine read_data_i1d_new
subroutine read_data_iscalar_new(filename,fieldname,data,domain,timelevel,append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  integer, intent(out) :: data     
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain 

  real :: r_data
  call read_data_scalar_new(filename,fieldname,r_data,domain,timelevel,append_pelist_name, no_domain)
  data = CEILING(r_data)
end subroutine read_data_iscalar_new
!=====================================================================================
subroutine read_data_3d_new(filename,fieldname,data,domain,timelevel,append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  real, dimension(:,:,:), intent(out) :: data ! 3 dimensional data    
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain 
  character(len=128) :: name
  character(len=256) :: fname
  integer :: unit, siz_in(4), siz(4), i, j, k
  integer :: nfile  ! index of the opened file in array files
  integer :: ndim, nvar, natt, ntime,var_dim, tlev=1
  integer :: index_field ! position of the fieldname in the list of variables
  integer :: iscomp, iecomp, jscomp, jecomp, cxsize, cysize,cxsize_max, cysize_max
  integer :: isdata, iedata, jsdata, jedata, dxsize, dysize,dxsize_max, dysize_max
  integer :: isglobal, ieglobal, jsglobal, jeglobal, gxsize, gysize,gxsize_max,gysize_max
  logical :: data_is_global
  logical :: file_opened, found, is_no_domain = .false.
  integer :: index_axis
  type(file_type), dimension(:), pointer  :: new_files_read
! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_3d_new):  module not initialized')  
  data_is_global=.false.
  file_opened=.false.
  if (PRESENT(timelevel)) then
     tlev = timelevel
  else
     tlev = 1
  endif
  fname = trim(filename)   
  if (PRESENT(append_pelist_name)) then
     if (append_pelist_name) then
        allocate(pelist(mpp_npes()))        
        call mpp_get_current_pelist(pelist,pelist_name)
        fname = trim(filename)//trim(pelist_name)
        deallocate(pelist)
     endif
  endif    

  is_no_domain = .false.
  if (PRESENT(no_domain)) THEN
     if(PRESENT(domain) .AND. no_domain) &
       call mpp_error(FATAL, 'fms_io(write_data_3d_new): no_domain cannot be .true. when optional argument domain is present.')
     is_no_domain = no_domain
  endif

  if (PRESENT(domain) ) then
     call mpp_get_compute_domain(domain,iscomp,iecomp,jscomp,jecomp,cxsize, & 
          cxsize_max,cysize,cysize_max)
     call mpp_get_data_domain(domain,isdata,iedata,jsdata,jedata,dxsize,dxsize_max,&
          dysize,dysize_max)
     call mpp_get_global_domain(domain,isglobal,ieglobal,jsglobal,jeglobal,gxsize, &
          gxsize_max,gysize,gysize_max)
     if (gxsize == size(data,1) .and. gysize == size(data,2)) data_is_global = .true.
  else  if (ASSOCIATED(Current_domain)  .AND. .NOT. is_no_domain ) then
     gxsize=ieg-isg+1
     gysize=jeg-jsg+1
     dxsize=ied-isd+1
     dysize=jed-jsd+1
     cxsize=ie-is+1
     cysize=je-js+1
     if (gxsize == size(data,1) .and. gysize == size(data,2)) data_is_global = .true.
  else 

     data_is_global =.true.
     gxsize = size(data,1)
     gysize = size(data,2)
     siz_in(3) = size(data,3)
  endif

  if (data_is_global .and. fset_r == MPP_MULTI) &
       call mpp_error(FATAL,'fms_io(read_data_3d_new): can not do global read on multi fileset ' &
       //'change fileset_reading to single')  
! Need to check if filename has been opened or not
  nfile = 0
  if(num_files_r ==0) then
     file_opened = .false.
  else
     file_opened=.false.
     do i=1,num_files_r
        if (files_read(i)%filename == trim(fname))  then
           nfile = i
           file_opened = .true.
           exit ! file is already opened
        endif
     enddo
  endif
  if (.not. file_opened) then !Need to open the file now     
     if (fset_r == MPP_MULTI .and. thread_r == MPP_SINGLE) then
        call mpp_error(FATAL,'fms_io(read_data_3d_new): single-threaded read from multi fileset not allowed' &
             //'change either threading_read or fileset_read')
     endif
     call mpp_open(unit,trim(fname),form=form,action=MPP_RDONLY,threading=thread_r, &
          fileset=fset_r)
! Increase num_files_r and set file_type 
     if(num_files_r == max_files_r) &  ! need to have bigger max_files_r
          call mpp_error(FATAL,'fms_io read_data: max_files_r exceeded, increase it via fms_io_nml')
     num_files_r=num_files_r + 1
     files_read(num_files_r)%filename = trim(fname)
     nfile = num_files_r
     files_read(nfile)%unit = unit
  else
     unit = files_read(nfile)%unit
  endif  
! Get info of this file and field   
  if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
     call mpp_get_info(unit,files_read(nfile)%ndim, &
          files_read(nfile)%nvar, files_read(nfile)%natt,files_read(nfile)%max_ntime)
     if(files_read(nfile)%max_ntime < 1)  files_read(nfile)%max_ntime = 1
     nvar = files_read(nfile)%nvar  
     ndim = files_read(nfile)%ndim
     if (nvar > max_fields) then
        write(error_msg,'(I3,"/",I3)') nvar,max_fields
        call mpp_error(FATAL,'fms_io(read_data_3d_new)1: max_fields too small needs increasing,nvar/max_fields=' &
             //trim(error_msg)//'in file'//trim(filename))
     endif
     call mpp_get_fields(unit,files_read(nfile)%fields(1:nvar))     
     call mpp_get_axes(unit,files_read(nfile)%axes(1:ndim))
     siz_in = 1
     index_field = -1
     found = .false.
     files_read(nfile)%is_dimvar(:) = .false.
     do i=1,files_read(nfile)%nvar
        call mpp_get_atts(files_read(nfile)%fields(i),name=name,ndim=var_dim,siz=siz_in)
        if(var_dim .lt.3) then
           do j=var_dim+1,3
              siz_in(j)=1
           enddo
        endif
        if (lowercase(trim(name)) == lowercase(trim(fieldname))) then ! found the variable
           index_field = i
           files_read(nfile)%fieldname(i) = fieldname
           files_read(nfile)%siz(i,:)  = siz_in
           files_read(nfile)%gsiz(i,:) = siz_in
           if (fset_r == MPP_SINGLE) then

              if (siz_in(1) /= gxsize .or. siz_in(2) /= gysize .or. siz_in(3) /= &
                   size(data,3)) then
                 PRINT *, gxsize, gysize, size(data, 3), siz_in(1), siz_in(2), siz_in(3)
                 call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                   ' in file '//trim(filename)//': field size mismatch 1')
              endif
           else if (fset_r == MPP_MULTI) then
              if (siz_in(1) /= dxsize .or. siz_in(2) /= dysize & 
                   .or. siz_in(3) /= size(data,3)) &
                   call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                   ' in file '//trim(filename)//': field size mismatch 2')
           endif
           found = .true.
           exit  !jump out of i loop
        endif
     enddo
     if(.not. found) then
        if (nvar+ndim > max_fields) then
           write(error_msg,'(I3,"/",I3)') nvar+ndim, max_fields 
           call  mpp_error(FATAL,'fms_io(read_data_3d_new)2: max_fields exceeded, needs increasing, nvar/max_fields=' &
                //trim(error_msg)//' in file '//trim(filename))             
        endif
        do i=1,files_read(nfile)%ndim
           call mpp_get_atts(files_read(nfile)%axes(i),name=name, len= siz_in(1)) 
           if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
              index_field = i+nvar
              files_read(nfile)%is_dimvar(index_field) = .true.
              files_read(nfile)%fieldname(index_field) = fieldname
              files_read(nfile)%siz(index_field,:)  = siz_in
              files_read(nfile)%gsiz(index_field,:) = siz_in
              if (fset_r == MPP_SINGLE) then
                 if (siz_in(1) /= gxsize) &
                      call mpp_error(FATAL,'fms_io(read_data_3d_new), field '//trim(fieldname)// &
                      ' in file '//trim(filename)//' field size mismatch 3')
              endif
              siz(1)= siz_in(1)
              found = .true.
              exit
           endif
        enddo
     endif

!     PRINT *, 'fms_io: ', tlev, fieldname, files_read(nfile)%max_ntime, data_is_global, index_field, files_read(nfile)%nvar 
     if(index_field <1) call mpp_error(FATAL, 'fms_io, read_data_3d_new: field '//trim(fieldname)// &
          ' NOT found in file '//trim(filename))
     if ( tlev < 1 .or. files_read(nfile)%max_ntime < tlev)  then
        write(error_msg,'(I5,"/",I5)') tlev, files_read(nfile)%max_ntime
        call mpp_error(FATAL,'fms_io(read_data_3d_new): time level out of range, time level/max_time_level=' &
             //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
     endif
     if (data_is_global) then
        if (files_read(nfile)%is_dimvar(index_field)) then
           index_axis = index_field - files_read(nfile)%nvar 
           call mpp_get_axis_data( files_read(nfile)%axes(index_axis),data(:,1,1))
        else
           call mpp_read(unit,files_read(nfile)%fields(index_field),data(:,:,:),tlev)
        endif
     else 
        if (files_read(nfile)%is_dimvar(index_field)) call mpp_error(FATAL,'fms_io(read_data_3d_new): domain is present' &
             //' but the variable is a dimension variable.  Remove domain flag for this var')
        if (PRESENT(domain)) then
           call mpp_read(unit,files_read(nfile)%fields(index_field),domain,data,tlev)
        else
           call mpp_read(unit, files_read(nfile)%fields(index_field),Current_domain,data,tlev)
        endif
     endif
  endif  
  return
end subroutine read_data_3d_new
!.............................................................. 
! </SUBROUTINE>
subroutine read_data_2d_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  real, dimension(:,:), intent(out) :: data     !2 dimensional data 
  real, dimension(size(data,1),size(data,2),1) :: data_3d
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain
  data_3d = 0.0

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain)
  data(:,:) = data_3d(:,:,1)
end subroutine read_data_2d_new
!.....................................................................
subroutine read_data_1d_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain)
  IMPLICIT NONE
  character(len=*), intent(in) :: filename, fieldname
  real, dimension(:), intent(out) :: data     !1 dimensional data 
  real, dimension(size(data,1),1,1) :: data_3d
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain
  
  data_3d = 0.0  

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain)

  data(:) = data_3d(:,1,1)  
end subroutine read_data_1d_new
!.....................................................................

subroutine read_data_scalar_new(filename,fieldname,data,domain,timelevel,&
     append_pelist_name, no_domain)
  IMPLICIT NONE
! this subroutine is for reading a single number
  character(len=*), intent(in) :: filename, fieldname
  real, intent(out) :: data     !zero dimension data 
  real, dimension(1,1,1) :: data_3d
  type(domain2d), intent(in), optional :: domain
  integer, intent(in) , optional :: timelevel
  logical, intent(in), optional :: append_pelist_name, no_domain

  data_3d = 0.0

  call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel,&
       append_pelist_name, no_domain)

  data = data_3d(1,1,1)

end subroutine read_data_scalar_new
!.....................................................................

function unique_axes(array,dom)
  IMPLICIT NONE
  integer, dimension(:), intent(in) :: array
  type(domain1d), dimension(:), intent(in), optional :: dom
  integer, dimension(max_axes) :: unique_axes
  integer :: i,n,j
  logical :: dup
  
  unique_axes=-1
  n=1  
  do i=1,size(array(:))
     dup = .false.
     if (array(i) == 0) exit 
     do j=i-1,1,-1
        if (array(i) == array(j)) then
           if (PRESENT(dom)) then
              if (dom(i) == dom(j)) then
                 dup = .true.
                 exit
              endif
           else
              dup = .true.
              exit
           endif
        endif
     enddo
     if (.not. dup) then
!          unique_axes(n) = array(i)
        unique_axes(n) = i
        n=n+1
     endif
  enddo
  
end function unique_axes

  !#######################################################################
  !#######################################################################
  !   --------- routines for reading distributed data ---------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !
  ! reading can be done either by all PEs (default) or by only the root PE
  ! this is controlled by namelist variable "read_all_pe".
  !#######################################################################

subroutine read_data_2d ( unit, data, end )

  integer, intent(in)                        :: unit
  real,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  real, dimension(isg:ieg,jsg:jeg) :: gdata
  integer :: len

  include "read_data_2d.inc"  
end subroutine read_data_2d

!#######################################################################

subroutine read_ldata_2d ( unit, data, end )

  integer, intent(in)                        :: unit
  logical, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end  
  logical, dimension(isg:ieg,jsg:jeg) :: gdata
  integer :: len

  include "read_data_2d.inc"
end subroutine read_ldata_2d
!#######################################################################

subroutine read_idata_2d ( unit, data, end )

  integer, intent(in)                        :: unit
  integer, intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end
  
  integer, dimension(isg:ieg,jsg:jeg) :: gdata
  integer :: len

  include "read_data_2d.inc"  
end subroutine read_idata_2d

!#######################################################################

subroutine read_cdata_2d ( unit, data, end )

  integer, intent(in)                        :: unit
  complex,    intent(out), dimension(isd:,jsd:) :: data
  logical, intent(out), optional             :: end
  complex, dimension(isg:ieg,jsg:jeg) :: gdata
  integer :: len

  include "read_data_2d.inc"
end subroutine read_cdata_2d

!#######################################################################

subroutine read_data_3d ( unit, data, end )

  integer, intent(in)                          :: unit
  real,    intent(out), dimension(isd:,jsd:,:) :: data
  logical, intent(out), optional               :: end  
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer :: len

  include "read_data_3d.inc"
end subroutine read_data_3d

!#######################################################################

subroutine read_cdata_3d ( unit, data, end )

  integer, intent(in)                          :: unit
  complex, intent(out), dimension(isd:,jsd:,:) :: data
  logical, intent(out), optional               :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
  integer :: len

  include "read_data_3d.inc"  
end subroutine read_cdata_3d

!#######################################################################

subroutine read_data_4d ( unit, data, end )

  integer, intent(in)                            :: unit
  real,    intent(out), dimension(isd:,jsd:,:,:) :: data
  logical, intent(out), optional                 :: end
  real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: len

! WARNING: memory usage with this routine could be costly
 
  include "read_data_4d.inc"  
end subroutine read_data_4d

!#######################################################################

subroutine read_cdata_4d ( unit, data, end )

  integer, intent(in)                            :: unit
  complex, intent(out), dimension(isd:,jsd:,:,:) :: data
  logical, intent(out), optional                 :: end
  complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
  integer :: len

! WARNING: memory usage with this routine could be costly
  include "read_data_4d.inc"
end subroutine read_cdata_4d

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

subroutine write_cdata_2d ( unit, data )

  integer, intent(in)                       :: unit
  complex, intent(in), dimension(isd:,jsd:) :: data
  complex, dimension(isg:ieg,jsg:jeg) :: gdata

  include "write_data.inc"
end subroutine write_cdata_2d

!#######################################################################

subroutine write_data_3d ( unit, data )

  integer, intent(in) :: unit
  real,    intent(in), dimension(isd:,jsd:,:) :: data
  real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
    
  include "write_data.inc"
end subroutine write_data_3d

!#######################################################################

subroutine write_cdata_3d ( unit, data )

  integer, intent(in) :: unit
  complex, intent(in), dimension(isd:,jsd:,:) :: data
  complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

  include "write_data.inc"
end subroutine write_cdata_3d

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
  integer :: i,j

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


  
end module fms_io_mod


