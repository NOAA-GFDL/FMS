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


module time_interp_external2_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov">M.J. Harrison</CONTACT>
!
!<REVIEWER EMAIL="hsimmons@iarc.uaf.edu">Harper Simmons</REVIEWER>
!
!<OVERVIEW>
! Perform I/O and time interpolation of external fields (contained in a file).
!</OVERVIEW>

!<DESCRIPTION>
! Perform I/O and time interpolation for external fields.
! Uses udunits library to calculate calendar dates and
! convert units.  Allows for reading data decomposed across
! model horizontal grid using optional domain2d argument
!
! data are defined over data domain for domain2d data
! (halo values are NOT updated by this module)
!
!</DESCRIPTION>
!
!<NAMELIST NAME="time_interp_external_nml">
! <DATA NAME="num_io_buffers" TYPE="integer">
! size of record dimension for internal buffer.  This is useful for tuning i/o performance
! particularly for large datasets (e.g. daily flux fields)
! </DATA>
!</NAMELIST>

  use platform_mod, only : DOUBLE_KIND => r8_kind
  use fms_mod, only : write_version_number
  use mpp_mod, only : mpp_error,FATAL,WARNING,mpp_pe, stdout, stdlog, NOTE
  use mpp_mod, only : input_nml_file
  use time_manager_mod, only : time_type, get_date, set_date, operator ( >= ) , operator ( + ) , days_in_month, &
                            operator( - ), operator ( / ) , days_in_year, increment_time, &
                            set_time, get_time, operator( > ), get_calendar_type, NO_CALENDAR
  use get_cal_time_mod, only : get_cal_time
  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_data_domain, &
       mpp_get_global_domain, NULL_DOMAIN2D
  use time_interp_mod, only : time_interp, time_interp_init
  use axis_utils2_mod, only : get_axis_cart, get_axis_modulo, get_axis_modulo_times
  use fms_mod, only : lowercase, check_nml_error
  use platform_mod, only: r8_kind
  use horiz_interp_mod, only : horiz_interp, horiz_interp_type
  use fms2_io_mod,      only : Valid_t, FmsNetcdfDomainFile_t, open_file, get_unlimited_dimension_name, &
                               variable_att_exists, FmsNetcdfFile_t, &
                               variable_exists, get_valid, get_variable_num_dimensions, read_data, &
                               is_valid, close_file, get_dimension_size, get_variable_dimension_names, &
                               get_variable_size, get_time_calendar, get_variable_missing, get_variable_units

  implicit none
  private

! Include variable "version" to be written to log file.
#include<file_version.h>

  integer, parameter, public  :: NO_REGION=0, INSIDE_REGION=1, OUTSIDE_REGION=2
  integer, parameter, private :: modulo_year= 0001
  integer, parameter, private :: LINEAR_TIME_INTERP = 1 ! not used currently
  integer, parameter, public  :: SUCCESS = 0, ERR_FIELD_NOT_FOUND = 1
  real,    parameter, private :: DEFAULT_MISSING_VALUE = -1e20
  integer,            private :: max_fields = 100, max_files= 40
  integer, private :: num_fields = 0, num_files=0
  ! denotes time intervals in file (interpreted from metadata)
  integer, private :: num_io_buffers = 2 ! set -1 to read all records from disk into memory
  logical, private :: module_initialized = .false.
  logical, private :: debug_this_module = .false.

  public init_external_field, time_interp_external, time_interp_external_init, &
       time_interp_external_exit, get_external_field_size, get_time_axis, get_external_field_missing
  public set_override_region, reset_src_data_region
  public get_external_fileobj

  private find_buf_index,&
         set_time_modulo

  type, private :: ext_fieldtype
     type(FmsNetcdfFile_t), pointer :: fileobj=>NULL() ! keep unit open when not reading all records
     character(len=128) :: name, units
     integer :: siz(4), ndim
     character(len=32) :: axisname(4)
     type(domain2d) :: domain
     type(time_type), dimension(:), pointer :: time =>NULL() ! midpoint of time interval
     type(time_type), dimension(:), pointer :: start_time =>NULL(), end_time =>NULL()
     type(time_type), dimension(:), pointer :: period =>NULL()
     logical :: modulo_time ! denote climatological time axis
     real, dimension(:,:,:,:), pointer :: data =>NULL() ! defined over data domain or global domain
     logical, dimension(:,:,:,:), pointer :: mask =>NULL() ! defined over data domain or global domain
     integer, dimension(:), pointer :: ibuf  =>NULL() ! record numbers associated with buffers
     real, dimension(:,:,:,:), pointer :: src_data  =>NULL() ! input data buffer
     type(valid_t) :: valid ! data validator
     integer :: nbuf
     logical :: domain_present
     real(DOUBLE_KIND) :: slope, intercept
     integer :: isc,iec,jsc,jec
     type(time_type) :: modulo_time_beg, modulo_time_end
     logical :: have_modulo_times, correct_leap_year_inconsistency
     integer :: region_type
     integer :: is_region, ie_region, js_region, je_region
     integer :: is_src, ie_src, js_src, je_src
     integer :: tdim
     integer :: numwindows
     logical, dimension(:,:), pointer :: need_compute=>NULL()
     real    :: missing ! missing value
  end type ext_fieldtype

  type, private :: filetype
     character(len=128) :: filename = ''
     type(FmsNetcdfFile_t), pointer :: fileobj => NULL()
  end type filetype

  interface time_interp_external
     module procedure time_interp_external_0d
     module procedure time_interp_external_2d
     module procedure time_interp_external_3d
  end interface

  integer :: outunit

  type(ext_fieldtype), save, private, pointer :: field(:) => NULL()
  type(filetype),      save, private, pointer :: opened_files(:) => NULL()
!Balaji: really should use field%missing
  integer, private, parameter :: dk = DOUBLE_KIND ! ensures that time_interp_missing is in range for mixed-mode
                                                  ! compiling
  real(DOUBLE_KIND), private, parameter :: time_interp_missing=-1e99_dk
  contains

! <SUBROUTINE NAME="time_interp_external_init">
!
! <DESCRIPTION>
! Initialize the time_interp_external module
! </DESCRIPTION>
!
    subroutine time_interp_external_init()

      integer :: ioun, io_status, logunit, ierr

      namelist /time_interp_external_nml/ num_io_buffers, debug_this_module, &
                                          max_fields, max_files

      ! open and read namelist

      if(module_initialized) return

      logunit = stdlog()
      outunit = stdout()
      call write_version_number("TIME_INTERP_EXTERNAL_MOD", version)

      read (input_nml_file, time_interp_external_nml, iostat=io_status)
      ierr = check_nml_error(io_status, 'time_interp_external_nml')

      write(logunit,time_interp_external_nml)
      call realloc_fields(max_fields)
      call realloc_files(max_files)

      module_initialized = .true.

      call time_interp_init()

      return

    end  subroutine time_interp_external_init
! </SUBROUTINE> NAME="time_interp_external_init"


!<FUNCTION NAME="init_external_field" TYPE="integer">
!
!<DESCRIPTION>
! initialize an external field.  Buffer "num_io_buffers" (default=2) in memory to reduce memory allocations.
! distributed reads are supported using the optional "domain" flag.
! Units conversion via the optional "desired_units" flag using udunits_mod.
!
! Return integer id of field for future calls to time_interp_external.
!
! </DESCRIPTION>
!
!<IN  NAME="file" TYPE="character(len=*)">
! filename
!</IN>
!<IN  NAME="fieldname" TYPE="character(len=*)">
! fieldname (in file)
!</IN>
!<IN NAME="format" TYPE="integer">
! mpp_io flag for format of file (optional). Currently only "MPP_NETCDF" supported
!</IN>
!<IN NAME="threading" TYPE="integer">
! mpp_io flag for threading (optional).  "MPP_SINGLE" means root pe reads global field and distributes to other PEs
! "MPP_MULTI" means all PEs read data
!</IN>
!<IN NAME="domain" TYPE="mpp_domains_mod:domain2d">
! domain flag (optional)
!</IN>
!<IN NAME="desired_units" TYPE="character(len=*)">
! Target units for data (optional), e.g. convert from deg_K to deg_C.
! Failure to convert using udunits will result in failure of this module.
!</IN>
!<IN NAME="verbose" TYPE="logical">
! verbose flag for debugging (optional).
!</IN>
!<INOUT NAME="axis_centers" TYPE="axistype" DIM="(4)">
! MPP_IO axistype array for grid centers ordered X-Y-Z-T (optional).
!</INOUT>
!<INOUT NAME="axis_sizes" TYPE="integer" DIM="(4)">
!  array of axis lengths ordered X-Y-Z-T (optional).
!</INOUT>


    function init_external_field(file,fieldname,domain,desired_units,&
         verbose,axis_names, axis_sizes,override,correct_leap_year_inconsistency,&
         permit_calendar_conversion,use_comp_domain,ierr, nwindows, ignore_axis_atts, ongrid )

      character(len=*), intent(in)            :: file,fieldname
      logical, intent(in), optional           :: verbose
      character(len=*), intent(in), optional  :: desired_units
      type(domain2d), intent(in), optional    :: domain
      integer, intent(inout), optional        :: axis_sizes(4)
      character(len=*), intent(out), optional :: axis_names(4)
      logical, intent(in), optional           :: override, correct_leap_year_inconsistency,&
                                                 permit_calendar_conversion,use_comp_domain
      integer,          intent(out), optional :: ierr
      integer,          intent(in),  optional :: nwindows
      logical, optional                       :: ignore_axis_atts
      logical, optional                       :: ongrid !< Optional flag indicating if the data is ongrid

      logical :: ongrid_local !< Flag indicating if the data is ongrid

      integer :: init_external_field

      real(DOUBLE_KIND) :: slope, intercept
      integer :: unit,ndim,nvar,natt,ntime,i,j
      integer :: iscomp,iecomp,jscomp,jecomp,isglobal,ieglobal,jsglobal,jeglobal
      integer :: isdata,iedata,jsdata,jedata, dxsize, dysize,dxsize_max,dysize_max
      logical :: verb, transpose_xy,use_comp_domain1
      real, dimension(:), allocatable :: tstamp, tstart, tend, tavg
      character(len=1) :: cart
      character(len=1), dimension(4) :: cart_dir
      character(len=128) :: units, fld_units
      character(len=128) :: name, msg, calendar_type, timebeg, timeend
      character(len=128) :: timename, timeunits
      character(len=128), allocatable :: axisname(:)
      integer,            allocatable :: axislen(:)
      integer :: siz(4), siz_in(4), gxsize, gysize,gxsize_max, gysize_max
      type(time_type) :: tdiff
      integer :: yr, mon, day, hr, minu, sec
      integer :: len, nfile, nfields_orig, nbuf, nx,ny
      integer :: numwindows
      logical :: ignore_axatts
      logical :: have_modulo_time
      type(FmsNetcdfFile_t), pointer :: fileobj=>NULL()

      if (.not. module_initialized) call mpp_error(FATAL,'Must call time_interp_external_init first')
      if(present(ierr)) ierr = SUCCESS
      ignore_axatts=.false.
      cart_dir(1)='X';cart_dir(2)='Y';cart_dir(3)='Z';cart_dir(4)='T'
      if(present(ignore_axis_atts)) ignore_axatts = ignore_axis_atts
      use_comp_domain1 = .false.
      if(PRESENT(use_comp_domain)) use_comp_domain1 = use_comp_domain
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.
      numwindows = 1
      if(present(nwindows)) numwindows = nwindows

      units = 'same'
      if (PRESENT(desired_units)) then
          units = desired_units
          call mpp_error(FATAL,'==> Unit conversion via time_interp_external &
               &has been temporarily deprecated.  Previous versions of&
               &this module used udunits_mod to perform unit conversion.&
               &  Udunits_mod is in the process of being replaced since &
               &there were portability issues associated with this code.&
               & Please remove the desired_units argument from calls to &
               &this routine.')
      endif
      nfile = 0
      do i=1,num_files
         if(trim(opened_files(i)%filename) == trim(file)) then
            nfile = i
            exit  ! file is already opened
         endif
      enddo
      if(nfile == 0) then
         num_files = num_files + 1
         if(num_files > max_files) then ! not enough space in the file table, reallocate it
            !--- z1l: For the case of multiple thread, realoc_files will cause memory leak.
            !---      If multiple threads are working on file A. One of the thread finished first and
            !---      begin to work on file B, the realloc_files will cause problem for
            !---      other threads are working on the file A.
            !  call realloc_files(2*size(opened_files))
            call mpp_error(FATAL, "time_interp_external: num_files is greater than max_files, "// &
                                      "increase time_interp_external_nml max_files")
         endif
         opened_files(num_files)%filename = trim(file)
         allocate(opened_files(num_files)%fileobj)
         fileobj => opened_files(num_files)%fileobj

         if(.not. open_file(opened_files(num_files)%fileobj, trim(file), 'read')) &
             call mpp_error(FATAL, 'time_interp_external_mod: Error in opening file '//trim(file))
      else
         fileobj => opened_files(nfile)%fileobj
      endif

      call get_unlimited_dimension_name(fileobj, timename)
      call get_dimension_size(fileobj, timename, ntime)

      if (ntime < 1) then
          write(msg,'(a15,a,a58)') 'external field ',trim(fieldname),&
           ' does not have an associated record dimension (REQUIRED) '
          call mpp_error(FATAL,trim(msg))
      endif

      !--- get time calendar_type
      call get_time_calendar(fileobj, timename, calendar_type)

      !--- get time units
      call get_variable_units(fileobj, timename, timeunits)

      !--- get timebeg and timeend
      have_modulo_time = get_axis_modulo_times(fileobj, timename, timebeg, timeend)

      allocate(tstamp(ntime),tstart(ntime),tend(ntime),tavg(ntime))
      call read_data(fileobj, timename, tstamp)
      transpose_xy = .false.
      isdata=1; iedata=1; jsdata=1; jedata=1
      gxsize=1; gysize=1
      siz_in = 1

      if (PRESENT(domain)) then
         call mpp_get_compute_domain(domain,iscomp,iecomp,jscomp,jecomp)
         nx = iecomp-iscomp+1; ny = jecomp-jscomp+1
         call mpp_get_data_domain(domain,isdata,iedata,jsdata,jedata,dxsize,dxsize_max,dysize,dysize_max)
         call mpp_get_global_domain(domain,isglobal,ieglobal,jsglobal,jeglobal,gxsize,gxsize_max,gysize,gysize_max)
         ongrid_local = .false.
         if (present(ongrid)) ongrid_local = ongrid
         !> If this is an ongrid case, set is[e]js[e]data to be equal to the compute domain.
         !! This is what it is used to allocate space for the data!
         if (ongrid_local) then
              isdata=iscomp
              iedata=iecomp
              jsdata=jscomp
              jedata=jecomp
         endif
      elseif(use_comp_domain1) then
         call mpp_error(FATAL,"init_external_field:"//&
              " use_comp_domain=true but domain is not present")
      endif

      init_external_field = -1
      nfields_orig = num_fields

      tavg = -1.0
      tstart = tstamp
      tend = tstamp
      if(variable_att_exists(fileobj, fieldname, 'time_avg_info')) then
        if(variable_exists(fileobj, 'average_T1')) call read_data(fileobj, 'average_T1', tstart)
        if(variable_exists(fileobj, 'average_T2')) call read_data(fileobj, 'average_T2', tend)
        if(variable_exists(fileobj, 'average_DT')) call read_data(fileobj, 'average_DT', tavg)
      endif

      if (.not. variable_exists(fileobj, fieldname) ) then
        if (present(ierr)) then
           ierr = ERR_FIELD_NOT_FOUND
           return
        else
           call mpp_error(FATAL,'external field "'//trim(fieldname)//'" not found in file "'//trim(file)//'"')
        endif
      endif

      num_fields = num_fields + 1
      if(num_fields > max_fields) then
         !--- z1l: For the case of multiple thread, realoc_fields will cause memory leak.
         !---      If multiple threads are working on field A. One of the thread finished first and
         !---      begin to work on field B, the realloc_files will cause problem for
         !---      other threads are working on the field A.
         !call realloc_fields(size(field)*2)
         call mpp_error(FATAL, "time_interp_external: num_fields is greater than max_fields, "// &
                                   "increase time_interp_external_nml max_fields")
      endif

      !--- get field attribute
      call get_variable_units(fileobj, fieldname, fld_units)

      init_external_field = num_fields
      field(num_fields)%fileobj => fileobj
      field(num_fields)%name = trim(fieldname)
      field(num_fields)%units = trim(fld_units)
      field(num_fields)%isc = 1
      field(num_fields)%iec = 1
      field(num_fields)%jsc = 1
      field(num_fields)%jec = 1
      field(num_fields)%region_type = NO_REGION
      field(num_fields)%is_region   = 0
      field(num_fields)%ie_region   = -1
      field(num_fields)%js_region   = 0
      field(num_fields)%je_region   = -1
      if (PRESENT(domain)) then
         field(num_fields)%domain_present = .true.
         field(num_fields)%domain = domain
         field(num_fields)%isc=iscomp;field(num_fields)%iec = iecomp
         field(num_fields)%jsc=jscomp;field(num_fields)%jec = jecomp
      else
         field(num_fields)%domain_present = .false.
      endif

      field(num_fields)%valid = get_valid(fileobj, fieldname)
      ndim = get_variable_num_dimensions(fileobj, fieldname)
      if (ndim > 4) call mpp_error(FATAL, &
           'invalid array rank <=4d fields supported')

      field(num_fields)%siz = 1
      field(num_fields)%ndim = ndim
      field(num_fields)%tdim = 4
      !--- get field missing value
      field(num_fields)%missing = get_variable_missing(fileobj, fieldname)

      allocate(axisname(ndim), axislen(ndim))

      call get_variable_dimension_names(fileobj, fieldname, axisname)
      call get_variable_size(fileobj, fieldname, axislen)
      do j=1,field(num_fields)%ndim
         call get_axis_cart(fileobj, axisname(j), cart)
         len = axislen(j)
         if (cart == 'N' .and. .not. ignore_axatts) then
            write(msg,'(a,"/",a)')  trim(file),trim(fieldname)
            call mpp_error(FATAL,'file/field '//trim(msg)// &
                 ' couldnt recognize axis atts in time_interp_external')
         else if (cart == 'N' .and. ignore_axatts) then
            cart = cart_dir(j)
         endif
         select case (cart)
         case ('X')
            if (j.eq.2) transpose_xy = .true.
            if (.not.PRESENT(domain) .and. .not.PRESENT(override)) then
               isdata=1;iedata=len
               iscomp=1;iecomp=len
               gxsize = len
               dxsize = len
               field(num_fields)%isc=iscomp;field(num_fields)%iec=iecomp
            elseif (PRESENT(override)) then
               gxsize = len
               if (PRESENT(axis_sizes)) axis_sizes(1) = len
            endif
            field(num_fields)%axisname(1) = axisname(j)
            if(use_comp_domain1) then
               field(num_fields)%siz(1) = nx
            else
               field(num_fields)%siz(1) = dxsize
            endif
            if (len /= gxsize) then
               write(msg,'(a,"/",a)')  trim(file),trim(fieldname)
               call mpp_error(FATAL,'time_interp_ext, file/field '//trim(msg)//' x dim doesnt match model')
            endif
         case ('Y')
            field(num_fields)%axisname(2) = axisname(j)
            if (.not.PRESENT(domain) .and. .not.PRESENT(override)) then
               jsdata=1;jedata=len
               jscomp=1;jecomp=len
               gysize = len
               dysize = len
               field(num_fields)%jsc=jscomp;field(num_fields)%jec=jecomp
            elseif (PRESENT(override)) then
               gysize = len
               if (PRESENT(axis_sizes)) axis_sizes(2) = len
            endif
            if(use_comp_domain1) then
               field(num_fields)%siz(2) = ny
            else
               field(num_fields)%siz(2) = dysize
            endif
            if (len /= gysize) then
               write(msg,'(a,"/",a)')  trim(file),trim(fieldname)
               call mpp_error(FATAL,'time_interp_ext, file/field '//trim(msg)//' y dim doesnt match model')
            endif
         case ('Z')
            field(num_fields)%axisname(3) = axisname(j)
            field(num_fields)%siz(3) = len
         case ('T')
            field(num_fields)%axisname(4) = axisname(j)
            field(num_fields)%siz(4) = ntime
            field(num_fields)%tdim   = j
         end select
      enddo
      siz = field(num_fields)%siz
      if(PRESENT(axis_names)) axis_names = field(num_fields)%axisname
      if (PRESENT(axis_sizes) .and. .not.PRESENT(override)) then
         axis_sizes = field(num_fields)%siz
      endif

      if (verb) write(outunit,'(a,4i6)') 'field x,y,z,t local size= ',siz
      if (verb) write(outunit,*) 'field contains data in units = ',trim(field(num_fields)%units)
      if (transpose_xy) call mpp_error(FATAL,'axis ordering not supported')
      if (num_io_buffers .le. 1) call mpp_error(FATAL,'time_interp_ext:num_io_buffers should be at least 2')
      nbuf = min(num_io_buffers,siz(4))

      field(num_fields)%numwindows = numwindows
      allocate(field(num_fields)%need_compute(nbuf, numwindows))
      field(num_fields)%need_compute = .true.

      allocate(field(num_fields)%data(isdata:iedata,jsdata:jedata,siz(3),nbuf),&
               field(num_fields)%mask(isdata:iedata,jsdata:jedata,siz(3),nbuf) )
      field(num_fields)%mask = .false.
      field(num_fields)%data = 0.0
         slope=1.0;intercept=0.0
!             if (units /= 'same') call convert_units(trim(field(num_fields)%units),trim(units),slope,intercept)
!             if (verb.and.units /= 'same') then
!                 write(outunit,*) 'attempting to convert data to units = ',trim(units)
!                 write(outunit,'(a,f8.3,a,f8.3)') 'factor = ',slope,' offset= ',intercept
!             endif
      field(num_fields)%slope = slope
      field(num_fields)%intercept = intercept
      allocate(field(num_fields)%ibuf(nbuf))
      field(num_fields)%ibuf = -1
      field(num_fields)%nbuf =  0 ! initialize buffer number so that first reading fills data(:,:,:,1)
      if(PRESENT(override)) then
         field(num_fields)%is_src = 1
         field(num_fields)%ie_src = gxsize
         field(num_fields)%js_src = 1
         field(num_fields)%je_src = gysize
         allocate(field(num_fields)%src_data(gxsize,gysize,siz(3),nbuf))
      else
         field(num_fields)%is_src = isdata
         field(num_fields)%ie_src = iedata
         field(num_fields)%js_src = jsdata
         field(num_fields)%je_src = jedata
         allocate(field(num_fields)%src_data(isdata:iedata,jsdata:jedata,siz(3),nbuf))
      endif

      allocate(field(num_fields)%time(ntime))
      allocate(field(num_fields)%period(ntime))
      allocate(field(num_fields)%start_time(ntime))
      allocate(field(num_fields)%end_time(ntime))

      do j=1,ntime
         field(num_fields)%time(j)       = get_cal_time(tstamp(j),trim(timeunits),trim(calendar_type),permit_calendar_conversion)
         field(num_fields)%start_time(j) = get_cal_time(tstart(j),trim(timeunits),trim(calendar_type),permit_calendar_conversion)
         field(num_fields)%end_time(j)   = get_cal_time(  tend(j),trim(timeunits),trim(calendar_type),permit_calendar_conversion)
      enddo

      if (field(num_fields)%modulo_time) then
         call set_time_modulo(field(num_fields)%Time)
         call set_time_modulo(field(num_fields)%start_time)
         call set_time_modulo(field(num_fields)%end_time)
      endif

      if(present(correct_leap_year_inconsistency)) then
        field(num_fields)%correct_leap_year_inconsistency = correct_leap_year_inconsistency
      else
        field(num_fields)%correct_leap_year_inconsistency = .false.
      endif

      if(have_modulo_time) then
         if(get_calendar_type() == NO_CALENDAR) then
            field(num_fields)%modulo_time_beg = set_time(timebeg)
            field(num_fields)%modulo_time_end = set_time(timeend)
         else
            field(num_fields)%modulo_time_beg = set_date(timebeg)
            field(num_fields)%modulo_time_end = set_date(timeend)
         endif
         field(num_fields)%have_modulo_times = .true.
      else
         field(num_fields)%have_modulo_times = .false.
      endif
      if(ntime == 1) then
         call mpp_error(NOTE, 'time_interp_external_mod: file '//trim(file)//'  has only one time level')
      else
         do j= 1, ntime
            field(num_fields)%period(j) = field(num_fields)%end_time(j)-field(num_fields)%start_time(j)
            if (field(num_fields)%period(j) > set_time(0,0)) then
               call get_time(field(num_fields)%period(j), sec, day)
               sec = sec/2+mod(day,2)*43200
               day = day/2
               field(num_fields)%time(j) = field(num_fields)%start_time(j)+&
                    set_time(sec,day)
            else
               if (j > 1 .and. j < ntime) then
                  tdiff = field(num_fields)%time(j+1) -  field(num_fields)%time(j-1)
                  call get_time(tdiff, sec, day)
                  sec = sec/2+mod(day,2)*43200
                  day = day/2
                  field(num_fields)%period(j) = set_time(sec,day)
                  sec = sec/2+mod(day,2)*43200
                  day = day/2
                  field(num_fields)%start_time(j) = field(num_fields)%time(j) - set_time(sec,day)
                  field(num_fields)%end_time(j) = field(num_fields)%time(j) + set_time(sec,day)
               elseif ( j == 1) then
                  tdiff = field(num_fields)%time(2) -  field(num_fields)%time(1)
                  call get_time(tdiff, sec, day)
                  field(num_fields)%period(j) = set_time(sec,day)
                  sec = sec/2+mod(day,2)*43200
                  day = day/2
                  field(num_fields)%start_time(j) = field(num_fields)%time(j) - set_time(sec,day)
                  field(num_fields)%end_time(j) = field(num_fields)%time(j) + set_time(sec,day)
               else
                  tdiff = field(num_fields)%time(ntime) -  field(num_fields)%time(ntime-1)
                  call get_time(tdiff, sec, day)
                  field(num_fields)%period(j) = set_time(sec,day)
                  sec = sec/2+mod(day,2)*43200
                  day = day/2
                  field(num_fields)%start_time(j) = field(num_fields)%time(j) - set_time(sec,day)
                  field(num_fields)%end_time(j) = field(num_fields)%time(j) + set_time(sec,day)
               endif
            endif
         enddo
      endif

      do j=1,ntime-1
         if (field(num_fields)%time(j) >= field(num_fields)%time(j+1)) then
            write(msg,'(A,i20)') "times not monotonically increasing. Filename: " &
                 //TRIM(file)//"  field:  "//TRIM(fieldname)//" timeslice: ", j
            call mpp_error(FATAL, TRIM(msg))
         endif
      enddo

      field(num_fields)%modulo_time = get_axis_modulo(fileobj, timename)

      if (verb) then
         if (field(num_fields)%modulo_time) write(outunit,*) 'data are being treated as modulo in time'
         do j= 1, ntime
            write(outunit,*) 'time index,  ', j
            call get_date(field(num_fields)%start_time(j),yr,mon,day,hr,minu,sec)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'start time: yyyy/mm/dd hh:mm:ss= ',yr,'/',mon,'/',day,hr,':',minu,':',sec
            call get_date(field(num_fields)%time(j),yr,mon,day,hr,minu,sec)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'mid time: yyyy/mm/dd hh:mm:ss= ',yr,'/',mon,'/',day,hr,':',minu,':',sec
            call get_date(field(num_fields)%end_time(j),yr,mon,day,hr,minu,sec)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'end time: yyyy/mm/dd hh:mm:ss= ',yr,'/',mon,'/',day,hr,':',minu,':',sec
         enddo
      end if

      deallocate(axisname, axislen)
      deallocate(tstamp, tstart, tend, tavg)

      return

    end function init_external_field

!</FUNCTION> NAME="init_external_field"


    subroutine time_interp_external_2d(index, time, data_in, interp, verbose,horz_interp, mask_out, &
               is_in, ie_in, js_in, je_in, window_id)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real, dimension(:,:), intent(inout) :: data_in
      integer, intent(in), optional :: interp
      logical, intent(in), optional :: verbose
      type(horiz_interp_type),intent(in), optional :: horz_interp
      logical, dimension(:,:), intent(out), optional :: mask_out ! set to true where output data is valid
      integer,                  intent(in), optional :: is_in, ie_in, js_in, je_in
      integer,                  intent(in), optional :: window_id

      real   , dimension(size(data_in,1), size(data_in,2), 1) :: data_out
      logical, dimension(size(data_in,1), size(data_in,2), 1) :: mask3d

      data_out(:,:,1) = data_in(:,:) ! fill initial values for the portions of array that are not touched by 3d routine
      call time_interp_external_3d(index, time, data_out, interp, verbose, horz_interp, mask3d, &
                                   is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      data_in(:,:) = data_out(:,:,1)
      if (PRESENT(mask_out)) mask_out(:,:) = mask3d(:,:,1)

      return
    end subroutine time_interp_external_2d


    function get_external_fileobj(filename, fileobj)
       character(len=*),             intent(in) :: filename
       type(FmsNetcdfFile_t), intent(out) :: fileobj
       logical                                  :: get_external_fileobj
       integer :: i

       get_external_fileobj = .false.
       do i=1,num_files
          if(trim(opened_files(i)%filename) == trim(filename)) then
            fileobj = opened_files(i)%fileobj
            get_external_fileobj = .true.
            exit  ! file is already opened
         endif
      enddo

    end function get_external_fileobj


!<SUBROUTINE NAME="time_interp_external" >
!
!<DESCRIPTION>
! Provide data from external file interpolated to current model time.
! Data may be local to current processor or global, depending on
! "init_external_field" flags.
!</DESCRIPTION>
!
!<IN NAME="index" TYPE="integer">
! index of external field from previous call to init_external_field
!</IN>
!<IN NAME="time" TYPE="time_manager_mod:time_type">
! target time for data
!</IN>
!<INOUT NAME="data" TYPE="real" DIM="(:,:),(:,:,:)">
! global or local data array
!</INOUT>
!<IN NAME="interp" TYPE="integer">
! time_interp_external defined interpolation method (optional).  Currently this module only supports
! LINEAR_TIME_INTERP.
!</IN>
!<IN NAME="verbose" TYPE="logical">
! verbose flag for debugging (optional).
!</IN>

    subroutine time_interp_external_3d(index, time, data, interp,verbose,horz_interp, mask_out, is_in, ie_in, js_in, je_in, window_id)

      integer,                    intent(in)           :: index
      type(time_type),            intent(in)           :: time
      real, dimension(:,:,:),  intent(inout)           :: data
      integer,                    intent(in), optional :: interp
      logical,                    intent(in), optional :: verbose
      type(horiz_interp_type),    intent(in), optional :: horz_interp
      logical, dimension(:,:,:), intent(out), optional :: mask_out ! set to true where output data is valid
      integer,                    intent(in), optional :: is_in, ie_in, js_in, je_in
      integer,                    intent(in), optional :: window_id

      integer :: nx, ny, nz, interp_method, t1, t2
      integer :: i1, i2, isc, iec, jsc, jec, mod_time
      integer :: yy, mm, dd, hh, min, ss
      character(len=256) :: err_msg, filename

      integer :: isw, iew, jsw, jew, nxw, nyw
          ! these are boundaries of the updated portion of the "data" argument
          ! they are calculated using sizes of the "data" and isc,iec,jsc,jsc
          ! fileds from respective input field, to center the updated portion
          ! in the output array

      real :: w1,w2
      logical :: verb
      character(len=16) :: message1, message2

      nx = size(data,1)
      ny = size(data,2)
      nz = size(data,3)

      interp_method = LINEAR_TIME_INTERP
      if (PRESENT(interp)) interp_method = interp
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.

      if (index < 1.or.index > num_fields) &
           call mpp_error(FATAL,'invalid index in call to time_interp_ext -- field was not initialized or failed to initialize')

      isc=field(index)%isc;iec=field(index)%iec
      jsc=field(index)%jsc;jec=field(index)%jec

      if( field(index)%numwindows == 1 ) then
         nxw = iec-isc+1
         nyw = jec-jsc+1
      else
         if( .not. present(is_in) .or. .not. present(ie_in) .or. .not. present(js_in) .or. .not. present(je_in) ) then
            call mpp_error(FATAL, 'time_interp_external: is_in, ie_in, js_in and je_in must be present ' // &
                                  'when numwindows > 1, field='//trim(field(index)%name))
         endif
         nxw = ie_in - is_in + 1
         nyw = je_in - js_in + 1
         isc = isc + is_in - 1
         iec = isc + ie_in - is_in
         jsc = jsc + js_in - 1
         jec = jsc + je_in - js_in
      endif

      isw = (nx-nxw)/2+1; iew = isw+nxw-1
      jsw = (ny-nyw)/2+1; jew = jsw+nyw-1

      if (nx < nxw .or. ny < nyw .or. nz < field(index)%siz(3)) then
         write(message1,'(i6,2i5)') nx,ny,nz
         call mpp_error(FATAL,'field '//trim(field(index)%name)//' Array size mismatch in time_interp_external.'// &
         ' Array "data" is too small. shape(data)='//message1)
      endif
      if(PRESENT(mask_out)) then
        if (size(mask_out,1) /= nx .or. size(mask_out,2) /= ny .or. size(mask_out,3) /= nz) then
          write(message1,'(i6,2i5)') nx,ny,nz
          write(message2,'(i6,2i5)') size(mask_out,1),size(mask_out,2),size(mask_out,3)
          call mpp_error(FATAL,'field '//trim(field(index)%name)//' array size mismatch in time_interp_external.'// &
          ' Shape of array "mask_out" does not match that of array "data".'// &
          ' shape(data)='//message1//' shape(mask_out)='//message2)
        endif
      endif

      if (field(index)%siz(4) == 1) then
         ! only one record in the file => time-independent field
         call load_record(field(index),1,horz_interp, is_in, ie_in ,js_in, je_in,window_id)
         i1 = find_buf_index(1,field(index)%ibuf)
         if( field(index)%region_type == NO_REGION ) then
            where(field(index)%mask(isc:iec,jsc:jec,:,i1))
               data(isw:iew,jsw:jew,:) = field(index)%data(isc:iec,jsc:jec,:,i1)
            elsewhere
!               data(isw:iew,jsw:jew,:) = time_interp_missing !field(index)%missing? Balaji
               data(isw:iew,jsw:jew,:) = field(index)%missing
            end where
         else
            where(field(index)%mask(isc:iec,jsc:jec,:,i1))
               data(isw:iew,jsw:jew,:) = field(index)%data(isc:iec,jsc:jec,:,i1)
            end where
         endif
         if(PRESENT(mask_out)) &
              mask_out(isw:iew,jsw:jew,:) = field(index)%mask(isc:iec,jsc:jec,:,i1)
      else
        if(field(index)%have_modulo_times) then
          call time_interp(time,field(index)%modulo_time_beg, field(index)%modulo_time_end, field(index)%time(:), &
                          w2, t1, t2, field(index)%correct_leap_year_inconsistency, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(field(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 1: "//trim(err_msg)//&
                    ",file="//trim(filename)//",field="//trim(field(index)%name) )
          endif
        else
          if(field(index)%modulo_time) then
            mod_time=1
          else
            mod_time=0
          endif
          call time_interp(time,field(index)%time(:),w2,t1,t2,modtime=mod_time, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(field(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 2: "//trim(err_msg)//&
                    ",file="//trim(filename)//",field="//trim(field(index)%name) )
          endif
        endif
         w1 = 1.0-w2
         if (verb) then
            call get_date(time,yy,mm,dd,hh,min,ss)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',min,':',ss
            write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
         endif

         call load_record(field(index),t1,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
         call load_record(field(index),t2,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
         i1 = find_buf_index(t1,field(index)%ibuf)
         i2 = find_buf_index(t2,field(index)%ibuf)
         if(i1<0.or.i2<0) &
              call mpp_error(FATAL,'time_interp_external : records were not loaded correctly in memory')

         if (verb) then
            write(outunit,*) 'ibuf= ',field(index)%ibuf
            write(outunit,*) 'i1,i2= ',i1, i2
         endif

         if( field(index)%region_type == NO_REGION ) then
            where(field(index)%mask(isc:iec,jsc:jec,:,i1).and.field(index)%mask(isc:iec,jsc:jec,:,i2))
               data(isw:iew,jsw:jew,:) = field(index)%data(isc:iec,jsc:jec,:,i1)*w1 + &
                    field(index)%data(isc:iec,jsc:jec,:,i2)*w2
            elsewhere
!               data(isw:iew,jsw:jew,:) = time_interp_missing !field(index)%missing? Balaji
               data(isw:iew,jsw:jew,:) = field(index)%missing
            end where
         else
            where(field(index)%mask(isc:iec,jsc:jec,:,i1).and.field(index)%mask(isc:iec,jsc:jec,:,i2))
               data(isw:iew,jsw:jew,:) = field(index)%data(isc:iec,jsc:jec,:,i1)*w1 + &
                    field(index)%data(isc:iec,jsc:jec,:,i2)*w2
            end where
         endif
         if(PRESENT(mask_out)) &
              mask_out(isw:iew,jsw:jew,:) = &
                                        field(index)%mask(isc:iec,jsc:jec,:,i1).and.&
                                        field(index)%mask(isc:iec,jsc:jec,:,i2)
      endif

    end subroutine time_interp_external_3d
!</SUBROUTINE> NAME="time_interp_external"

    subroutine time_interp_external_0d(index, time, data, verbose)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real, intent(inout) :: data
      logical, intent(in), optional :: verbose

      integer :: t1, t2
      integer :: i1, i2, mod_time
      integer :: yy, mm, dd, hh, min, ss
      character(len=256) :: err_msg, filename

      real :: w1,w2
      logical :: verb

      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.

      if (index < 1.or.index > num_fields) &
           call mpp_error(FATAL,'invalid index in call to time_interp_ext -- field was not initialized or failed to initialize')

      if (field(index)%siz(4) == 1) then
         ! only one record in the file => time-independent field
         call load_record_0d(field(index),1)
         i1 = find_buf_index(1,field(index)%ibuf)
         data = field(index)%data(1,1,1,i1)
      else
        if(field(index)%have_modulo_times) then
          call time_interp(time,field(index)%modulo_time_beg, field(index)%modulo_time_end, field(index)%time(:), &
                          w2, t1, t2, field(index)%correct_leap_year_inconsistency, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(field(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 3:"//trim(err_msg)//&
                    ",file="//trim(filename)//",field="//trim(field(index)%name) )
          endif
        else
          if(field(index)%modulo_time) then
            mod_time=1
          else
            mod_time=0
          endif
          call time_interp(time,field(index)%time(:),w2,t1,t2,modtime=mod_time, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(field(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 4:"//trim(err_msg)// &
                    ",file="//trim(filename)//",field="//trim(field(index)%name) )
          endif
        endif
         w1 = 1.0-w2
         if (verb) then
            call get_date(time,yy,mm,dd,hh,min,ss)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',min,':',ss
            write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
         endif
         call load_record_0d(field(index),t1)
         call load_record_0d(field(index),t2)
         i1 = find_buf_index(t1,field(index)%ibuf)
         i2 = find_buf_index(t2,field(index)%ibuf)

         if(i1<0.or.i2<0) &
              call mpp_error(FATAL,'time_interp_external : records were not loaded correctly in memory')
         data = field(index)%data(1,1,1,i1)*w1 + field(index)%data(1,1,1,i2)*w2
         if (verb) then
            write(outunit,*) 'ibuf= ',field(index)%ibuf
            write(outunit,*) 'i1,i2= ',i1, i2
         endif
      endif

    end subroutine time_interp_external_0d

    subroutine set_time_modulo(Time)

      type(time_type), intent(inout), dimension(:) :: Time

      integer :: ntime, n
      integer :: yr, mon, dy, hr, minu, sec

      ntime = size(Time(:))

      do n = 1, ntime
         call get_date(Time(n), yr, mon, dy, hr, minu, sec)
         yr = modulo_year
         Time(n) = set_date(yr, mon, dy, hr, minu, sec)
      enddo


    end subroutine set_time_modulo

! ============================================================================
! load specified record from file
subroutine load_record(field, rec, interp, is_in, ie_in, js_in, je_in, window_id_in)
  type(ext_fieldtype),     intent(inout)        :: field
  integer            ,     intent(in)           :: rec    ! record number
  type(horiz_interp_type), intent(in), optional :: interp
  integer,                 intent(in), optional :: is_in, ie_in, js_in, je_in
  integer,                 intent(in), optional :: window_id_in

  ! ---- local vars
  integer :: ib ! index in the array of input buffers
  integer :: isw,iew,jsw,jew ! boundaries of the domain on each window
  integer :: is_region, ie_region, js_region, je_region, i, j, n
  integer :: start(4), nread(4)
  logical :: need_compute
  real    :: mask_in(size(field%src_data,1),size(field%src_data,2),size(field%src_data,3))
  real, allocatable :: mask_out(:,:,:)
  integer :: window_id

  window_id = 1
  if( PRESENT(window_id_in) ) window_id = window_id_in
  need_compute = .true.

!$OMP CRITICAL
  ib = find_buf_index(rec,field%ibuf)

  if(ib>0) then
     !--- do nothing
     need_compute = .false.
  else
     ! calculate current buffer number in round-robin fasion
     field%nbuf = field%nbuf + 1
     if(field%nbuf > size(field%data,4).or.field%nbuf <= 0) field%nbuf = 1
     ib = field%nbuf
     field%ibuf(ib) = rec
     field%need_compute(ib,:) = .true.

     if (debug_this_module) write(outunit,*) 'reading record without domain for field ',trim(field%name)
     start = 1; nread = 1
     start(1) = field%is_src; nread(1) = field%ie_src - field%is_src + 1
     start(2) = field%js_src; nread(2) = field%je_src - field%js_src + 1
     start(3) = 1;            nread(3) = size(field%src_data,3)
     start(field%tdim) = rec; nread(field%tdim) = 1
     call read_data(field%fileobj,field%name,field%src_data(:,:,:,ib),corner=start,edge_lengths=nread)
  endif
!$OMP END CRITICAL
  isw=field%isc;iew=field%iec
  jsw=field%jsc;jew=field%jec

  if( field%numwindows > 1) then
     if( .NOT. PRESENT(is_in) .OR. .NOT. PRESENT(ie_in) .OR. .NOT. PRESENT(js_in) .OR. .NOT. PRESENT(je_in) ) then
        call mpp_error(FATAL, 'time_interp_external(load_record): is_in, ie_in, js_in, je_in must be present when numwindows>1')
     endif
     isw = isw + is_in - 1
     iew = isw + ie_in - is_in
     jsw = jsw + js_in - 1
     jew = jsw + je_in - js_in
  endif

  ! interpolate to target grid

  need_compute = field%need_compute(ib, window_id)
  if(need_compute) then
     if(PRESENT(interp)) then
        is_region = field%is_region; ie_region = field%ie_region
        js_region = field%js_region; je_region = field%je_region
        mask_in = 0.0
        where (is_valid(field%src_data(:,:,:,ib), field%valid)) mask_in = 1.0
        if ( field%region_type .NE. NO_REGION ) then
           if( ANY(mask_in == 0.0) ) then
              call mpp_error(FATAL, "time_interp_external: mask_in should be all 1 when region_type is not NO_REGION")
           endif
           if( field%region_type == OUTSIDE_REGION) then
              do j = js_region, je_region
                 do i = is_region, ie_region
                    mask_in(i,j,:) = 0.0
                 enddo
              enddo
           else  ! field%region_choice == INSIDE_REGION
              do j = 1, size(mask_in,2)
                 do i = 1, size(mask_in,1)
                    if( j<js_region .OR. j>je_region .OR. i<is_region .OR. i>ie_region ) mask_in(i,j,:) = 0.0
                 enddo
              enddo
           endif
        endif
        allocate(mask_out(isw:iew,jsw:jew, size(field%src_data,3)))
        call horiz_interp(interp,field%src_data(:,:,:,ib),field%data(isw:iew,jsw:jew,:,ib), &
             mask_in=mask_in, &
             mask_out=mask_out)

        field%mask(isw:iew,jsw:jew,:,ib) = mask_out(isw:iew,jsw:jew,:) > 0
        deallocate(mask_out)
        field%need_compute(ib, window_id) = .false.
     else
        if ( field%region_type .NE. NO_REGION ) then
           call mpp_error(FATAL, "time_interp_external: region_type should be NO_REGION when interp is not present")
        endif
        field%data(isw:iew,jsw:jew,:,ib) = field%src_data(isw:iew,jsw:jew,:,ib)
        field%mask(isw:iew,jsw:jew,:,ib) = is_valid(field%data(isw:iew,jsw:jew,:,ib),field%valid)
     endif
     ! convert units
     where(field%mask(isw:iew,jsw:jew,:,ib)) field%data(isw:iew,jsw:jew,:,ib) = &
          field%data(isw:iew,jsw:jew,:,ib)*field%slope + field%intercept
  endif

end subroutine load_record


subroutine load_record_0d(field, rec)
  type(ext_fieldtype),     intent(inout)        :: field
  integer            ,     intent(in)           :: rec    ! record number
  ! ---- local vars
  integer :: ib ! index in the array of input buffers
  integer :: start(4), nread(4)

  ib = find_buf_index(rec,field%ibuf)

  if(ib>0) then
     return
  else
     ! calculate current buffer number in round-robin fasion
     field%nbuf = field%nbuf + 1
     if(field%nbuf > size(field%data,4).or.field%nbuf <= 0) field%nbuf = 1
     ib = field%nbuf
     field%ibuf(ib) = rec

     if (debug_this_module) write(outunit,*) 'reading record without domain for field ',trim(field%name)
     start = 1; nread = 1
     start(3) = 1;            nread(3) = size(field%src_data,3)
     start(field%tdim) = rec; nread(field%tdim) = 1
     call read_data(field%fileobj,field%name,field%src_data(:,:,:,ib),corner=start,edge_lengths=nread)
     if ( field%region_type .NE. NO_REGION ) then
        call mpp_error(FATAL, "time_interp_external: region_type should be NO_REGION when field is scalar")
     endif
     field%data(1,1,:,ib) = field%src_data(1,1,:,ib)
     field%mask(1,1,:,ib) = is_valid(field%data(1,1,:,ib),field%valid)
     ! convert units
     where(field%mask(1,1,:,ib)) field%data(1,1,:,ib) = &
          field%data(1,1,:,ib)*field%slope + field%intercept
  endif

end subroutine load_record_0d

! ============================================================================
subroutine reset_src_data_region(index, is, ie, js, je)
   integer, intent(in) :: index
   integer, intent(in) :: is, ie, js, je
   integer             :: nk, nbuf

   if( is == field(index)%is_src .AND. ie == field(index)%ie_src .AND. &
       js == field(index)%js_src .AND. ie == field(index)%je_src ) return

   if( .NOT. ASSOCIATED(field(index)%src_data) ) call mpp_error(FATAL, &
       "time_interp_external: field(index)%src_data is not associated")
   nk = size(field(index)%src_data,3)
   nbuf = size(field(index)%src_data,4)
   deallocate(field(index)%src_data)
   allocate(field(index)%src_data(is:ie,js:je,nk,nbuf))
   field(index)%is_src = is
   field(index)%ie_src = ie
   field(index)%js_src = js
   field(index)%je_src = je


end subroutine reset_src_data_region

! ============================================================================
subroutine set_override_region(index, region_type, is_region, ie_region, js_region, je_region)
   integer, intent(in) :: index, region_type
   integer, intent(in) :: is_region, ie_region, js_region, je_region

   field(index)%region_type = region_type
   field(index)%is_region   = is_region
   field(index)%ie_region   = ie_region
   field(index)%js_region   = js_region
   field(index)%je_region   = je_region

   return

end subroutine set_override_region

! ============================================================================
! reallocates array of fields, increasing its size
subroutine realloc_files(n)
  integer, intent(in) :: n ! new size

  type(filetype), pointer :: ptr(:)
  integer :: i

  if (associated(opened_files)) then
     if (n <= size(opened_files)) return ! do nothing, if requested size no more than current
  endif

  allocate(ptr(n))
  do i = 1, size(ptr)
     ptr(i)%filename = ''
     if(Associated(ptr(i)%fileobj)) then
        call close_file(ptr(i)%fileobj)
        deallocate(ptr(i)%fileobj)
     endif
  enddo

  if (associated(opened_files))then
     ptr(1:size(opened_files)) = opened_files(:)
     deallocate(opened_files)
  endif
  opened_files => ptr

end subroutine realloc_files

! ============================================================================
! reallocates array of fields,increasing its size
subroutine realloc_fields(n)
  integer, intent(in) :: n ! new size

  type(ext_fieldtype), pointer :: ptr(:)
  integer :: i, ier

  if (associated(field)) then
     if (n <= size(field)) return ! do nothing if requested size no more then current
  endif

  allocate(ptr(n))
  do i=1,size(ptr)
     ptr(i)%fileobj => NULL()
     ptr(i)%name=''
     ptr(i)%units=''
     ptr(i)%siz=-1
     ptr(i)%ndim=-1
     ptr(i)%domain = NULL_DOMAIN2D
     if (ASSOCIATED(ptr(i)%time))       DEALLOCATE(ptr(i)%time, stat=ier)
     if (ASSOCIATED(ptr(i)%start_time)) DEALLOCATE(ptr(i)%start_time, stat=ier)
     if (ASSOCIATED(ptr(i)%end_time))   DEALLOCATE(ptr(i)%end_time, stat=ier)
     if (ASSOCIATED(ptr(i)%period)) DEALLOCATE(ptr(i)%period, stat=ier)
     ptr(i)%modulo_time=.false.
     if (ASSOCIATED(ptr(i)%data)) DEALLOCATE(ptr(i)%data, stat=ier)
     if (ASSOCIATED(ptr(i)%ibuf)) DEALLOCATE(ptr(i)%ibuf, stat=ier)
     if (ASSOCIATED(ptr(i)%src_data)) DEALLOCATE(ptr(i)%src_data, stat=ier)
     ptr(i)%nbuf=-1
     ptr(i)%domain_present=.false.
     ptr(i)%slope=1.0
     ptr(i)%intercept=0.0
     ptr(i)%isc=-1;ptr(i)%iec=-1
     ptr(i)%jsc=-1;ptr(i)%jec=-1
  enddo
  if (associated(field)) then
     ptr(1:size(field)) = field(:)
     deallocate(field)
  endif
  field=>ptr

end subroutine realloc_fields


    function find_buf_index(indx,buf)
      integer :: indx
      integer, dimension(:) :: buf
      integer :: find_buf_index

      integer :: nbuf, i

      nbuf = size(buf(:))

      find_buf_index = -1

      do i=1,nbuf
         if (buf(i) == indx) then
            find_buf_index = i
            exit
         endif
      enddo

    end function find_buf_index

!<FUNCTION NAME="get_external_field_size" TYPE="integer" DIM="(4)">
!
!<DESCRIPTION>
! return size of field after call to init_external_field.
! Ordering is X/Y/Z/T.
! This call only makes sense for non-distributed reads.
!</DESCRIPTION>
!
!<IN NAME="index" TYPE="integer">
! returned from previous call to init_external_field.
!</IN>

    function get_external_field_size(index)

      integer :: index
      integer :: get_external_field_size(4)

      if (index .lt. 1 .or. index .gt. num_fields) &
           call mpp_error(FATAL,'invalid index in call to get_external_field_size')


      get_external_field_size(1) = field(index)%siz(1)
      get_external_field_size(2) = field(index)%siz(2)
      get_external_field_size(3) = field(index)%siz(3)
      get_external_field_size(4) = field(index)%siz(4)

    end function get_external_field_size
!</FUNCTION> NAME="get_external_field_size"


!<FUNCTION NAME="get_external_field_missing" TYPE="real">
!
!<DESCRIPTION>
! return missing value
!</DESCRIPTION>
!
!<IN NAME="index" TYPE="integer">
! returned from previous call to init_external_field.
!</IN>

    function get_external_field_missing(index)

      integer :: index
      real :: missing
      real :: get_external_field_missing

      if (index .lt. 1 .or. index .gt. num_fields) &
           call mpp_error(FATAL,'invalid index in call to get_external_field_size')


      get_external_field_missing = field(index)%missing

    end function get_external_field_missing
!</FUNCTION> NAME="get_external_field_missing"

! ===========================================================================
subroutine get_time_axis(index, time)
  integer        , intent(in)  :: index   ! field id
  type(time_type), intent(out) :: time(:) ! array of time values to be filled

  integer :: n ! size of the data to be assigned

  if (index < 1.or.index > num_fields) &
       call mpp_error(FATAL,'invalid index in call to get_time_axis')

  n = min(size(time),size(field(index)%time))

  time(1:n) = field(index)%time(1:n)
end subroutine

!<SUBROUTINE NAME="time_interp_external_exit">
!
!<DESCRIPTION>
! exit time_interp_external_mod.  Close all open files and
! release storage
!</DESCRIPTION>

    subroutine time_interp_external_exit()

      integer :: i,j
!
! release storage arrays
!
      do i=1,num_fields
         deallocate(field(i)%time,field(i)%start_time,field(i)%end_time,&
              field(i)%period,field(i)%data,field(i)%mask,field(i)%ibuf)
         if (ASSOCIATED(field(i)%src_data)) deallocate(field(i)%src_data)
         field(i)%domain = NULL_DOMAIN2D
         field(i)%nbuf = 0
         field(i)%slope = 0.
         field(i)%intercept = 0.
      enddo

      !-- close all the files opended
      do i = 1, num_files
         call close_file(opened_files(i)%fileobj)
         deallocate(opened_files(i)%fileobj)
      enddo

      deallocate(field)
      deallocate(opened_files)

      num_fields = 0

      module_initialized = .false.

    end subroutine time_interp_external_exit
!</SUBROUTINE> NAME="time_interp_external_exit"

end module time_interp_external2_mod
