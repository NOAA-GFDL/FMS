#include  <os.h>

module time_interp_external_mod
!
!<CONTACT EMAIL="mjh@gfdl.noaa.gov">M.J. Harrison</CONTACT>
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

  use mpp_mod, only : mpp_error,FATAL,mpp_pe, stdout, stdlog
  use mpp_io_mod, only : mpp_open, mpp_get_atts, mpp_get_info, MPP_NETCDF, MPP_MULTI, MPP_SINGLE,&
       mpp_get_times, MPP_RDONLY, MPP_ASCII, default_axis,axistype,fieldtype,atttype, &
       mpp_get_axes, mpp_get_fields, mpp_read, default_field, mpp_close
  use time_manager_mod, only : time_type, get_date, set_date, operator ( >= ) , operator ( + ) , days_in_month, &
                              operator ( / ) , days_in_year, increment_time
  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain, NULL_DOMAIN2D
!  use time_interp_mod, only : get_time_weights, LINEAR_TIME_INTERP, modulo_year
  use time_interp_mod, only : time_interp
  use udunits_mod, only : get_cal_time, convert_units 
  use axis_utils_mod, only : get_axis_cart, get_axis_modulo
  use fms_mod, only : lowercase

  implicit none

  character(len=128), private :: version='CVS $Id: time_interp_external.F90,v 1.2 2002/07/16 22:57:03 fms Exp $'
  character(len=128), private :: tagname='Tag $Name: havana $'

  integer, parameter, private :: max_fields = 100, modulo_year= 0001
  integer, parameter, private :: LINEAR_TIME_INTERP = 1 ! not used currently
  integer, private :: num_fields = 0
  ! denotes time intervals in file (interpreted from metadata)
  integer, parameter, private :: FIELD_UNDEF = -1, FIELD_ANNUAL = 1, FIELD_MONTHLY = 2, FIELD_DAILY = 3, FIELD_HOURLY = 4
  integer, private :: num_io_buffers = -1 ! set -1 to read all records from disk into memory 
  logical, private :: time_interp_ext_initialized = .false.

  public init_external_field, time_interp_external, time_interp_external_init, time_interp_external_exit, get_external_field_size

  private set_monthly_time, set_annual_time, set_daily_time, set_hourly_time, find_buf_index

  type, private :: ext_fieldtype
     integer :: unit ! keep unit open when not reading all records
     character(len=128) :: name, units
     integer :: siz(4), tavg_type, ndim
     type(domain2d) :: domain
     type(axistype) :: axes(4)
     type(time_type), dimension(:), pointer :: time ! midpoint of time interval
     type(time_type), dimension(:), pointer :: start_time, end_time
     type(fieldtype) :: field ! mpp_io type
     real, pointer, dimension(:) :: period ! unit of days
     logical :: modulo_time ! denote climatological time axis
     real, dimension(:,:,:,:), pointer :: data ! defined over data domain or global domain
     integer, dimension(:), pointer :: ibuf
     real, dimension(:,:), pointer :: buf2d
     real, dimension(:,:,:), pointer :: buf3d
     integer :: nbuf
     logical :: domain_present
     real(DOUBLE_KIND) :: slope, intercept
     integer :: isc,iec,jsc,jec
  end type ext_fieldtype

  interface time_interp_external
     module procedure time_interp_external_2d
     module procedure time_interp_external_3d
  end interface

  type(ext_fieldtype), private :: field(max_fields)
  
  contains

! <SUBROUTINE NAME="time_interp_external_init">
!
! <DESCRIPTION>
! Initialize the time_interp_external module
! </DESCRIPTION>
!
    subroutine time_interp_external_init()

      integer :: ioun, io_status, i

      namelist /time_interp_external_nml/ num_io_buffers

      ! open and read namelist

      write(stdlog(),'(/a/)') version
      write(stdlog(),'(/a/)') tagname

      call mpp_open(ioun,'input.nml',action=MPP_RDONLY,form=MPP_ASCII)
      read(ioun,time_interp_external_nml,iostat=io_status)
      write(stdlog(),time_interp_external_nml)
      if (io_status .gt. 0) then
         call mpp_error(FATAL,'=>Error reading time_interp_external_nml')
      endif
      call mpp_close(ioun)

      do i=1,max_fields
         field(i)%unit=-1
         field(i)%name=''
         field(i)%units=''
         field(i)%siz=-1
         field(i)%tavg_type=-1
         field(i)%ndim=-1
         field(i)%domain = NULL_DOMAIN2D
         field(i)%axes(:) = default_axis
         if (ASSOCIATED(field(i)%time)) NULLIFY(field(i)%time)
         if (ASSOCIATED(field(i)%start_time)) NULLIFY(field(i)%start_time)
         if (ASSOCIATED(field(i)%end_time)) NULLIFY(field(i)%end_time)
         field(i)%field = default_field
         if (ASSOCIATED(field(i)%period)) NULLIFY(field(i)%period)
         field(i)%modulo_time=.false.
         if (ASSOCIATED(field(i)%data)) NULLIFY(field(i)%data)
         if (ASSOCIATED(field(i)%ibuf)) NULLIFY(field(i)%ibuf)
         if (ASSOCIATED(field(i)%buf2d)) NULLIFY(field(i)%buf2d)
         if (ASSOCIATED(field(i)%buf3d)) NULLIFY(field(i)%buf3d)
         field(i)%nbuf=-1
         field(i)%domain_present=.false.
         field(i)%slope=1.0
         field(i)%intercept=0.0
         field(i)%isc=-1;field(i)%iec=-1
         field(i)%jsc=-1;field(i)%jec=-1
      enddo

      time_interp_ext_initialized = .true.


      return
      
    end  subroutine time_interp_external_init
! </SUBROUTINE> NAME="time_interp_external_init"


!<FUNCTION NAME="init_external_field" TYPE="integer">
!
!<DESCRIPTION>
! initialize an external field.  Buffer entire field to memory (default) or
! store "num_io_buffers" in memory to reduce memory allocations. 
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
!<INOUT NAME="axis_centers" TYPE="axistype" DIM=(4)>
! MPP_IO axistype array for grid centers ordered X-Y-Z-T (optional).
!</INOUT>
!<INOUT NAME="axis_sizes" TYPE="integer" DIM=(4)>
!  array of axis lengths ordered X-Y-Z-T (optional).
!</INOUT>


    function init_external_field(file,fieldname,format,threading,domain,desired_units,verbose,axis_centers,axis_sizes)
      
      character(len=*), intent(in) :: file,fieldname
      integer, intent(in), optional :: format, threading
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: desired_units
      type(domain2d), intent(in), optional :: domain
      type(axistype), intent(inout), optional :: axis_centers(4)
      integer, intent(inout), optional :: axis_sizes(4)

      integer :: init_external_field
      
      type(fieldtype), dimension(:), allocatable :: flds
      type(axistype), dimension(:), allocatable :: axes, fld_axes
      type(axistype) :: time_axis
      type(atttype), allocatable, dimension(:) :: global_atts, atts
      
      real(DOUBLE_KIND) :: slope, intercept
      integer :: form, thread, fset, unit,ndim,nvar,natt,ntime,i,j
      integer :: iscomp,iecomp,jscomp,jecomp,isglobal,ieglobal,jsglobal,jeglobal
      integer :: isdata,iedata,jsdata,jedata, dxsize, dysize,dxsize_max,dysize_max
      logical :: verb, transpose_xy
      real, dimension(:), allocatable :: tstamp
      character(len=1) :: cart
      character(len=128) :: units, fld_units
      character(len=128) :: name
      integer :: siz(4), siz_in(4), gxsize, gysize,gxsize_max, gysize_max
      type(time_type) :: tdiff
      integer :: yr, mon, day, hr, minu, sec
      integer :: yr2, mon2, day2, hr2, min2, sec2
      integer :: len, natts

      if (.not. time_interp_ext_initialized) call mpp_error(FATAL,'Must call time_interp_external_init first')

      form=MPP_NETCDF
      if (PRESENT(format)) form = format
      thread = MPP_MULTI
      if (PRESENT(threading)) thread = threading
      fset = MPP_SINGLE
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      units = 'same'
      if (PRESENT(desired_units)) units = desired_units
      
      call mpp_open(unit,trim(file),MPP_RDONLY,form,threading=thread,&
           fileset=fset)
      call mpp_get_info(unit,ndim,nvar,natt,ntime)
      allocate(global_atts(natt))
      call mpp_get_atts(unit, global_atts)
      allocate(axes(ndim))
      call mpp_get_axes(unit, axes, time_axis)
      allocate(flds(nvar))
      call mpp_get_fields(unit,flds)
      allocate(tstamp(ntime))
      call mpp_get_times(unit,tstamp)

      transpose_xy = .false.

      if (PRESENT(domain)) then
          call mpp_get_compute_domain(domain,iscomp,iecomp,jscomp,jecomp)
          call mpp_get_data_domain(domain,isdata,iedata,jsdata,jedata,dxsize,dxsize_max,dysize,dysize_max)
          call mpp_get_global_domain(domain,isglobal,ieglobal,jsglobal,jeglobal,gxsize,gxsize_max,gysize,gysize_max)
      endif
      
      init_external_field = -1
      
      do i=1,nvar
         call mpp_get_atts(flds(i),name=name,units=fld_units,ndim=ndim,siz=siz_in)
         if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
             if (verb) write(stdout(),*) 'found field ',trim(fieldname), ' in file !!'
             num_fields = num_fields + 1
             init_external_field = num_fields
             field(num_fields)%unit = unit
             field(num_fields)%name = trim(name)
             field(num_fields)%units = trim(fld_units)
             field(num_fields)%field = flds(i)
             if (PRESENT(domain)) then
                field(num_fields)%domain_present = .true.
                field(num_fields)%domain = domain
                field(num_fields)%isc=iscomp;field(num_fields)%iec = iecomp
                field(num_fields)%jsc=jscomp;field(num_fields)%jec = jecomp
             else
                field(num_fields)%domain_present = .false.
             endif
             allocate(fld_axes(ndim))
             call mpp_get_atts(flds(i),axes=fld_axes)
             if (ndim < 2 .or. ndim > 4) call mpp_error(FATAL, &
                  'invalid array rank only 2-4d fields supported')
             field(num_fields)%siz = 1
             field(num_fields)%ndim = ndim
             do j=1,field(num_fields)%ndim
                cart = 'N'
                call get_axis_cart(fld_axes(j), cart)
                call mpp_get_atts(fld_axes(j),len=len)
                if (cart == 'N') call mpp_error(FATAL,'couldnt recognize axis atts in time_interp_external')
                select case (cart)
                case ('X')
                    if (j.eq.2) transpose_xy = .true.
                    if (.not.PRESENT(domain)) then
                       isdata=1;iedata=len
                       iscomp=1;iecomp=len
                       gxsize = len
                       dxsize = len
                       field(num_fields)%isc=iscomp;field(num_fields)%iec=iecomp
                    endif
                    field(num_fields)%axes(1) = fld_axes(j)
                    field(num_fields)%siz(1) = dxsize
                    if (len /= gxsize) call mpp_error(FATAL,'x dim doesnt match model')
                case ('Y')
                    field(num_fields)%axes(2) = fld_axes(j)
                    if (.not.PRESENT(domain)) then
                       jsdata=1;jedata=len
                       jscomp=1;jecomp=len
                       gysize = len
                       dysize = len
                       field(num_fields)%jsc=jscomp;field(num_fields)%jec=jecomp
                    endif
                    field(num_fields)%siz(2) = dysize
                    if (len /= gysize) call mpp_error(FATAL,'y dim doesnt match model')                        
                case ('Z')
                    field(num_fields)%axes(3) = fld_axes(j)
                    field(num_fields)%siz(3) = siz_in(3)
                case ('T')
                    field(num_fields)%axes(4) = fld_axes(j)
                    field(num_fields)%siz(4) = ntime
                end select
             enddo
             siz = field(num_fields)%siz

             if (PRESENT(axis_centers)) then
                axis_centers = field(num_fields)%axes
             endif

             if (PRESENT(axis_sizes)) then
                axis_sizes = field(num_fields)%siz
             endif

             deallocate(fld_axes)
             if (verb) write(stdout(),'(a,4i6)') 'field x,y,z,t local size= ',siz
             if (verb) write(stdout(),*) 'field contains data in units = ',trim(field(num_fields)%units)
             if (transpose_xy) call mpp_error(FATAL,'axis ordering not supported')
             if (num_io_buffers == -1) then
                field(num_fields)%nbuf = siz(4)
                allocate(field(num_fields)%data(isdata:iedata,jsdata:jedata,siz(3),siz(4)))
             else
                if (num_io_buffers .le. 1) call mpp_error(FATAL,'num_io_buffers should be at least 2')
                field(num_fields)%nbuf = min(num_io_buffers,siz(4))
                allocate(field(num_fields)%data(isdata:iedata,jsdata:jedata,siz(3),field(num_fields)%nbuf))
             endif
             slope=1.0;intercept=0.0
             if (units /= 'same') call convert_units(trim(field(num_fields)%units),trim(units),slope,intercept)
             if (verb.and.units /= 'same') then
                 write(stdout(),*) 'attempting to convert data to units = ',trim(units)
                 write(stdout(),'(a,f8.3,a,f8.3)') 'factor = ',slope,' offset= ',intercept
             endif
             field(num_fields)%slope = slope
             field(num_fields)%intercept = intercept
             allocate(field(num_fields)%ibuf(field(num_fields)%nbuf))
             select case (field(num_fields)%ndim)
             case (4)
                allocate(field(num_fields)%buf3d(isdata:iedata,jsdata:jedata,siz_in(3)))
                do j=1,field(num_fields)%nbuf
                   if (PRESENT(domain)) then
                      call mpp_read(unit,flds(i),domain,field(num_fields)%buf3d,j)
                   else
                      call mpp_read(unit,flds(i),field(num_fields)%buf3d,j)
                   endif
                   field(num_fields)%data(iscomp:iecomp,jscomp:jecomp,:,j) = &
                        field(num_fields)%buf3d(iscomp:iecomp,jscomp:jecomp,:)*slope + intercept 
                   field(num_fields)%ibuf(j) = j
                end do
             case (3)
                allocate(field(num_fields)%buf2d(isdata:iedata,jsdata:jedata))
                do j=1,field(num_fields)%nbuf
                   if (PRESENT(domain)) then
                      call mpp_read(unit,flds(i),domain,field(num_fields)%buf2d,j)
                   else
                      call mpp_read(unit,flds(i),field(num_fields)%buf2d,j)
                   endif
                   field(num_fields)%data(iscomp:iecomp,jscomp:jecomp,1,j) = &
                        field(num_fields)%buf2d(iscomp:iecomp,jscomp:jecomp)*slope + intercept
                   field(num_fields)%ibuf(j) = j
                enddo
             end select
             allocate(field(num_fields)%time(ntime))
             call mpp_get_atts(time_axis,units=units)
             do j=1,ntime
                call get_cal_time(tstamp(j),trim(units),yr,mon,day,hr,minu,sec)
                field(num_fields)%time(j) = set_date(yr,mon,day,hr,minu,sec)
             enddo
             do j=1,ntime-1
                if (field(num_fields)%time(j) >= field(num_fields)%time(j+1)) &
                     call mpp_error(FATAL,'times not monotonically increasing')
             enddo
             
             ! try to determine whether the data are hourly, daily, monthly or annual mean values
             ! were not taking time values to represent the exact midpoint of the
             ! time intervals since this usually isnt what the times are intended to
             ! represent.  Only checking the first 2 records.  This is not a very good test but
             ! should be sufficient for most cases
             call get_date(field(num_fields)%time(1),yr,mon,day,hr,minu,sec)
             field(num_fields)%tavg_type = FIELD_UNDEF
             if (ntime .gt. 1) then
                call get_date(field(num_fields)%time(2),yr2,mon2,day2,hr2,min2,sec2)
                if (yr2-yr == 1 .and. mon /= 12) field(num_fields)%tavg_type = FIELD_ANNUAL
                if (yr2-yr == 1 .and. mon == 12 .and. mon2 == 1) field(num_fields)%tavg_type = FIELD_MONTHLY
                if (yr2-yr == 0 .and. mon2-mon == 1) field(num_fields)%tavg_type = FIELD_MONTHLY
                if (yr2-yr == 0 .and. mon2-mon == 0 .and. day2-day == 1) field(num_fields)%tavg_type = FIELD_DAILY
                if (yr2-yr == 0 .and. mon2-mon == 0 .and. day2-day == 0 .and. hr2-hr == 1) &
                     field(num_fields)%tavg_type = FIELD_HOURLY
                ! a modest attempt at more robust time frequency detection. In the case where the 2nd record crosses a
                ! month boundary, check the third record. If it doesn't cross a month boundary then the field
                ! is probably daily resolution.
                if (field(num_fields)%tavg_type .eq. FIELD_MONTHLY) then
                   if (ntime .gt. 2) then
                      call get_date(field(num_fields)%time(3),yr,mon,day,hr,minu,sec)
                      if (yr-yr2 == 0 .and. mon-mon2 == 0 .and. day-day2 == 1) field(num_fields)%tavg_type = FIELD_DAILY
                   endif
                endif
             else
                field(num_fields)%tavg_type = FIELD_ANNUAL
             endif
             if (field(num_fields)%tavg_type == FIELD_UNDEF) call mpp_error(FATAL,'error determining time averaging interval')
             field(num_fields)%modulo_time = get_axis_modulo(time_axis)
             allocate(field(num_fields)%period(ntime))
             allocate(field(num_fields)%start_time(ntime))
             allocate(field(num_fields)%end_time(ntime))
             ! determine time boundaries based on assumed averaging interval.
             ! we are not using this information now, but storing it for future use
             select case (field(num_fields)%tavg_type)
             case (FIELD_ANNUAL)
                call set_annual_time(field(num_fields)%time,field(num_fields)%start_time,field(num_fields)%end_time,&
                     field(num_fields)%period,field(num_fields)%modulo_time,verbose=verb)
             case (FIELD_MONTHLY)
                call set_monthly_time(field(num_fields)%time,field(num_fields)%start_time,field(num_fields)%end_time,&
                  field(num_fields)%period,field(num_fields)%modulo_time,verbose=verb)
             case (FIELD_DAILY)
                call set_daily_time(field(num_fields)%time,field(num_fields)%start_time,field(num_fields)%end_time,&
                     field(num_fields)%period,field(num_fields)%modulo_time,verbose=verb)
             case (FIELD_HOURLY)
                call set_hourly_time(field(num_fields)%time,field(num_fields)%start_time,field(num_fields)%end_time,&
                     field(num_fields)%period,field(num_fields)%modulo_time,verbose=verb)
             case default
                call mpp_error(FATAL,'not able to parse time ')    
             end select
             if (verb) then
                if (field(num_fields)%modulo_time) write(stdout(),*) 'data are being treated as modulo in time'
                select case (field(num_fields)%tavg_type)
                case (FIELD_ANNUAL)
                   write(stdout(),*) 'data are being interpreted as annual means'
                case (FIELD_MONTHLY)
                   write(stdout(),*) 'data are being interpreted as monthly means'
                case (FIELD_DAILY)
                   write(stdout(),*) 'data are being interpreted as daily means'
                case (FIELD_HOURLY)
                   write(stdout(),*) 'data are being interpreted as hourly means'
                end select
             endif
          end if
       enddo
      
      if (field(num_fields)%name == 'none') call mpp_error(FATAL,'external field field not found')



      deallocate(global_atts)
      deallocate(axes)
      deallocate(flds)

      return
      
    end function init_external_field
    
!</FUNCTION> NAME="init_external_field"


    subroutine time_interp_external_2d(index, time, data_in, interp, verbose)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real, dimension(:,:), intent(inout) :: data_in
      integer, intent(in), optional :: interp
      logical, intent(in), optional :: verbose
      
      integer :: t1, t2
      real, dimension(size(data_in,1), size(data_in,2), 1) :: data_out
      integer :: interp_out = LINEAR_TIME_INTERP
      logical :: verbose_out = .false.

      if (PRESENT(interp)) interp_out = interp
      if (PRESENT(verbose)) verbose_out = verbose

      call time_interp_external_3d(index, time, data_out, interp_out, verbose_out)

      data_in(:,:) = data_out(:,:,1)

      return
    end subroutine time_interp_external_2d

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

    subroutine time_interp_external_3d(index, time, data, interp,verbose)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real, dimension(:,:,:), intent(inout) :: data
      integer, intent(in), optional :: interp
      logical, intent(in), optional :: verbose
      
      integer :: nx, ny, nz, nt, interp_method, t1, t2, i
      integer :: yr1, mon1, day1, hr1, min1, sec1
      integer :: yr2, mon2, day2, hr2, min2, sec2
      integer :: i1, i2, isc, iec, jsc, jec, mod_time
      integer :: yy, mm, dd, hh, min, ss
      real :: w1,w2
      logical :: verb
      type(time_type) :: time_out

      nx = size(data,1)
      ny = size(data,2)
      nz = size(data,3)

      interp_method = LINEAR_TIME_INTERP
      if (PRESENT(interp)) interp_method = interp
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      
      if (nx.ne.field(index)%siz(1).or.ny.ne.field(index)%siz(2).or.nz.ne.field(index)%siz(3)) &
           call mpp_error(FATAL,'array size mismatch in time_interp_external')

      if (index > num_fields) call mpp_error(FATAL,'index exceeds available datasets')

      
      isc=field(index)%isc;iec=field(index)%iec
      jsc=field(index)%jsc;jec=field(index)%jec

      
      if (field(index)%tavg_type == FIELD_ANNUAL .and. field(index)%siz(4) == 1) then
         data = field(index)%data(:,:,:,1)
      else
!         call get_time_weights(time, mid_times,t1,t2,w1,w2,field(index)%modulo_time,verbose=verb)
         mod_time=0;if(field(index)%modulo_time) mod_time=1
         call time_interp(time,field(index)%time(:),w2,t1,t2,modtime=mod_time)
         w1 = 1.0-w2
         if (verb) then
            call get_date(time,yy,mm,dd,hh,min,ss)
            write(stdout(),'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',min,':',ss
            write(stdout(),*) 't1, t2, w1, w2= ', t1, t2, w1, w2
         endif
         i1 = find_buf_index(t1,field(index)%ibuf)
         i2 = find_buf_index(t2,field(index)%ibuf)

         if (verb) then
            write(stdout(),*) 'ibuf= ',field(index)%ibuf
            write(stdout(),*) 'i1,i2= ',i1, i2
         endif

         do while (i1 == -1 .or. i2 == -1)
            if (i1 == -1) then
               ! shift data array and index array
               ! this may be expensive for large arrays since it involves array copies
               field(index)%data = eoshift(field(index)%data,1,dim=4)
               field(index)%ibuf = eoshift(field(index)%ibuf,1)
               if (field(index)%domain_present) then
                  select case (field(index)%ndim)
                  case (3)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%domain,field(index)%buf2d,t1)
                     field(index)%data(isc:iec,jsc:jec,1,field(index)%nbuf) = &
                          field(index)%buf2d(isc:iec,jsc:jec)*field(index)%slope + field(index)%intercept
                  case (4)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%domain,field(index)%buf3d,t1)
                     field(index)%data(isc:iec,jsc:jec,:,field(index)%nbuf) = &
                          field(index)%buf3d(isc:iec,jsc:jec,:)*field(index)%slope + field(index)%intercept
                  end select
               else
                  select case (field(index)%ndim)
                  case (3)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%buf2d,t1)
                     field(index)%data(isc:iec,jsc:jec,1,field(index)%nbuf) = &
                          field(index)%buf2d(isc:iec,jsc:jec)*field(index)%slope + field(index)%intercept
                  case (4)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%buf3d,t1)
                     field(index)%data(isc:iec,jsc:jec,:,field(index)%nbuf) = &
                          field(index)%buf3d(isc:iec,jsc:jec,:)*field(index)%slope + field(index)%intercept
                  end select
               endif
               field(index)%ibuf(field(index)%nbuf) = t1
               i1 = find_buf_index(t1,field(index)%ibuf)
               i2 = find_buf_index(t2,field(index)%ibuf)           
               if (verb) then
                  write(stdout(),*) 'ibuf= ',field(index)%ibuf
                  write(stdout(),*) 'i1,i2= ',i1, i2
               endif
            endif
            if (i2 == -1) then
               ! shift data array and index array
               ! this may be expensive for large arrays since it involves array copies
               field(index)%data = eoshift(field(index)%data,1,dim=4)
               field(index)%ibuf = eoshift(field(index)%ibuf,1)
               if (field(index)%domain_present) then
                  select case(field(index)%ndim)
                  case (3)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%domain,field(index)%buf2d,t2)
                     field(index)%data(isc:iec,jsc:jec,1,field(index)%nbuf) = &
                          field(index)%buf2d(isc:iec,jsc:jec)*field(index)%slope + field(index)%intercept
                  case (4)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%domain,field(index)%buf3d,t2)
                     field(index)%data(isc:iec,jsc:jec,:,field(index)%nbuf) = &
                          field(index)%buf3d(isc:iec,jsc:jec,:)*field(index)%slope + field(index)%intercept
                  case default
                     call mpp_error(FATAL,'invalid array rank')
                  end select
               else
                  select case(field(index)%ndim)
                  case (3)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%buf2d,t2)
                     field(index)%data(isc:iec,jsc:jec,1,field(index)%nbuf) = &
                          field(index)%buf2d(isc:iec,jsc:jec)*field(index)%slope + field(index)%intercept
                  case (4)
                     call mpp_read(field(index)%unit,field(index)%field,field(index)%buf3d,t2)
                     field(index)%data(isc:iec,jsc:jec,:,field(index)%nbuf) = &
                          field(index)%buf3d(isc:iec,jsc:jec,:)*field(index)%slope + field(index)%intercept
                  case default
                     call mpp_error(FATAL,'invalid array rank')
                  end select
               endif
               field(index)%ibuf(field(index)%nbuf) = t2
               i1 = find_buf_index(t1,field(index)%ibuf)
               i2 = find_buf_index(t2,field(index)%ibuf)           
               if (verb) then
                  write(stdout(),*) 'ibuf= ',field(index)%ibuf
                  write(stdout(),*) 'i1,i2= ',i1, i2
               endif
            endif
         end do
         data = field(index)%data(:,:,:,i1)*w1 + (field(index)%data(:,:,:,i2))*w2
      endif

      return
      
99    format (a,i4,'/',i2,'/',i2,'/',1x,i2,':',i2,':',i2)

    end subroutine time_interp_external_3d
!</SUBROUTINE> NAME="time_interp_external"
    
    subroutine set_monthly_time(time,start_time,end_time,period,modulo,verbose)
! USER BEWARE :::  incoming time needs to lie within the intended month
! i.e. 2/1/1900 will be interpreted as February data
      type(time_type), intent(inout), dimension(:) :: time, start_time, end_time
      real, intent(out), dimension(:) :: period
      logical, intent(in), optional :: modulo,verbose

      integer :: yr,mon,day,hr,minu,sec,monp1,ndays,ntime,i
      logical :: mod_time, leap_yr,verb

      ntime = size(time)
      mod_time = .false.
      if (PRESENT(modulo)) mod_time = modulo
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      
      do i=1,ntime
         call get_date(time(i),yr,mon,day,hr,minu,sec)
         if (verb .and. mpp_pe() == 0) write(stdout(),99) 'time from dataset yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         if (mod_time) yr = modulo_year ! set year to climatological time if modulo axis
         if (verb .and. mod_time)write(stdout(),*) 'adjusting year to climatological year'
         time(i) = set_date(yr,mon,day,hr,0,0)         
! position time at midpoint of month
         ndays = days_in_month(time(i))
         period(i) = ndays
         day = ndays/2+1
         if(mod(ndays,2) /= 0) then
             hr = 12
         else
             hr = 0
         endif
         time(i) = set_date(yr,mon,day,hr,0,0)         
         start_time(i) = set_date(yr,mon,1,0,0,0)
         monp1 = mod(mon,12) + 1
         end_time(i) = set_date(yr,monp1,1,0,0,0)
         if (verb) then
             call get_date(time(i),yr,mon,day,hr,minu,sec)
             if (verb .and. mpp_pe() == 0) write(*,99) 'adjusted time yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         endif         
      enddo
99    format (a,i4,'/',i2,'/',i2,'/',1x,i2,':',i2,':',i2)
      
      return
    end subroutine set_monthly_time

    subroutine set_annual_time(time,start_time,end_time,period,modulo,verbose)

      type(time_type), intent(inout), dimension(:) :: time, start_time, end_time
      real, intent(out), dimension(:) :: period
      logical, intent(in), optional :: modulo, verbose

      integer :: yr,mon,day,hr,minu,sec,monp1,ndays,ntime,i
      type(time_type) :: year_mid, year_beg, year_end
      logical :: mod_time, leap_yr, verb

      ntime = size(time)
      mod_time = .false.
      if (PRESENT(modulo)) mod_time = modulo
      verb = .false.
      if (PRESENT(verbose)) verb = verbose
      
      do i=1,ntime
         call get_date(time(i),yr,mon,day,hr,minu,sec)
         if (verb .and. mpp_pe() == 0) write(stdout(),99) 'time from dataset yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         if (mod_time) yr = modulo_year ! set year to climatological time if modulo axis
         if (verb .and. mod_time .and. mpp_pe() == 0) write(stdout(),*) 'adjusting year to climatological year'
         year_beg =  set_date(yr,1,1)
         year_end =  set_date(yr+1,1,1)
         time(i) = (year_beg + year_end) / 2
         if (verb .and. mpp_pe() == 0) then
             call get_date(time(i),yr,mon,day,hr,minu,sec)
             write(*,99) 'adjusted time yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         endif
         start_time(i) = year_beg
         end_time(i) = year_end
         if (mod_time) then
             period(i) = 365.25
         else
             period(i) = days_in_year(time(i))
         endif
      enddo
           
99    format (a,i4,'/',i2,'/',i2,'/',1x,i2,':',i2,':',i2) 
      return
    end subroutine set_annual_time

    subroutine set_daily_time(time,start_time,end_time,period,modulo,verbose)
! USER BEWARE :: day is assumed to lie within (not inclusive ) intended time period 
      type(time_type), intent(inout), dimension(:) :: time, start_time, end_time
      real, intent(out), dimension(:) :: period
      logical, intent(in), optional :: modulo,verbose

      integer :: yr,mon,day,hr,minu,sec,monp1,ndays,ntime,i
      logical :: mod_time, leap_yr,verb
      character(len=32) :: msg

      ntime = size(time)
      mod_time = .false.
      if (PRESENT(modulo)) mod_time = modulo
      verb=.false.
      if (PRESENT(verbose)) verb = verbose
      
      do i=1,ntime
         call get_date(time(i),yr,mon,day,hr,minu,sec)
         if (verb .and. mpp_pe() ==0) write(stdout(),99) 'time from dataset yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         if (mod_time) yr = modulo_year ! set year to climatological time if modulo axis
         if (verb.and.mod_time .and. mpp_pe() == 0) write(stdout(),*) 'adjusting year to climatological year'

! position time at midpoint of day
         hr = 12
         period(i) = 1
         time(i) = set_date(yr,mon,day,hr,0,0)
         start_time(i) = set_date(yr,mon,day,0,0,0)
         end_time(i) = increment_time(start_time(i),0,1)
         if (verb .and. mpp_pe() ==0) then
             call get_date(time(i),yr,mon,day,hr,minu,sec)
             write(*,99) 'adjusted time yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         endif         
      enddo
99    format (a,i4,'/',i2,'/',i2,'/',1x,i2,':',i2,':',i2)       
      return
    end subroutine set_daily_time

     subroutine set_hourly_time(time,start_time,end_time,period,modulo,verbose)
! USER BEWARE :: hour is assumed to lie within (not inclusive ) intended time period 
      type(time_type), intent(inout), dimension(:) :: time, start_time, end_time
      real, intent(out), dimension(:) :: period
      logical, intent(in), optional :: modulo,verbose

      integer :: yr,mon,day,hr,minu,sec,monp1,ndays,ntime,i
      logical :: mod_time, leap_yr,verb

      ntime = size(time)
      mod_time = .false.
      if (PRESENT(modulo)) mod_time = modulo
      verb=.false.
      if (PRESENT(verbose)) verb = verbose

      do i=1,ntime
         call get_date(time(i),yr,mon,day,hr,minu,sec)
         if (verb .and. mpp_pe() == 0) write(stdout(),99) 'time from dataset yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         if (mod_time) yr = modulo_year ! set year to climatological time if modulo axis
         if (verb.and.mod_time.and.mpp_pe()==0) write(stdout(),*) 'adjusting year to climatological year'
! position time on 1/2 hour
         minu = 30
         period(i) = 1.0/24.0
         time(i) = set_date(yr,mon,day,hr,minu,0)
         start_time(i) = set_date(yr,mon,day,hr,0,0)
         end_time(i) = increment_time(start_time(i),3600,0)
         if (verb.and.mpp_pe()==0) then
             call get_date(time(i),yr,mon,day,hr,minu,sec)
             write(*,99) 'adjusted time yyyy/mm/dd hh:mm:ss= ',yr,mon,day,hr,minu,sec
         endif      
      enddo
99    format (a,i4,'/',i2,'/',i2,'/',1x,i2,':',i2,':',i2)       
      return
    end subroutine set_hourly_time


    function find_buf_index(indx,buf)
      integer :: indx
      integer, dimension(:) :: buf
      integer :: find_buf_index

      integer :: nbuf, i
      
      nbuf = size(buf)

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
      
      if (index .lt. 1 .or. index .gt. num_fields) call mpp_error(FATAL,'invalid index in call to get_external_field_size')


      get_external_field_size(1) = field(index)%siz(1)
      get_external_field_size(2) = field(index)%siz(2)
      get_external_field_size(3) = field(index)%siz(3)
      get_external_field_size(4) = field(index)%siz(4)

    end function get_external_field_size
!</FUNCTION> NAME="get_external_field"

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
              field(i)%period,field(i)%data,field(i)%ibuf)
         if (ASSOCIATED(field(i)%buf2d)) deallocate(field(i)%buf2d)
         if (ASSOCIATED(field(i)%buf3d)) deallocate(field(i)%buf3d)
         do j=1,4
            field(i)%axes(j) = default_axis
         enddo
         field(i)%domain = NULL_DOMAIN2D
         field(i)%field = default_field
         field(i)%nbuf = 0
         field(i)%slope = 0.
         field(i)%intercept = 0.
      enddo
      
      num_fields = 0

      time_interp_ext_initialized = .false.

    end subroutine time_interp_external_exit
!</SUBROUTINE> NAME="time_interp_external_exit"

end module time_interp_external_mod

#ifdef test_time_interp_external

program test_time_interp_ext

use mpp_mod, only : mpp_init, mpp_exit, mpp_npes, stdout, stdlog, FATAL, mpp_error
use mpp_io_mod, only : mpp_io_init, mpp_io_exit, mpp_open, MPP_RDONLY, MPP_ASCII, mpp_close
use mpp_domains_mod, only : mpp_domains_init, domain2d, mpp_define_layout, mpp_define_domains,&
     mpp_global_sum, mpp_global_max, mpp_global_min, BITWISE_EXACT_SUM, mpp_get_compute_domain, &
     mpp_domains_set_stack_size
use time_interp_external_mod, only : time_interp_external, time_interp_external_init,&
     time_interp_external_exit, time_interp_external, init_external_field, get_External_field_size
use time_manager_mod, only : get_date, set_date, time_manager_init, set_calendar_type, JULIAN, time_type, increment_time

implicit none



integer :: id, i, io_status, unit
character(len=128) :: filename, fieldname
type(time_type) :: time
real, allocatable, dimension(:,:,:) :: data_d, data_g
type(domain2d) :: domain
integer :: layout(2), fld_size(4)
integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
integer :: yy, mm, dd, hh, mm, ss
real :: sm,mx,mn
integer :: ntime=12,year0=1991,month0=1,day0=1,days_inc=31

namelist /test_time_interp_ext_nml/ filename, fieldname,ntime,year0,month0,&
     day0,days_inc


call mpp_init
call mpp_io_init
call mpp_domains_init
call time_interp_external_init
call time_manager_init
call set_calendar_type(JULIAN)


call mpp_open(unit,'input.nml',action=MPP_RDONLY,form=MPP_ASCII)
read(unit,test_time_interp_ext_nml,iostat=io_status)
write(stdlog(),test_time_interp_ext_nml)
if (io_status .gt. 0) then
   call mpp_error(FATAL,'=>Error reading test_time_interp_ext_nml')
endif
call mpp_close(unit)

write(stdout(),*) 'INTERPOLATING NON DECOMPOSED FIELDS'
write(stdout(),*) '======================================'

call time_interp_external_init

id = init_external_field(filename,fieldname,verbose=.true.)

fld_size = get_external_field_size(id)

allocate(data_g(fld_size(1),fld_size(2),fld_size(3)))

time = set_date(year0,month0,day0,0,0,0)

do i=1,ntime
   call time_interp_external(id,time,data_g,verbose=.true.)
   sm = sum(data_g)
   mn = minval(data_g)
   mx = maxval(data_g)
   write(stdout(),*) 'sum= ', sm
   write(stdout(),*) 'max= ', mx
   write(stdout(),*) 'min= ', mn
   time = increment_time(time,0,days_inc)
enddo

call mpp_define_layout((/1,fld_size(1),1,fld_size(2)/),mpp_npes(),layout)
call mpp_define_domains((/1,fld_size(1),1,fld_size(2)/),layout,domain)
call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
call mpp_get_compute_domain(domain,isd,ied,jsd,jed)

call mpp_domains_set_stack_size(fld_size(1)*fld_size(2)*min(fld_size(3),1)*2)
allocate(data_d(isd:ied,jsd:jed,fld_size(3)))


write(stdout(),*) 'INTERPOLATING DOMAIN DECOMPOSED FIELDS'
write(stdout(),*) '======================================'

id = init_external_field(filename,fieldname,domain=domain,verbose=.true.)

time = set_date(year0,month0,day0)

do i=1,ntime
   call time_interp_external(id,time,data_d,verbose=.true.)
   sm = mpp_global_sum(domain,data_d,flags=BITWISE_EXACT_SUM)
   mx = mpp_global_max(domain,data_d)
   mn = mpp_global_min(domain,data_d)
   write(stdout(),*) 'global sum= ', sm
   write(stdout(),*) 'global max= ', mx
   write(stdout(),*) 'global min= ', mn
   time = increment_time(time,0,days_inc)
enddo

call time_interp_external_exit


call mpp_io_exit
call mpp_exit
stop

end program test_time_interp_ext
#endif
    
  

      







