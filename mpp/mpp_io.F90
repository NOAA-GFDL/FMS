!-----------------------------------------------------------------------
!                 Parallel I/O for message-passing codes
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_io.html"/>

! <OVERVIEW>
!   <TT>mpp_io_mod</TT>, is a set of simple calls for parallel I/O on
!   distributed systems. It is geared toward the writing of data in netCDF
!   format. It requires the modules <LINK
!   SRC="mpp_domains.html">mpp_domains_mod</LINK> and <LINK
!   SRC="mpp.html">mpp_mod</LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   In massively parallel environments, an often difficult problem is
!   the reading and writing of data to files on disk. MPI-IO and MPI-2 IO
!   are moving toward providing this capability, but are currently not
!   widely implemented. Further, it is a rather abstruse
!   API. <TT>mpp_io_mod</TT> is an attempt at a simple API encompassing a
!   certain variety of the I/O tasks that will be required. It does not
!   attempt to be an all-encompassing standard such as MPI, however, it
!   can be implemented in MPI if so desired. It is equally simple to add
!   parallel I/O capability to <TT>mpp_io_mod</TT> based on vendor-specific
!   APIs while providing a layer of insulation for user codes.
!   
!   The <TT>mpp_io_mod</TT> parallel I/O API built on top of the <LINK
!   SRC="mpp_domains.html">mpp_domains_mod</LINK> and <LINK
!   SRC="mpp.html">mpp_mod</LINK> API for domain decomposition and
!   message passing. Features of <TT>mpp_io_mod</TT> include:
!   
!    1) Simple, minimal API, with free access to underlying API for more
!   complicated stuff.<BR/>
!    2) Self-describing files: comprehensive header information
!   (metadata) in the file itself.<BR/>
!    3) Strong focus on performance of parallel write: the climate models
!   for which it is designed typically read a minimal amount of data
!   (typically at the beginning of the run); but on the other hand, tend
!   to write copious amounts of data during the run. An interface for
!   reading is also supplied, but its performance has not yet been optimized.<BR/>
!    4) Integrated netCDF capability: <LINK SRC
!   ="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</LINK> is a
!   data format widely used in the climate/weather modeling
!   community. netCDF is considered the principal medium of data storage
!   for <TT>mpp_io_mod</TT>. But I provide a raw unformatted
!   fortran I/O capability in case netCDF is not an option, either due to
!   unavailability, inappropriateness, or poor performance.<BR/>
!    5) May require off-line post-processing: a tool for this purpose,
!   <TT>mppnccombine</TT>, is available. GFDL users may use
!   <TT>~hnv/pub/mppnccombine</TT>. Outside users may obtain the
!   source <LINK SRC
!   ="ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</LINK>.  It
!   can be compiled on any C compiler and linked with the netCDF
!   library. The program is free and is covered by the <LINK SRC
!   ="ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</LINK>.
!   
!   The internal representation of the data being written out is
!   assumed be the default real type, which can be 4 or 8-byte. Time data
!   is always written as 8-bytes to avoid overflow on climatic time scales
!   in units of seconds.
!   
!   <LINK SRC="modes"></LINK><H4>I/O modes in <TT>mpp_io_mod</TT></H4>
!   
!   The I/O activity critical to performance in the models for which
!   <TT>mpp_io_mod</TT> is designed is typically the writing of large
!   datasets on a model grid volume produced at intervals during
!   a run. Consider a 3D grid volume, where model arrays are stored as
!   <TT>(i,j,k)</TT>. The domain decomposition is typically along
!   <TT>i</TT> or <TT>j</TT>: thus to store data to disk as a global
!   volume, the distributed chunks of data have to be seen as
!   non-contiguous. If we attempt to have all PEs write this data into a
!   single file, performance can be seriously compromised because of the
!   data reordering that will be required. Possible options are to have
!   one PE acquire all the data and write it out, or to have all the PEs
!   write independent files, which are recombined offline. These three
!   modes of operation are described in the <TT>mpp_io_mod</TT> terminology
!   in terms of two parameters, <I>threading</I> and <I>fileset</I>,
!   as follows:
!   
!   <I>Single-threaded I/O:</I> a single PE acquires all the data
!   and writes it out.<BR/>
!   <I>Multi-threaded, single-fileset I/O:</I> many PEs write to a
!   single file.<BR/>
!    <I>Multi-threaded, multi-fileset I/O:</I> many PEs write to
!   independent files. This is also called <I>distributed I/O</I>.
!   
!   The middle option is the most difficult to achieve performance. The
!   choice of one of these modes is made when a file is opened for I/O, in
!   <LINK SRC="#mpp_open">mpp_open</LINK>.
!   
!   <LINK name="metadata"></LINK><H4>Metadata in <TT>mpp_io_mod</TT></H4>
!   
!   A requirement of the design of <TT>mpp_io_mod</TT> is that the file must
!   be entirely self-describing: comprehensive header information
!   describing its contents is present in the header of every file. The
!   header information follows the model of netCDF. Variables in the file
!   are divided into <I>axes</I> and <I>fields</I>. An axis describes a
!   co-ordinate variable, e.g <TT>x,y,z,t</TT>. A field consists of data in
!   the space described by the axes. An axis is described in
!   <TT>mpp_io_mod</TT> using the defined type <TT>axistype</TT>:
!   
!   <PRE>
!   type, public :: axistype
!      sequence
!      character(len=128) :: name
!      character(len=128) :: units
!      character(len=256) :: longname
!      character(len=8) :: cartesian
!      integer :: len
!      integer :: sense           !+/-1, depth or height?
!      type(domain1D), pointer :: domain
!      real, dimension(:), pointer :: data
!      integer :: id, did
!      integer :: type  ! external NetCDF type format for axis data
!      integer :: natt
!      type(atttype), pointer :: Att(:) ! axis attributes
!   end type axistype
!   </PRE>
!   
!   A field is described using the type <TT>fieldtype</TT>:
!   
!   <PRE>
!   type, public :: fieldtype
!      sequence
!      character(len=128) :: name
!      character(len=128) :: units
!      character(len=256) :: longname
!      real :: min, max, missing, fill, scale, add
!      integer :: pack
!      type(axistype), dimension(:), pointer :: axes
!      integer, dimension(:), pointer :: size
!      integer :: time_axis_index
!      integer :: id
!      integer :: type ! external NetCDF format for field data
!      integer :: natt, ndim
!      type(atttype), pointer :: Att(:) ! field metadata
!   end type fieldtype
!   </PRE>
!   
!   An attribute (global, field or axis) is described using the <TT>atttype</TT>:
!   
!   <PRE>
!   type, public :: atttype
!      sequence
!      integer :: type, len
!      character(len=128) :: name
!      character(len=256)  :: catt
!      real(FLOAT_KIND), pointer :: fatt(:)
!   end type atttype
!   </PRE>
!   
!   <LINK name="packing"></LINK>This default set of field attributes corresponds
!   closely to various conventions established for netCDF files. The
!   <TT>pack</TT> attribute of a field defines whether or not a
!   field is to be packed on output. Allowed values of
!   <TT>pack</TT> are 1,2,4 and 8. The value of
!   <TT>pack</TT> is the number of variables written into 8
!   bytes. In typical use, we write 4-byte reals to netCDF output; thus
!   the default value of <TT>pack</TT> is 2. For
!   <TT>pack</TT> = 4 or 8, packing uses a simple-minded linear
!   scaling scheme using the <TT>scale</TT> and <TT>add</TT>
!   attributes. There is thus likely to be a significant loss of dynamic
!   range with packing. When a field is declared to be packed, the
!   <TT>missing</TT> and <TT>fill</TT> attributes, if
!   supplied, are packed also.
!   
!   Please note that the pack values are the same even if the default
!   real is 4 bytes, i.e <TT>PACK=1</TT> still follows the definition
!   above and writes out 8 bytes.
!   
!   A set of <I>attributes</I> for each variable is also available. The
!   variable definitions and attribute information is written/read by calling
!   <LINK SRC="#mpp_write_meta">mpp_write_meta</LINK> or <LINK SRC="#mpp_read_meta">mpp_read_meta</LINK>. A typical calling
!   sequence for writing data might be:
!   
!   <PRE>
!   ...
!     type(domain2D), dimension(:), allocatable, target :: domain
!     type(fieldtype) :: field
!     type(axistype) :: x, y, z, t
!   ...
!     call mpp_define_domains( (/1,nx,1,ny/), domain )
!     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
!                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
!   ...
!     call mpp_write_meta( unit, x, 'X', 'km', 'X distance', &
!          domain=domain(pe)%x, data=(/(float(i),i=1,nx)/) )
!     call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', &
!          domain=domain(pe)%y, data=(/(float(i),i=1,ny)/) )
!     call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', &
!          data=(/(float(i),i=1,nz)/) )
!     call mpp_write_meta( unit, t, 'Time', 'second', 'Time' )
!   
!     call mpp_write_meta( unit, field, (/x,y,z,t/), 'a', '(m/s)', AAA', &
!          missing=-1e36 )
!   ...
!     call mpp_write( unit, x )
!     call mpp_write( unit, y )
!     call mpp_write( unit, z )
!   ...
!   </PRE>
!   
!   In this example, <TT>x</TT> and <TT>y</TT> have been
!   declared as distributed axes, since a domain decomposition has been
!   associated. <TT>z</TT> and <TT>t</TT> are undistributed
!   axes. <TT>t</TT> is known to be a <I>record</I> axis (netCDF
!   terminology) since we do not allocate the <TT>data</TT> element
!   of the <TT>axistype</TT>. <I>Only one record axis may be
!   associated with a file.</I> The call to <LINK
!   SRC="#mpp_write_meta">mpp_write_meta</LINK> initializes
!   the axes, and associates a unique variable ID with each axis. The call
!   to <TT>mpp_write_meta</TT> with argument <TT>field</TT>
!   declared <TT>field</TT> to be a 4D variable that is a function
!   of <TT>(x,y,z,t)</TT>, and a unique variable ID is associated
!   with it. A 3D field will be written at each call to
!   <TT>mpp_write(field)</TT>.
!   
!   The data to any variable, including axes, is written by
!   <TT>mpp_write</TT>.
!   
!   Any additional attributes of variables can be added through
!   subsequent <TT>mpp_write_meta</TT> calls, using the variable ID as a
!   handle. <I>Global</I> attributes, associated with the dataset as a
!   whole, can also be written thus. See the <LINK
!   SRC="#mpp_write_meta">mpp_write_meta</LINK> call syntax below
!   for further details.
!   
!   You cannot interleave calls to <TT>mpp_write</TT> and
!   <TT>mpp_write_meta</TT>: the first call to
!   <TT>mpp_write</TT> implies that metadata specification is
!   complete.
!   
!   A typical calling sequence for reading data might be:
!   
!   <PRE>
!   ...
!     integer :: unit, natt, nvar, ntime
!     type(domain2D), dimension(:), allocatable, target :: domain
!     type(fieldtype), allocatable, dimension(:) :: fields
!     type(atttype), allocatable, dimension(:) :: global_atts
!     real, allocatable, dimension(:) :: times
!   ...
!     call mpp_define_domains( (/1,nx,1,ny/), domain )
!   
!     call mpp_read_meta(unit)
!     call mpp_get_info(unit,natt,nvar,ntime)
!     allocate(global_atts(natt))
!     call mpp_get_atts(unit,global_atts)
!     allocate(fields(nvar))
!     call mpp_get_vars(unit, fields)
!     allocate(times(ntime))
!     call mpp_get_times(unit, times)
!   
!     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
!                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
!   ...
!     do i=1, nvar
!       if (fields(i)%name == 'a')  call mpp_read(unit,fields(i),domain(pe), a,
!                                                 tindex)
!     enddo
!   ...
!   </PRE>
!   
!   In this example, the data are distributed as in the previous
!   example. The call to <LINK
!   SRC="#mpp_read_meta">mpp_read_meta</LINK> initializes
!   all of the metadata associated with the file, including global
!   attributes, variable attributes and non-record dimension data. The
!   call to <TT>mpp_get_info</TT> returns the number of global
!   attributes (<TT>natt</TT>), variables (<TT>nvar</TT>) and
!   time levels (<TT>ntime</TT>) associated with the file
!   identified by a unique ID (<TT>unit</TT>).
!   <TT>mpp_get_atts</TT> returns all global attributes for
!   the file in the derived type <TT>atttype(natt)</TT>.
!   <TT>mpp_get_vars</TT> returns variable types
!   (<TT>fieldtype(nvar)</TT>).  Since the record dimension data are not allocated for calls to <LINK SRC="#mpp_write">mpp_write</LINK>, a separate call to  <TT>mpp_get_times</TT> is required to access record dimension data.  Subsequent calls to
!   <TT>mpp_read</TT> return the field data arrays corresponding to
!   the fieldtype.  The <TT>domain</TT> type is an optional
!   argument.  If <TT>domain</TT> is omitted, the incoming field
!   array should be dimensioned for the global domain, otherwise, the
!   field data is assigned to the computational domain of a local array.
!   
!   <I>Multi-fileset</I> reads are not supported with <TT>mpp_read</TT>.

! </DESCRIPTION>

module mpp_io_mod

use mpp_data_mod,       only : default_field, default_axis, default_att
use mpp_datatype_mod,   only : axistype, atttype, fieldtype, validtype
use mpp_parameter_mod,  only : MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR, MPP_ASCII
use mpp_parameter_mod,  only : MPP_IEEE32, MPP_NATIVE, MPP_NETCDF, MPP_SEQUENTIAL
use mpp_parameter_mod,  only : MPP_DIRECT, MPP_SINGLE, MPP_MULTI, MPP_DELETE, MPP_COLLECT
use mpp_io_util_mod,    only : mpp_get_iospec, mpp_get_id, mpp_get_ncid, mpp_get_unit_range, mpp_is_valid
use mpp_io_util_mod,    only : mpp_set_unit_range, mpp_get_info, mpp_get_atts, mpp_get_fields
use mpp_io_util_mod,    only : mpp_get_times, mpp_get_axes, mpp_get_recdimid, mpp_get_axis_data
use mpp_io_util_mod,    only : mpp_io_set_stack_size, mpp_get_field_index, mpp_get_axis_index
use mpp_io_misc_mod,    only : mpp_io_init, mpp_io_exit, netcdf_err, mpp_flush
use mpp_io_write_mod,   only : mpp_write, mpp_write_meta, mpp_copy_meta, mpp_modify_meta
use mpp_io_read_mod,    only : mpp_read, mpp_read_meta, mpp_get_tavg_info
use mpp_io_connect_mod, only : mpp_open, mpp_close

implicit none
private

  character(len=128) :: version= &
       '$Id: mpp_io.F90,v 12.0 2005/04/14 17:58:20 fms Exp $'
  character(len=128) :: tagname= &
       '$Name: lima $'

  !--- public parameters  -----------------------------------------------
  public :: MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR, MPP_ASCII, MPP_IEEE32
  public :: MPP_NATIVE, MPP_NETCDF, MPP_SEQUENTIAL, MPP_DIRECT, MPP_SINGLE
  public :: MPP_MULTI, MPP_DELETE, MPP_COLLECT

  !--- public data type ------------------------------------------------
  public :: axistype, atttype, fieldtype, validtype

  !--- public data -----------------------------------------------------
  public :: default_field, default_axis, default_att
    
  !--- public interface from mpp_io_util_mod ----------------------
  public :: mpp_get_iospec, mpp_get_id, mpp_get_ncid, mpp_get_unit_range, mpp_is_valid
  public :: mpp_set_unit_range, mpp_get_info, mpp_get_atts, mpp_get_fields
  public :: mpp_get_times, mpp_get_axes, mpp_get_recdimid, mpp_get_axis_data
  public :: mpp_io_set_stack_size, mpp_get_field_index, mpp_get_axis_index

  !--- public interface from mpp_io_misc_mod ----------------------
  public :: mpp_io_init, mpp_io_exit, netcdf_err, mpp_flush

  !--- public interface from mpp_io_write_mod ---------------------
  public :: mpp_write, mpp_write_meta, mpp_copy_meta, mpp_modify_meta

  !--- public interface from mpp_io_read_mod ---------------------
  public :: mpp_read, mpp_read_meta, mpp_get_tavg_info

  !--- public interface from mpp_io_switch_mod ---------------------
  public :: mpp_open, mpp_close

end module mpp_io_mod

#ifdef test_mpp_io
program mpp_io_test
#include <fms_platform.h>

  use mpp_mod,         only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error
  use mpp_mod,         only : FATAL, mpp_chksum, MPP_DEBUG, mpp_set_stack_size
  use mpp_mod,         only : mpp_broadcast, mpp_sync, mpp_exit
  use mpp_domains_mod, only : mpp_define_domains, mpp_domains_set_stack_size, domain1D
  use mpp_domains_mod, only : domain2D, mpp_define_layout, mpp_get_domain_components
  use mpp_domains_mod, only : mpp_get_data_domain, mpp_get_compute_domain, mpp_domains_exit
  use mpp_io_mod,      only : mpp_io_init, mpp_write_meta, axistype, fieldtype, atttype
  use mpp_io_mod,      only : MPP_RDONLY, mpp_open, MPP_OVERWR, MPP_ASCII, MPP_SINGLE
  use mpp_io_mod,      only : MPP_NETCDF, MPP_MULTI, mpp_get_atts, mpp_write, mpp_close
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_get_fields, mpp_get_times
  use mpp_io_mod,      only : mpp_read, mpp_io_exit

  implicit none

#ifdef use_netCDF
#include <netcdf.inc>
#endif

  !--- namelist definition
  integer           :: nx=128, ny=128, nz=40, nt=2
  integer           :: halo=2, stackmax=1500000, stackmaxd=500000
  logical           :: debug=.FALSE.  
  character(len=64) :: file='test', iospec='-F cachea' 
  namelist / mpp_io_nml / nx, ny, nz, nt, halo, stackmax, stackmaxd, debug, file, iospec

  integer        :: pe, npes
  type(domain2D) :: domain

  integer            :: is, ie, js, je, isd, ied, jsd, jed
  integer            :: tk, tk0, tks_per_sec
  integer            :: i,j,k, unit=7, layout(2)
  logical            :: opened
  character(len=64)  :: varname
  integer            :: ndim, nvar, natt, ntime
  real(DOUBLE_KIND)  :: time
  type(axistype)     :: x, y, z, t
  type(fieldtype)    :: f
  type(domain1D)     :: xdom, ydom
  integer(LONG_KIND) :: rchk, chk
  type(atttype),          allocatable :: atts(:)
  type(fieldtype),        allocatable :: vars(:)
  type(axistype),         allocatable :: axes(:)
  real(DOUBLE_KIND),      allocatable :: tstamp(:)
  real, dimension(:,:,:), allocatable :: data, gdata, rdata

  call mpp_init()
  pe = mpp_pe()
  npes = mpp_npes()

!possibly open a file called mpp_io.nml
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='mpp_io.nml', err=10 )
  read( unit,mpp_io_nml )
  close(unit)
10 continue

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_io_init(MPP_DEBUG)
  else
      call mpp_io_init()
  end if
  call mpp_set_stack_size(stackmax)
  call mpp_domains_set_stack_size(stackmaxd)

  if( pe.EQ.mpp_root_pe() )then
      print '(a,6i4)', 'npes, nx, ny, nz, nt, halo=', npes, nx, ny, nz, nt, halo
      print *, 'Using NEW domaintypes and calls...'
  end if
!define global data array
  allocate( gdata(nx,ny,nz) )
  if( pe.EQ.mpp_root_pe() )then
!      call random_number(gdata) )
!fill in global array: with k.iiijjj
      gdata = 0.
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               gdata(i,j,k) = k + i*1e-3 + j*1e-6
            end do
         end do
      end do
  end if
  call mpp_broadcast( gdata, size(gdata), mpp_root_pe() )

!define domain decomposition
  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo )
  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
  call mpp_get_domain_components( domain, xdom, ydom )
  allocate( data(isd:ied,jsd:jed,nz) )
  data(is:ie,js:je,:) = gdata(is:ie,js:je,:)

!tests
  write( file,'(a,i3.3)' )trim(file), npes

!sequential write: single-threaded formatted: only if small
  if( nx*ny*nz*nt.LT.1000 )then
      if( pe.EQ.mpp_root_pe() )print *, 'sequential write: single-threaded formatted'
!here the only test is a successful write: please look at test.txt for verification.
      call mpp_open( unit, trim(file)//'s.txt', action=MPP_OVERWR, form=MPP_ASCII, threading=MPP_SINGLE )
      call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nx)/) )
      call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,ny)/) )
      call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
      call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
      call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data' )
      call mpp_write( unit, x )
      call mpp_write( unit, y )
      call mpp_write( unit, z )
      do i = 0,nt-1
         time = i*10.
         call mpp_write( unit, f, domain, data, time )
      end do
      call mpp_close(unit)
  end if

!netCDF distributed write
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF distributed write'
  call mpp_open( unit, trim(file)//'d', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nx)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,ny)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data' )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time )
  end do
  call mpp_close(unit)
  
!netCDF single-threaded write
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF single-threaded write'
  call mpp_open( unit, trim(file)//'s', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nx)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,ny)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time )
  end do
  call mpp_close(unit)

!netCDF multi-threaded read
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF multi-threaded read'
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(file)//'s', action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  allocate( atts(natt) )
  allocate( axes(ndim) )
  allocate( vars(nvar) )
  allocate( tstamp(ntime) )
  call mpp_get_atts ( unit, atts(:) )
  call mpp_get_axes ( unit, axes(:) )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_times( unit, tstamp(:) )

  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  allocate( rdata(is:ie,js:je,nz) )
  call mpp_read( unit, vars(1), domain, rdata, 1 )
  rchk = mpp_chksum(rdata(is:ie,js:je,:))
  chk  = mpp_chksum( data(is:ie,js:je,:))
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', 'checksum=', rchk, chk
  if( rchk.NE.chk )call mpp_error( FATAL, 'Checksum error on multi-threaded netCDF read.' )

  call mpp_io_exit()
  call mpp_domains_exit()
  call mpp_exit()

end program mpp_io_test

#endif
