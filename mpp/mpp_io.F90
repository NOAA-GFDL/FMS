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

!these are used to determine hardware/OS/compiler

#if defined(_CRAY) || defined(__sgi)
#define SGICRAY
#endif

#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(__sgi)
#define SGICRAY_MPP
#endif

!machines that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)
#define use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#endif

module mpp_io_mod
  use mpp_mod
  use mpp_domains_mod
  implicit none
  private
  character(len=256), private :: version='$Id: mpp_io.F90,v 5.5 2000/07/28 20:17:19 fms Exp $' !RCS ID

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, parameter, private :: stdin=100, stdout=101, stderr=102
#else
  integer, parameter, private :: stdin=5, stdout=6, stderr=0
#endif
  integer, private :: pe, npes

  type, public :: axistype
     sequence
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     character(len=8) :: cartesian
     integer :: sense           !+/-1, depth or height?
     type(domain1D), pointer :: domain !if pointer is associated, it is a distributed data axis
     real, pointer :: data(:)   !axis values (not used if time axis)
     integer :: id, did         !id is the "variable ID", did is the "dimension ID": netCDF requires 2 IDs for axes
  end type axistype

  type, public :: fieldtype
     sequence
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     real :: min, max, missing, fill, scale, add
     integer :: pack
     type(axistype), pointer :: axes(:) !axes associated with field
!size, time_axis_index redundantly hold info already contained in axes
!it's clunky and inelegant, but required so that axes can be shared among multiple files
     integer, dimension(:), pointer :: size
     integer :: time_axis_index
     integer :: id
  end type fieldtype
  type(fieldtype), public :: default_field !provided to users with default components

  type, public :: filetype
     sequence
     character(len=256) :: name
     integer :: action, format, access, threading, fileset, record, ncid
     logical :: opened, initialized, nohdrs
     integer :: time_level
     real :: time
     integer :: id              !variable ID of time axis associated with file (only one time axis per file)
  end type filetype

!action on open
  integer, parameter, public :: MPP_WRONLY=100, MPP_RDONLY=101, MPP_APPEND=102, MPP_OVERWR=103
!format
  integer, parameter, public :: MPP_ASCII=200,  MPP_IEEE32=201, MPP_NATIVE=202, MPP_NETCDF=203
!access
  integer, parameter, public :: MPP_SEQUENTIAL=300, MPP_DIRECT=301
!threading, fileset
  integer, parameter, public :: MPP_SINGLE=400, MPP_MULTI=401
!action on close
  integer, parameter, public :: MPP_DELETE=501, MPP_COLLECT=502

  type(filetype), private, allocatable :: mpp_file(:)
  integer, private :: records_per_pe
  integer, private :: maxunits, unit_begin, unit_end
  integer, private :: varnum=0
  integer, private :: error
  character(len=256) :: text
!null unit: returned by PEs not participating in IO after a collective call
  integer, parameter, private :: NULLUNIT=-1
  real, parameter, private :: NULLTIME=-1.
  logical, private :: verbose=.FALSE., debug=.FALSE., mpp_io_initialized=.FALSE.

  interface mpp_write_meta
     module procedure mpp_write_meta
     module procedure mpp_write_meta_scalar_r
     module procedure mpp_write_meta_scalar_i
     module procedure mpp_write_meta_axis
     module procedure mpp_write_meta_field
     module procedure mpp_write_meta_global
     module procedure mpp_write_meta_global_scalar_r
     module procedure mpp_write_meta_global_scalar_i
  end interface

  interface mpp_write
     module procedure mpp_write_3D_dist2D
     module procedure mpp_write_2D_dist2D
     module procedure mpp_write_2D_dist1D
     module procedure mpp_write_1D_dist1D
     module procedure mpp_write_3D
     module procedure mpp_write_2D
     module procedure mpp_write_1D
     module procedure mpp_write_0D
     module procedure mpp_write_axis
  end interface

  public :: mpp_close, mpp_flush, mpp_get_iospec, mpp_get_ncid, mpp_get_unit_range, mpp_io_init, mpp_io_exit, mpp_open, &
            mpp_set_unit_range, mpp_write, mpp_write_meta

#ifdef use_netCDF
#include <netcdf.inc>
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!               mpp_io_init: initialize parallel I/O                   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_io_init( flags, maxunit )
      integer, intent(in), optional :: flags, maxunit
!initialize IO package: initialize mpp_file array, set valid range of units for fortran IO

      if( mpp_io_initialized )return
      if( KIND(0.).NE.DOUBLE_KIND )call mpp_error( FATAL, 'MPP_IO: default reals must be 8-byte: recompile.' )
      call mpp_init(flags)           !if mpp_init has been called, this call will merely return
      pe = mpp_pe()
      npes = mpp_npes()
      call mpp_domains_init(flags)

      maxunits = 64
      if( PRESENT(maxunit) )maxunits = maxunit
      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if
!initialize default_field
      default_field%name = 'noname'
      default_field%units = 'nounits'
      default_field%longname = 'noname'
!largest possible 4-byte reals
      default_field%min = -huge(1._4)
      default_field%max =  huge(1._4)
      default_field%missing = -1e36
      default_field%fill = -1e36
      default_field%scale = 1.
      default_field%add = 0.
      default_field%pack = 1
      default_field%time_axis_index = -1 !this value will never match any index
      
!up to MAXUNITS fortran units and MAXUNITS netCDF units are supported
!file attributes (opened, format, access, threading, fileset) are saved against the unit number
!external handles to netCDF units are saved from maxunits+1:2*maxunits
      allocate( mpp_file(NULLUNIT:2*maxunits) ) !starts at NULLUNIT=-1, used by non-participant PEs in single-threaded I/O
      mpp_file(:)%name   = ' '
      mpp_file(:)%action    = -1
      mpp_file(:)%format    = -1
      mpp_file(:)%threading = -1
      mpp_file(:)%fileset   = -1
      mpp_file(:)%record    = -1
      mpp_file(:)%ncid      = -1
      mpp_file(:)%opened = .FALSE.
      mpp_file(:)%initialized = .FALSE.
      mpp_file(:)%time_level = 0
      mpp_file(:)%time = NULLTIME
      mpp_file(:)%id = -1
!NULLUNIT "file" is always single-threaded, open and initialized (to pass checks in mpp_write)
      mpp_file(NULLUNIT)%threading = MPP_SINGLE
      mpp_file(NULLUNIT)%opened = .TRUE.
      mpp_file(NULLUNIT)%initialized = .TRUE.

!set range of allowed fortran unit numbers: could be compiler-dependent (should not overlap stdin/out/err)
      call mpp_set_unit_range( 7, maxunits )

      if( pe.EQ.0 )then
          write( stdout,'(/a)' )'MPP_IO module '//trim(version)
#ifdef use_netCDF
          text = NF_INQ_LIBVERS()
          write( stdout,* )'Using netCDF library version '//trim(text)
#endif
      endif

#ifdef CRAYPVP
!we require every file to be assigned threadwise: PVPs default to global, and are reset here
      call ASSIGN( 'assign -P thread p:%', error )
#endif

      call mpp_sync()
      mpp_io_initialized = .TRUE.
      return
    end subroutine mpp_io_init

    subroutine mpp_io_exit()
      integer :: unit

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_IO_EXIT: must first call mpp_io_init.' )
!close all open fortran units
      do unit = unit_begin,unit_end
         if( mpp_file(unit)%opened )close(unit)
      end do
#ifdef use_netCDF
!close all open netCDF units
      do unit = maxunits+1,2*maxunits
         if( mpp_file(unit)%opened )error = NF_CLOSE(mpp_file(unit)%ncid)
      end do
#endif

      deallocate(mpp_file)
      mpp_io_initialized = .FALSE.
      return
    end subroutine mpp_io_exit
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                !
!                   OPENING AND CLOSING FILES: mpp_open() and mpp_close()                        !
!                                                                                                !
! mpp_open( unit, file, action, form, access, threading, fileset, iospec, nohdrs, recl, pelist ) !
!      integer, intent(out) :: unit                                                              !
!      character(len=*), intent(in) :: file                                                      !
!      integer, intent(in), optional :: action, form, access, threading, fileset, recl           !
!      character(len=*), intent(in), optional :: iospec                                          !
!      logical, intent(in), optional :: nohdrs                                                   !
!      integer, optional, intent(in) :: pelist(:) !default ALL                                   !
!                                                                                                !
!  unit is intent(OUT): always _returned_by_ mpp_open()                                          !
!  file is the filename: REQUIRED                                                                !
!    we append .nc to filename if it is a netCDF file                                            !
!    we append .<pppp> to filename if fileset is private (pppp is PE number)                     !
!  iospec is a system hint for I/O organization, e.g assign(1) on SGI/Cray systems.              !
!  if nohdrs is .TRUE. headers are not written on non-netCDF writes.                             !
!  nohdrs has no effect when action=MPP_RDONLY|MPP_APPEND or when form=MPP_NETCDF                !
! FLAGS:                                                                                         !
!    action is one of MPP_RDONLY, MPP_APPEND or MPP_WRONLY                                       !
!    form is one of MPP_ASCII:  formatted read/write                                             !
!                   MPP_NATIVE: unformatted read/write with no conversion                        !
!                   MPP_IEEE32: unformatted read/write with conversion to IEEE32                 !
!                   MPP_NETCDF: unformatted read/write with conversion to netCDF                 !
!    access is one of MPP_SEQUENTIAL or MPP_DIRECT (ignored for netCDF)                          !
!      RECL argument is REQUIRED for direct access IO                                            !
!    threading is one of MPP_SINGLE or MPP_MULTI                                                 !
!      single-threaded IO in a multi-PE run is done by PE0                                       !
!    fileset is one of MPP_MULTI and MPP_SINGLE                                                  !
!      fileset is only used for multi-threaded I/O                                               !
!      if all I/O PEs in <pelist> use a single fileset, they write to the same file              !
!      if all I/O PEs in <pelist> use a multi  fileset, they each write an independent file      !
!  recl is the record length in bytes                                                            !
!  pelist is the list of I/O PEs (currently ALL)                                                 !
!                                                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_open( unit, file, action, form, access, threading, fileset, iospec, nohdrs, recl, pelist )
      integer, intent(out) :: unit
      character(len=*), intent(in) :: file
      integer, intent(in), optional :: action, form, access, threading, fileset, recl
      character(len=*), intent(in), optional :: iospec
      logical, intent(in), optional :: nohdrs
      integer, intent(in), optional :: pelist(:) !default ALL

      character(len=16) :: act, acc, for, pos
      integer :: action_flag, form_flag, access_flag, threading_flag, fileset_flag
      logical :: exists
      character(len=64) :: filespec
      type(axistype) :: unlim    !used by netCDF with mpp_append

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_OPEN: must first call mpp_io_init.' )
!set flags
      action_flag = MPP_WRONLY        !default
      if( PRESENT(action) )action_flag = action
      form_flag = MPP_ASCII
      if( PRESENT(form) )form_flag = form
#ifndef use_netCDF
      if( form_flag.EQ.MPP_NETCDF ) &
           call mpp_error( FATAL, 'MPP_OPEN: To open a file with form=MPP_NETCDF, you must compile mpp_io with -Duse_netCDF.' )
#endif
      access_flag = MPP_SEQUENTIAL
      if( PRESENT(access) )access_flag = access
      threading_flag = MPP_SINGLE
      if( npes.GT.1 .AND. PRESENT(threading) )threading_flag = threading
      fileset_flag = MPP_MULTI
      if( PRESENT(fileset) )fileset_flag = fileset
      if( threading_flag.EQ.MPP_SINGLE )fileset_flag = MPP_SINGLE

!get a unit number
      if( threading_flag.EQ.MPP_SINGLE )then
          if( pe.NE.0 .AND. action_flag.NE.MPP_RDONLY )then
              unit = NULLUNIT           !PEs not participating in IO from this mpp_open() will return this value for unit
              return
          end if
      end if
      if( form_flag.EQ.MPP_NETCDF )then
          do unit = maxunits+1,2*maxunits
             if( .NOT.mpp_file(unit)%opened )exit
          end do
          if( unit.GT.2*maxunits )call mpp_error( FATAL, 'MPP_OPEN: too many open netCDF files.' )
      else
          do unit = unit_begin, unit_end
             inquire( unit,OPENED=mpp_file(unit)%opened )
             if( .NOT.mpp_file(unit)%opened )exit
          end do
          if( unit.GT.unit_end )call mpp_error( FATAL, 'MPP_OPEN: no available units.' )
      end if

!get a filename
      text = file
      if( form_flag.EQ.MPP_NETCDF )text = trim(file)//'.nc'
      if( fileset_flag.EQ.MPP_MULTI )write( text,'(a,i4.4)' )trim(text)//'.', pe
      mpp_file(unit)%name = text
      if( verbose )print '(a,2i3,x,a,5i5)', 'MPP_OPEN: PE, unit, filename, action, format, access, threading, fileset=', &
           pe, unit, trim(mpp_file(unit)%name), action_flag, form_flag, access_flag, threading_flag, fileset_flag

!action: read, write, overwrite, append: act and pos are ignored by netCDF
      if( action_flag.EQ.MPP_RDONLY )then
          act = 'READ'
          pos = 'REWIND'
          if( form_flag.EQ.MPP_NETCDF )call mpp_error( FATAL, 'MPP_OPEN: only writes are currently supported with netCDF.' )
      else if( action_flag.EQ.MPP_WRONLY .OR. action_flag.EQ.MPP_OVERWR )then
          act = 'WRITE'
          pos = 'REWIND'
      else if( action_flag.EQ.MPP_APPEND )then
          act = 'WRITE'
          pos = 'APPEND'
      else
          call mpp_error( FATAL, 'MPP_OPEN: action must be one of MPP_WRONLY, MPP_APPEND or MPP_RDONLY.' )
      end if

!access: sequential or direct: ignored by netCDF
      if( form_flag.NE.MPP_NETCDF )then
          if( access_flag.EQ.MPP_SEQUENTIAL )then
              acc = 'SEQUENTIAL'
          else if( access_flag.EQ.MPP_DIRECT )then
              acc = 'DIRECT'
              if( form_flag.EQ.MPP_ASCII )call mpp_error( FATAL, 'MPP_OPEN: formatted direct access I/O is prohibited.' )
              if( .NOT.PRESENT(recl) ) &
                   call mpp_error( FATAL, 'MPP_OPEN: recl (record length in bytes) must be specified with access=MPP_DIRECT.' )
              mpp_file(unit)%record = 1
              records_per_pe = 1 !each PE writes 1 record per mpp_write
          else
              call mpp_error( FATAL, 'MPP_OPEN: access must be one of MPP_SEQUENTIAL or MPP_DIRECT.' )
          end if
      end if

!threading: SINGLE or MULTI
      if( threading_flag.EQ.MPP_MULTI )then
!fileset: MULTI or SINGLE (only for multi-threaded I/O
          if( fileset_flag.EQ.MPP_SINGLE )then
              if( form_flag.EQ.MPP_NETCDF ) &
                   call mpp_error( FATAL, 'MPP_OPEN: we currently do not support single-file multi-threaded netCDF I/O.' )
#ifdef _CRAYT3E
              call ASSIGN( 'assign -I -F global.privpos f:'//trim(mpp_file(unit)%name), error )
#endif
          else if( fileset_flag.NE.MPP_MULTI )then
              call mpp_error( FATAL, 'MPP_OPEN: fileset must be one of MPP_MULTI or MPP_SINGLE.' )
          end if
      else if( threading_flag.NE.MPP_SINGLE )then
          call mpp_error( FATAL, 'MPP_OPEN: threading must be one of MPP_SINGLE or MPP_MULTI.' )
      end if

!apply I/O specs before opening the file
!note that -P refers to the scope of a fortran unit, which is always thread-private even if file is shared
#ifdef CRAYPVP
      call ASSIGN( 'assign -I -P thread  f:'//trim(mpp_file(unit)%name), error )
#endif
#ifdef _CRAYT3E
      call ASSIGN( 'assign -I -P private f:'//trim(mpp_file(unit)%name), error )
#endif
      if( PRESENT(iospec) )then
!iospec provides hints to the system on how to organize I/O
!on Cray systems this is done through 'assign', see assign(1) and assign(3F)
!on other systems this will be expanded as needed
!no error checks here on whether the supplied iospec is valid
#ifdef SGICRAY
          call ASSIGN( 'assign -I '//trim(iospec)//' f:'//trim(mpp_file(unit)%name), error )
          if( form_flag.EQ.MPP_NETCDF )then
!for netCDF on SGI/Cray systems we pass it to the environment variable NETCDF_XFFIOSPEC
!ideally we should parse iospec, pass the argument of -F to NETCDF_FFIOSPEC, and the rest to NETCDF_XFFIOSPEC
!maybe I'll get around to it someday
!PXFSETENV is a POSIX-standard routine for setting environment variables from fortran
              call PXFSETENV( 'NETCDF_XFFIOSPEC', 0, trim(iospec), 0, 1, error )
          end if
#endif
      end if

!open the file as specified above for various formats
      if( form_flag.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          if( action_flag.EQ.MPP_WRONLY )then
              error = NF_CREATE( trim(mpp_file(unit)%name), NF_NOCLOBBER, mpp_file(unit)%ncid ); call netcdf_err(error)
              if( verbose )print '(a,i3,i16)', 'MPP_OPEN: new netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_OVERWR )then
              error = NF_CREATE( trim(mpp_file(unit)%name), NF_CLOBBER,   mpp_file(unit)%ncid ); call netcdf_err(error)
              action_flag = MPP_WRONLY !after setting clobber, there is no further distinction btwn MPP_WRONLY and MPP_OVERWR
              if( verbose )print '(a,i3,i16)', 'MPP_OPEN: overwrite netCDF file: pe, ncid=', pe, mpp_file(unit)%ncid
          else if( action_flag.EQ.MPP_APPEND )then
              error = NF_OPEN( trim(mpp_file(unit)%name), NF_WRITE, mpp_file(unit)%ncid ); call netcdf_err(error)
!get the current time level of the file: writes to this file will be at next time level
              error = NF_INQ_UNLIMDIM( mpp_file(unit)%ncid, unlim%did )
              if( error.EQ.NF_NOERR )then
                  error = NF_INQ_DIM( mpp_file(unit)%ncid, unlim%did, unlim%name, mpp_file(unit)%time_level )
                  call netcdf_err(error)
                  error = NF_INQ_VARID( mpp_file(unit)%ncid, unlim%name, mpp_file(unit)%id ); call netcdf_err(error)
              end if
              if( verbose )print '(a,i3,i16,i4)', 'MPP_OPEN: append to existing netCDF file: pe, ncid, time_axis_id=',&
                   pe, mpp_file(unit)%ncid, mpp_file(unit)%id
          end if
          mpp_file(unit)%opened = .TRUE.
#endif
      else
!format: ascii, native, or IEEE 32 bit
          if( form_flag.EQ.MPP_ASCII )then
              for = 'FORMATTED'
          else if( form_flag.EQ.MPP_IEEE32 )then
              for = 'UNFORMATTED'
!assign -N is currently unsupported on SGI
#ifdef _CRAY
              call ASSIGN( 'assign -I -N ieee_32 f:'//trim(mpp_file(unit)%name), error )
#endif
          else if( form_flag.EQ.MPP_NATIVE )then
              for = 'UNFORMATTED'
          else
              call mpp_error( FATAL, 'MPP_OPEN: form must be one of MPP_ASCII, MPP_NATIVE, MPP_IEEE32 or MPP_NETCDF.' )
          end if
          inquire( file=trim(mpp_file(unit)%name), EXIST=exists )
          if( exists .AND. action_flag.EQ.MPP_WRONLY ) &
               call mpp_error( WARNING, 'MPP_OPEN: File '//trim(mpp_file(unit)%name)//' opened WRONLY already exists!' )
          if( action_flag.EQ.MPP_OVERWR )action_flag = MPP_WRONLY
!perform the OPEN here
          if( PRESENT(recl) )then
              if( verbose )print '(2(x,a,i3),5(x,a),a,i8)', 'MPP_OPEN: PE=', pe, &
                   'unit=', unit, trim(mpp_file(unit)%name), 'attributes=', trim(acc), trim(for), trim(act), ' RECL=', recl
              open( unit, file=trim(mpp_file(unit)%name), access=acc, form=for, action=act, recl=recl )
          else
              if( verbose )print '(2(x,a,i3),6(x,a))',      'MPP_OPEN: PE=', pe, &
                   'unit=', unit, trim(mpp_file(unit)%name), 'attributes=', trim(acc), trim(for), trim(pos), trim(act)
              open( unit, file=trim(mpp_file(unit)%name), access=acc, form=for, action=act, position=pos )
          end if
!check if OPEN worked
          inquire( unit,OPENED=mpp_file(unit)%opened )
          if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_OPEN: error in OPEN() statement.' )
      end if
      mpp_file(unit)%action = action_flag
      mpp_file(unit)%format = form_flag
      mpp_file(unit)%access = access_flag
      mpp_file(unit)%threading = threading_flag
      mpp_file(unit)%fileset = fileset_flag
      if( PRESENT(nohdrs) )mpp_file(unit)%nohdrs = nohdrs

      if( action_flag.EQ.MPP_WRONLY )then
          if( form_flag.NE.MPP_NETCDF .AND. access_flag.EQ.MPP_DIRECT )call mpp_write_meta( unit, 'record_length', ival=recl )
!actual file name
          call mpp_write_meta( unit, 'filename', cval=mpp_file(unit)%name )
!MPP_IO package version
          call mpp_write_meta( unit, 'MPP_IO_VERSION', cval=trim(version) )
!filecount for multifileset
          if( threading_flag.EQ.MPP_MULTI .AND. fileset_flag.EQ.MPP_MULTI ) &
               call mpp_write_meta( unit, 'NumFilesInSet', ival=npes )
      end if

      return
    end subroutine mpp_open

    subroutine mpp_close( unit, action )
      integer, intent(in) :: unit
      integer, intent(in), optional :: action
      character(len=8) :: status
      logical :: collect

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_CLOSE: must first call mpp_io_init.' )
      if( unit.EQ.NULLUNIT )return !nothing was actually opened on this unit

!action on close
      status = 'KEEP'
!collect is supposed to launch the post-processing collector tool for multi-fileset
      collect = .FALSE.
      if( PRESENT(action) )then
          if( action.EQ.MPP_DELETE )then
              status = 'DELETE'
          else if( action.EQ.MPP_COLLECT )then
              collect = .FALSE.         !should be TRUE but this is not yet ready
              call mpp_error( WARNING, 'MPP_CLOSE: the COLLECT operation is not yet implemented.' )
          else
              call mpp_error( FATAL, 'MPP_CLOSE: action must be one of MPP_DELETE or MPP_COLLECT.' )
          end if
      end if
      if( mpp_file(unit)%fileset.NE.MPP_MULTI )collect = .FALSE.
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          error = NF_CLOSE(mpp_file(unit)%ncid); call netcdf_err(error)
#endif
      else
          close(unit,status=status)
      end if
#ifdef SGICRAY
!this line deleted: since the FILENV is a shared file, this might cause a problem in
! multi-threaded I/O if one PE does assign -R before another one has opened it.
!      call ASSIGN( 'assign -R f:'//trim(mpp_file(unit)%name), error )
#endif
      mpp_file(unit)%name = ' '
      mpp_file(unit)%action    = -1
      mpp_file(unit)%format    = -1
      mpp_file(unit)%access    = -1
      mpp_file(unit)%threading = -1
      mpp_file(unit)%fileset   = -1
      mpp_file(unit)%record    = -1
      mpp_file(unit)%ncid      = -1
      mpp_file(unit)%opened = .FALSE.
      mpp_file(unit)%initialized = .FALSE.
      mpp_file(unit)%id = -1
      mpp_file(unit)%time_level = 0
      mpp_file(unit)%time = NULLTIME
      return
    end subroutine mpp_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                 !
!                             MPP_WRITE_META                                      !
!                                                                                 !
!This series of routines is used to describe the contents of the file             !
!being written on <unit>. Each file can contain any number of fields,             !
!which can be functions of 0-3 spatial axes and 0-1 time axes. Axis               !
!descriptors are stored in the <axistype> structure and field                     !
!descriptors in the <fieldtype> structure.                                        !
!                                                                                 !
!  type, public :: axistype                                                       !
!     sequence                                                                    !
!     character(len=128) :: name                                                  !
!     character(len=128) :: units                                                 !
!     character(len=256) :: longname                                              !
!     integer :: sense           !+/-1, depth or height?                          !
!     type(domain1D), pointer :: domain                                           !
!     real, pointer :: data(:) !axis values (not used if time axis)               !
!     integer :: id                                                               !
!  end type axistype                                                              !
!                                                                                 !
!  type, public :: fieldtype                                                      !
!     sequence                                                                    !
!     character(len=128) :: name                                                  !
!     character(len=128) :: units                                                 !
!     character(len=256) :: longname                                              !
!     real :: min, max, missing, fill, scale, add                                 !
!     type(axistype), pointer :: axis(:)                                          !
!     integer :: id                                                               !
!  end type fieldtype                                                             !
!                                                                                 !
!The metadata contained in the type is always written for each axis and           !
!field. Any other metadata one wishes to attach to an axis or field               !
!can subsequently be passed to mpp_write_meta using the ID, as shown below.       !
!                                                                                 !
!mpp_write_meta can take several forms:                                           !
!                                                                                 !
!  mpp_write_meta( unit, name, rval=rval, pack=pack )                             !
!  mpp_write_meta( unit, name, ival=ival )                                        !
!  mpp_write_meta( unit, name, cval=cval )                                        !
!      integer, intent(in) :: unit                                                !
!      character(len=*), intent(in) :: name                                       !
!      real, intent(in), optional :: rval(:)                                      !
!      integer, intent(in), optional :: ival(:)                                   !
!      character(len=*), intent(in), optional :: cval                             !
!                                                                                 !
!    This form defines global metadata associated with the file as a              !
!    whole. The attribute is named <name> and can take on a real, integer         !
!    or character value. <rval> and <ival> can be scalar or 1D arrays.            !
!                                                                                 !
!  mpp_write_meta( unit, id, name, rval=rval, pack=pack )                         !
!  mpp_write_meta( unit, id, name, ival=ival )                                    !
!  mpp_write_meta( unit, id, name, cval=cval )                                    !
!      integer, intent(in) :: unit, id                                            !
!      character(len=*), intent(in) :: name                                       !
!      real, intent(in), optional :: rval(:)                                      !
!      integer, intent(in), optional :: ival(:)                                   !
!      character(len=*), intent(in), optional :: cval                             !
!                                                                                 !
!    This form defines metadata associated with a previously defined              !
!    axis or field, identified to mpp_write_meta by its unique ID <id>.           !
!    The attribute is named <name> and can take on a real, integer                !
!    or character value. <rval> and <ival> can be scalar or 1D arrays.            !
!    This need not be called for attributes already contained in                  !
!    the type.                                                                    !
!                                                                                 !
!    PACK can take values 1,2,4,8. This only has meaning when writing             !
!    floating point numbers. The value of PACK defines the number of words that   !
!    are written into 8 bytes. For pack=4 and pack=8, an integer value is         !
!    written: rval is assumed to have been scaled to the appropriate dynamic      !
!    range.                                                                       !
!    PACK currently only works for netCDF files, and is ignored otherwise.        !
!                                                                                 !
!   subroutine mpp_write_meta_axis( unit, axis, name, units, longname, &          !
!        cartesian, sense, domain, data )                                         !
!     integer, intent(in) :: unit                                                 !
!     type(axistype), intent(inout) :: axis                                       !
!     character(len=*), intent(in) :: name, units, longname                       !
!     character(len=*), intent(in), optional :: cartesian                         !
!     integer, intent(in), optional :: sense                                      !
!     type(domain1D), intent(in), optional, target :: domain                      !
!     real, intent(in), optional :: data(:)                                       !
!                                                                                 !
!    This form defines a time or space axis. Metadata corresponding to the type   !
!    above are written to the file on <unit>. A unique ID for subsequent          !
!    references to this axis is returned in axis%id. If the <domain>              !
!    element is present, this is recognized as a distributed data axis            !
!    and domain decomposition information is also written if required (the        !
!    domain decomposition info is required for multi-fileset multi-threaded       !
!    I/O). If the <data> element is allocated, it is considered to be a space     !
!    axis, otherwise it is a time axis with an unlimited dimension. Only one      !
!    time axis is allowed per file.                                               !
!                                                                                 !
!   subroutine mpp_write_meta_field( unit, field, axes, name, units, longname, &  !
!        min, max, missing, fill, scale, add, pack )                              !
!     integer, intent(in) :: unit                                                 !
!     type(fieldtype), intent(out) :: field                                       !
!     type(axistype), intent(in) :: axes(:)                                       !
!     character(len=*), intent(in) :: name, units, longname                       !
!     real, intent(in), optional :: min, max, missing, fill, scale, add           !
!     integer, intent(in), optional :: pack                                       !
!                                                                                 !
!    This form defines a field. Metadata corresponding to the type                !
!    above are written to the file on <unit>. A unique ID for subsequent          !
!    references to this field is returned in field%id. At least one axis          !
!    must be associated, 0D variables are not considered. mpp_write_meta          !
!    must previously have been called on all axes associated with this            !
!    field.                                                                       !
!                                                                                 !
! The mpp_write_meta package also includes subroutines write_attribute and        !
! write_attribute_netcdf, that are private to this module.                        !
!                                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_write_meta_global( unit, name, rval, ival, cval, pack )
!writes a global metadata attribute to unit <unit>
!attribute <name> can be an real, integer or character
!one and only one of rval, ival, and cval should be present
!the first found will be used
!for a non-netCDF file, it is encoded into a string "GLOBAL <name> <val>"
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack

      if( .NOT.mpp_io_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset.EQ.MPP_SINGLE   .AND. pe.NE.0 )return
      if( mpp_file(unit)%action.NE.MPP_WRONLY )return !no writing metadata on APPEND
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          call write_attribute_netcdf( unit, NF_GLOBAL, name, rval, ival, cval, pack )
#endif
      else
          call write_attribute( unit, 'GLOBAL '//trim(name), rval, ival, cval, pack )
      end if

      return
    end subroutine mpp_write_meta_global

!versions of above to support <rval> and <ival> as scalars (because of f90 strict rank matching)
    subroutine mpp_write_meta_global_scalar_r( unit, name, rval, pack )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, intent(in) :: rval
      integer, intent(in), optional :: pack

      call mpp_write_meta_global( unit, name, rval=(/rval/), pack=pack )
      return
    end subroutine mpp_write_meta_global_scalar_r

    subroutine mpp_write_meta_global_scalar_i( unit, name, ival )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: ival

      call mpp_write_meta_global( unit, name, ival=(/ival/) )
      return
    end subroutine mpp_write_meta_global_scalar_i

    subroutine mpp_write_meta( unit, id, name, rval, ival, cval, pack )
!writes a metadata attribute for variable <id> to unit <unit>
!attribute <name> can be an real, integer or character
!one and only one of rval, ival, and cval should be present
!the first found will be used
!for a non-netCDF file, it is encoded into a string "<id> <name> <val>"
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack

      if( .NOT.mpp_io_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset.EQ.MPP_SINGLE   .AND. pe.NE.0 )return
      if( mpp_file(unit)%action.NE.MPP_WRONLY )return !no writing metadata on APPEND
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
          call write_attribute_netcdf( unit, id, name, rval, ival, cval, pack )
      else
          write( text, '(a,i4,a)' )'VARIABLE ', id, ' '//name
          call write_attribute( unit, trim(text), rval, ival, cval, pack )
      end if

      return
    end subroutine mpp_write_meta

!versions of above to support <rval> and <ival> as scalar (because of f90 strict rank matching)
    subroutine mpp_write_meta_scalar_r( unit, id, name, rval, pack )
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      real, intent(in) :: rval
      integer, intent(in), optional :: pack

      call mpp_write_meta( unit, id, name, rval=(/rval/), pack=pack )
      return
    end subroutine mpp_write_meta_scalar_r

    subroutine mpp_write_meta_scalar_i( unit, id, name, ival )
      integer, intent(in) :: unit, id
      character(len=*), intent(in) :: name
      integer, intent(in) :: ival

      call mpp_write_meta( unit, id, name, ival=(/ival/) )
      return
    end subroutine mpp_write_meta_scalar_i

    subroutine mpp_write_meta_axis( unit, axis, name, units, longname, cartesian, sense, domain, data )
!load the values in an axistype (still need to call mpp_write)
!write metadata attributes for axis
!it is declared intent(inout) so you can nullify pointers in the incoming object if needed
!the f90 standard doesn't guarantee that intent(out) on a type guarantees that its pointer components will be unassociated
      integer, intent(in) :: unit
      type(axistype), intent(inout) :: axis
      character(len=*), intent(in) :: name, units, longname
      character(len=*), intent(in), optional :: cartesian
      integer, intent(in), optional :: sense
      type(domain1D), intent(in), optional, target :: domain
      real, intent(in), optional :: data(:)
      character(len=256) :: text

      if( .NOT.mpp_io_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset.EQ.MPP_SINGLE   .AND. pe.NE.0 )return
      if( mpp_file(unit)%action.NE.MPP_WRONLY )return !no writing metadata on APPEND
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

!pre-existing pointers need to be nullified
      if( ASSOCIATED(axis%domain) )NULLIFY(axis%domain)
      if( ASSOCIATED(axis%data)   )NULLIFY(axis%data)
!load axistype
      axis%name     = name
      axis%units    = units
      axis%longname = longname
      if( PRESENT(cartesian) )axis%cartesian = cartesian
      if( PRESENT(sense)     )axis%sense     = sense
      if( PRESENT(domain)    )axis%domain => domain
      if( PRESENT(data) )then
          if( PRESENT(domain) )then
              if( size(data).NE.domain%global%size ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: size(data).NE.domain%global%size.' )
              allocate( axis%data(domain%global%start_index:domain%global%end_index) )
          else
              allocate( axis%data(size(data)) )
          end if
          axis%data = data
      end if

!write metadata
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
!write axis def
          if( ASSOCIATED(axis%data) )then !space axis
              if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. ASSOCIATED(axis%domain) )then
                  error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, axis%domain%compute%size, axis%did )
              else
                  error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, size(axis%data),          axis%did )
              end if
              call netcdf_err(error)
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_FLOAT, 1, axis%did, axis%id ); call netcdf_err(error)
          else                            !time axis
              if( mpp_file(unit)%id.NE.-1 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: There is already a time axis for this file.' )
              error = NF_DEF_DIM( mpp_file(unit)%ncid, axis%name, NF_UNLIMITED, axis%did ); call netcdf_err(error)
              error = NF_DEF_VAR( mpp_file(unit)%ncid, axis%name, NF_DOUBLE, 1, axis%did, axis%id ); call netcdf_err(error)
              mpp_file(unit)%id = axis%id !file ID is the same as time axis varID
          end if
#endif
      else
          varnum = varnum + 1
          axis%id = varnum
          axis%did = varnum
!write axis def
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' name'
          call write_attribute( unit, trim(text), cval=axis%name )
          write( text, '(a,i4,a)' )'AXIS ', axis%id, ' size'
          if( ASSOCIATED(axis%data) )then !space axis
              if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. ASSOCIATED(axis%domain) )then
                  call write_attribute( unit, trim(text), ival=(/axis%domain%compute%size/) )
              else
                  call write_attribute( unit, trim(text), ival=(/size(axis%data)/) )
              end if
          else                            !time axis
              if( mpp_file(unit)%id.NE.-1 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_AXIS: There is already a time axis for this file.' )
              call write_attribute( unit, trim(text), ival=(/0/) ) !a size of 0 indicates time axis
              mpp_file(unit)%id = axis%id
          end if
      end if
!write axis attributes
      call mpp_write_meta( unit, axis%id, 'long_name', cval=axis%longname )
      call mpp_write_meta( unit, axis%id, 'units',     cval=axis%units    )
      if( PRESENT(cartesian) )call mpp_write_meta( unit, axis%id, 'cartesian_axis', cval=axis%cartesian )
      if( PRESENT(sense) )then
          if( sense.EQ.-1 )then
              call mpp_write_meta( unit, axis%id, 'positive', cval='down' )
          else if( sense.EQ.1 )then
              call mpp_write_meta( unit, axis%id, 'positive', cval='up' )
          end if
!silently ignore values of sense other than +/-1.
      end if
      if( mpp_file(unit)%threading.EQ.MPP_MULTI .AND. mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. ASSOCIATED(axis%domain) )then
          call mpp_write_meta( unit, axis%id, 'domain_decomposition',                 &
               ival=(/ axis%domain%global%start_index,  axis%domain%global%end_index, &
                       axis%domain%compute%start_index, axis%domain%compute%end_index /) )
      end if
      if( verbose )print '(a,2i3,x,a,2i3)', 'MPP_WRITE_META: Wrote axis metadata, pe, unit, axis%name, axis%id, axis%did=', &
           pe, unit, trim(axis%name), axis%id, axis%did 

      return
    end subroutine mpp_write_meta_axis

    subroutine mpp_write_meta_field( unit, field, axes, name, units, longname, min, max, missing, fill, scale, add, pack )
!define field: must have already called mpp_write_meta(axis) for each axis
      integer, intent(in) :: unit
      type(fieldtype), intent(out) :: field
      type(axistype), intent(in) :: axes(:)
      character(len=*), intent(in) :: name, units, longname
      real, intent(in), optional :: min, max, missing, fill, scale, add
      integer, intent(in), optional :: pack
!this array is required because of f77 binding on netCDF interface
      integer, allocatable :: axis_id(:)
      real :: a, b
      integer :: i

      if( .NOT.mpp_io_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset.EQ.MPP_SINGLE   .AND. pe.NE.0 )return
      if( mpp_file(unit)%action.NE.MPP_WRONLY )return !no writing metadata on APPEND
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

!pre-existing pointers need to be nullified
      if( ASSOCIATED(field%axes) )NULLIFY(field%axes)
!fill in field metadata
      field%name = name
      field%units = units
      field%longname = longname
      allocate( field%axes(size(axes)) )
      field%axes = axes
      field%time_axis_index = -1 !this value will never match any axis index
!size is buffer area for the corresponding axis info: it is required to buffer this info in the fieldtype
!because axis might be reused in different files
      allocate( field%size(size(axes)) )
      do i = 1,size(axes)
         if( ASSOCIATED(axes(i)%data) )then !space axis
             field%size(i) = size(axes(i)%data)
         else               !time
             field%size(i) = 1
             field%time_axis_index = i
         end if
      end do
!attributes
      if( PRESENT(min) )field%min = min
      if( PRESENT(max) )field%max = max
      if( PRESENT(missing) )field%missing = missing
      if( PRESENT(fill) )field%fill = fill
      if( PRESENT(scale) )field%scale = scale
      if( PRESENT(add) )field%add = add
      
!pack is currently used only for netCDF
      field%pack = 2        !default write 32-bit floats
      if( PRESENT(pack) )then
          if( pack.NE.1 .AND. pack.NE.2 )then
              if( pack.NE.4 .AND. pack.NE.8 ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: only legal packing values are 1,2,4,8.' )
              if( .NOT.PRESENT(scale) .OR. .NOT.PRESENT(add) ) &
                   call mpp_error( FATAL, 'MPP_WRITE_META_FIELD: scale and add must be supplied when pack=4 or 8.' )
          end if
          field%pack = pack
      end if
      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          allocate( axis_id(size(field%axes)) )
          do i = 1,size(field%axes)
             axis_id(i) = field%axes(i)%did
          end do
!write field def
          if( field%pack.EQ.1 )then
              error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_DOUBLE, size(field%axes), axis_id, field%id )
          else if( field%pack.EQ.2 )then
              error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_FLOAT,  size(field%axes), axis_id, field%id )
          else if( field%pack.EQ.4 )then
              error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_SHORT,  size(field%axes), axis_id, field%id )
          else if( field%pack.EQ.8 )then
              error = NF_DEF_VAR( mpp_file(unit)%ncid, field%name, NF_BYTE,   size(field%axes), axis_id, field%id )
          end if
          call netcdf_err(error)
#endif
      else
          varnum = varnum + 1
          field%id = varnum
          if( PRESENT(pack) )call mpp_error( WARNING, 'MPP_WRITE_META: Packing is currently available only on netCDF files.' )
!write field def
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' name'
          call write_attribute( unit, trim(text), cval=field%name )
          write( text, '(a,i4,a)' )'FIELD ', field%id, ' axes'
          call write_attribute( unit, trim(text), ival=field%axes(:)%did )
      end if
!write field attributes: these names follow netCDF conventions
      call mpp_write_meta( unit, field%id, 'long_name', cval=field%longname )
      call mpp_write_meta( unit, field%id, 'units',     cval=field%units    )
!all real attributes must be written as packed
      if( PRESENT(min) .AND. PRESENT(max) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/min,max/), pack=pack )
          else
              a = nint((min-add)/scale)
              b = nint((max-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_range', rval=(/a,  b  /), pack=pack )
          end if
      else if( PRESENT(min) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_min', rval=field%min, pack=pack )
          else
              a = nint((min-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_min', rval=a, pack=pack )
          end if
      else if( PRESENT(max) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'valid_max', rval=field%max, pack=pack )
          else
              a = nint((max-add)/scale)
              call mpp_write_meta( unit, field%id, 'valid_max', rval=a, pack=pack )
          end if
      end if
      if( PRESENT(missing) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, 'missing_value', rval=field%missing, pack=pack )
          else
              a = nint((missing-add)/scale)
              call mpp_write_meta( unit, field%id, 'missing_value', rval=a, pack=pack )
          end if
      end if
      if( PRESENT(fill) )then
          if( field%pack.EQ.1 .OR. field%pack.EQ.2 )then
              call mpp_write_meta( unit, field%id, '_FillValue', rval=field%missing, pack=pack )
          else
              a = nint((fill-add)/scale)
              call mpp_write_meta( unit, field%id, '_FillValue', rval=a, pack=pack )
          end if
      end if
      if( field%pack.NE.1 .AND. field%pack.NE.2 )then
          call mpp_write_meta( unit, field%id, 'packing', ival=field%pack )
          if( PRESENT(scale) )call mpp_write_meta( unit, field%id, 'scale_factor',  rval=field%scale )
          if( PRESENT(add)   )call mpp_write_meta( unit, field%id, 'add_offset',    rval=field%add   )
      end if
      if( verbose )print '(a,2i3,x,a,i3)', 'MPP_WRITE_META: Wrote field metadata: pe, unit, field%name, field%id=', &
           pe, unit, trim(field%name), field%id 

      return
    end subroutine mpp_write_meta_field

    subroutine write_attribute( unit, name, rval, ival, cval, pack )
!called to write metadata for non-netCDF I/O
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, intent(in), optional :: rval(:)
      integer, intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
!pack is currently ignored in this routine: only used by netCDF I/O
      integer, intent(in), optional :: pack

      if( mpp_file(unit)%nohdrs )return
!encode text string
      if( PRESENT(rval) )then
          write( text,* )trim(name)//'=', rval
      else if( PRESENT(ival) )then
          write( text,* )trim(name)//'=', ival
      else if( PRESENT(cval) )then
          text = ' '//trim(name)//'='//trim(cval)
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE: one of rval, ival, cval must be present.' )
      end if
      if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
          write( unit,fmt='(a)' )trim(text)//char(10)
      else                      !MPP_IEEE32 or MPP_NATIVE
          if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
              write(unit)trim(text)//char(10)
          else                  !MPP_DIRECT
              write( unit,rec=mpp_file(unit)%record )trim(text)//char(10)
              if( verbose )print '(a,i3,a,i3)', 'WRITE_ATTRIBUTE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              mpp_file(unit)%record = mpp_file(unit)%record + 1
          end if
      end if
      return
    end subroutine write_attribute

    subroutine write_attribute_netcdf( unit, id, name, rval, ival, cval, pack )
!called to write metadata for netCDF I/O
      integer, intent(in) :: unit
      integer, intent(in) :: id
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack
      integer :: lenc
      integer, allocatable :: rval_i(:)
#ifdef use_netCDF
      if( PRESENT(rval) )then
!pack is only meaningful for FP numbers
          if( PRESENT(pack) )then
              if( pack.EQ.1 )then
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_DOUBLE, size(rval), rval ); call netcdf_err(error)
              else if( pack.EQ.2 )then
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval), rval ); call netcdf_err(error)
              else if( pack.EQ.4 )then
                  allocate( rval_i(size(rval)) )
                  rval_i = rval
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_SHORT,  size(rval_i), rval ); call netcdf_err(error)
                  deallocate(rval_i)
              else if( pack.EQ.8 )then
                  allocate( rval_i(size(rval)) )
                  rval_i = rval
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_BYTE,   size(rval_i), rval ); call netcdf_err(error)
                  deallocate(rval_i)
              else
                  call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: only legal packing values are 1,2,4,8.' )
              end if
          else
!default is to write FLOATs (32-bit)
              error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT,  size(rval), rval ); call netcdf_err(error)
          end if
      else if( PRESENT(ival) )then
          error = NF_PUT_ATT_INT( mpp_file(unit)%ncid, id, name, NF_INT, size(ival), ival ); call netcdf_err(error)
      else if( present(cval) )then
          error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, name, len_trim(cval), cval ); call netcdf_err(error)
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: one of rval, ival, cval must be present.' )
      end if
#endif use_netCDF 
      return
    end subroutine write_attribute_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                             MPP_WRITE                                !
!                                                                      !
! mpp_write is used to write data to the file on <unit> using the      !
! file parameters supplied by mpp_open(). Axis and field definitions   !
! must have previously been written to the file using mpp_write_meta.  !
!                                                                      !
! mpp_write can take 2 forms, one for distributed data and one for     !
! non-distributed data. Distributed data refer to arrays whose two     !
! fastest-varying indices are domain-decomposed. Distributed data      !
! must be 2D or 3D (in space). Non-distributed data can be 0-3D.       !
!                                                                      !
! In all calls to mpp_write, tstamp is an optional argument. It is to  !
! be omitted if the field was defined not to be a function of time.    !
! Results are unpredictable if the argument is supplied for a time-    !
! independent field, or omitted for a time-dependent field. Repeated   !
! writes of a time-independent field are also not recommended. One     !
! time level of one field is written per call.                         !
!                                                                      !
!                                                                      !
! For non-distributed data, use                                        !
!                                                                      !
!  mpp_write( unit, field, data, tstamp )                              !
!     integer, intent(in) :: unit                                      !
!     type(fieldtype), intent(in) :: field                             !
!     real, optional :: tstamp                                         !
!     data is real and can be scalar or of rank 1-3.                   !
!                                                                      !
! For distributed data, use                                            !
!                                                                      !
!  mpp_write( unit, field, domain, data, tstamp )                      !
!     integer, intent(in) :: unit                                      !
!     type(fieldtype), intent(in) :: field                             !
!     type(domain2D), intent(in) :: domain                             !
!     real, optional :: tstamp                                         !
!     data is real and can be of rank 2 or 3.                          !
!                                                                      !
!  mpp_write( unit, axis )                                             !
!     integer, intent(in) :: unit                                      !
!     type(axistype), intent(in) :: axis                               !
!                                                                      !
! This call writes the actual co-ordinate values along each space      !
! axis. It must be called once for each space axis after all other     !
! metadata has been written.                                           !
!                                                                      !
! The mpp_write package also includes the routine write_record which   !
! performs the actual write. This routine is private to this module.   !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_write_3D_dist2D( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      real, intent(inout) :: data(domain%x%data%start_index:,domain%y%data%start_index:,:)
      real, intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data
      real :: cdata(domain%x%compute%start_index:domain%x%compute%end_index, &
                    domain%y%compute%start_index:domain%y%compute%end_index,size(data,3))
!global_domain and gdata are used to globalize data for multi-PE single-threaded I/O
!      type(domain2D), allocatable :: global_domain(:)
      real, allocatable :: gdata(:,:,:)

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( domain%x%data%is_global .AND. domain%y%data%is_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(domain%x%global%start_index:domain%x%global%end_index, &
                              domain%y%global%start_index:domain%y%global%end_index,size(data,3)) )
              call mpp_get_global( domain, data, gdata )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else
!store compute domain as contiguous data and pass to write_record
          cdata(:,:,:) = data(domain%x%compute%start_index:domain%x%compute%end_index, &
                              domain%y%compute%start_index:domain%y%compute%end_index,:)
          call write_record( unit, field, size(cdata), cdata, tstamp, domain )
      end if

      return
    end subroutine mpp_write_3D_dist2D

    subroutine mpp_write_2D_dist2D( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      real, intent(inout) :: data(:,:)
      real, intent(in), optional :: tstamp
      real :: data_3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data_3D )
      ptr = LOC(data)
#else
      data_3D(:,:,1) = data(:,:)
#endif
      call mpp_write_3D_dist2D( unit, field, domain, data_3D, tstamp )
      return
    end subroutine mpp_write_2D_dist2D

    subroutine mpp_write_2D_dist1D( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      real, intent(inout) :: data(domain%data%start_index:,:)
      real, intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data
      real :: cdata(domain%compute%start_index:domain%compute%end_index,size(data,2))
!global_domain and gdata are used to globalize data for multi-PE single-threaded I/O
      type(domain1D), allocatable :: global_domain(:)
      type(domain2D) :: write_domain(1)
      real, allocatable :: gdata(:,:)

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( domain%data%is_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( global_domain(0:domain%ndomains-1) )
              call mpp_define_domains( (/domain%global%start_index,domain%global%end_index/), &
                   global_domain, flags=GLOBAL_DATA_DOMAIN, pelist=domain%pelist, extent=domain%sizelist )
              
              allocate( gdata(domain%global%start_index:domain%global%end_index,size(data,2)) )
              gdata(domain%compute%start_index:domain%compute%end_index,:) = &
               data(domain%compute%start_index:domain%compute%end_index,:)
              call mpp_set_halo_size(domain%compute%max_size*size(data,2))
              call mpp_update_domains( gdata, global_domain(pe) )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else
!store compute domain as contiguous data and pass to write_record
          cdata(:,:) = data(domain%compute%start_index:domain%compute%end_index,:)
!write_domain is a fake 2D domain for passing to write_record
!its x axis is <domain>
!its y axis is the global undistributed second axis
          write_domain(1)%x = domain
          call mpp_define_domains( (/1,size(data,2)/), write_domain%y, flags=GLOBAL_COMPUTE_DOMAIN )
          call write_record( unit, field, size(cdata), cdata, tstamp, write_domain(1) )
      end if

      return
    end subroutine mpp_write_2D_dist1D

    subroutine mpp_write_1D_dist1D( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      real, intent(inout) :: data(:)
      real, intent(in), optional :: tstamp
      real :: data2D(size(data,1),1)

      data2D(:,1) = data(:)
      call mpp_write_2D_dist1D( unit, field, domain, data2D, tstamp )
      return
    end subroutine mpp_write_1D_dist1D

    subroutine mpp_write_3D( unit, field, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      real, intent(in) :: data(:,:,:)
      real, intent(in), optional :: tstamp

      call write_record( unit, field, size(data), data, tstamp )
    end subroutine mpp_write_3D

    subroutine mpp_write_2D( unit, field, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      real, intent(in) :: data(:,:)
      real, intent(in), optional :: tstamp

      call write_record( unit, field, size(data), data, tstamp )
    end subroutine mpp_write_2D

    subroutine mpp_write_1D( unit, field, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      real, intent(in) :: data(:)
      real, intent(in), optional :: tstamp

      call write_record( unit, field, size(data), data, tstamp )
    end subroutine mpp_write_1D

    subroutine mpp_write_0D( unit, field, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      real, intent(in) :: data
      real, intent(in), optional :: tstamp

      call write_record( unit, field, 1, (/data/), tstamp )
    end subroutine mpp_write_0D

    subroutine mpp_write_axis( unit, axis )
      integer, intent(in) :: unit
      type(axistype), intent(in) :: axis
      type(fieldtype) :: field

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset  .EQ.MPP_SINGLE .AND. pe.NE.0 )return
!we convert axis to type(fieldtype) in order to call write_record
      field = default_field
      allocate( field%axes(1) )
      field%axes(1) = axis
      allocate( field%size(1) )
      field%id = axis%id
      if( mpp_file(unit)%fileset.EQ.MPP_MULTI .AND. ASSOCIATED(axis%domain) )then
          field%size(1) = axis%domain%compute%size
          call write_record( unit, field, axis%domain%compute%size, axis%data(axis%domain%compute%start_index:) )
      else
          field%size(1) = size(axis%data)
          call write_record( unit, field, size(axis%data),          axis%data )
      end if
      return
    end subroutine mpp_write_axis

    subroutine write_record( unit, field, nwords, data, time_in, domain )
!routine that is finally called by all mpp_write routines to perform the write
!a non-netCDF record contains:
!      field ID
!      a set of 4 coordinates (is:ie,js:je) giving the data subdomain
!      a timelevel and a timestamp (=NULLTIME if field is static)
!      3D real data (stored as 1D)
!if you are using direct access I/O, the RECL argument to OPEN must be large enough for the above
!in a global direct access file, record position on PE is given by %record.

!Treatment of timestamp:
!   We assume that static fields have been passed without a timestamp.
!   Here that is converted into a timestamp of NULLTIME.
!   For non-netCDF fields, field is treated no differently, but is written
!   with a timestamp of NULLTIME. There is no check in the code to prevent
!   the user from repeatedly writing a static field.

      integer, intent(in) :: unit, nwords
      type(fieldtype), intent(in) :: field
      real, intent(in) :: data(nwords)
      real, intent(in), optional :: time_in
      type(domain2D), intent(in), optional :: domain
      integer, dimension(size(field%axes)) :: start, axsiz
      real :: time
      integer :: time_level
      logical :: newtime
      integer :: subdomain(4)
      integer :: packed_data(nwords)
#ifdef __sgi
      real(FLOAT_KIND) :: data_r4(nwords)
#endif
      integer :: i

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return
      if( mpp_file(unit)%fileset  .EQ.MPP_SINGLE .AND. pe.NE.0 )return

      if( .NOT.mpp_file(unit)%initialized )then
!this is the first call to mpp_write
!we now declare the file to be initialized
!if this is netCDF we switch file from DEFINE mode to DATA mode
          if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
!NOFILL is probably required for parallel: any circumstances in which not advisable?
              error = NF_SET_FILL( mpp_file(unit)%ncid, NF_NOFILL, i ); call netcdf_err(error)
              if( mpp_file(unit)%action.EQ.MPP_WRONLY )error = NF_ENDDEF(mpp_file(unit)%ncid); call netcdf_err(error)
#endif
          else
              call mpp_write_meta( unit, 'END', cval='metadata' )
          end if
          mpp_file(unit)%initialized = .TRUE.
          if( verbose )print '(a,i3,a)', 'MPP_WRITE: PE=', pe, ' initialized file '//trim(mpp_file(unit)%name)//'.'
      end if

!initialize time: by default assume NULLTIME
      time = NULLTIME
      time_level = -1
      newtime = .FALSE.
      if( PRESENT(time_in) )time = time_in
!increment time level if new time
      if( time.GT.mpp_file(unit)%time+EPSILON(time) )then !new time
          mpp_file(unit)%time_level = mpp_file(unit)%time_level + 1
          mpp_file(unit)%time = time
          newtime = .TRUE.
      end if
      if( verbose )print '(a,2i3,2i5,es13.5)', 'MPP_WRITE: PE, unit, %id, %time_level, %time=',&
           pe, unit, mpp_file(unit)%id, mpp_file(unit)%time_level, mpp_file(unit)%time

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
!define netCDF data block to be written:
!  time axis: START = time level
!             AXSIZ = 1
!  space axis: if there is no domain info
!              START = 1
!              AXSIZ = field%size(axis)
!          if there IS domain info:
!              start of domain is compute%start_index for multi-file I/O
!                                 global%start_index for all other cases
!              this number must be converted to 1 for NF_PUT_VAR
!                  (netCDF fortran calls are with reference to 1),
!          So, START = compute%start_index - <start of domain> + 1
!              AXSIZ = usually compute%size
!          However, if compute%start_index-compute%end_index+1.NE.compute%size,
!              we assume that the call is passing a subdomain.
!              To pass a subdomain, you must pass a domain2D object that satisfies the following:
!                  global%start_index must contain the <start of domain> as defined above;
!                  the data domain and compute domain must refer to the subdomain being passed.
!              In this case, START = compute%start_index - <start of domain> + 1
!                            AXSIZ = compute%start_index - compute%end_index + 1
! NOTE: passing of subdomains will fail for multi-PE single-threaded I/O,
!       since that attempts to gather all data on PE 0.
          start = 1
          do i = 1,size(field%axes)
             axsiz(i) = field%size(i)
             if( i.EQ.field%time_axis_index )start(i) = mpp_file(unit)%time_level
          end do
          if( PRESENT(domain) )then
              axsiz(1) = domain%x%compute%size
              axsiz(2) = domain%y%compute%size
              if( npes.GT.1 .AND. mpp_file(unit)%fileset.EQ.MPP_SINGLE )then
                  start(1) = domain%x%compute%start_index - domain%x%global%start_index + 1
                  start(2) = domain%y%compute%start_index - domain%y%global%start_index + 1
              else
                  if( domain%x%compute%end_index-domain%x%compute%start_index+1.NE.domain%x%compute%size )then
                      start(1) = domain%x%compute%start_index - domain%x%global%start_index + 1
                      axsiz(1) = domain%x%compute%end_index - domain%x%compute%start_index + 1
                  end if
                  if( domain%y%compute%end_index-domain%y%compute%start_index+1.NE.domain%y%compute%size )then
                      start(2) = domain%y%compute%start_index - domain%y%global%start_index + 1
                      axsiz(2) = domain%y%compute%end_index - domain%y%compute%start_index + 1
                  end if
              end if
          end if
          if( debug )print '(a,2i3,12i4)', 'WRITE_RECORD: PE, unit, start, axsiz=', pe, unit, start, axsiz
#ifdef use_netCDF
!write time information if new time
          if( newtime )then
              error = NF_PUT_VAR1_DOUBLE( mpp_file(unit)%ncid, mpp_file(unit)%id, mpp_file(unit)%time_level, time )
              call netcdf_err(error)
          end if
          if( field%pack.LE.2 )then
              error = NF_PUT_VARA_DOUBLE( mpp_file(unit)%ncid, field%id, start, axsiz, data        ); call netcdf_err(error)
          else              !convert to integer using scale and add: no error check on packed data representation
              packed_data = nint((data-field%add)/field%scale)
              error = NF_PUT_VARA_INT   ( mpp_file(unit)%ncid, field%id, start, axsiz, packed_data )
          end if
#endif
      else                      !non-netCDF
!subdomain contains (/is,ie,js,je/)
          if( PRESENT(domain) )then
              subdomain = (/ domain%x%compute%start_index, domain%x%compute%end_index, &
                             domain%y%compute%start_index, domain%y%compute%end_index /)
          else
              subdomain = -1    ! -1 means use global value from axis metadata
          end if
          if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
              write( unit,* )field%id, subdomain, time_level, time, data
          else                      !MPP_IEEE32 or MPP_NATIVE
              if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
#ifdef __sgi
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write(unit)field%id, subdomain, time_level, time, data_r4
                  else
                      write(unit)field%id, subdomain, time_level, time, data
                  end if
#else
                  write(unit)field%id, subdomain, time_level, time, data
#endif
              else                  !MPP_DIRECT
#ifdef __sgi
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data_r4
                  else
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
                  end if
#else
                  write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
#endif
                  if( debug )print '(a,i3,a,i3)', 'MPP_WRITE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              end if
          end if
      end if

!recompute current record for direct access I/O
      if( mpp_file(unit)%access.EQ.MPP_DIRECT )then
          if( mpp_file(unit)%fileset.EQ.MPP_SINGLE )then
!assumes all PEs participate in I/O: modify later
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe*npes
          else
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe
          end if
      end if

      return
    end subroutine write_record

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!    mpp_get_filespec, mpp_set_filespec: OS-dependent filespecs        !
!         on SGICRAY this is currently done through assign(3F).        !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_flush(unit)
!flush the output on a unit, syncing with disk
      integer, intent(in) :: unit

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_FLUSH: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_FLUSH: invalid unit number.' )
      if( .NOT.mpp_file(unit)%initialized )call mpp_error( FATAL, 'MPP_FLUSH: cannot flush a file during writing of metadata.' )
      if( mpp_file(unit)%threading.EQ.MPP_SINGLE .AND. pe.NE.0 )return

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          error = NF_SYNC(mpp_file(unit)%ncid); call netcdf_err(error)
#endif
      else
          call FLUSH(unit)
      end if
      return
    end subroutine mpp_flush

    subroutine mpp_get_iospec( unit, iospec )
      integer, intent(in) :: unit
      character(len=*), intent(out) :: iospec

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_GET_IOSPEC: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_GET_IOSPEC: invalid unit number.' )
#ifdef SGICRAY
!currently will write to stdout: don't know how to trap and return as string to iospec
      call ASSIGN( 'assign -V f:'//trim(mpp_file(unit)%name), error )
#endif
      return
    end subroutine mpp_get_iospec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!         netCDF-specific routines: mpp_get_ncid, netcdf_error         !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function mpp_get_ncid(unit)
      integer :: mpp_get_ncid
      integer, intent(in) :: unit

      mpp_get_ncid = mpp_file(unit)%ncid
      return
    end function mpp_get_ncid

    subroutine netcdf_err(err)
      integer, intent(in) :: err
      character(len=80) :: errmsg
      integer :: unit

#ifdef use_netCDF
      if( err.EQ.NF_NOERR )return
      errmsg = NF_STRERROR(err)
      call mpp_io_exit()        !make sure you close all open files
      call mpp_error( FATAL, 'NETCDF ERROR: '//trim(errmsg) )
#endif
      return
    end subroutine netcdf_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!       minor routines: mpp_get_unit_range, mpp_set_unit_range         !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_get_unit_range( unit_begin_out, unit_end_out )
      integer, intent(out) ::      unit_begin_out, unit_end_out

      unit_begin_out = unit_begin; unit_end_out = unit_end
      return
    end subroutine mpp_get_unit_range

    subroutine mpp_set_unit_range( unit_begin_in, unit_end_in )
      integer, intent(in) ::       unit_begin_in, unit_end_in

      if( unit_begin_in.GT.unit_end_in )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_begin_in.GT.unit_end_in.' )
      if( unit_begin_in.LT.0           )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_begin_in.LT.0.' )
      if( unit_end_in  .GT.maxunits    )call mpp_error( FATAL, 'MPP_SET_UNIT_RANGE: unit_end_in.GT.maxunits.' )
      unit_begin = unit_begin_in; unit_end = unit_end_in
      return
    end subroutine mpp_set_unit_range

end module mpp_io_mod
