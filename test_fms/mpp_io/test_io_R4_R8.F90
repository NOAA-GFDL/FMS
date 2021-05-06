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
!> @file
!> @brief unit test for mpp_write and mpp_read
!> @email gfdl.climate.model.info@noaa.gov
!> @description Tests mpp_write and mpp_read for reads/writes
!>  with mixed precision reals on non-mosaic files
program test_io_R4_R8

  use platform_mod,    only : r4_kind, r8_kind, i8_kind
  use mpp_mod,         only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sync_self
  use mpp_mod,         only : FATAL, NOTE, mpp_chksum, MPP_DEBUG, mpp_set_stack_size, MPP_CLOCK_SYNC
  use mpp_mod,         only : mpp_sync, mpp_exit, mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_domains_mod, only : mpp_define_domains, mpp_domains_set_stack_size, domain1D, mpp_get_global_domain
  use mpp_domains_mod, only : domain2D, mpp_define_layout, mpp_get_domain_components, mpp_define_mosaic
  use mpp_domains_mod, only : mpp_get_memory_domain, mpp_get_compute_domain, mpp_domains_exit
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_define_io_domain, mpp_deallocate_domain
  use mpp_io_mod,      only : mpp_io_init, mpp_write_meta, axistype, fieldtype, atttype
  use mpp_io_mod,      only : MPP_RDONLY, mpp_open, MPP_OVERWR, MPP_ASCII, MPP_SINGLE
  use mpp_io_mod,      only : MPP_NETCDF, MPP_MULTI, mpp_get_atts, mpp_write, mpp_close
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_get_fields, mpp_get_times
  use mpp_io_mod,      only : mpp_read, mpp_io_exit, MPP_APPEND

#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#endif

  implicit none

#ifdef use_netCDF
#include <netcdf.inc>
#endif

  !--- namelist definition
  integer           :: nx=360, ny=200, nz=50, nt=2
  integer           :: halo=2, stackmax=1500000, stackmaxd=2000000
  logical           :: debug=.FALSE.
  character(len=64) :: file='test', iospec='-F cachea'
  integer           :: layout(2) = (/1,1/)
  integer           :: ntiles_x=1, ntiles_y=1  ! total number of tiles will be ntiles_x*ntiles_y,
                                               ! the grid size for each tile will be (nx/ntiles_x, ny/ntiles_y)
                                               ! set ntiles > 1 to test the efficiency of mpp_io.
  integer           :: io_layout(2) = (/1,1/)  ! set io_layout to divide each tile into io_layout(1)*io_layout(2)
                                               ! group and write out data from the root pe of each group.
  integer           :: pack_size = 1

  namelist / test_mpp_io_R8_nml / nx, ny, nz, nt, halo, stackmax, stackmaxd, debug, file, iospec, &
                               ntiles_x, ntiles_y, layout, io_layout

  integer        :: pe, npes, io_status
  type(domain2D) :: domain

  integer            :: tks_per_sec
  integer            :: i,j,k, unit=7
  integer            :: id_single_tile_mult_file
  integer            :: id_mult_tile, id_single_tile_with_group, id_mult_tile_with_group
  logical            :: opened
  character(len=64)  :: varname

  real(r4_kind)               :: time4
  real(r8_kind)              :: time8
  type(axistype)     :: x, y, z, t
  type(fieldtype)    :: f
  type(domain1D)     :: xdom, ydom
  integer(i8_kind) :: rchk, chk
  real(r8_kind)                  :: doubledata = 0.0
  real(r8_kind)                  :: realarray(4)

! initialize modules and set up namelist
  call mpp_init()
  pe = mpp_pe()
  npes = mpp_npes()

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, test_mpp_io_R8_nml, iostat=io_status)
#else
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, file='input.nml', iostat=io_status)
  read( unit,test_mpp_io_R8_nml, iostat=io_status )
  close(unit)
#endif

  if (io_status > 0) then
      call mpp_error(FATAL,'=>test_mpp_io_R8: Error reading input.nml')
  endif

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_io_init(MPP_DEBUG)
  else
      call mpp_io_init()
  end if
  call mpp_set_stack_size(stackmax)
  call mpp_domains_set_stack_size(stackmaxd)

  write( file,'(a,i3.3)' )trim(file), npes

! determine the pack_size
  pack_size = size(transfer(doubledata, realarray))
  if( pack_size .NE. 1 .AND. pack_size .NE. 2) call mpp_error(FATAL,'test_mpp_io_R8: pack_size should be 1 or 2')

! test read/writes for different symmetries and positions and

  if(ntiles_x == 1 .and. ntiles_y == 1 .and. io_layout(1) == 1 .and. io_layout(2) == 1) then
     call test_netcdf_io_R4('Simple')
     call test_netcdf_io_R4('Symmetry_T_cell')
     call test_netcdf_io_R4('Symmetry_E_cell')
     call test_netcdf_io_R4('Symmetry_N_cell')
     call test_netcdf_io_R4('Symmetry_C_cell')
     call test_netcdf_io_R4('Symmetry_T_cell_memory')
     call test_netcdf_io_R4('Symmetry_E_cell_memory')
     call test_netcdf_io_R4('Symmetry_N_cell_memory')
     call test_netcdf_io_R4('Symmetry_C_cell_memory')
     call test_netcdf_io_R8('Simple')
     call test_netcdf_io_R8('Symmetry_T_cell')
     call test_netcdf_io_R8('Symmetry_E_cell')
     call test_netcdf_io_R8('Symmetry_N_cell')
     call test_netcdf_io_R8('Symmetry_C_cell')
     call test_netcdf_io_R8('Symmetry_T_cell_memory')
     call test_netcdf_io_R8('Symmetry_E_cell_memory')
     call test_netcdf_io_R8('Symmetry_N_cell_memory')
     call test_netcdf_io_R8('Symmetry_C_cell_memory')
  else
     call mpp_error( FATAL, 'test_mpp_io_R8: Invalid nml parameters for non-mosaic files')
  endif

  call mpp_io_exit()
  call mpp_domains_exit()
  call mpp_exit()

  contains

  !------------------------------------------------------------------
  !> Tests mpp_ reads and writes on netcdf files for both 32 and 64-bit reals
  subroutine test_netcdf_io_R4(type)
  character(len=*), intent(in) :: type
  integer :: ndim, nvar, natt, ntime
  integer :: is, ie, js, je, isd, ied, jsd, jed, ism, iem, jsm, jem
  integer :: position, msize(2), ioff, joff, nxg, nyg
  logical :: symmetry
  type(atttype),          allocatable :: atts(:)
  type(fieldtype),        allocatable :: vars(:)
  type(axistype),         allocatable :: axes(:)
  real(r8_kind),                   allocatable :: tstamp(:)
  real(r4_kind), dimension(:,:,:), allocatable  :: data4, gdata4, rdata4
  !--- determine the shift and symmetry according to type,
  select case(type)
  case('Simple')
     position = CENTER; symmetry = .false.
  case('Symmetry_T_cell', 'Symmetry_T_cell_memory')
     position = CENTER; symmetry = .true.
  case('Symmetry_E_cell', 'Symmetry_E_cell_memory')
     position = EAST;   symmetry = .true.
  case('Symmetry_N_cell', 'Symmetry_N_cell_memory')
     position = NORTH;  symmetry = .true.
  case('Symmetry_C_cell', 'Symmetry_C_cell_memory')
     position = CORNER; symmetry = .true.
  case default
     call mpp_error(FATAL, "type = "//type//" is not a valid test type")
  end select

!define domain decomposition
  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  if(index(type,"memory") == 0) then
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry )
  else  ! on memory domain
     msize(1) = nx/layout(1) + 2*halo + 2
     msize(2) = ny/layout(2) + 2*halo + 2
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry, &
                              memory_size = msize )
  end if

  call mpp_get_compute_domain( domain, is,  ie,  js,  je, position=position  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, position=position )
  call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
  call mpp_get_global_domain ( domain, xsize=nxg, ysize=nyg, position=position )
  call mpp_get_domain_components( domain, xdom, ydom )

!define global data arrays
  allocate( gdata4(nxg,nyg,nz) )
  gdata4 = 0.
  do k = 1,nz
     do j = 1,nyg
        do i = 1,nxg
           gdata4(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do

  ioff = ism - isd; joff = jsm - jsd
  allocate( data4(ism:iem,jsm:jem,nz) )
  data4 = 0
  data4(is+ioff:ie+ioff,js+joff:je+joff,:) = gdata4(is:ie,js:je,:)

! test single thread ascii writes
 if( nx*ny*nz*nt.LT.1000 .AND. index(type,"memory") .NE. 0) then
  if( index(type,"memory") .NE. 0 )then
      call mpp_open( unit, trim(file)//'sR4.txt', action=MPP_OVERWR, form=MPP_ASCII, threading=MPP_SINGLE )
      call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
      call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
      call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
      call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
      call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data' )
      call mpp_write( unit, x )
      call mpp_write( unit, y )
      call mpp_write( unit, z )
      do i = 0,nt-1
         time4 = i*10.
         call mpp_write( unit, f, domain, data4, time4)
      end do
      call mpp_close(unit)
  end if
 end if
!netCDF distributed write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'dR4', action=MPP_OVERWR, &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time4 = i*10.
     call mpp_write( unit, f, domain, data4, time4 )
  end do
  call mpp_close(unit)

!netCDF single-threaded write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR4', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )

  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )

  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )

  do i = 0,nt-1
     time4 = i*10.
     call mpp_write( unit, f, domain, data4, time4)
  end do
  call mpp_close(unit)
!reopen and test appending write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR4', action=MPP_APPEND, form=MPP_NETCDF, threading=MPP_SINGLE )
  call mpp_write( unit, f, domain, data4, time4)
  call mpp_close( unit )

! clear out for reads
  allocate( rdata4(is:ie,js:je,nz) )

!netCDF multi-threaded read
  call mpp_sync()
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR4', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
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
  call mpp_read( unit, vars(1), domain, rdata4, 1 )
  ! compare read and stored chksums
  rchk = mpp_chksum(rdata4(is:ie,js:je,:))
  chk  = mpp_chksum( data4(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//' R4: single-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'R4 Checksum error on multi-threaded/single-fileset netCDF read for type ' &
               //trim(type) )
  end if
  call mpp_close(unit)
  deallocate( atts, axes, vars, tstamp )

!netCDF distributed read
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(type)//"_"//trim(file)//'dR4', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
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
  rdata4 = 0

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )

  call mpp_read( unit, vars(1), domain, rdata4, 1 )
  ! compare read and stored chksums
  rchk = mpp_chksum(rdata4(is:ie,js:je,:))
  chk  = mpp_chksum( data4(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//' R4: multi-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'R4 Checksum error on multi-threaded/multi-fileset netCDF read for type ' &
           //trim(type) )
  end if
  deallocate( atts, axes, vars, tstamp )
  deallocate( rdata4, gdata4, data4)

  end subroutine test_netcdf_io_R4

  subroutine test_netcdf_io_R8(type)
  character(len=*), intent(in) :: type
  integer :: ndim, nvar, natt, ntime
  integer :: is, ie, js, je, isd, ied, jsd, jed, ism, iem, jsm, jem
  integer :: position, msize(2), ioff, joff, nxg, nyg
  logical :: symmetry
  type(atttype),          allocatable :: atts(:)
  type(fieldtype),        allocatable :: vars(:)
  type(axistype),         allocatable :: axes(:)
  real(r8_kind),                   allocatable :: tstamp8(:)
  real(r8_kind), dimension(:,:,:), allocatable :: data8, gdata8, rdata8
  !--- determine the shift and symmetry according to type,
  select case(type)
  case('Simple')
     position = CENTER; symmetry = .false.
  case('Symmetry_T_cell', 'Symmetry_T_cell_memory')
     position = CENTER; symmetry = .true.
  case('Symmetry_E_cell', 'Symmetry_E_cell_memory')
     position = EAST;   symmetry = .true.
  case('Symmetry_N_cell', 'Symmetry_N_cell_memory')
     position = NORTH;  symmetry = .true.
  case('Symmetry_C_cell', 'Symmetry_C_cell_memory')
     position = CORNER; symmetry = .true.
  case default
     call mpp_error(FATAL, "type = "//type//" is not a valid test type")
  end select

!define domain decomposition
  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  if(index(type,"memory") == 0) then
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry )
  else  ! on memory domain
     msize(1) = nx/layout(1) + 2*halo + 2
     msize(2) = ny/layout(2) + 2*halo + 2
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, symmetry = symmetry, &
                              memory_size = msize )
  end if

  call mpp_get_compute_domain( domain, is,  ie,  js,  je, position=position  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, position=position )
  call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )
  call mpp_get_global_domain ( domain, xsize=nxg, ysize=nyg, position=position )
  call mpp_get_domain_components( domain, xdom, ydom )

!define global data arrays for 64 bit reals
  allocate( gdata8(nxg,nyg,nz) )
  gdata8 = 0.
  do k = 1,nz
     do j = 1,nyg
        do i = 1,nxg
           gdata8(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do

  ioff = ism - isd; joff = jsm - jsd
  allocate( data8(ism:iem,jsm:jem,nz) )
  data8 = 0
  data8(is+ioff:ie+ioff,js+joff:je+joff,:) = gdata8(is:ie,js:je,:)


!sequential write: single-threaded formatted: only if small
! test ascii writes
  if( nx*ny*nz*nt.LT.1000 .AND. index(type,"memory") .NE. 0 )then
!here the only test is a successful write: please look at test.txt for verification.
      call mpp_open( unit, trim(file)//'sR8.txt', action=MPP_OVERWR, form=MPP_ASCII, threading=MPP_SINGLE )
      call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
      call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
      call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance',              data=(/(i-1.,i=1,nz)/) )
      call mpp_write_meta( unit, t, 'T', 'sec', 'Time' )
      call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data' )
      call mpp_write( unit, x )
      call mpp_write( unit, y )
      call mpp_write( unit, z )
      do i = 0,nt-1
         time8 = i*10.
         call mpp_write( unit, f, domain, data8, time8)
      end do
      call mpp_close(unit)
  end if

!netCDF distributed write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'dR8', action=MPP_OVERWR, &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time8 = i*10.
     call mpp_write( unit, f, domain, data8, time8 )
  end do
  call mpp_close(unit)

!netCDF single-threaded write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR8', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nxg)/) )

  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nyg)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )

  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )

  do i = 0,nt-1
     time8 = i*10.
     call mpp_write( unit, f, domain, data8, time8)
  end do
  call mpp_close(unit)
  allocate( rdata8(is:ie,js:je,nz) )

!reopen and test appending write
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR8', action=MPP_APPEND, form=MPP_NETCDF, threading=MPP_SINGLE )
  call mpp_write( unit, f, domain, data8, time8)
  call mpp_close( unit )

!netCDF multi-threaded read
  call mpp_sync()
  call mpp_open( unit, trim(type)//"_"//trim(file)//'sR8', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  allocate( atts(natt) )
  allocate( axes(ndim) )
  allocate( vars(nvar) )
  allocate( tstamp8(ntime) )
  call mpp_get_atts ( unit, atts(:) )
  call mpp_get_axes ( unit, axes(:) )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_times( unit, tstamp8(:) )

  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  call mpp_read( unit, vars(1), domain, rdata8, 1 )
  rchk = mpp_chksum(rdata8(is:ie,js:je,:))
  chk  = mpp_chksum( data8(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//' R8: single-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'R8 Checksum error on multi-threaded/single-fileset netCDF read for type ' &
               //trim(type) )
  end if

  call mpp_close(unit)
  deallocate( atts, axes, vars, tstamp8 )

!netCDF distributed read
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(type)//"_"//trim(file)//'dR8', action=MPP_RDONLY,  &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  allocate( atts(natt) )
  allocate( axes(ndim) )
  allocate( vars(nvar) )
  allocate( tstamp8(ntime) )
  call mpp_get_atts ( unit, atts(:) )
  call mpp_get_axes ( unit, axes(:) )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_times( unit, tstamp8(:) )

  call mpp_get_atts(vars(1),name=varname)
  rdata8 = 0

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )

  call mpp_read( unit, vars(1), domain, rdata8, 1 )

  rchk = mpp_chksum(rdata8(is:ie,js:je,:))
  chk  = mpp_chksum( data8(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//' R8: multi-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'R8 Checksum error on multi-threaded/multi-fileset netCDF read for type ' &
           //trim(type) )
  end if

  deallocate( atts, axes, vars, tstamp8 )
  deallocate( rdata8, gdata8, data8)

  end subroutine test_netcdf_io_R8

end program test_io_R4_R8
