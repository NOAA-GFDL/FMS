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
!> @author Ryan Mulhall
!> @brief Unit test for mpp_write/mpp_read on mosaics
!> @email gfdl.climate.model.info@noaa.gov
!> @description Performs reads and writes on mosaic files using mpp_write
!> and mpp_read using 32 and 64 bit reals
program test_io_mosaic_R4_R8

  use platform_mod
  use mpp_mod,         only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sync_self
  use mpp_mod,         only : FATAL, NOTE, mpp_chksum, MPP_DEBUG, mpp_set_stack_size, MPP_CLOCK_SYNC
  use mpp_mod,         only : mpp_sync, mpp_exit, mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init_test_full_init
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
  character(len=64) :: file='test_mosaic', iospec='-F cachea'
  integer           :: layout(2) = (/1,1/)
  integer           :: ntiles_x=3, ntiles_y=4  ! total number of tiles will be ntiles_x*ntiles_y,
                                               ! the grid size for each tile will be (nx/ntiles_x, ny/ntiles_y)
                                               ! set ntiles > 1 to test the efficiency of mpp_io.
  integer           :: io_layout(2) = (/1,1/)  ! set io_layout to divide each tile into io_layout(1)*io_layout(2)
                                               ! group and write out data from the root pe of each group.
  integer           :: pack_size = 1

  namelist / test_io_mosaic_nml / nx, ny, nz, nt, halo, stackmax, stackmaxd, debug, file, iospec, &
                               ntiles_x, ntiles_y, layout, io_layout

  integer        :: pe, npes, io_status
  type(domain2D) :: domain
  integer            :: tks_per_sec
  integer            :: i,j,k, unit=7
  logical            :: opened
  character(len=64)  :: varname
  real(r8_kind)      :: time8
  real(r4_kind)      :: time4
  type(axistype)     :: x, y, z, t
  type(fieldtype)    :: f
  type(domain1D)     :: xdom, ydom
  integer(i8_kind)   :: rchk, chk
  real(r8_kind)      :: doubledata = 0.0
  real               :: realarray(4)
  integer            :: ierr
! initialize mpp and read input namelist
  call mpp_init(test_level=mpp_init_test_full_init)
  pe = mpp_pe()
  npes = mpp_npes()
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, file='input.nml', iostat=io_status)
  read( unit,test_io_mosaic_nml, iostat=io_status )
  close(unit)

      if (io_status > 0) then
         call mpp_error(FATAL,'=>test_io_mosaic_R4_R8: Error reading input.nml')
      endif

  call SYSTEM_CLOCK( count_rate=tks_per_sec )
  if( debug )then
      call mpp_io_init(MPP_DEBUG)
  else
      call mpp_io_init()
  end if
  call mpp_set_stack_size(stackmax)
  call mpp_domains_set_stack_size(stackmaxd)

  if( pe.EQ.mpp_root_pe() )then
      print '(a,6i6)', 'npes, nx, ny, nz, nt, halo=', npes, nx, ny, nz, nt, halo
      print *, 'Using NEW domaintypes and calls...'
  end if

  write( file,'(a,i3.3)' )trim(file), npes
! set layouts to test with mosaic tiles
  io_layout(1) = 3
  io_layout(2) = 2
  layout(1) = 3
  layout(2) = 4
! determine the pack_size
  pack_size = size(transfer(doubledata, realarray))
  if( pack_size .NE. 1 .AND. pack_size .NE. 2) call mpp_error(FATAL,'test_io_mosaic_R4_R8: pack_size should be 1 or 2')

  ! test different mosaic reads with r4
  call test_netcdf_io_mosaic_R4('Single_tile_mult_file_R4', layout, 1, 1, (/1,1/) )
  call test_netcdf_io_mosaic_R4('Single_tile_with_group_R4', layout, 1, 1, io_layout)
  call test_netcdf_io_mosaic_R4('Mult_tile_R4', layout, io_layout(1), io_layout(2), (/1,1/))
  call test_netcdf_io_mosaic_R4('Mult_tile_with_group_R4', layout, ntiles_x, ntiles_y, io_layout)
  ! test different mosaic reads with r8
  call test_netcdf_io_mosaic_R8('Single_tile_mult_file_R8', layout, 1, 1, (/1,1/) )
  call test_netcdf_io_mosaic_R8('Single_tile_with_group_R8', layout, 1, 1, io_layout)
  call test_netcdf_io_mosaic_R8('Mult_tile_R8', layout, io_layout(1), io_layout(2), (/1,1/))
  call test_netcdf_io_mosaic_R8('Mult_tile_with_group_R8', layout, ntiles_x, ntiles_y, io_layout)

  call mpp_io_exit()
  call mpp_domains_exit()
  call MPI_FINALIZE(ierr)

  contains

  !------------------------------------------------------------------
  subroutine test_netcdf_io_mosaic_R4(type, layout, ntiles_x, ntiles_y, io_layout)
  character(len=*), intent(in) :: type
  integer,          intent(in) :: layout(:)
  integer,          intent(in) :: io_layout(:)
  integer,          intent(in) :: ntiles_x, ntiles_y

  integer                              :: ndim, nvar, natt, ntime
  integer                              :: isc, iec, jsc, jec, nlon, nlat, n, i, j
  integer                              :: my_tile, ncontacts, npes_per_tile, ntiles
  integer, dimension(:),   allocatable :: tile1, istart1, iend1, jstart1, jend1
  integer, dimension(:),   allocatable :: tile2, istart2, iend2, jstart2, jend2
  integer, dimension(:),   allocatable :: pe_start, pe_end
  integer, dimension(:,:), allocatable :: layout2D, global_indices
  character(len=64)                    :: output_file
  logical                              :: is_root_pe
  real(r4_kind), dimension(:,:,:), allocatable  :: data, rdata
  type(fieldtype), save                :: vars(1)

  npes = mpp_npes()

  ncontacts = 0
  ntiles = ntiles_x*ntiles_y

  npes_per_tile = npes/ntiles
  my_tile       = mpp_pe()/npes_per_tile + 1
  is_root_pe = .false.
  if(mpp_pe() == (my_tile-1)*npes_per_tile ) is_root_pe = .true.

  allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
  !--- for simplify purpose, we assume all the tiles have the same size.
  do n = 1, ntiles
     pe_start(n) = (n-1)*npes_per_tile
     pe_end(n)   = n*npes_per_tile-1
  end do
  if(ntiles>1) then
     nlon = nx/ntiles_x
     nlat = ny/ntiles_y
  else
     nlon = nx
     nlat = ny
  endif

  do n = 1, ntiles
     global_indices(:,n) = (/1,nlon,1,nlat/)
     layout2D(1,n)         = layout(1)/ntiles_x
     layout2D(2,n)         = layout(2)/ntiles_y
  end do

  call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, ncontacts, tile1, tile2, &
                         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, pe_start, pe_end, &
                         name = type)
  call mpp_get_compute_domain( domain, isc,  iec,  jsc,  jec  )
  call mpp_get_domain_components(domain, xdom, ydom)
  allocate( data (isc:iec,jsc:jec,nz) )
  allocate( rdata(isc:iec,jsc:jec,nz) )
  do k = 1,nz
     do j = jsc, jec
        do i = isc, iec
           data(i,j,k)  = k + i*1e-3 + j*1e-6
        enddo
     enddo
  enddo

  ! open with netcdf distribute write if ntiles = 1, otherwise single-thread write
  output_file = type
  select case(type)
  case("Single_tile_single_file_R4")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE )
  case("Single_tile_mult_file_R4")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile_R4")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group_R4")
     call mpp_define_io_domain(domain, io_layout)
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case("Mult_tile_with_group_R4")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)

  case default
     call mpp_error(FATAL, "program test_io_mosaic_R4_R8: invaid value of type="//type)
  end select
  ! write data
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nlon)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nlat)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time4 = i*10.
     call mpp_write( unit, f, domain, data, time4 )
  end do
  call mpp_close(unit)

  call mpp_sync()               !wait for previous write to complete
  ! reopen file and check results
  select case(type)
  case("Single_tile_single_file_R4")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  case("Single_tile_mult_file_R4")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile_R4")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, &
         fileset=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group_R4", "Mult_tile_with_group_R4")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case default
     call mpp_error(FATAL, "program test_io_mosaic_R4_R8: invaid value of type="//type)
  end select

  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  do i = 0,nt-1
     call mpp_read( unit, vars(1), domain, rdata, 1 )
  enddo
  ! compare read and stored data to validate successful write/read
  rchk = mpp_chksum(rdata)
  chk  = mpp_chksum( data)
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on netCDF read for type ' &
               //trim(type) )
  end if

  deallocate( rdata, data)
  call mpp_deallocate_domain(domain)

  end subroutine test_netcdf_io_mosaic_R4

  subroutine test_netcdf_io_mosaic_R8(type, layout, ntiles_x, ntiles_y, io_layout)
  character(len=*), intent(in) :: type
  integer,          intent(in) :: layout(:)
  integer,          intent(in) :: io_layout(:)
  integer,          intent(in) :: ntiles_x, ntiles_y

  integer                              :: ndim, nvar, natt, ntime
  integer                              :: isc, iec, jsc, jec, nlon, nlat, n, i, j
  integer                              :: my_tile, ncontacts, npes_per_tile, ntiles
  integer, dimension(:),   allocatable :: tile1, istart1, iend1, jstart1, jend1
  integer, dimension(:),   allocatable :: tile2, istart2, iend2, jstart2, jend2
  integer, dimension(:),   allocatable :: pe_start, pe_end
  integer, dimension(:,:), allocatable :: layout2D, global_indices
  character(len=64)                    :: output_file
  logical                              :: is_root_pe
  real(r8_kind), dimension(:,:,:), allocatable  :: data, rdata
  type(fieldtype), save                :: vars(1)

  npes = mpp_npes()

  ncontacts = 0
  ntiles = ntiles_x*ntiles_y

  npes_per_tile = npes/ntiles
  my_tile       = mpp_pe()/npes_per_tile + 1
  is_root_pe = .false.
  if(mpp_pe() == (my_tile-1)*npes_per_tile ) is_root_pe = .true.

  allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
  !--- for simplify purpose, we assume all the tiles have the same size.
  do n = 1, ntiles
     pe_start(n) = (n-1)*npes_per_tile
     pe_end(n)   = n*npes_per_tile-1
  end do
  if(ntiles>1) then
     nlon = nx/ntiles_x
     nlat = ny/ntiles_y
  else
     nlon = nx
     nlat = ny
  endif

  do n = 1, ntiles
     global_indices(:,n) = (/1,nlon,1,nlat/)
     layout2D(1,n)         = layout(1)/ntiles_x
     layout2D(2,n)         = layout(2)/ntiles_y
  end do

  call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, ncontacts, tile1, tile2, &
                         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, pe_start, pe_end, &
                         name = type)
  call mpp_get_compute_domain( domain, isc,  iec,  jsc,  jec  )
  call mpp_get_domain_components(domain, xdom, ydom)
  allocate( data (isc:iec,jsc:jec,nz) )
  allocate( rdata(isc:iec,jsc:jec,nz) )
  do k = 1,nz
     do j = jsc, jec
        do i = isc, iec
           data(i,j,k)  = k + i*1e-3 + j*1e-6
        enddo
     enddo
  enddo

  ! open with netcdf distribute write if ntiles = 1, otherwise single-thread write
  output_file = type
  select case(type)
  case("Single_tile_single_file_R8")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE )
  case("Single_tile_mult_file_R8")
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile_R8")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group_R8")
     call mpp_define_io_domain(domain, io_layout)
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case("Mult_tile_with_group_R8")
     write(output_file, '(a,I4.4)') type//'.tile', my_tile
     call mpp_open( unit, output_file, action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)

  case default
     call mpp_error(FATAL, "program test_io_mosaic_R4_R8: invaid value of type="//type)
  end select
  ! write data
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nlon)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,nlat)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=pack_size )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time8 = i*10.
     call mpp_write( unit, f, domain, data, time8 )
  end do
  call mpp_close(unit)

  call mpp_sync()               !wait for previous write to complete
  ! reopen file and check results
  select case(type)
  case("Single_tile_single_file_R8")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
  case("Single_tile_mult_file_R8")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  case("Mult_tile_R8")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, &
         fileset=MPP_SINGLE, is_root_pe=is_root_pe )
  case("Single_tile_with_group_R8", "Mult_tile_with_group_R8")
     call mpp_open( unit, output_file, action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI, domain=domain)
  case default
     call mpp_error(FATAL, "program test_io_mosaic_R4_R8: invaid value of type="//type)
  end select

  call mpp_get_info( unit, ndim, nvar, natt, ntime )
  call mpp_get_fields ( unit, vars(:) )
  call mpp_get_atts(vars(1),name=varname)

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  do i = 0,nt-1
     call mpp_read( unit, vars(1), domain, rdata, 1 )
  enddo
  ! compare read and stored data to validate successful write/read
  rchk = mpp_chksum(rdata)
  chk  = mpp_chksum( data)
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on netCDF read for type ' &
               //trim(type) )
  end if


  deallocate( rdata, data)
  call mpp_deallocate_domain(domain)

  end subroutine test_netcdf_io_mosaic_R8

end program test_io_mosaic_R4_R8
