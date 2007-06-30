#ifdef test_mpp_io
program test
#include <fms_platform.h>

  use mpp_mod,         only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error
  use mpp_mod,         only : FATAL, NOTE, mpp_chksum, MPP_DEBUG, mpp_set_stack_size
  use mpp_mod,         only : mpp_sync, mpp_exit
  use mpp_domains_mod, only : mpp_define_domains, mpp_domains_set_stack_size, domain1D
  use mpp_domains_mod, only : domain2D, mpp_define_layout, mpp_get_domain_components
  use mpp_domains_mod, only : mpp_get_memory_domain, mpp_get_compute_domain, mpp_domains_exit
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER, mpp_get_data_domain
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
  namelist / test_mpp_io_nml / nx, ny, nz, nt, halo, stackmax, stackmaxd, debug, file, iospec

  integer        :: pe, npes
  type(domain2D) :: domain

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

  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='input.nml', err=10 )
  read( unit,test_mpp_io_nml )
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

  write( file,'(a,i3.3)' )trim(file), npes

  call test_netcdf_io('Simple')
  call test_netcdf_io('Symmetry_T_cell')
  call test_netcdf_io('Symmetry_E_cell')
  call test_netcdf_io('Symmetry_N_cell')
  call test_netcdf_io('Symmetry_C_cell')
  call test_netcdf_io('Symmetry_T_cell_memory')
  call test_netcdf_io('Symmetry_E_cell_memory')
  call test_netcdf_io('Symmetry_N_cell_memory')
  call test_netcdf_io('Symmetry_C_cell_memory')
  call mpp_io_exit()
  call mpp_domains_exit()
  call mpp_exit()

  contains

  !------------------------------------------------------------------

  subroutine test_netcdf_io(type)
  character(len=*), intent(in) :: type

  integer :: is, ie, js, je, isd, ied, jsd, jed, ism, iem, jsm, jem
  integer :: ishift, jshift, msize(2), ioff, joff
  logical :: symmetry

  !--- determine the shift and symmetry according to type, 
  select case(type)
  case('Simple')
     ishift = 0; jshift = 0; symmetry = .false.
  case('Symmetry_T_cell', 'Symmetry_T_cell_memory')
     ishift = 0; jshift = 0; symmetry = .true.
  case('Symmetry_E_cell', 'Symmetry_E_cell_memory')
     ishift = 1; jshift = 0; symmetry = .true.
  case('Symmetry_N_cell', 'Symmetry_N_cell_memory')
     ishift = 0; jshift = 1; symmetry = .true.
  case('Symmetry_C_cell', 'Symmetry_C_cell_memory')
     ishift = 1; jshift = 1; symmetry = .true.
  case default
     call mpp_error(FATAL, "type = "//type//" is not a valid test type")
  end select

!define global data array
  allocate( gdata(nx+ishift,ny+jshift,nz) )
  gdata = 0.
  do k = 1,nz
     do j = 1,ny+jshift
        do i = 1,nx+ishift
           gdata(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do

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

  call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
  call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
  call mpp_get_memory_domain ( domain, ism, iem, jsm, jem )
  call mpp_get_domain_components( domain, xdom, ydom )
  ioff = ism - isd; joff = jsm - jsd
  ie  = ie +ishift; je  = je +jshift
  ied = ied+ishift; jed = jed+jshift  
  iem = iem+ishift; jem = jem+jshift 
  allocate( data(ism:iem,jsm:jem,nz) )
  data(is+ioff:ie+ioff,js+joff:je+joff,:) = gdata(is:ie,js:je,:)

!tests

!sequential write: single-threaded formatted: only if small
  if( nx*ny*nz*nt.LT.1000 .AND. index(type,"memory") .NE. 0 )then
      if( pe.EQ.mpp_root_pe() )print *, 'sequential write: single-threaded formatted'
!here the only test is a successful write: please look at test.txt for verification.
      call mpp_open( unit, trim(file)//'s.txt', action=MPP_OVERWR, form=MPP_ASCII, threading=MPP_SINGLE )
      call mpp_write_meta( unit, x, 'X', 'km', 'X distance', domain=xdom, data=(/(i-1.,i=1,nx+ishift)/) )
      call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', domain=ydom, data=(/(i-1.,i=1,ny+jshift)/) )
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
  call mpp_open( unit, trim(type)//"_"//trim(file)//'d', action=MPP_OVERWR, &
                 form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_MULTI )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nx+ishift)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,ny+jshift)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z', data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )
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
  call mpp_open( unit, trim(type)//"_"//trim(file)//'s', action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )
  call mpp_write_meta( unit, x, 'X', 'km', 'X distance', 'X', domain=xdom, data=(/(i-1.,i=1,nx+ishift)/) )
  call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', 'Y', domain=ydom, data=(/(i-1.,i=1,ny+jshift)/) )
  call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', 'Z',              data=(/(i-1.,i=1,nz)/) )
  call mpp_write_meta( unit, t, 'T', 'sec', 'Time', 'T' )
  call mpp_write_meta( unit, f, (/x,y,z,t/), 'Data', 'metres', 'Random data', pack=1 )
  call mpp_write( unit, x )
  call mpp_write( unit, y )
  call mpp_write( unit, z )
  do i = 0,nt-1
     time = i*10.
     call mpp_write( unit, f, domain, data, time)
  end do
  call mpp_close(unit)

  allocate( rdata(is:ie,js:je,nz) )
!netCDF multi-threaded read
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF multi-threaded read'
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(type)//"_"//trim(file)//'s', action=MPP_RDONLY,  &
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
  call mpp_read( unit, vars(1), domain, rdata, 1 )
  rchk = mpp_chksum(rdata(is:ie,js:je,:))
  chk  = mpp_chksum( data(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': single-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on multi-threaded/single-fileset netCDF read for type ' &
               //trim(type) )
  end if

  deallocate( atts, axes, vars, tstamp )

!netCDF distributed read
  if( pe.EQ.mpp_root_pe() )print *, 'netCDF multi-threaded read'
  call mpp_sync()               !wait for previous write to complete
  call mpp_open( unit, trim(type)//"_"//trim(file)//'d', action=MPP_RDONLY,  &
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
  rdata = 0

  if( varname.NE.'Data' )call mpp_error( FATAL, 'File being read is not the expected one.' )
  call mpp_read( unit, vars(1), domain, rdata, 1 )
  rchk = mpp_chksum(rdata(is:ie,js:je,:))
  chk  = mpp_chksum( data(is+ioff:ie+ioff,js+joff:je+joff,:))
  if( pe.EQ.mpp_root_pe() )print '(a,2z18)', trim(type)//' checksum=', rchk, chk
  if( rchk == chk ) then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(type)//': multi-fileset: data comparison are OK.' )
  else
      call mpp_error( FATAL, 'Checksum error on multi-threaded/multi-fileset netCDF read for type ' &
           //trim(type) )
  end if

  deallocate( atts, axes, vars, tstamp )

  deallocate( rdata, gdata, data)

  end subroutine test_netcdf_io

end program test

#endif
