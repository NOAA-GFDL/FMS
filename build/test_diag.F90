program test_diag

use mpp_mod
use time_manager_mod
use diag_manager_mod!, only:register_diag_field,send_data
use mpp_domains_mod
 USE constants_mod
 USE diag_axis_mod
  USE fms_mod
  USE fms_io_mod
implicit none

integer :: nx = 128
integer :: ny = 128
integer :: nz = 40

real :: x(128)
real :: y(128)

real, dimension (128,128) :: twod
real, dimension (:),allocatable :: oned
real :: zerod

integer :: oned_axes(1)
integer :: twod_axes(2)


integer :: ii,diag_field_id,i,ix,iy

logical :: diag_result, tf

character (len=100) :: err
type diag_info
 character (len=15) :: modd
 character (len=15) :: varname
 character (len=100) :: longname
 integer  :: id
endtype diag_info
type(time_type) :: time
type(diag_info) :: rh
type(diag_info) :: shflx 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer :: unit=7
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened

  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  integer :: wide_halo_x = 0, wide_halo_y = 0
  integer :: nx_cubic = 0, ny_cubic = 0
  logical :: test_performance = .false.
  logical :: test_interface = .true.
  logical :: test_nest_domain = .false.
  logical :: test_edge_update = .false.
  logical :: test_group = .false.
  logical :: test_cubic_grid_redistribute = .false.
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
  logical :: test_get_nbr = .FALSE.
  logical :: test_boundary = .false.
  logical :: test_global_sum = .false.
  integer :: ensemble_size
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  integer :: num_iter = 1
  integer :: num_fields = 4


  !--- namelist variable for nest domain
  integer :: tile_fine   = 1
  integer :: tile_coarse = 1
  integer :: istart_fine = 0, iend_fine = -1, jstart_fine = 0, jend_fine = -1
  integer :: istart_coarse = 0, iend_coarse = -1, jstart_coarse = 0, jend_coarse = -1
  integer :: npes_coarse = 0
  integer :: npes_fine   = 0
  integer :: extra_halo = 0
  logical :: mix_2D_3D = .false.
  logical :: test_subset = .false.
  logical :: test_unstruct = .false.
  integer :: nthreads = 1
  integer ::  j, k
!  integer :: layout(2)
  integer :: id
  integer :: outunit, errunit, io_status
  integer :: get_cpu_affinity, base_cpu, omp_get_num_threads, omp_get_thread_num



    character(len=40) :: type='Cubic-Grid'

    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        ::  l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem, lsg, leg

    integer, allocatable, dimension(:)       :: pe_start, pe_end, npts_tile, grid_index, ntiles_grid
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: x1, x2, g1, g2
    real,    allocatable, dimension(:,:,:)   :: a1, a2, gdata
    real,    allocatable, dimension(:,:)     :: rmask
    real,    allocatable, dimension(:)       :: frac_crit
    logical, allocatable, dimension(:,:,:)   :: lmask
    integer, allocatable, dimension(:)       :: isl, iel, jsl, jel
    logical            :: cubic_grid
    character(len=3)   :: text
    integer            :: nx_save, ny_save, tile
    integer            :: ntotal_land, istart, iend, pos
    integer :: pe, npes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: nlon1, nlat1, nlon2, nlat2
  INTEGER, DIMENSION(2) :: layout = (/0,0/)
  INTEGER :: test_number=1
  INTEGER :: nlon=18, nlat=18, nlev=2
  INTEGER :: io_layout(2) = (/0,0/)
  INTEGER :: nstep = 2
  TYPE(time_type) ::  Time_step, Time_end, Time_start, Run_length
  LOGICAL :: used, test_successful
  CHARACTER(len=256) :: err_msg

  INTEGER :: nyc1, jsw, jew, isw, iew
  INTEGER :: numthreads=1, ny_per_thread, idthread
  INTEGER :: months=0, days=0, dt_step=0
integer :: ierr,nml_unit, log_unit, out_unit 
 INTEGER :: id_dat2
  NAMELIST /test_diag_manager_nml/ layout, test_number, nlon, nlat, nlev, io_layout, numthreads, &
                                   dt_step, months, days

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialize stuff
  CALL fms_init
  log_unit = stdlog()
  out_unit = stdout()
  CALL constants_init
  CALL set_calendar_type(JULIAN)
#ifdef INTERNAL_FILE_NML
  READ (input_nml_file, NML=test_diag_manager_nml, IOSTAT=ierr)
#else
  IF ( file_exist('input.nml') ) THEN
     nml_unit = open_namelist_file()
     READ(nml_unit, nml=test_diag_manager_nml, iostat=ierr)
     CALL close_file(nml_unit)
  ELSE
     ! Set ierr to an arbitrary positive number if input.nml does not exist.
     ierr = 100
  END IF
#endif
  WRITE (out_unit,test_diag_manager_nml)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET UP THE UNSTRUCTURED GRID
call mpp_init()

  pe = mpp_pe()
  npes = mpp_npes()

    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
           goto 450
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
           goto 450
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from test_unstruct_update ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_unstruct_update: npes should be multiple of ntiles No test is done for '//trim(type))
           goto 450
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
       allocate(frac_crit(ntiles))
       frac_crit(1) = 0.3; frac_crit(2) = 0.1; frac_crit(3) = 0.6
       frac_crit(4) = 0.2; frac_crit(5) = 0.4; frac_crit(6) = 0.5

    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    !--- define domain
    if( cubic_grid ) then
       call define_cubic_mosaic(type, SG_domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( SG_domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( SG_domain, isd, ied, jsd, jed )

    allocate(lmask(nx,ny,ntiles))
    allocate(npts_tile(ntiles))
    lmask = .false.
    if(mpp_pe() == mpp_root_pe() ) then
       allocate(rmask(nx,ny))
       !--- construct gmask.
       do n = 1, ntiles
          call random_number(rmask)
          do j = 1, ny
             do i = 1, nx
                if(rmask(i,j) > frac_crit(n)) then
                   lmask(i,j,n) = .true.
                endif
             enddo
          enddo
          npts_tile(n) = count(lmask(:,:,n))
       enddo
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
       l = 0
       allocate(isl(0:mpp_npes()-1), iel(0:mpp_npes()-1))
       allocate(jsl(0:mpp_npes()-1), jel(0:mpp_npes()-1))
       call mpp_get_compute_domains(SG_domain,xbegin=isl,xend=iel,ybegin=jsl,yend=jel)

       do n = 1, ntiles
          do j = 1, ny
             do i = 1, nx
                if(lmask(i,j,n)) then
                   l = l + 1
                   grid_index(l) = (j-1)*nx+i
                endif
             enddo
          enddo
       enddo
       deallocate(rmask, isl, iel, jsl, jel)
    endif
    call mpp_broadcast(npts_tile, ntiles, mpp_root_pe())
    if(mpp_pe() .NE. mpp_root_pe()) then
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
    endif
    call mpp_broadcast(grid_index, ntotal_land, mpp_root_pe())

    allocate(ntiles_grid(ntotal_land))
    ntiles_grid = 1
   !--- define the unstructured grid domain
    call mpp_define_unstruct_domain(UG_domain, SG_domain, npts_tile, ntiles_grid, mpp_npes(), 1, grid_index, name="LAND unstruct")
    call mpp_get_UG_compute_domain(UG_domain, istart, iend)

    !--- figure out lmask according to grid_index
    pos = 0
    do n = 1, ntiles
       do l = 1, npts_tile(n)
          pos = pos + 1
          j = (grid_index(pos)-1)/nx + 1
          i = mod((grid_index(pos)-1),nx) + 1
          lmask(i,j,n) = .true.
       enddo
    enddo

    !--- set up data
    allocate(gdata(nx,ny,ntiles))
    gdata = -999
    do n = 1, ntiles
       do j = 1, ny
          do i = 1, nx
             if(lmask(i,j,n)) then
                gdata(i,j,n) = n*1.e+3 + i + j*1.e-3
             endif
          end do
       end do
    end do

    !--- test the 2-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,1), a2(isc:iec,jsc:jec,1 ) )

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
          a1(i,j,1) = gdata(i,j,tile)
       enddo
    enddo
    a2 = -999
    write(mpp_pe()+1000,*) "npts_tile = "
    write(mpp_pe()+1000,*) npts_tile
    write(mpp_pe()+1000,*) "a1 = ", isc, iec, jsc, jec
    do j = jsc, jec
       write(mpp_pe()+1000,*) a1(:,j,1)
    enddo

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       x2(l,1) = gdata(i,j,tile)
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))

    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))


    deallocate(a1,a2,x1,x2)

    !--- test the 3-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,nz), a2(isc:iec,jsc:jec,nz ) )

    tile = mpp_pe()/npes_per_tile + 1
    do k = 1, nz
       do j = jsc, jec
          do i = isc, iec
             a1(i,j,k) = gdata(i,j,tile)
             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
          enddo
       enddo
    enddo
    a2 = -999

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       do k = 1, nz
          x2(l,k) = gdata(i,j,tile) + k*1.e-6
       enddo
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1, x1)

    write(mpp_pe()+1000,*) "x1 = ", istart, iend
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)


    deallocate(a1,a2,x1,x2)

    !--- test the 2-D data is on data domain
    allocate( a1(isd:ied, jsd:jed,1), a2(isd:ied,jsd:jed,1 ) )
    a1 = -999; a2 = -999

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
          a1(i,j,1) = gdata(i,j,tile)
       enddo
    enddo
    a2 = -999
    write(mpp_pe()+1000,*) "npts_tile = "
    write(mpp_pe()+1000,*) npts_tile

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       x2(l,1) = gdata(i,j,tile)
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))

    write(mpp_pe()+1000,*) "x1 = ", istart, iend
    write(mpp_pe()+1000,*) x1
    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))


    deallocate(a1,a2,x1,x2)

    !--- test the 3-D data is on computing domain
    allocate( a1(isd:ied, jsd:jed,nz), a2(isd:ied,jsd:jed,nz ) )
    a1 = -999; a2 = -999

    tile = mpp_pe()/npes_per_tile + 1
    do k = 1, nz
       do j = jsc, jec
          do i = isc, iec
             a1(i,j,k) = gdata(i,j,tile)
             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
          enddo
       enddo
    enddo
    a2 = -999
    write(mpp_pe()+1000,*) "npts_tile = "
    write(mpp_pe()+1000,*) npts_tile
    do j = jsc, jec
       write(mpp_pe()+1000,*) a1(:,j,1)
    enddo

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       do k = 1, nz
          x2(l,k) = gdata(i,j,tile) + k*1.e-6
       enddo
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1, x1)

    write(mpp_pe()+1000,*) "x1 = ", istart, iend
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)


    deallocate(a1,a2,x1,x2)









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

450 write(6,*) "Looks like we made it.  Look how far we've come, Baby"
zerod = 123.4
!oned = 111.1
twod = 222.2
do i=1,nx
x(i) = real(ii)
y(i) = real(ii)*2.0
enddo

!type diag_info
! character (len=15) :: modd
! character (len=15) :: varname
! character (len=100) :: longname
! integer  :: id
!endtype diag_info
 rh%modd="flux"
 rh%varname="rh_ref"
 rh%longname="Relativve Humidity"

 shflx%modd = "flux"
 shflx%varname = "shflx"
 shflx%longname = "Suface heat flux"


write (6,*) "Beginning test."
write(6,*) "MPP_INIT CALLED"
 call mpp_init()

if ( mpp_pe() == mpp_root_pe() ) write (6,*) "init diag"
 time =  set_time(0, 0, 1, err)
write (6,*) err
 if ( mpp_pe() == mpp_root_pe() ) write (6,*) tf
 call diag_manager_init()
!  if ( mpp_pe() == mpp_root_pe() ) 
!write (6,*) "init x-axis"
  !     INTEGER FUNCTION diag_axis_init(name, data, units, cart_name, long_name,
  !           direction, set_name, edges, Domain, Domain2, aux, tile_count)
! ix = diag_axis_init("lon",x,"num","X","x-direction")
!
! iy = diag_axis_init("lat",y,"num","Y","y-direction")
!
 iy = diag_axis_init("uns",y,"num","U","unstrctured grid",domainU=UG_domain)
! oned_axes(1)=ix
! twod_axes(1)=ix
! twod_axes(1)=iy


! if ( mpp_pe() == mpp_root_pe() ) 
!write (6,*) "ix = ", ix
write (6,*) "iy = ", iy

  id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/iy/), Time, 'sample data', 'K')
  ix      = register_diag_field('test_diag_manager_mod', 'dat2_rms', (/iy/), Time, 'sample data', 'K')
allocate ( oned(size(y,dim=1)) )
oned=111.1
 time =  set_time(0, 0, 1, err)
  used = send_data(id_dat2, oned, Time, err_msg=err_msg)
  used = send_data(ix, oned, Time, err_msg=err_msg)
 time =  set_time(1, 0, 1, err)
  used = send_data(id_dat2, oned, Time, err_msg=err_msg)
  used = send_data(ix, oned, Time, err_msg=err_msg)

 write (6,*) "SEND_DATA_ER",err_msg,used

!  rh%id = register_diag_field ("test_diag_manager_mod", "dat2") !> Returns the field index to be used in subsequent calls to send_data
                                              !! Is used as diag_field_id
!  shflx%id = register_diag_field (trim(shflx%modd) , trim(shflx%varname), oned_axes,time)

!! For our example here, rh_ref will be a scalar
! diag_result = send_data(diag_field_id,zerod,time)
!  diag_result = send_data(diag_field_id,zerod)

 if ( mpp_pe() == mpp_root_pe() ) write (6,*) "The result of register (field index) is ",id_dat2
! ii =  get_diag_field_id("flux", "shflx")
!  if ( mpp_pe() == mpp_root_pe() ) write (6,*) "The result of get_diag_field_id is ",ii


! write(6,*) send_data(rh%id,zerod)
    call diag_manager_end (Time)
    call fms_end()
write (6,*) "End test."
end program test_diag
