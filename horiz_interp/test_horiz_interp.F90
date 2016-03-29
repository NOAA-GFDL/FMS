#ifdef TEST_HORIZ_INTERP
!z1l: currently only test bilinear interpolation.

program test

  use mpp_mod,          only : mpp_error, mpp_pe,  mpp_npes, mpp_root_pe
  use mpp_mod,          only : FATAL, stdout, stdlog, mpp_chksum
  use mpp_io_mod,       only : mpp_open, mpp_close, mpp_read
  use mpp_io_mod,       only : axistype, fieldtype
  use mpp_io_mod,       only : mpp_get_info, mpp_get_fields, mpp_get_times
  use mpp_io_mod,       only : mpp_get_axes, mpp_get_axis_data, mpp_get_atts
  use mpp_io_mod,       only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
  use mpp_io_mod,       only : mpp_get_att_name, mpp_get_att_char, mpp_get_att_type, mpp_get_att_real
  use mpp_io_mod,       only : mpp_write_meta, axistype, fieldtype, mpp_write, mpp_close 
  use mpp_domains_mod,  only : mpp_update_domains, mpp_define_domains, domain1d
  use mpp_domains_mod,  only : domain2d, mpp_define_layout, mpp_get_compute_domain
  use mpp_domains_mod,  only : mpp_get_domain_components
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type
  use horiz_interp_mod, only : horiz_interp_init
  use axis_utils_mod,   only : get_axis_cart
  use fms_io_mod,       only : read_data, write_data
  use fms_io_mod,       only : field_size, fms_io_exit
  use fms_mod,          only : fms_init, fms_end, open_namelist_file, close_file, file_exist, field_exist
  use fms_mod,          only : check_nml_error, write_version_number, lowercase
  use constants_mod,    only : constants_init, PI
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type

implicit none

  character(len=256) :: src_file = ""
  character(len=256) :: dst_grid = "INPUT/grid_spec.nc"
  character(len=256) :: field_name = ""
  character(len=256) :: dst_file = "output.nc"
  logical            :: new_missing_handle = .false.
  character(len=256) :: interp_method = "bilinear"
  logical            :: use_2d_version = .false.
  logical            :: write_remap_index = .false.
  integer            :: layout(2) = (/1,1/)

  real, allocatable, dimension(:)            :: x_src, y_src
  real, allocatable, dimension(:,:)          :: x_src_2d, y_src_2d
  real, allocatable, dimension(:,:)          :: x_dst, y_dst
  type(axistype), allocatable, dimension(:)  :: axes
  type(fieldtype), allocatable, dimension(:) :: fields
  type(domain2d)                             :: Domain
  integer :: unit, ierr, io, src_unit, src_field_index, nk
  integer :: nx_src, ny_src, nx_dst, ny_dst, is, ie, js, je
  real    :: missing_value
  integer :: n, ntimes

  namelist /test_horiz_interp_nml/ src_file, field_name, dst_file, dst_grid, new_missing_handle, &
           interp_method, use_2d_version, layout, write_remap_index

  call fms_init
  call horiz_interp_init
  call constants_init

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, test_horiz_interp_nml, iostat=io)
  ierr = check_nml_error(io, 'test_horiz_interp_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr=1
     do while (ierr /= 0)
        read(unit, nml=test_horiz_interp_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'test_horiz_interp_nml')
     enddo
10   call close_file (unit)
  endif
#endif


  if( .not. file_exist(src_file) ) call mpp_error(FATAL, &
       "test_horiz_interp: src_file = "//trim(src_file)//" does not exist")
  if( .not. field_exist(src_file, field_name) ) call mpp_error(FATAL, &
       "test_horiz_interp: field_name = "//trim(field_name)//" does not exist in file "//trim(src_file) )


  ! reading the grid information and missing value from src_file
  call mpp_open(src_unit, trim(src_file), &
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
  call read_src_file()

  ! reading the dst_grid file
  call read_dst_grid()

  !--- currently only test for first time level. The following will read the input data,
  !--- do the remapping and write out data
  call process_data()
  
  call mpp_close(src_unit)

  call fms_io_exit

  call fms_end()

contains


  ! read the input data, do the remapping and write out data
  subroutine process_data()
     real, allocatable, dimension(:,:,:) :: src_data, dst_data
     type(horiz_interp_type)             :: Interp
     type(axistype)     :: xaxis, yaxis, zaxis, taxis
     type(domain1D)     :: xdom, ydom
     type(fieldtype)    :: field, field_istart,field_iend,field_jstart,field_jend
     real    :: D2R
     integer :: k, i

     call mpp_get_domain_components( domain, xdom, ydom )
     !--- set up meta data
!     call mpp_open(unit,dst_file,action=MPP_OVERWR,form=MPP_NETCDF,threading=MPP_MULTI, fileset=MPP_MULTI)
     call mpp_open(unit,dst_file,action=MPP_OVERWR,form=MPP_NETCDF,threading=MPP_MULTI, fileset=MPP_MULTI)
     call mpp_write_meta( unit, xaxis, 'lon', 'none', 'X index', 'X', domain=xdom, data=(/(i-1.,i=1,nx_dst)/) )
     call mpp_write_meta( unit, yaxis, 'lat', 'none', 'Y index', 'Y', domain=ydom, data=(/(i-1.,i=1,ny_dst)/) )
     call mpp_write_meta( unit, zaxis, 'level', 'none', 'Z index', 'Z', data=(/(i-1.,i=1,nk)/) )
     call mpp_write_meta( unit, taxis, 'time', 'none', 'T index', 'T' )
     call mpp_write_meta( unit, field, (/xaxis, yaxis, zaxis, taxis/), field_name, 'none', 'none', missing=missing_value)
     if(write_remap_index) then
        call mpp_write_meta( unit, field_istart, (/xaxis, yaxis/), "istart", 'none', 'none')
        call mpp_write_meta( unit, field_iend, (/xaxis, yaxis/), "iend", 'none', 'none')     
        call mpp_write_meta( unit, field_jstart, (/xaxis, yaxis/), "jstart", 'none', 'none')
        call mpp_write_meta( unit, field_jend, (/xaxis, yaxis/), "jend", 'none', 'none')    
     endif

     call mpp_write( unit, xaxis )
     call mpp_write( unit, yaxis )
     call mpp_write( unit, zaxis )


     D2R = PI/180.

     allocate(src_data(nx_src,ny_src,nk))
     allocate(dst_data(is:ie,js:je,nk))

     if(trim(interp_method) == "bilinear" .and. use_2d_version) then
        write(stdout(),*) "use 2-D version of bilinear interpolation"
        call horiz_interp_new(Interp, x_src_2d*D2R, y_src_2d*D2R, x_dst*D2R, y_dst*D2R, &
          interp_method = trim(interp_method) )    
     else
        write(stdout(),*) "use 1-D version of interpolation"
        call horiz_interp_new(Interp, x_src*D2R, y_src*D2R, x_dst*D2R, y_dst*D2R, &
          interp_method = trim(interp_method), grid_at_center = .true. )    
     endif

     if(write_remap_index) then
        dst_data(:,:,1) = interp%i_lon(:,:,1)
        call mpp_write(unit, field_istart, domain, dst_data(:,:,1))
        dst_data(:,:,1) = interp%i_lon(:,:,2)
        call mpp_write(unit, field_iend, domain, dst_data(:,:,1))
        dst_data(:,:,1) = interp%j_lat(:,:,1)
        call mpp_write(unit, field_jstart, domain, dst_data(:,:,1))
        dst_data(:,:,1) = interp%j_lat(:,:,2)
        call mpp_write(unit, field_jend, domain, dst_data(:,:,1))
     endif

     dst_data = 0
     do n = 1, ntimes
        call mpp_read(src_unit,fields(src_field_index),src_data, tindex=n)

        do k = 1, nk
           call horiz_interp(interp, src_data(:,:,k), dst_data(:,:,k), &
                missing_value=missing_value, new_missing_handle=new_missing_handle)
        enddo
        call mpp_write(unit, field, domain, dst_data, n*1.0)
        write(stdout(),*) "chksum at time = ", n, " = ", mpp_chksum(dst_data)
     enddo

     call mpp_close(unit)
     deallocate(src_data, dst_data)
  
  end subroutine process_data


  subroutine read_dst_grid
    integer :: start(4), nread(4)
    character(len=256) :: tile_file
    integer :: i, j, siz(4)
    real, allocatable :: tmp(:,:)

    if (field_exist(dst_grid, "gridfiles" )) then
       call read_data(dst_grid, "gridfiles", tile_file)
       tile_file = 'INPUT/'//trim(tile_file)
    else
       call mpp_error(FATAL, "test_horiz_interp: field gridfiles does not exist in file "//trim(dst_grid) )
    endif

    call field_size(tile_file, "x", siz)
    nx_dst = (siz(1)-1)/2
    ny_dst = (siz(2)-1)/2

    if(layout(1)*layout(2) .NE. mpp_npes() ) then
       call mpp_define_layout( (/1,nx_dst,1,ny_dst/), mpp_npes(), layout )
    end if

    call mpp_define_domains( (/1,nx_dst,1,ny_dst/), layout, Domain, name='test_data_override')

    call mpp_get_compute_domain(Domain, is, ie, js, je)

    allocate(tmp(2*is-1:2*ie+1,2*js-1:2*je+1))
    allocate(x_dst(is:ie,js:je), y_dst(is:ie,js:je))
 
    start = 1; nread = 1
    start(1) = 2*is-1; nread(1) = 2*(ie-is+1)+1
    start(2) = 2*js-1; nread(2) = 2*(je-js+1)+1

    call read_data(tile_file, "x", tmp, start, nread, domain) 
    do j = js, je
       do i = is, ie
          x_dst(i,j) = tmp(2*i,2*j)
       enddo
    enddo
    call read_data(tile_file, "y", tmp, start, nread, domain)   
    do j = js, je
       do i = is, ie
          y_dst(i,j) = tmp(2*i,2*j)
       enddo
    enddo

    deallocate(tmp)

  end subroutine read_dst_grid

  subroutine read_src_file

    integer                                    :: ndim, nvar, natt, n
    integer                                    :: nt, i, j, k, jj, len1
    character(len=1)                           :: cart
    character(len=32)                          :: name

    call mpp_get_info(src_unit, ndim, nvar, natt, ntimes)
    
    allocate(fields(nvar))
    call mpp_get_fields(src_unit, fields)
    src_field_index = 0
    do i=1,nvar
       call mpp_get_atts(fields(i),name=name)
       if (lowercase(trim(field_name)) == lowercase(trim(name))) then
             src_field_index = i
       endif
    enddo
    if(src_field_index == 0) call mpp_error(FATAL, 'test_horiz_interp: field '&
            //trim(field_name)//' is not in the file '//trim(src_file) )
    !--- get the src grid
    call mpp_get_atts(fields(src_field_index),ndim=ndim)
    allocate(axes(ndim))
    call mpp_get_atts(fields(src_field_index),axes=axes)
    nx_src=0; ny_src=0; nk=1
    do j=1,ndim
       call mpp_get_atts(axes(j),len=len1)
       call get_axis_cart(axes(j),cart)
       select case (cart)
       case ('X')
          nx_src = len1
          allocate(x_src(nx_src))
          call mpp_get_axis_data(axes(j),x_src)
       case('Y')
          ny_src = len1
          allocate(y_src(ny_src))
          call mpp_get_axis_data(axes(j),y_src)
       case('Z')
          nk = len1
       end select
    enddo

    if(nx_src==0) call mpp_error(FATAL,'test_horiz_interp: file ' &
         //trim(src_file)//' does not contain axis with cartesian attributes = "X" ')
    if(ny_src==0) call mpp_error(FATAL,'test_horiz_interp: file '&
         //trim(src_file)//' does not contain axis with cartesian attributes = "Y" ')

    allocate(x_src_2d(nx_src,ny_src), y_src_2d(nx_src,ny_src))
    do i = 1, nx_src
       x_src_2d(i,:) = x_src(i)
    enddo
    do i = 1, ny_src
       y_src_2d(:,i) = y_src(i)
    enddo


    !--- get the missing value
    call mpp_get_atts(fields(src_field_index),missing=missing_value)

  end subroutine read_src_file


end program test

#else
module null_test_horiz_interp
end module  

#endif 
