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
  use mpp_domains_mod,  only : mpp_get_domain_components, mpp_define_mosaic, mpp_define_io_domain
  use mosaic_mod,       only : get_mosaic_ntiles
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type
  use horiz_interp_mod, only : horiz_interp_init
  use axis_utils_mod,   only : get_axis_cart
  use fms_io_mod,       only : read_data, write_data
  use fms_io_mod,       only : field_size, fms_io_exit, get_mosaic_tile_file
  use fms_mod,          only : fms_init, fms_end, open_namelist_file, close_file, file_exist, field_exist
  use fms_mod,          only : check_nml_error, write_version_number, lowercase
  use constants_mod,    only : constants_init, PI
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type
  use axis_utils_mod,    only : get_axis_bounds
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
  integer            :: io_layout(2) = (/1,1/)

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
     real, allocatable :: mask_src(:,:)
     real    :: D2R
     integer :: k, i, j
     character(len=128) :: dst_file2

     call get_mosaic_tile_file(dst_file, dst_file2, .false., domain=domain)

     call mpp_get_domain_components( domain, xdom, ydom )
     !--- set up meta data
!     call mpp_open(unit,dst_file,action=MPP_OVERWR,form=MPP_NETCDF,threading=MPP_MULTI, fileset=MPP_MULTI)
     call mpp_open(unit,dst_file2,action=MPP_OVERWR,form=MPP_NETCDF,domain=domain)
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

     if(trim(interp_method) == 'conservative' ) then
        ! get the mask_in
        call mpp_read(src_unit,fields(src_field_index),src_data, tindex=1)
        allocate(mask_src(nx_src,ny_src))
        mask_src = 1.0
        do j = 1, ny_src
           do i = 1, nx_src
              if(src_data(i,j,1) == missing_value) mask_src(i,j) = 0.0
           enddo
        enddo

        write(stdout(),*) "use 2-D version of conservative interpolation"
        call horiz_interp_new(Interp, x_src*D2R, y_src*D2R, x_dst*D2R, y_dst*D2R, &
          interp_method = trim(interp_method), mask_in=mask_src )
        deallocate(mask_src)
     else if(trim(interp_method) == "bilinear" .and. use_2d_version) then
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
    integer :: i, j, siz(4), ntiles, mytile, npes_per_tile
    real, allocatable :: tmp(:,:)

    ntiles = get_mosaic_ntiles(dst_grid)
    if(ntiles .NE. 1 .and. ntiles .NE. 6) call mpp_error(FATAL, "test_horiz_interp: ntiles should be 1 or 6")
    npes_per_tile = mpp_npes()/ntiles
    mytile = mpp_pe()/npes_per_tile + 1

    if (field_exist(dst_grid, "gridfiles" )) then
       call read_data(dst_grid, "gridfiles", tile_file, level=mytile)
       tile_file = 'INPUT/'//trim(tile_file)
    else
       call mpp_error(FATAL, "test_horiz_interp: field gridfiles does not exist in file "//trim(dst_grid) )
    endif

    call field_size(tile_file, "x", siz)
    nx_dst = (siz(1)-1)/2
    ny_dst = (siz(2)-1)/2

    if(layout(1)*layout(2)*ntiles .NE. mpp_npes() ) then
       call mpp_define_layout( (/1,nx_dst,1,ny_dst/), mpp_npes()/ntiles, layout )
    end if

    if(ntiles == 1) then
       call mpp_define_domains( (/1,nx_dst,1,ny_dst/), layout, Domain, name='test_data_override')
    else
       call define_cubic_mosaic(domain, nx_dst, ny_dst, layout)
    endif
    call mpp_define_io_domain(domain, io_layout)

    call mpp_get_compute_domain(Domain, is, ie, js, je)

    allocate(tmp(2*is-1:2*ie+1,2*js-1:2*je+1))

    start = 1; nread = 1
    start(1) = 2*is-1; nread(1) = 2*(ie-is+1)+1
    start(2) = 2*js-1; nread(2) = 2*(je-js+1)+1

    call read_data(tile_file, "x", tmp, start, nread, domain)
    if(trim(interp_method) == 'conservative' ) then
       allocate(x_dst(is:ie+1,js:je+1), y_dst(is:ie+1,js:je+1))
       do j = js, je+1
          do i = is, ie+1
             x_dst(i,j) = tmp(2*i-1,2*j-1)
          enddo
       enddo
    else
       allocate(x_dst(is:ie,js:je), y_dst(is:ie,js:je))
       do j = js, je
          do i = is, ie
             x_dst(i,j) = tmp(2*i,2*j)
          enddo
       enddo
    endif
    call read_data(tile_file, "y", tmp, start, nread, domain)

    if(trim(interp_method) == 'conservative' ) then
       do j = js, je+1
          do i = is, ie+1
             y_dst(i,j) = tmp(2*i-1,2*j-1)
          enddo
       enddo
    else
       do j = js, je
          do i = is, ie
             y_dst(i,j) = tmp(2*i,2*j)
          enddo
       enddo
    endif

    deallocate(tmp)

  end subroutine read_dst_grid

  !--- define mosaic domain for cubic grid
  subroutine define_cubic_mosaic(domain, ni, nj, layout)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: ni, nj
    integer   :: pe_start(6), pe_end(6)
    integer   :: global_indices(4,6), layout2d(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer :: p, npes_per_tile, i

    ntiles = 6
    num_contact = 12
    p = 0
    npes_per_tile = mpp_npes()/ntiles
    do i = 1, 6
       layout2d(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
       pe_start(i) = p
       p = p + npes_per_tile
       pe_end(i) = p-1
    enddo


    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj
    call mpp_define_mosaic(global_indices, layout2d, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., name = 'cubic mosaic'  )

    return

  end subroutine define_cubic_mosaic


  subroutine read_src_file

    integer                                    :: ndim, nvar, natt, n
    integer                                    :: nt, i, j, k, jj, len1
    character(len=1)                           :: cart
    character(len=32)                          :: name
    type(axistype) :: axis_bounds(2)

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
          if(trim(interp_method) == 'conservative' ) then
             allocate(x_src(nx_src+1))
             call get_axis_bounds(axes(j),axis_bounds(1), axes)
             call mpp_get_axis_data(axis_bounds(1), x_src)
          else
             allocate(x_src(nx_src))
             call mpp_get_axis_data(axes(j),x_src)
          endif
       case('Y')
          ny_src = len1
          if(trim(interp_method) == 'conservative' ) then
             allocate(y_src(ny_src+1))
             call get_axis_bounds(axes(j),axis_bounds(2), axes)
             call mpp_get_axis_data(axis_bounds(2), y_src)
          else
             allocate(y_src(ny_src))
             call mpp_get_axis_data(axes(j),y_src)
          endif
       case('Z')
          nk = len1
       end select
    enddo

    if(nx_src==0) call mpp_error(FATAL,'test_horiz_interp: file ' &
         //trim(src_file)//' does not contain axis with cartesian attributes = "X" ')
    if(ny_src==0) call mpp_error(FATAL,'test_horiz_interp: file '&
         //trim(src_file)//' does not contain axis with cartesian attributes = "Y" ')

    if(trim(interp_method) .ne. 'conservative') then
       allocate(x_src_2d(nx_src,ny_src), y_src_2d(nx_src,ny_src))
       do i = 1, nx_src
          x_src_2d(i,:) = x_src(i)
       enddo
       do i = 1, ny_src
          y_src_2d(:,i) = y_src(i)
       enddo
    endif

    !--- get the missing value
    call mpp_get_atts(fields(src_field_index),missing=missing_value)


  end subroutine read_src_file


end program test

#else
module null_test_horiz_interp
end module

#endif
