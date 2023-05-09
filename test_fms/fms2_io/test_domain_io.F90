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

program test_domain_read
  use   mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_domains, mpp_define_io_domain, &
                               mpp_get_compute_domain, mpp_get_data_domain, domain2d
  use   mpp_mod,         only: mpp_chksum, mpp_pe, mpp_root_pe, mpp_error, FATAL
  use   fms2_io_mod,     only: open_file, register_axis, register_variable_attribute, close_file, &
                               FmsNetcdfDomainFile_t, write_data, register_field, read_data, &
                               parse_mask_table
  use   fms_mod,         only: fms_init, fms_end
  use   platform_mod,    only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none

  type data_type
    real(kind=r4_kind),    allocatable   :: var_r4(:,:,:,:,:)
    real(kind=r8_kind),    allocatable   :: var_r8(:,:,:,:,:)
    integer(kind=i4_kind), allocatable   :: var_i4(:,:,:,:,:)
    integer(kind=i8_kind), allocatable   :: var_i8(:,:,:,:,:)
  end type

  integer, dimension(2)                 :: layout = (/2,3/) !< Domain layout
  integer                               :: ndim1            !< Number of points in dim1
  integer                               :: ndim2            !< Number of points in dim2
  integer                               :: ndim3            !< Number of points in dim3
  integer                               :: ndim4            !< Number of points in dim4
  integer                               :: ndim5            !< Number of points in dim5
  type(domain2d)                        :: Domain           !< Domain with mask table
  type(FmsNetcdfDomainFile_t)           :: fileobj          !< fms2io fileobj for domain decomposed
  character(len=6), dimension(5)        :: names            !< Dimensions names
  type(data_type)                       :: var_data_in      !< Variable data written in
  type(data_type)                       :: var_data_out     !< Variable data read in

  call fms_init

  ndim1 = 96
  ndim2 = 96
  ndim3 = 2
  ndim4 = 2
  ndim5 = 2

  call mpp_domains_set_stack_size(17280000)
  call mpp_define_domains( (/1,ndim1,1,ndim2/), layout, Domain, xhalo=3, yhalo=2)
  call mpp_define_io_domain(Domain, (/1,1/))

  if (open_file(fileobj, "test_domain_write.nc", "overwrite", domain, nc_format="netcdf4")) then
    names(1) = "dim1"
    names(2) = "dim2"
    names(3) = "dim3"
    names(4) = "dim4"
    names(5) = "dim5"

    call register_axis(fileobj, "dim1", "x")
    call register_axis(fileobj, "dim2", "y")
    call register_axis(fileobj, "dim3", ndim3)
    call register_axis(fileobj, "dim4", ndim4)
    call register_axis(fileobj, "dim5", ndim5)

    call register_field_wrapper(fileobj, "var2", names, 2)
    call register_field_wrapper(fileobj, "var3", names, 3)
    call register_field_wrapper(fileobj, "var4", names, 4)
    call register_field_wrapper(fileobj, "var5", names, 5)

    call var_data_alloc(var_data_in, Domain, ndim3, ndim4, ndim5)
    call var_data_init(var_data_in)
    call var_data_set(var_data_in, Domain, ndim3, ndim4, ndim5)

    call write_data_wrapper(fileobj, "r4", var_data_in%var_r4)
    call write_data_wrapper(fileobj, "r8", var_data_in%var_r8)
    call write_data_wrapper(fileobj, "i4", var_data_in%var_i4)
    call write_data_wrapper(fileobj, "i8", var_data_in%var_i8)

    call close_file(fileobj)
  endif

  if (open_file(fileobj, "test_domain_write.nc", "read", domain, nc_format="netcdf4")) then
    call register_axis(fileobj, "dim1", "x")
    call register_axis(fileobj, "dim2", "y")

    call var_data_alloc(var_data_out, Domain, ndim3, ndim4, ndim5)
    call read_data_wrapper(fileobj, "var2", 2, var_data_out, var_data_in)
    call read_data_wrapper(fileobj, "var3", 3, var_data_out, var_data_in)
    call read_data_wrapper(fileobj, "var4", 4, var_data_out, var_data_in)
    call read_data_wrapper(fileobj, "var5", 5, var_data_out, var_data_in)

    call close_file(fileobj)
  endif
  call fms_end

  contains

  subroutine register_field_wrapper(fileobj, var_name, dimension_names, ndim)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj            !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_name           !< Name of the variable
    character(len=*),            intent(in)    :: dimension_names(:) !< dimension names
    integer,                     intent(in)    :: ndim               !< Number of dimension

    call register_field(fileobj, trim(var_name)//"_r8", "double", names(1:ndim))
    call register_field(fileobj, trim(var_name)//"_r4", "float",  names(1:ndim))
    call register_field(fileobj, trim(var_name)//"_i8", "int",    names(1:ndim))
    call register_field(fileobj, trim(var_name)//"_i4", "int64",  names(1:ndim))
  end subroutine register_field_wrapper

  subroutine var_data_alloc(var_data, Domain, dim3, dim4, dim5)
    type(data_type),             intent(inout)    :: var_data         !< Variable data
    type(domain2d),              intent(in)       :: Domain           !< Domain with mask table
    integer,                     intent(in)       :: dim3             !< Size of dim3
    integer,                     intent(in)       :: dim4             !< Size of dim4
    integer,                     intent(in)       :: dim5             !< Size of dim5

    integer :: is !< Starting x index
    integer :: ie !< Ending x index
    integer :: js !< Starting y index
    integer :: je !< Ending y index

    call mpp_get_data_domain(Domain, is, ie, js, je) !< This includes halos (but halos will not be written!)

    allocate(var_data%var_r4(is:ie, js:je, dim3, dim4, dim5))
    allocate(var_data%var_r8(is:ie, js:je, dim3, dim4, dim5))
    allocate(var_data%var_i4(is:ie, js:je, dim3, dim4, dim5))
    allocate(var_data%var_i8(is:ie, js:je, dim3, dim4, dim5))
  end subroutine var_data_alloc

  subroutine var_data_init(var_data)
    type(data_type),             intent(inout)    :: var_data         !< Variable data

    var_data%var_r4 = real(-999.99, kind=r4_kind)
    var_data%var_r8 = real(-999.99, kind=r8_kind)
    var_data%var_i4 = int(-999, kind=i4_kind)
    var_data%var_i8 = int(-999, kind=i8_kind)
  end subroutine var_data_init

  subroutine var_data_set(var_data, domain, dim3, dim4, dim5)
    type(data_type),             intent(inout)    :: var_data         !< Variable data
    type(domain2d),              intent(in)       :: Domain           !< Domain with mask table
    integer,                     intent(in)       :: dim3             !< Size of dim3
    integer,                     intent(in)       :: dim4             !< Size of dim4
    integer,                     intent(in)       :: dim5             !< Size of dim5

    integer :: i, j, k, l, m

    integer :: is !< Starting x index
    integer :: ie !< Ending x index
    integer :: js !< Starting y index
    integer :: je !< Ending y index
    integer :: var_count

    call mpp_get_compute_domain(Domain, is, ie, js, je) !< This does not include halos!

    var_count = 0
    do i = is, ie
      do j = js, je
        do k = 1, dim3
          do l = 1, dim4
            do m = 1, dim5
              var_count = var_count + 1
              var_data%var_r4(i,j,k,l,m) = real(var_count, kind=r4_kind)
              var_data%var_r8(i,j,k,l,m) = real(var_count, kind=r8_kind)
              var_data%var_i4(i,j,k,l,m) = int(var_count, kind=i4_kind)
              var_data%var_i8(i,j,k,l,m) = int(var_count, kind=i8_kind)
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine

  subroutine write_data_wrapper(fileobj, var_kind, var_data)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj             !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_kind            !< The kind of the variable
    class(*),                    intent(in)    :: var_data(:,:,:,:,:) !< Variable data

    call write_data(fileobj, "var2_"//trim(var_kind), var_data(:,:,1,1,1))
    call write_data(fileobj, "var3_"//trim(var_kind), var_data(:,:,:,1,1))
    call write_data(fileobj, "var4_"//trim(var_kind), var_data(:,:,:,:,1))
    call write_data(fileobj, "var5_"//trim(var_kind), var_data(:,:,:,:,:))

  end subroutine

  subroutine read_data_wrapper(fileobj, var_name, dim, var_data, ref_data)
    type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj             !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_name            !< The kind of the variable
    integer,                     intent(in)    :: dim
    type(data_type),             intent(inout) :: var_data
    type(data_type),             intent(inout) :: ref_data

    select case(dim)
    case(2)
      call var_data_init(var_data)
      call read_data(fileobj, trim(var_name)//"_r4", var_data%var_r4(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,1,1,1)), mpp_chksum(ref_data%var_r4(:,:,1,1,1)), "var2_r4")

      call read_data(fileobj, trim(var_name)//"_r8", var_data%var_r8(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,1,1,1)), mpp_chksum(ref_data%var_r8(:,:,1,1,1)), "var2_r8")

      call read_data(fileobj, trim(var_name)//"_i4", var_data%var_i4(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,1,1,1)), mpp_chksum(ref_data%var_i4(:,:,1,1,1)), "var2_i4")

      call read_data(fileobj, trim(var_name)//"_i8", var_data%var_i8(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,1,1,1)), mpp_chksum(ref_data%var_i8(:,:,1,1,1)), "var2_i8")
    case(3)
      call var_data_init(var_data)
      call read_data(fileobj, trim(var_name)//"_r4", var_data%var_r4(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,1,1)), mpp_chksum(ref_data%var_r4(:,:,:,1,1)), "var3_r4")

      call read_data(fileobj, trim(var_name)//"_r8", var_data%var_r8(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,1,1)), mpp_chksum(ref_data%var_r8(:,:,:,1,1)), "var3_r8")

      call read_data(fileobj, trim(var_name)//"_i4", var_data%var_i4(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,1,1)), mpp_chksum(ref_data%var_i4(:,:,:,1,1)), "var3_i4")

      call read_data(fileobj, trim(var_name)//"_i8", var_data%var_i8(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,1,1)), mpp_chksum(ref_data%var_i8(:,:,:,1,1)), "var3_i8")
    case(4)
      call var_data_init(var_data)
      call read_data(fileobj, trim(var_name)//"_r4", var_data%var_r4(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,:,1)), mpp_chksum(ref_data%var_r4(:,:,:,:,1)), "var4_r4")

      call read_data(fileobj, trim(var_name)//"_r8", var_data%var_r8(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,:,1)), mpp_chksum(ref_data%var_r8(:,:,:,:,1)), "var4_r8")

      call read_data(fileobj, trim(var_name)//"_i4", var_data%var_i4(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,1)), mpp_chksum(ref_data%var_i4(:,:,:,:,1)), "var4_i4")

      call read_data(fileobj, trim(var_name)//"_i8", var_data%var_i8(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,1)), mpp_chksum(ref_data%var_i8(:,:,:,:,1)), "var4_i8")
    case(5)
      call var_data_init(var_data)
      call read_data(fileobj, trim(var_name)//"_r4", var_data%var_r4(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,:,:)), mpp_chksum(ref_data%var_r4(:,:,:,:,:)), "var5_r4")

      call read_data(fileobj, trim(var_name)//"_r8", var_data%var_r8(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,:,:)), mpp_chksum(ref_data%var_r8(:,:,:,:,:)), "var5_r8")

      call read_data(fileobj, trim(var_name)//"_i4", var_data%var_i4(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,:)), mpp_chksum(ref_data%var_i4(:,:,:,:,:)), "var5_i4")

      call read_data(fileobj, trim(var_name)//"_i8", var_data%var_i8(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,:)), mpp_chksum(ref_data%var_i8(:,:,:,:,:)), "var5_i8")
    end select

  end subroutine read_data_wrapper

  subroutine compare_var_data(check_sum_in, check_sum_ref, varname)
    integer(kind=i8_kind), intent(in) :: check_sum_in
    integer(kind=i8_kind), intent(in) :: check_sum_ref
    character(len=*), intent(in) :: varname

   if (check_sum_ref .ne. check_sum_in) call mpp_error(FATAL, &
     "Checksums do not match for variable: "//trim(varname))
  end subroutine compare_var_data

end program test_domain_read
