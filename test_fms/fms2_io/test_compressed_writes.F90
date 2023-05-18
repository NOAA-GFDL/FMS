program test_compressed_writes
  use fms_mod, only: fms_init, fms_end
  use fms2_io_mod
  use mpp_mod
  use   platform_mod,    only: r4_kind, r8_kind, i4_kind, i8_kind

  implicit none

  type data_type
    real(kind=r4_kind),    allocatable   :: var_r4(:,:,:,:,:)
    real(kind=r8_kind),    allocatable   :: var_r8(:,:,:,:,:)
    integer(kind=i4_kind), allocatable   :: var_i4(:,:,:,:,:)
    integer(kind=i8_kind), allocatable   :: var_i8(:,:,:,:,:)
  end type

  type(FmsNetcdfFile_t)                 :: fileobj             !< fms2io fileobj for domain decomposed
  character(len=6), dimension(5)        :: names               !< Dimensions names
  type(data_type)                       :: var_data_in         !< Variable data written in
  type(data_type)                       :: var_data_out         !< Variable data read
  type(data_type)                       :: var_data_ref         !< Variable data read
  integer :: ndim2 = 2
  integer :: ndim3 = 3
  integer :: ndim4 = 4
  integer :: ndim5 = 1
  integer, allocatable :: pes(:)

  call fms_init

  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  if (open_file(fileobj, "test_compressed_writes.nc", "overwrite", nc_format="netcdf4", pelist=pes)) then

    names(1) = "c_xy"
    names(2) = "dim2"
    names(3) = "dim3"
    names(4) = "dim4"
    names(5) = "dim5"

    call register_axis(fileobj, "c_xy", mpp_pe()+1, is_compressed=.true.)
    call register_axis(fileobj, "dim2", ndim2)
    call register_axis(fileobj, "dim3", ndim3)
    call register_axis(fileobj, "dim4", ndim4)
    call register_axis(fileobj, "dim5", ndim5)

    call register_field_wrapper(fileobj, "var1", names, 1)
    call register_field_wrapper(fileobj, "var2", names, 2)
    call register_field_wrapper(fileobj, "var3", names, 3)
    call register_field_wrapper(fileobj, "var4", names, 4)
    call register_field_wrapper(fileobj, "var5", names, 5)

    call var_data_alloc(var_data_in, ndim2, ndim3, ndim4, ndim5, mpp_pe()+1)
    call var_data_set(var_data_in, mpp_pe())

    call write_data_wrapper(fileobj, "r4", var_data_in%var_r4)
    call write_data_wrapper(fileobj, "r8", var_data_in%var_r8)
    call write_data_wrapper(fileobj, "i4", var_data_in%var_i4)
    call write_data_wrapper(fileobj, "i8", var_data_in%var_i8)

    call close_file(fileobj)
  endif

  !< Now check answers
  if (mpp_pe() .eq. mpp_root_pe()) then
    !call mpp_set_current_pelist((/mpp_pe()/))
    if (open_file(fileobj, "test_compressed_writes.nc", "read", nc_format="netcdf4")) then
      call var_data_alloc(var_data_out, ndim2, ndim3, ndim4, ndim5, sum(pes)+mpp_npes())
      call var_data_alloc(var_data_ref, ndim2, ndim3, ndim4, ndim5, sum(pes)+mpp_npes())
      call var_data_set_ref(var_data_ref)

      call read_data_wrapper(fileobj, "var2", 1, var_data_out, var_data_ref)
      call read_data_wrapper(fileobj, "var2", 2, var_data_out, var_data_ref)
      call read_data_wrapper(fileobj, "var3", 3, var_data_out, var_data_ref)
      call read_data_wrapper(fileobj, "var4", 4, var_data_out, var_data_ref)
      call read_data_wrapper(fileobj, "var5", 5, var_data_out, var_data_ref)

      call close_file(fileobj)
    endif
  endif
  call fms_end

  contains

  !> @brief registers all of the possible variable types for a given
  !! number of dimensions
  subroutine register_field_wrapper(fileob, var_name, dimension_names, ndim)
    type(FmsNetcdfFile_t),       intent(inout) :: fileob             !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_name           !< Name of the variable
    character(len=*),            intent(in)    :: dimension_names(:) !< dimension names
    integer,                     intent(in)    :: ndim               !< Number of dimension

    call register_field(fileob, trim(var_name)//"_r8", "double", names(1:ndim))
    call register_field(fileob, trim(var_name)//"_r4", "float",  names(1:ndim))
    call register_field(fileob, trim(var_name)//"_i8", "int64",  names(1:ndim))
    call register_field(fileob, trim(var_name)//"_i4", "int",    names(1:ndim))
  end subroutine register_field_wrapper

!> @brief Allocates the variable to be the size of data compute domain for x and y
  !! and for a given size for the 3rd 4th and 5th dimension
  subroutine var_data_alloc(var_data, dim2, dim3, dim4, dim5, compressed_dim)
    type(data_type),             intent(inout)    :: var_data         !< Variable data
    integer,                     intent(in)       :: dim2             !< Size of dim3
    integer,                     intent(in)       :: dim3             !< Size of dim3
    integer,                     intent(in)       :: dim4             !< Size of dim4
    integer,                     intent(in)       :: dim5             !< Size of dim5
    integer,                     intent(in)       :: compressed_dim   !< Size of the compressed dimension

    allocate(var_data%var_r4(compressed_dim, dim2, dim3, dim4, dim5))
    allocate(var_data%var_r8(compressed_dim, dim2, dim3, dim4, dim5))
    allocate(var_data%var_i4(compressed_dim, dim2, dim3, dim4, dim5))
    allocate(var_data%var_i8(compressed_dim, dim2, dim3, dim4, dim5))
  end subroutine var_data_alloc


  !> @brief Sets the data to some value
  subroutine var_data_set(var_data, var_value)
    type(data_type),             intent(inout)    :: var_data         !< Variable data
    integer,                     intent(in)       :: var_value        !< Value to set the data as

    var_data%var_r4 = real(var_value, kind=r4_kind)
    var_data%var_r8 = real(var_value, kind=r8_kind)
    var_data%var_i4 = int(var_value, kind=i4_kind)
    var_data%var_i8 = int(var_value, kind=i8_kind)
  end subroutine var_data_set

  !> @brief Writes the data for a give variable type
  subroutine write_data_wrapper(fileob, var_kind, var_data)
    type(FmsNetcdfFile_t),       intent(inout) :: fileob              !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_kind            !< The kind of the variable
    class(*),                    intent(in)    :: var_data(:,:,:,:,:) !< Variable data

    call write_data(fileob, "var1_"//trim(var_kind), var_data(:,1,1,1,1))
    call write_data(fileob, "var2_"//trim(var_kind), var_data(:,:,1,1,1))
    call write_data(fileob, "var3_"//trim(var_kind), var_data(:,:,:,1,1))
    call write_data(fileob, "var4_"//trim(var_kind), var_data(:,:,:,:,1))
    call write_data(fileob, "var5_"//trim(var_kind), var_data(:,:,:,:,:))

  end subroutine

  subroutine read_data_wrapper(fileob, var_name, dim, var_data, ref_data)
    type(FmsNetcdfFile_t),        intent(inout):: fileob              !< fms2io fileobj for domain decomposed
    character(len=*),            intent(in)    :: var_name            !< The kind of the variable
    integer,                     intent(in)    :: dim
    type(data_type),             intent(inout) :: var_data
    type(data_type),             intent(in)    :: ref_data

    integer :: i,j
    select case(dim)
    case(1)
      call var_data_set(var_data, -999)
      call read_data(fileob, trim(var_name)//"_r4", var_data%var_r4(:,1,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,1,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r4(:,1,1,1,1), (/mpp_pe()/)), "var2_r4")

      call read_data(fileob, trim(var_name)//"_r8", var_data%var_r8(:,1,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,1,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r8(:,1,1,1,1), (/mpp_pe()/)), "var2_r8")

      call read_data(fileob, trim(var_name)//"_i4", var_data%var_i4(:,1,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,1,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i4(:,1,1,1,1), (/mpp_pe()/)), "var2_i4")

      call read_data(fileob, trim(var_name)//"_i8", var_data%var_i8(:,1,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i8(:,1,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i8(:,1,1,1,1), (/mpp_pe()/)), "var2_i8")
    case(2)
      call var_data_set(var_data, -999)
      call read_data(fileob, trim(var_name)//"_r4", var_data%var_r4(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r4(:,:,1,1,1), (/mpp_pe()/)), "var2_r4")

      call read_data(fileob, trim(var_name)//"_r8", var_data%var_r8(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r8(:,:,1,1,1), (/mpp_pe()/)), "var2_r8")

      call read_data(fileob, trim(var_name)//"_i4", var_data%var_i4(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i4(:,:,1,1,1), (/mpp_pe()/)), "var2_i4")

      call read_data(fileob, trim(var_name)//"_i8", var_data%var_i8(:,:,1,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i8(:,:,1,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i8(:,:,1,1,1), (/mpp_pe()/)), "var2_i8")
    case(3)
      call var_data_set(var_data, -999)
      call read_data(fileob, trim(var_name)//"_r4", var_data%var_r4(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r4(:,:,:,1,1), (/mpp_pe()/)), "var3_r4")

      call read_data(fileob, trim(var_name)//"_r8", var_data%var_r8(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r8(:,:,:,1,1), (/mpp_pe()/)), "var3_r8")

      call read_data(fileob, trim(var_name)//"_i4", var_data%var_i4(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i4(:,:,:,1,1), (/mpp_pe()/)), "var3_i4")

      call read_data(fileob, trim(var_name)//"_i8", var_data%var_i8(:,:,:,1,1))
      call compare_var_data(mpp_chksum(var_data%var_i8(:,:,:,1,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i8(:,:,:,1,1), (/mpp_pe()/)), "var3_i8")
    case(4)
      call var_data_set(var_data, -999)
      call read_data(fileob, trim(var_name)//"_r4", var_data%var_r4(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,:,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r4(:,:,:,:,1), (/mpp_pe()/)), "var4_r4")

      call read_data(fileob, trim(var_name)//"_r8", var_data%var_r8(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,:,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r8(:,:,:,:,1), (/mpp_pe()/)), "var4_r8")

      call read_data(fileob, trim(var_name)//"_i4", var_data%var_i4(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i4(:,:,:,:,1), (/mpp_pe()/)), "var4_i4")

      call read_data(fileob, trim(var_name)//"_i8", var_data%var_i8(:,:,:,:,1))
      call compare_var_data(mpp_chksum(var_data%var_i8(:,:,:,:,1), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i8(:,:,:,:,1), (/mpp_pe()/)), "var4_i8")
    case(5)
      call var_data_set(var_data, -999)
      call read_data(fileob, trim(var_name)//"_r4", var_data%var_r4(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_r4(:,:,:,:,:), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r4(:,:,:,:,:), (/mpp_pe()/)), "var5_r4")

      call read_data(fileob, trim(var_name)//"_r8", var_data%var_r8(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_r8(:,:,:,:,:), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_r8(:,:,:,:,:), (/mpp_pe()/)), "var5_r8")

      call read_data(fileob, trim(var_name)//"_i4", var_data%var_i4(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_i4(:,:,:,:,:), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i4(:,:,:,:,:), (/mpp_pe()/)), "var5_i4")

      call read_data(fileob, trim(var_name)//"_i8", var_data%var_i8(:,:,:,:,:))
      call compare_var_data(mpp_chksum(var_data%var_i8(:,:,:,:,:), (/mpp_pe()/)), &
        mpp_chksum(ref_data%var_i8(:,:,:,:,:), (/mpp_pe()/)), "var5_i8")
    end select
  end subroutine

  subroutine var_data_set_ref(var_data)
    type(data_type),             intent(inout) :: var_data
    integer :: starting_index, ending_index, i,j

    ending_index = 0
    do i = 1, size(pes)
        starting_index = ending_index + 1
        ending_index = starting_index + pes(i)

        var_data%var_r4(starting_index:ending_index,:,:,:,:) = real(pes(i), kind=r4_kind)
        var_data%var_r8(starting_index:ending_index,:,:,:,:) = real(pes(i), kind=r8_kind)
        var_data%var_i4(starting_index:ending_index,:,:,:,:) = int(pes(i), kind=i4_kind)
        var_data%var_i8(starting_index:ending_index,:,:,:,:) = int(pes(i), kind=i8_kind)

    enddo
  end subroutine

  !> @brief Compares two checksums and crashes if they are not the same
  subroutine compare_var_data(check_sum_in, check_sum_ref, varname)
    integer(kind=i8_kind), intent(in) :: check_sum_in
    integer(kind=i8_kind), intent(in) :: check_sum_ref
    character(len=*), intent(in) :: varname

   if (check_sum_ref .ne. check_sum_in) call mpp_error(FATAL, &
     "Checksums do not match for variable: "//trim(varname))
  end subroutine compare_var_data

end program test_compressed_writes
