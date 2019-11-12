subroutine unpack_data_0d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_0d


subroutine unpack_data_1d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), dimension(:), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_1d

subroutine unpack_data_2d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), dimension(:,:), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_2d

subroutine unpack_data_3d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), dimension(:,:,:), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_3d

subroutine unpack_data_4d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), dimension(:,:,:,:), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_4d

subroutine unpack_data_5d(fileobj, varid, var_data)
  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: varid
  class(*), dimension(:,:,:,:,:), intent(inout) :: var_data !< Array that the data
                                                 !! will be read into.

  integer(kind=int32) :: buf_int32
  integer(kind=int64) :: buf_int64
  real(kind=real32) :: buf_real32
  real(kind=real64) :: buf_real64
  integer :: err

! 		First, check if scale factor exists, if it does assign it to a buffer and multiply 
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "scale_factor")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_int64 )
			var_data = var_data * buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data * buf_real64
      class default
    end select
  endif

! 		Then, check if add_offset exists, if it does assign it to a buffer and add
! 		it to the data

  if (attribute_exists(fileobj%ncid, varid, "add_offset")) then
    select type(var_data)
      type is (integer(kind=int32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_int32
      type is (integer(kind=int64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_int64 )
			var_data = var_data + buf_int64
    		call check_netcdf_code(err)
      type is (real(kind=real32))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real32 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real32
      type is (real(kind=real64))
	      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_real64 )
    		call check_netcdf_code(err)
			var_data = var_data + buf_real64
      class default
    end select
  endif

end subroutine unpack_data_5d
