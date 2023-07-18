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
program test_packed_reads
  use fms_mod,      only: fms_init, fms_end
  use platform_mod, only: r8_kind, r4_kind, i8_kind
  use mpp_mod,      only: mpp_error, FATAL, mpp_chksum
  use fms2_io_mod,  only: open_file, close_file, FmsNetcdfFile_t, read_data
  use netcdf

  implicit none

  real(kind=r8_kind)    :: var_r8(10, 1, 1, 1, 1)      !< original r8 data
  real(kind=r8_kind)    :: var_r8_out(10, 1, 1, 1, 1)  !< expected r8 data
  real(kind=r8_kind)    :: var_r8_out2(10, 1, 1, 1, 1) !< r8 data from file
  real(kind=r8_kind)    :: scale_factor_r8             !< r8 scale factor
  real(kind=r8_kind)    :: add_offset_r8               !< r8 offset

  real(kind=r4_kind)    :: var_r4(10, 1, 1, 1, 1)      !< original r4 data
  real(kind=r4_kind)    :: var_r4_out(10, 1, 1, 1, 1)  !< expected r4 data
  real(kind=r4_kind)    :: var_r4_out2(10, 1, 1, 1, 1) !< r4 data from file
  real(kind=r4_kind)    :: scale_factor_r4             !< r4 scale factor
  real(kind=r4_kind)    :: add_offset_r4               !< r4 offset

  type(FmsNetcdfFile_t) :: fileobj !< Fms2_io fileobj

  integer(2)            :: packed_data_r8(10, 1, 1, 1, 1) !< Packed data calculated from r8
  integer(2)            :: packed_data_r4(10, 1, 1, 1, 1) !< Packed data calculated from r4

  integer :: ncid                                                   !< netcdf file id
  integer :: varid1_r8, varid2_r8, varid3_r8, varid4_r8, varid5_r8  !< Variable ids
  integer :: varid1_r4, varid2_r4, varid3_r4, varid4_r4, varid5_r4  !< variable ids
  integer :: dimid1, dimid2, dimid3, dimid4, dimid5                 !< Dimension ids
  integer :: dim(5)                                                 !< Array of dimension ids
  integer :: i                                                      !< For do loops
  integer :: err                                                    !< netcdf error code

  call fms_init()

  do i = 1, size(var_r4, 1)
    var_r4(i,:,:,:,:) = 90_r4_kind - real(i, kind=r4_kind)/13_r4_kind + real(i, kind=r4_kind)
    var_r8(i,:,:,:,:) = 90_r8_kind - real(i, kind=r8_kind)/13_r8_kind + real(i, kind=r8_kind)
  enddo

  scale_factor_r4 = real((maxval(var_r4) - minval(var_r4)) / (2**(16-1)), kind=r4_kind)
  add_offset_r4 = real(minval(var_r4) + 2 ** (16 - 1) * scale_factor_r4, kind=r4_kind)
  packed_data_r4 = nint((var_r4 - add_offset_r4) / scale_factor_r4, kind=2)
  var_r4_out = packed_data_r4*scale_factor_r4 + add_offset_r4

  scale_factor_r8 = real((maxval(var_r8) - minval(var_r8)) / (2**(16-1)), kind=r8_kind)
  add_offset_r8 = real(minval(var_r8) + 2 ** (16 - 1) * scale_factor_r8, kind=r8_kind)
  packed_data_r8 = nint((var_r8 - add_offset_r8) / scale_factor_r8, kind=2)
  var_r8_out = packed_data_r8*scale_factor_r8 + add_offset_r8

  err = nf90_create("short_file.nc", ior(nf90_clobber, nf90_64bit_offset), ncid)
  err = nf90_def_dim(ncid, 'dim1', 10, dimid1)
  err = nf90_def_dim(ncid, 'dim2', 1, dimid2)
  err = nf90_def_dim(ncid, 'dim3', 1, dimid3)
  err = nf90_def_dim(ncid, 'dim4', 1, dimid4)
  err = nf90_def_dim(ncid, 'dim5', 1, dimid5)

  dim = (/dimid1, dimid2, dimid3, dimid4, dimid5/)

  call write_var_metadata(ncid, 'var_1d', dim(1:1), scale_factor_r4, scale_factor_r8, &
    add_offset_r4, add_offset_r8, varid1_r8, varid1_r4)
  call write_var_metadata(ncid, 'var_2d', dim(1:2), scale_factor_r4, scale_factor_r8, &
    add_offset_r4, add_offset_r8, varid2_r8, varid2_r4)
  call write_var_metadata(ncid, 'var_3d', dim(1:3), scale_factor_r4, scale_factor_r8, &
    add_offset_r4, add_offset_r8, varid3_r8, varid3_r4)
  call write_var_metadata(ncid, 'var_4d', dim(1:4), scale_factor_r4, scale_factor_r8, &
    add_offset_r4, add_offset_r8, varid4_r8, varid4_r4)
  call write_var_metadata(ncid, 'var_5d', dim(1:5), scale_factor_r4, scale_factor_r8, &
    add_offset_r4, add_offset_r8, varid5_r8, varid5_r4)

  call check(nf90_enddef(ncid))

  call check(nf90_put_var(ncid, varid1_r8, packed_data_r8(:,1,1,1,1)))
  call check(nf90_put_var(ncid, varid2_r8, packed_data_r8(:,:,1,1,1)))
  call check(nf90_put_var(ncid, varid3_r8, packed_data_r8(:,:,:,1,1)))
  call check(nf90_put_var(ncid, varid4_r8, packed_data_r8(:,:,:,:,1)))
  call check(nf90_put_var(ncid, varid5_r8, packed_data_r8(:,:,:,:,:)))

  call check(nf90_put_var(ncid, varid1_r4, packed_data_r4(:,1,1,1,1)))
  call check(nf90_put_var(ncid, varid2_r4, packed_data_r4(:,:,1,1,1)))
  call check(nf90_put_var(ncid, varid3_r4, packed_data_r4(:,:,:,1,1)))
  call check(nf90_put_var(ncid, varid4_r4, packed_data_r4(:,:,:,:,1)))
  call check(nf90_put_var(ncid, varid5_r4, packed_data_r4(:,:,:,:,:)))

  call check(nf90_close(ncid))

  if (open_file(fileobj, "short_file.nc", "read")) then
    var_r8_out2 = -999_r8_kind
    call read_data(fileobj, "var_1d_r8", var_r8_out2(:,1,1,1,1))
    call compare_data(mpp_chksum(var_r8_out2(:,1,1,1,1)), mpp_chksum(var_r8_out(:,1,1,1,1)), "var_1d_r8")

    var_r8_out2 = -999_r8_kind
    call read_data(fileobj, "var_2d_r8", var_r8_out2(:,:,1,1,1))
    call compare_data(mpp_chksum(var_r8_out2(:,:,1,1,1)), mpp_chksum(var_r8_out(:,:,1,1,1)), "var_1d_r8")

    var_r8_out2 = -999_r8_kind
    call read_data(fileobj, "var_3d_r8", var_r8_out2(:,:,:,1,1))
    call compare_data(mpp_chksum(var_r8_out2(:,:,:,1,1)), mpp_chksum(var_r8_out(:,:,:,1,1)), "var_1d_r8")

    var_r8_out2 = -999_r8_kind
    call read_data(fileobj, "var_4d_r8", var_r8_out2(:,:,:,:,1))
    call compare_data(mpp_chksum(var_r8_out2(:,:,:,:,1)), mpp_chksum(var_r8_out(:,:,:,:,1)), "var_1d_r8")

    var_r8_out2 = -999_r8_kind
    call read_data(fileobj, "var_5d_r8", var_r8_out2(:,:,:,:,:))
    call compare_data(mpp_chksum(var_r8_out2(:,:,:,:,:)), mpp_chksum(var_r8_out(:,:,:,:,:)), "var_1d_r8")

    var_r4_out2 = -999_r4_kind
    call read_data(fileobj, "var_1d_r4", var_r4_out2(:,1,1,1,1))
    call compare_data(mpp_chksum(var_r4_out2(:,1,1,1,1)), mpp_chksum(var_r4_out(:,1,1,1,1)), "var_1d_r4")

    var_r4_out2 = -999_r4_kind
    call read_data(fileobj, "var_2d_r4", var_r4_out2(:,:,1,1,1))
    call compare_data(mpp_chksum(var_r4_out2(:,:,1,1,1)), mpp_chksum(var_r4_out(:,:,1,1,1)), "var_1d_r4")

    var_r4_out2 = -999_r4_kind
    call read_data(fileobj, "var_3d_r4", var_r4_out2(:,:,:,1,1))
    call compare_data(mpp_chksum(var_r4_out2(:,:,:,1,1)), mpp_chksum(var_r4_out(:,:,:,1,1)), "var_1d_r4")

    var_r4_out2 = -999_r4_kind
    call read_data(fileobj, "var_4d_r4", var_r4_out2(:,:,:,:,1))
    call compare_data(mpp_chksum(var_r4_out2(:,:,:,:,1)), mpp_chksum(var_r4_out(:,:,:,:,1)), "var_1d_r4")

    var_r4_out2 = -999_r4_kind
    call read_data(fileobj, "var_5d_r4", var_r4_out2(:,:,:,:,:))
    call compare_data(mpp_chksum(var_r4_out2(:,:,:,:,:)), mpp_chksum(var_r4_out(:,:,:,:,:)), "var_1d_r4")


    call close_file(fileobj)
  endif

  call fms_end()

  contains

  !> @brief Write out the variable data
  subroutine write_var_metadata(fileid, varname, dimids, sfactor_r4, sfactor_r8, &
    offset_r4, offset_r8, varid_r4, varid_r8)
    integer,               intent(in)  :: fileid     !< netcdf file id
    character(len=*),      intent(in)  :: varname    !< variable name
    integer,               intent(in)  :: dimids(:)  !< array of the dimension ids
    real(kind=r4_kind),    intent(in)  :: sfactor_r4 !< Scale factor in r4 precision
    real(kind=r8_kind),    intent(in)  :: sfactor_r8 !< Scale factor in r8 precision
    real(kind=r4_kind),    intent(in)  :: offset_r4  !< offset in r4 precision
    real(kind=r8_kind),    intent(in)  :: offset_r8  !< offset in r8 precision
    integer,               intent(out) :: varid_r4   !< variable id for the r4 variable
    integer,               intent(out) :: varid_r8   !< variable id for the r8 variable

    call check(nf90_def_var(ncid, trim(varname)//"_r4", nf90_short, dimids, varid_r4))
    call check(nf90_put_att(ncid, varid_r4, "scale_factor", sfactor_r4))
    call check(nf90_put_att(ncid, varid_r4, "add_offset", offset_r4))

    call check(nf90_def_var(ncid, trim(varname)//"_r8", nf90_short, dimids, varid_r8))
    call check(nf90_put_att(ncid, varid_r8, "scale_factor", sfactor_r8))
    call check(nf90_put_att(ncid, varid_r8, "add_offset", offset_r8))

  end subroutine write_var_metadata

  !> @brief Compare two checksums
  subroutine compare_data(checksum_in, checksum_out, varname)
    integer(kind=i8_kind), intent(in) :: checksum_in    !< The data checksum
    integer(kind=i8_kind), intent(in) :: checksum_out   !< The reference checksum
    character(len=*),      intent(in) :: varname        !< The variable's name (for error messages)

   if (checksum_in .ne. checksum_out) call mpp_error(FATAL, &
     "Checksums do not match for variable: "//trim(varname))
  end subroutine

  !> @brief Check the netcdf error code
  subroutine check(status)
    integer, intent ( in) :: status !< netcdf error code

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program test_packed_reads
