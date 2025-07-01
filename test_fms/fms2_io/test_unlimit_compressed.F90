!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief  This programs tests the register_unlimited_compressed_axis feature in
!! fms2io
program test_unlimit_compressed

use fms2_io_mod, only: open_file, register_unlimited_compressed_axis, register_restart_field, &
                       write_restart, read_restart, close_file, FmsNetcdfDomainFile_t, &
                       get_dimension_size
use mpp_mod, only: mpp_pe, mpp_error, FATAL, mpp_sum
use mpp_domains_mod, only: mpp_define_domains, mpp_define_io_domain, domain2d
use fms_mod, only: fms_init, fms_end

implicit none

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io fileobj
type(domain2d)              :: Domain         !< Domain
integer, allocatable        :: vdata(:)       !< Variable data
integer                     :: dimsize        !< Size of the dimension

call fms_init()

call mpp_define_domains( (/1,96,1,96/), (/1,6/), Domain)
call mpp_define_io_domain(Domain, (/1,2/))

allocate(vdata(mpp_pe())) !< The size of vdata each different for each PE
vdata = mpp_pe()

!< Writes
if (open_file(fileobj, "filename.nc", "overwrite", Domain, is_restart=.true.)) then
  call register_unlimited_compressed_axis(fileobj, "i", mpp_pe())
  call register_restart_field(fileobj, "var", vdata, (/"i"/))
  call write_restart(fileobj)

  call close_file(fileobj)
endif

deallocate(vdata)

!< Reads
if (open_file(fileobj, "filename.nc", "read", Domain, is_restart=.true.)) then
  call get_dimension_size(fileobj, "i", dimsize)

  allocate(vdata(dimsize))

  call register_restart_field(fileobj, "var", vdata, (/"i"/))
  call read_restart(fileobj)
  call close_file(fileobj)
endif

!< Check if it worked:
select case (mpp_pe())
case (0, 1, 2)
  !< PE 0, 1, and 2 are going to be reading the .0001 file and the size of i for
  !! that file is 3
  if (dimsize .ne. 3) call mpp_error(FATAL, "dimsize is not the correct the size")
case (3, 4, 5)
  !< PE 3, 4, and 5 are going to be reading the .0002 file and the size of i for
  !! that file is 12
  if (dimsize .ne. 12) call mpp_error(FATAL, "dimsize is not the correct the size")
end select

call fms_end()

end program test_unlimit_compressed
