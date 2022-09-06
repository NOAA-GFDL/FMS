program test

use fms2_io_mod
use mpp_mod
use mpp_domains_mod
use fms_mod, only: fms_init, fms_end

implicit none

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io fileobj
type(domain2d)              :: Domain         !< Domain
integer, allocatable :: vdata(:)
integer :: dimsize

call fms_init()

call mpp_define_domains( (/1,96,1,96/), (/1,6/), Domain)
call mpp_define_io_domain(Domain, (/1,2/))

allocate(vdata(mpp_pe())) !< this is allcate(vdata(nbergs)) in icebergs
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
write (mpp_pe()+100, *) dimsize
write (mpp_pe()+100, *) vdata

call fms_end()

end program
