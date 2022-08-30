program test

use fms2_io_mod
use mpp_mod
use mpp_domains_mod
use fms_mod, only: fms_init, fms_end

implicit none

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io fileobj
type(domain2d)              :: Domain         !< Domain
integer, allocatable :: vdata(:)

call fms_init()

call mpp_define_domains( (/1,96,1,96/), (/1,6/), Domain)
call mpp_define_io_domain(Domain, (/1,2/))

allocate(vdata(mpp_pe()))
vdata = mpp_pe()

if (open_file(fileobj, "filename.nc", "overwrite", Domain, is_restart=.true.)) then
  call register_axis(fileobj, "i", mpp_pe(), is_compressed=.true.)
  call register_restart_field(fileobj, "var", vdata, (/"i"/))
  call write_restart(fileobj)

  call close_file(fileobj)
endif

call fms_end()

end program
