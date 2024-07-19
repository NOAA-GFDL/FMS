program test_data_override_weights

use fms_mod, only: fms_init, fms_end
use data_override_mod
use time_manager_mod
use mpp_mod
use mpp_domains_mod
use platform_mod

implicit none

integer, parameter       :: lkind = DO_TEST_KIND_
integer                  :: nlat, nlon           !< Number of longitudes and latitudes
integer                  :: is, ie, js, je       !< The starting and ending indices of the data_domain
type(domain2d)           :: domain               !< Domain
real(lkind), allocatable :: sst_obs(:,:)         !< Data overriden using data_override
real(lkind), allocatable :: sst_obs_weights(:,:) !< Data overriden using the weight files in data_override
type(time_type)          :: Time                 !< The time of the simulation
logical                  :: success              !< .True. if the data_override was successful
integer                  :: i, j                 !< For do loops

call fms_init()
call set_calendar_type(NOLEAP)

nlon = 1440
nlat = 1080

!< Create a domain nlonXnlat with mask
call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), (/1,1/), Domain, name='test_data_override_weights')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, is, ie, js, je)

call data_override_init(ice_domain_in=domain, mode=lkind)
Time = set_date(1994,6,6,0,0,0)
allocate(sst_obs(is:ie,js:je))
sst_obs = 999.
call data_override('ICE','sst_obs',sst_obs, Time, override=success)
if (.not. success) call mpp_error(FATAL, "Data override failed")

allocate(sst_obs_weights(is:ie,js:je))
sst_obs_weights = 999.
call data_override('ICE','sst_obs_weights',sst_obs_weights, Time, override=success)
if (.not. success) call mpp_error(FATAL, "Data override failed for when using weights")

do i = is, ie
  do j = js, je
    if (abs(sst_obs_weights(i,j) - sst_obs(i,j)) > 1e-06) then
      print *, "i = ", i, " j = ", j, " ", sst_obs_weights(i,j), " vs ", sst_obs(i,j)
      print *, "diff =", abs(sst_obs_weights(i,j) - sst_obs(i,j))
      call mpp_error(FATAL, "The data is not the expected result")
    endif
  enddo
enddo

call fms_end()

end program test_data_override_weights