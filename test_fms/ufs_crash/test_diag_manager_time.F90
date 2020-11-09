program test_diag_manager

use   mpp_domains_mod
use   diag_manager_mod
use   fms_mod
use  time_manager_mod, only: time_type, set_calendar_type, set_date, NOLEAP, JULIAN, operator(+), set_time, print_time

implicit none

type(time_type)                   :: Time
integer, dimension(2)             :: layout = (/1,1/)
integer :: nlon, nlat, nz
type(domain2d)                    :: Domain
real, dimension(:), allocatable :: x, y, z
integer :: i, j
integer :: is, ie, js, je
real, allocatable, dimension(:,:,:) :: sst, ice
integer :: id_x, id_y, id_z, id_sst, id_ice
integer :: used

call fms_init
call set_calendar_type(JULIAN)
call diag_manager_init

nlon = 20
nlat = 20
nz = 5

call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_diag_manager')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_compute_domain(Domain, is, ie, js, je)

! Set up the data
allocate(x(nlon), y(nlat), z(nz))
allocate(sst(is:ie,js:je,1:nz), ice(is:ie,js:je,1:nz))

do i=1,nlon
  x(i) = i
enddo
do j=1,nlat
  y(j) = j
enddo
do i=1,nz
   z(i) = i
enddo

sst = 666.66
ice = 619.0

! Set up the intial time
Time = set_date(2,1,1,0,0,0)

! Register the diags
id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)
id_z  = diag_axis_init('z',  z,  'point_Z', 'z', long_name='point_Z')
id_sst = register_diag_field  ('test_diag_manager_mod', 'sst', (/id_x,id_y,id_z/), Time, 'SST', 'K')

! Send the axis data
used = send_data(id_x, x, Time)
used = send_data(id_y, y, Time)
used = send_data(id_z, z, Time)

! Increase the time and send data
do i=1,20
Time = set_date(2,1,i,0,0,0)
if(id_sst > 0) used = send_data(id_sst, sst, Time)
enddo

call diag_manager_end(Time)
call fms_end

end program test_diag_manager
