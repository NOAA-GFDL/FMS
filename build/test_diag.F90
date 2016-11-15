program test_diag

use mpp_mod
use time_manager_mod
use diag_manager_mod!, only:register_diag_field,send_data

implicit none

integer, parameter :: nx = 10
integer, parameter :: ny = 5

real :: x(nx)
real :: y(ny)

real, dimension (nx,ny) :: twod
real, dimension (nx) :: oned
real :: zerod

integer :: oned_axes(1)
integer :: twod_axes(2)


integer :: ii,diag_field_id,i,ix,iy

logical :: diag_result, tf

character (len=100) :: err
type diag_info
 character (len=15) :: modd
 character (len=15) :: varname
 character (len=100) :: longname
 integer  :: id
endtype diag_info
type(time_type) :: time
type(diag_info) :: rh
type(diag_info) :: shflx 
!type time_type
!   private
!   integer:: seconds
!   integer:: days
!   integer:: ticks
!   integer:: dummy ! added as a workaround bug on IRIX64 (AP)
!end type time_type

!time%seconds=3600
!time%days=1
!time%ticks=60
!time%dummy = 0

zerod = 123.4
oned = 111.1
twod = 222.2
do i=1,nx
x(i) = real(ii)
y(i) = real(ii)*2.0
enddo

!type diag_info
! character (len=15) :: modd
! character (len=15) :: varname
! character (len=100) :: longname
! integer  :: id
!endtype diag_info
 rh%modd="flux"
 rh%varname="rh_ref"
 rh%longname="Relativve Humidity"

 shflx%modd = "flux"
 shflx%varname = "shflx"
 shflx%longname = "Suface heat flux"


write (6,*) "Beginning test."

 call mpp_init()
 if ( mpp_pe() == mpp_root_pe() ) write (6,*) "init diag"

 time =  set_time(0, 0, 1, err)
write (6,*) err
 if ( mpp_pe() == mpp_root_pe() ) write (6,*) tf

 call diag_manager_init()
!  if ( mpp_pe() == mpp_root_pe() ) 
write (6,*) "init x-axis"
  !     INTEGER FUNCTION diag_axis_init(name, data, units, cart_name, long_name,
  !           direction, set_name, edges, Domain, Domain2, aux, tile_count)
 ix = diag_axis_init("lon",x,"num","X","x-direction")

 iy = diag_axis_init("lat",y,"num","Y","y-direction")


 oned_axes(1)=ix
 twod_axes(1)=ix
! twod_axes(1)=iy


! if ( mpp_pe() == mpp_root_pe() ) 
write (6,*) "ix = ", ix
write (6,*) "iy = ", iy


  rh%id = register_diag_field ("flux", "rh_ref") !> Returns the field index to be used in subsequent calls to send_data
                                              !! Is used as diag_field_id
!  shflx%id = register_diag_field (trim(shflx%modd) , trim(shflx%varname), oned_axes,time)

!! For our example here, rh_ref will be a scalar
! diag_result = send_data(diag_field_id,zerod,time)
!  diag_result = send_data(diag_field_id,zerod)

 if ( mpp_pe() == mpp_root_pe() ) write (6,*) "The result of register (field index) is ",rh%id
! ii =  get_diag_field_id("flux", "shflx")
!  if ( mpp_pe() == mpp_root_pe() ) write (6,*) "The result of get_diag_field_id is ",ii


! write(6,*) send_data(rh%id,zerod)

write (6,*) "End test."
end program test_diag
