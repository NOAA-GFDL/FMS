module spherical_regrid_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov">
!M.J. Harrison
!</CONTACT>
!
!<REVIEWER EMAIL="Anand.Gnanadesikan@noaa.gov">
!Anand Gnandasekian
!</REVIEWER>
!

!<OVERVIEW>
! Map between logically rectangular grids using inverse great circle 
! weights.
!</OVERVIEW>

!<DESCRIPTION>
! A maximum range of influence (max_dist) is selected 
! along with the number of nearest neighbors (num_nbrs)
!
! Currently only supports a single regridding map 
! at a time.  Subsequent calls to spherical_regrid_init
! will erase previous mapping information
!
! Optional masking information for the source grid is used to mask destination
! grid points.  If the nearest neighbor from the source grid is masked, the
! destination point is masked, otherwise if destination point has at least one
! valid point in it's region of influence, it obtains a valid value.
!
! Schematic:
!    The following example illustrates a destination grid location (+) with 
!    a (R)adius of influence (in radians) denoted by (=).  Valid source grid 
!    locations (o) which fall within the radius of influence of the destination 
!    point are used in the mapping.  Masked points (x) do not contribute to the mapping
!    but are included as nearest neighbors. In this case, 4 valid source grid points fall
!    within the radius of influence.   
!
!<PRE>   
!               o  o   o   o   o   o   o
!                      =========
!                     =         = 
!               o  o = o   o   o = o   o
!                   =        R    = 
!                   =      +----->=
!               o  o   o   x   x   x   x
!                     =         =
!                      =========
!               o  o   o   x   x   x   x
!
!</PRE>
!
!</DESCRIPTION>
!

use mpp_mod, only : mpp_error, NOTE, FATAL, mpp_pe,mpp_root_pe

implicit none

private

real, parameter :: max_dist_default = 0.17  ! radians
real, parameter :: epsln=1.e-10, pi=3.14159265358979323846, large=1.e20, pih=pi/2.,pi2=pi*2.
real, parameter :: deg2rad=pi/180.
integer, public :: num_neighbors, map_dest_size, map_src_size, map_dest_xsize, map_dest_ysize, map_src_xsize, map_src_ysize
integer, dimension(:,:), allocatable, public :: map_src_add
logical, dimension(:), allocatable :: lmask_src, lmask_dest
real, dimension(:,:), allocatable :: map_src_dist
real, dimension(:,:), allocatable, public :: map_src_wgt
real :: max_src_dist = max_dist_default
logical :: spherical_regrid_initialized = .false., radial_src_search = .true., src_is_modulo = .true.

real, dimension(:), allocatable :: dest_1d
real, dimension(:), allocatable ::  src_1d

public spherical_regrid_init, regrid, spherical_distance, spherical_regrid_exit

interface spherical_regrid_init
   module procedure spherical_regrid_init_1d
   module procedure spherical_regrid_init_1d_src
   module procedure spherical_regrid_init_1d_dest
   module procedure spherical_regrid_init_2d
end interface


contains

!
!<SUBROUTINE NAME="regrid">
!
!<DESCRIPTION>
! map an array from src grid to dest grid
!</DESCRIPTION>
!
!<IN NAME="src" TYPE="real, pointer" DIM="(:,:)">
! data on source grid
!</IN>
!<INOUT NAME="dest" TYPE="real, pointer" DIM="(:,:)">
! data on destination grid
!</INOUT>
subroutine regrid(dest,src)

real, dimension(:,:), pointer :: dest
real, dimension(:,:), pointer :: src

integer :: n1,n2, i, j


n1 = size(dest,1)*size(dest,2)
n2 = size(src,1)*size(src,2)

if (n1 /= map_dest_size .or. n2 /= map_src_size) call mpp_error(FATAL,'=> incorrect array size in call to remap')


if (.NOT. allocated(dest_1d)) allocate(dest_1d(map_dest_size))
if (.NOT. allocated(src_1d)) allocate(src_1d(map_src_size))

src_1d = reshape(src,(/map_src_size/))
dest_1d=0.0


do i=1,map_dest_size
   if (lmask_dest(i)) then
      do j=1,num_neighbors
         if (map_src_add(i,j) /= 0) then
            dest_1d(i) = dest_1d(i)+src_1d(map_src_add(i,j))*map_src_wgt(i,j)
         endif
      enddo
   endif
enddo

dest = reshape(dest_1d,(/map_dest_xsize,map_dest_ysize/))

return

end subroutine regrid
!</SUBROUTINE> NAME="regrid"


subroutine spherical_regrid_init_1d(x_dest,y_dest,x_src,y_src,num_nbrs,mask_src,max_dist,src_modulo)

! x_dest,y_dest = destination grid lon,lat
! x_src,y_src = source grid lon,lat

real, dimension(:), pointer :: x_dest, x_src, y_dest, y_src
integer, intent(in) :: num_nbrs
real, dimension(:,:), optional, pointer :: mask_src
real, intent(in), optional :: max_dist
logical, intent(in), optional :: src_modulo

real, dimension(:,:), pointer :: x_dest_out, y_dest_out
real, dimension(:,:), pointer :: x_src_out, y_src_out

integer :: nx_src,ny_src,nx_dest,ny_dest,i,j,len

nx_src = size(x_src(:))
ny_src = size(y_src(:))
nx_dest = size(x_dest(:))
ny_dest = size(y_dest(:))

allocate(x_dest_out(nx_dest,ny_dest))
allocate(y_dest_out(nx_dest,ny_dest))
allocate(x_src_out(nx_src,ny_src))
allocate(y_src_out(nx_src,ny_src))

do j=1,ny_src
   x_src_out(:,j) = x_src(:)
enddo

do i=1,nx_src
   y_src_out(i,:) = y_src(:)
enddo

do j=1,ny_dest
   x_dest_out(:,j) = x_dest(:)
enddo

do i=1,nx_dest
   y_dest_out(i,:) = y_dest(:)
enddo

call spherical_regrid_init_2d(x_dest_out,y_dest_out,x_src_out,y_src_out,num_nbrs,mask_src,max_dist,src_modulo)

deallocate(x_dest_out, y_dest_out, x_src_out, y_src_out)

return

end subroutine spherical_regrid_init_1d

subroutine spherical_regrid_init_1d_src(x_dest,y_dest,x_src,y_src,num_nbrs,mask_src,max_dist,src_modulo)

! x_dest,y_dest = destination grid lon,lat
! x_src,y_src = source grid lon,lat

real, dimension(:,:), pointer :: x_dest, y_dest
real, dimension(:), pointer :: x_src, y_src
integer, intent(in) :: num_nbrs
real, dimension(:,:), optional, pointer :: mask_src
real, intent(in), optional :: max_dist
logical, intent(in), optional :: src_modulo

real, dimension(:,:), pointer :: x_src_out, y_src_out

integer :: nx_src,ny_src,nx_dest, ny_dest, i,j,len

nx_src = size(x_src(:))
ny_src = size(y_src(:))
nx_dest = size(x_dest,1)
ny_dest = size(x_dest,2)

allocate(x_src_out(nx_src,ny_src))
allocate(y_src_out(nx_src,ny_src))

do j=1,ny_src
   x_src_out(:,j) = x_src(:)
enddo

do i=1,nx_src
   y_src_out(i,:) = y_src(:)
enddo

call spherical_regrid_init_2d(x_dest,y_dest,x_src_out,y_src_out,num_nbrs,mask_src,max_dist,src_modulo)

deallocate(x_src_out, y_src_out)

return

end subroutine spherical_regrid_init_1d_src

subroutine spherical_regrid_init_1d_dest(x_dest,y_dest,x_src,y_src,num_nbrs,mask_src,max_dist,src_modulo)

! x_dest,y_dest = destination grid lon,lat
! x_src,y_src = source grid lon,lat

real, dimension(:), pointer :: x_dest, y_dest
real, dimension(:,:), pointer :: x_src, y_src
integer, intent(in) :: num_nbrs
real, dimension(:,:), optional, pointer :: mask_src
real, intent(in), optional :: max_dist
logical, intent(in), optional :: src_modulo

real, dimension(:,:), pointer :: x_dest_out, y_dest_out

integer :: nx_dest,ny_dest, nx_src, ny_src,i,j, len

nx_src = size(x_src,1)
ny_src = size(x_src,2)
nx_dest = size(x_dest(:))
ny_dest = size(y_dest(:))

allocate(x_dest_out(nx_dest,ny_dest))
allocate(y_dest_out(nx_dest,ny_dest))

do j=1,ny_dest
   x_dest_out(:,j) = x_dest(:)
enddo

do i=1,nx_dest
   y_dest_out(i,:) = y_dest(:)
enddo

call spherical_regrid_init_2d(x_dest_out,y_dest_out,x_src,y_src,num_nbrs,mask_src,max_dist,src_modulo)

deallocate(x_dest_out, y_dest_out)

return

end subroutine spherical_regrid_init_1d_dest

!<SUBROUTINE NAME="spherical_regrid_init">
!
!<DESCRIPTION>
! initialize spherical regrid information.  Identify nearest points and calculate 
! regriding weights.
!</DESCRIPTION>
!
!<IN NAME="x_dest" TYPE="real, pointer" DIM="(:),(:,:)" UNITS="degrees">
! longitudes of destination grid 
!</IN>
!<IN NAME="y_dest" TYPE="real, pointer" DIM="(:),(:,:)" UNITS="degrees">
! latitudes of destination grid
!</IN>
!<IN NAME="x_src" TYPE="real, pointer" DIM="(:),(:,:)" UNITS="degrees">
! longitudes of source grid 
!</IN>
!<IN NAME="y_src" TYPE="real, pointer" DIM="(:),(:,:)" UNITS="degrees">
! latitudes of source grid
!</IN>
!<IN NAME="num_nbrs" TYPE="integer">
! maximum number of neighbors for regridding
!</IN>
!<IN NAME="mask_src" TYPE="real, pointer" DIM="(:,:)">
! real mask for source grid (optional)
!</IN>
!<IN NAME="max_dist" TYPE="real" UNITS="radians">
! radius of influence around destination grid points (optional)
!</IN>
!<IN NAME="src_modulo" TYPE="logical">
! .T. if longitudes in source grid are modulo (optional)
!</IN>

subroutine spherical_regrid_init_2d(x_dest,y_dest,x_src,y_src,num_nbrs,mask_src,max_dist,src_modulo)

! x_dest,y_dest = destination grid lon,lat
! x_src,y_src = source grid lon,lat

real, dimension(:,:), pointer :: x_dest, x_src, y_dest, y_src
integer, intent(in) :: num_nbrs
real, dimension(:,:), optional, pointer :: mask_src
real, optional :: max_dist
logical, intent(in), optional :: src_modulo

integer :: i, j, n, m, nx_dest, ny_dest, nx_src, ny_src, k, num
integer :: bound, bound_start, bound_end, i0, j0, i_left, i_right
real, dimension(:), allocatable :: theta_dest,phi_dest
real, dimension(:), allocatable :: theta_src,phi_src
real :: min_theta_dest, max_theta_dest, min_phi_dest, max_phi_dest
real :: min_theta_src, max_theta_src, min_phi_src, max_phi_src
logical :: continue_search, found_neighbors, continue_radial_search, result, debug
real :: d, sum, dtheta, dphi,nearest,res
integer :: step, i_nearest, len, step_size

map_dest_xsize=size(x_dest,1);map_dest_ysize=size(x_dest,2)
map_src_xsize=size(x_src,1);map_src_ysize=size(x_src,2)
map_dest_size = map_dest_xsize*map_dest_ysize
map_src_size = map_src_xsize*map_src_ysize


if (num_nbrs <= 0) call mpp_error(FATAL,'num_neighbors must be > 0')

num_neighbors = num_nbrs

max_src_dist = max_dist_default
if (PRESENT(max_dist)) max_src_dist = max_dist

if (PRESENT(src_modulo)) src_is_modulo = src_modulo

allocate(theta_dest(map_dest_size))
allocate(phi_dest(map_dest_size))
allocate(theta_src(map_src_size))
allocate(phi_src(map_src_size))

if (map_dest_size /= size(y_dest,1)*size(y_dest,2) .or. map_src_size /= size(y_src,1)*size(y_src,2)) &
     call mpp_error(FATAL, '=> grids not conformable')

if (spherical_regrid_initialized) then
    call mpp_error(NOTE,'recalculating regrid weights ...')
    deallocate(map_src_add,map_src_dist, &
         map_src_wgt,lmask_src, lmask_dest)
endif


allocate(map_src_add(map_dest_size,num_neighbors))
allocate(map_src_dist(map_dest_size,num_neighbors))
allocate(map_src_wgt(map_dest_size,num_neighbors))
allocate(lmask_src(map_src_size),lmask_dest(map_dest_size))

lmask_src=.true.
if (PRESENT(mask_src)) then
   if (size(mask_src,1)*size(mask_src,2) /= map_src_size) call mpp_error(FATAL,'mask not conformable')
   n=0
   do j=1,map_src_ysize
      do i=1,map_src_xsize
         n=n+1
         if (mask_src(i,j) == 0.) lmask_src(n) = .false.
      enddo
   enddo
endif

lmask_dest = .true.

theta_dest = reshape(x_dest,(/map_dest_size/))*deg2rad
theta_src = reshape(x_src,(/map_src_size/))*deg2rad
phi_dest = reshape(y_dest,(/map_dest_size/))*deg2rad
phi_src = reshape(y_src,(/map_src_size/))*deg2rad

map_src_add = 0
map_src_dist = large
map_src_wgt = 0.0

min_theta_dest=pi2;max_theta_dest=0.;min_phi_dest=pi;max_phi_dest=-pi
min_theta_src=pi2;max_theta_src=0.;min_phi_src=pi;max_phi_src=-pi

where(theta_dest<0.0)  theta_dest = theta_dest+pi2
where(theta_dest>pi2)  theta_dest = theta_dest-pi2
where(theta_src<0.0)  theta_src = theta_src+pi2
where(theta_src>pi2)  theta_src = theta_src-pi2

where(phi_dest < -pih) phi_dest = -pih
where(phi_dest > pih)  phi_dest =  pih
where(phi_src < -pih) phi_src = -pih
where(phi_src > pih)  phi_src =  pih    

do j=1,map_dest_size
   min_theta_dest = min(min_theta_dest,theta_dest(j))
   max_theta_dest = max(max_theta_dest,theta_dest(j))
   min_phi_dest = min(min_phi_dest,phi_dest(j))
   max_phi_dest = max(max_phi_dest,phi_dest(j))
enddo

do i=1,map_src_size
   min_theta_src = min(min_theta_src,theta_src(i))
   max_theta_src = max(max_theta_src,theta_src(i))
   min_phi_src = min(min_phi_src,phi_src(i))
   max_phi_src = max(max_phi_src,phi_src(i))
enddo

if (min_phi_dest < min_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
if (max_phi_dest > max_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
if (min_theta_dest < min_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
if (max_theta_dest > max_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'

do j=1,map_dest_size
   found_neighbors=.false.
   continue_search=.true.
   step = 1
   step_size = 100
   nearest = 1.e3
   do while (continue_search .and. step_size > 0)
      do while (step <= map_src_size .and. continue_search)
         ! count land points as nearest neighbors
         d = spherical_distance(theta_dest(j),phi_dest(j),theta_src(step),phi_src(step))
         if (d <= max_src_dist) then
            found_neighbors = update_dest_neighbors(j,step,d)
            if (found_neighbors) then
               n = 0
               i0 = mod(step,map_src_xsize)
               if (i0 == 0) i0 = map_src_xsize
               res = float(step)/float(map_src_xsize)
               j0 = ceiling(res)
               continue_radial_search = .true.
               do while (continue_radial_search)
                  n = n+1 ! radial counter
                  ! left boundary 
                  i_left = i0-n
                  if (i_left <= 0) then
                     if (src_is_modulo) then
                        i_left = map_src_xsize + i_left
                     else
                        i_left = 1
                     endif
                  endif
                  bound_start = max(j0-n-1,0)*map_src_xsize + i_left
                  bound_end = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_left
                  if (bound_end > map_src_size) call mpp_error(FATAL)
                  bound = bound_start
                  continue_radial_search = .false.
                  do while (bound <= bound_end)
                     d = spherical_distance(theta_dest(j),phi_dest(j),theta_src(bound),phi_src(bound))
                     result = update_dest_neighbors(j,bound,d)
                     bound = bound + map_src_xsize
                     if (result) continue_radial_search = .true.
                     if (result) found_neighbors = .true.
                  enddo
                  ! right boundary 
                  i_right = i0+n
                  if (i_right > map_src_xsize) then
                     if (src_is_modulo) then
                        i_right = i_right - map_src_xsize
                     else
                        i_right = map_src_xsize
                     endif
                  endif
                  bound_start = max(j0-n-1,0)*map_src_xsize + i_right
                  bound_end = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_right
                  bound = bound_start
                  if (bound_end > map_src_size) call mpp_error(FATAL)
                  do while (bound <= bound_end)
                     d = spherical_distance(theta_dest(j),phi_dest(j),theta_src(bound),phi_src(bound))
                     result = update_dest_neighbors(j,bound,d)
                     bound = bound + map_src_xsize
                     if (result) continue_radial_search = .true.
                     if (result) found_neighbors = .true.
                  enddo
                  ! bottom boundary 
                  bound_start = max(j0-n-1,0)*map_src_xsize + i_left 
                  bound_end =  max(j0-n-1,0)*map_src_xsize  + i_right 
                  if (bound_start > bound_end) then
                     bound_start = max(j0-n-1,0)*map_src_xsize + 1
                     bound_end = max(j0-n,1)*map_src_xsize
                  endif
                  bound = bound_start
                  if (bound_end > map_src_size) call mpp_error(FATAL)
                  do while (bound <= bound_end)
                     d = spherical_distance(theta_dest(j),phi_dest(j),theta_src(bound),phi_src(bound))
                     result = update_dest_neighbors(j,bound,d)
                     bound = bound + 1
                     if (result) continue_radial_search = .true.
                     if (result) found_neighbors = .true.
                  enddo
                  ! top boundary 
                  bound_start = min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_left
                  bound_end =  min(j0+n-1,map_src_ysize-1)*map_src_xsize + i_right
                  if (bound_start > bound_end) then
                     bound_start = min(j0+n-1,map_src_ysize-1)*map_src_xsize + 1
                     bound_end = min(j0+n,map_src_ysize-1)*map_src_xsize
                  endif
                  bound = bound_start
                  if (bound_end > map_src_size) call mpp_error(FATAL)
                  do while (bound <= bound_end)
                     d = spherical_distance(theta_dest(j),phi_dest(j),theta_src(bound),phi_src(bound))
                     result = update_dest_neighbors(j,bound,d)
                     bound = bound + 1
                     if (result) continue_radial_search = .true.
                     if (result) found_neighbors = .true.
                  enddo
               enddo
               continue_search = .false. ! stop looking
            endif
         endif
         step=step+step_size
      enddo ! search loop
      step = 1
      step_size = step_size/2
   enddo

   if (.not.found_neighbors) then ! no neighbors found 
      lmask_dest(j) = .false.
   endif

enddo ! map dest


do j=1,map_dest_size
   sum=0.0
! neighbors are sorted nearest to farthest
! check nearest to see if it is a land point
   if (map_src_add(j,1) == 0) then
      lmask_dest(j) = .false.
   else if (.not.lmask_src(map_src_add(j,1))) then
      lmask_dest(j) = .false. 
!
! compare first 2 nearest neighbors -- if they are nearly
! equidistant then use this mask for robustness
!
      if (map_src_add(j,2) /= 0 .and. abs(map_src_dist(j,2) - map_src_dist(j,1)) < epsln) then
         if (.not.lmask_src(map_src_add(j,2))) then
            lmask_dest(j) = .false.
         endif
      endif

   endif

 
   do n=1,num_neighbors
      if (map_src_dist(j,n) <= epsln) then
          map_src_wgt(j,n) = large
          sum = sum+large
      elseif (map_src_dist(j,n) <= max_src_dist .and. map_src_dist(j,n) > epsln .and. lmask_src(map_src_add(j,n))) then
          map_src_wgt(j,n) = 1.0/map_src_dist(j,n)
          sum = sum+map_src_wgt(j,n)
      else
          map_src_wgt(j,n)=0.0
      endif
   enddo
   if (sum > epsln) then
      map_src_wgt(j,:) = map_src_wgt(j,:)/sum
   else
      lmask_dest(j) = .false.
   endif
enddo
   
deallocate(theta_dest,theta_src,phi_dest,phi_src)

return

spherical_regrid_initialized = .true.

end subroutine spherical_regrid_init_2d
!</SUBROUTINE> NAME="spherical_regrid_init"

function spherical_distance(theta1,phi1,theta2,phi2)

real, intent(in) :: theta1, phi1, theta2, phi2
real :: spherical_distance

real :: r1(3), r2(3), cross(3), s, dot, ang

! this is a simple, enough way to calculate distance on the sphere
! first, construct cartesian vectors r1 and r2
! then calculate the cross-product which is proportional to the area
! between the 2 vectors.  The angular distance is arcsin of the 
! distancealong the sphere
!
! theta is longitude and phi is latitude
!


r1(1) = cos(theta1)*cos(phi1);r1(2)=sin(theta1)*cos(phi1);r1(3)=sin(phi1)
r2(1) = cos(theta2)*cos(phi2);r2(2)=sin(theta2)*cos(phi2);r2(3)=sin(phi2)

cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

s = sqrt(cross(1)**2.+cross(2)**2.+cross(3)**2.)

s = min(s,1.0-epsln)

dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

if (dot > 0) then
    ang = asin(s)
else if (dot < 0) then
    ang = pi - asin(s)
else
    ang = pi/2.
endif

spherical_distance = abs(ang) ! in radians

return

end function spherical_distance

function update_dest_neighbors(dest_add,src_add,d)

integer, intent(in) :: dest_add, src_add
real, intent(in) :: d
logical :: update_dest_neighbors

integer :: n,m

update_dest_neighbors = .false.

if (d .le. max_src_dist) then
   NLOOP : do n=1,num_neighbors
      DIST_CHK : if (d .lt. map_src_dist(dest_add,n)) then
         if (n > 1 .and. src_add == map_src_add(dest_add,n-1)) exit NLOOP
         do m=num_neighbors,n+1,-1
            map_src_add(dest_add,m) = map_src_add(dest_add,m-1)
            map_src_dist(dest_add,m) = map_src_dist(dest_add,m-1)
         enddo
         map_src_add(dest_add,n) = src_add
         map_src_dist(dest_add,n) = d
         update_dest_neighbors = .true.
         exit NLOOP ! n loop
      endif DIST_CHK
   end do NLOOP
endif

return

end function update_dest_neighbors

!<SUBROUTINE NAME="spherical_regrid_exit">
!
!<DESCRIPTION>
! deallocate storage for spherical_regrid_mod.
!</DESCRIPTION>
subroutine spherical_regrid_exit()

deallocate(map_src_add,map_src_dist,map_src_wgt,lmask_src,&
         lmask_dest)
spherical_regrid_initialized= .false.

return

end subroutine spherical_regrid_exit
!</SUBROUTINE> NAME="spherical_regrid_exit"

end module spherical_regrid_mod

#ifdef test_
program test

use spherical_regrid_mod

implicit none

real :: x1,y1,x2,y2
real, parameter :: deg2rad = 0.01745
integer :: i

do  
   print *, 'type positions x1,y1,x2,y2'
   read(5,*,end=99,err=99) x1,y1,x2,y2
!   x1=x1*deg2rad;x2=x2*deg2rad;y1=y1*deg2rad;y2=y2*deg2rad
   print *, 'x1,y1= ',x1,y1
   print *, 'x2,y2= ',x2,y2
   print *, spherical_distance(x1,y1,x2,y2)
end do

99 continue

stop

end program test
#endif




