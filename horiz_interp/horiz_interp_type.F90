module horiz_interp_type_mod
! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     define derived data type that contains indices and weights used for subsequent 
!      interpolations.
! </OVERVIEW>

! <DESCRIPTION>
!     define derived data type that contains indices and weights used for subsequent 
!      interpolations.
! </DESCRIPTION>


use mpp_mod, only : mpp_send, mpp_recv


implicit none
private


public :: horiz_interp_type, stats


!<PUBLICTYPE >
 type horiz_interp_type
   real,    dimension(:,:), pointer   :: faci =>NULL(), facj =>NULL()   !weights for conservative scheme
   integer, dimension(:,:), pointer   :: ilon =>NULL(), jlat =>NULL()   !indices for conservative scheme
   real,    dimension(:,:), pointer   :: area_src =>NULL()              !area of the source grid
   real,    dimension(:,:), pointer   :: area_dst =>NULL()              !area of the destination grid
   real,    dimension(:,:,:), pointer :: wti =>NULL(),wtj =>NULL()      !weights for bilinear interpolation 
   integer, dimension(:,:,:), pointer :: i_lon =>NULL(), j_lat =>NULL() !indices for bilinear interpolation 
                                                                        !and spherical regrid
   real,    dimension(:,:,:), pointer :: src_dist =>NULL()              !distance between destination grid and 
                                                                        !neighbor source grid.
   logical, dimension(:,:), pointer   :: found_neighbors =>NULL()       !indicate whether destination grid 
                                                                        !has some source grid around it.
   real                               :: max_src_dist
   integer, dimension(:,:), pointer   :: num_found => NULL()
   integer                            :: nlon_src, nlat_src !size of source grid
   integer                            :: nlon_dst, nlat_dst !size of destination grid
   integer                            :: interp_method      !interpolation method.
                                                            !=1, conservative scheme
                                                            !=2, bilinear interpolation
                                                            !=3, spherical regrid
 end type
!</PUBLICTYPE>

contains

!#######################################################################
!---This statistics is for bilinear interpolation and spherical regrid.
 subroutine stats ( dat, low, high, avg, miss, missing_value, mask )
 real,    intent(in)  :: dat(:,:)
 real,    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real, intent(in), optional :: missing_value
 real,    intent(in), optional :: mask(:,:)

 real :: dsum, npts, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2)

   dsum = 0.0
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5)
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0
   
   npts = size(dat(:,:)) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
      ! Force use of "scalar", integer pointer mpp interface
         call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo         
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      ! Force use of "scalar", integer pointer mpp interface
      call mpp_send(buffer_real(1),plen=3,to_pe=root_pe)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int(1), plen=2, to_pe=root_pe)
    endif

    return

 end subroutine stats




end module horiz_interp_type_mod
