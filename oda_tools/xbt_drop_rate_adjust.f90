module xbt_adjust

  use oda_types_mod, only : ocean_profile_type, missing_value

  implicit none
  
  real, parameter :: s1=0.30731408,s2=6.707e-9,s3=-8.1899e-5,sa=3.227,&
                     sb=-2.17e-4
  real, parameter :: t1=0.29585798,t2=1.002e-9,t3=-3.1658e-5,ta=3.426,&
       tb=-4.7e-4

contains
  
subroutine xbt_drop_rate_adjust(station)


  type(ocean_profile_type), intent(inout) :: station


  integer :: k
  real :: dpth_orig, dpth_new,tdrop
  integer :: fix_depth
  
  if (.not. station%accepted) return

  fix_depth = int(station%fix_depth)
  select case(fix_depth)
  case(-1)
      return
  case(0)
      return      
  case(1)
! use Hanawa et al (1994) drop rate correction
      do k=1,station%levels
         dpth_orig = station%depth(k)
         if (dpth_orig .ne. missing_value) then
             dpth_new = (1.0417*dpth_orig) - (75.096*(1.0-((1.0-(0.0002063*dpth_orig)))**0.5))
             station%depth(k)=dpth_new
         endif
      enddo
      station%fix_depth=-1.0
  case (2)
! use Kizu et al (2005) correction
      do k=1,station%levels
         dpth_orig = station%depth(k)
         if (dpth_orig .ne. missing_value) then
             if (dpth_orig .le. 250.0) then
                 dpth_new = dpth_orig*0.9572
             else if (dpth_orig .le. 500.) then
                 dpth_new = dpth_orig*0.9565
             else if (dpth_orig .le. 750.0) then
                 dpth_new = dpth_orig*0.9558
             else if (dpth_orig .le. 1000.) then
                 dpth_new = dpth_orig*0.9550
             else if (dpth_orig .le. 1250.0) then
                 dpth_new = dpth_orig*0.9542
             else if (dpth_orig .le. 1500.0) then
                 dpth_new = dpth_orig*0.9533
             else
                 dpth_new = dpth_orig*0.9524
             endif
             station%depth(k)=dpth_new       
         endif
      enddo
      station%fix_depth=-1.0      
  case(103)
      do k=1,station%levels      
         dpth_orig = station%depth(k)
         if (dpth_orig .ne. missing_value) then
             tdrop=(s1*dpth_orig + s2) - s3
             dpth_new = sa*tdrop + sb*tdrop*tdrop
             station%depth(k)=dpth_new             
         endif
      enddo
      station%fix_depth=-1.0      
  case(104)
      do k=1,station%levels
         dpth_orig = station%depth(k)
         if (dpth_orig .ne. missing_value) then
             tdrop=(t1*dpth_orig + t2) - t3
             dpth_new = ta*tdrop + tb*tdrop*tdrop
             station%depth(k)=dpth_new             
         endif
      enddo
      station%fix_depth=-1.0      
  end select

  return
end subroutine xbt_drop_rate_adjust

end module xbt_adjust
