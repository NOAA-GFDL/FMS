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
!> @defgroup monin_obukhov_mod monin_obukhov_mod
!> @ingroup monin_obukhov
!> @brief Routines for computing surface drag coefficients
!! from data at the lowest model level
!! and for computing the profile of fields
!! between the lowest model level and the ground
!! using Monin-Obukhov scaling


subroutine MO_DRAG_1D_(pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, &
     u_star, b_star, avail)

real(kind=FMS_MO_KIND_), intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
real(kind=FMS_MO_KIND_), intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
logical, optional,       intent(in),    dimension(:) :: avail

logical            :: lavail, avail_dummy(1)
integer            :: n, ier

integer,                 parameter :: max_iter = 20
real(kind=FMS_MO_KIND_), parameter :: error=1.e-04, zeta_min=1.e-06, small=1.e-04

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_drag_1d in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

n      = size(pt)
lavail = .false.
if(present(avail)) lavail = .true.


if(lavail) then
if (count(avail) .eq. 0) return
call monin_obukhov_drag_1d(real(grav,kind=FMS_MO_KIND_), real(vonkarm,kind=FMS_MO_KIND_),               &
   & error, zeta_min, max_iter, real(small,kind=FMS_MO_KIND_),                         &
   & neutral, stable_option, new_mo_option, real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &!miz
   & real(drag_min_heat,kind=FMS_MO_KIND_), real(drag_min_moist,kind=FMS_MO_KIND_), &
   real(drag_min_mom,kind=FMS_MO_KIND_),              &
   & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
   & drag_q, u_star, b_star, lavail, avail, ier)
else
call monin_obukhov_drag_1d(real(grav,kind=FMS_MO_KIND_), real(vonkarm,kind=FMS_MO_KIND_),               &
   & error, zeta_min, max_iter, real(small,kind=FMS_MO_KIND_),                         &
   & neutral, stable_option, new_mo_option, real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &!miz
   & real(drag_min_heat,kind=FMS_MO_KIND_), real(drag_min_moist,kind=FMS_MO_KIND_), &
   real(drag_min_mom,kind=FMS_MO_KIND_),             &
   & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
   & drag_q, u_star, b_star, lavail, avail_dummy, ier)
endif

end subroutine MO_DRAG_1D_


!=======================================================================

subroutine MO_PROFILE_1D_(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_t, del_q, avail)

real(kind=FMS_MO_KIND_), intent(in)                :: zref, zref_t
real(kind=FMS_MO_KIND_), intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: del_m, del_t, del_q
logical, optional,       intent(in) , dimension(:) :: avail

logical                            :: dummy_avail(1)
integer                            :: n, ier

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_profile_1d in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

n = size(z)
if(present(avail)) then

if (count(avail) .eq. 0) return

call monin_obukhov_profile_1d(real(vonkarm,kind=FMS_MO_KIND_), &
   & neutral, stable_option, new_mo_option, real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &
   & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
   & del_m, del_t, del_q, .true., avail, ier)

else

call monin_obukhov_profile_1d(real(vonkarm,kind=FMS_MO_KIND_), &
   & neutral, stable_option, new_mo_option, real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &
   & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
   & del_m, del_t, del_q, .false., dummy_avail, ier)

endif

end subroutine MO_PROFILE_1D_

!=======================================================================

subroutine STABLE_MIX_3D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:,:,:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:,:)  :: mix

integer :: n, m, i, j, k, ier
real(kind=FMS_MO_KIND_), dimension(:), allocatable :: rich1d
real(kind=FMS_MO_KIND_), dimension(:), allocatable :: mix1d

if(.not.module_is_initialized) call error_mesg('stable_mix_3d in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

n = size(rich,1)*size(rich,2)*size(rich,3)

allocate(rich1d(n), mix1d(n))
rich1d = 0.0
mix1d = 0.0

do k = 1, size(rich,3)
        do j = 1, size(rich,2)
                do i = 1, size(rich,1)
                        rich1d(m) = rich(i, j, k)
                        mix1d(m) = mix(i, j, k)
                        m = m+1
                end do
        end do 
end do

call monin_obukhov_stable_mix(stable_option, real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &
& n, rich1d, mix1d, ier)

end subroutine STABLE_MIX_3D_

!=======================================================================

subroutine MO_DIFF_2D_N_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:,:) :: z
real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:)   :: u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:,:) :: k_m, k_h

integer                               :: ni, nj, nk, ier
real(kind=FMS_MO_KIND_), parameter    :: ustar_min = 1.e-10

if(.not.module_is_initialized) call error_mesg('mo_diff_2d_n in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(real(vonkarm,kind=FMS_MO_KIND_),                     &
     & ustar_min,                                     &
     & neutral, stable_option, new_mo_option,real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

end subroutine MO_DIFF_2D_N_

!=======================================================================
! The following routines are used by the public interfaces above
!=======================================================================

subroutine SOLVE_ZETA_(rich, z, z0, zt, zq, f_m, f_t, f_q, mask)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:) :: rich, z, z0, zt, zq
logical,                 intent(in) , dimension(:) :: mask
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: f_m, f_t, f_q


real(kind=FMS_MO_KIND_), parameter :: error    = 1.e-04
real(kind=FMS_MO_KIND_), parameter :: zeta_min = 1.e-06
integer,                 parameter :: max_iter = 20

real(kind=FMS_MO_KIND_)    :: max_cor
integer                    :: iter

real(kind=FMS_MO_KIND_), dimension(size(rich(:))) ::   &
     d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
     ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
     phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
     zeta_0, zeta_t, zeta_q, df_m, df_t

logical,                 dimension(size(rich(:))) :: mask_1


z_z0 = z/z0
z_zt = z/zt
z_zq = z/zq
ln_z_z0 = log(z_z0)
ln_z_zt = log(z_zt)
ln_z_zq = log(z_zq)

corr = real(0.0,kind=FMS_MO_KIND_)
mask_1 = mask

! initial guess

where(mask_1)
zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
elsewhere
zeta = real(0.0,kind=FMS_MO_KIND_)
end where

where (mask_1 .and. rich >= real(0.0,kind=FMS_MO_KIND_))
zeta = zeta/(real(1.0,kind=FMS_MO_KIND_) - rich/real(rich_crit,kind=FMS_MO_KIND_))
end where

iter_loop: do iter = 1, max_iter

where (mask_1 .and. abs(zeta).lt.zeta_min)
zeta = real(0.0,kind=FMS_MO_KIND_)
f_m = ln_z_z0
f_t = ln_z_zt
f_q = ln_z_zq
mask_1 = .false.  ! don't do any more calculations at these pts
end where

where (mask_1)
rzeta  = real(1.0,kind=FMS_MO_KIND_)/zeta
zeta_0 = zeta/z_z0
zeta_t = zeta/z_zt
zeta_q = zeta/z_zq
elsewhere
zeta_0 = real(0.0,kind=FMS_MO_KIND_)
zeta_t = real(0.0,kind=FMS_MO_KIND_)
zeta_q = real(0.0,kind=FMS_MO_KIND_)
end where

call mo_derivative_m(phi_m  , zeta  , mask_1)
call mo_derivative_m(phi_m_0, zeta_0, mask_1)
call mo_derivative_t(phi_t  , zeta  , mask_1)
call mo_derivative_t(phi_t_0, zeta_t, mask_1)

call mo_integral(f_m, zeta, zeta_0, ln_z_z0, mask_1)
call mo_integral(f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1)

where (mask_1)
df_m  = (phi_m - phi_m_0)*rzeta
df_t  = (phi_t - phi_t_0)*rzeta
rich_1 = zeta*f_t/(f_m*f_m)
d_rich = rich_1*( rzeta +  df_t/f_t - real(2.0,kind=FMS_MO_KIND_) *df_m/f_m)
correction = (rich - rich_1)/d_rich
corr = min(abs(correction),abs(correction/zeta))
 ! the criterion corr < error seems to work ok, but is a bit arbitrary
 !  when zeta is small the tolerance is reduced
end where

max_cor= maxval(corr)

if(max_cor > error) then
mask_1 = mask_1 .and. (corr > error)
  ! change the mask so computation proceeds only on non-converged points
where(mask_1)
 zeta = zeta + correction
end where
cycle iter_loop
else
return
end if

end do iter_loop

call error_mesg ('solve_zeta in monin_obukhov_mod',  &
            'surface drag iteration did not converge', FATAL)

end subroutine SOLVE_ZETA_

!=======================================================================

subroutine MO_DERIVATIVE_M_(phi_m, zeta, mask)

! the differential similarity function for momentum

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: phi_m
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta
logical,                 intent(in),  dimension(:) :: mask

logical,                 dimension(size(zeta(:))) :: stable, unstable
real(kind=FMS_MO_KIND_), dimension(size(zeta(:))) :: x

stable   = mask .and. zeta >= real(0.0,kind=FMS_MO_KIND_)
unstable = mask .and. zeta <  real(0.0,kind=FMS_MO_KIND_)

where (unstable)
x     = (1 - real(16.0,kind=FMS_MO_KIND_)*zeta  )**(real(-0.5,kind=FMS_MO_KIND_))
phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
end where

if(stable_option == 1) then

where (stable)
phi_m = real(1.0,kind=FMS_MO_KIND_) + zeta * (real(5.0,kind=FMS_MO_KIND_) &
       + real(b_stab,kind=FMS_MO_KIND_)*zeta)/(real(1.0,kind=FMS_MO_KIND_) + zeta)
end where

else if(stable_option == 2) then

where (stable .and. zeta < real(zeta_trans,kind=FMS_MO_KIND_))
phi_m = 1 + real(5.0,kind=FMS_MO_KIND_)*zeta
end where
where (stable .and. zeta >= real(zeta_trans,kind=FMS_MO_KIND_))
phi_m = real(lambda,kind=FMS_MO_KIND_) + real(b_stab,kind=FMS_MO_KIND_)*zeta
end where

endif

return
end subroutine MO_DERIVATIVE_M_

!=======================================================================

subroutine MO_DERIVATIVE_T_(phi_t, zeta, mask)

! the differential similarity function for buoyancy and tracers

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: phi_t
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta
logical,                 intent(in),  dimension(:) :: mask

logical, dimension(size(zeta(:))) :: stable, unstable

stable   = mask .and. zeta >= real(0.0,kind=FMS_MO_KIND_)
unstable = mask .and. zeta <  real(0.0,kind=FMS_MO_KIND_)

where (unstable)
phi_t = (1 - real(16.0,kind=FMS_MO_KIND_)*zeta)**(real(-0.5,kind=FMS_MO_KIND_))
end where

if(stable_option == 1) then

where (stable)
phi_t = real(1.0,kind=FMS_MO_KIND_) + zeta*(real(5.0,kind=FMS_MO_KIND_) &
        + real(b_stab,kind=FMS_MO_KIND_)*zeta)/(real(1.0,kind=FMS_MO_KIND_) + zeta)
end where

else if(stable_option == 2) then

where (stable .and. zeta < real(zeta_trans,kind=FMS_MO_KIND_))
phi_t = 1 + real(5.0,kind=FMS_MO_KIND_)*zeta
end where
where (stable .and. zeta >= real(zeta_trans,kind=FMS_MO_KIND_))
phi_t = real(lambda,kind=FMS_MO_KIND_) + real(b_stab,kind=FMS_MO_KIND_)*zeta
end where

endif

return
end subroutine MO_DERIVATIVE_T_

!=======================================================================

subroutine MO_INTEGRAL_TQ_(psi_t, psi_q, zeta, zeta_t, zeta_q, &
                      ln_z_zt, ln_z_zq, mask)

! the integral similarity function for moisture and tracers

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: psi_t, psi_q
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
logical,                 intent(in),  dimension(:) :: mask

real(kind=FMS_MO_KIND_), dimension(size(zeta(:))) :: x, x_t, x_q

logical,                 dimension(size(zeta(:))) :: stable, unstable, &
                                             weakly_stable, strongly_stable

stable   = mask .and. zeta >= real(0.0,kind=FMS_MO_KIND_)
unstable = mask .and. zeta <  real(0.0,kind=FMS_MO_KIND_)

where(unstable)

x     = sqrt(1 - real(16.0,kind=FMS_MO_KIND_)*zeta)
x_t   = sqrt(1 - real(16.0,kind=FMS_MO_KIND_)*zeta_t)
x_q   = sqrt(1 - real(16.0,kind=FMS_MO_KIND_)*zeta_q)

psi_t = ln_z_zt - real(2.0,kind=FMS_MO_KIND_)*log( (real(1.0,kind=FMS_MO_KIND_) + x) &
        / (real(1.0,kind=FMS_MO_KIND_) + x_t) )
psi_q = ln_z_zq - real(2.0,kind=FMS_MO_KIND_)*log( (real(1.0,kind=FMS_MO_KIND_) + x) &
        / (real(1.0,kind=FMS_MO_KIND_) + x_q) )

end where

if( stable_option == 1) then

where (stable)

psi_t = ln_z_zt + (real(5.0,kind=FMS_MO_KIND_) - real(b_stab,kind=FMS_MO_KIND_))*log((real(1.0,kind=FMS_MO_KIND_) &
        + zeta)/(real(1.0,kind=FMS_MO_KIND_) + zeta_t)) + real(b_stab,kind=FMS_MO_KIND_)*(zeta - zeta_t)

psi_q = ln_z_zq + (real(5.0,kind=FMS_MO_KIND_) - real(b_stab,kind=FMS_MO_KIND_))*log((real(1.0,kind=FMS_MO_KIND_) &
        + zeta)/(real(1.0,kind=FMS_MO_KIND_) + zeta_q)) + real(b_stab,kind=FMS_MO_KIND_)*(zeta - zeta_q)

end where

else if (stable_option == 2) then

weakly_stable   = stable .and. zeta <= real(zeta_trans,kind=FMS_MO_KIND_)
strongly_stable = stable .and. zeta >  real(zeta_trans,kind=FMS_MO_KIND_)

where (weakly_stable)
psi_t = ln_z_zt + real(5.0,kind=FMS_MO_KIND_)*(zeta - zeta_t)
psi_q = ln_z_zq + real(5.0,kind=FMS_MO_KIND_)*(zeta - zeta_q)
end where

where(strongly_stable)
x = (real(lambda,kind=FMS_MO_KIND_) - real(1.0,kind=FMS_MO_KIND_))*log(zeta/real(zeta_trans,kind=FMS_MO_KIND_)) &
     + real(b_stab,kind=FMS_MO_KIND_)*(zeta - real(zeta_trans,kind=FMS_MO_KIND_))
endwhere

where (strongly_stable .and. zeta_t <= real(zeta_trans,kind=FMS_MO_KIND_))
psi_t = ln_z_zt + x + real(5.0,kind=FMS_MO_KIND_)*(real(zeta_trans,kind=FMS_MO_KIND_) - zeta_t)
end where
where (strongly_stable .and. zeta_t > real(zeta_trans,kind=FMS_MO_KIND_))
psi_t = real(lambda,kind=FMS_MO_KIND_)*ln_z_zt + real(b_stab,kind=FMS_MO_KIND_)*(zeta  - zeta_t)
endwhere

where (strongly_stable .and. zeta_q <= real(zeta_trans,kind=FMS_MO_KIND_))
psi_q = ln_z_zq + x + real(5.0,kind=FMS_MO_KIND_)*(real(zeta_trans,kind=FMS_MO_KIND_) - zeta_q)
end where
where (strongly_stable .and. zeta_q > real(zeta_trans,kind=FMS_MO_KIND_))
psi_q = real(lambda,kind=FMS_MO_KIND_)*ln_z_zq + real(b_stab,kind=FMS_MO_KIND_)*(zeta  - zeta_q)
endwhere

end if

return
end subroutine MO_INTEGRAL_TQ_

!=======================================================================

subroutine MO_INTEGRAL_M_(psi_m, zeta, zeta_0, ln_z_z0, mask)

!  the integral similarity function for momentum

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: psi_m
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta, zeta_0, ln_z_z0
logical,                 intent(in),  dimension(:) :: mask

real(kind=FMS_MO_KIND_), dimension(size(zeta(:)))  :: x, x_0, x1, x1_0, num, denom, y

logical,                 dimension(size(zeta(:)))  :: stable, unstable, &
                             weakly_stable, strongly_stable

stable   = mask .and. zeta >= real(0.0,kind=FMS_MO_KIND_)
unstable = mask .and. zeta <  real(0.0,kind=FMS_MO_KIND_)

where(unstable)

x     = sqrt(1 - real(16.0,kind=FMS_MO_KIND_)*zeta)
x_0   = sqrt(1 - real(16.0,kind=FMS_MO_KIND_)*zeta_0)

x      = sqrt(x)
x_0    = sqrt(x_0)

x1     = real(1.0,kind=FMS_MO_KIND_) + x
x1_0   = real(1.0,kind=FMS_MO_KIND_) + x_0

num    = x1*x1*(real(1.0,kind=FMS_MO_KIND_) + x*x)
denom  = x1_0*x1_0*(real(1.0,kind=FMS_MO_KIND_) + x_0*x_0)
y      = atan(x) - atan(x_0)
psi_m  = ln_z_z0 - log(num/denom) + 2*y

end where

if( stable_option == 1) then

where (stable)
psi_m = ln_z_z0 + (real(5.0,kind=FMS_MO_KIND_) - real(b_stab,kind=FMS_MO_KIND_))*log((real(1.0,kind=FMS_MO_KIND_) &
        + zeta)/(real(1.0,kind=FMS_MO_KIND_) + zeta_0)) + real(b_stab,kind=FMS_MO_KIND_)*(zeta - zeta_0)
end where

else if (stable_option == 2) then

weakly_stable   = stable .and. zeta <= real(zeta_trans,kind=FMS_MO_KIND_)
strongly_stable = stable .and. zeta >  real(zeta_trans,kind=FMS_MO_KIND_)

where (weakly_stable)
psi_m = ln_z_z0 + real(5.0,kind=FMS_MO_KIND_)*(zeta - zeta_0)
end where

where(strongly_stable)
x = (real(lambda,kind=FMS_MO_KIND_) - real(1.0,kind=FMS_MO_KIND_))*log(zeta/real(zeta_trans,kind=FMS_MO_KIND_)) & 
    + real(b_stab,kind=FMS_MO_KIND_)*(zeta - real(zeta_trans,kind=FMS_MO_KIND_))
endwhere

where (strongly_stable .and. zeta_0 <= real(zeta_trans,kind=FMS_MO_KIND_))
psi_m = ln_z_z0 + x + real(5.0,kind=FMS_MO_KIND_)*(real(zeta_trans,kind=FMS_MO_KIND_)- zeta_0)
end where
where (strongly_stable .and. zeta_0 > real(zeta_trans,kind=FMS_MO_KIND_))
psi_m = real(lambda,kind=FMS_MO_KIND_)*ln_z_z0 + real(b_stab,kind=FMS_MO_KIND_)*(zeta  - zeta_0)
endwhere

end if

return
end subroutine MO_INTEGRAL_M_


!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================


subroutine MO_DRAG_2D_(pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=FMS_MO_KIND_), intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
real(kind=FMS_MO_KIND_), intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real(kind=FMS_MO_KIND_), intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
call mo_drag(pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), &
              speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
              u_star(:,j), b_star(:,j))
end do


return
end subroutine MO_DRAG_2D_

!=======================================================================
subroutine MO_DRAG_0D_(pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=FMS_MO_KIND_), intent(in)    :: z, speed, pt, pt0, z0, zt, zq
real(kind=FMS_MO_KIND_), intent(out)   :: drag_m, drag_t, drag_q, u_star, b_star

real(kind=FMS_MO_KIND_), dimension(1)  :: pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                 drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1

pt_1   (1) = pt
pt0_1  (1) = pt0
z_1    (1) = z
z0_1   (1) = z0
zt_1   (1) = zt
zq_1   (1) = zq
speed_1(1) = speed

call mo_drag(pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
            drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1)

drag_m = drag_m_1(1)
drag_t = drag_t_1(1)
drag_q = drag_q_1(1)
u_star = u_star_1(1)
b_star = b_star_1(1)

return
end subroutine MO_DRAG_0D_
!=======================================================================

subroutine MO_PROFILE_2D_(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_h, del_q)

real(kind=FMS_MO_KIND_), intent(in)                  :: zref, zref_t
real(kind=FMS_MO_KIND_), intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:) :: del_m, del_h, del_q

integer :: j

do j = 1, size(z,2)
call mo_profile(zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                 zq(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                 del_m(:,j), del_h (:,j), del_q (:,j))
enddo

return
end subroutine MO_PROFILE_2D_

!=======================================================================

subroutine MO_PROFILE_0D_(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_h, del_q)

real(kind=FMS_MO_KIND_), intent(in)   :: zref, zref_t
real(kind=FMS_MO_KIND_), intent(in)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_), intent(out)  :: del_m, del_h, del_q

real(kind=FMS_MO_KIND_), dimension(1) :: z_1, z0_1, zt_1, zq_1, u_star_1, b_star_1, q_star_1, &
                                          del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile(zref, zref_t, z_1, z0_1, zt_1, zq_1, &
               u_star_1, b_star_1, q_star_1,        &
               del_m_1, del_h_1, del_q_1)

del_m = del_m_1(1)
del_h = del_h_1(1)
del_q = del_q_1(1)


return
end subroutine MO_PROFILE_0D_

!=======================================================================

subroutine MO_PROFILE_1D_N_(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_t, del_q, avail)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:)   :: zref
real(kind=FMS_MO_KIND_), intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:) :: del_m, del_t, del_q
logical, optional,       intent(in) , dimension(:)   :: avail

integer :: k

do k = 1, size(zref(:))
if(present(avail)) then
call mo_profile(zref(k), zref(k), z, z0, zt, zq, &
  u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
else
 call mo_profile(zref(k), zref(k), z, z0, zt, zq, &
  u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k))
endif
enddo

return
end subroutine MO_PROFILE_1D_N_

!=======================================================================

subroutine MO_PROFILE_0D_N_(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_t, del_q)

real(kind=FMS_MO_KIND_),    intent(in),  dimension(:) :: zref
real(kind=FMS_MO_KIND_),    intent(in)                :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_),    intent(out), dimension(:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
call mo_profile(zref(k), zref(k), z, z0, zt, zq, &
  u_star, b_star, q_star, del_m(k), del_t(k), del_q(k))
enddo

return
end subroutine MO_PROFILE_0D_N_

!=======================================================================

subroutine MO_PROFILE_2D_N_(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                    del_m, del_t, del_q)

real(kind=FMS_MO_KIND_),    intent(in),  dimension(:)     :: zref
real(kind=FMS_MO_KIND_),    intent(in),  dimension(:,:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_),    intent(out), dimension(:,:,:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
call mo_profile(zref(k), zref(k), z, z0, zt, zq, &
  u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

return
end subroutine MO_PROFILE_2D_N_

!=======================================================================

subroutine MO_DIFF_2D_1_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:) :: z, u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:) :: k_m, k_h

real(kind=FMS_MO_KIND_), dimension(size(z,1),size(z,2),1) :: z_n, k_m_n, k_h_n

z_n(:,:,1) = z

call mo_diff(z_n, u_star, b_star, k_m_n, k_h_n)

k_m = k_m_n(:,:,1)
k_h = k_h_n(:,:,1)

return
end subroutine MO_DIFF_2D_1_


!=======================================================================

subroutine MO_DIFF_1D_1_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: z, u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: k_m, k_h

real(kind=FMS_MO_KIND_), dimension(size(z),1,1) :: z_n, k_m_n, k_h_n
real(kind=FMS_MO_KIND_), dimension(size(z),1)   :: u_star_n, b_star_n

z_n   (:,1,1) = z
u_star_n(:,1) = u_star
b_star_n(:,1) = b_star

call mo_diff(z_n, u_star_n, b_star_n, k_m_n, k_h_n)

k_m = k_m_n(:,1,1)
k_h = k_h_n(:,1,1)

return
end subroutine MO_DIFF_1D_1_

!=======================================================================

subroutine MO_DIFF_1D_N_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:) :: z
real(kind=FMS_MO_KIND_), intent(in),  dimension(:)   :: u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:) :: k_m, k_h

real(kind=FMS_MO_KIND_), dimension(size(z,1),1)            :: u_star2, b_star2
real(kind=FMS_MO_KIND_), dimension(size(z,1),1, size(z,2)) :: z2, k_m2, k_h2

integer :: n

do n = 1, size(z,2)
z2   (:,1,n) = z(:,n)
enddo
u_star2(:,1) = u_star
b_star2(:,1) = b_star

call mo_diff(z2, u_star2, b_star2, k_m2, k_h2)

do n = 1, size(z,2)
k_m(:,n) = k_m2(:,1,n)
k_h(:,n) = k_h2(:,1,n)
enddo

return
end subroutine MO_DIFF_1D_N_

!=======================================================================

subroutine MO_DIFF_0D_1_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in)  :: z, u_star, b_star
real(kind=FMS_MO_KIND_), intent(out) :: k_m, k_h

integer                                   :: ni, nj, nk, ier
real(kind=FMS_MO_KIND_), parameter        :: ustar_min = 1.e-10
real(kind=FMS_MO_KIND_), dimension(1,1,1) :: z_a, k_m_a, k_h_a
real(kind=FMS_MO_KIND_), dimension(1,1)   :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_1 in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = 1
z_a(1,1,1) = z
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm,kind=FMS_MO_KIND_), ustar_min, neutral, stable_option, new_mo_option,  &
                  real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), ni, nj, nk, z_a,                            &!miz
                  u_star_a, b_star_a, k_m_a, k_h_a, ier)
k_m = k_m_a(1,1,1)
k_h = k_h_a(1,1,1)
end subroutine MO_DIFF_0D_1_

!=======================================================================

subroutine MO_DIFF_0D_N_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: z
real(kind=FMS_MO_KIND_), intent(in)                :: u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: k_m, k_h

integer                                         :: ni, nj, nk, ier
real(kind=FMS_MO_KIND_), parameter              :: ustar_min = 1.e-10
real(kind=FMS_MO_KIND_), dimension(1,1,size(z)) :: z_a, k_m_a, k_h_a
real(kind=FMS_MO_KIND_), dimension(1,1)         :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_n in monin_obukhov_mod', &
'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = size(z(:))
z_a(1,1,:) = z(:)
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm,kind=FMS_MO_KIND_), ustar_min, neutral, stable_option,new_mo_option, &
                        real(rich_crit,kind=FMS_MO_KIND_), real(zeta_trans,kind=FMS_MO_KIND_), ni, nj, nk, z_a,                    &!miz
                        u_star_a, b_star_a, k_m_a, k_h_a, ier)
k_m(:) = k_m_a(1,1,:)
k_h(:) = k_h_a(1,1,:)
end subroutine MO_DIFF_0D_N_

!=======================================================================

subroutine STABLE_MIX_2D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:,:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:)  :: mix

real(kind=FMS_MO_KIND_), dimension(size(rich,1),size(rich,2),1) :: rich_3d, mix_3d

rich_3d(:,:,1) = rich

call stable_mix(rich_3d, mix_3d)

mix = mix_3d(:,:,1)

return
end subroutine STABLE_MIX_2D_


!=======================================================================

subroutine STABLE_MIX_1D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:)  :: mix

real(kind=FMS_MO_KIND_), dimension(size(rich),1,1)  :: rich_3d, mix_3d

rich_3d(:,1,1) = rich

call stable_mix(rich_3d, mix_3d)

mix = mix_3d(:,1,1)

return
end subroutine STABLE_MIX_1D_

!=======================================================================

subroutine STABLE_MIX_0D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in)  :: rich
real(kind=FMS_MO_KIND_), intent(out) :: mix

real(kind=FMS_MO_KIND_), dimension(1,1,1) :: rich_3d, mix_3d

rich_3d(1,1,1) = rich

call stable_mix(rich_3d, mix_3d)

mix = mix_3d(1,1,1)

return
end subroutine STABLE_MIX_0D_
!=======================================================================

!> @}
! close documentation grouping