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
!> @defgroup monin_obukhov_inter monin_obukhov_inter
!> @ingroup monin_obukhov
!> @brief Utility routines to be used in @ref monin_obukhov_mod

!> @file
!> @brief File for @ref monin_obukhov_inter

!> @addtogroup monin_obukhov_inter
!> @{
module monin_obukhov_inter
implicit none
private


public :: monin_obukhov_diff
public :: monin_obukhov_drag_1d
public :: monin_obukhov_solve_zeta
public :: monin_obukhov_derivative_t
public :: monin_obukhov_derivative_m
public :: monin_obukhov_profile_1d
public :: monin_obukhov_integral_m
public :: monin_obukhov_integral_tq
public :: monin_obukhov_stable_mix


contains


pure subroutine monin_obukhov_diff(vonkarm,                &
     & ustar_min,                                     &
     & neutral, stable_option,new_mo_option,rich_crit, zeta_trans, &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

  real   , intent(in   )                        :: vonkarm
  real   , intent(in   )                        :: ustar_min !< = 1.e-10
  logical, intent(in   )                        :: neutral
  integer, intent(in   )                        :: stable_option
  logical, intent(in   )                        :: new_mo_option !miz
  real   , intent(in   )                        :: rich_crit, zeta_trans
  integer, intent(in   )                        :: ni, nj, nk
  real   , intent(in   ), dimension(ni, nj, nk) :: z
  real   , intent(in   ), dimension(ni, nj)     :: u_star, b_star
  real   , intent(  out), dimension(ni, nj, nk) :: k_m, k_h
  integer, intent(  out)                        :: ier

  real , dimension(ni, nj) :: phi_m, phi_h, zeta, uss
  integer :: j, k
  logical, dimension(ni) :: mask

  ier = 0

  mask = .true.
  uss = max(u_star, ustar_min)

  if(neutral) then
     do k = 1, size(z,3)
        k_m(:,:,k) = vonkarm *uss*z(:,:,k)
        k_h(:,:,k) = k_m(:,:,k)
     end do
  else
     do k = 1, size(z,3)
        zeta = - vonkarm * b_star*z(:,:,k)/(uss*uss)
        do j = 1, size(z,2)
           call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
                & ni, phi_m(:,j), zeta(:,j), mask, ier)
           call monin_obukhov_derivative_t(stable_option, new_mo_option,rich_crit, zeta_trans, &
                & ni, phi_h(:,j), zeta(:,j), mask, ier)
        enddo
        k_m(:,:,k) = vonkarm * uss*z(:,:,k)/phi_m
        k_h(:,:,k) = vonkarm * uss*z(:,:,k)/phi_h
     end do
  endif

end subroutine monin_obukhov_diff


pure subroutine monin_obukhov_drag_1d(grav, vonkarm,               &
     & error, zeta_min, max_iter, small,                         &
     & neutral, stable_option, new_mo_option, rich_crit, zeta_trans,&
     & drag_min_heat, drag_min_moist, drag_min_mom,              &
     & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
     & drag_q, u_star, b_star, lavail, avail, ier)

  real   , intent(in   )                :: grav
  real   , intent(in   )                :: vonkarm
  real   , intent(in   )                :: error    !< = 1.e-04
  real   , intent(in   )                :: zeta_min !< = 1.e-06
  integer, intent(in   )                :: max_iter !< = 20
  real   , intent(in   )                :: small    !< = 1.e-04
  logical, intent(in   )                :: neutral
  integer, intent(in   )                :: stable_option
  logical, intent(in   )                :: new_mo_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  real   , intent(in   )                :: drag_min_heat, drag_min_moist, drag_min_mom
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: pt, pt0, z, z0, zt, zq, speed
  real   , intent(inout), dimension(n)  :: drag_m, drag_t, drag_q, u_star, b_star
  logical, intent(in   )                :: lavail !< whether to use provided mask or not
  logical, intent(in   ), dimension(n)  :: avail  !< provided mask
  integer, intent(out  )                :: ier

  real   , dimension(n) :: rich, fm, ft, fq, zz
  logical, dimension(n) :: mask, mask_1, mask_2
  real   , dimension(n) :: delta_b !!, us, bs, qs
  real                  :: r_crit, sqrt_drag_min_heat
  real                  :: sqrt_drag_min_moist, sqrt_drag_min_mom
  real                  :: us, bs, qs
  integer               :: i

  ier = 0
  r_crit = 0.95*rich_crit  ! convergence can get slow if one is
                           ! close to rich_crit
  sqrt_drag_min_heat = 0.0
  if(drag_min_heat.ne.0.0) sqrt_drag_min_heat = sqrt(drag_min_heat)
  sqrt_drag_min_moist = 0.0
  if(drag_min_moist.ne.0.0) sqrt_drag_min_moist = sqrt(drag_min_moist)
  sqrt_drag_min_mom = 0.0
  if(drag_min_mom.ne.0.0) sqrt_drag_min_mom = sqrt(drag_min_mom)

  mask = .true.
  if(lavail) mask = avail

  where(mask)
     delta_b = grav*(pt0 - pt)/pt0
     rich    = - z*delta_b/(speed*speed + small)
     zz      = max(z,z0,zt,zq)
  elsewhere
     rich = 0.0
  end where

  if(neutral) then

     do i = 1, n
        if(mask(i)) then
           fm(i)   = log(zz(i)/z0(i))
           ft(i)   = log(zz(i)/zt(i))
           fq(i)   = log(zz(i)/zq(i))
           us   = vonkarm/fm(i)
           bs   = vonkarm/ft(i)
           qs   = vonkarm/fq(i)
           drag_m(i)    = us*us
           drag_t(i)    = us*bs
           drag_q(i)    = us*qs
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        end if
     enddo

  else

     mask_1 = mask .and. rich <  r_crit
     mask_2 = mask .and. rich >= r_crit

     do i = 1, n
        if(mask_2(i)) then
           drag_m(i)   = drag_min_mom
           drag_t(i)   = drag_min_heat
           drag_q(i)   = drag_min_moist
           us       = sqrt_drag_min_mom
           bs       = sqrt_drag_min_heat
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        end if
     enddo

     call monin_obukhov_solve_zeta (error, zeta_min, max_iter, small, &
          & stable_option, new_mo_option, rich_crit, zeta_trans,      &
          & n, rich, zz, z0, zt, zq, fm, ft, fq, mask_1, ier)

     do i = 1, n
        if(mask_1(i)) then
           us   = max(vonkarm/fm(i), sqrt_drag_min_mom)
           bs   = max(vonkarm/ft(i), sqrt_drag_min_heat)
           qs   = max(vonkarm/fq(i), sqrt_drag_min_moist)
           drag_m(i)   = us*us
           drag_t(i)   = us*bs
           drag_q(i)   = us*qs
           u_star(i)   = us*speed(i)
           b_star(i)   = bs*delta_b(i)
        endif
     enddo

  end if

end subroutine monin_obukhov_drag_1d


pure subroutine monin_obukhov_solve_zeta(error, zeta_min, max_iter, small,  &
     & stable_option, new_mo_option, rich_crit, zeta_trans,        & !miz
     & n, rich, z, z0, zt, zq, f_m, f_t, f_q, mask, ier)

  real   , intent(in   )                :: error    !< = 1.e-04
  real   , intent(in   )                :: zeta_min !< = 1.e-06
  integer, intent(in   )                :: max_iter !< = 20
  real   , intent(in   )                :: small    !< = 1.e-04
  integer, intent(in   )                :: stable_option
  logical, intent(in   )                :: new_mo_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(in   ), dimension(n)  :: rich, z, z0, zt, zq
  logical, intent(in   ), dimension(n)  :: mask
  real   , intent(  out), dimension(n)  :: f_m, f_t, f_q
  integer, intent(  out)                :: ier

  real    :: max_cor
  integer :: iter
  real, dimension(n) ::   &
       d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
       ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
       phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
       zeta_0, zeta_t, zeta_q, df_m, df_t
  logical, dimension(n) :: mask_1

  ier = 0

  z_z0 = z/z0
  z_zt = z/zt
  z_zq = z/zq
  ln_z_z0 = log(z_z0)
  ln_z_zt = log(z_zt)
  ln_z_zq = log(z_zq)

  corr = 0.0
  mask_1 = mask

  ! initial guess

  zeta = 0.0
  where(mask_1)
     zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
  end where

  where (mask_1 .and. rich >= 0.0)
     zeta = zeta/(1.0 - rich/rich_crit)
  end where

  iter_loop: do iter = 1, max_iter

     where (mask_1 .and. abs(zeta).lt.zeta_min)
        zeta = 0.0
        f_m = ln_z_z0
        f_t = ln_z_zt
        f_q = ln_z_zq
        mask_1 = .false.  ! don't do any more calculations at these pts
     end where


     zeta_0 = 0.0
     zeta_t = 0.0
     zeta_q = 0.0
     where (mask_1)
        rzeta  = 1.0/zeta
        zeta_0 = zeta/z_z0
        zeta_t = zeta/z_zt
        zeta_q = zeta/z_zq
     end where

     call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m  , zeta  , mask_1, ier)
     call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
          & n, phi_m_0, zeta_0,  mask_1, ier)
     call monin_obukhov_derivative_t(stable_option, new_mo_option,rich_crit, zeta_trans, &
          & n, phi_t  , zeta  , mask_1, ier)
     call monin_obukhov_derivative_t(stable_option, new_mo_option,rich_crit, zeta_trans, &
          & n, phi_t_0, zeta_t, mask_1, ier)

     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m, zeta, zeta_0, ln_z_z0, mask_1, ier)
     call monin_obukhov_integral_tq(stable_option, new_mo_option, rich_crit, zeta_trans, &
          & n, f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1, ier)

     where (mask_1)
        df_m  = (phi_m - phi_m_0)*rzeta
        df_t  = (phi_t - phi_t_0)*rzeta
        rich_1 = zeta*f_t/(f_m*f_m)
        d_rich = rich_1*( rzeta +  df_t/f_t - 2.0 *df_m/f_m)
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

  ier = 1 ! surface drag iteration did not converge

end subroutine monin_obukhov_solve_zeta


!> The differential similarity function for buoyancy and tracers
! seems to be the same as monin_obukhov_derivative_m?
pure subroutine monin_obukhov_derivative_t(stable_option,new_mo_option,rich_crit, zeta_trans, &
     & n, phi_t, zeta, mask, ier)

  integer, intent(in   )                :: stable_option
  logical, intent(in   )                :: new_mo_option !miz
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: phi_t
  real   , intent(in   ), dimension(n)  :: zeta
  logical, intent(in   ), dimension(n)  :: mask
  integer, intent(  out)                :: ier

  logical, dimension(n) :: stable, unstable
  real                  :: b_stab, lambda

  ier = 0
  b_stab     = 1.0/rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

!miz: modified to include new monin-obukhov option
  if (new_mo_option) then
     where (unstable)
        phi_t = (1 - 16.0*zeta)**(-1./3.)
     end where
  else
  where (unstable)
     phi_t = (1 - 16.0*zeta)**(-0.5)
  end where
  end if
!miz

  if(stable_option == 1) then

     where (stable)
        phi_t = 1.0 + zeta*(5.0 + b_stab*zeta)/(1.0 + zeta)
     end where

  else if(stable_option == 2) then

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_t = 1 + 5.0*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_t = lambda + b_stab*zeta
     end where

  endif

end subroutine monin_obukhov_derivative_t


! the differential similarity function for momentum
pure subroutine monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
     & n, phi_m, zeta, mask, ier)

  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(  out), dimension(n)  :: phi_m
  real   , intent(in   ), dimension(n)  :: zeta
  logical, intent(in   ), dimension(n)  :: mask
  integer, intent(out  )                :: ier

  logical, dimension(n) :: stable, unstable
  real   , dimension(n) :: x
  real                  :: b_stab, lambda

  ier = 0
  b_stab     = 1.0/rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     x     = (1 - 16.0*zeta  )**(-0.5)
     phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
  end where

  if(stable_option == 1) then

     where (stable)
        phi_m = 1.0 + zeta  *(5.0 + b_stab*zeta)/(1.0 + zeta)
     end where

  else if(stable_option == 2) then

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_m = 1 + 5.0*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_m = lambda + b_stab*zeta
     end where

  endif

end subroutine monin_obukhov_derivative_m


pure subroutine monin_obukhov_profile_1d(vonkarm, &
     & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &
     & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
     & del_m, del_t, del_q, lavail, avail, ier)

  real   , intent(in   )                :: vonkarm
  logical, intent(in   )                :: neutral
  integer, intent(in   )                :: stable_option
  logical, intent(in   )                :: new_mo_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real,    intent(in   )                :: zref, zref_t
  real,    intent(in   ), dimension(n)  :: z, z0, zt, zq, u_star, b_star, q_star
  real,    intent(  out), dimension(n)  :: del_m, del_t, del_q
  logical, intent(in   )                :: lavail !< whether to use provided mask or not
  logical, intent(in   ), dimension(n)  :: avail  !< provided mask
  integer, intent(out  )                :: ier

  real, dimension(n) :: zeta, zeta_0, zeta_t, zeta_q, zeta_ref, zeta_ref_t, &
       ln_z_z0, ln_z_zt, ln_z_zq, ln_z_zref, ln_z_zref_t,  &
       f_m_ref, f_m, f_t_ref, f_t, f_q_ref, f_q,           &
       mo_length_inv
  logical, dimension(n) :: mask

  ier = 0

  mask = .true.
  if(lavail) mask = avail

  del_m = 0.0  ! zero output arrays
  del_t = 0.0
  del_q = 0.0

  where(mask)
     ln_z_z0     = log(z/z0)
     ln_z_zt     = log(z/zt)
     ln_z_zq     = log(z/zq)
     ln_z_zref   = log(z/zref)
     ln_z_zref_t = log(z/zref_t)
  endwhere

  if(neutral) then

     where(mask)
        del_m = 1.0 - ln_z_zref  /ln_z_z0
        del_t = 1.0 - ln_z_zref_t/ln_z_zt
        del_q = 1.0 - ln_z_zref_t/ln_z_zq
     endwhere

  else

     where(mask .and. u_star > 0.0)
        mo_length_inv = - vonkarm * b_star/(u_star*u_star)
        zeta       = z     *mo_length_inv
        zeta_0     = z0    *mo_length_inv
        zeta_t     = zt    *mo_length_inv
        zeta_q     = zq    *mo_length_inv
        zeta_ref   = zref  *mo_length_inv
        zeta_ref_t = zref_t*mo_length_inv
     endwhere

     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m,     zeta, zeta_0,   ln_z_z0,   mask, ier)
     call monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
          & n, f_m_ref, zeta, zeta_ref, ln_z_zref, mask, ier)

     call monin_obukhov_integral_tq(stable_option, new_mo_option, rich_crit, zeta_trans, &
          & n, f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask, ier)
     call monin_obukhov_integral_tq(stable_option, new_mo_option, rich_crit, zeta_trans, &
          & n, f_t_ref, f_q_ref, zeta, zeta_ref_t, zeta_ref_t, ln_z_zref_t, ln_z_zref_t,  mask, ier)

     where(mask)
        del_m = 1.0 - f_m_ref/f_m
        del_t = 1.0 - f_t_ref/f_t
        del_q = 1.0 - f_q_ref/f_q
     endwhere

  end if


end subroutine monin_obukhov_profile_1d


!> The integral similarity function for momentum
pure subroutine monin_obukhov_integral_m(stable_option, rich_crit, zeta_trans, &
     & n, psi_m, zeta, zeta_0, ln_z_z0, mask, ier)

  integer, intent(in   )                :: stable_option
  real   , intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(inout), dimension(n)  :: psi_m
  real   , intent(in)   , dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical, intent(in)   , dimension(n)  :: mask
  integer, intent(out)                  :: ier

  real                   :: b_stab, lambda
  real, dimension(n) :: x, x_0, x1, x1_0, num, denom, y
  logical, dimension(n) :: stable, unstable, &
       weakly_stable, strongly_stable

  ier = 0

  b_stab     = 1.0/rich_crit

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where(unstable)

     x     = sqrt(1 - 16.0*zeta)
     x_0   = sqrt(1 - 16.0*zeta_0)

     x      = sqrt(x)
     x_0    = sqrt(x_0)

     x1     = 1.0 + x
     x1_0   = 1.0 + x_0

     num    = x1*x1*(1.0 + x*x)
     denom  = x1_0*x1_0*(1.0 + x_0*x_0)
     y      = atan(x) - atan(x_0)
     psi_m  = ln_z_z0 - log(num/denom) + 2*y

  end where

  if( stable_option == 1) then

     where (stable)
        psi_m = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
             + b_stab*(zeta - zeta_0)
     end where

  else if (stable_option == 2) then

     lambda = 1.0 + (5.0 - b_stab)*zeta_trans

     weakly_stable   = stable .and. zeta <= zeta_trans
     strongly_stable = stable .and. zeta >  zeta_trans

     where (weakly_stable)
        psi_m = ln_z_z0 + 5.0*(zeta - zeta_0)
     end where

     where(strongly_stable)
        x = (lambda - 1.0)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
     endwhere

     where (strongly_stable .and. zeta_0 <= zeta_trans)
        psi_m = ln_z_z0 + x + 5.0*(zeta_trans - zeta_0)
     end where
     where (strongly_stable .and. zeta_0 > zeta_trans)
        psi_m = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
     endwhere

  end if

end subroutine monin_obukhov_integral_m


!> The integral similarity function for moisture and tracers
pure subroutine monin_obukhov_integral_tq(stable_option, new_mo_option, rich_crit, zeta_trans, &
     & n, psi_t, psi_q, zeta, zeta_t, zeta_q, &
     & ln_z_zt, ln_z_zq, mask, ier)

  integer, intent(in   )                :: stable_option
  logical, intent(in   )                :: new_mo_option !miz
  real,    intent(in   )                :: rich_crit, zeta_trans
  integer, intent(in   )                :: n
  real   , intent(inout), dimension(n)  :: psi_t, psi_q
  real   , intent(in)   , dimension(n)  :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  logical, intent(in)   , dimension(n)  :: mask
  integer, intent(  out)                :: ier

  real, dimension(n)     :: x, x_t, x_q
  logical, dimension(n)  :: stable, unstable, &
                                  weakly_stable, strongly_stable
  real                   :: b_stab, lambda
  real                   :: s3 !miz

  ier = 0

  b_stab     = 1.0/rich_crit

stable   = mask .and. zeta >= 0.0
unstable = mask .and. zeta <  0.0

!miz: modified to include a new monin-obukhov option
if (new_mo_option) then
 s3 = sqrt(3.0)
 where(unstable)
  x     = (1 - 16.0*zeta)**(1./3.)
  x_t   = (1 - 16.0*zeta_t)**(1./3.)
  x_q   = (1 - 16.0*zeta_q)**(1./3.)

  psi_t = ln_z_zt - 1.5*log((x**2+x+1)/(x_t**2 + x_t + 1)) + s3*(atan((2*x+1)/s3) - atan((2*x_t + 1)/s3))
  psi_q = ln_z_zq - 1.5*log((x**2+x+1)/(x_q**2 + x_q + 1)) + s3*(atan((2*x+1)/s3) - atan((2*x_q + 1)/s3))
  end where
else

where(unstable)

  x     = sqrt(1 - 16.0*zeta)
  x_t   = sqrt(1 - 16.0*zeta_t)
  x_q   = sqrt(1 - 16.0*zeta_q)

  psi_t = ln_z_zt - 2.0*log( (1.0 + x)/(1.0 + x_t) )
  psi_q = ln_z_zq - 2.0*log( (1.0 + x)/(1.0 + x_q) )

end where
end if
!miz

if( stable_option == 1) then

  where (stable)

    psi_t = ln_z_zt + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_t)) &
       + b_stab*(zeta - zeta_t)
    psi_q = ln_z_zq + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_q)) &
       + b_stab*(zeta - zeta_q)

  end where

else if (stable_option == 2) then

   lambda = 1.0 + (5.0 - b_stab)*zeta_trans

  weakly_stable   = stable .and. zeta <= zeta_trans
  strongly_stable = stable .and. zeta >  zeta_trans

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0*(zeta - zeta_t)
    psi_q = ln_z_zq + 5.0*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (lambda - 1.0)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
  endwhere

  where (strongly_stable .and. zeta_t <= zeta_trans)
    psi_t = ln_z_zt + x + 5.0*(zeta_trans - zeta_t)
  end where
  where (strongly_stable .and. zeta_t > zeta_trans)
    psi_t = lambda*ln_z_zt + b_stab*(zeta  - zeta_t)
  endwhere

  where (strongly_stable .and. zeta_q <= zeta_trans)
    psi_q = ln_z_zq + x + 5.0*(zeta_trans - zeta_q)
  end where
  where (strongly_stable .and. zeta_q > zeta_trans)
    psi_q = lambda*ln_z_zq + b_stab*(zeta  - zeta_q)
  endwhere

end if

end subroutine monin_obukhov_integral_tq


pure subroutine monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
     &                              n, rich, mix, ier)

  integer, intent(in   )                 :: stable_option
  real   , intent(in   )                 :: rich_crit, zeta_trans
  integer, intent(in   )                 :: n
  real   , intent(in   ), dimension(n)   :: rich
  real   , intent(  out), dimension(n)   :: mix
  integer, intent(  out)                 :: ier

  real               :: r, a, b, c, zeta, phi
  real               :: b_stab, rich_trans, lambda
  integer            :: i

  ier = 0

mix = 0.0
b_stab     = 1.0/rich_crit
rich_trans = zeta_trans/(1.0 + 5.0*zeta_trans)

if(stable_option == 1) then

     c = - 1.0
     do i = 1, n
        if(rich(i) > 0.0 .and. rich(i) < rich_crit) then
           r = 1.0/rich(i)
           a = r - b_stab
           b = r - (1.0 + 5.0)
           zeta = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a)
           phi = 1.0 + b_stab*zeta + (5.0 - b_stab)*zeta/(1.0 + zeta)
           mix(i) = 1./(phi*phi)
     endif
  end do

else if(stable_option == 2) then

  lambda = 1.0 + (5.0 - b_stab)*zeta_trans

  where(rich > 0.0 .and. rich <= rich_trans)
    mix = (1.0 - 5.0*rich)**2
  end where
  where(rich > rich_trans .and. rich < rich_crit)
    mix = ((1.0 - b_stab*rich)/lambda)**2
  end where

end if

end subroutine monin_obukhov_stable_mix

end module monin_obukhov_inter
!> @}
! close documentation grouping
