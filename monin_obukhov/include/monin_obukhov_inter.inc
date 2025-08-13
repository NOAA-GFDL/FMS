!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @addtogroup monin_obukhov_inter
!> @{


pure subroutine MONIN_OBUKHOV_DIFF_(vonkarm,                &
     & ustar_min,                                     &
     & neutral, stable_option,new_mo_option,rich_crit, zeta_trans, &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

  real(kind=FMS_MO_KIND_), intent(in)                         :: vonkarm
  real(kind=FMS_MO_KIND_), intent(in)                         :: ustar_min !< = 1.0E-10
  logical, intent(in)                                         :: neutral
  integer, intent(in)                                         :: stable_option
  logical, intent(in)                                         :: new_mo_option !miz
  real(kind=FMS_MO_KIND_), intent(in)                         :: rich_crit, zeta_trans
  integer, intent(in)                                         :: ni, nj, nk
  real(kind=FMS_MO_KIND_), intent(in), dimension(ni, nj, nk)  :: z
  real(kind=FMS_MO_KIND_), intent(in), dimension(ni, nj)      :: u_star, b_star
  real(kind=FMS_MO_KIND_), intent(out), dimension(ni, nj, nk) :: k_m, k_h
  integer, intent(out)                                        :: ier

  real(kind=FMS_MO_KIND_),             dimension(ni, nj)      :: phi_m, phi_h, zeta, uss
  integer                :: j, k
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

end subroutine MONIN_OBUKHOV_DIFF_


pure subroutine MONIN_OBUKHOV_DRAG_1D_(grav, vonkarm,               &
     & error, zeta_min, max_iter, small,                         &
     & neutral, stable_option, new_mo_option, rich_crit, zeta_trans,&
     & drag_min_heat, drag_min_moist, drag_min_mom,              &
     & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
     & drag_q, u_star, b_star, lavail, avail, ier)

  real(kind=FMS_MO_KIND_), intent(in)                  :: grav
  real(kind=FMS_MO_KIND_), intent(in)                  :: vonkarm
  real(kind=FMS_MO_KIND_), intent(in)                  :: error    !< = 1.0E-04
  real(kind=FMS_MO_KIND_), intent(in)                  :: zeta_min !< = 1.0E-06
  integer, intent(in)                                  :: max_iter !< = 20
  real(kind=FMS_MO_KIND_), intent(in)                  :: small    !< = 1.0E-04
  logical, intent(in)                                  :: neutral
  integer, intent(in)                                  :: stable_option
  logical, intent(in)                                  :: new_mo_option
  real(kind=FMS_MO_KIND_), intent(in)                  :: rich_crit, zeta_trans
  real(kind=FMS_MO_KIND_), intent(in)                  :: drag_min_heat, drag_min_moist, drag_min_mom
  integer, intent(in)                                  :: n
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)    :: pt, pt0, z, z0, zt, zq, speed
  real(kind=FMS_MO_KIND_), intent(inout), dimension(n) :: drag_m, drag_t, drag_q, u_star, b_star
  logical, intent(in)                :: lavail !< whether to use provided mask or not
  logical, intent(in), dimension(n)  :: avail  !< provided mask
  integer, intent(out)               :: ier

  real(kind=FMS_MO_KIND_), dimension(n) :: rich, fm, ft, fq, zz
  logical, dimension(n)                 :: mask, mask_1, mask_2
  real(kind=FMS_MO_KIND_), dimension(n) :: delta_b !!, us, bs, qs
  real(kind=FMS_MO_KIND_)               :: r_crit, sqrt_drag_min_heat
  real(kind=FMS_MO_KIND_)               :: sqrt_drag_min_moist, sqrt_drag_min_mom
  real(kind=FMS_MO_KIND_)               :: us, bs, qs
  integer                               :: i
  integer, parameter                    :: lkind = FMS_MO_KIND_

  ier = 0
  r_crit = 0.95_lkind*rich_crit ! convergence can get slow if one is
                           ! close to rich_crit
  sqrt_drag_min_heat = 0.0_lkind
  if(drag_min_heat.ne.0.0_lkind)  sqrt_drag_min_heat  = sqrt(drag_min_heat)
  sqrt_drag_min_moist = 0.0_lkind
  if(drag_min_moist.ne.0.0_lkind) sqrt_drag_min_moist = sqrt(drag_min_moist)
  sqrt_drag_min_mom = 0.0_lkind
  if(drag_min_mom.ne.0.0_lkind)   sqrt_drag_min_mom   = sqrt(drag_min_mom)

  mask = .true.
  if(lavail) mask = avail

  where(mask)
     delta_b = grav*(pt0 - pt)/pt0
     rich    = - z*delta_b/(speed*speed + small)
     zz      = max(z,z0,zt,zq)
  elsewhere
     rich = 0.0_lkind
  end where

  if(neutral) then

     do i = 1, n
        if(mask(i)) then
           fm(i)     = log(zz(i)/z0(i))
           ft(i)     = log(zz(i)/zt(i))
           fq(i)     = log(zz(i)/zq(i))
           us        = vonkarm/fm(i)
           bs        = vonkarm/ft(i)
           qs        = vonkarm/fq(i)
           drag_m(i) = us*us
           drag_t(i) = us*bs
           drag_q(i) = us*qs
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        end if
     enddo

  else

     mask_1 = mask .and. rich <  r_crit
     mask_2 = mask .and. rich >= r_crit

     do i = 1, n
        if(mask_2(i)) then
           drag_m(i) = drag_min_mom
           drag_t(i) = drag_min_heat
           drag_q(i) = drag_min_moist
           us        = sqrt_drag_min_mom
           bs        = sqrt_drag_min_heat
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        end if
     enddo

     call monin_obukhov_solve_zeta (error, zeta_min, max_iter, small, &
          & stable_option, new_mo_option, rich_crit, zeta_trans,      &
          & n, rich, zz, z0, zt, zq, fm, ft, fq, mask_1, ier)

     do i = 1, n
        if(mask_1(i)) then
           us        = max(vonkarm/fm(i), sqrt_drag_min_mom)
           bs        = max(vonkarm/ft(i), sqrt_drag_min_heat)
           qs        = max(vonkarm/fq(i), sqrt_drag_min_moist)
           drag_m(i) = us*us
           drag_t(i) = us*bs
           drag_q(i) = us*qs
           u_star(i) = us*speed(i)
           b_star(i) = bs*delta_b(i)
        endif
     enddo

  end if

end subroutine MONIN_OBUKHOV_DRAG_1D_


pure subroutine MONIN_OBUKHOV_SOLVE_ZETA_(error, zeta_min, max_iter, small,  &
     & stable_option, new_mo_option, rich_crit, zeta_trans,        & !miz
     & n, rich, z, z0, zt, zq, f_m, f_t, f_q, mask, ier)

  real(kind=FMS_MO_KIND_), intent(in)                :: error    !< = 1.0E-04
  real(kind=FMS_MO_KIND_), intent(in)                :: zeta_min !< = 1.0E-06
  integer,                 intent(in)                :: max_iter !< = 20
  real(kind=FMS_MO_KIND_), intent(in)                :: small    !< = 1.0E-04
  integer,                 intent(in)                :: stable_option
  logical,                 intent(in)                :: new_mo_option
  real(kind=FMS_MO_KIND_), intent(in)                :: rich_crit, zeta_trans
  integer,                 intent(in)                :: n
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)  :: rich, z, z0, zt, zq
  logical,                 intent(in), dimension(n)  :: mask
  real(kind=FMS_MO_KIND_), intent(out), dimension(n) :: f_m, f_t, f_q
  integer,                 intent(out)               :: ier

  real(kind=FMS_MO_KIND_)               :: max_cor
  integer                               :: iter
  real(kind=FMS_MO_KIND_), dimension(n) ::   &
       d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
       ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
       phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
       zeta_0, zeta_t, zeta_q, df_m, df_t
  logical, dimension(n)                 :: mask_1
  integer, parameter                    :: lkind = FMS_MO_KIND_

  ier = 0

  z_z0    = z/z0
  z_zt    = z/zt
  z_zq    = z/zq
  ln_z_z0 = log(z_z0)
  ln_z_zt = log(z_zt)
  ln_z_zq = log(z_zq)

  corr = 0.0_lkind
  mask_1 = mask

  ! initial guess

  zeta = 0.0_lkind
  where(mask_1)
     zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
  end where

  where (mask_1 .and. rich >= 0.0_lkind)
     zeta = zeta/(1.0_lkind - rich/rich_crit)
  end where

  iter_loop: do iter = 1, max_iter

     where (mask_1 .and. abs(zeta).lt.zeta_min)
        zeta = 0.0_lkind
        f_m = ln_z_z0
        f_t = ln_z_zt
        f_q = ln_z_zq
        mask_1 = .false.  ! don't do any more calculations at these pts
     end where


     zeta_0 = 0.0_lkind
     zeta_t = 0.0_lkind
     zeta_q = 0.0_lkind
     where (mask_1)
        rzeta  = 1.0_lkind/zeta
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
        df_m       = (phi_m - phi_m_0)*rzeta
        df_t       = (phi_t - phi_t_0)*rzeta
        rich_1     = zeta*f_t/(f_m*f_m)
        d_rich     = rich_1*( rzeta +  df_t/f_t - 2.0_lkind *df_m/f_m)
        correction = (rich - rich_1)/d_rich
        corr       = min(abs(correction),abs(correction/zeta))
        ! the criterion corr < error seems to work ok, but is a bit arbitrary
        !  when zeta is small the tolerance is reduced
     end where

     max_cor = maxval(corr)

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

end subroutine MONIN_OBUKHOV_SOLVE_ZETA_


!> The differential similarity function for buoyancy and tracers
! seems to be the same as monin_obukhov_derivative_m?
pure subroutine MONIN_OBUKHOV_DERIVATIVE_T_(stable_option,new_mo_option,rich_crit, zeta_trans, &
     & n, phi_t, zeta, mask, ier)

  integer,                 intent(in)                :: stable_option
  logical,                 intent(in)                :: new_mo_option !miz
  real(kind=FMS_MO_KIND_), intent(in)                :: rich_crit, zeta_trans
  integer,                 intent(in)                :: n
  real(kind=FMS_MO_KIND_), intent(out), dimension(n) :: phi_t
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)  :: zeta
  logical,                 intent(in), dimension(n)  :: mask
  integer,                 intent(out)               :: ier

  logical, dimension(n)                              :: stable, unstable
  real(kind=FMS_MO_KIND_)                            :: b_stab, lambda
  integer, parameter                                 :: lkind = FMS_MO_KIND_

  ier      = 0
  b_stab   = 1.0_lkind/rich_crit

  stable   = mask .and. zeta >= 0.0_lkind
  unstable = mask .and. zeta <  0.0_lkind

!miz: modified to include new monin-obukhov option
  if (new_mo_option) then
     where (unstable)
        phi_t = (1.0_lkind - 16.0_lkind*zeta)**(-1.0_lkind/3.0_lkind)
     end where
  else
  where (unstable)
     phi_t = (1.0_lkind - 16.0_lkind*zeta)**(-0.5_lkind)
  end where
  end if
!miz

  if(stable_option == 1) then

     where (stable)
        phi_t = 1.0_lkind + zeta*(5.0_lkind + b_stab*zeta)/(1.0_lkind + zeta)
     end where

  else if(stable_option == 2) then

     lambda = 1.0_lkind + (5.0_lkind - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_t = 1.0_lkind + 5.0_lkind*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_t = lambda + b_stab*zeta
     end where

  endif

end subroutine MONIN_OBUKHOV_DERIVATIVE_T_


! the differential similarity function for momentum
pure subroutine MONIN_OBUKHOV_DERIVATIVE_M_(stable_option, rich_crit, zeta_trans, &
     & n, phi_m, zeta, mask, ier)

  integer,                 intent(in)                :: stable_option
  real(kind=FMS_MO_KIND_), intent(in)                :: rich_crit, zeta_trans
  integer,                 intent(in)                :: n
  real(kind=FMS_MO_KIND_), intent(out), dimension(n) :: phi_m
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)  :: zeta
  logical,                 intent(in), dimension(n)  :: mask
  integer,                 intent(out)               :: ier

  logical, dimension(n)                 :: stable, unstable
  real(kind=FMS_MO_KIND_), dimension(n) :: x
  real(kind=FMS_MO_KIND_)               :: b_stab, lambda
  integer, parameter                    :: lkind = FMS_MO_KIND_

  ier        = 0
  b_stab     = 1.0_lkind/rich_crit

  stable   = mask .and. zeta >= 0.0_lkind
  unstable = mask .and. zeta <  0.0_lkind

  where (unstable)
     x     = (1.0_lkind - 16.0_lkind*zeta  )**(-0.5_lkind)
     phi_m = sqrt(x) ! phi_m = (1 - 16.0*zeta)**(-0.25)
  end where

  if(stable_option == 1) then

     where (stable)
        phi_m = 1.0_lkind + zeta  *(5.0_lkind + b_stab*zeta)/(1.0_lkind + zeta)
     end where

  else if(stable_option == 2) then

     lambda = 1.0_lkind + (5.0_lkind - b_stab)*zeta_trans

     where (stable .and. zeta < zeta_trans)
        phi_m = 1.0_lkind + 5.0_lkind*zeta
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_m = lambda + b_stab*zeta
     end where

  endif

end subroutine MONIN_OBUKHOV_DERIVATIVE_M_


pure subroutine MONIN_OBUKHOV_PROFILE_1D_(vonkarm, &
     & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &
     & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
     & del_m, del_t, del_q, lavail, avail, ier)

  real(kind=FMS_MO_KIND_), intent(in)                 :: vonkarm
  logical,                 intent(in)                 :: neutral
  integer,                 intent(in)                 :: stable_option
  logical,                 intent(in)                 :: new_mo_option
  real(kind=FMS_MO_KIND_), intent(in)                 :: rich_crit, zeta_trans
  integer,                 intent(in)                 :: n
  real(kind=FMS_MO_KIND_), intent(in)                 :: zref, zref_t
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)   :: z, z0, zt, zq, u_star, b_star, q_star
  real(kind=FMS_MO_KIND_), intent(out), dimension(n)  :: del_m, del_t, del_q
  logical,                 intent(in)                 :: lavail !< whether to use provided mask or not
  logical,                 intent(in), dimension(n)   :: avail  !< provided mask
  integer,                 intent(out)                :: ier

  real(kind=FMS_MO_KIND_), dimension(n) :: zeta, zeta_0, zeta_t, zeta_q, zeta_ref, zeta_ref_t, &
       ln_z_z0, ln_z_zt, ln_z_zq, ln_z_zref, ln_z_zref_t,  &
       f_m_ref, f_m, f_t_ref, f_t, f_q_ref, f_q,           &
       mo_length_inv
  logical, dimension(n)                 :: mask
  integer, parameter                    :: lkind = FMS_MO_KIND_

  ier = 0

  mask = .true.
  if(lavail) mask = avail

  del_m = 0.0_lkind  ! zero output arrays
  del_t = 0.0_lkind
  del_q = 0.0_lkind

  where(mask)
     ln_z_z0     = log(z/z0)
     ln_z_zt     = log(z/zt)
     ln_z_zq     = log(z/zq)
     ln_z_zref   = log(z/zref)
     ln_z_zref_t = log(z/zref_t)
  endwhere

  if(neutral) then

     where(mask)
        del_m = 1.0_lkind - ln_z_zref  /ln_z_z0
        del_t = 1.0_lkind - ln_z_zref_t/ln_z_zt
        del_q = 1.0_lkind - ln_z_zref_t/ln_z_zq
     endwhere

  else

     where(mask .and. u_star > 0.0_lkind)
        mo_length_inv = - vonkarm * b_star/(u_star*u_star)
        zeta          = z     *mo_length_inv
        zeta_0        = z0    *mo_length_inv
        zeta_t        = zt    *mo_length_inv
        zeta_q        = zq    *mo_length_inv
        zeta_ref      = zref  *mo_length_inv
        zeta_ref_t    = zref_t*mo_length_inv
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
        del_m = 1.0_lkind - f_m_ref/f_m
        del_t = 1.0_lkind - f_t_ref/f_t
        del_q = 1.0_lkind - f_q_ref/f_q
     endwhere

  end if


end subroutine MONIN_OBUKHOV_PROFILE_1D_


!> The integral similarity function for momentum
pure subroutine MONIN_OBUKHOV_INTEGRAL_M_(stable_option, rich_crit, zeta_trans, &
     & n, psi_m, zeta, zeta_0, ln_z_z0, mask, ier)

  integer,                 intent(in)                   :: stable_option
  real(kind=FMS_MO_KIND_), intent(in)                   :: rich_crit, zeta_trans
  integer,                 intent(in)                   :: n
  real(kind=FMS_MO_KIND_), intent(inout), dimension(n)  :: psi_m
  real(kind=FMS_MO_KIND_), intent(in)   , dimension(n)  :: zeta, zeta_0, ln_z_z0
  logical,                 intent(in)   , dimension(n)  :: mask
  integer,                 intent(out)                  :: ier

  real(kind=FMS_MO_KIND_)                   :: b_stab, lambda
  real(kind=FMS_MO_KIND_), dimension(n)     :: x, x_0, x1, x1_0, num, denom, y
  logical, dimension(n)                     :: stable, unstable, &
                                               weakly_stable, strongly_stable
  integer, parameter                        :: lkind = FMS_MO_KIND_

  ier      = 0

  b_stab   = 1.0_lkind/rich_crit

  stable   = mask .and. zeta >= 0.0_lkind
  unstable = mask .and. zeta <  0.0_lkind

  where(unstable)

     x      = sqrt(1.0_lkind - 16.0_lkind*zeta)
     x_0    = sqrt(1.0_lkind - 16.0_lkind*zeta_0)

     x      = sqrt(x)
     x_0    = sqrt(x_0)

     x1     = 1.0_lkind + x
     x1_0   = 1.0_lkind + x_0

     num    = x1*x1*(1.0_lkind + x*x)
     denom  = x1_0*x1_0*(1.0_lkind + x_0*x_0)
     y      = atan(x) - atan(x_0)
     psi_m  = ln_z_z0 - log(num/denom) + 2.0_lkind*y

  end where

  if( stable_option == 1) then

     where (stable)
        psi_m = ln_z_z0 + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_0)) &
             + b_stab*(zeta - zeta_0)
     end where

  else if (stable_option == 2) then

     lambda = 1.0_lkind + (5.0_lkind - b_stab)*zeta_trans

     weakly_stable   = stable .and. zeta <= zeta_trans
     strongly_stable = stable .and. zeta >  zeta_trans

     where (weakly_stable)
        psi_m = ln_z_z0 + 5.0_lkind*(zeta - zeta_0)
     end where

     where(strongly_stable)
        x = (lambda - 1.0_lkind)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
     endwhere

     where (strongly_stable .and. zeta_0 <= zeta_trans)
        psi_m = ln_z_z0 + x + 5.0_lkind*(zeta_trans - zeta_0)
     end where
     where (strongly_stable .and. zeta_0 > zeta_trans)
        psi_m = lambda*ln_z_z0 + b_stab*(zeta  - zeta_0)
     endwhere

  end if

end subroutine MONIN_OBUKHOV_INTEGRAL_M_


!> The integral similarity function for moisture and tracers
pure subroutine MONIN_OBUKHOV_INTEGRAL_TQ_(stable_option, new_mo_option, rich_crit, zeta_trans, &
     & n, psi_t, psi_q, zeta, zeta_t, zeta_q, &
     & ln_z_zt, ln_z_zq, mask, ier)

  integer,                 intent(in)                  :: stable_option
  logical,                 intent(in)                  :: new_mo_option !miz
  real(kind=FMS_MO_KIND_), intent(in)                  :: rich_crit, zeta_trans
  integer,                 intent(in)                  :: n
  real(kind=FMS_MO_KIND_), intent(inout), dimension(n) :: psi_t, psi_q
  real(kind=FMS_MO_KIND_), intent(in)   , dimension(n) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
  logical,                 intent(in)   , dimension(n) :: mask
  integer,                 intent(out)                 :: ier

  real(kind=FMS_MO_KIND_), dimension(n)     :: x, x_t, x_q
  logical, dimension(n)                     :: stable, unstable, &
                                               weakly_stable, strongly_stable
  real(kind=FMS_MO_KIND_)                   :: b_stab, lambda
  real(kind=FMS_MO_KIND_)                   :: s3 !miz
  integer, parameter                        :: lkind = FMS_MO_KIND_

  ier = 0

  b_stab     = 1.0_lkind/rich_crit

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

!miz: modified to include a new monin-obukhov option
if (new_mo_option) then
 s3 = sqrt(3.0_lkind)
 where(unstable)
  x     = (1.0_lkind - 16.0_lkind*zeta)**(1.0_lkind/3.0_lkind)
  x_t   = (1.0_lkind - 16.0_lkind*zeta_t)**(1.0_lkind/3.0_lkind)
  x_q   = (1.0_lkind - 16.0_lkind*zeta_q)**(1.0_lkind/3.0_lkind)

  psi_t = ln_z_zt - 1.5_lkind*log((x**2+x+1.0_lkind)/(x_t**2 + x_t + 1.0_lkind)) + s3 * &
          (atan((2.0_lkind*x+1.0_lkind)/s3) - atan((2.0_lkind*x_t + 1.0_lkind)/s3))
  psi_q = ln_z_zq - 1.5_lkind*log((x**2+x+1.0_lkind)/(x_q**2 + x_q + 1.0_lkind)) + s3 * &
          (atan((2.0_lkind*x+1.0_lkind)/s3) - atan((2.0_lkind*x_q + 1.0_lkind)/s3))
  end where
else

where(unstable)

  x     = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_t   = sqrt(1.0_lkind - 16.0_lkind*zeta_t)
  x_q   = sqrt(1.0_lkind - 16.0_lkind*zeta_q)

  psi_t = ln_z_zt - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_t) )
  psi_q = ln_z_zq - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_q) )

end where
end if
!miz

if( stable_option == 1) then

  where (stable)

    psi_t = ln_z_zt + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_t)) &
       + b_stab*(zeta - zeta_t)
    psi_q = ln_z_zq + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_q)) &
       + b_stab*(zeta - zeta_q)

  end where

else if (stable_option == 2) then

   lambda = 1.0_lkind + (5.0_lkind - b_stab)*zeta_trans

  weakly_stable   = stable .and. zeta <= zeta_trans
  strongly_stable = stable .and. zeta >  zeta_trans

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0_lkind*(zeta - zeta_t)
    psi_q = ln_z_zq + 5.0_lkind*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (lambda - 1.0_lkind)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans)
  endwhere

  where (strongly_stable .and. zeta_t <= zeta_trans)
    psi_t = ln_z_zt + x + 5.0_lkind*(zeta_trans - zeta_t)
  end where
  where (strongly_stable .and. zeta_t > zeta_trans)
    psi_t = lambda*ln_z_zt + b_stab*(zeta  - zeta_t)
  endwhere

  where (strongly_stable .and. zeta_q <= zeta_trans)
    psi_q = ln_z_zq + x + 5.0_lkind*(zeta_trans - zeta_q)
  end where
  where (strongly_stable .and. zeta_q > zeta_trans)
    psi_q = lambda*ln_z_zq + b_stab*(zeta  - zeta_q)
  endwhere

end if

end subroutine MONIN_OBUKHOV_INTEGRAL_TQ_


pure subroutine MONIN_OBUKHOV_STABLE_MIX_(stable_option, rich_crit, zeta_trans, &
     &                              n, rich, mix, ier)

  integer,                 intent(in)                :: stable_option
  real(kind=FMS_MO_KIND_), intent(in)                :: rich_crit, zeta_trans
  integer,                 intent(in)                :: n
  real(kind=FMS_MO_KIND_), intent(in), dimension(n)  :: rich
  real(kind=FMS_MO_KIND_), intent(out), dimension(n) :: mix
  integer,                 intent(out)               :: ier

  real(kind=FMS_MO_KIND_)              :: r, a, b, c, zeta, phi
  real(kind=FMS_MO_KIND_)              :: b_stab, rich_trans, lambda
  integer                              :: i
  integer, parameter                   :: lkind = FMS_MO_KIND_

  ier = 0

mix        = 0.0_lkind
b_stab     = 1.0_lkind/rich_crit
rich_trans = zeta_trans/(1.0_lkind + 5.0_lkind*zeta_trans)

if(stable_option == 1) then

     c = - 1.0_lkind
     do i = 1, n
        if(rich(i) > 0.0_lkind .and. rich(i) < rich_crit) then
           r      = 1.0_lkind/rich(i)
           a      = r - b_stab
           b      = r - (1.0_lkind + 5.0_lkind)
           zeta   = (-b + sqrt(b*b - 4.0_lkind*a*c))/(2.0_lkind*a)
           phi    = 1.0_lkind + b_stab*zeta + (5.0_lkind - b_stab)*zeta/(1.0_lkind + zeta)
           mix(i) = 1.0_lkind/(phi*phi)
     endif
  end do

else if(stable_option == 2) then

  lambda = 1.0_lkind + (5.0_lkind - b_stab)*zeta_trans

  where(rich > 0.0_lkind .and. rich <= rich_trans)
    mix = (1.0_lkind - 5.0_lkind*rich)**2
  end where
  where(rich > rich_trans .and. rich < rich_crit)
    mix = ((1.0_lkind - b_stab*rich)/lambda)**2
  end where

end if

end subroutine MONIN_OBUKHOV_STABLE_MIX_

!> @}
! close documentation grouping
