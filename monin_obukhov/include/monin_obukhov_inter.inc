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

!> @addtogroup monin_obukhov_inter
!> @{


pure subroutine MONIN_OBUKHOV_DIFF_(vonkarm,                &
     & ustar_min,                                     &
     & neutral, stable_option,new_mo_option,rich_crit, zeta_trans, &
     & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier)

  real(kind=FMS_MO_KIND_), intent(in)                         :: vonkarm
  real(kind=FMS_MO_KIND_), intent(in)                         :: ustar_min !< = 1.e-10
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
        k_m(:,:,k) = real(vonkarm *uss*z(:,:,k), FMS_MO_KIND_)
        k_h(:,:,k) = real(k_m(:,:,k), FMS_MO_KIND_)
     end do
  else
     do k = 1, size(z,3)
        zeta = real(- vonkarm * b_star*z(:,:,k)/(uss*uss), FMS_MO_KIND_)
        do j = 1, size(z,2)
           call monin_obukhov_derivative_m(stable_option, rich_crit, zeta_trans, &
                & ni, phi_m(:,j), zeta(:,j), mask, ier)
           call monin_obukhov_derivative_t(stable_option, new_mo_option,rich_crit, zeta_trans, &
                & ni, phi_h(:,j), zeta(:,j), mask, ier)
        enddo
        k_m(:,:,k) = real(vonkarm * uss*z(:,:,k)/phi_m, FMS_MO_KIND_)
        k_h(:,:,k) = real(vonkarm * uss*z(:,:,k)/phi_h, FMS_MO_KIND_)
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
  real(kind=FMS_MO_KIND_), intent(in)                  :: error    !< = 1.e-04
  real(kind=FMS_MO_KIND_), intent(in)                  :: zeta_min !< = 1.e-06
  integer, intent(in)                                  :: max_iter !< = 20
  real(kind=FMS_MO_KIND_), intent(in)                  :: small    !< = 1.e-04
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
  r_crit = real(0.95_lkind*rich_crit, FMS_MO_KIND_)  ! convergence can get slow if one is
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
     delta_b = real(grav*(pt0 - pt)/pt0,               FMS_MO_KIND_)
     rich    = real(- z*delta_b/(speed*speed + small), FMS_MO_KIND_)
     zz      = real(max(z,z0,zt,zq),                   FMS_MO_KIND_)
  elsewhere
     rich = 0.0_lkind
  end where

  if(neutral) then

     do i = 1, n
        if(mask(i)) then
           fm(i)     = real(log(zz(i)/z0(i)), FMS_MO_KIND_)
           ft(i)     = real(log(zz(i)/zt(i)), FMS_MO_KIND_)
           fq(i)     = real(log(zz(i)/zq(i)), FMS_MO_KIND_)
           us        = real(vonkarm/fm(i),    FMS_MO_KIND_)
           bs        = real(vonkarm/ft(i),    FMS_MO_KIND_)
           qs        = real(vonkarm/fq(i),    FMS_MO_KIND_)
           drag_m(i) = real(us*us,            FMS_MO_KIND_)
           drag_t(i) = real(us*bs,            FMS_MO_KIND_)
           drag_q(i) = real(us*qs,            FMS_MO_KIND_)
           u_star(i) = real(us*speed(i),      FMS_MO_KIND_)
           b_star(i) = real(bs*delta_b(i),    FMS_MO_KIND_)
        end if
     enddo

  else

     mask_1 = mask .and. rich <  r_crit
     mask_2 = mask .and. rich >= r_crit

     do i = 1, n
        if(mask_2(i)) then
           drag_m(i) = real(drag_min_mom,       FMS_MO_KIND_)
           drag_t(i) = real(drag_min_heat,      FMS_MO_KIND_)
           drag_q(i) = real(drag_min_moist,     FMS_MO_KIND_)
           us        = real(sqrt_drag_min_mom,  FMS_MO_KIND_)
           bs        = real(sqrt_drag_min_heat, FMS_MO_KIND_)
           u_star(i) = real(us*speed(i),        FMS_MO_KIND_)
           b_star(i) = real(bs*delta_b(i),      FMS_MO_KIND_)
        end if
     enddo

     call monin_obukhov_solve_zeta (error, zeta_min, max_iter, small, &
          & stable_option, new_mo_option, rich_crit, zeta_trans,      &
          & n, rich, zz, z0, zt, zq, fm, ft, fq, mask_1, ier)

     do i = 1, n
        if(mask_1(i)) then
           us        = real(max(vonkarm/fm(i), sqrt_drag_min_mom),   FMS_MO_KIND_)
           bs        = real(max(vonkarm/ft(i), sqrt_drag_min_heat),  FMS_MO_KIND_)
           qs        = real(max(vonkarm/fq(i), sqrt_drag_min_moist), FMS_MO_KIND_)
           drag_m(i) = real(us*us,                                   FMS_MO_KIND_)
           drag_t(i) = real(us*bs,                                   FMS_MO_KIND_)
           drag_q(i) = real(us*qs,                                   FMS_MO_KIND_)
           u_star(i) = real(us*speed(i),                             FMS_MO_KIND_)
           b_star(i) = real(bs*delta_b(i),                           FMS_MO_KIND_)
        endif
     enddo

  end if

end subroutine MONIN_OBUKHOV_DRAG_1D_


pure subroutine MONIN_OBUKHOV_SOLVE_ZETA_(error, zeta_min, max_iter, small,  &
     & stable_option, new_mo_option, rich_crit, zeta_trans,        & !miz
     & n, rich, z, z0, zt, zq, f_m, f_t, f_q, mask, ier)

  real(kind=FMS_MO_KIND_), intent(in)                :: error    !< = 1.e-04
  real(kind=FMS_MO_KIND_), intent(in)                :: zeta_min !< = 1.e-06
  integer,                 intent(in)                :: max_iter !< = 20
  real(kind=FMS_MO_KIND_), intent(in)                :: small    !< = 1.e-04
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

  z_z0    = real(z/z0,      FMS_MO_KIND_)
  z_zt    = real(z/zt,      FMS_MO_KIND_)
  z_zq    = real(z/zq,      FMS_MO_KIND_)
  ln_z_z0 = real(log(z_z0), FMS_MO_KIND_)
  ln_z_zt = real(log(z_zt), FMS_MO_KIND_)
  ln_z_zq = real(log(z_zq), FMS_MO_KIND_)

  corr = 0.0_lkind
  mask_1 = mask

  ! initial guess

  zeta = 0.0_lkind
  where(mask_1)
     zeta = real(rich*ln_z_z0*ln_z_z0/ln_z_zt, FMS_MO_KIND_)
  end where

  where (mask_1 .and. rich >= 0.0)
     zeta = real(zeta/(1.0 - rich/rich_crit), FMS_MO_KIND_)
  end where

  iter_loop: do iter = 1, max_iter

     where (mask_1 .and. abs(zeta).lt.zeta_min)
        zeta = 0.0_lkind
        f_m = real(ln_z_z0, FMS_MO_KIND_)
        f_t = real(ln_z_zt, FMS_MO_KIND_)
        f_q = real(ln_z_zq, FMS_MO_KIND_)
        mask_1 = .false.  ! don't do any more calculations at these pts
     end where


     zeta_0 = 0.0_lkind
     zeta_t = 0.0_lkind
     zeta_q = 0.0_lkind
     where (mask_1)
        rzeta  = real(1.0_lkind/zeta, FMS_MO_KIND_)
        zeta_0 = real(zeta/z_z0,      FMS_MO_KIND_)
        zeta_t = real(zeta/z_zt,      FMS_MO_KIND_)
        zeta_q = real(zeta/z_zq,      FMS_MO_KIND_)
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
        df_m       = real((phi_m - phi_m_0)*rzeta,                           FMS_MO_KIND_)
        df_t       = real((phi_t - phi_t_0)*rzeta,                           FMS_MO_KIND_)
        rich_1     = real(zeta*f_t/(f_m*f_m),                                FMS_MO_KIND_)
        d_rich     = real(rich_1*( rzeta +  df_t/f_t - 2.0_lkind *df_m/f_m), FMS_MO_KIND_)
        correction = real((rich - rich_1)/d_rich,                            FMS_MO_KIND_)
        corr       = real(min(abs(correction),abs(correction/zeta)),         FMS_MO_KIND_)
        ! the criterion corr < error seems to work ok, but is a bit arbitrary
        !  when zeta is small the tolerance is reduced
     end where

     max_cor = real(maxval(corr), FMS_MO_KIND_)

     if(max_cor > error) then
        mask_1 = mask_1 .and. (corr > error)
        ! change the mask so computation proceeds only on non-converged points
        where(mask_1)
           zeta = real(zeta + correction, FMS_MO_KIND_)
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
  b_stab   = real(1.0_lkind/rich_crit, FMS_MO_KIND_)

  stable   = mask .and. zeta >= 0.0_lkind
  unstable = mask .and. zeta <  0.0_lkind

!miz: modified to include new monin-obukhov option
  if (new_mo_option) then
     where (unstable)
        phi_t = real((1 - 16.0_lkind*zeta)**(-1.0_lkind/3.0_lkind), FMS_MO_KIND_)
     end where
  else
  where (unstable)
     phi_t = real((1 - 16.0_lkind*zeta)**(-0.5_lkind), FMS_MO_KIND_)
  end where
  end if
!miz

  if(stable_option == 1) then

     where (stable)
        phi_t = real(1.0_lkind + zeta*(5.0_lkind + b_stab*zeta)/(1.0_lkind + zeta), FMS_MO_KIND_)
     end where

  else if(stable_option == 2) then

     lambda = real(1.0_lkind + (5.0_lkind - b_stab)*zeta_trans, FMS_MO_KIND_)

     where (stable .and. zeta < zeta_trans)
        phi_t = real(1 + 5.0_lkind*zeta, FMS_MO_KIND_)
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_t = real(lambda + b_stab*zeta, FMS_MO_KIND_)
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
  b_stab     = real(1.0_lkind/rich_crit, FMS_MO_KIND_)

  stable   = mask .and. zeta >= 0.0
  unstable = mask .and. zeta <  0.0

  where (unstable)
     x     = real((1 - 16.0_lkind*zeta  )**(-0.5_lkind), FMS_MO_KIND_)
     phi_m = real(sqrt(x),                               FMS_MO_KIND_) ! phi_m = (1 - 16.0*zeta)**(-0.25)
  end where

  if(stable_option == 1) then

     where (stable)
        phi_m = real(1.0_lkind + zeta  *(5.0_lkind + b_stab*zeta)/(1.0 + zeta), FMS_MO_KIND_)
     end where

  else if(stable_option == 2) then

     lambda = real(1.0_lkind + (5.0_lkind - b_stab)*zeta_trans, FMS_MO_KIND_)

     where (stable .and. zeta < zeta_trans)
        phi_m = real(1 + 5.0_lkind*zeta, FMS_MO_KIND_)
     end where
     where (stable .and. zeta >= zeta_trans)
        phi_m = real(lambda + b_stab*zeta, FMS_MO_KIND_)
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
     ln_z_z0     = real(log(z/z0),     FMS_MO_KIND_)
     ln_z_zt     = real(log(z/zt),     FMS_MO_KIND_)
     ln_z_zq     = real(log(z/zq),     FMS_MO_KIND_)
     ln_z_zref   = real(log(z/zref),   FMS_MO_KIND_)
     ln_z_zref_t = real(log(z/zref_t), FMS_MO_KIND_)
  endwhere

  if(neutral) then

     where(mask)
        del_m = real(1.0_lkind - ln_z_zref  /ln_z_z0, FMS_MO_KIND_)
        del_t = real(1.0_lkind - ln_z_zref_t/ln_z_zt, FMS_MO_KIND_)
        del_q = real(1.0_lkind - ln_z_zref_t/ln_z_zq, FMS_MO_KIND_)
     endwhere

  else

     where(mask .and. u_star > 0.0_lkind)
        mo_length_inv = real(- vonkarm * b_star/(u_star*u_star), FMS_MO_KIND_)
        zeta          = real(z     *mo_length_inv, FMS_MO_KIND_)
        zeta_0        = real(z0    *mo_length_inv, FMS_MO_KIND_)
        zeta_t        = real(zt    *mo_length_inv, FMS_MO_KIND_)
        zeta_q        = real(zq    *mo_length_inv, FMS_MO_KIND_)
        zeta_ref      = real(zref  *mo_length_inv, FMS_MO_KIND_)
        zeta_ref_t    = real(zref_t*mo_length_inv, FMS_MO_KIND_)
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
        del_m = real(1.0_lkind - f_m_ref/f_m, FMS_MO_KIND_)
        del_t = real(1.0_lkind - f_t_ref/f_t, FMS_MO_KIND_)
        del_q = real(1.0_lkind - f_q_ref/f_q, FMS_MO_KIND_)
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

  b_stab   = real(1.0_lkind/rich_crit, FMS_MO_KIND_)

  stable   = mask .and. zeta >= 0.0_lkind
  unstable = mask .and. zeta <  0.0_lkind

  where(unstable)

     x      = real(sqrt(1 - 16.0_lkind*zeta),       FMS_MO_KIND_)
     x_0    = real(sqrt(1 - 16.0_lkind*zeta_0),     FMS_MO_KIND_)

     x      = real(sqrt(x),                         FMS_MO_KIND_)
     x_0    = real(sqrt(x_0),                       FMS_MO_KIND_)

     x1     = real(1.0_lkind + x,                   FMS_MO_KIND_)
     x1_0   = real(1.0_lkind + x_0,                 FMS_MO_KIND_)

     num    = real(x1*x1*(1.0_lkind + x*x),         FMS_MO_KIND_)
     denom  = real(x1_0*x1_0*(1.0_lkind + x_0*x_0), FMS_MO_KIND_)
     y      = real(atan(x) - atan(x_0),             FMS_MO_KIND_)
     psi_m  = real(ln_z_z0 - log(num/denom) + 2*y,  FMS_MO_KIND_)

  end where

  if( stable_option == 1) then

     where (stable)
        psi_m = real(ln_z_z0 + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_0)) &
             + b_stab*(zeta - zeta_0), FMS_MO_KIND_)
     end where

  else if (stable_option == 2) then

     lambda = real(1.0_lkind + (5.0_lkind - b_stab)*zeta_trans, FMS_MO_KIND_)

     weakly_stable   = stable .and. zeta <= zeta_trans
     strongly_stable = stable .and. zeta >  zeta_trans

     where (weakly_stable)
        psi_m = real(ln_z_z0 + 5.0_lkind*(zeta - zeta_0), FMS_MO_KIND_)
     end where

     where(strongly_stable)
        x = real((lambda - 1.0_lkind)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans), FMS_MO_KIND_)
     endwhere

     where (strongly_stable .and. zeta_0 <= zeta_trans)
        psi_m = real(ln_z_z0 + x + 5.0_lkind*(zeta_trans - zeta_0), FMS_MO_KIND_)
     end where
     where (strongly_stable .and. zeta_0 > zeta_trans)
        psi_m = real(lambda*ln_z_z0 + b_stab*(zeta  - zeta_0), FMS_MO_KIND_)
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

  b_stab     = real(1.0_lkind/rich_crit, FMS_MO_KIND_)

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

!miz: modified to include a new monin-obukhov option
if (new_mo_option) then
 s3 = sqrt(3.0_lkind)
 where(unstable)
  x     = real((1 - 16.0_lkind*zeta)**(1.0_lkind/3.0_lkind),   FMS_MO_KIND_)
  x_t   = real((1 - 16.0_lkind*zeta_t)**(1.0_lkind/3.0_lkind), FMS_MO_KIND_)
  x_q   = real((1 - 16.0*zeta_q)**(1.0_lkind/3.0_lkind),       FMS_MO_KIND_)

  psi_t = real(ln_z_zt - 1.5_lkind*log((x**2+x+1)/(x_t**2 + x_t + 1)) + s3*(atan((2*x+1)/s3) - &
               atan((2*x_t + 1)/s3)),                          FMS_MO_KIND_)
  psi_q = real(ln_z_zq - 1.5_lkind*log((x**2+x+1)/(x_q**2 + x_q + 1)) + s3*(atan((2*x+1)/s3) - &
              atan((2*x_q + 1)/s3)),                           FMS_MO_KIND_)
  end where
else

where(unstable)

  x     = real(sqrt(1 - 16.0_lkind*zeta),                                    FMS_MO_KIND_)
  x_t   = real(sqrt(1 - 16.0_lkind*zeta_t),                                  FMS_MO_KIND_)
  x_q   = real(sqrt(1 - 16.0_lkind*zeta_q),                                  FMS_MO_KIND_)

  psi_t = real(ln_z_zt - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_t) ), FMS_MO_KIND_)
  psi_q = real(ln_z_zq - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_q) ), FMS_MO_KIND_)

end where
end if
!miz

if( stable_option == 1) then

  where (stable)

    psi_t = real(ln_z_zt + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_t)) &
       + b_stab*(zeta - zeta_t), FMS_MO_KIND_)
    psi_q = real(ln_z_zq + (5.0_lkind - b_stab)*log((1.0_lkind + zeta)/(1.0_lkind + zeta_q)) &
       + b_stab*(zeta - zeta_q), FMS_MO_KIND_)

  end where

else if (stable_option == 2) then

   lambda = real(1.0_lkind + (5.0_lkind - b_stab)*zeta_trans, FMS_MO_KIND_)

  weakly_stable   = stable .and. zeta <= zeta_trans
  strongly_stable = stable .and. zeta >  zeta_trans

  where (weakly_stable)
    psi_t = real(ln_z_zt + 5.0_lkind*(zeta - zeta_t), FMS_MO_KIND_)
    psi_q = real(ln_z_zq + 5.0_lkind*(zeta - zeta_q), FMS_MO_KIND_)
  end where

  where(strongly_stable)
    x = real((lambda - 1.0_lkind)*log(zeta/zeta_trans) + b_stab*(zeta - zeta_trans), FMS_MO_KIND_)
  endwhere

  where (strongly_stable .and. zeta_t <= zeta_trans)
    psi_t = real(ln_z_zt + x + 5.0_lkind*(zeta_trans - zeta_t), FMS_MO_KIND_)
  end where
  where (strongly_stable .and. zeta_t > zeta_trans)
    psi_t = real(lambda*ln_z_zt + b_stab*(zeta  - zeta_t),      FMS_MO_KIND_)
  endwhere

  where (strongly_stable .and. zeta_q <= zeta_trans)
    psi_q = real(ln_z_zq + x + 5.0_lkind*(zeta_trans - zeta_q), FMS_MO_KIND_)
  end where
  where (strongly_stable .and. zeta_q > zeta_trans)
    psi_q = real(lambda*ln_z_zq + b_stab*(zeta  - zeta_q),      FMS_MO_KIND_)
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
b_stab     = real(1.0_lkind/rich_crit,                           FMS_MO_KIND_)
rich_trans = real(zeta_trans/(1.0_lkind + 5.0_lkind*zeta_trans), FMS_MO_KIND_)

if(stable_option == 1) then

     c = - 1.0_lkind
     do i = 1, n
        if(rich(i) > 0.0_lkind .and. rich(i) < rich_crit) then
           r      = real(1.0_lkind/rich(i),                                                      FMS_MO_KIND_)
           a      = real(r - b_stab,                                                             FMS_MO_KIND_)
           b      = real(r - (1.0_lkind + 5.0_lkind),                                            FMS_MO_KIND_)
           zeta   = real((-b + sqrt(b*b - 4.0_lkind*a*c))/(2.0_lkind*a),                         FMS_MO_KIND_)
           phi    = real(1.0_lkind + b_stab*zeta + (5.0_lkind - b_stab)*zeta/(1.0_lkind + zeta), FMS_MO_KIND_)
           mix(i) = real(1.0_lkind/(phi*phi),                                                    FMS_MO_KIND_)
     endif
  end do

else if(stable_option == 2) then

  lambda = real(1.0_lkind + (5.0_lkind - b_stab)*zeta_trans, FMS_MO_KIND_)

  where(rich > 0.0_lkind .and. rich <= rich_trans)
    mix = real((1.0_lkind - 5.0_lkind*rich)**2,              FMS_MO_KIND_)
  end where
  where(rich > rich_trans .and. rich < rich_crit)
    mix = real(((1.0_lkind - b_stab*rich)/lambda)**2,        FMS_MO_KIND_)
  end where

end if

end subroutine MONIN_OBUKHOV_STABLE_MIX_

!> @}
! close documentation grouping
