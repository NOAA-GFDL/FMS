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


subroutine MO_DRAG_1D_ &
         (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, &
          u_star, b_star, avail)

real(kind=FMS_MO_KIND_), intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
real(kind=FMS_MO_KIND_), intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
logical, intent(in), optional, dimension(:)          :: avail

logical                                              :: lavail, avail_dummy(1)
integer                                              :: n, ier

integer, parameter                                   :: max_iter = 20
integer, parameter                                   :: lkind    = FMS_MO_KIND_
real(kind=FMS_MO_KIND_), parameter                   :: error    = 1.0E-04_lkind, &
                                                        zeta_min = 1.0E-06_lkind, &
                                                        small    = 1.0E-04_lkind

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_drag_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n      = size(pt)
lavail = .false.
if(present(avail)) lavail = .true.


if(lavail) then
   if (count(avail) .eq. 0) return
   call monin_obukhov_drag_1d(real(grav, FMS_MO_KIND_), real(vonkarm, FMS_MO_KIND_),    &
        & error, zeta_min, max_iter, real(small, FMS_MO_KIND_), neutral, stable_option, &
        & new_mo_option, real(rich_crit, FMS_MO_KIND_), real(zeta_trans, FMS_MO_KIND_), &!miz
        & real(drag_min_heat, FMS_MO_KIND_), real(drag_min_moist, FMS_MO_KIND_),        &
        & real(drag_min_mom, FMS_MO_KIND_), n, pt, pt0, z, z0, zt,                      &
        & zq, speed, drag_m, drag_t, drag_q, u_star, b_star, lavail, avail, ier)
else
call monin_obukhov_drag_1d(real(grav, FMS_MO_KIND_), real(vonkarm, FMS_MO_KIND_),             &
        & error, zeta_min, max_iter, real(small, FMS_MO_KIND_), neutral, stable_option,       &
        & new_mo_option, real(rich_crit, FMS_MO_KIND_), real(zeta_trans, FMS_MO_KIND_),       &!miz
        & real(drag_min_heat, FMS_MO_KIND_), real(drag_min_moist, FMS_MO_KIND_),              &
        & real(drag_min_mom, FMS_MO_KIND_), n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, &
        & drag_q, u_star, b_star, lavail, avail_dummy, ier)
endif

end subroutine MO_DRAG_1D_


!=======================================================================

subroutine MO_PROFILE_1D_(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real(kind=FMS_MO_KIND_),    intent(in)                :: zref, zref_t
real(kind=FMS_MO_KIND_),    intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_),    intent(out), dimension(:) :: del_m, del_t, del_q
logical, intent(in) , optional, dimension(:)          :: avail

logical                                               :: dummy_avail(1)
integer                                               :: n, ier

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_profile_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(z)
if(present(avail)) then

   if (count(avail) .eq. 0) return

   call monin_obukhov_profile_1d(real(vonkarm, FMS_MO_KIND_), &
        & neutral, stable_option, new_mo_option, real(rich_crit, FMS_MO_KIND_),   &
        & real(zeta_trans, FMS_MO_KIND_), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .true., avail, ier)

else

   call monin_obukhov_profile_1d(real(vonkarm, FMS_MO_KIND_), &
        & neutral, stable_option, new_mo_option, real(rich_crit, FMS_MO_KIND_),   &
        & real(zeta_trans, FMS_MO_KIND_), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .false., dummy_avail, ier)

endif

end subroutine MO_PROFILE_1D_

!=======================================================================

subroutine STABLE_MIX_3D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:,:,:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: n3 !< Size of dimension 3 of mix and rich
integer :: i, j !< Loop indices

n2 = size(mix, 2)
n3 = size(mix, 3)

do j=1, n3
  do i=1, n2
    call stable_mix(rich(:, i, j), mix(:, i, j))
  enddo
enddo

end subroutine STABLE_MIX_3D_

!=======================================================================

subroutine MO_DIFF_2D_N_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:,:) :: z
real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:)   :: u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:,:) :: k_m, k_h

integer                                                :: ni, nj, nk, ier
integer, parameter                                     :: lkind = FMS_MO_KIND_
real(kind=FMS_MO_KIND_), parameter                     :: ustar_min = 1.0E-10_lkind

if(.not.module_is_initialized) call error_mesg('mo_diff_2d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(real(vonkarm, FMS_MO_KIND_), ustar_min, neutral,                      &
                             & stable_option, new_mo_option, real(rich_crit, FMS_MO_KIND_),   &
                             & real(zeta_trans, FMS_MO_KIND_), ni, nj, nk, z, u_star, b_star, &
                             & k_m, k_h, ier)

end subroutine MO_DIFF_2D_N_

!=======================================================================
! The following routines are used by the public interfaces above
!=======================================================================

subroutine SOLVE_ZETA_(rich, z, z0, zt, zq, f_m, f_t, f_q, mask)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:) :: rich, z, z0, zt, zq
logical, intent(in) , dimension(:)                 :: mask
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: f_m, f_t, f_q

integer, parameter                                 :: lkind    = FMS_MO_KIND_
real(kind=FMS_MO_KIND_), parameter                 :: error    = 1.0E-04_lkind
real(kind=FMS_MO_KIND_), parameter                 :: zeta_min = 1.0E-06_lkind
integer, parameter                                 :: max_iter = 20

real(kind=FMS_MO_KIND_)                            :: max_cor
integer                                            :: iter

real(kind=FMS_MO_KIND_), dimension(size(rich(:)))  ::   &
          d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
          ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
          phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
          zeta_0, zeta_t, zeta_q, df_m, df_t

logical, dimension(size(rich(:)))                  :: mask_1

z_z0    = z/z0
z_zt    = z/zt
z_zq    = z/zq
ln_z_z0 = log(z_z0)
ln_z_zt = log(z_zt)
ln_z_zq = log(z_zq)

corr = 0.0_lkind
mask_1 = mask

! initial guess

where(mask_1)
  zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
elsewhere
  zeta = 0.0_lkind
end where

where (mask_1 .and. rich >= 0.0_lkind)
  zeta = zeta/(1.0_lkind - rich/real(rich_crit, FMS_MO_KIND_))
end where

iter_loop: do iter = 1, max_iter

  where (mask_1 .and. abs(zeta).lt.zeta_min)
    zeta   = 0.0_lkind
    f_m    = ln_z_z0
    f_t    = ln_z_zt
    f_q    = ln_z_zq
    mask_1 = .false.  ! don't do any more calculations at these pts
  end where

  where (mask_1)
    rzeta  = 1.0_lkind/zeta
    zeta_0 = zeta/z_z0
    zeta_t = zeta/z_zt
    zeta_q = zeta/z_zq
  elsewhere
    zeta_0 = 0.0_lkind
    zeta_t = 0.0_lkind
    zeta_q = 0.0_lkind
  end where

  call mo_derivative_m(phi_m  , zeta  , mask_1)
  call mo_derivative_m(phi_m_0, zeta_0, mask_1)
  call mo_derivative_t(phi_t  , zeta  , mask_1)
  call mo_derivative_t(phi_t_0, zeta_t, mask_1)

  call mo_integral_m(f_m, zeta, zeta_0, ln_z_z0, mask_1)
  call mo_integral_tq(f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1)

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

call error_mesg ('solve_zeta in monin_obukhov_mod',  &
                 'surface drag iteration did not converge', FATAL)

end subroutine SOLVE_ZETA_

!=======================================================================

subroutine MO_DERIVATIVE_M_(phi_m, zeta, mask)

! the differential similarity function for momentum

real(kind=FMS_MO_KIND_), intent(out),  dimension(:) :: phi_m
real(kind=FMS_MO_KIND_), intent(in),   dimension(:) :: zeta
logical,                 intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
real(kind=FMS_MO_KIND_), dimension(size(zeta(:)))   :: x
integer, parameter                                  :: lkind = FMS_MO_KIND_

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  x     = (1.0_lkind - 16.0_lkind*zeta  )**(-0.5_lkind)
  phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
end where

if(stable_option == 1) then

  where (stable)
    phi_m = 1.0_lkind + zeta*(5.0_lkind + real(b_stab, FMS_MO_KIND_) &
            *zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,FMS_MO_KIND_))
    phi_m = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,FMS_MO_KIND_))
    phi_m = real(lambda, FMS_MO_KIND_) + real(b_stab, FMS_MO_KIND_)*zeta
  end where

endif

return
end subroutine MO_DERIVATIVE_M_

!=======================================================================

subroutine MO_DERIVATIVE_T_(phi_t, zeta, mask)

! the differential similarity function for buoyancy and tracers

real(kind=FMS_MO_KIND_), intent(out),  dimension(:) :: phi_t
real(kind=FMS_MO_KIND_), intent(in),   dimension(:) :: zeta
logical ,                intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
integer, parameter                                  :: lkind = FMS_MO_KIND_

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  phi_t = (1.0_lkind - 16.0_lkind*zeta)**(-0.5_lkind)
end where

if(stable_option == 1) then

  where (stable)
    phi_t = 1.0_lkind + zeta * (5.0_lkind + real(b_stab, FMS_MO_KIND_)&
            * zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,FMS_MO_KIND_))
    phi_t = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,FMS_MO_KIND_))
    phi_t = real(lambda, FMS_MO_KIND_) + real(b_stab, FMS_MO_KIND_)*zeta
  end where

endif

return
end subroutine MO_DERIVATIVE_T_

!=======================================================================

subroutine MO_INTEGRAL_TQ_ (psi_t, psi_q, zeta, zeta_t, zeta_q, &
                           ln_z_zt, ln_z_zq, mask)

! the integral similarity function for moisture and tracers

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: psi_t, psi_q
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
logical , intent(in),  dimension(:)                :: mask

real(kind=FMS_MO_KIND_), dimension(size(zeta(:)))  :: x, x_t, x_q

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = FMS_MO_KIND_

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where(unstable)

  x     = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_t   = sqrt(1.0_lkind - 16.0_lkind*zeta_t)
  x_q   = sqrt(1.0_lkind - 16.0_lkind*zeta_q)

  psi_t = ln_z_zt - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_t) )
  psi_q = ln_z_zq - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_q) )

end where

if( stable_option == 1) then

  where (stable)

    psi_t = ln_z_zt + (5.0_lkind - real(b_stab, FMS_MO_KIND_)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_t)) &
            + real(b_stab, FMS_MO_KIND_)*(zeta - zeta_t)
    psi_q = ln_z_zq + (5.0_lkind - real(b_stab, FMS_MO_KIND_)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_q)) &
            + real(b_stab, FMS_MO_KIND_)*(zeta - zeta_q)

  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,FMS_MO_KIND_)
  strongly_stable = stable .and. zeta >  real(zeta_trans,FMS_MO_KIND_)

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0_lkind*(zeta - zeta_t)
    psi_q = ln_z_zq + 5.0_lkind*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (real(lambda, FMS_MO_KIND_) - 1.0_lkind)*log(zeta/real(zeta_trans, FMS_MO_KIND_)) + &
         real(b_stab, FMS_MO_KIND_)*(zeta - real(zeta_trans, FMS_MO_KIND_))
  endwhere

  where (strongly_stable .and. zeta_t <= real(zeta_trans,FMS_MO_KIND_))
    psi_t = ln_z_zt + x + 5.0_lkind * (real(zeta_trans, FMS_MO_KIND_) - zeta_t)
  end where

  where (strongly_stable .and. zeta_t > real(zeta_trans,FMS_MO_KIND_))
    psi_t = real(lambda, FMS_MO_KIND_)* ln_z_zt &
            + real(b_stab, FMS_MO_KIND_)*(zeta - zeta_t)
  endwhere

  where (strongly_stable .and. zeta_q <= real(zeta_trans,FMS_MO_KIND_))
    psi_q = ln_z_zq + x + 5.0_lkind &
            *(real(zeta_trans, FMS_MO_KIND_) - zeta_q)
  end where

  where (strongly_stable .and. zeta_q > real(zeta_trans,FMS_MO_KIND_))
    psi_q = real(lambda, FMS_MO_KIND_)*ln_z_zq + real(b_stab, FMS_MO_KIND_) &
            * (zeta - zeta_q)
  endwhere

end if

return
end subroutine MO_INTEGRAL_TQ_

!=======================================================================

subroutine MO_INTEGRAL_M_ (psi_m, zeta, zeta_0, ln_z_z0, mask)

!  the integral similarity function for momentum

real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: psi_m
real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: zeta, zeta_0, ln_z_z0
logical,                 intent(in),  dimension(:) :: mask

real(kind=FMS_MO_KIND_), dimension(size(zeta(:)))  :: x, x_0, x1, x1_0, num, denom, y

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = FMS_MO_KIND_

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
    psi_m = ln_z_z0 + (5.0_lkind - real(b_stab, FMS_MO_KIND_)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_0)) &
            + real(b_stab, FMS_MO_KIND_)*(zeta - zeta_0)
  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,FMS_MO_KIND_)
  strongly_stable = stable .and. zeta >  real(zeta_trans,FMS_MO_KIND_)

  where (weakly_stable)
    psi_m = ln_z_z0 + 5.0_lkind*(zeta - zeta_0)
  end where

  where(strongly_stable)
    x = (real(lambda, FMS_MO_KIND_) - 1.0_lkind)*log(zeta/real(zeta_trans, FMS_MO_KIND_)) + &
        real(b_stab, FMS_MO_KIND_)*(zeta - real(zeta_trans, FMS_MO_KIND_))
  endwhere

  where (strongly_stable .and. zeta_0 <= real(zeta_trans,FMS_MO_KIND_))
    psi_m = ln_z_z0 + x + 5.0_lkind &
            *(real(zeta_trans, FMS_MO_KIND_) - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > real(zeta_trans,FMS_MO_KIND_))
    psi_m = real(lambda, FMS_MO_KIND_)*ln_z_z0 + real(b_stab, FMS_MO_KIND_) &
            *(zeta - zeta_0)
  endwhere

end if

return
end subroutine MO_INTEGRAL_M_


!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================


subroutine MO_DRAG_2D_ &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=FMS_MO_KIND_), intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
real(kind=FMS_MO_KIND_), intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real(kind=FMS_MO_KIND_), intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
  call mo_drag (pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), &
                   speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
                   u_star(:,j), b_star(:,j))
end do


return
end subroutine MO_DRAG_2D_

!=======================================================================
subroutine MO_DRAG_0D_ &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

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

call mo_drag (pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
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
  call mo_profile (zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                      zq(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                      del_m(:,j), del_h (:,j), del_q (:,j))
enddo

return
end subroutine MO_PROFILE_2D_

!=======================================================================

subroutine MO_PROFILE_0D_(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real(kind=FMS_MO_KIND_), intent(in)  :: zref, zref_t
real(kind=FMS_MO_KIND_), intent(in)  :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_), intent(out) :: del_m, del_h, del_q

real(kind=FMS_MO_KIND_), dimension(1) :: z_1, z0_1, zt_1, zq_1, u_star_1, b_star_1, q_star_1, &
                      del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile (zref, zref_t, z_1, z0_1, zt_1, zq_1, &
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

real(kind=FMS_MO_KIND_),    intent(in),  dimension(:)   :: zref
real(kind=FMS_MO_KIND_),    intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=FMS_MO_KIND_),    intent(out), dimension(:,:) :: del_m, del_t, del_q
logical, intent(in) , optional,          dimension(:)   :: avail

integer :: k

do k = 1, size(zref(:))
  if(present(avail)) then
    call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
  else
      call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
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
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
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
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

return
end subroutine MO_PROFILE_2D_N_

!=======================================================================

subroutine MO_DIFF_2D_1_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:,:)      :: z, u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:)      :: k_m, k_h

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

real(kind=FMS_MO_KIND_), dimension(size(z),1,1)    :: z_n, k_m_n, k_h_n
real(kind=FMS_MO_KIND_), dimension(size(z),1)      :: u_star_n, b_star_n

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

real(kind=FMS_MO_KIND_), intent(in)       :: z, u_star, b_star
real(kind=FMS_MO_KIND_), intent(out)      :: k_m, k_h

integer                                   :: ni, nj, nk, ier
integer, parameter                        :: lkind = FMS_MO_KIND_
real(kind=FMS_MO_KIND_), parameter        :: ustar_min = 1.0E-10_lkind
real(kind=FMS_MO_KIND_), dimension(1,1,1) :: z_a, k_m_a, k_h_a
real(kind=FMS_MO_KIND_), dimension(1,1)   :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_1 in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = 1
z_a(1,1,1)    = z
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, FMS_MO_KIND_), ustar_min, neutral,               &
                        & stable_option, new_mo_option, real(rich_crit, FMS_MO_KIND_), &
                        & real(zeta_trans, FMS_MO_KIND_), ni, nj, nk, z_a, u_star_a,   & !miz
                        & b_star_a, k_m_a, k_h_a, ier)
k_m = k_m_a(1,1,1)
k_h = k_h_a(1,1,1)

end subroutine MO_DIFF_0D_1_

!=======================================================================

subroutine MO_DIFF_0D_N_(z, u_star, b_star, k_m, k_h)

real(kind=FMS_MO_KIND_), intent(in),  dimension(:) :: z
real(kind=FMS_MO_KIND_), intent(in)                :: u_star, b_star
real(kind=FMS_MO_KIND_), intent(out), dimension(:) :: k_m, k_h

integer                                            :: ni, nj, nk, ier
integer, parameter                                 :: lkind = FMS_MO_KIND_
real(kind=FMS_MO_KIND_), parameter                 :: ustar_min = 1.0E-10_lkind
real(kind=FMS_MO_KIND_), dimension(1,1,size(z))    :: z_a, k_m_a, k_h_a
real(kind=FMS_MO_KIND_), dimension(1,1)            :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = size(z(:))
z_a(1,1,:)    = z(:)
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, FMS_MO_KIND_), ustar_min, neutral,              &
                       & stable_option, new_mo_option, real(rich_crit, FMS_MO_KIND_), &
                       & real(zeta_trans, FMS_MO_KIND_), ni, nj, nk, z_a, u_star_a,   &
                       & b_star_a, k_m_a, k_h_a, ier)
k_m(:) = k_m_a(1,1,:)
k_h(:) = k_h_a(1,1,:)
end subroutine MO_DIFF_0D_N_

!=======================================================================

subroutine STABLE_MIX_2D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:,:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: i !< Loop index

n2 = size(mix, 2)

do i=1, n2
  call stable_mix(rich(:, i), mix(:, i))
enddo

end subroutine STABLE_MIX_2D_


!=======================================================================

subroutine STABLE_MIX_1D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in) , dimension(:)  :: rich
real(kind=FMS_MO_KIND_), intent(out), dimension(:)  :: mix
integer :: n !< Size of mix and rich
integer :: ierr !< Error code set by monin_obukhov_stable_mix

if (.not.module_is_initialized) call error_mesg('stable_mix in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(mix)

call monin_obukhov_stable_mix(stable_option, real(rich_crit,FMS_MO_KIND_), &
                              & real(zeta_trans,FMS_MO_KIND_), n, rich, mix, ierr)

end subroutine STABLE_MIX_1D_

!=======================================================================

subroutine STABLE_MIX_0D_(rich, mix)

real(kind=FMS_MO_KIND_), intent(in)       :: rich
real(kind=FMS_MO_KIND_), intent(out)      :: mix
real(kind=FMS_MO_KIND_), dimension(1)     :: mix_1d !< Representation of mix as a dimension(1) array

call stable_mix([rich], mix_1d)

mix = mix_1d(1)

end subroutine STABLE_MIX_0D_
!=======================================================================

!> @}
! close documentation grouping
