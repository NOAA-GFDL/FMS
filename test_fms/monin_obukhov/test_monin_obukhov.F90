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

program test_monin_obukhov

  use monin_obukhov_inter, only: monin_obukhov_drag_1d, monin_obukhov_stable_mix, monin_obukhov_diff, monin_obukhov_profile_1d
  use mpp_mod, only : mpp_error, FATAL, stdout, mpp_init, mpp_exit

  implicit none
  integer, parameter :: i8 = selected_int_kind(18)
  integer(i8)        :: ier_tot, ier

  real    :: grav, vonkarm, error, zeta_min, small, ustar_min
  real    :: zref, zref_t
  integer :: max_iter

  real    :: rich_crit, zeta_trans
  real    :: drag_min_heat, drag_min_moist, drag_min_mom
  logical :: neutral
  integer :: stable_option
  logical :: new_mo_option

  call mpp_init()

  grav          = 9.80
  vonkarm       = 0.4
  error         = 1.0e-4
  zeta_min      = 1.0e-6
  max_iter      = 20
  small         = 1.0e-4
  neutral       = .false.
  stable_option = 1
  new_mo_option = .false.
  rich_crit     =10.0
  zeta_trans    = 0.5
  drag_min_heat  = 1.0e-5
  drag_min_moist= 1.0e-5
  drag_min_mom  = 1.0e-5
  ustar_min     = 1.e-10


  zref   = 10.
  zref_t = 2.


  ier_tot = 0
  call test_drag
  print *,'test_drag                    ier = ', ier
  ier_tot = ier_tot + ier

  call test_stable_mix
  print *,'test_stable_mix              ier = ', ier
  ier_tot = ier_tot + ier

  call test_diff
  print *,'test_diff                    ier = ', ier
  ier_tot = ier_tot + ier

  call test_profile
  print *,'test_profile                 ier = ', ier
  ier_tot = ier_tot + ier

  if(ier_tot/=0) then
      print *, 'Test_monin_obukhov: result different from expected result'
  else
     print *, 'No error detected.'
  endif

  call mpp_exit()

  CONTAINS

    subroutine test_drag

      integer(i8)        :: w

      integer :: i, ier_l, CHKSUM_DRAG
      integer, parameter :: n = 5
      logical :: avail(n), lavail

      real, dimension(n) :: pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star

      ! potential temperature
      pt     = (/ 268.559120403867, 269.799228886728, 277.443023238556, 295.79192777341, 293.268717243262 /)
      pt0    = (/ 273.42369841804 , 272.551410044203, 278.638168565727, 298.133068766049, 292.898163706587/)
      z      = (/ 29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475, 33.2184943356517/)
      z0     = (/ 5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05, 0.07/)
      zt     = (/ 3.69403636275411e-05, 0.0001, 1.01735489109205e-05, 7.63933834969505e-05, 0.00947346982656289/)
      zq     = (/ 5.72575636226887e-05, 0.0001, 5.72575636226887e-05, 5.72575636226887e-05, 5.72575636226887e-05/)
      speed  = (/ 2.9693638452068, 2.43308757772094, 5.69418282305367, 9.5608693754561, 4.35302260074334/)
      lavail = .true.
      avail  = (/.true., .true., .true., .true., .true. /)

      drag_m = 0
      drag_t = 0
      drag_q = 0
      u_star = 0
      b_star = 0

      call monin_obukhov_drag_1d(grav, vonkarm,               &
           & error, zeta_min, max_iter, small,                         &
           & neutral, stable_option, new_mo_option, rich_crit, zeta_trans,&
           & drag_min_heat, drag_min_moist, drag_min_mom,              &
           & n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t,         &
           & drag_q, u_star, b_star, lavail, avail, ier_l)

      ! check sum results
      w = 0
      w = w + transfer(sum(drag_m), w)
      w = w + transfer(sum(drag_t), w)
      w = w + transfer(sum(drag_q), w)
      w = w + transfer(sum(u_star), w)
      w = w + transfer(sum(b_star), w)

      ! plug in check sum here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95)
#define CHKSUM_DRAG 4466746452959549648
#endif
#if defined(_PGF95)
#define CHKSUM_DRAG 4466746452959549650
#endif


      print *,'chksum test_drag      : ', w, ' ref ', CHKSUM_DRAG
      ier = CHKSUM_DRAG - w

    end subroutine test_drag

    subroutine test_stable_mix

      integer(i8)        :: w

      integer, parameter      :: n = 5
      real   , dimension(n)   :: rich
      real   , dimension(n)   :: mix
      integer                 :: ier_l, CHKSUM_STABLE_MIX


      stable_option = 1
      rich_crit     = 10.0
      zeta_trans    =  0.5

      rich = (/1650.92431853365, 1650.9256285137, 77.7636819036559, 1.92806556391324, 0.414767442012442/)


      call monin_obukhov_stable_mix(stable_option, rich_crit, zeta_trans, &
           &                              n, rich, mix, ier_l)

      w = transfer( sum(mix) , w)

      ! plug in check sum here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95)
#define CHKSUM_STABLE_MIX 4590035772608644256
#endif
#if defined(_PGF95)
#define CHKSUM_STABLE_MIX 4590035772608644258
#endif

      print *,'chksum test_stable_mix: ', w, ' ref ', CHKSUM_STABLE_MIX
      ier = CHKSUM_STABLE_MIX - w

    end subroutine test_stable_mix

    !========================================================================

    subroutine test_diff

      integer(i8)        :: w

      integer, parameter             :: ni=1, nj=1, nk=1
      real   , dimension(ni, nj, nk) :: z
      real   , dimension(ni, nj)     :: u_star, b_star
      real   , dimension(ni, nj, nk) :: k_m, k_h
      integer                        :: ier_l, CHKSUM_DIFF

      z      = 19.9982554527751
      u_star = 0.129638955971075
      b_star = 0.000991799765557209

      call monin_obukhov_diff(vonkarm,                        &
           & ustar_min,                                     &
           & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &!miz
           & ni, nj, nk, z, u_star, b_star, k_m, k_h, ier_l)

      w = 0
      w = w + transfer( sum(k_m) , w)
      w = w + transfer( sum(k_h) , w)

      ! plug check sum in here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95) || defined(_PGF95)
#define CHKSUM_DIFF -9222066590093362639
#endif

      print *,'chksum test_diff      : ', w, ' ref ', CHKSUM_DIFF

      ier = CHKSUM_DIFF - w

    end subroutine test_diff

    !========================================================================

    subroutine test_profile

      integer(i8)        :: w

      integer, parameter :: n = 5
      integer            :: ier_l, CHKSUM_PROFILE

      logical :: avail(n)

      real, dimension(n) :: z, z0, zt, zq, u_star, b_star, q_star
      real, dimension(n) :: del_m, del_t, del_q

      z      = (/ 29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475, 33.2184943356517 /)
      z0     = (/ 5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05, 0.07/)
      zt     = (/ 3.69403636275411e-05, 0.0001, 1.01735489109205e-05, 7.63933834969505e-05, 0.00947346982656289/)
      zq     = (/ 5.72575636226887e-05, 0.0001, 5.72575636226887e-05, 5.72575636226887e-05, 5.72575636226887e-05/)
      u_star = (/ 0.109462510724615, 0.0932942802513508, 0.223232887323184, 0.290918439028557, 0.260087579361467/)
      b_star = (/ 0.00690834676781433, 0.00428178089592372, 0.00121229800895103, 0.00262353784027441, -0.000570314880866852/)
      q_star = (/ 0.000110861442197537, 9.44983279664197e-05, 4.17643828631936e-05, 0.000133135421415819, 9.36317815993945e-06/)

      avail = (/ .true., .true.,.true.,.true.,.true. /)

      call monin_obukhov_profile_1d(vonkarm, &
           & neutral, stable_option, new_mo_option, rich_crit, zeta_trans, &
           & n, zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
           & del_m, del_t, del_q, .true., avail, ier_l)

      ! check sum results
      w = 0
      w = w + transfer(sum(del_m), w)
      w = w + transfer(sum(del_t), w)
      w = w + transfer(sum(del_q), w)

      ! plug check sum in here>>>
#if defined(__INTEL_COMPILER) || defined(_LF95)
#define CHKSUM_PROFILE -4596910845317820786
#endif
#if defined(_PGF95)
#define CHKSUM_PROFILE -4596910845317820785
#endif

      print *,'chksum test_profile   : ', w, ' ref ', CHKSUM_PROFILE

      ier = CHKSUM_PROFILE - w

    end subroutine test_profile


end program test_monin_obukhov
