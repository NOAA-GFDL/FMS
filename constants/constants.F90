
module constants_mod

!-----------------------------------------------------
!
!  Defines useful constants for Earth in mks units.
!
!-----------------------------------------------------

use fms_mod, only: write_version_number

implicit none
private

character(len=128) :: version='$Id: constants.F90,v 1.3 2002/07/16 22:54:42 fms Exp $'
character(len=128) :: tagname='$Name: havana $'
logical :: do_log = .true.
logical :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------
!------------ physical constants ---------------

real, public, parameter :: RADIUS = 6376.e3      !  radius of the earth (meters)
real, public, parameter :: OMEGA  = 7.292e-5     !  rotation rate of planet (1/sec)
real, public, parameter :: GRAV   = 9.80         !  acceleration due to gravity (m/s2)
real, public, parameter :: RDGAS  = 287.04       !  gas constant for dry air (J/Kg/deg)
real, public, parameter :: KAPPA  = 2./7.        !  RDGAS / CP
real, public, parameter :: cp     = RDGAS/KAPPA  !  spec heat cap of dry air (J/kg/deg)

!------------ water vapor constants ---------------

real, public, parameter :: RVGAS = 461.50        !  gas constant for water vapor (J/Kg/deg)
real, public, parameter :: DENS_H2O = 1000.      !  density of liquid water (Kg/m3)
real, public, parameter :: HLV = 2.500e6         !  latent heat of evaporation (J/Kg)
real, public, parameter :: HLF = 3.34e5          !  latent heat of fusion (J/Kg)
real, public, parameter :: HLS = 2.834e6         !  latent heat of sublimation (J/Kg)
real, public, parameter :: TFREEZE = 273.16      !  temp where fresh water freezes (deg K)

!------------ miscellaneous constants ---------------

real, public, parameter :: STEFAN  =  5.6734e-8  !  Stefan-Boltzmann constant (W/m2/deg4)
real, public, parameter :: VONKARM =  0.40       !  Von Karman constant
real, public, parameter :: PI      =  3.14159265358979323846 ! is it enough?
!-----------------------------------------------------------------------

public constants_init

contains

! optional initialization routine
! only purpose is to write version to log file

subroutine constants_init

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

  if (.not.do_log) return
  call write_version_number (version,tagname)
  do_log = .false.

end subroutine constants_init

!-----------------------------------------------------------------------

end module constants_mod

