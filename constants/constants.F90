!! \brief Defines useful constants for Earth. Constants are defined as real
!! parameters. Constants are accessed through the "use" statement.
!!
!! \author Bruce Wyman <Bruce.Wyman@noaa.gov>
!!
!! \link http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/ \endlink
module constants_mod

implicit none
private

! Include variable "version" to be written to log file.
#include<file_version.h>
!-----------------------------------------------------------------------
! version is public so that write_version_number can be called for constants_mod
! by fms_init
public :: version

real :: realnumber !< dummy variable to use in HUGE initializations

real, public, parameter :: RADIUS = 6371.0e3 !< Radius of the Earth [m]
real, public, parameter :: OMEGA  = 7.292e-5 !< Rotation rate of the Earth [1/s]
real, public, parameter :: GRAV   = 9.80 !< Acceleration due to gravity [m/s^2]
real, public, parameter :: RDGAS  = 287.04 !< Gas constant for dry air [J/kg/deg]
real, public, parameter :: KAPPA  = 2./7. !< RDGAS / CP_AIR [dimensionless]
real, public, parameter :: CP_AIR = RDGAS/KAPPA !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
real, public, parameter :: CP_OCEAN = 3989.24495292815 !< Specific heat capacity taken from McDougall (2002) "Potential Enthalpy ..." [J/kg/deg]
real, public, parameter :: RHO0    = 1.035e3 !< Average density of sea water [kg/m^3]
real, public, parameter :: RHO0R   = 1.0/RHO0 !< Reciprocal of average density of sea water [m^3/kg]
real, public, parameter :: RHO_CP  = RHO0*CP_OCEAN !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]

real, public, parameter :: ES0 = 1.0 !< Humidity factor. Controls the humidity content of the atmosphere through
                                     !! the Saturation Vapour Pressure expression when using DO_SIMPLE. [dimensionless]
real, public, parameter :: RVGAS = 461.50 !< Gas constant for water vapor [J/kg/deg]
real, public, parameter :: CP_VAPOR = 4.0*RVGAS !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]
real, public, parameter :: DENS_H2O = 1000. !< Density of liquid water [kg/m^3]
real, public, parameter :: HLV = 2.500e6 !< Latent heat of evaporation [J/kg]
real, public, parameter :: HLF = 3.34e5 !< Latent heat of fusion [J/kg]
real, public, parameter :: HLS = HLV + HLF !< Latent heat of sublimation [J/kg]
real, public, parameter :: TFREEZE = 273.16 !< Freezing temperature of fresh water [K]

real, public, parameter :: WTMAIR = 2.896440E+01 !< Molecular weight of air [AMU]
real, public, parameter :: WTMH2O = WTMAIR*(RDGAS/RVGAS) !pjp OK to change value because not used yet. !< Molecular weight of water [AMU]
!real, public, parameter :: WTMO3  = 47.99820E+01
real, public, parameter :: WTMOZONE =  47.99820 !< Molecular weight of ozone [AMU]
real, public, parameter :: WTMC     =  12.00000 !< Molecular weight of carbon [AMU]
real, public, parameter :: WTMCO2   =  44.00995 !< Molecular weight of carbon dioxide [AMU]
real, public, parameter :: WTMO2    =  31.9988 !< Molecular weight of molecular oxygen [AMU]
real, public, parameter :: WTMCFC11 = 137.3681 !< Molecular weight of CFC-11 (CCl3F) [AMU]
real, public, parameter :: WTMCFC12 = 120.9135 !< Molecular weight of CFC-21 (CCl2F2) [AMU]
real, public, parameter :: DIFFAC = 1.660000E+00 !< Diffusivity factor [dimensionless]
real, public, parameter :: SECONDS_PER_DAY  = 8.640000E+04 !< Seconds in a day [s]
real, public, parameter :: SECONDS_PER_HOUR = 3600. !< Seconds in an hour [s]
real, public, parameter :: SECONDS_PER_MINUTE=60. !< Seconds in a minute [s]
real, public, parameter :: AVOGNO = 6.023000E+23 !< Avogadro's number [atoms/mole]
real, public, parameter :: PSTD   = 1.013250E+06 !< Mean sea level pressure [dynes/cm^2]
real, public, parameter :: PSTD_MKS    = 101325.0 !< Mean sea level pressure [N/m^2]
!real, public, parameter :: REARTH  = 6.356766E+08 !pjp Not used anywhere.

real, public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to heating rate in degrees per day [deg sec/(cm day)]
real, public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to heating rate in degrees per day [deg sec/(m day)]
real, public, parameter :: O2MIXRAT    = 2.0953E-01 !< Mixing ratio of molecular oxygen in air [dimensionless]
real, public, parameter :: RHOAIR      = 1.292269 !< Reference atmospheric density [kg/m^3]
real, public, parameter :: ALOGMIN     = -50.0 !< Minimum value allowed as argument to log function [N/A]

real, public, parameter :: STEFAN  = 5.6734e-8 !< Stefan-Boltzmann constant [W/m^2/deg^4]
real, public, parameter :: VONKARM = 0.40 !< Von Karman constant [dimensionless]
real, public, parameter :: PI_8    = 3.14159265358979323846d0 !< Ratio of circle circumference to diameter [N/A] (REAL(KIND=8))
real, public, parameter :: PI      = 3.14159265358979323846 !< Ratio of circle circumference to diameter [N/A]
real, public, parameter :: RAD_TO_DEG=180./PI !< Degrees per radian [deg/rad]
real, public, parameter :: DEG_TO_RAD=PI/180. !< Radians per degree [rad/deg]
real, public, parameter :: RADIAN  = RAD_TO_DEG !< Equal to RAD_TO_DEG. Named RADIAN for backward compatability. [rad/deg]
real, public, parameter :: C2DBARS = 1.e-4 !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
real, public, parameter :: KELVIN  = 273.15 !< Degrees Kelvin at zero Celsius [K]
real, public, parameter :: EPSLN   = 1.0e-40 !< A small number to prevent divide by zero exceptions [N/A]
!-----------------------------------------------------------------------
public :: constants_init

contains

!> \brief dummy routine.
subroutine constants_init

end subroutine constants_init

end module constants_mod

!   <FUTURE>
!   1.  Renaming of constants.
!   </FUTURE>
!   <FUTURE>
!   2.  Additional constants.
!   </FUTURE>
!   <NOTE>
!    Constants have been declared as type REAL, PARAMETER.
!
!    The value a constant can not be changed in a users program.
!    New constants can be defined in terms of values from the
!    constants module using a parameter statement.<br><br>
!
!    The name given to a particular constant may be changed.<br><br>
!
!    Constants can be used on the right side on an assignment statement
!    (their value can not be reassigned).
!
!
!<TESTPROGRAM NAME="EXAMPLE">
!<PRE>
!    use constants_mod, only:  TFREEZE, grav_new => GRAV
!    real, parameter :: grav_inv = 1.0 / grav_new
!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!    geopotential(:,:) = height(:,:) * grav_new
!</PRE>
!</TESTPROGRAM>
!   </NOTE>
