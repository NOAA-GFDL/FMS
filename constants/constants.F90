
module constants_mod

!-----------------------------------------------------------------------

!------------ physical constants ---------------

real, public, parameter :: radius = 6376.e3
real, public, parameter :: omega  = 7.292e-5
real, public, parameter :: grav   = 9.80
real, public, parameter :: rdgas  = 287.04
real, public, parameter :: kappa  = 2./7.
real, public, parameter :: cp     = rdgas/kappa
real, public, parameter :: p00    = 1000.e2

!------------ water vapor constants ---------------

real, public, parameter :: rvgas = 461.50
real, public, parameter :: dens_h2o = 1000.
real, public, parameter :: hlv = 2.500e6
real, public, parameter :: hlf = 3.34e5
real, public, parameter :: hls = 2.834e6
real, public, parameter :: tfreeze = 273.16

!------------ miscellaneous constants ---------------

real, public, parameter :: stefan  =  5.6734e-8
real, public, parameter :: vonkarm =  0.40

!! real, private, parameter :: mv = 18.016, md = 28.966
!! real, public, parameter :: d622  = rdgas/rvgas
!! real, public, parameter :: d378  = 1.0-d622, d608 = d378/d622

!-----------------------------------------------------------------------

end module constants_mod

