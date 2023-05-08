program test_diurnal_solar

  use fms_mod,           only: fms_init, fms_end, error_mesg
  use mpp_mod,           only : mpp_error, FATAL, stdout, mpp_init, mpp_exit
  use astronomy_mod,     only: astronomy_init, diurnal_solar
  use time_manager_mod,  only: JULIAN, set_calendar_type
  use platform_mod,      only: r4_kind, r8_kind
  use constants_mod,     only: PI

  implicit none

  ! Module variables
  real   :: ecc   = 0.01671_r8_kind!< Eccentricity of Earth's orbit [dimensionless]
  real   :: obliq = 23.439_r8_kind  !< Obliquity [degrees]
  real   :: per   = 102.932_r8_kind !< Longitude of perihelion with respect to autumnal equinox in NH [degrees]
  integer         :: num_angles = 3600
  real, parameter :: twopi = 2.0 * PI
  real, parameter :: deg_to_rad = twopi / 360.0 
  real, dimension(:), allocatable :: orb_angle



  call fms_init()
  call set_calendar_type(JULIAN)
  call astronomy_init
  
  call test_diurnal_solar_0d
  call fms_end

  contains

  subroutine test_diurnal_solar_0d
    real :: lat, lon, gmt, time_since_ae
    real :: cosz, fracday, rrsun
    real :: ref_angle, ref_dec, ref_rrsun
    real :: ang, dec

    !testing at the equator, at 12:00PM, during the atumnal equinox
    lat = 0.0
    lon = 0.0
    gmt = PI
    time_since_ae = 0.0

    call diurnal_solar(lat, lon, gmt, time_since_ae, cosz, fracday, rrsun)
    print *, "ang=", ang, "dec=", dec, "rrsun=", rrsun

    !expecting 1.0 for cosz because the sun is directly overehead during the equinox at the equator
    if (cosz .ne. 1.0 .or. fracday .ne. 1.0) call mpp_error(FATAL, &
    "test_diurnial_solar: incorrect scalar lat and lon values")

    !find reference value for rrsun at an angle of 0.0 for autumnal equinox
    ref_angle = angle(time_since_ae)
    ref_dec   = declination(ref_angle)
    ref_rrsun = r_inv_squared(ref_angle)
    print *, ref_angle, ref_dec, ref_rrsun
    !if (rrsun .ne.  ref_rrsun) call mpp_error(FATAL, &
    !"test_diurnal_solar: incorrect rrsun value")

  end subroutine test_diurnal_solar_0d

  !functions need to calculate reference values
  function angle(t)

    real, intent(in) :: t
    real :: angle
    real :: norm_time
    real :: x 
    integer :: int, int_1

    allocate (orb_angle(0:num_angles))

    norm_time = t * float(num_angles) / twopi
    int       = floor(norm_time)
    int       = modulo(int,num_angles)
    int_1     = int+1
    x         = norm_time - floor(norm_time)
    angle    = (1.0 - x) * orb_angle(int) + x * orb_angle(int_1)
    angle    = modulo(angle, twopi)

  end function angle

  function declination(ang)

    real, intent(in) :: ang
    real :: declination
    real :: rad_obliq
    real :: sin_dec

    rad_obliq = obliq * deg_to_rad
    sin_dec = - sin(rad_obliq) * sin(ang)
    declination = asin(sin_dec)

  end function declination

  function r_inv_squared(ang)

    real, intent(in) :: ang
    real :: r_inv_squared
    real :: r
    real :: rad_per

    rad_per  = per * deg_to_rad
    r        = 1.0 - ecc**2 / (1.0 + ecc * cos(ang - rad_per))
    r_inv_squared = r**(-2)

  end function r_inv_squared

  
end program test_diurnal_solar