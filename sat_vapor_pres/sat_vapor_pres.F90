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
!> @defgroup sat_vapor_pres_mod sat_vapor_pres_mod
!> @ingroup sat_vapor_pres
!> @brief Routines for computing the saturation vapor pressure (es),
!! the specific humidity (qs) and vapor mixing ratio (mrs)
!> Given a specified relative humidity, calculates es, qs, and mrs, as well as their
!! derivatives with respect to temperature, and also includes routines
!! to initialize the look-up table.
!! This module contains routines for determining the saturation vapor
!! pressure (<TT>ES</TT>) from lookup tables constructed using equations given
!! in the Smithsonian tables.  The <TT>ES</TT> lookup tables are valid between
!! -160C and +100C (approx 113K to 373K).
!!
!! The values of <TT>ES</TT> are computed over ice from -160C to -20C,
!! over water from 0C to 100C, and a blended value (over water and ice)
!! from -20C to 0C.
!!
!! Routines are also included to calculate the saturation specific
!! humidity and saturation mixing ratio for vapor, and their deriv-
!! atives with respect to temperature.  By default, the values returned
!! are those at saturation; optionally, values of q and mr at a spec-
!! ified relative humidity may instead be returned. Two forms are
!! available; the approximate form that has been traditionally used in
!! GCMs, and an exact form provided by SJ Lin in which saturation is
!! reached while maintaining constant pressure and temperature.
!!
!! This version was written for non-vector machines.
!! See the <LINK SRC="#NOTES">notes</LINK> section for details on vectorization.
!!
!!    arguments
!!    ---------
!!      temp    intent in       temperature in degrees kelvin
!!      es      intent out      saturation vapor pressure in Pascals
!!      des     intent out      derivative of saturation vapor pressure
!!                              with respect to temperature
!!                              (Pascals/degree)
!!      press   intent in       atmospheric pressure in Pascals
!!      qs      intent out      specific humidity at relative humidity hc
!!                              (kg(vapor) / kg(moist air)
!!      mrs     intent out      mixing ratio at relative humidity hc
!!                              (kg(vapor) / kg(dry air)
!!
!!   optional arguments
!!   ------------------
!!      q       intent in       vapor specific humidity
!!                              (kg(vapor) / kg(moist air)
!!      hc      intent in       relative humidity at which output
!!                              fields are desired: default is 100 %
!!      dqsdT   intent out      derivative of saturation specific
!!                              humidity with respect to temperature
!!                              (kg(vapor) / kg(moist air) /degree)
!!      mr      intent in       vapor mixing ratio
!!                              (kg(vapor) / kg(dry air)
!!      dmrsdT  intent out      derivative of saturation mixing ratio
!!                              with respect to temperature
!!                              (kg(vapor) / kg(dry air) /degree)
!!      esat    intent out      saturation vapor pressure
!!                              (Pascals)
!!      err_msg intent out      character string to hold error message
!!      es_over_liq
!!              intent  in      use es table wrt liquid only
!!
!! Example Usages:
!!
!!              call lookup_es  (temp, es, err_msg)
!!
!!              call lookup_des (temp, des, err_msg)
!!
!!              call lookup_es_des (temp, es, des, err_msg)
!!
!!              call lookup_es2 (temp, es, err_msg)
!!
!!              call lookup_des2 (temp, des, err_msg)
!!
!!              call lookup_es2_des2 (temp, es, des, err_msg)
!!
!!              call compute_qs (temp, press, qs, q, hc, dqsdT, esat,
!!                               err_msg, es_over_liq)
!!
!!              call compute_mrs (temp, press, mrs, mr, hc, dmrsdT, esat,
!!                                err_msg, es_over_liq)

module sat_vapor_pres_mod

!-----------------------------------------------------------------------
!
!
!    arguments
!    ---------
!      temp    intent in       temperature in degrees kelvin
!      es      intent out      saturation vapor pressure in Pascals
!      des     intent out      derivative of saturation vapor pressure
!                              with respect to temperature
!                              (Pascals/degree)
!      press   intent in       atmospheric pressure in Pascals
!      qs      intent out      specific humidity at relative humidity hc
!                              (kg(vapor) / kg(moist air)
!      mrs     intent out      mixing ratio at relative humidity hc
!                              (kg(vapor) / kg(dry air)
!
!   optional arguments
!   ------------------
!      q       intent in       vapor specific humidity
!                              (kg(vapor) / kg(moist air)
!      hc      intent in       relative humidity at which output
!                              fields are desired: default is 100 %
!      dqsdT   intent out      derivative of saturation specific
!                              humidity with respect to temperature
!                              (kg(vapor) / kg(moist air) /degree)
!      mr      intent in       vapor mixing ratio
!                              (kg(vapor) / kg(dry air)
!      dmrsdT  intent out      derivative of saturation mixing ratio
!                              with respect to temperature
!                              (kg(vapor) / kg(dry air) /degree)
!      esat    intent out      saturation vapor pressure
!                              (Pascals)
!      err_msg intent out      character string to hold error message
!      es_over_liq
!              intent  in      use es table wrt liquid only
!
!-----------------------------------------------------------------------

! <CONTACT EMAIL="Bruce.Wyman@noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Routines for determining the saturation vapor pressure
!   (<TT>ES</TT>), saturation vapor specific humidity and saturation
!   vapor mixing ratio, and their derivatives with respect to
!   temperature.
! </OVERVIEW>

! <DESCRIPTION>
!   This module contains routines for determining the saturation vapor
!   pressure (<TT>ES</TT>) from lookup tables constructed using equations given
!   in the Smithsonian tables.  The <TT>ES</TT> lookup tables are valid between
!   -160C and +100C (approx 113K to 373K).

!   The values of <TT>ES</TT> are computed over ice from -160C to -20C,
!   over water from 0C to 100C, and a blended value (over water and ice)
!   from -20C to 0C.

!   Routines are also included to calculate the saturation specific
!   humidity and saturation mixing ratio for vapor, and their deriv-
!   atives with respect to temperature.  By default, the values returned
!   are those at saturation; optionally, values of q and mr at a spec-
!   ified relative humidity may instead be returned. Two forms are
!   available; the approximate form that has been traditionally used in
!   GCMs, and an exact form provided by SJ Lin in which saturation is
!   reached while maintaining constant pressure and temperature.

!   This version was written for non-vector machines.
!   See the <LINK SRC="#NOTES">notes</LINK> section for details on vectorization.

! </DESCRIPTION>

! <PUBLIC>
!   Description summarizing public interface.
! </PUBLIC>

 use         constants_mod, only:  TFREEZE, RDGAS, RVGAS, HLV, ES0
 use        fms_mod, only:  write_version_number, stdout, stdlog, mpp_pe, mpp_root_pe, &
                            mpp_error, FATAL, fms_error_handler,   &
                            error_mesg, check_nml_error
 use        mpp_mod, only: input_nml_file
 use  sat_vapor_pres_k_mod, only:  sat_vapor_pres_init_k, lookup_es_k, &
                                   lookup_des_k, lookup_es_des_k, &
                                   lookup_es2_k,  &
                                   lookup_des2_k, lookup_es2_des2_k, &
                                   lookup_es3_k,  &
                                   lookup_des3_k, lookup_es3_des3_k, &
                                   compute_qs_k, compute_mrs_k
 use platform_mod, only: r4_kind, r8_kind

implicit none
private

 public :: lookup_es, lookup_des, sat_vapor_pres_init
 public :: lookup_es2, lookup_des2, lookup_es2_des2
 public :: lookup_es3, lookup_des3, lookup_es3_des3
 public :: lookup_es_des, compute_qs, compute_mrs
!public :: compute_es
 public :: escomp, descomp ! for backward compatibility
                           ! use lookup_es, lookup_des instead

!-----------------------------------------------------------------------

! <INTERFACE NAME="lookup_es">

!   <OVERVIEW>
!     For the given temperatures, returns the saturation vapor pressures.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the given temperatures these routines return the
!     saturation vapor pressure (esat). The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call lookup_es( temp, esat, err_msg )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure in pascals.
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="err_msg" UNITS="      " TYPE="character">
!     Character string containing error message to be returned to
!     calling routine.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 !> @brief For the given temperatures, returns the saturation vapor pressures
 !!
 !> For the given temperatures these routines return the saturation vapor pressure(esat).
 !! The return values are derived from lookup tables.
 !! Example usage:
 !! @code{.F90} call lookup_es( temp, esat, err_msg ) @endcode
 !!
 !! @param temp Temperature in degrees Kelvin.
 !! @param esat Saturation vapor pressure in pascals.
 !!             May be a scalar, 1d, 2d, or 3d array
 !!             Must have the same order and size as temp.
 !! @param err_msg Character string containing error message to be returned to
 !!                calling routine.
 !! @throws FATAL table overflow, nbad=##
 !!     Temperature(s) provided to the saturation vapor pressure lookup
 !!          are outside the valid range of the lookup table (-160 to 100 deg C).
 !!          This may be due to a numerical instability in the model.
 !!          Information should have been printed to standard output to help
 !!          determine where the instability may have occurred.
 !!          If the lookup table needs a larger temperature range,
 !!          then parameters in the module header must be modified.
 !> @ingroup sat_vapor_pres_mod
 interface lookup_es
   module procedure lookup_es_0d_r4, lookup_es_0d_r8
   module procedure lookup_es_1d_r4, lookup_es_1d_r8
   module procedure lookup_es_2d_r4, lookup_es_2d_r8
   module procedure lookup_es_3d_r4, lookup_es_3d_r8
 end interface lookup_es
 !> Provided for backward compatibility (to be removed soon)
 !> @ingroup sat_vapor_pres_mod
 interface escomp
   module procedure lookup_es_0d_r4, lookup_es_0d_r8
   module procedure lookup_es_1d_r4, lookup_es_1d_r8
   module procedure lookup_es_2d_r4, lookup_es_2d_r8
   module procedure lookup_es_3d_r4, lookup_es_3d_r8
 end interface escomp
! </INTERFACE>
!-----------------------------------------------------------------------
! <INTERFACE NAME="lookup_des">

!   <OVERVIEW>
!     For the given temperatures, returns the derivative of saturation vapor pressure
!     with respect to temperature.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the given temperatures these routines return the derivative of esat w.r.t.
!     temperature (desat). The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call lookup_des( temp, desat )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Derivative of saturation vapor pressure w.r.t. temperature
!                 in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="err_msg" UNITS="      " TYPE="character">
!     Character string containing error message to be returned to
!     calling routine.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 !> For the given temperatures, returns the derivative of saturation vapor pressure
 !! with respect to temperature.
 !!
 !! For the given temperatures these routines return the derivtive of esat w.r.t. temperature
 !! (desat). The return values are derived from lookup tables.
 !!
 !! @param [in] temp Temperature in degrees kelvin
 !! @param [out] desat Derivative of saturation vapor pressure w.r.t. temperature
 !!                 in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
 !!                 Must have the same order and size as temp.
 !! @param [out] err_msg Character string containing error message to be returned to
 !!     calling routine.
 !!
 !! @error FATAL table overflow, nbad=##
 !!        Temperature(s) provided to the saturation vapor pressure lookup
 !!        are outside the valid range of the lookup table (-160 to 100 deg C).
 !!        This may be due to a numerical instability in the model.
 !!        Information should have been printed to standard output to help
 !!        determine where the instability may have occurred.
 !!        If the lookup table needs a larger temperature range,
 !!        then parameters in the module header must be modified.
 !!
 !! <br>Example usage:
 !! @code{.F90} call lookup_des( temp, desat) @endcode
 !> @ingroup sat_vapor_pres_mod
 interface lookup_des
   module procedure lookup_des_0d_r4, lookup_des_0d_r8
   module procedure lookup_des_1d_r4, lookup_des_1d_r8
   module procedure lookup_des_2d_r4, lookup_des_2d_r8
   module procedure lookup_des_3d_r4, lookup_des_3d_r8
 end interface lookup_des
! </INTERFACE>
 !> Provided for backward compatibility (to be removed soon)
 !> @ingroup sat_vapor_pres_mod
 interface descomp
   module procedure lookup_des_0d_r4, lookup_des_0d_r8
   module procedure lookup_des_1d_r4, lookup_des_1d_r8
   module procedure lookup_des_2d_r4, lookup_des_2d_r8
   module procedure lookup_des_3d_r4, lookup_des_3d_r8
 end interface descomp

!-----------------------------------------------------------------------

! <INTERFACE NAME="lookup_es_des">

!   <OVERVIEW>
!     For the given temperatures, returns the saturation vapor pressure
!     and the derivative of saturation vapor pressure with respect to
!     temperature.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the given temperatures these routines return the
!     saturation vapor pressure (esat) and the derivative of esat w.r.t
!     temperature (desat). The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call lookup_es_des( temp, esat, desat, err_msg )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure in pascals.
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="desat" UNITS="pascal/ degree" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Derivative of saturation vapor pressure w.r.t. temperature
!                 in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="err_msg" UNITS="      " TYPE="character">
!     Character string containing error message to be returned to
!     calling routine.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 !> @brief For the given temperatures, returns the saturation vapor pressure
 !! and the derivative of saturation vapor pressure with respect to
 !! temperature.
 !!
 !> For the given temperatures these routines return the
 !! saturation vapor pressure (esat) and the derivative of esat w.r.t
 !! temperature (desat). The return values are derived from
 !! lookup tables (see notes below).
 !!
 !! <br>Example usage:
 !! @code{.F90} call lookup_es_des( temp, esat, desat, err_msg ) @endcode
 !!
 !! @param temp Temperature in degrees Kelvin.
 !! @param [out] esat Saturation vapor pressure in pascals. May be a scalar, 1d, 2d, or 3d array.
 !!              Must have the same order and size as temp.
 !! @param [out] desat Derivative of saturation vapor pressure w.r.t. temperature
 !!                    in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
 !!                    Must have the same order and size as temp.
 !! @param [out] err_msg Character string containing error message to be returned to
 !!                      calling routine.
 !! @error FATAL table overflow, nbad=##
 !! Temperature(s) provided to the saturation vapor pressure lookup
 !! are outside the valid range of the lookup table (-160 to 100 deg C).
 !! This may be due to a numerical instability in the model.
 !! Information should have been printed to standard output to help
 !! determine where the instability may have occurred.
 !! If the lookup table needs a larger temperature range,
 !! then parameters in the module header must be modified.
 !> @ingroup sat_vapor_pres_mod
 interface lookup_es_des
   module procedure lookup_es_des_0d_r4, lookup_es_des_0d_r8
   module procedure lookup_es_des_1d_r4, lookup_es_des_1d_r8
   module procedure lookup_es_des_2d_r4, lookup_es_des_2d_r8
   module procedure lookup_es_des_3d_r4, lookup_es_des_3d_r8
 end interface lookup_es_des

 !> @ingroup sat_vapor_pres_mod
 interface lookup_es2
   module procedure lookup_es2_0d_r4, lookup_es2_0d_r8
   module procedure lookup_es2_1d_r4, lookup_es2_1d_r8
   module procedure lookup_es2_2d_r4, lookup_es2_2d_r8
   module procedure lookup_es2_3d_r4, lookup_es2_3d_r8
 end interface lookup_es2

 !> @ingroup sat_vapor_pres_mod
 interface lookup_des2
   module procedure lookup_des2_0d_r4, lookup_des2_0d_r8
   module procedure lookup_des2_1d_r4, lookup_des2_1d_r8
   module procedure lookup_des2_2d_r4, lookup_des2_2d_r8
   module procedure lookup_des2_3d_r4, lookup_des2_3d_r8
 end interface lookup_des2

 !> @ingroup sat_vapor_pres_mod
 interface lookup_es2_des2
   module procedure lookup_es2_des2_0d_r4, lookup_es2_des2_0d_r8
   module procedure lookup_es2_des2_1d_r4, lookup_es2_des2_1d_r8
   module procedure lookup_es2_des2_2d_r4, lookup_es2_des2_2d_r8
   module procedure lookup_es2_des2_3d_r4, lookup_es2_des2_3d_r8
 end interface lookup_es2_des2

 !> @ingroup sat_vapor_pres_mod
 interface lookup_es3
   module procedure lookup_es3_0d_r4, lookup_es3_0d_r8
   module procedure lookup_es3_1d_r4, lookup_es3_1d_r8
   module procedure lookup_es3_2d_r4, lookup_es3_2d_r8
   module procedure lookup_es3_3d_r4, lookup_es3_3d_r8
 end interface lookup_es3

 !> @ingroup sat_vapor_pres_mod
 interface lookup_des3
   module procedure lookup_des3_0d_r4, lookup_des3_0d_r8
   module procedure lookup_des3_1d_r4, lookup_des3_1d_r8
   module procedure lookup_des3_2d_r4, lookup_des3_2d_r8
   module procedure lookup_des3_3d_r4, lookup_des3_3d_r8
 end interface lookup_des3

 !> @ingroup sat_vapor_pres_mod
 interface lookup_es3_des3
   module procedure lookup_es3_des3_0d_r4, lookup_es3_des3_0d_r8
   module procedure lookup_es3_des3_1d_r4, lookup_es3_des3_1d_r8
   module procedure lookup_es3_des3_2d_r4, lookup_es3_des3_2d_r8
   module procedure lookup_es3_des3_3d_r4, lookup_es3_des3_3d_r8
 end interface lookup_es3_des3

!-----------------------------------------------------------------------

! <INTERFACE NAME="compute_qs">

!   <OVERVIEW>
!     For the given temperatures, pressures and optionally vapor
!     specific humidity, returns the specific humidity at saturation
!     (optionally at relative humidity hc instead of at saturation) and
!     optionally the derivative of saturation specific humidity w.r.t.
!     temperature, and the saturation vapor pressure.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the input temperature and pressure these routines return the
!     specific humidity (qsat) at saturation (unless optional argument
!     hc is used to specify the relative humidity at which qsat should
!     apply) and, if desired, the derivative of qsat w.r.t temperature
!     (dqsdT) and / or the saturation vapor pressure (esat). If the
!     optional input argument specific humidity (q) is present, the
!     exact expression for qs is used; if q is not present the tradit-
!     ional form (valid at saturation) is used. if the optional qsat
!     derivative argument is present, the derivative of qsat w.r.t.
!     temperature will also be returned, defined consistent with the
!     expression used for qsat. The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call compute_qs( temp, press, qsat, q, hc, dqsdT, esat, err_msg )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Air pressure in Pascals.
!   </IN>
!   <OUT NAME="qsat" UNITS="kg(vapor) / kg(moist air)" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Specific humidity in kg (vapor) / kg (moist air)
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>
!   <IN NAME="q" UNIT="kg(vapor) / kg (moist air)" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Vapor specific humidity in kg (vapor) / kg (moist air).
!     If present, exact formulation for qsat and dqsdT will be used.
!   </IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)">
!     Relative humidity at which output variables are desired.
!     If not present, values will apply at saturation.
!   </IN>
!   <OUT NAME="dqsdT" UNITS="kg(vapor) / kg(moist air) / degree" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Derivative of saturation specific humidity w.r.t. temperature
!                 in kg(vapor) / kg(moist air) / degree. May be a
!                 scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure. May be a scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="err_msg" UNITS="      " TYPE="character">
!     Character string containing error message to be returned to
!     calling routine.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 !> @brief For the given temperatures, pressures and optionally vapor
 !! specific humidity, returns the specific humidity at saturation
 !! (optionally at relative humidity hc instead of at saturation) and
 !! optionally the derivative of saturation specific humidity w.r.t.
 !! temperature, and the saturation vapor pressure.
 !!
 !! For the input temperature and pressure these routines return the
 !! specific humidity (qsat) at saturation (unless optional argument
 !! hc is used to specify the relative humidity at which qsat should
 !! apply) and, if desired, the derivative of qsat w.r.t temperature
 !! (dqsdT) and / or the saturation vapor pressure (esat). If the
 !! optional input argument specific humidity (q) is present, the
 !! exact expression for qs is used; if q is not present the tradit-
 !! ional form (valid at saturation) is used. if the optional qsat
 !! derivative argument is present, the derivative of qsat w.r.t.
 !! temperature will also be returned, defined consistent with the
 !! expression used for qsat. The return values are derived from
 !! lookup tables (see notes below).
 !!
 !! Example usage:
 !! @code{.F90} call compute_qs( temp, press, qsat, q, hc, dqsdT, esat, err_msg ) @endcode
 !!
 !> @ingroup sat_vapor_pres_mod
 interface compute_qs
   module procedure compute_qs_0d_r4, compute_qs_0d_r8
   module procedure compute_qs_1d_r4, compute_qs_1d_r8
   module procedure compute_qs_2d_r4, compute_qs_2d_r8
   module procedure compute_qs_3d_r4, compute_qs_3d_r8
 end interface compute_qs

!-----------------------------------------------------------------------

! <INTERFACE NAME="compute_mrs">

!   <OVERVIEW>
!     For the given temperatures, pressures and optionally vapor
!     mixing ratio, returns the  vapor mixing ratio at saturation
!     (optionally at relative humidity hc instead of at saturation) and
!     optionally the derivative of saturation vapor mixing ratio w.r.t.
!     temperature, and the saturation vapor pressure.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the input temperature and pressure these routines return the
!     vapor mixing ratio (mrsat) at saturation (unless optional argument
!     hc is used to specify the relative humidity at which mrsat should
!     apply) and, if desired, the derivative of mrsat w.r.t temperature
!     (dmrsdT) and / or the saturation vapor pressure (esat). If the
!     optional input argument specific humidity (mr) is present, the
!     exact expression for mrs is used; if qr is not present the tradit-
!     ional form (valid at saturation) is used. if the optional mrsat
!     derivative argument is present, the derivative of mrsat w.r.t.
!     temperature will also be returned, defined consistent with the
!     expression used for mrsat. The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call compute_mrs( temp, press, mrsat, mr, hc, dmrsdT, esat,
!                       err_msg )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <IN NAME="press" UNIT="Pascals" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Air pressure in Pascals.
!   </IN>
!   <OUT NAME="mrsat" UNITS="kg(vapor) / kg (dry air)" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Vapor mixing ratio in kg (vapor) / kg (dry air)
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>
!   <IN NAME="mr" UNIT="kg(vapor) / kg (dry air)" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Vapor mixing ratio in kg (vapor) / kg (dry air).
!     If present, exact formulation for mrsat and dmrsdT will be used.
!   </IN>
!   <IN NAME="hc" UNIT="fraction" TYPE="real" DIM="(scalar)">
!     Relative humidity at which output variables are desired.
!     If not present, values will apply at saturation.
!   </IN>
!   <OUT NAME="dmrsdT" UNITS="kg(vapor) / kg(dry air) / degree" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Derivative of saturation vapor mixing ratio w.r.t. temperature
!                 in kg(vapor) / kg(dry air) / degree. May be a
!                 scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="esat" UNITS="Pascals" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure. May be a scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>
!   <OUT NAME="err_msg" UNITS="      " TYPE="character">
!     Character string containing error message to be returned to
!     calling routine.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 !> For the given temperatures, pressures and optionally vapor
 !! mixing ratio, returns the  vapor mixing ratio at saturation
 !! (optionally at relative humidity hc instead of at saturation) and
 !! optionally the derivative of saturation vapor mixing ratio w.r.t.
 !! temperature, and the saturation vapor pressure.
 !!
 !! For the input temperature and pressure these routines return the
 !! vapor mixing ratio (mrsat) at saturation (unless optional argument
 !! hc is used to specify the relative humidity at which mrsat should
 !! apply) and, if desired, the derivative of mrsat w.r.t temperature
 !! (dmrsdT) and / or the saturation vapor pressure (esat). If the
 !! optional input argument specific humidity (mr) is present, the
 !! exact expression for mrs is used; if qr is not present the tradit-
 !! ional form (valid at saturation) is used. if the optional mrsat
 !! derivative argument is present, the derivative of mrsat w.r.t.
 !! temperature will also be returned, defined consistent with the
 !! expression used for mrsat. The return values are derived from
 !! lookup tables (see notes below).
 !!
 !! <br>Example usage:
 !! @code{.F90} call compute_mrs( temp, press, mrsat, mr, hc, dmrsdT, esat,
 !!                       err_msg ) @endcode
 !> @ingroup sat_vapor_pres_mod
 interface compute_mrs
   module procedure compute_mrs_0d_r4, compute_mrs_0d_r8
   module procedure compute_mrs_1d_r4, compute_mrs_1d_r8
   module procedure compute_mrs_2d_r4, compute_mrs_2d_r8
   module procedure compute_mrs_3d_r4, compute_mrs_3d_r8
 end interface compute_mrs

!-----------------------------------------------------------------------
! <INTERFACE NAME="compute_es">

!   <OVERVIEW>
!     For the given temperatures, computes the saturation vapor pressures.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Computes saturation vapor pressure for the given temperature using
!     the equations given in the Smithsonian Meteorological Tables.
!     Between -20C and 0C a blended value over ice and water is returned.
!   </DESCRIPTION>
!   <TEMPLATE>
!     es = compute_es ( temp )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure in pascals.
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>

!interface compute_es
!  module procedure compute_es_0d, compute_es_1d, compute_es_2d, compute_es_3d
!end interface
! </INTERFACE>
!-----------------------------------------------------------------------
 !> @ingroup sat_vapor_pres_mod
 interface check_1d
    module procedure check_1d_r4, check_1d_r8
 end interface check_1d

 interface check_2d
    module procedure check_2d_r4, check_2d_r8
 end interface check_2d

 !> @ingroup sat_vapor_pres_mod
 interface temp_check
   module procedure temp_check_1d_r4, temp_check_1d_r8
   module procedure temp_check_2d_r4, temp_check_2d_r8
   module procedure temp_check_3d_r4, temp_check_3d_r8
 end interface temp_check

 !> @ingroup sat_vapor_pres_mod
 interface show_all_bad
   module procedure show_all_bad_0d_r4, show_all_bad_0d_r8
   module procedure show_all_bad_1d_r4, show_all_bad_1d_r8
   module procedure show_all_bad_2d_r4, show_all_bad_2d_r8
   module procedure show_all_bad_3d_r4, show_all_bad_3d_r8
 end interface show_all_bad

!> @addtogroup sat_vapor_pres_mod
!> @{
!-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>

 logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------
!  parameters for use in computing qs and mrs

 real(r8_kind), parameter :: EPSILO = real(RDGAS,r8_kind)/real(RVGAS, r8_kind)
 real(r8_kind), parameter :: ZVIR = real(RVGAS,r8_kind)/real(RDGAS,r8_kind) - 1.0_r8_kind

!-----------------------------------------------------------------------
!  parameters for table size and resolution

 integer, public :: tcmin = -160  ! minimum temperature (degC) in lookup table
 integer, public :: tcmax =  100  ! maximum temperature (degC) in lookup table
 integer :: esres =  10   ! table resolution (increments per degree)
 integer :: nsize  ! (tcmax-tcmin)*esres+1    !  lookup table size
 integer :: nlim   ! nsize-1

 integer :: stdoutunit=0
!-----------------------------------------------------------------------
!  variables needed by temp_check
 real(r8_kind) :: tmin, dtinv, teps

! The default values below preserve the behavior of omsk and earlier revisions.
 logical :: show_bad_value_count_by_slice=.true.
 logical :: show_all_bad_values=.false.
 logical :: use_exact_qs = .false.
 logical :: do_simple             =.false.
 logical :: construct_table_wrt_liq = .false.
 logical :: construct_table_wrt_liq_and_ice = .false.

 namelist / sat_vapor_pres_nml / show_bad_value_count_by_slice, show_all_bad_values, &
                                 use_exact_qs, do_simple, &
                                 construct_table_wrt_liq, &
                                 construct_table_wrt_liq_and_ice

contains

 subroutine sat_vapor_pres_init(err_msg)

!  =================================================================
!  +                                                               +
!  +             construction of the es table                      +
!  +                                                               +
!  + this table is constructed from es equations from the          +
!  + smithsonian tables.  the es input is computed from values     +
!  + (in one-tenth of a degree increments) of es over ice          +
!  + from -153c to 0c and values of es over water from 0c to 102c. +
!  + output table contains these data interleaved with their       +
!  + derivatives with respect to temperature except between -20c   +
!  + and 0c where blended (over water and over ice) es values and  +
!  + derivatives are calculated.                                   +
!  +   note: all es computation is done in pascals                 +
!  =================================================================

  character(len=*), intent(out), optional :: err_msg
  character(len=128) :: err_msg_local
  integer :: unit, ierr, io

! return silently if this routine has already been called
  if (module_is_initialized) return

!---- read namelist input ----
  read (input_nml_file, sat_vapor_pres_nml, iostat=io)
  ierr = check_nml_error(io,'sat_vapor_pres_nml')

! write version number and namelist to log file
  call write_version_number("SAT_VAPOR_PRES_MOD", version)
  unit = stdlog()
  stdoutunit = stdout()
  if (mpp_pe() == mpp_root_pe()) write (unit, nml=sat_vapor_pres_nml)

  if(do_simple) then
    tcmin = -173
    tcmax =  350
  endif
  nsize = (tcmax-tcmin)*esres+1
  nlim  = nsize-1
  call sat_vapor_pres_init_k(nsize, real(tcmin,r8_kind), real(tcmax,r8_kind), &
                             real(TFREEZE,r8_kind), real(HLV,r8_kind),&
                             real(RVGAS,r8_kind), real(ES0,r8_kind), err_msg_local, use_exact_qs, do_simple,&
                             construct_table_wrt_liq, &
                             construct_table_wrt_liq_and_ice, &
                             teps, tmin, dtinv)
  if ( err_msg_local == '' ) then
     if(present(err_msg)) err_msg = ''
  else
     if(fms_error_handler('lookup_es',err_msg_local,err_msg)) return
  endif

  module_is_initialized = .true.

end subroutine sat_vapor_pres_init

#include "sat_vapor_pres_r4.fh"
#include "sat_vapor_pres_r8.fh"

!#######################################################################
!#######################################################################
!-------------------------------------------------------------------
!                Computation of the es values
!
!   Saturation vapor pressure (es) values are computed from
!   equations in the Smithsonian meteorological tables page 350.
!   For temperatures < 0C, sat vapor pres is computed over ice.
!   For temperatures > -20C, sat vapor pres is computed over water.
!   Between -20C and 0C the returned value is blended (over water
!   and over ice).  All sat vapor pres values are returned in pascals.
!
!   Reference:  Smithsonian meteorological tables, page 350.
!-------------------------------------------------------------------

! <FUNCTION NAME="compute_es_1d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
! </FUNCTION>
!function compute_es_1d (tem) result (es)
!real, intent(in) :: tem(:)
!real :: es(size(tem,1))

!es = compute_es_k(tem, TFREEZE)

!end function compute_es_1d
!--------------------------------------------------------

! <FUNCTION NAME="compute_es_0d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
! </FUNCTION>
!function compute_es_0d (tem) result (es)
!real, intent(in) :: tem
!real :: es
!real, dimension(1) :: tem1, es1

!  tem1(1) = tem
!  es1 = compute_es_1d (tem1)
!  es = es1(1)

!end function compute_es_0d

!--------------------------

! <FUNCTION NAME="compute_es_2d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
! </FUNCTION>
!function compute_es_2d (tem) result (es)
!real, intent(in) :: tem(:,:)
!real, dimension(size(tem,1),size(tem,2)) :: es
!integer :: j

!   do j = 1, size(tem,2)
!     es(:,j) = compute_es_1d (tem(:,j))
!   enddo

!end function compute_es_2d

!--------------------------
! <FUNCTION NAME="compute_es_3d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
! </FUNCTION>
!function compute_es_3d (tem) result (es)
!real, intent(in) :: tem(:,:,:)
!real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: es
!integer :: j, k

!   do k = 1, size(tem,3)
!   do j = 1, size(tem,2)
!     es(:,j,k) = compute_es_1d (tem(:,j,k))
!   enddo
!   enddo

!end function compute_es_3d

!#######################################################################
end module sat_vapor_pres_mod
!#######################################################################

! <INFO>

!   <REFERENCE>
!     Smithsonian Meteorological Tables Page 350.
!   </REFERENCE>

!   <BUG>
!     No error checking is done to make sure that the size of the
!     input and output fields match.
!   </BUG>

!   <NOTE>
!     1. <B>Vectorization</B><BR/>
!        To create a vector version the lookup routines need to be modified.
!    The local variables: tmp, del, ind, should be changed to arrays
!    with the same size and order as input array temp.
!
!     2. <B>Construction of the <TT>ES</TT> tables</B><BR/>
!         The tables are constructed using the saturation vapor pressure (<TT>ES</TT>)
!    equations in the Smithsonian tables. The tables are valid between
!    -160C to +100C with increments at 1/10 degree. Between -160C and -20C
!    values of <TT>ES</TT> over ice are used, between 0C and 100C values of<TT> ES</TT>
!    over water are used, between -20C and 0C blended values of <TT>ES</TT>
!    (over water and over ice) are used.
!
!    There are three tables constructed: <TT>ES</TT>, first derivative
!       (<TT>ES'</TT>), and
!    second derivative (<TT>ES</TT>'').  The ES table is constructed directly from
!    the equations in the Smithsonian tables. The <TT>ES</TT>' table is constructed
!    by bracketing temperature values at +/- 0.01 degrees. The <TT>ES</TT>'' table
!    is estimated by using centered differencing of the <TT>ES</TT>' table.
!
!     3. <B>Determination of <TT>es</TT> and <TT>es'</TT> from lookup tables</B><BR/>
!         Values of the saturation vapor pressure (<TT>es</TT>) and the
!    derivative (<TT>es'</TT>) are determined at temperature (T) from the lookup
!    tables (<TT>ES</TT>, <TT>ES'</TT>, <TT>ES''</TT>)
!    using the following formula.
!<PRE>
!    es (T) = ES(t) + ES'(t) * dt + 0.5 * ES''(t) * dt**2
!    es'(T) = ES'(t) + ES''(t) * dt
!
!    where     t = lookup table temperature closest to T
!             dt = T - t
!</PRE>
!
!     4. Internal (private) parameters<BR/>
!       These parameters can be modified to increase/decrease the size/range
!    of the lookup tables.
!<PRE>
!!    tcmin   The minimum temperature (in deg C) in the lookup tables.
!!              [integer, default: tcmin = -160]
!!
!!    tcmax   The maximum temperature (in deg C) in the lookup tables.
!!              [integer, default: tcmin = +100]
!!</PRE>
!!   </NOTE>
!
!!   <TESTPROGRAM NAME="test_sat_vapor_pres">
!<PRE>
!use sat_vapor_pres_mod
!implicit none
!
!integer, parameter :: ipts=500, jpts=100, kpts=50, nloop=1
!real, dimension(ipts,jpts,kpts) :: t,es,esn,des,desn
!integer :: n
!
!! generate temperatures between 120K and 340K
!  call random_number (t)
!  t = 130. + t * 200.
!
!! initialize the tables (optional)
!  call sat_vapor_pres_init
!
!! compute actual es and "almost" actual des
!   es = compute_es  (t)
!  des = compute_des (t)
!
!do n = 1, nloop
!! es and des
!  call lookup_es  (t, esn)
!  call lookup_des (t,desn)
!enddo
!
!! terminate, print deviation from actual
!  print *, 'size=',ipts,jpts,kpts,nloop
!  print *, 'err es  = ', sum((esn-es)**2)
!  print *, 'err des = ', sum((desn-des)**2)
!
!contains
!
!!----------------------------------
!! routine to estimate derivative
!
! function compute_des (tem) result (des)
! real, intent(in) :: tem(:,:,:)
! real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: des,esp,esm
! real, parameter :: tdel = .01
!    esp = compute_es (tem+tdel)
!    esm = compute_es (tem-tdel)
!    des = (esp-esm)/(2*tdel)
! end function compute_des
!!----------------------------------
!
!end program test_sat_vapor_pres
!</PRE>
!   </TESTPROGRAM>
! </INFO>
!> @}
! close documentation grouping
