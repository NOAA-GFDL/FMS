              module sat_vapor_pres_params_mod

character(len=128)   :: version = '$Id: sat_vapor_pres_params.F90,v 13.0 2006/03/28 21:42:57 fms Exp $'
public

! minimum temperature (degC) in lookup table:
integer, parameter :: TCMIN = -160 

! maximum temperature (degC) in lookup table
integer, parameter :: TCMAX =  100 

! table resolution (increments per degree)
integer, parameter :: ESRES =  10  

! lookup table size
integer, parameter :: NSIZE = (tcmax-tcmin)*esres+1   

integer, parameter :: NLIM  = nsize-1
real, parameter    :: TFREEZE = 273.16


              end module sat_vapor_pres_params_mod
