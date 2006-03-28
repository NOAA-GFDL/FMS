
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

!  lookup table limits in degK:
real, parameter    :: TMIN = real(TCMIN) + TFREEZE  
real, parameter    :: TMAX = real(TCMAX) + TFREEZE  

real, parameter    :: DTRES = 1./real(ESRES)
real, parameter    :: DTINV = real(ESRES)
real, parameter    :: TEPS = 1./real(2*ESRES)

