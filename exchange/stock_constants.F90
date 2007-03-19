module stock_constants_mod

  character(len=128), parameter :: version = '$Id: stock_constants.F90,v 14.0 2007/03/15 22:39:07 fms Exp $'

  integer,           parameter                :: ISTOCK_WATER=1, ISTOCK_HEAT=2
  character(len=16), parameter, dimension(2)  :: STOCK_NAMES=(/'water', 'heat '/)
  integer,           parameter                :: ISTOCK_TOP=1, ISTOCK_BOTTOM=2, ISTOCK_SIDE=3

end module stock_constants_mod
