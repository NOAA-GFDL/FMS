module stock_constants_mod

  character(len=128), parameter :: version = '$Id: stock_constants.F90,v 1.1.2.7 2006/11/13 21:02:41 ap Exp $'

  integer,           parameter                :: ISTOCK_WATER=1, ISTOCK_HEAT=2
  character(len=16), parameter, dimension(2)  :: STOCK_NAMES=(/'water', 'heat '/)
  integer,           parameter                :: ISTOCK_TOP=1, ISTOCK_BOTTOM=2, ISTOCK_SIDE=3

end module stock_constants_mod
