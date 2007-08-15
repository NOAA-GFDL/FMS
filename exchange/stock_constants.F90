module stock_constants_mod

  character(len=128), parameter :: version = '$Id: stock_constants.F90,v 15.0 2007/08/14 04:13:43 fms Exp $'

  integer,           parameter                :: NELEMS=2
  integer,           parameter                :: ISTOCK_WATER=1, ISTOCK_HEAT=2
  character(len=5) , parameter, dimension(2)  :: STOCK_NAMES=(/'water', 'heat '/)
  character(len=12), parameter, dimension(2)  :: STOCK_UNITS=(/'[Kg/m^2]    ','[Joules/m^2]'/)
  integer,           parameter                :: ISTOCK_TOP=1, ISTOCK_BOTTOM=2, ISTOCK_SIDE=3
  integer                                     :: stocks_file
end module stock_constants_mod
