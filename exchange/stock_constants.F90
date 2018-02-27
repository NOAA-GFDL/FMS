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

module stock_constants_mod

  use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_sum
  use fms_mod, only : write_version_number
  use time_manager_mod, only : time_type, get_time
  use time_manager_mod, only : operator(+), operator(-)
  use diag_manager_mod, only : register_diag_field,send_data

  implicit none

  ! Include variable "version" to be written to log file.
#include<file_version.h>

  integer,public,    parameter                :: NELEMS=3
  integer,           parameter                :: NELEMS_report=3
  integer,public,    parameter                :: ISTOCK_WATER=1, ISTOCK_HEAT=2, ISTOCK_SALT=3
  integer,public,    parameter                :: ISTOCK_TOP=1, ISTOCK_BOTTOM=2, ISTOCK_SIDE=3
  integer,public                              :: stocks_file
  ! Stock related stuff
  ! Shallow (no constructor) data structures holding the starting stock values (per PE) and
  ! flux integrated increments at present time.

  integer, parameter :: NSIDES  = 3         ! top, bottom, side
 
  type stock_type ! per PE values
     real  :: q_start = 0.0    ! total stocks at start time
     real  :: q_now   = 0.0    ! total stocks at time t    
     
     ! The dq's below are the stocks increments at the present time
     ! delta_t * surf integr of flux
     ! one for each side (ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE)
     real  :: dq(NSIDES)    = 0.0    ! stock increments at present time on the Ice   grid     
     real  :: dq_IN(NSIDES) = 0.0    ! stock increments at present time on the Ocean grid       
  end type stock_type

  type(stock_type), save, dimension(NELEMS) :: Atm_stock, Ocn_stock, Lnd_stock, Ice_stock
  type(time_type), save :: init_time

  public stocks_report
  public stocks_report_init
  public stocks_set_init_time

  integer,private,   parameter                :: NCOMPS=4
  integer,private,   parameter                :: ISTOCK_ATM=1,ISTOCK_LND=2,ISTOCK_ICE=3,ISTOCK_OCN=4
  character(len=3) , parameter, dimension(NCOMPS)  :: COMP_NAMES=(/'ATM', 'LND', 'ICE', 'OCN'/)


  character(len=5) , parameter, dimension(NELEMS)  :: STOCK_NAMES=(/'water', 'heat ', 'salt '/)
  character(len=12), parameter, dimension(NELEMS)  :: STOCK_UNITS=(/'[Kg]    ','[Joules]','[Kg]    '/)


contains


    subroutine stocks_report_init(Time)
    type(time_type)               , intent(in) :: Time

    character(len=80) :: formatString,space
    integer :: i,s
    real, dimension(NELEMS) :: val_atm, val_lnd, val_ice, val_ocn

    ! Write the version of this file to the log file
    call write_version_number('STOCK_CONSTANTS_MOD', version)

    do i = 1, NELEMS_report
       val_atm(i) = Atm_stock(i)%q_start
       val_lnd(i) = Lnd_stock(i)%q_start
       val_ice(i) = Ice_stock(i)%q_start
       val_ocn(i) = Ocn_stock(i)%q_start
       call mpp_sum(val_atm(i))
       call mpp_sum(val_lnd(i))
       call mpp_sum(val_ice(i))
       call mpp_sum(val_ocn(i))
    enddo



    if(mpp_pe() == mpp_root_pe()) then
!       earth_area = 4.*PI*Radius**2

       write(stocks_file,*) '================Stocks Report Guide====================================='
       write(stocks_file,*) ' '
       write(stocks_file,*) 'S(t) = Total amount     of a tracer in the component model at time t.'
       write(stocks_file,*) '       Calculated via the component model itself.'
       write(stocks_file,*) ' '
       write(stocks_file,*) 'F(t) = Cumulative input of a tracer to the component model at time t.'
       write(stocks_file,*) '       Calculated via interchange of fluxes with other component models.'
       write(stocks_file,*) ' '
       write(stocks_file,*) 'S(t) - S(0) = Cumulative increase of the component stocks at time t'
       write(stocks_file,*) '              Calculated by the component itself.'
       write(stocks_file,*) ' '
       write(stocks_file,*) 'In a conserving component F(t)=S(t)-S(0) to within numerical accuracy.'
       write(stocks_file,*) ' '
       write(stocks_file,*) 'Component Model refers to one of OCN, ATM, LND or ICE'
       write(stocks_file,*) ''
       write(stocks_file,*) 'NOTE: When use_lag_fluxes=.true. is used in coupler, the ocean stocks '
       write(stocks_file,*) '      calculations are in error by an order which scales as the inverse'
       write(stocks_file,*) '      of the number of time steps.'
       write(stocks_file,*) ' '
       write(stocks_file,*) '======================================================================='       


       write(stocks_file,*) '======================Initial Stock S(0)==============================='       
!The following produces  formatString='(5x,a,a,12x,a,a, 9x)' but is general to handle more elements         
       formatString= '(5x'
       do i=1,NELEMS_report
          s = 25-len_trim(STOCK_NAMES(i))-len_trim(STOCK_UNITS(i))
          write(space,'(i2)') s 
          formatString= trim(formatString)//',a,a,'//trim(space)
          formatString= trim(formatString)//trim('x') 
       enddo
       formatString= trim(formatString)//')'
       
       write(stocks_file,formatString) (trim(STOCK_NAMES(i)),trim(STOCK_UNITS(i)), i=1,NELEMS_report)

!The following produces  formatString=' (a,x,es22.15,3x,es22.15,3x)' but is general to handle more elements
       formatString= '(a,x'
       do i=1,NELEMS_report
          write(space,'(i2)') s 
          formatString= trim(formatString)//',es22.15,3x'
       enddo
       formatString= trim(formatString)//')'
       
      
       write(stocks_file,formatString) 'ATM', (val_atm(i), i=1,NELEMS_report) 
       write(stocks_file,formatString) 'LND', (val_lnd(i), i=1,NELEMS_report)
       write(stocks_file,formatString) 'ICE', (val_ice(i), i=1,NELEMS_report)
       write(stocks_file,formatString) 'OCN', (val_ocn(i), i=1,NELEMS_report)

       write(stocks_file,*) '========================================================================'       
       write(stocks_file,'(a)'  ) ' '!blank line

    end if
    
    call stocks_set_init_time(Time)

  end subroutine stocks_report_init


  subroutine stocks_report(Time)
    type(time_type)               , intent(in) :: Time

    type(time_type) :: timeSinceStart
    type(stock_type) :: stck
    real, dimension(NCOMPS) :: f_value, f_ice_grid, f_ocn_grid, f_ocn_btf, q_start, q_now,c_value
    character(len=80) :: formatString
    integer :: iday0, isec0, iday, isec, hours
    real    :: days
    integer :: diagID , comp,elem,i
    integer, parameter :: initID = -2 ! initial value for diag IDs. Must not be equal to the value 
    ! that register_diag_field returns when it can't register the filed -- otherwise the registration 
    ! is attempted every time this subroutine is called

    integer, dimension(NCOMPS,NELEMS), save :: f_valueDiagID = initID
    integer, dimension(NCOMPS,NELEMS), save :: c_valueDiagID = initID
    integer, dimension(NCOMPS,NELEMS), save :: fmc_valueDiagID = initID
    integer, dimension(NCOMPS,NELEMS), save :: f_lostDiagID = initID

    real :: diagField
    logical :: used
    character(len=30) :: field_name, units

    if(mpp_pe()==mpp_root_pe()) then
       call get_time(init_time, isec0, iday0)
       call get_time(Time, isec, iday)
       
       hours = iday*24 + isec/3600 - iday0*24 - isec0/3600
       days  = hours/24.  
       write(stocks_file,*) '==============================================='
       write(stocks_file,'(a,f12.3)') 't = TimeSinceStart[days]= ',days 
       write(stocks_file,*) '==============================================='
    endif

    do elem = 1,NELEMS_report

       do comp = 1,NCOMPS

          if(comp == ISTOCK_ATM) stck = Atm_stock(elem)
          if(comp == ISTOCK_LND) stck = Lnd_stock(elem)
          if(comp == ISTOCK_ICE) stck = Ice_stock(elem)
          if(comp == ISTOCK_OCN) stck = Ocn_stock(elem)


          f_ice_grid(comp) = sum(stck%dq)
          f_ocn_grid(comp) = sum(stck%dq_IN)
          f_ocn_btf(comp)  = stck%dq_IN( ISTOCK_BOTTOM )

          q_start(comp) = stck%q_start
          q_now(comp)   = stck%q_now 

          call mpp_sum(f_ice_grid(comp))
          call mpp_sum(f_ocn_grid(comp))
          call mpp_sum(f_ocn_btf(comp))
          call mpp_sum(q_start(comp))
          call mpp_sum(q_now(comp))

          c_value(comp) = q_now(comp) - q_start(comp)

          if(mpp_pe() == mpp_root_pe()) then

             if(f_valueDiagID(comp,elem) == initID) then
                field_name = trim(COMP_NAMES(comp)) // trim(STOCK_NAMES(elem))
                field_name  = trim(field_name) // 'StocksChange_Flux'
                units = trim(STOCK_UNITS(elem))
                f_valueDiagID(comp,elem) = register_diag_field('stock_print', field_name, Time, &
                     units=units)
             endif

             if(c_valueDiagID(comp,elem) == initID) then
                field_name = trim(COMP_NAMES(comp)) // trim(STOCK_NAMES(elem))
                field_name = trim(field_name) // 'StocksChange_Comp'
                units = trim(STOCK_UNITS(elem))
                c_valueDiagID(comp,elem) = register_diag_field('stock_print', field_name, Time, &
                     units=units)
             endif

             if(fmc_valueDiagID(comp,elem) == initID) then
                field_name = trim(COMP_NAMES(comp)) // trim(STOCK_NAMES(elem))
                field_name = trim(field_name) // 'StocksChange_Diff'
                units = trim(STOCK_UNITS(elem))
                fmc_valueDiagID(comp,elem) = register_diag_field('stock_print', field_name, Time, &
                     units=units)
             endif

             f_value(comp) = f_ice_grid(comp)

             if(comp == ISTOCK_OCN) then

                f_value(comp) = f_ocn_grid(comp)

                if(f_lostDiagID(comp,elem) == initID) then
                   field_name = trim(COMP_NAMES(comp)) // trim(STOCK_NAMES(elem))
                   field_name = trim(field_name) // 'StocksExchangeLost'
                   units = trim(STOCK_UNITS(elem))
                   f_lostDiagID(comp,elem) = register_diag_field('stock_print', field_name, Time, &
                        units=units)
                endif

                DiagID=f_lostDiagID(comp,elem)
                diagField = f_ice_grid(comp) - f_ocn_grid(comp)
                if (DiagID > 0)  used = send_data(DiagID, diagField, Time)

             endif


             DiagID=f_valueDiagID(comp,elem)
             diagField = f_value(comp)
             if (DiagID > 0)  used = send_data(DiagID, diagField, Time)
             DiagID=c_valueDiagID(comp,elem)
             diagField = c_value(comp)
             if (DiagID > 0)  used = send_data(DiagID, diagField, Time)
             DiagID=fmc_valueDiagID(comp,elem)
             diagField = f_value(comp)-c_value(comp)
             if (DiagID > 0)  used = send_data(DiagID, diagField, Time)


             !             formatString = '(a,a,a,i16,2x,es22.15,2x,es22.15,2x,es22.15,2x,es22.15,2x,es22.15,2x,es22.15)'
             !
             !             write(stocks_file,formatString) trim(COMP_NAMES(comp)),STOCK_NAMES(elem),STOCK_UNITS(elem) &
             !                  ,hours, q_now, q_now-q_start, f_value, f_value - (q_now - q_start), (f_value - (q_now - q_start))/q_start


          endif
       enddo


       if(mpp_pe()==mpp_root_pe()) then
!          write(stocks_file,'(a)'  ) ' '!blank line
!          write(stocks_file,'(a,f12.3)') 't = TimeSinceStart[days]= ',days 
!          write(stocks_file,'(a)'  )   ' '!blank line
!          write(stocks_file,'(a,30x,a,20x,a,20x,a,20x,a)') 'Component ','ATM','LND','ICE','OCN'
!          write(stocks_file,'(55x,a,20x,a,20x,a,20x,a)')  'ATM','LND','ICE','OCN'
!          write(stocks_file,'(a,f12.3,12x,a,20x,a,20x,a,20x,a)') 't = TimeSinceStart[days]= ',days,'ATM','LND','ICE','OCN'

          write(stocks_file,'(a,a,40x,a,20x,a,20x,a,20x,a)') 'Stocks of ',trim(STOCK_NAMES(elem)),'ATM','LND','ICE','OCN'
          formatString = '(a,a,2x,es22.15,2x,es22.15,2x,es22.15,2x,es22.15)'

          write(stocks_file,formatString) 'Total =S(t)               ',STOCK_UNITS(elem),&
               ( q_now(i), i=1,NCOMPS)
          write(stocks_file,formatString) 'Change=S(t)-S(0)          ',STOCK_UNITS(elem),&
               ( q_now(i)-q_start(i), i=1,NCOMPS)
          write(stocks_file,formatString) 'Input =F(t)               ',STOCK_UNITS(elem),&
               ( f_value(i), i=1,NCOMPS)
          write(stocks_file,formatString) 'Diff  =F(t) - (S(t)-S(0)) ',STOCK_UNITS(elem),&
               ( f_value(i) - c_value(i), i=1,NCOMPS)                           
          write(stocks_file,formatString) 'Error =Diff/S(0)          ','[NonDim]    ', &
               ((f_value(i) - c_value(i))/(1+q_start(i)), i=1,NCOMPS)  !added 1 to avoid div by zero. Assuming q_start large          

          write(stocks_file,'(a)'  ) ' '!blank line
          formatString = '(a,a,a,6x,es22.15)'
          write(stocks_file,formatString) 'Lost Stocks in the exchange between Ice and Ocean ',trim(STOCK_NAMES(elem)),trim(STOCK_UNITS(elem)),  &
               f_ice_grid(ISTOCK_OCN) - f_ocn_grid(ISTOCK_OCN) + f_ocn_btf(ISTOCK_OCN)

          write(stocks_file,'(a)') ' ' !blank line  
          write(stocks_file,'(a)') ' ' !blank line

       endif
    enddo

  end subroutine stocks_report

  subroutine stocks_set_init_time(Time)
    type(time_type)     , intent(in) :: Time
    init_time = Time
    
  end subroutine stocks_set_init_time

end module stock_constants_mod
