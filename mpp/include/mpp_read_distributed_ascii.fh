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
subroutine MPP_READ_DISTRIBUTED_ASCII_1D_ (unit,fmt,ssize,data,iostat)
  integer, intent(in)               :: unit
  character(*), intent(in)          :: fmt
  integer, intent(in)               :: ssize
  MPP_TYPE_, dimension(:), intent(inout) :: data
  integer, intent(out)              :: iostat

  integer, allocatable :: pelist(:)
  logical              :: is_ioroot=.false.

  if(.not.module_is_initialized) call mpp_error(FATAL,'MPP_READ_DISTRIBUTED_ASCII_1D_:  module not initialized')

  iostat = 0
  call mpp_dist_io_pelist(ssize,pelist)  ! ALLOCATE and create pelist if size of group > 1
  if(.not. ALLOCATED(pelist)) &
           call mpp_error(FATAL,'MPP_READ_DISTRIBUTED_ASCII_1D_:: pelist allocation failed')
  is_ioroot = mpp_is_dist_ioroot(ssize)
  if(is_ioroot) then
    if(trim(fmt)=='*')then
       read(unit,*,iostat=iostat) data
    else
       read(unit,fmt=trim(fmt),iostat=iostat) data
    endif
    if(iostat /= 0) return  ! Calling routine must handle error
  endif

  call mpp_broadcast(data,size(data),pelist(1),pelist)
  deallocate(pelist)  ! Don't forget to deallocate pelist
end subroutine MPP_READ_DISTRIBUTED_ASCII_1D_
