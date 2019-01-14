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

function MPP_CHKSUM_INT_( var, pelist, mask_val )
  integer(LONG_KIND) :: MPP_CHKSUM_INT_
      MPP_TYPE_, intent(in) :: var MPP_RANK_
      integer, optional :: pelist(:)
  MPP_TYPE_, intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     MPP_CHKSUM_INT_ = sum( INT( PACK(var,var/=mask_val),LONG_KIND) )
  else
     MPP_CHKSUM_INT_ = sum(INT(var,LONG_KIND))
  end if

      call mpp_sum( MPP_CHKSUM_INT_, pelist )
      return

    end function MPP_CHKSUM_INT_


!Handles real mask for easier implimentation
! until exists full integer vartypes...
function MPP_CHKSUM_INT_RMASK_( var, pelist, mask_val )
  integer(LONG_KIND) :: MPP_CHKSUM_INT_RMASK_
  MPP_TYPE_, intent(in) :: var MPP_RANK_
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(INT_KIND)::i4tmp(2)=0
  real(FLOAT_KIND)::r4tmp(2)=0
  integer(LONG_KIND) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4,tmpStr5
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT     
  !!! Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val,INT_KIND) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT        
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, _FillValue, pack and mask_val. // &
            & Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. Continuing by using the default MPP_FILL_INT. // &
            & THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  MPP_CHKSUM_INT_RMASK_ = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function MPP_CHKSUM_INT_RMASK_
