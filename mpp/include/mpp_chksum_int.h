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
  integer(KIND(var))::tmpVarP
  integer(KIND(mask_val)) :: tmpFullP,tmpZeroMasked


  if ( KIND(mask_val) == KIND(var)) then
     !same numBytes
     !cast to MPP_TYPE_
     tmpVarP = TRANSFER(mask_val , tmpVarP)
  else if (KIND(mask_val) /=  KIND(var) ) then
     !check if still safe to duck type to lower precision
     tmpFullP = TRANSFER(mask_val,tmpFullP ) ! transfer bits to int mold of same numBytes
     tmpZeroMasked = IBITS(tmpFullP,0,BIT_SIZE(tmpVarP))
     if (tmpFullP /= tmpZeroMasked ) then
        !as an int, mask_val is actually /using/ higher bits than var, and so not a valid mask...
        call mpp_error(FATAL, "mpp_chksum_int.h was called with real mask_val, and mask_val can not be safely cast to int type of var (nonzero high bits).")
     end if
     tmpVarP = INT(tmpFullP, KIND(var) ) ! attempt cast of mask_val as int down to the int precision of var, could do transfer...
  end if

  MPP_CHKSUM_INT_RMASK_ = mpp_chksum(var,pelist,mask_val=tmpVarP)

  return

end function MPP_CHKSUM_INT_RMASK_
