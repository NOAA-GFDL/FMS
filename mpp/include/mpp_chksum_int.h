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
  real(KIND(var)) :: tmpVarP
  real(KIND(mask_val)) :: tmpMaskP
  integer(KIND(var)) :: mask_val_bitcastToVarKind
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=256) :: errStr

  if ( KIND(mask_val) == KIND(var)) then !sameKind
     !cast to MPP_TYPE_
     mask_val_bitcastToVarKind = TRANSFER(mask_val , mask_val_bitcastToVarKind)
  else ! KIND of mask_val and var are different
    !check to see if can lossless cast to the precision KIND of var and back
    tmpVarP  = REAL( mask_val , KIND=KIND(var) ) ! cast to precision of var
    tmpMaskP = REAL( tmpVarP  , KIND=KIND(mask_val) ) !cast back to precision of mask_val

    if (  tmpMaskP == mask_val ) then !safeCast
       !We can use cast of mask_val
       mask_val_bitcastToVarKind = TRANSFER( tmpVarP , mask_val_bitcastToVarKind )
    else 
      ! construct detailed errStr
      errStr = "mpp_chksum: mpp_chksum_i" 
      write(unit=tmpStr1,fmt="(I1)") KIND(var) 
      write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
      errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL(" 
      write(unit=tmpstr3,fmt="(I1)") KIND(mask_val) 
      errStr = errStr // tmpstr3 // ") mask_val="
      write(unit=tmpstr4,fmt=*) mask_val
      errStr = errStr // trim(tmpstr4) // ", mask cannot be safely cast. Check _FillValue and mask_val. Hint: Try MPP_FILL_{INT,FLOAT,DOUBLE}."
      call mpp_error(FATAL, trim(errStr) )
    end if ! safeCast
  end if ! sameKind

  ! proceed as integer checksum using bitcast mask_val
  MPP_CHKSUM_INT_RMASK_ = mpp_chksum(var,pelist,mask_val=mask_val_bitcastToVarKind)

  return

end function MPP_CHKSUM_INT_RMASK_
