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
  integer(KIND(var))::itmpVarP(2)
  real(KIND(var))::rtmpVarP(2)
  !!integer(KIND(mask_val)) :: tmpFullP,tmpZeroMasked
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=256) :: errStr


  ! typically should catch real "default_fill" .OR. mask_val == MPP_FILL_FLOAT
  if (mask_val == MPP_FILL_DOUBLE ) then
     !use MPP_FILL_INT
     imask_val = MPP_FILL_INT
  else if ( CEILING(mask_val,INT_KIND) == MPP_FILL_INT ) then
     ! NETCDF fill values proveide CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
     imask_val = MPP_FILL_INT
  else if ( KIND(mask_val) == KIND(var) ) then
     ! same size see if previously ducktyped integere via transfer
     itmpVarP(1) = TRANSFER(mask_val , itmpVarP(1))
     if ( itmpVarP(1) == MPP_FILL_INT ) imask_val = MPP_FILL_INT
  else if ( KIND(mask_val)==DOUBLE_KIND .AND. KIND(var)==INT_KIND ) then
     itmpVarP = TRANSFER(mask_val, itmpVarP)
     ! did we pack int fill
     if ( ANY(itmpVarP == MPP_FILL_INT ) .OR. ANY( itmpVarP == CEILING(MPP_FILL_DOUBLE,LONG_KIND) ) ) then
        imask_val = MPP_FILL_INT
     else ! did we pack real fill
        rtmpVarP = TRANSFER(mask_val, rtmpVarP)
        if ( ANY(rtmpVarP == MPP_FILL_DOUBLE ) .OR. ANY(CEILING(rtmpVarP,INT_KIND) == MPP_FILL_INT) ) then
           imask_val = MPP_FILL_INT
        end if
     end if
! here be dragons, we should not be doing either of these things
  else if ( KIND(mask_val)==FLOAT_KIND .AND. KIND(var)==LONG_KIND ) then
     call mpp_error(NOTE, "Your have called mpp_chksum with a var(LONG_KIND) and mask_val(FLOAT_KIND) // &
          and are not using the supported MPP_FILL_ parameters. Taking CEILING(mask_val,LONG_KING),// & 
          please confirm this is what you want and use MPP_FILL_ to make this message go away.")
     imask_val = CEILING(mask_val,KIND=LONG_KIND)
  else 
     ! construct detailed errStr
      errStr = "mpp_chksum: mpp_chksum_i" 
      write(unit=tmpStr1,fmt="(I1)") KIND(var) 
      write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
      errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL(" 
      write(unit=tmpstr3,fmt="(I1)") KIND(mask_val) 
      errStr = errStr // tmpstr3 // ") mask_val="
      write(unit=tmpstr4,fmt=*) mask_val
      errStr = errStr // trim(tmpstr4) // ", andthis is not a supported cast. Check _FillValue and mask_val. Hint: Try MPP_FILL_{INT,FLOAT,DOUBLE}. Defaulting to MPP_FILL_INT."
      call mpp_error(NOTE, trim(errStr) ) 

      imask_val = MPP_FILL_INT
   end if

  MPP_CHKSUM_INT_RMASK_ = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function MPP_CHKSUM_INT_RMASK_
