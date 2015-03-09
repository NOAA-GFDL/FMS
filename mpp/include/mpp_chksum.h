function MPP_CHKSUM_( var, pelist , mask_val)
!mold is a dummy array to be used by TRANSFER()
!must be same TYPE as result
!result is LONG_KIND, which will actually be int ifdef no_8byte_integers
  !optional mask_val is masked away in checksum_int.h function via PACK()
  integer(LONG_KIND) :: MPP_CHKSUM_
  integer(LONG_KIND) :: mold(1)
      MPP_TYPE_, intent(in) :: var MPP_RANK_
      integer, intent(in), optional :: pelist(:)
  MPP_TYPE_, intent(in),optional :: mask_val

  if ( PRESENT(mask_val) ) then
     MPP_CHKSUM_ = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      MPP_CHKSUM_ = mpp_chksum( TRANSFER(var,mold), pelist )
  end if
  
      return
    end function MPP_CHKSUM_
