function MPP_CHKSUM_( var, pelist, mask_val )
!mold is a dummy array to be used by TRANSFER()
!must be same TYPE as result
!result is LONG_KIND, which will actually be int ifdef no_8byte_integers
  !mold and mask_val must be same numBytes, otherwise undefined behavior
      integer(LONG_KIND) :: MPP_CHKSUM_
      MPP_TYPE_, intent(in) :: var
      integer, intent(in), optional :: pelist(:)
      integer(LONG_KIND) :: mold(1)
  MPP_TYPE_, intent(in), optional :: mask_val
      pointer( p, mold )

      p = LOC(var)

  if ( PRESENT(mask_val) ) then
     MPP_CHKSUM_ = mpp_chksum( mold, pelist, TRANSFER(mask_val, mold(1)) )
  else
      MPP_CHKSUM_ = mpp_chksum( mold, pelist )
  end if
      return
    end function MPP_CHKSUM_
