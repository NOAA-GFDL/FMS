
! -*-f90-*-
! $Id: mpp_sum.inc,v 13.0 2006/03/28 21:41:40 fms Exp $

!#######################################################################

    subroutine MPP_SUM_SCALAR_( a, pelist )
!sums array a when only first element is passed: this routine just converts to a call to MPP_SUM_
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_( b, 1, pelist )
      a = b(1)
      return
    end subroutine MPP_SUM_SCALAR_

!#######################################################################
    subroutine MPP_SUM_2D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_2D_

!#######################################################################
    subroutine MPP_SUM_3D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_3D_

!#######################################################################
    subroutine MPP_SUM_4D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_4D_

!#######################################################################
    subroutine MPP_SUM_5D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_5D_
