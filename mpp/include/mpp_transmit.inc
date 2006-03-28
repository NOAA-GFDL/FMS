
! -*-f90-*-
! $Id: mpp_transmit.inc,v 13.0 2006/03/28 21:41:46 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_TRANSMIT_SCALAR_( put_data, to_pe, get_data, from_pe, plen, glen)
      integer, intent(in) :: to_pe, from_pe
      MPP_TYPE_, intent(in)  :: put_data
      MPP_TYPE_, intent(out) :: get_data
      integer, optional, intent(in) :: plen, glen
      integer                       :: put_len, get_len
      MPP_TYPE_ :: put_data1D(1), get_data1D(1)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call MPP_TRANSMIT_ ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      put_data1D(1) = put_data
      call MPP_TRANSMIT_ ( put_data1D, 1, to_pe, get_data1D, 1, from_pe )
      get_data = get_data1D(1)
#endif
      return
    end subroutine MPP_TRANSMIT_SCALAR_

    subroutine MPP_TRANSMIT_2D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_2D_

    subroutine MPP_TRANSMIT_3D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_3D_

    subroutine MPP_TRANSMIT_4D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_4D_

    subroutine MPP_TRANSMIT_5D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_5D_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_RECV_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(*)
      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_

    subroutine MPP_SEND_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(*)
      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_

    subroutine MPP_RECV_SCALAR_( get_data, from_pe, glen )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      MPP_TYPE_, intent(out) :: get_data
      integer, optional, intent(in) :: glen
      integer                       :: get_len
      MPP_TYPE_ :: get_data1D(1)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, get_data1D )
      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe )
#else
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, 1, from_pe )
      get_data = get_data1D(1)
#endif
    end subroutine MPP_RECV_SCALAR_

    subroutine MPP_SEND_SCALAR_( put_data, to_pe, plen)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      MPP_TYPE_, intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer                       :: put_len
      MPP_TYPE_ :: put_data1D(1)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE )
#else
      put_data1D(1) = put_data
      call mpp_transmit( put_data1D, 1, to_pe, dummy, 1, NULL_PE )
#endif
    end subroutine MPP_SEND_SCALAR_

    subroutine MPP_RECV_2D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:)
      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_2D_

    subroutine MPP_SEND_2D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:)
      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_2D_

    subroutine MPP_RECV_3D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_3D_

    subroutine MPP_SEND_3D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:)
      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_3D_

    subroutine MPP_RECV_4D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_4D_

    subroutine MPP_SEND_4D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_4D_

    subroutine MPP_RECV_5D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_5D_

    subroutine MPP_SEND_5D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_5D_
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_SCALAR_( data, from_pe, pelist )
      MPP_TYPE_, intent(inout) :: data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(1)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )

      ptr = LOC(data)
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )
#else
      data1D(1) = data
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )
      data = data1D(1)
#endif
      return
    end subroutine MPP_BROADCAST_SCALAR_
    
    subroutine MPP_BROADCAST_2D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_2D_
    
    subroutine MPP_BROADCAST_3D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_3D_
    
    subroutine MPP_BROADCAST_4D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_4D_
    
    subroutine MPP_BROADCAST_5D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_5D_
