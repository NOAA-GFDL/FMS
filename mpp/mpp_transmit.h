    subroutine MPP_TRANSMIT_( put_data, put_len, to_pe, get_data, get_len, from_pe )
!a message-passing routine intended to be reminiscent equally of both MPI and SHMEM

!put_data and get_data are contiguous MPP_TYPE_ arrays

!at each call, your put_data array is put to   to_pe's get_data
!              your get_data array is got from from_pe's put_data
!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get

!special PE designations:
!      NULL_PE: to disable a put or a get (e.g at boundaries)
!      ANY_PE:  if remote PE for the put or get is to be unspecific
!      ALL_PES: broadcast and collect operations (collect not yet implemented)

!ideally we would not pass length, but this f77-style call performs better (arrays passed by address, not descriptor)
!further, this permits <length> contiguous words from an array of any rank to be passed (avoiding f90 rank conformance check)

!caller is responsible for completion checks (mpp_sync_self) before and after

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(*)
      MPP_TYPE_, intent(out) :: get_data(*)
#ifdef use_libSMA
      integer :: np
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      MPP_TYPE_ :: remote_data(get_len)
      pointer( ptr_remote_data, remote_data )
      MPP_TYPE_ :: broadcast_data(get_len)
      pointer( ptr, broadcast_data )
      integer :: words
      character(len=8) :: text
#endif

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return
      
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
#ifdef use_libSMA
!send data pointer to to_pe
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( status(to_pe), MPP_WAIT )
          status(to_pe) = MPP_WAIT !prohibit puts to to_pe until it has retrieved this message
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          data_loc = LOC(put_data)
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INTEGER_PUT( mpp_from_pe, pe, 1, to_pe )
          call SHMEM_PUT8( remote_data_loc(pe), data_loc, 1, to_pe )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_SEND, put_len*MPP_TYPE_BYTELEN_ )
#endif use_libSMA
#ifdef use_libMPI
!use non-blocking sends
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_ISEND( put_data, put_len, MPI_TYPE_, to_pe, tag, MPI_COMM_WORLD, request(to_pe), error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_SEND, put_len*MPP_TYPE_BYTELEN_ )
#endif

      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len )call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
#ifdef use_libSMA
          ptr = LOC(mpp_stack)
          words = size(broadcast_data)*size(transfer(broadcast_data(1),word))
          if( words.GT.mpp_stack_size )then
              write( text, '(i8)' )words
              call mpp_error( FATAL, 'MPP_TRANSMIT user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
          end if
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( npes.GT.1 )then
              broadcast_data(1:get_len) = put_data(1:get_len)
              call mpp_sync()
#ifdef _CRAYT90
              call SHMEM_UDCFLUSH !invalidate data cache
#endif
              call SHMEM_BROADCAST_( broadcast_data, broadcast_data, get_len, from_pe, 0,0,npes, sync )
              call mpp_sync()
              get_data(1:get_len) = broadcast_data(1:get_len)
          end if
          if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, get_len*MPP_TYPE_BYTELEN_ )
#endif
!SHMEM_BROADCAST does not broadcast to itself:
! need copy if get_data and put_data are not the same.
!contrariwise MPI_BCAST only copies itself:
! need copy prior to calling MPI_BCAST.
!thus this copy operation is placed inbetween.
          if( pe.EQ.from_pe )then
#ifdef use_CRI_pointers
!dir$ IVDEP
              if( LOC(get_data).NE.LOC(put_data) ) &
#endif
                   get_data(1:put_len) = put_data(1:put_len)
          end if
#ifdef use_libMPI
          if( npes.GT.1 )then
             if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
             call MPI_BCAST( get_data, put_len, MPI_TYPE_, from_pe, MPI_COMM_WORLD, error )
             if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, put_len*MPP_TYPE_BYTELEN_ )
          endif
#endif
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets
#ifdef use_libSMA
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe along with to_pe=ANY_PE.' )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_GET_( get_data, put_data, get_len, from_pe )
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe that you have retrieved this message
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          return
#endif
#ifdef use_libMPI
!...but you cannot have a pure get with MPI
          call mpp_error( FATAL, 'MPP_TRANSMIT: you cannot transmit to ANY_PE using MPI.' )
#endif

      else if( to_pe.NE.NULL_PE )then	!no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( debug )write( stderr,* )'pe, from_pe, remote_data_loc(from_pe)=', pe, from_pe, remote_data_loc(from_pe)
          call SHMEM_INT8_WAIT( remote_data_loc(from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(from_pe)
          remote_data_loc(from_pe) = MPP_WAIT !reset
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe we have retrieved the location
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
!dir$ IVDEP
          if( LOC(get_data).NE.LOC(remote_data) )get_data(1:get_len) = remote_data(1:get_len)
#else
          call SHMEM_GET_( get_data, remote_data, get_len, from_pe )
#endif
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#elif  use_libMPI
!receive from from_pe
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_RECV( get_data, get_len, MPI_TYPE_, from_pe, MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#else !neither use_libSMA nor use_libMPI
          if( pe.EQ.from_pe )then
#ifdef use_CRI_pointers
!dir$ IVDEP
              if( LOC(get_data).NE.LOC(put_data) ) &
#endif
                   get_data(1:put_len) = put_data(1:put_len)
          end if
#endif

      else if( from_pe.EQ.ANY_PE )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
!since we don't know which PE is sending us data, we wait for remote PE to send us its ID
!this is only required for !CRAYPVP  && !sgi_mipspro, but is done there too, so that we can send put_is_done back.
          call shmem_integer_wait( mpp_from_pe, ANY_PE )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( remote_data_loc(mpp_from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(mpp_from_pe)
          remote_data_loc(mpp_from_pe) = MPP_WAIT !reset
          call SHMEM_PUT8( status(pe), MPP_READY, 1, mpp_from_pe ) !tell mpp_from_pe we have retrieved the location
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
!dir$ IVDEP
          if( LOC(get_data).NE.LOC(remote_data) )get_data(1:get_len) = remote_data(1:get_len)
#else
          call SHMEM_GET_( get_data, remote_data, get_len, mpp_from_pe )
#endif
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          mpp_from_pe = ANY_PE   !reset
#endif use_libSMA
#ifdef use_libMPI
!receive from MPI_ANY_SOURCE
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_RECV( get_data, get_len, MPI_TYPE_, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#endif

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if
      return
    end subroutine MPP_TRANSMIT_

    subroutine MPP_TRANSMIT_SCALAR_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data
      MPP_TYPE_, intent(out) :: get_data
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
#else
      call mpp_error( FATAL, 'MPP_TRANSMIT_SCALAR_ requires Cray pointers.' )
#endif
      call MPP_TRANSMIT_ ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      return
    end subroutine MPP_TRANSMIT_SCALAR_

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

    subroutine MPP_RECV_SCALAR_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data
      MPP_TYPE_ :: get_data1D(get_len)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, get_data1D )
      ptr = LOC(get_data)
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe )
#else
      call mpp_error( FATAL, 'MPP_RECV_SCALAR_ currently requires CRI pointers.' )
#endif
    end subroutine MPP_RECV_SCALAR_

    subroutine MPP_SEND_SCALAR_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data
      MPP_TYPE_ :: put_data1D(put_len)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE )
#else
      call mpp_error( FATAL, 'MPP_SEND_SCALAR_ currently requires CRI pointers.' )
#endif
    end subroutine MPP_SEND_SCALAR_
    
