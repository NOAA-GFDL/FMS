    subroutine MPP_SUM_( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
!we are using f77-style call: array passed by address and not descriptor; further, the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a(*)
      integer :: peset(3)
#ifdef use_libSMA
!first <length> words are array, rest are pWrk
      MPP_TYPE_ :: work(length+length/2+1+SHMEM_REDUCE_MIN_WRKDATA_SIZE)
      pointer( ptr, work )
      integer :: words
      character(len=8) :: text
#else
      MPP_TYPE_ :: work(length)
#endif

      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
      end if

      if( PRESENT(pelist) )then
          call make_pe_set(pelist,peset)
      else
#ifdef use_libSMA
          peset(1) = 0
          peset(2) = 0
          peset(3) = npes
#endif
#ifdef use_libMPI
          peset(1) = MPI_COMM_WORLD
#endif
      end if
      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
!allocate space from the stack for pwrk and b
      ptr = LOC(mpp_stack)
      words = size(work)*size(transfer(work(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_SUM user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      work(1:length) = a(1:length)
      call SHMEM_BARRIER_ALL()
      call SHMEM_SUM_( work, work, length, peset(1), peset(2), peset(3), work(length+1), sync )
#endif use_libSMA
#ifdef use_libMPI
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, length, MPI_TYPE_, MPI_SUM, peset(1), error )
#endif
      a(1:length) = work(1:length)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_SUM_

    subroutine MPP_SUM_SCALAR_( a, length, pelist )
!sums array a when only first element is passed: this routine just converts to a call to MPP_SUM_
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a
      MPP_TYPE_ :: b(length)
#ifdef use_CRI_pointers
      pointer( ptr, b )
      ptr = loc(a)
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_( b, length, pelist )
#else
      call mpp_error( FATAL, 'MPP_SUM_SCALAR_: currently requires CRI pointers.' )
#endif
      return
    end subroutine MPP_SUM_SCALAR_
