    subroutine MPP_REDUCE_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: peset(3)
#ifdef use_libSMA
!work holds pWrk array + 1 word for symmetric copy of a
      MPP_TYPE_ :: work(SHMEM_REDUCE_MIN_WRKDATA_SIZE+1)
      pointer( ptr, work )
      integer :: words
      character(len=8) :: text
#endif use_libSMA
#ifdef use_libMPI
      MPP_TYPE_ :: work
#endif
      if( .NOT.mpp_initialized )call mpp_error( FATAL, 'MPP_REDUCE: You must first call mpp_init.' )
      if( npes.EQ.1 )return

      if( PRESENT(pelist) )then
          if( size(pelist).EQ.1 )return
          call make_pe_set( pelist, peset )
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
          call mpp_error( FATAL, 'MPP_REDUCE user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      
      work(1) = a
      call SHMEM_REDUCE_( work, work, 1, peset(1), peset(2), peset(3), work(2), sync )
      call SHMEM_BARRIER_ALL()
      a = work(1)
#endif use_libSMA
#ifdef use_libMPI
      if( verbose )call mpp_error( NOTE, 'MPP_REDUCE_: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, 1, MPI_TYPE_, MPI_REDUCE_, peset(1), error )
      a = work
#endif
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_REDUCE_
