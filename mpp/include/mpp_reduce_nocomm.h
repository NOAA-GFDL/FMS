    subroutine MPP_REDUCE_0D_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine MPP_REDUCE_0D_

    subroutine MPP_REDUCE_1D_( a, length, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine MPP_REDUCE_1D_

