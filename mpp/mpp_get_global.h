    subroutine MPP_GET_GLOBAL_( domain, local_field, global_field )
!given an array defined on the data domain <domain>, constructs a global field
!USE WITH CARE! a global 3D array could occupy a lot of memory!
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local_field(domain%x%data%start_index:,domain%y%data%start_index: MPP_EXTRA_INDICES_ )
      MPP_TYPE_, intent(out) :: global_field(domain%x%global%start_index:,domain%y%global%start_index: MPP_EXTRA_INDICES_ )
!      type(domain2D), allocatable, target :: global_domain(:)
      integer :: i
#ifdef use_CRI_pointers
      MPP_TYPE_ :: global_field_loc
      pointer( ptr, global_field_loc )

      ptr = LOC(global_field)
#else
      call mpp_error( FATAL, 'MPP_GET_GLOBAL_: currently requires Cray pointers.' )
#endif
      if( .NOT.mpp_domains_initialized )call mpp_error( FATAL, 'MPP_GET_GLOBAL_3D_: You must first call mpp_domains_init.' )
      if( size(local_field,1).NE.domain%x%data%size .OR. &
          size(local_field,2).NE.domain%y%data%size .OR. &
          size(global_field,1).NE.domain%x%global%size .OR. &
          size(global_field,2).NE.domain%y%global%size ) &
           call mpp_error( FATAL, 'MPP_GET_GLOBAL_: argument mismatch, check domain and array sizes.' )
!copy compute domain from local to global array
      global_field = 0.
      global_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index MPP_EXTRA_INDICES_ ) = &
       local_field(domain%x%compute%start_index:domain%x%compute%end_index, &
                   domain%y%compute%start_index:domain%y%compute%end_index MPP_EXTRA_INDICES_ )
      if( npes.EQ.1 )return

!      call mpp_sum( global_field(domain%x%global%start_index,domain%y%global%start_index,1), size(global_field) )
      call mpp_sum( global_field_loc, size(global_field) )

      return
    end subroutine MPP_GET_GLOBAL_
