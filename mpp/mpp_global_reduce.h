  function MPP_GLOBAL_REDUCE_3D_( domain, field, position )
    MPP_TYPE_: MPP_GLOBAL_REDUCE_
    type(domain2D), intent(in) :: domain
    integer, dimension(3), intent(out) :: position
    MPP_TYPE_, intent(in)  ::  field(domain%x%data%begin:,domain%y%data%begin:,:)

    if( size(field,1).NE.domain%x%data%size .OR. size(field,2).NE.domain%y%data%size ) &
         call mpp_error( FATAL, 'MPP_GLOBAL_REDUCE: incoming array does not match domain.' )

    call mpp_global_field( domain, field2D, global2D, flags )
    MPP_GLOBAL_REDUCE_ = reduce(global2D)

    return
  end function MPP_GLOBAL_REDUCE_
