  function MPP_GLOBAL_SUM_( domain, field, flags )
    MPP_TYPE_: MPP_GLOBAL_SUM_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in)  ::  field(domain%x%data%begin:,domain%y%data%begin: MPP_EXTRA_INDICES )
    integer, intent(in), optional :: flags
    MPP_TYPE_ :: field2D(domain%x%data%begin:domain%x%data%end,domain%y%data%begin:domain%y%data%end)
    MPP_TYPE_ :: global2D(domain%x%global%begin:domain%x%global%end,domain%y%global%begin:domain%y%global%end)
    integer :: i,j

    if( size(field,1).NE.domain%x%data%size .OR. size(field,2).NE.domain%y%data%size ) &
         call mpp_error( FATAL, 'MPP_GLOBAL_SUM: incoming array does not match domain.' )

    do j = domain%y%compute%begin, domain%y%compute%end
       do i = domain%x%compute%begin, domain%x%compute%end
          field2D(i,j) = sum( field(i,j MPP_EXTRA_INDICES) )
       end do
    end do

    call mpp_global_field( domain, field2D, global2D, flags )
    MPP_GLOBAL_SUM_ = sum(global2D)

    return
  end function MPP_GLOBAL_SUM_
