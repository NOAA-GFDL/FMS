  function MPP_GLOBAL_SUM_( domain, field, flags, new )
    MPP_TYPE_ :: MPP_GLOBAL_SUM_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,: MPP_EXTRA_INDICES_ )
    integer, intent(in), optional :: flags
    logical, intent(in), optional :: new

!    MPP_TYPE_, allocatable, dimension(:),save :: field1D
    MPP_TYPE_, dimension(:,:), allocatable :: field2D
!z1l    MPP_TYPE_, dimension(domain%x%compute%begin:domain%x%compute%end+1,domain%y%compute%begin:domain%y%compute%end+1) :: field2D
!    pointer(ptr_field2D,field2D)
    MPP_TYPE_, allocatable, dimension(:,:) :: field2Dold, global2D
    integer :: i,j, ioff,joff, ishift, jshift, isc, iec, jsc, jec, is, js
    integer :: gxsize, gysize       
!    logical :: use_new
!    integer, save :: f1D_len=0
!    integer(LONG_KIND), save :: f_addr
!    type(domain2D), save :: domain_prev
    integer :: global_flag

    global_flag = NON_BITWISE_EXACT_SUM
    if(present(flags)) global_flag = flags

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_global_domain(domain, xsize = gxsize, ysize = gysize)
    call mpp_get_global_shift(domain, size(field,1), size(field,2), ishift, jshift)

    if( size(field,1).EQ.domain%x%compute%size+ishift .AND. size(field,2).EQ.domain%y%compute%size+jshift )then
!field is on compute domain
        ioff = -domain%x%compute%begin + 1 
        joff = -domain%y%compute%begin + 1
    else if( size(field,1).EQ.domain%x%data%size+ishift .AND. size(field,2).EQ.domain%y%data%size+jshift )then
!field is on data domain
        ioff = -domain%x%data%begin + 1
        joff = -domain%y%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: incoming field array must match either compute domain or data domain.' )
    end if

    !--- Add extra point and/or shift the index of data if needed.
    is = isc
    if(ishift .NE. 0) then
       iec = iec + ishift
       if( domain%x%cyclic ) then
          isc = isc + ishift
          is  = isc
       else if(domain%x%compute%begin == domain%x%global%begin) then
          gxsize = gxsize+ishift
          is = isc;
       else
          gxsize = gxsize+ishift
          is = isc + ishift
       endif
    endif
    js = jsc
    if(jshift .NE. 0) then
       jec = jec + jshift
       if( domain%y%cyclic ) then
          jsc = jsc + jshift
          js  = jsc
       else if(domain%y%compute%begin == domain%y%global%begin) then
          gysize = gysize+jshift
          js = jsc
       else
          gysize = gysize+jshift
          js = jsc + jshift
       endif
    endif

    if( global_flag == BITWISE_EXACT_SUM )then
        !this is bitwise exact across different PE counts.

        !--- z1l: comments out the following because cray pointer will cause problem when
        !---      implementing symmetric domain.
!!$        use_new=.false.; if(PRESENT(new))use_new=new
!!$        if(use_new)then
!!$          
!!$
!!$          if(f1D_len<(iec-isc+1)*(jec-jsc+1) )then
!!$            if(ALLOCATED(field1D))then
!!$              call mpp_global_field_free_comm(domain_prev,f_addr,ksize=1)
!!$              deallocate(field1D)
!!$            endif
!!$            f1D_len = (iec-isc+1)*(jec-jsc+1)
!!$            allocate(field1D(f1D_len))
!!$            f_addr=LOC(field1D); domain_prev%id = domain%id
!!$          endif
!!$          ptr_field2D = f_addr
!!$        else
!!$          allocate( field2Dold(isc:iec+1, jsc:jec+1) )
!!$          ptr_field2D = LOC(field2Dold)
!!$        endif
        allocate( field2D (isc:iec,jsc:jec) )
        allocate( global2D( gxsize, gysize ) )
        global2D = 0.
        do j = jsc, jec
           do i = isc, iec
               field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff MPP_EXTRA_INDICES_) )
           end do
        end do

        call mpp_global_field( domain, field2D(isc:iec,jsc:jec), global2D, new=new )
        MPP_GLOBAL_SUM_ = sum(global2D)
        deallocate(global2D); if(allocated(field2Dold))deallocate(field2Dold)
    else
        !this is not bitwise-exact across different PE counts
        MPP_GLOBAL_SUM_ = sum( field(is+ioff:iec+ioff, js+joff:jec+joff MPP_EXTRA_INDICES_) )
        call mpp_sum( MPP_GLOBAL_SUM_, domain%list(:)%pe )
    end if

    return
  end function MPP_GLOBAL_SUM_
