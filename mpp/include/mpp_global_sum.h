  function MPP_GLOBAL_SUM_( domain, field, flags, new )
    MPP_TYPE_ :: MPP_GLOBAL_SUM_
    type(domain2D), intent(in) :: domain
    MPP_TYPE_, intent(in) :: field(:,: MPP_EXTRA_INDICES_ )
    integer, intent(in), optional :: flags
    logical, intent(in), optional :: new
    MPP_TYPE_, allocatable, dimension(:),save :: field1D
    MPP_TYPE_, dimension(domain%x%compute%begin:domain%x%compute%end,domain%y%compute%begin:domain%y%compute%end) :: field2D
    pointer(ptr_field2D,field2D)
    MPP_TYPE_, allocatable, dimension(:,:) :: field2Dold, global2D
    integer :: i,j, ioff,joff
    logical :: use_new
    integer, save :: f1D_len=0
    integer(LONG_KIND), save :: f_addr
    type(domain2D), save :: domain_prev


    if( size(field,1).EQ.domain%x%compute%size .AND. size(field,2).EQ.domain%y%compute%size )then
!field is on compute domain
        ioff = -domain%x%compute%begin + 1
        joff = -domain%y%compute%begin + 1
    else if( size(field,1).EQ.domain%x%data%size .AND. size(field,2).EQ.domain%y%data%size )then
!field is on data domain
        ioff = -domain%x%data%begin + 1
        joff = -domain%y%data%begin + 1
    else
        call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: incoming field array must match either compute domain or data domain.' )
    end if
    if( PRESENT(flags) )then
        if( flags.NE.BITWISE_EXACT_SUM )call mpp_error( FATAL, 'MPP_GLOBAL_SUM_: only valid flag is BITWISE_EXACT_SUM.' )
!this is bitwise exact across different PE counts.

        use_new=.false.; if(PRESENT(new))use_new=new
        if(use_new)then
          if(f1D_len<(domain%x%compute%end-domain%x%compute%begin+1)*(domain%y%compute%end-domain%y%compute%begin+1))then
            if(ALLOCATED(field1D))then
              call mpp_global_field_free_comm(domain_prev,f_addr,ksize=1)
              deallocate(field1D)
            endif
            f1D_len = (domain%x%compute%end-domain%x%compute%begin+1)*(domain%y%compute%end-domain%y%compute%begin+1)
            allocate(field1D(f1D_len))
            f_addr=LOC(field1D); domain_prev%id = domain%id
          endif
          ptr_field2D = f_addr
        else
          allocate( field2Dold(domain%x%compute%begin:domain%x%compute%end,domain%y%compute%begin:domain%y%compute%end) )
          ptr_field2D = LOC(field2Dold)
        endif

        allocate( global2D(domain%x%global%size,domain%y%global%size) )
        do j = domain%y%compute%begin, domain%y%compute%end
           do i = domain%x%compute%begin, domain%x%compute%end
              field2D(i,j) = sum( field(i+ioff:i+ioff,j+joff:j+joff MPP_EXTRA_INDICES_) )
           end do
        end do

        call mpp_global_field( domain, field2D, global2D, new=new )
        MPP_GLOBAL_SUM_ = sum(global2D)
        deallocate(global2D); if(allocated(field2Dold))deallocate(field2Dold)
    else
!this is not bitwise-exact across different PE counts
        MPP_GLOBAL_SUM_ = sum( field(domain%x%compute%begin+ioff:domain%x%compute%end+ioff, &
                                     domain%y%compute%begin+joff:domain%y%compute%end+joff MPP_EXTRA_INDICES_) )
        call mpp_sum( MPP_GLOBAL_SUM_, domain%list(:)%pe )
    end if

    return
  end function MPP_GLOBAL_SUM_
