    subroutine MPP_READ_2DDECOMP_2D_( unit, field, domain, data, tindex )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_READ_2DDECOMP_2D_

    subroutine MPP_READ_2DDECOMP_3D_( unit, field, domain, data, tindex )
!mpp_read reads <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(domain%x%data%start_index:,domain%y%data%start_index:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_, allocatable :: cdata(:,:,:)
      MPP_TYPE_, allocatable :: gdata(:)
      integer :: len, lenx,leny,lenz,i,j,k,n

      if (.NOT. present(tindex) .AND. mpp_file(unit)%time_level .ne. -1) &
      call mpp_error(FATAL, 'MPP_READ: need to specify a time level for data with time axis')

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_READ: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_READ: invalid unit number.' )

      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          lenx=(domain%x%global%end_index-domain%x%global%start_index+1)
          leny=(domain%y%global%end_index-domain%y%global%start_index+1)
          lenz=size(data,3)
          len=lenx*leny*lenz
          allocate( gdata(len))          
! read field on pe 0 and pass to all pes
          if( pe.EQ.0 ) call read_record( unit, field, len, gdata, tindex )
! broadcasting global array, this can be expensive!          
          call mpp_transmit( gdata, len, ALL_PES, gdata, len,0)
!          
          do k=1,size(data,3)
             do j=domain%y%compute%start_index,domain%y%compute%end_index
                do i=domain%x%compute%start_index,domain%x%compute%end_index
                    n=(i-domain%x%global%start_index +1)+&
                      (j-domain%y%global%start_index)*lenx+&
                      (k-1)*lenx*leny
                    data(i,j,k)=gdata(n)
                enddo
             enddo
          enddo
      else
! for uniprocessor or multithreaded read
! read compute domain as contiguous data

          allocate( cdata(domain%x%compute%size,domain%y%compute%size,size(data,3)) )
          call read_record(unit,field,size(cdata),cdata,tindex,domain)

          data(domain%x%compute%start_index:domain%x%compute%end_index,&
               domain%y%compute%start_index:domain%y%compute%end_index,:)&
               =cdata(:,:,:) 
      end if

      return
    end subroutine MPP_READ_2DDECOMP_3D_

    subroutine MPP_READ_2DDECOMP_4D_( unit, field, domain, data, tindex )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),size(data,3)*size(data,4))
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_READ_2DDECOMP_4D_

    subroutine MPP_READ_2DDECOMP_5D_( unit, field, domain, data, tindex )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      integer, intent(in), optional :: tindex
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),size(data,3)*size(data,4)*size(data,5))
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_READ_2DDECOMP_5D_
