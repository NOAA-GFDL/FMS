    subroutine MPP_WRITE_1DDECOMP_1D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:)
      real, intent(in), optional :: tstamp
#ifdef use_CRI_pointers
      MPP_TYPE_ :: data2D(size(data,1),1)
      pointer( ptr, data2D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data2D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE: requires Cray pointers.' _
#endif
      return
    end subroutine MPP_WRITE_1DDECOMP_1D_

    subroutine MPP_WRITE_1DDECOMP_2D_( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(domain%data%start_index:,:)
      real, intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data
!gdata for global single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:) :: cdata, gdata
      type(domain2D) :: write_domain(1)

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( domain%data%is_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(domain%global%start_index:domain%global%end_index,size(data,2)) )
              gdata = 0.
              gdata(domain%compute%start_index:domain%compute%end_index,:) = &
               data(domain%compute%start_index:domain%compute%end_index,:)
              call mpp_sum( gdata(domain%compute%start_index,1), size(gdata), domain%pelist )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(domain%compute%size,size(data,2)) )
          cdata(:,:) = data(domain%compute%start_index:domain%compute%end_index,:)
!write_domain is a fake 2D domain for passing to write_record
!its x axis is <domain>
!its y axis is the global undistributed second axis
          write_domain(1)%x = domain
          call mpp_define_domains( (/1,size(data,2)/), write_domain%y, flags=GLOBAL_COMPUTE_DOMAIN )
          call write_record( unit, field, size(cdata), cdata, tstamp, write_domain(1) )
      end if

      return
    end subroutine MPP_WRITE_1DDECOMP_2D_

    subroutine MPP_WRITE_1DDECOMP_3D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      real, intent(in), optional :: tstamp
#ifdef use_CRI_pointers
      MPP_TYPE_ :: data2D(size(data,1),size(data,2)*size(data,3))
      pointer( ptr, data2D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data2D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE: requires Cray pointers.' _
#endif
      return
    end subroutine MPP_WRITE_1DDECOMP_3D_

    subroutine MPP_WRITE_1DDECOMP_4D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      real, intent(in), optional :: tstamp
#ifdef use_CRI_pointers
      MPP_TYPE_ :: data2D(size(data,1),size(data,2)*size(data,3)*size(data,4))
      pointer( ptr, data2D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data2D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE: requires Cray pointers.' _
#endif
      return
    end subroutine MPP_WRITE_1DDECOMP_4D_

    subroutine MPP_WRITE_1DDECOMP_5D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      real, intent(in), optional :: tstamp
#ifdef use_CRI_pointers
      MPP_TYPE_ :: data2D(size(data,1),size(data,2)*size(data,3)*size(data,4)*size(data,5))
      pointer( ptr, data2D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data2D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE: requires Cray pointers.' _
#endif
      return
    end subroutine MPP_WRITE_1DDECOMP_5D_

    
