    subroutine MPP_WRITE_2DDECOMP_2D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      real, intent(in), optional :: tstamp
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data3D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_WRITE_2DDECOMP_2D_

    subroutine MPP_WRITE_2DDECOMP_3D_( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(domain%x%data%start_index:,domain%y%data%start_index:,:)
      real, intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data
!global_domain and gdata are used to globalize data for multi-PE single-threaded I/O
!      type(domain2D), allocatable :: global_domain(:)
      MPP_TYPE_, allocatable, dimension(:,:,:) :: cdata, gdata

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( domain%x%data%is_global .AND. domain%y%data%is_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(domain%x%global%start_index:domain%x%global%end_index, &
                              domain%y%global%start_index:domain%y%global%end_index,size(data,3)) )
              call mpp_get_global( domain, data, gdata )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(domain%x%compute%size,domain%y%compute%size,size(data,3)) )
          cdata(:,:,:) = data(domain%x%compute%start_index:domain%x%compute%end_index, &
                              domain%y%compute%start_index:domain%y%compute%end_index,:)
          call write_record( unit, field, size(cdata), cdata, tstamp, domain )
      end if

      return
    end subroutine MPP_WRITE_2DDECOMP_3D_

    subroutine MPP_WRITE_2DDECOMP_4D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      real, intent(in), optional :: tstamp
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),size(data,3)*size(data,4))
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data3D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_WRITE_2DDECOMP_4D_

    subroutine MPP_WRITE_2DDECOMP_5D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      real, intent(in), optional :: tstamp
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),size(data,3)*size(data,4)*size(data,5))
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_write( unit, field, domain, data3D, tstamp )
#else
      call mpp_error( FATAL, 'MPP_WRITE_2DDECOMP_2D_: requires Cray pointers.' )
#endif
      return
    end subroutine MPP_WRITE_2DDECOMP_5D_
