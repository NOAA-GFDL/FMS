    subroutine MPP_WRITE_1DDECOMP_1D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(inout) :: domain
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
      type(domain1D), intent(inout) :: domain !must have intent(out) as well because active domain might be reset
      MPP_TYPE_, intent(inout) :: data(:,:)
      real, intent(in), optional :: tstamp
!cdata is used to store compute domain as contiguous data (NEW: receives either on compute or data domain)
!gdata for global single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:) :: cdata, gdata
      integer :: is, ie, isd, ied, isg, ieg
      logical :: data_has_halo, halo_is_global

      if( .NOT.mpp_io_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie  )
      call mpp_get_data_domain   ( domain, isd, ied, is_global=halo_is_global )
      call mpp_get_global_domain ( domain, isg, ieg )
      data_has_halo = size(data,1).NE.ie-is+1
      if( data_has_halo .AND. size(data,1).NE.ied-isd+1 ) &
           call mpp_error( FATAL, 'MPP_WRITE: data must be either on compute domain or data domain.' )
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halo_is_global )then
              call mpp_update_domains( data, domain )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(isg:ieg,size(data,2)) )
              call mpp_global_field( domain, data, gdata )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata), gdata, tstamp )
          end if
      else if( data_has_halo )then
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(is:ie,size(data,2)) )
          cdata(:,:) = data(is-isd+1:ie-isd+1,:)
!problem here? no domain1D version of write_record
!          call write_record( unit, field, size(cdata), cdata, tstamp, domain )
          call mpp_error( FATAL, 'MPP_WRITE: no domain1D version of write_record available.' )
      else
!data is already contiguous
          call write_record( unit, field, size(data), data, tstamp )
      end if

      return
    end subroutine MPP_WRITE_1DDECOMP_2D_

    subroutine MPP_WRITE_1DDECOMP_3D_( unit, field, domain, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain1D), intent(inout) :: domain
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
      type(domain1D), intent(inout) :: domain
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
      type(domain1D), intent(inout) :: domain
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

    
