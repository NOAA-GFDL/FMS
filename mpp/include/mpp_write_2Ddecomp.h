    subroutine MPP_WRITE_2DDECOMP_2D_( unit, field, domain, data, tstamp)
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      real(DOUBLE_KIND), intent(in), optional :: tstamp

      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
#ifdef use_CRI_pointers
      pointer( ptr, data3D )
      ptr = LOC(data)
#else
      data3D = RESHAPE( data, SHAPE(data3D) )
#endif
      call mpp_write( unit, field, domain, data3D, tstamp )
      return
    end subroutine MPP_WRITE_2DDECOMP_2D_

    subroutine MPP_WRITE_2DDECOMP_3D_( unit, field, domain, data, tstamp )
!mpp_write writes <data> which has the domain decomposition <domain>
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(inout) :: domain !must have intent(out) as well because active domain might be reset
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      real(DOUBLE_KIND), intent(in), optional :: tstamp

!cdata is used to store compute domain as contiguous data
!gdata is used to globalize data for multi-PE single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:,:) :: cdata, gdata
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg
      integer :: ishift, jshift, position

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, y_is_global=y_is_global )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg )
      position = CENTER
      ishift = 0; jshift = 0
      if(mpp_domain_is_symmetry(domain)) then
         ishift = field%axes(1)%shift; jshift = field%axes(2)%shift
         if(ishift == 1 .AND. jshift == 1) then
            position = CORNER
         else if(ishift == 1) then
            position = EAST
         else if(jshift == 1) then
            position = NORTH
         endif
      endif

      if( size(data,1).EQ.ie-is+1+ishift .AND. size(data,2).EQ.je-js+1+jshift )then
          data_has_halos = .FALSE.
      else if( size(data,1).EQ.ied-isd+1+ishift .AND. size(data,2).EQ.jed-jsd+1+jshift )then
          data_has_halos = .TRUE.
      else
          write( stderr(),'(a,10i5)' )'MPP_WRITE_2DDECOMP fails on field '//trim(field%name)// &
               ': is,ie+ishift,js,je+jshift, isd,ied+ishift,jsd,jed+jshift, size(data,1), size(data,2)=', &
               is,ie,js,je, isd,ied,jsd,jed, size(data,1), size(data,2)
          call mpp_error( FATAL, 'MPP_WRITE: data must be either on compute domain or data domain.' )
      end if
      if(ie == ieg) ie  = ie  + ishift
      if(je == jeg) je  = je  + jshift
      ieg = ieg + ishift; jeg = jeg + jshift

      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then
              call mpp_update_domains( data, domain, position = position )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(data(:,:,:)), data, tstamp )
          else
!put field onto global domain
              allocate( gdata(isg:ieg,jsg:jeg,size(data,3)) )
              call mpp_global_field( domain, data, gdata )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if( pe.NE.0 )return
              call write_record( unit, field, size(gdata(:,:,:)), gdata, tstamp )
          end if
      else if( data_has_halos )then
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(is:ie,js:je,size(data,3)) )
          cdata(:,:,:) = data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:)
          call write_record( unit, field, size(cdata(:,:,:)), cdata, tstamp, domain )
      else
!data is already contiguous
          call write_record( unit, field, size(data(:,:,:)), data, tstamp, domain )
      end if

      return
    end subroutine MPP_WRITE_2DDECOMP_3D_
