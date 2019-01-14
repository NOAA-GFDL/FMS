!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
    subroutine WRITE_RECORD_( unit, field, nwords, data, time_in, domain, tile_count)
!routine that is finally called by all mpp_write routines to perform the write
!a non-netCDF record contains:
!      field ID
!      a set of 4 coordinates (is:ie,js:je) giving the data subdomain
!      a timelevel and a timestamp (=NULLTIME if field is static)
!      3D real data (stored as 1D)
!if you are using direct access I/O, the RECL argument to OPEN must be large enough for the above
!in a global direct access file, record position on PE is given by %record.

!Treatment of timestamp:
!   We assume that static fields have been passed without a timestamp.
!   Here that is converted into a timestamp of NULLTIME.
!   For non-netCDF fields, field is treated no differently, but is written
!   with a timestamp of NULLTIME. There is no check in the code to prevent
!   the user from repeatedly writing a static field.

      integer,           intent(in)           :: unit, nwords
      type(fieldtype),   intent(in)           :: field
      MPP_TYPE_,         intent(in)           :: data(nwords)
      MPP_TYPE_,         intent(in), optional :: time_in
      type(domain2D),    intent(in), optional :: domain
      integer,           intent(in), optional :: tile_count
      integer, dimension(size(field%axes(:))) :: start, axsiz
      real(DOUBLE_KIND) :: time
      integer :: time_level
      logical :: newtime
      integer :: subdomain(4)
      integer :: packed_data(nwords)
      integer :: i, is, ie, js, je

      real(FLOAT_KIND) :: data_r4(nwords)
      pointer( ptr1, data_r4)
      pointer( ptr2, packed_data)

      if (mpp_io_stack_size < nwords) call mpp_io_set_stack_size(nwords)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%write_on_this_pe) return
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )
      if( .NOT.mpp_file(unit)%initialized )then
!this is the first call to mpp_write
!we now declare the file to be initialized
!if this is netCDF we switch file from DEFINE mode to DATA mode
          if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
!NOFILL is probably required for parallel: any circumstances in which not advisable?
              error = NF_SET_FILL( mpp_file(unit)%ncid, NF_NOFILL, i ); call netcdf_err( error, mpp_file(unit) )
              if( mpp_file(unit)%action.EQ.MPP_WRONLY )then
                 if(header_buffer_val>0) then
                    error = NF__ENDDEF(mpp_file(unit)%ncid,header_buffer_val,4,0,4)
                 else
                    error = NF_ENDDEF(mpp_file(unit)%ncid)
                 endif
              endif
              call netcdf_err( error, mpp_file(unit) )
#endif
          else
              call mpp_write_meta( unit, 'END', cval='metadata' )
          end if
          mpp_file(unit)%initialized = .TRUE.
          if( verbose )print '(a,i6,a)', 'MPP_WRITE: PE=', pe, ' initialized file '//trim(mpp_file(unit)%name)//'.'
      end if

!initialize time: by default assume NULLTIME
      time = NULLTIME
      time_level = -1
      newtime = .FALSE.
      if( PRESENT(time_in) )time = time_in
!increment time level if new time
      if( time.GT.mpp_file(unit)%time+EPSILON(time) )then !new time
          mpp_file(unit)%time_level = mpp_file(unit)%time_level + 1
          mpp_file(unit)%time = time
          newtime = .TRUE.
      end if
      if( verbose )print '(a,2i6,2i5,es13.5)', 'MPP_WRITE: PE, unit, %id, %time_level, %time=',&
           pe, unit, mpp_file(unit)%id, mpp_file(unit)%time_level, mpp_file(unit)%time

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
          ptr2 = LOC(mpp_io_stack(1))
!define netCDF data block to be written:
!  time axis: START = time level
!             AXSIZ = 1
!  space axis: if there is no domain info
!              START = 1
!              AXSIZ = field%size(axis)
!          if there IS domain info:
!              start of domain is compute%start_index for multi-file I/O
!                                 global%start_index for all other cases
!              this number must be converted to 1 for NF_PUT_VAR
!                  (netCDF fortran calls are with reference to 1),
!          So, START = compute%start_index - <start of domain> + 1
!              AXSIZ = usually compute%size
!          However, if compute%start_index-compute%end_index+1.NE.compute%size,
!              we assume that the call is passing a subdomain.
!              To pass a subdomain, you must pass a domain2D object that satisfies the following:
!                  global%start_index must contain the <start of domain> as defined above;
!                  the data domain and compute domain must refer to the subdomain being passed.
!              In this case, START = compute%start_index - <start of domain> + 1
!                            AXSIZ = compute%start_index - compute%end_index + 1
! NOTE: passing of subdomains will fail for multi-PE single-threaded I/O,
!       since that attempts to gather all data on PE 0.
          start = 1
          do i = 1,size(field%axes(:))
             axsiz(i) = field%size(i)
             if( i.EQ.field%time_axis_index )start(i) = mpp_file(unit)%time_level
             start(i) = max(start(i),1)
          end do

          if( debug )print '(a,2i6,12i6)', 'WRITE_RECORD: PE, unit, start, axsiz=', pe, unit, start, axsiz
#ifdef use_netCDF
!write time information if new time
          if( newtime )then
              if( KIND(time).EQ.DOUBLE_KIND )then
                  error = NF_PUT_VAR1_DOUBLE( mpp_file(unit)%ncid, mpp_file(unit)%id, mpp_file(unit)%time_level, time )
              else if( KIND(time).EQ.FLOAT_KIND )then
                  error = NF_PUT_VAR1_REAL  ( mpp_file(unit)%ncid, mpp_file(unit)%id, mpp_file(unit)%time_level, time )
              end if
          end if
          if( field%pack == 0 )then
              packed_data = CEILING(data)
              error = NF_PUT_VARA_INT   ( mpp_file(unit)%ncid, field%id, start, axsiz, packed_data )
          elseif( field%pack.GT.0 .and. field%pack.LE.2 )then
              if( KIND(data).EQ.DOUBLE_KIND )then
                  error = NF_PUT_VARA_DOUBLE( mpp_file(unit)%ncid, field%id, start, axsiz, data )
              else if( KIND(data).EQ.FLOAT_KIND )then
                  error = NF_PUT_VARA_REAL  ( mpp_file(unit)%ncid, field%id, start, axsiz, data )
              end if
          else              !convert to integer using scale and add: no error check on packed data representation
              packed_data = nint((data-field%add)/field%scale)
              error = NF_PUT_VARA_INT   ( mpp_file(unit)%ncid, field%id, start, axsiz, packed_data )
          end if
          call netcdf_err( error, mpp_file(unit), field=field )
#endif
      else                      !non-netCDF
          ptr1 = LOC(mpp_io_stack(1))
!subdomain contains (/is,ie,js,je/)
          if( PRESENT(domain) )then
              call mpp_get_compute_domain(domain, is, ie, js, je)
              subdomain(:) = (/ is, ie, js, je /)
          else
              subdomain(:) = -1    ! -1 means use global value from axis metadata
          end if
          if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
              write( unit,* )field%id, subdomain, time_level, time, data
          else                      !MPP_IEEE32 or MPP_NATIVE
              if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
#ifdef __sgi
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write(unit)field%id, subdomain, time_level, time, data_r4
                  else
                      write(unit)field%id, subdomain, time_level, time, data
                  end if
#else
                  write(unit)field%id, subdomain, time_level, time, data
#endif
              else                  !MPP_DIRECT
#ifdef __sgi
                  if( mpp_file(unit)%format.EQ.MPP_IEEE32 )then
                      data_r4 = data !IEEE conversion layer on SGI until assign -N ieee_32 is supported
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data_r4
                  else
                      write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
                  end if
#else
                  write( unit, rec=mpp_file(unit)%record )field%id, subdomain, time_level, time, data
#endif
                  if( debug )print '(a,i6,a,i6)', 'MPP_WRITE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              end if
          end if
      end if

!recompute current record for direct access I/O
      if( mpp_file(unit)%access.EQ.MPP_DIRECT )then
          if( mpp_file(unit)%fileset.EQ.MPP_SINGLE )then
!assumes all PEs participate in I/O: modify later
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe*npes
          else
              mpp_file(unit)%record = mpp_file(unit)%record + records_per_pe
          end if
      end if

      return
    end subroutine WRITE_RECORD_

    subroutine MPP_WRITE_2DDECOMP_2D_( unit, field, domain, data, tstamp, tile_count, default_data)
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      MPP_TYPE_,              intent(in), optional :: tstamp
      integer,           intent(in), optional :: tile_count
      MPP_TYPE_,         intent(in), optional :: default_data

      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
      pointer( ptr, data3D )
      ptr = LOC(data)

      call mpp_write( unit, field, domain, data3D, tstamp, tile_count, default_data)
      return
    end subroutine MPP_WRITE_2DDECOMP_2D_

    subroutine MPP_WRITE_2DDECOMP_3D_( unit, field, domain, data, tstamp, tile_count, default_data)
!mpp_write writes <data> which has the domain decomposition <domain>
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(in)           :: field
      type(domain2D),    intent(inout)        :: domain 
      MPP_TYPE_,         intent(inout)        :: data(:,:,:)
      MPP_TYPE_,         intent(in), optional :: tstamp
      integer,           intent(in), optional :: tile_count
      MPP_TYPE_,         intent(in), optional :: default_data

!cdata is used to store compute domain as contiguous data
!gdata is used to globalize data for multi-PE single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:,:) :: cdata, gdata
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ism, iem, jsm, jem
      integer :: position, errunit
      type(domain2d), pointer :: io_domain=>NULL()

      call mpp_clock_begin(mpp_write_clock)

      errunit = stderr()
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      position = field%position

      call mpp_get_compute_domain( domain, is,  ie,  js,  je, tile_count=tile_count, position=position  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, &
                                   y_is_global=y_is_global, tile_count=tile_count, position=position )
      call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )

      if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
          data_has_halos = .FALSE.
      else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then
          data_has_halos = .TRUE.
      else
          write( errunit,'(a,10i5)' )'MPP_WRITE_2DDECOMP fails on field '//trim(field%name)// &
               ': is,ie,js,je, ism,iem,jsm,jem, size(data,1), size(data,2)=', &
               is,ie,js,je, ism,iem,jsm,jem, size(data,1), size(data,2)
          call mpp_error( FATAL, 'MPP_WRITE: data must be either on compute domain or data domain.' )
      end if
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then
              call mpp_update_domains( data, domain, position = position )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(data(:,:,:)), data, tstamp)
              endif
          else
!put field onto global domain
              call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
              if(mpp_file(unit)%write_on_this_pe .OR. .NOT. global_field_on_root_pe) then
                  allocate( gdata(isg:ieg,jsg:jeg,size(data,3)) )
              else
                  allocate( gdata(1,1,1))
              endif
              if(global_field_on_root_pe) then
                 call mpp_global_field( domain, data, gdata, position = position, &
                                        flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY,   &
                                        default_data=default_data)
              else
                 call mpp_global_field( domain, data, gdata, position = position, &
                                        default_data=default_data)
              endif
!all non-0 PEs have passed their data to PE 0 and may now exit
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(gdata(:,:,:)), gdata, tstamp)
              endif
              deallocate(gdata)
          end if
      else if(mpp_file(unit)%io_domain_exist ) then
          if( halos_are_global )then
              call mpp_update_domains( data, domain, position = position )
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(data(:,:,:)), data, tstamp)
              endif
          else
              io_domain=>mpp_get_io_domain(mpp_file(unit)%domain) 
              call mpp_get_global_domain ( io_domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
              if(mpp_file(unit)%write_on_this_pe .OR. .NOT. global_field_on_root_pe) then
                 allocate( gdata(isg:ieg,jsg:jeg,size(data,3)) )
              else
                 allocate( gdata(1,1,1))
              endif
              if(global_field_on_root_pe) then
                 call mpp_global_field( io_domain, data, gdata, position = position, &
                                        flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY,      &
                                        default_data=default_data)
              else
                 call mpp_global_field( io_domain, data, gdata, position = position, &
                                        default_data=default_data)
              endif
              io_domain => NULL()
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(gdata(:,:,:)), gdata, tstamp)
              endif
              deallocate( gdata )
          endif
      else if( data_has_halos )then
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(is:ie,js:je,size(data,3)) )
          cdata(:,:,:) = data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:)
          call WRITE_RECORD_( unit, field, size(cdata(:,:,:)), cdata, tstamp, domain, tile_count ) 
      else
!data is already contiguous
          call WRITE_RECORD_( unit, field, size(data(:,:,:)), data, tstamp, domain, tile_count )
      end if

      call mpp_clock_end(mpp_write_clock)

      return
    end subroutine MPP_WRITE_2DDECOMP_3D_

    subroutine MPP_WRITE_2DDECOMP_4D_( unit, field, domain, data, tstamp, tile_count, default_data)
!mpp_write writes <data> which has the domain decomposition <domain>
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(in)           :: field
      type(domain2D),    intent(inout)        :: domain 
      MPP_TYPE_,         intent(inout)        :: data(:,:,:,:)
      MPP_TYPE_,         intent(in), optional :: tstamp
      integer,           intent(in), optional :: tile_count
      MPP_TYPE_,         intent(in), optional :: default_data

!cdata is used to store compute domain as contiguous data
!gdata is used to globalize data for multi-PE single-threaded I/O
      MPP_TYPE_, allocatable, dimension(:,:,:,:) :: cdata, gdata
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ism, iem, jsm, jem
      integer :: position, errunit
      type(domain2d), pointer :: io_domain=>NULL()

      errunit = stderr()
      call mpp_clock_begin(mpp_write_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_WRITE: invalid unit number.' )

      position = field%position

      call mpp_get_compute_domain( domain, is,  ie,  js,  je, tile_count=tile_count, position=position  )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, &
                                   y_is_global=y_is_global, tile_count=tile_count, position=position )
      call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=position )

      if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
          data_has_halos = .FALSE.
      else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then
          data_has_halos = .TRUE.
      else
          write( errunit,'(a,10i5)' )'MPP_WRITE_2DDECOMP fails on field '//trim(field%name)// &
               ': is,ie,js,je, ism,iem,jsm,jem, size(data,1), size(data,2)=', &
               is,ie,js,je, ism,iem,jsm,jem, size(data,1), size(data,2)
          call mpp_error( FATAL, 'MPP_WRITE: data must be either on compute domain or data domain.' )
      end if
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then
              call mpp_update_domains( data, domain, position = position )
!all non-0 PEs have passed their data to PE 0 and may now exit
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(data(:,:,:,:)), data, tstamp)
              endif
          else
!put field onto global domain
              call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
              if(mpp_file(unit)%write_on_this_pe .OR. .NOT. global_field_on_root_pe) then
                  allocate( gdata(isg:ieg,jsg:jeg,size(data,3),size(data,4)) )
              else
                  allocate( gdata(1,1,1,1))
              endif
              if(global_field_on_root_pe) then
                 call mpp_global_field( domain, data, gdata, position = position, &
                                        flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY, &
                                        default_data=default_data)
              else
                 call mpp_global_field( domain, data, gdata, position = position, &
                                        default_data=default_data)
              endif
!all non-0 PEs have passed their data to PE 0 and may now exit
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(gdata(:,:,:,:)), gdata, tstamp)
              endif
              deallocate(gdata)
          end if
      else if(mpp_file(unit)%io_domain_exist ) then
          if( halos_are_global )then
              if(npes .GT. 1) call mpp_update_domains( data, domain, position = position )
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(data(:,:,:,:)), data, tstamp)
              endif
          else
              io_domain=>mpp_get_io_domain(mpp_file(unit)%domain) 
              call mpp_get_global_domain ( io_domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
              if(mpp_file(unit)%write_on_this_pe .OR. .NOT. global_field_on_root_pe) then
                 allocate( gdata(isg:ieg,jsg:jeg,size(data,3),size(data,4)) )
              else
                 allocate( gdata(1,1,1,1))
              endif
              if(global_field_on_root_pe) then
                 call mpp_global_field( io_domain, data, gdata, position = position, &
                                        flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY,      &
                                        default_data=default_data)
              else
                 call mpp_global_field( io_domain, data, gdata, position = position, &
                                        default_data=default_data)
              endif
              io_domain => NULL()
              if(mpp_file(unit)%write_on_this_pe ) then
                 call WRITE_RECORD_( unit, field, size(gdata(:,:,:,:)), gdata, tstamp)
              endif
              deallocate( gdata )
          endif
      else if( data_has_halos )then
!store compute domain as contiguous data and pass to write_record
          allocate( cdata(is:ie,js:je,size(data,3),size(data,4)) )
          cdata(:,:,:,:) = data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:,:)
          call WRITE_RECORD_( unit, field, size(cdata(:,:,:,:)), cdata, tstamp, domain, tile_count ) 
      else
!data is already contiguous
          call WRITE_RECORD_( unit, field, size(data(:,:,:,:)), data, tstamp, domain, tile_count )
      end if

      call mpp_clock_end(mpp_write_clock)

      return
    end subroutine MPP_WRITE_2DDECOMP_4D_

