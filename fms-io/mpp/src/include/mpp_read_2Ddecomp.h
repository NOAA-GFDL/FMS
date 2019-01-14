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
    subroutine READ_RECORD_CORE_(unit, field, nwords, data, start, axsiz)
      integer,         intent(in)    :: unit
      type(fieldtype), intent(in)    :: field
      integer,         intent(in)    :: nwords
      MPP_TYPE_,        intent(inout) :: data(nwords)
      integer,         intent(in)    :: start(:), axsiz(:)

      integer(SHORT_KIND) :: i2vals(nwords)
!rab used in conjunction with transfer intrinsic to determine size of a variable
      integer(KIND=1) :: one_byte(8)
      integer         :: word_sz
!#ifdef __sgi
      integer(INT_KIND) :: ivals(nwords)
      real(FLOAT_KIND) :: rvals(nwords)
!#else
!      integer :: ivals(nwords)
!      real :: rvals(nwords)
!#endif

      real(DOUBLE_KIND) :: r8vals(nwords)
      pointer( ptr1, i2vals )
      pointer( ptr2, ivals )
      pointer( ptr3, rvals )
      pointer( ptr4, r8vals )
      if (mpp_io_stack_size < nwords) call mpp_io_set_stack_size(nwords)

#ifdef use_netCDF
      word_sz = size(transfer(data(1),one_byte))

          select case (field%type)
             case(NF_BYTE)
! use type conversion
                call mpp_error( FATAL, 'MPP_READ: does not support NF_BYTE packing' )
             case(NF_SHORT)
                ptr1 = LOC(mpp_io_stack(1))
                error = NF_GET_VARA_INT2  ( mpp_file(unit)%ncid, field%id, start, axsiz, i2vals )
                call netcdf_err( error, mpp_file(unit), field=field )
                if(field%scale == 1.0 .and. field%add == 0.0) then
                   data(:)=i2vals(:)
                else
                   data(:)=i2vals(:)*field%scale + field%add
                end if
             case(NF_INT)

                ptr2 = LOC(mpp_io_stack(1))

                error = NF_GET_VARA_INT   ( mpp_file(unit)%ncid, field%id, start, axsiz, ivals  )
                call netcdf_err( error, mpp_file(unit), field=field )
                if(field%scale == 1.0 .and. field%add == 0.0) then
                   data(:)=ivals(:)
                else
                   data(:)=ivals(:)*field%scale + field%add
                end if
             case(NF_FLOAT)
                ptr3 = LOC(mpp_io_stack(1))
                if (size(transfer(rvals(1),one_byte)) .eq. word_sz) then
                  error = NF_GET_VARA_REAL  ( mpp_file(unit)%ncid, field%id, start, axsiz, data  )
                  call netcdf_err( error, mpp_file(unit), field=field )
                  if(field%scale /= 1.0 .or. field%add /= 0.0) then
                     data(:)=data(:)*field%scale + field%add
                  end if
                else
                  error = NF_GET_VARA_REAL  ( mpp_file(unit)%ncid, field%id, start, axsiz, rvals  )
                  call netcdf_err( error, mpp_file(unit), field=field )
                  if(field%scale == 1.0 .and. field%add == 0.0) then
                     data(:)=rvals(:)
                  else
                     data(:)=rvals(:)*field%scale + field%add
                  end if
                end if
             case(NF_DOUBLE)
                ptr4 = LOC(mpp_io_stack(1))
                if (size(transfer(r8vals(1),one_byte)) .eq. word_sz) then
                  error = NF_GET_VARA_DOUBLE( mpp_file(unit)%ncid, field%id, start, axsiz, data )
                  call netcdf_err( error, mpp_file(unit), field=field )
                  if(field%scale /= 1.0 .or. field%add /= 0.0) then
                     data(:)=data(:)*field%scale + field%add
                  end if
                else
                  error = NF_GET_VARA_DOUBLE( mpp_file(unit)%ncid, field%id, start, axsiz, r8vals )
                  call netcdf_err( error, mpp_file(unit), field=field )
                  if(field%scale == 1.0 .and. field%add == 0.0) then
                     data(:)=r8vals(:)
                  else
                     data(:)=r8vals(:)*field%scale + field%add
                  end if
                end if
             case default
                call mpp_error( FATAL, 'MPP_READ: invalid pack value' )
          end select
#else
      call mpp_error( FATAL, 'MPP_READ currently requires use_netCDF option' )
#endif

    end subroutine READ_RECORD_CORE_


    subroutine READ_RECORD_( unit, field, nwords, data, time_level, domain, position, tile_count, start_in, axsiz_in )
!routine that is finally called by all mpp_read routines to perform the read
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

      integer,         intent(in)             :: unit, nwords
      type(fieldtype), intent(in)             :: field
      MPP_TYPE_,      intent(inout)           :: data(nwords)
      integer,        intent(in),    optional :: time_level
      type(domain2D), intent(in),    optional :: domain
      integer,        intent(in),    optional :: position, tile_count
      integer,        intent(in),    optional :: start_in(:), axsiz_in(:)
      integer, dimension(size(field%axes(:))) :: start, axsiz

      integer :: tlevel !,subdomain(4)


      integer :: i, error, is, ie, js, je, isg, ieg, jsg, jeg
      type(domain2d), pointer :: io_domain=>NULL()

      if (.not.PRESENT(time_level)) then
          tlevel = 0
      else
          tlevel = time_level
      endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'READ_RECORD: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'READ_RECORD: invalid unit number.' )
      if( .NOT.mpp_file(unit)%read_on_this_pe )return

      if( .NOT.mpp_file(unit)%initialized ) call mpp_error( FATAL, 'MPP_READ: must first call mpp_read_meta.' )
      if( mpp_file(unit)%format .NE. MPP_NETCDF ) call mpp_error( FATAL, 'Currently dont support non-NetCDF mpp read' )

      if (.not.PRESENT(time_level)) then
          tlevel = 0
      else
          tlevel = time_level
      endif

      if( verbose )print '(a,2i6,2i5)', 'MPP_READ: PE, unit, %id, %time_level =',&
           pe, unit, mpp_file(unit)%id, tlevel
      if( PRESENT(start_in) .AND. PRESENT(axsiz_in) ) then
         if(size(start(:)) > size(start_in(:)) )call mpp_error( FATAL, 'MPP_READ: size(start_in) < size(start)')
         if(size(axsiz(:)) > size(axsiz_in(:)) )call mpp_error( FATAL, 'MPP_READ: size(axsiz_in) < size(axsiz)')
         start(:) = start_in(1:size(start(:)))
         axsiz(:) = axsiz_in(1:size(axsiz(:)))
      else
!define netCDF data block to be read:
!  time axis: START = time level
!             AXSIZ = 1
!  space axis: if there is no domain info
!              START = 1
!              AXSIZ = field%size(axis)
!          if there IS domain info:
!              start of domain is compute%start_index for multi-file I/O
!                                 global%start_index for all other cases
!              this number must be converted to 1 for NF_GET_VAR
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
             if( i .EQ. field%time_axis_index )start(i) = tlevel
          end do
          if( PRESENT(domain) )then
              call mpp_get_compute_domain( domain, is,  ie,  js,  je, tile_count=tile_count, position=position  )
              call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
              axsiz(1) = ie-is+1
              axsiz(2) = je-js+1
              if( mpp_file(unit)%fileset.EQ.MPP_SINGLE )then
                 if( npes.GT.1 )then
                    start(1) = is - isg + 1
                    start(2) = js - jsg + 1
                 else   !--- z1l fix a problem related obc when npes = 1
                    if( ie-is+1.NE.ieg-isg+1 )then
                       start(1) = is - isg + 1
                       axsiz(1) = ie - is + 1
                    end if
                    if( je-js+1.NE.jeg-jsg+1 )then
                       start(2) = js - jsg + 1
                       axsiz(2) = je - js + 1
                    end if
                 end if
              else if( mpp_file(unit)%io_domain_exist ) then
                 io_domain=>mpp_get_io_domain(domain)
                 call mpp_get_compute_domain( io_domain, is,  ie,  js,  je, tile_count=tile_count, position=position  )
                 call mpp_get_global_domain ( io_domain, isg, ieg, jsg, jeg, tile_count=tile_count, position=position )
                 start(1) = is - isg + 1
                 start(2) = js - jsg + 1
                 io_domain => NULL()
              end if
          end if
      endif
      if( verbose )print '(a,2i6,i6,12i4)', 'READ_RECORD: PE, unit, nwords, start, axsiz=', pe, unit, nwords, start, axsiz

      call READ_RECORD_CORE_(unit, field, nwords, data, start, axsiz)

      return
    end subroutine READ_RECORD_

    subroutine MPP_READ_2DDECOMP_2D_( unit, field, domain, data, tindex, tile_count )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in), optional :: tindex, tile_count
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),1)
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex, tile_count)
      return
    end subroutine MPP_READ_2DDECOMP_2D_

    subroutine MPP_READ_2DDECOMP_3D_( unit, field, domain, data, tindex, tile_count )
!mpp_read reads <data> which has the domain decomposition <domain>
      integer,           intent(in) :: unit
      type(fieldtype),   intent(in) :: field
      type(domain2D),    intent(in) :: domain
      MPP_TYPE_,      intent(inout) :: data(:,:,:)
      integer, intent(in), optional :: tindex, tile_count

      MPP_TYPE_, allocatable :: cdata(:,:,:)
      MPP_TYPE_, allocatable :: gdata(:)
      integer :: len, lenx,leny,lenz,i,j,k,n
!NEW: data may be on compute OR data domain
      logical :: data_has_halos, halos_are_global, x_is_global, y_is_global
      integer :: is, ie, js, je, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ism, iem, jsm, jem
      integer :: ioff, joff, position

      call mpp_clock_begin(mpp_read_clock)
      
      if (.NOT. present(tindex) .AND. mpp_file(unit)%time_level .ne. -1) &
      call mpp_error(FATAL, 'MPP_READ: need to specify a time level for data with time axis')

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_READ: invalid unit number.' )

      call mpp_get_compute_domain( domain, is,  ie,  js,  je, tile_count=tile_count )
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, x_is_global=x_is_global, &
                                   y_is_global=y_is_global, tile_count=tile_count )
      call mpp_get_memory_domain ( domain, ism, iem, jsm, jem )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, tile_count=tile_count )

      ! when domain is symmetry, extra point is needed for some data on x/y direction
      position = CENTER
      if(mpp_domain_is_symmetry(domain)) then
         if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 ) then  ! CENTER
            data_has_halos = .FALSE.
         else if( size(data,1).EQ.ie-is+2 .AND. size(data,2).EQ.je-js+1 ) then ! EAST
            data_has_halos = .FALSE.
            position = EAST
            ie = ie + 1
         else if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+2 ) then ! NORTH
            position = NORTH
            data_has_halos = .FALSE.
            je = je + 1
         else if( size(data,1).EQ.ie-is+2 .AND. size(data,2).EQ.je-js+2 ) then ! CORNER
            position = CORNER
            data_has_halos = .FALSE.
            ie = ie + 1;   je = je + 1
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then ! CENTER
            data_has_halos = .TRUE.
         else if( size(data,1).EQ.iem-ism+2 .AND. size(data,2).EQ.jem-jsm+1 )then ! EAST
            position = EAST
            data_has_halos = .TRUE.
            ie = ie + 1
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+2 )then ! NORTH
            position = NORTH
            data_has_halos = .TRUE.
            je = je + 1
         else if( size(data,1).EQ.iem-ism+2 .AND. size(data,2).EQ.jem-jsm+2 )then ! CORNER
            position = CORNER
            data_has_halos = .TRUE.
            ie = ie + 1;  je = je + 1
         else
            call mpp_error( FATAL, 'MPP_READ: when domain is symmetry, data must be either on ' &
                      //'compute domain or data domain with the consideration of shifting.' )
         end if
      else
         if( size(data,1).EQ.ie-is+1 .AND. size(data,2).EQ.je-js+1 )then
            data_has_halos = .FALSE.
         else if( size(data,1).EQ.iem-ism+1 .AND. size(data,2).EQ.jem-jsm+1 )then
            data_has_halos = .TRUE.
         else
            call mpp_error( FATAL, 'MPP_READ: data must be either on compute domain or data domain.' )
         end if
      endif
      halos_are_global = x_is_global .AND. y_is_global
      if( npes.GT.1 .AND. mpp_file(unit)%threading.EQ.MPP_SINGLE )then
          if( halos_are_global )then !you can read directly into data array
              if( pe.EQ.0 )call READ_RECORD_( unit, field, size(data(:,:,:)), data, tindex )
          else
              lenx=size(data,1)
              leny=size(data,2)
              lenz=size(data,3)
              len=lenx*leny*lenz
              allocate(gdata(len))          
! read field on pe 0 and pass to all pes
              if( pe.EQ.0 ) call READ_RECORD_( unit, field, len, gdata, tindex )
! broadcasting global array, this can be expensive!          
              call mpp_transmit( put_data=gdata(1), plen=len, to_pe=ALL_PES, &
                                 get_data=gdata(1), glen=len, from_pe=0 )
              ioff = is; joff = js
              if( data_has_halos )then
                  ioff = isd; joff = jsd
              end if
              do k=1,size(data,3)
                 do j=js,je
                    do i=is,ie
                       n=(i-isg+1) + (j-jsg)*lenx + (k-1)*lenx*leny
                       data(i-ioff+1,j-joff+1,k)=gdata(n)
                    enddo
                 enddo
              enddo
              deallocate(gdata)
              call mpp_sync_self() ! ensure MPI_ISEND is done.
          end if
      else if( data_has_halos )then
! for uniprocessor or multithreaded read
! read compute domain as contiguous data

          allocate( cdata(is:ie,js:je,size(data,3)) )
          call READ_RECORD_(unit,field,size(cdata(:,:,:)),cdata,tindex,domain,position,tile_count)

          data(is-isd+1:ie-isd+1,js-jsd+1:je-jsd+1,:) = cdata(:,:,:)
          deallocate(cdata)
      else
          call READ_RECORD_(unit,field,size(data(:,:,:)),data,tindex,domain,position,tile_count)
      end if

      call mpp_clock_end(mpp_read_clock)

      return
    end subroutine MPP_READ_2DDECOMP_3D_


    subroutine MPP_READ_2DDECOMP_4D_( unit, field, domain, data, tindex, tile_count )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      integer, intent(in), optional :: tindex, tile_count
      MPP_TYPE_ :: data3D(size(data,1),size(data,2),size(data,3)*size(data,4))
      pointer( ptr, data3D )
      ptr = LOC(data)
      call mpp_read( unit, field, domain, data3D, tindex, tile_count)
      return
    end subroutine MPP_READ_2DDECOMP_4D_

