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
    subroutine MPP_READ_COMPRESSED_1D_(unit, field, domain, data, tindex)
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      type(domain2D),  intent(in) :: domain
      MPP_TYPE_, intent(inout) :: data(:)
      integer,           intent(in), optional :: tindex

      MPP_TYPE_ :: data2D(size(data,1),1)
      pointer( ptr, data2D )
      ptr = LOC(data)

      call mpp_read(unit, field, domain, data2D, tindex)
      return
    end subroutine MPP_READ_COMPRESSED_1D_

    subroutine MPP_READ_COMPRESSED_2D_(unit, field, domain, data, tindex, start, nread, threading)
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(in)           :: field
      type(domain2D),    intent(in)           :: domain
      MPP_TYPE_,         intent(inout)        :: data(:,:)
      integer,           intent(in), optional :: tindex
      integer,           intent(in), optional :: start(:), nread(:)
      integer,           intent(in), optional :: threading

      integer, allocatable :: pelist(:)
      integer :: npes, p, threading_flag
      type(domain2d), pointer :: io_domain=>NULL()
      logical :: compute_chksum,print_compressed_chksum
      integer(LONG_KIND) ::chk

      call mpp_clock_begin(mpp_read_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: invalid unit number.' )

      print_compressed_chksum = .FALSE.

    if(size(data) > 0) then
      data = 0 !! zero out data so other tiles do not contribute junk to chksum
      threading_flag = MPP_SINGLE
      if( PRESENT(threading) )threading_flag = threading
      if( threading_flag == MPP_MULTI ) then
	 call read_record(unit,field,size(data(:,:)),data,tindex,start_in=start, axsiz_in=nread)
      else if( threading_flag == MPP_SINGLE ) then

	 io_domain=>mpp_get_io_domain(domain)
	 if(.not. ASSOCIATED(io_domain)) call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: io_domain must be defined.' )
	 npes = mpp_get_domain_npes(io_domain)
	 allocate(pelist(npes))
	 call mpp_get_pelist(io_domain,pelist)

	 if(mpp_pe() == pelist(1)) call read_record(unit,field,size(data(:,:)),data,tindex,start_in=start, axsiz_in=nread)

	 !--- z1l replace mpp_broadcast with mpp_send/mpp_recv to avoid hang in calling MPI_COMM_CREATE
	 !---     because size(pelist) might be different for different rank.
	 !--- prepost receive
	 if( mpp_pe() == pelist(1) ) then
	    do p = 2, npes
	       call mpp_send(data(1,1), plen=size(data(:,:)), to_pe=pelist(p), tag=COMM_TAG_1)
	    enddo
	    call mpp_sync_self()
	 else
	    call mpp_recv(data(1,1), glen=size(data(:,:)), from_pe=pelist(1), block=.false., tag=COMM_TAG_1)
	    call mpp_sync_self(check=EVENT_RECV)
	 endif

	 deallocate(pelist)
      else
	 call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: threading should be MPP_SINGLE or MPP_MULTI')
      endif
    endif

      compute_chksum = .FALSE.
      if (ANY(field%checksum /= default_field%checksum) ) compute_chksum = .TRUE.

      if (compute_chksum) then
#ifdef use_netCDF
	 if (field%type==NF_INT) then
            if (field%fill == MPP_FILL_DOUBLE .or. field%fill == real(MPP_FILL_INT) ) then
               chk = mpp_chksum( ceiling(data), mask_val=MPP_FILL_INT )
            else
	       call mpp_error(NOTE,"During mpp_io(mpp_read_compressed_2d) int field "//trim(field%name)// &
			      " found fill. Icebergs, or code using defaults can safely ignore. "// &
			      " If manually overriding compressed restart fills, confirm this is what you want.")
	       chk = mpp_chksum( ceiling(data), mask_val=field%fill)
	    end if
	 else !!real data
	    chk = mpp_chksum(data,mask_val=field%fill)
	 end if
#endif
	 !!compare
         if ( print_compressed_chksum) then
            if ( mpp_pe() == mpp_root_pe() ) then 
               print '(A,Z16)', "mpp_read_compressed_2d chksum: "//trim(field%name)//" = ", chk
               !! discuss making fatal after testing/review to match other routines.
               !! Need to do some nword-counting and digging with pjp
               !! this should be if ( chk /= field%checksum ) as it was tested @ulm_201505..
               if ( MOD(chk, field%checksum(1)) /= 0 ) then
                  print '(A,Z16)', "File stored checksum: "//trim(field%name)//" = ", field%checksum(1)
                  call mpp_error(NOTE,"mpp_read_compressed_2d chksum: "//trim(field%name)//" failed!")
               end if
            endif
         end if
      end if

      call mpp_clock_end(mpp_read_clock)
      return
    end subroutine MPP_READ_COMPRESSED_2D_

    subroutine MPP_READ_COMPRESSED_3D_(unit, field, domain, data, tindex, start, nread, threading)
      integer,		 intent(in)	      :: unit
      type(fieldtype),	 intent(in)	      :: field
      type(domain2D),	 intent(in)	      :: domain
      MPP_TYPE_,	 intent(inout)	      :: data(:,:,:)
      integer,		 intent(in), optional :: tindex
      integer,		 intent(in), optional :: start(:), nread(:)
      integer,		 intent(in), optional :: threading

      integer, allocatable :: pelist(:)
      integer :: npes, p, threading_flag
      type(domain2d), pointer :: io_domain=>NULL()
      logical :: compute_chksum,print_compressed_chksum
      integer(LONG_KIND) ::chk

      call mpp_clock_begin(mpp_read_clock)

      data = 0 !! zero out data so other tiles do not contribute junk to chksum

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_3D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_3D_: invalid unit number.' )

      print_compressed_chksum = .FALSE.
      threading_flag = MPP_SINGLE
      if( PRESENT(threading) )threading_flag = threading
      if( threading_flag == MPP_MULTI ) then
	 call read_record(unit,field,size(data(:,:,:)),data,tindex,start_in=start, axsiz_in=nread)
      else if( threading_flag == MPP_SINGLE ) then

	 io_domain=>mpp_get_io_domain(domain)
	 if(.not. ASSOCIATED(io_domain)) call mpp_error( FATAL, 'MPP_READ_COMPRESSED_3D_: io_domain must be defined.' )
	 npes = mpp_get_domain_npes(io_domain)
	 allocate(pelist(npes))
	 call mpp_get_pelist(io_domain,pelist)

	 if(mpp_pe() == pelist(1)) call read_record(unit,field,size(data(:,:,:)),data,tindex,start_in=start, axsiz_in=nread)

	 !--- z1l replace mpp_broadcast with mpp_send/mpp_recv to avoid hang in calling MPI_COMM_CREATE
	 !---	  because size(pelist) might be different for different rank.
	 !--- prepost receive
	 if( mpp_pe() == pelist(1) ) then
	    do p = 2, npes
	       call mpp_send(data(1,1,1), plen=size(data(:,:,:)), to_pe=pelist(p), tag=COMM_TAG_1)
	    enddo
	    call mpp_sync_self()
	 else
	    call mpp_recv(data(1,1,1), glen=size(data(:,:,:)), from_pe=pelist(1), block=.false., tag=COMM_TAG_1)
	    call mpp_sync_self(check=EVENT_RECV)
	 endif

	 deallocate(pelist)
      else
	 call mpp_error( FATAL, 'MPP_READ_COMPRESSED_3D_: threading should be MPP_SINGLE or MPP_MULTI')
      endif

      compute_chksum = .FALSE.
      if (ANY(field%checksum /= default_field%checksum) ) compute_chksum = .TRUE.

      if (compute_chksum) then
#ifdef use_netCDF
	 if (field%type==NF_INT) then
	    if (field%fill == MPP_FILL_DOUBLE .or. field%fill == real(MPP_FILL_INT) ) then
               chk = mpp_chksum( ceiling(data), mask_val=MPP_FILL_INT )
            else 
	       call mpp_error(NOTE,"During mpp_io(mpp_read_compressed_3d) int field "//trim(field%name)// &
			      " found fill. Icebergs, or code using defaults can safely ignore. "// &
			      " If manually overriding compressed restart fills, confirm this is what you want.")
	       chk = mpp_chksum( ceiling(data), mask_val=field%fill)
	    end if
	 else !!real
	    chk = mpp_chksum(data,mask_val=field%fill)
	 end if
#endif
	 !!compare
         if ( print_compressed_chksum) then
            if ( mpp_pe() == mpp_root_pe() ) then 
               print '(A,Z16)', "mpp_read_compressed_3d chksum: "//trim(field%name)//" = ", chk
               !! discuss making fatal after testing/review to match other routines.
               !! Need to do some nword-counting and digging with pjp
               !! this should be if ( chk /= field%checksum ) as it was tested @ulm_201505..
               if ( MOD(chk, field%checksum(1)) /= 0 ) then
                  print '(A,Z16)', "File stored checksum: "//trim(field%name)//" = ", field%checksum(1)
                  call mpp_error(NOTE,"mpp_read_compressed_3d chksum: "//trim(field%name)//" failed!")
               end if
            endif
         end if
      end if

      call mpp_clock_end(mpp_read_clock)
      return
    end subroutine MPP_READ_COMPRESSED_3D_
