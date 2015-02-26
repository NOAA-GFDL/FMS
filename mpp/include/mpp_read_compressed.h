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

      call mpp_clock_begin(mpp_read_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: invalid unit number.' )

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
      call mpp_clock_end(mpp_read_clock)
      return
    end subroutine MPP_READ_COMPRESSED_2D_
