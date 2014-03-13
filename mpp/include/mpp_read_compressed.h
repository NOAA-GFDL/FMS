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

    subroutine MPP_READ_COMPRESSED_2D_(unit, field, domain, data, tindex)
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(in)           :: field
      type(domain2D),    intent(in)           :: domain
      MPP_TYPE_,         intent(inout)        :: data(:,:)
      integer,              intent(in), optional :: tindex

      integer, allocatable :: pelist(:)
      integer :: npes
      type(domain2d), pointer :: io_domain=>NULL()

      call mpp_clock_begin(mpp_read_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: invalid unit number.' )

      io_domain=>mpp_get_io_domain(domain) 
      if(.not. ASSOCIATED(io_domain)) call mpp_error( FATAL, 'MPP_READ_COMPRESSED_2D_: io_domain must be defined.' )
      npes = mpp_get_domain_npes(io_domain)
      allocate(pelist(npes))
      call mpp_get_pelist(io_domain,pelist)

      if(mpp_pe() == pelist(1)) call read_record(unit,field,size(data(:,:)),data,tindex)

      call mpp_broadcast(data,size(data(:,:)),pelist(1),pelist)
      deallocate(pelist)
      call mpp_clock_end(mpp_read_clock)
      return
    end subroutine MPP_READ_COMPRESSED_2D_
