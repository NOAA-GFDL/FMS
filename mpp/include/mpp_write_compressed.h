    subroutine MPP_WRITE_COMPRESSED_1D_(unit, field, domain, data, nelems_io, tstamp, default_data)
      integer, intent(in) :: unit
      type(fieldtype), intent(inout) :: field
      type(domain2D), intent(inout) :: domain
      MPP_TYPE_, intent(inout) :: data(:)
      integer, intent(in) :: nelems_io(:)  ! number of compressed elements
      real,              intent(in), optional :: tstamp
      MPP_TYPE_,         intent(in), optional :: default_data

      MPP_TYPE_ :: data2D(size(data,1),1)
      pointer( ptr, data2D )
      ptr = LOC(data)

      call mpp_write_compressed(unit, field, domain, data2D, nelems_io, tstamp, default_data)
      return
    end subroutine MPP_WRITE_COMPRESSED_1D_

    subroutine MPP_WRITE_COMPRESSED_2D_(unit, field, domain, data, nelems_io, tstamp, default_data)
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(inout)        :: field
      type(domain2D),    intent(inout)        :: domain
      MPP_TYPE_,         intent(inout)        :: data(:,:)
      integer,           intent(in)           :: nelems_io(:)  ! number of compressed elements from each
                                                               ! member of the io_domain. It MUST have the
                                                               ! same order as the io_domain pelist.
      real,              intent(in), optional :: tstamp
      MPP_TYPE_,         intent(in), optional :: default_data

!cdata is used to store the io-domain compressed data
      MPP_TYPE_, allocatable, dimension(:,:) :: cdata
      MPP_TYPE_, allocatable, dimension(:,:) :: sbuff,rbuff
      MPP_TYPE_                            :: fill

      MPP_TYPE_ :: sbuff1D(size(data))
      MPP_TYPE_ :: rbuff1D(size(data,2)*sum(nelems_io(:)))
      pointer(sptr,sbuff1D); pointer(rptr,rbuff1D)
      integer, allocatable :: pelist(:)
      integer, allocatable :: nz_gather(:)
      integer :: i,j,nz,nelems,mynelems,idx,npes
      type(domain2d), pointer :: io_domain=>NULL()

      call mpp_clock_begin(mpp_write_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE_COMPRESSED_2D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_WRITE_COMPRESSED_2D_: invalid unit number.' )

      fill = 0
      if(PRESENT(default_data)) fill = default_data

      io_domain=>mpp_get_io_domain(domain) 
      if (.not. ASSOCIATED(io_domain)) call mpp_error( FATAL, 'MPP_WRITE_COMPRESSED_2D_: io_domain must be defined.' )
      npes = mpp_get_domain_npes(io_domain)
      allocate(pelist(npes))
      call mpp_get_pelist(io_domain,pelist)

      mynelems = size(data,1)
      nz = size(data,2)
      nelems = sum(nelems_io(:))
      ! Check that nz is consistent across all PEs in io_domain
      allocate(nz_gather(npes))
      call mpp_gather((/nz/), nz_gather, pelist)
      if ( mpp_file(unit)%write_on_this_pe.and.maxloc(nz_gather,1).ne.minloc(nz_gather,1) ) then
         call mpp_error( FATAL, 'MPP_WRITE_COMPRESSED_2D_: size(data,2) must be consistent across all PEs in io_domain' )
      end if
      deallocate(nz_gather)

      if(mpp_file(unit)%write_on_this_pe ) allocate(rbuff(nz,nelems))
      allocate(sbuff(nz,mynelems))

      ! this matrix inversion makes for easy gather to the IO root
      ! and a clear, concise unpack
      do j=1,mynelems
        do i=1,nz
          sbuff(i,j) = data(j,i) 
      enddo; enddo

   !  Note that the gatherV implied here is asymmetric; only root needs to know the vector of recv size
      sptr = LOC(sbuff); rptr = LOC(rbuff)
      call mpp_gather(sbuff1D,size(sbuff),rbuff1D,nz*nelems_io(:),pelist)

      if(mpp_file(unit)%write_on_this_pe ) then
         allocate(cdata(nelems,nz))
         cdata = fill
         do j=1,nz
           do i=1,nelems
             cdata(i,j) = rbuff(j,i)
         enddo; enddo
         ! cludge for now; need resizing accessor
         field%size(1) = nelems
         call write_record( unit, field, nelems*nz, cdata, tstamp)
         deallocate(rbuff,cdata)
      endif

      deallocate(sbuff,pelist)

      call mpp_clock_end(mpp_write_clock)

      return
    end subroutine MPP_WRITE_COMPRESSED_2D_
