!----------
!ug support

!------------------------------------------------------------------------------
!>Write data for a 1D field associated with an unstructured mpp domain to a
!!restart file.
subroutine mpp_io_unstructured_write_r_1D(funit, &
                                          field, &
                                          domain, &
                                          fdata, &
                                          nelems_io_in, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)                       :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)                      :: field        !<A field whose data will be written
    type(domainUG),intent(inout)                       :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:),intent(inout)                    :: fdata        !<The data that will be written to the file.
    integer(INT_KIND),dimension(:),intent(in),optional :: nelems_io_in !<Number of grid points in the unstructured dimension for each rank (correct sizes
                                                                       !!only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional                           :: tstamp       !<A time value.
    real,intent(in), optional                          :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill           !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND),dimension(:),allocatable :: nelems_io      !<Number of grid points in the unstructured dimension for each rank.
    integer(INT_KIND)                          :: compute_size   !<Size of the unstructured compute domain for the current rank.
    integer(INT_KIND)                          :: nelems         !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:),allocatable              :: cdata          !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: i              !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_1D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_1D:" &
                       //" the inputted file unit is already in use.")
    endif

   !Set the fill value for the field.
    fill = 0.0
    if (present(default_data)) then
        fill = default_data
    endif

   !Point to the I/O domain associated with the inputted unstructured mpp
   !domain.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)

   !Get the pelist associated with the I/O domain.
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)

   !Make sure that only the root rank of the pelist will write to the file.
   !This check is needed because data is only gathered on the lowest rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1) .and. .not. &
        mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_1D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_1D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
    endif

   !For the 1D unstructured case, data is assumed to be of the form
   !fdata = fdata(unstructured).  The number of data elements in the
   !unstructured dimension (size(fdata,1)) may differ between ranks.
   !If not passed in, the number of data elements in the unstructured
   !dimension must be gathered on the root rank of the pelist.  The number
   !data elements in the unstructured dimension should be equal to the size
   !of the unstructured computed domain.
    if (present(nelems_io_in)) then
        allocate(nelems_io(size(nelems_io_in)))
        nelems_io = nelems_io_in
    else
        allocate(nelems_io(io_domain_npes))
        nelems_io = 0
        call mpp_get_UG_compute_domain(domain, &
                                       size=compute_size)
        if (compute_size .ne. size(fdata)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_1D:" &
                           //" the size of the data does not match the" &
                           //" the size of the unstructured compute domain.")
        endif
        call mpp_gather((/compute_size/), &
                        nelems_io, &
                        pelist)
    endif

   !Allocate an array which will be used to gather the data to be written
   !onto the root rank of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems))
    else
        allocate(rbuff(1))
    endif

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(fdata, &
                    size(fdata), &
                    rbuff, &
                    nelems_io, &
                    pelist)

   !Write out the data to the file.  This is only done by the root rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems))
        cdata = fill
        do i = 1,nelems
            cdata(i) = rbuff(i)
        enddo
        field%size(1) = nelems
        call write_record(funit, &
                          field, &
                          nelems, &
                          cdata, &
                          tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(rbuff)
    deallocate(pelist)
    deallocate(nelems_io)

   !Stop the mpp timer.
    call mpp_clock_end(mpp_write_clock)

    return
end subroutine mpp_io_unstructured_write_r_1D

!------------------------------------------------------------------------------
!>Write data for a 2D field associated with an unstructured mpp domain to a
!!restart file.
subroutine mpp_io_unstructured_write_r_2D(funit, &
                                          field, &
                                          domain, &
                                          fdata, &
                                          nelems_io_in, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)             :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)            :: field        !<A field whose data will be written
    type(domainUG),intent(inout)             :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:,:),intent(inout)        :: fdata        !<The data that will be written to the file.
    integer,dimension(:),intent(in),optional :: nelems_io_in !<Number of grid points in the unstructured dimension for each rank (correct
                                                             !!sizes only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional                 :: tstamp       !<A time value.
    real,intent(in), optional                :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill             !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain        !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes   !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist           !<A pelist.
    integer(INT_KIND),dimension(:),allocatable :: nelems_io        !<Number of grid points in the unstructured dimension for each rank.
    integer(INT_KIND)                          :: compute_size     !<Size of the unstructured compute domain for the current rank.
    integer(INT_KIND),dimension(:),allocatable :: dim_size_gather  !<Array used to check that size(fdata,2) (non-unstructured dimension) is the same for all ranks in an I/O domain.
    integer(INT_KIND)                          :: size_fdata_dim_2 !<Number of data points in the non-unstructured dimension (size(fdata,2)).
    integer(INT_KIND)                          :: mynelems         !<Number of data points in the unstructured dimension (size(fdata,1)).
    real,dimension(:),allocatable              :: sbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    integer(INT_KIND)                          :: nelems           !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:,:),allocatable            :: cdata            !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: i                !<Loop variable.
    integer(INT_KIND)                          :: j                !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_2D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_2D:" &
                       //" the inputted file unit is already in use.")
    endif

   !Set the fill value for the field.
    fill = 0.0
    if (present(default_data)) then
        fill = default_data
    endif

   !Point to the I/O domain associated with the inputted unstructured mpp
   !domain.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)

   !Get the pelist associated with the I/O domain.
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)

   !Make sure that only the root rank of the pelist will write to the file.
   !This check is needed because data is only gathered on the lowest rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1) .and. .not. &
        mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_2D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_2D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
    endif

   !For the 2D unstructured case, data is assumed to be of the form
   !fdata = fdata(unstructured,z).  The number of data elements in the
   !unstructured dimension (size(fdata,1)) may differ between ranks.
   !If not passed in, the number of data elements in the unstructured
   !dimension must be gathered on the root rank of the pelist.  The number
   !data elements in the unstructured dimension should be equal to the size
   !of the unstructured computed domain.
    if (present(nelems_io_in)) then
        allocate(nelems_io(size(nelems_io_in)))
        nelems_io = nelems_io_in
    else
        allocate(nelems_io(io_domain_npes))
        nelems_io = 0
        call mpp_get_UG_compute_domain(domain, &
                                       size=compute_size)
        if (compute_size .ne. size(fdata,1)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_2D:" &
                           //" the size of the data does not match the" &
                           //" the size of the unstructured compute domain.")
        endif
        call mpp_gather((/compute_size/), &
                        nelems_io, &
                        pelist)
    endif

   !For the non-unstructured dimension, the number of data elements is
   !required to be the same for all ranks.  Perform a gather to check this.
    allocate(dim_size_gather(io_domain_npes))
    dim_size_gather = 0
    size_fdata_dim_2 = size(fdata,2)
    call mpp_gather((/size_fdata_dim_2/), &
                    dim_size_gather, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxloc(dim_size_gather,1) .ne. minloc(dim_size_gather,1)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_2D:" &
                           //" all ranks in the io_domain must have the same" &
                           //" number of data elements in the z-dimension" &
                           //" (size(data,2)).")
        endif
    endif
    deallocate(dim_size_gather)

   !Allocate arrays which will be used to gather the data to be written
   !onto the root rank of the pelist.
    mynelems = size(fdata,1)
    allocate(sbuff(mynelems*size_fdata_dim_2))
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems*size_fdata_dim_2))
    else
        allocate(rbuff(1))
    endif

   !Load the data into the sbuff array.  The data is transposed so that the
   !gather may be performed more easily.
    do j = 1,mynelems
        do i = 1,size_fdata_dim_2
            sbuff((j-1)*size_fdata_dim_2+i) = fdata(j,i)
        enddo
    enddo

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(sbuff, &
                    size(sbuff), &
                    rbuff, &
                    size_fdata_dim_2*nelems_io, &
                    pelist)

   !Write out the data to the file.  This is only done by the root rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems,size_fdata_dim_2))
        cdata = fill
        do j = 1,size_fdata_dim_2
            do i = 1,nelems
                cdata(i,j) = rbuff((i-1)*size_fdata_dim_2+j)
            enddo
        enddo
        field%size(1) = nelems
        call write_record(funit, &
                          field, &
                          nelems*size_fdata_dim_2, &
                          cdata, &
                          tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(sbuff)
    deallocate(rbuff)
    deallocate(pelist)
    deallocate(nelems_io)

   !Stop the mpp timer.
    call mpp_clock_end(mpp_write_clock)

    return
end subroutine mpp_io_unstructured_write_r_2D

!------------------------------------------------------------------------------
!>Write data for a 3D field associated with an unstructured mpp domain to a
!!restart file.
subroutine mpp_io_unstructured_write_r_3D(funit, &
                                          field, &
                                          domain, &
                                          fdata, &
                                          nelems_io_in, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)             :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)            :: field        !<A field whose data will be written
    type(domainUG),intent(inout)             :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:,:,:),intent(inout)      :: fdata        !<The data that will be written to the file.
    integer,dimension(:),intent(in),optional :: nelems_io_in !<Number of grid points in the unstructured dimension for each rank (correct
                                                             !!sizes only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional                 :: tstamp       !<A time value.
    real,intent(in), optional                :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill             !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain        !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes   !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist           !<A pelist.
    integer(INT_KIND),dimension(:),allocatable :: nelems_io        !<Number of grid points in the unstructured dimension for each rank.
    integer(INT_KIND)                          :: compute_size     !<Size of the unstructured compute domain for the current rank.
    integer(INT_KIND),dimension(:),allocatable :: dim_size_gather  !<Array used to check that size(fdata,2) and size(fdata,3) (non-unstructured dimensions) are the same for
                                                                   !!all ranks in an I/O domain.
    integer(INT_KIND)                          :: size_fdata_dim_2 !<Number of data points in a non-unstructured dimension (size(fdata,2)).
    integer(INT_KIND)                          :: size_fdata_dim_3 !<Number of data points in a non-unstructured dimension (size(fdata,3)).
    integer(INT_KIND)                          :: mynelems         !<Number of data points in the unstructured dimension (size(fdata,1)).
    real,dimension(:),allocatable              :: sbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    integer(INT_KIND)                          :: nelems           !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:,:,:),allocatable          :: cdata            !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: i                !<Loop variable.
    integer(INT_KIND)                          :: j                !<Loop variable.
    integer(INT_KIND)                          :: k                !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_3D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_3D:" &
                       //" the inputted file unit is already in use.")
    endif

   !Set the fill value for the field.
    fill = 0.0
    if (present(default_data)) then
        fill = default_data
    endif

   !Point to the I/O domain associated with the inputted unstructured mpp
   !domain.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)

   !Get the pelist associated with the I/O domain.
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)

   !Make sure that only the root rank of the pelist will write to the file.
   !This check is needed because data is only gathered on the lowest rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1) .and. .not. &
        mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_3D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_3D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
    endif

   !For the 3D unstructured case, data is assumed to be of the form
   !fdata = fdata(unstructured,z,cc).  The number of data elements in the
   !unstructured dimension (size(fdata,1)) may differ between ranks.
   !If not passed in, the number of data elements in the unstructured
   !dimension must be gathered on the root rank of the pelist.  The number
   !data elements in the unstructured dimension should be equal to the size
   !of the unstructured computed domain.
    if (present(nelems_io_in)) then
        allocate(nelems_io(size(nelems_io_in)))
        nelems_io = nelems_io_in
    else
        allocate(nelems_io(io_domain_npes))
        nelems_io = 0
        call mpp_get_UG_compute_domain(domain, &
                                       size=compute_size)
        if (compute_size .ne. size(fdata,1)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_3D:" &
                           //" the size of the data does not match the" &
                           //" the size of the unstructured compute domain.")
        endif
        call mpp_gather((/compute_size/), &
                        nelems_io, &
                        pelist)
    endif

   !The number of data elements in the non-unstructured dimensions are
   !required to be the same for all ranks.  Perform gathers to check this.
    allocate(dim_size_gather(io_domain_npes))
    dim_size_gather = 0
    size_fdata_dim_2 = size(fdata,2)
    call mpp_gather((/size_fdata_dim_2/), &
                    dim_size_gather, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxloc(dim_size_gather,1) .ne. minloc(dim_size_gather,1)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_3D:" &
                           //" all ranks in the io_domain must have the same" &
                           //" number of data elements in the z-dimension" &
                           //" (size(data,2)).")
        endif
    endif
    dim_size_gather = 0
    size_fdata_dim_3 = size(fdata,3)
    call mpp_gather((/size_fdata_dim_3/), &
                    dim_size_gather, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxloc(dim_size_gather,1) .ne. minloc(dim_size_gather,1)) then
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_write_compressed_r_3D:" &
                           //" all ranks in the io_domain must have the same" &
                           //" number of data elements in the cc-dimension" &
                           //" (size(data,3)).")
        endif
    endif
    deallocate(dim_size_gather)

   !Allocate arrays which will be used to gather the data to be written
   !onto the root rank of the pelist.
    mynelems = size(fdata,1)
    allocate(sbuff(mynelems*size_fdata_dim_2*size_fdata_dim_3))
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems*size_fdata_dim_2*size_fdata_dim_3))
    else
        allocate(rbuff(1))
    endif

   !Load the data into the sbuff array.  The data is transposed so that the
   !gather may be performed more easily.
    do k = 1,mynelems
        do j = 1,size_fdata_dim_2
            do i = 1,size_fdata_dim_3
                sbuff((k-1)*size_fdata_dim_2*size_fdata_dim_3+(j-1)*size_fdata_dim_3+i) = fdata(k,j,i)
            enddo
        enddo
    enddo

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(sbuff, &
                    size(sbuff), &
                    rbuff, &
                    nelems_io*size_fdata_dim_2*size_fdata_dim_3, &
                    pelist)

   !Write out the data to the file.  This is only done by the root rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems,size_fdata_dim_2,size_fdata_dim_3))
        cdata = fill
        do k = 1,size_fdata_dim_3
            do j = 1,size_fdata_dim_2
                do i = 1,nelems
                    cdata(i,j,k) = rbuff((i-1)*size_fdata_dim_2*size_fdata_dim_3+(j-1)*size_fdata_dim_3+k)
                enddo
            enddo
        enddo
        field%size(1) = nelems
        call write_record(funit, &
                          field, &
                          nelems*size_fdata_dim_2*size_fdata_dim_3, &
                          cdata, &
                          tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(sbuff)
    deallocate(rbuff)
    deallocate(pelist)
    deallocate(nelems_io)

   !Stop the mpp timer.
    call mpp_clock_end(mpp_write_clock)

    return
end subroutine mpp_io_unstructured_write_r_3D

!------------------------------------------------------------------------------

!----------
