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

!----------
!ug support

!------------------------------------------------------------------------------
!>Write data for a 1D field associated with an unstructured mpp domain to a
!!restart file.
subroutine mpp_io_unstructured_write_r_1D(funit, &
                                          field, &
                                          domain, &
                                          fdata, &
                                          nelems_io, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)    :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)   :: field        !<A field whose data will be written.
    type(domainUG),intent(inout)    :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:),intent(inout) :: fdata        !<The data that will be written to the file.
    integer,dimension(:),intent(in) :: nelems_io    !<Number of grid points in the compressed dimension for each rank (correct
                                                    !!sizes only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional        :: tstamp       !<A time value.
    real,intent(in), optional       :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill           !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: nelems         !<Total number of data points (sum(nelems_io)) to be written by the root rank of the pelist.
    real,dimension(:),allocatable              :: rbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:),allocatable              :: cdata          !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: i              !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_1D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_1D:" &
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
                       "mpp_io_unstructured_write_r_1D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_1D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
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
        call write_record_default(funit, &
                                  field, &
                                  nelems, &
                                  cdata, &
                                  tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(rbuff)
    deallocate(pelist)

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
                                          nelems_io, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)      :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)     :: field        !<A field whose data will be written
    type(domainUG),intent(inout)      :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:,:),intent(inout) :: fdata        !<The data that will be written to the file.
    integer,dimension(:),intent(in)   :: nelems_io    !<Number of grid points in the compressed dimension for each rank (correct
                                                      !!sizes only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional          :: tstamp       !<A time value.
    real,intent(in), optional         :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill           !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: dim_size_1     !<Number of data points in the first dimension (size(fdata,1)).
    integer(INT_KIND)                          :: dim_size_2     !<Number of data points in the second dimension (size(fdata,2)).
    real,dimension(:),allocatable              :: sbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    integer(INT_KIND)                          :: nelems         !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:,:),allocatable            :: cdata          !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: offset_r       !<Offset for the rbuff array used to reorder the data before the netCDF write.
    integer(INT_KIND)                          :: offset_c       !<Offset for the cdata array used to reorder the data before the netCDF write.
    integer(INT_KIND)                          :: i              !<Loop variable.
    integer(INT_KIND)                          :: j              !<Loop variable.
    integer(INT_KIND)                          :: k              !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_2D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_2D:" &
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
                       "mpp_io_unstructured_write_r_2D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_2D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
    endif

   !Load the data elements for each rank into a one dimensional array, which
   !will be used to gather the data onto the root rank of the pelist.
    allocate(sbuff(size(fdata)))
    dim_size_1 = size(fdata,1)
    dim_size_2 = size(fdata,2)
    do j = 1,dim_size_2
        do i = 1,dim_size_1
            sbuff((j-1)*dim_size_1+i) = fdata(i,j)
        enddo
    enddo

   !Allocate an array which will be used to gather the data to be written
   !onto the root rank of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems*dim_size_2))
    else
        allocate(rbuff(1))
    endif

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(sbuff, &
                    size(sbuff), &
                    rbuff, &
                    nelems_io*dim_size_2, &
                    pelist)

   !Reorder the gather data so that is of the form (nelems,dim_size_2).  Write
   !out the data to the file.  This is only done by the root rank of the
   !pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems,dim_size_2))
        cdata = fill
        do j = 1,dim_size_2
            offset_c = 0
            do k = 1,io_domain_npes
                if (k .gt. 1) then
                    offset_r = (j-1)*nelems_io(k) + dim_size_2*(sum(nelems_io(1:k-1)))
                else
                    offset_r = (j-1)*nelems_io(k)
                endif
                do i = 1,nelems_io(k)
                    cdata(i+offset_c,j) = rbuff(i+offset_r)
                enddo
                offset_c = offset_c + nelems_io(k)
            enddo
        enddo
        field%size(1) = nelems
        call write_record_default(funit, &
                                  field, &
                                  nelems*dim_size_2, &
                                  cdata, &
                                  tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(sbuff)
    deallocate(rbuff)
    deallocate(pelist)

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
                                          nelems_io, &
                                          tstamp, &
                                          default_data)

   !Inputs/outputs
    integer(INT_KIND),intent(in)        :: funit        !<A file unit for the to which the data will be written.
    type(fieldtype),intent(inout)       :: field        !<A field whose data will be written
    type(domainUG),intent(inout)        :: domain       !<An unstructured mpp domain associatd with the inputted file.
    real,dimension(:,:,:),intent(inout) :: fdata        !<The data that will be written to the file.
    integer,dimension(:),intent(in)     :: nelems_io    !<Number of grid points in the compressed dimension for each rank (correct
                                                        !!sizes only exist for the root rank of I/O domain pelist.)
    real,intent(in),optional            :: tstamp       !<A time value.
    real,intent(in), optional           :: default_data !<Fill value for the inputted field.

   !Local variables
    real                                       :: fill           !<Fill value for the inputted field.  This defaults to zero.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: dim_size_1     !<Number of data points in the first dimension (size(fdata,1)).
    integer(INT_KIND)                          :: dim_size_2     !<Number of data points in the second dimension (size(fdata,2)).
    integer(INT_KIND)                          :: dim_size_3     !<Number of data points in the second dimension (size(fdata,3)).
    real,dimension(:),allocatable              :: sbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    integer(INT_KIND)                          :: nelems         !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff          !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:,:,:),allocatable          :: cdata          !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: offset_r       !<Offset for the rbuff array used to reorder the data before the netCDF write.
    integer(INT_KIND)                          :: offset_c       !<Offset for the cdata array used to reorder the data before the netCDF write.
    integer(INT_KIND)                          :: i              !<Loop variable.
    integer(INT_KIND)                          :: j              !<Loop variable.
    integer(INT_KIND)                          :: k              !<Loop variable.
    integer(INT_KIND)                          :: m              !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_3D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_3D:" &
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
                       "mpp_io_unstructured_write_r_3D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_r_3D:" &
                       //" the non-root ranks of the pelist are not" &
                       //" allowed to perform the write.")
    endif

   !Load the data elements for each rank into a one dimensional array, which
   !will be used to gather the data onto the root rank of the pelist.
    allocate(sbuff(size(fdata)))
    dim_size_1 = size(fdata,1)
    dim_size_2 = size(fdata,2)
    dim_size_3 = size(fdata,3)
    do k = 1,dim_size_3
        do j = 1,dim_size_2
            do i = 1,dim_size_1
                sbuff((k-1)*dim_size_2*dim_size_1+(j-1)*dim_size_1+i) = fdata(i,j,k)
            enddo
        enddo
    enddo

   !Allocate an array which will be used to gather the data to be written
   !onto the root rank of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems*dim_size_2*dim_size_3))
    else
        allocate(rbuff(1))
    endif

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(sbuff, &
                    size(sbuff), &
                    rbuff, &
                    nelems_io*dim_size_2*dim_size_3, &
                    pelist)

   !Reorder the gather data so that is of the form (nelems,dim_size_2).  Write
   !out the data to the file.  This is only done by the root rank of the
   !pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems,dim_size_2,dim_size_3))
        cdata = fill
        do m = 1,dim_size_3
            do j = 1,dim_size_2
                offset_c = 0
                do k = 1,io_domain_npes
                    if (k .gt. 1) then
                        offset_r = (m-1)*dim_size_2*nelems_io(k) + &
                                   (j-1)*nelems_io(k) + &
                                   dim_size_2*dim_size_3*(sum(nelems_io(1:k-1)))
                    else
                        offset_r = (m-1)*dim_size_2*nelems_io(k) + &
                                   (j-1)*nelems_io(k)
                    endif
                    do i = 1,nelems_io(k)
                        cdata(i+offset_c,j,m) = rbuff(i+offset_r)
                    enddo
                    offset_c = offset_c + nelems_io(k)
                enddo
            enddo
        enddo
        field%size(1) = nelems
        call write_record_default(funit, &
                                  field, &
                                  nelems*dim_size_2*dim_size_3, &
                                  cdata, &
                                  tstamp)
        deallocate(cdata)
    endif

   !Deallocate local allocatables.
    deallocate(sbuff)
    deallocate(rbuff)
    deallocate(pelist)

   !Stop the mpp timer.
    call mpp_clock_end(mpp_write_clock)

    return
end subroutine mpp_io_unstructured_write_r_3D

!------------------------------------------------------------------------------
!>Write data for a 4D field associated with an unstructured mpp domain to a
!!restart file.
subroutine mpp_io_unstructured_write_r_4D(funit, &
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
    real,dimension(:,:,:,:),intent(inout)    :: fdata        !<The data that will be written to the file.
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
    integer(INT_KIND)                          :: size_fdata_dim_4 !<Number of data points in a non-unstructured dimension (size(fdata,3)).
    integer(INT_KIND)                          :: mynelems         !<Number of data points in the unstructured dimension (size(fdata,1)).
    real,dimension(:),allocatable              :: sbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    integer(INT_KIND)                          :: nelems           !<Total number of unstructured data points (sum(nelems_io) for the root rank of the pelist).
    real,dimension(:),allocatable              :: rbuff            !<Buffer used to gather the data onto the root rank of the pelist.
    real,dimension(:,:,:,:),allocatable        :: cdata            !<Array used to write the data to the file after the gather is performed.
    integer(INT_KIND)                          :: i                !<Loop variable.
    integer(INT_KIND)                          :: j                !<Loop variable.
    integer(INT_KIND)                          :: k                !<Loop variable.
    integer(INT_KIND)                          :: n                !<Loop variable.

   !Start the mpp timer.
   !mpp_write_clock is a module variable.
    call mpp_clock_begin(mpp_write_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_4D:" &
                       //" you must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_4D:" &
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
                       "mpp_io_unstructured_write_compressed_r_4D:" &
                       //" the root rank of the pelist must be allowed" &
                       //" to perform the write.")
    endif
    if (mpp_pe() .ne. pelist(1) .and. mpp_file(funit)%write_on_this_pe) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_write_compressed_r_4D:" &
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
        call mpp_get_UG_compute_domains(io_domain, &
                                       size=nelems_io)
    endif

   !The number of data elements in the non-unstructured dimensions are
   !required to be the same for all ranks.  Perform gathers to check this.
    size_fdata_dim_2 = size(fdata,2)
    size_fdata_dim_3 = size(fdata,3)
    size_fdata_dim_4 = size(fdata,4)

   !Allocate arrays which will be used to gather the data to be written
   !onto the root rank of the pelist.
    mynelems = size(fdata,1)
    allocate(sbuff(mynelems*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4))
    if (mpp_pe() .eq. pelist(1)) then
        nelems = sum(nelems_io)
        allocate(rbuff(nelems*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4))
    else
        allocate(rbuff(1))
    endif

   !Load the data into the sbuff array.  The data is transposed so that the
   !gather may be performed more easily.
    do k = 1,mynelems
        do j = 1,size_fdata_dim_2
            do i = 1,size_fdata_dim_3
                do n = 1,size_fdata_dim_4
                    sbuff((k-1)*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4 &
                          + (j-1)*size_fdata_dim_3*size_fdata_dim_4 &
                          + (i-1)*size_fdata_dim_4 + n) = fdata(k,j,i,n)
                enddo
            enddo
        enddo
    enddo

   !Perform the gather of data onto the root rank (pelist(1)).
    call mpp_gather(sbuff, &
                    size(sbuff), &
                    rbuff, &
                    nelems_io*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4, &
                    pelist)

   !Write out the data to the file.  This is only done by the root rank
   !of the pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(cdata(nelems,size_fdata_dim_2,size_fdata_dim_3,size_fdata_dim_4))
        cdata = fill
        do n = 1,size_fdata_dim_4
            do k = 1,size_fdata_dim_3
                do j = 1,size_fdata_dim_2
                    do i = 1,nelems
                        cdata(i,j,k,n) = rbuff((i-1)*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4 &
                                               + (j-1)*size_fdata_dim_3*size_fdata_dim_4 &
                                               + (k-1)*size_fdata_dim_4 + n)
                    enddo
                enddo
            enddo
        enddo
        field%size(1) = nelems
        call write_record_default(funit, &
                                  field, &
                                  nelems*size_fdata_dim_2*size_fdata_dim_3*size_fdata_dim_4, &
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
end subroutine mpp_io_unstructured_write_r_4D

!------------------------------------------------------------------------------

!----------
