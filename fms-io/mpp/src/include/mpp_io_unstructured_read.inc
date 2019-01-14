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
!>Read in one-dimensional data for a field associated with an unstructured
!!mpp domain.
subroutine mpp_io_unstructured_read_r_1D(funit, &
                                         field, &
                                         domain, &
                                         fdata, &
                                         tindex, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/outputs
    integer(INT_KIND),intent(in)                       :: funit     !<A file unit returned by mpp_open.
    type(fieldtype),intent(in)                         :: field     !<A field whose data will be read in from the file.
    type(domainUG),intent(in)                          :: domain    !<An unstructured mpp domain.
    real,dimension(:),intent(inout)                    :: fdata     !<The data that will be read in from the file.
    integer(INT_KIND),intent(in),optional              :: tindex    !<Time level index for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: start     !<Corner indices for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: nread     !<Edge lengths for a NetCDF file.
    integer(INT_KIND),intent(in),optional              :: threading !<Flag telling whether one or multiple ranks will read the file.

   !Local variables
    integer(INT_KIND)                          :: threading_flag !<Flag telling whether one or multiple ranks will read the file.  This defaults to MPP_SINGLE.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: p              !<Loop variable.
    logical(INT_KIND)                          :: compute_chksum !<Flag telling whether or not a check-sum of the read-in data is calculated.
    integer(LONG_KIND)                         :: chk            !<Calculated check-sum for the read in data.

   !Start the mpp timer.
   !mpp_read_clock is a module variable.
    call mpp_clock_begin(mpp_read_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_1D:" &
                       //" you must must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_1D:" &
                       //" the inputted file unit is already in use.")
    endif

   !If the data array has more than zero elements, then read in the data.
    if (size(fdata) .gt. 0) then

       !Initialize the data to zero.
        fdata = 0

       !Get the value for the "threading" flag.
        threading_flag = MPP_SINGLE
        if (present(threading)) then
            threading_flag = threading
        endif

       !Read in the data.
        if (threading_flag .eq. MPP_MULTI) then

           !For the multi-rank case, directly read in the data.
            call read_record(funit, &
                             field, &
                             size(fdata), &
                             fdata, &
                             tindex, &
                             start_in=start, &
                             axsiz_in=nread)
        elseif (threading_flag .eq. MPP_SINGLE) then

           !For the single-rank, first point to the I/O domain associated with
           !the inputted unstructured mpp domain.
            io_domain => null()
            io_domain => mpp_get_UG_io_domain(domain)

           !Get the pelist associated with the I/O domain.
            io_domain_npes = mpp_get_UG_domain_npes(io_domain)
            allocate(pelist(io_domain_npes))
            call mpp_get_UG_domain_pelist(io_domain, &
                                          pelist)
            io_domain => null()

           !Let only the root rank of the pelist read in the data.
            if (mpp_pe() .eq. pelist(1)) then
                call read_record(funit, &
                                 field, &
                                 size(fdata), &
                                 fdata, &
                                 tindex, &
                                 start_in=start, &
                                 axsiz_in=nread)
            endif

           !Send the data from the root rank to the rest of the ranks on the
           !pelist.
            if (mpp_pe() .eq. pelist(1)) then
                do p = 2,io_domain_npes
                    call mpp_send(fdata, &
                                  size(fdata), &
                                  pelist(p), &
                                  tag=COMM_TAG_1)
                enddo
                call mpp_sync_self()
            else
                call mpp_recv(fdata, &
                              size(fdata), &
                              pelist(1), &
                              block=.false., &
                              tag=COMM_TAG_1)
                call mpp_sync_self(check=EVENT_RECV)
            endif
            deallocate(pelist)
        else
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_read_r_1D:" &
                           //" threading should be MPP_SINGLE or MPP_MULTI")
        endif
    endif

   !Decided whether or not to compute a check-sum of the read-in data.  The
   !check-sum is calculated if the inputted field's checksum values are not
   !equal to the default checksum value for a field.
    compute_chksum = .false.
    if (any(field%checksum .ne. default_field%checksum)) then
        compute_chksum = .true.
    endif

   !If necessary, compute a check-sum of the read-in data.
    if (compute_chksum) then
#ifdef use_netCDF
        if (field%type .eq. NF_INT) then
            if (field%fill .eq. MPP_FILL_DOUBLE .or. field%fill .eq. &
                real(MPP_FILL_INT)) then
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=MPP_FILL_INT)
            else
                call mpp_error(NOTE, &
                               "mpp_io_unstructured_read_r_1D:" &
                               //" int field "//trim(field%name) &
                               //" found fill. Icebergs, or code using" &
                               //" defaults can safely ignore." &
                               //" If manually overriding compressed" &
                               //" restart fills, confirm this is what you" &
                               //" want.")
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=field%fill)
            endif
        else
            chk = mpp_chksum(fdata, &
                             mask_val=field%fill)
        endif
#endif
       !Print out the computed check-sum for the field.  This feature is
       !currently turned off.  Uncomment the following lines to turn it
       !back on.
!       if (mpp_pe() .eq. mpp_root_pe()) then
!           write(stdout(),'(A,Z16)') "mpp_read_compressed_2d chksum: " &
!                                     //trim(field%name)//" = ",chk
!           if (mod(chk,field%checksum(1)) .ne. 0) then
!               write(stdout(),'(A,Z16)') "File stored checksum: " &
!                                         //trim(field%name)//" = ", &
!                                         field%checksum(1)
!               call mpp_error(NOTE, &
!                              "mpp_io_unstructured_read_r_1D: " &
!                              //trim(field%name)//" failed!")
!           endif
!       endif
    endif

   !Stop the mpp timer.
    call mpp_clock_end(mpp_read_clock)

    return
end subroutine mpp_io_unstructured_read_r_1D

!------------------------------------------------------------------------------
!>Read in two-dimensional data for a field associated with an unstructured
!!mpp domain.
subroutine mpp_io_unstructured_read_r_2D(funit, &
                                         field, &
                                         domain, &
                                         fdata, &
                                         tindex, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/outputs
    integer(INT_KIND),intent(in)                       :: funit     !<A file unit returned by mpp_open.
    type(fieldtype),intent(in)                         :: field     !<A field whose data will be read in from the file.
    type(domainUG),intent(in)                          :: domain    !<An unstructured mpp domain.
    real,dimension(:,:),intent(inout)                  :: fdata     !<The data that will be read in from the file.
    integer(INT_KIND),intent(in),optional              :: tindex    !<Time level index for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: start     !<Corner indices for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: nread     !<Edge lengths for a NetCDF file.
    integer(INT_KIND),intent(in),optional              :: threading !<Flag telling whether one or multiple ranks will read the file.

   !Local variables
    integer(INT_KIND)                          :: threading_flag !<Flag telling whether one or multiple ranks will read the file.  This defaults to MPP_SINGLE.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: p              !<Loop variable.
    logical(INT_KIND)                          :: compute_chksum !<Flag telling whether or not a check-sum of the read-in data is calculated.
    integer(LONG_KIND)                         :: chk            !<Calculated check-sum for the read in data.

   !Start the mpp timer.
   !mpp_read_clock is a module variable.
    call mpp_clock_begin(mpp_read_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_2D:" &
                       //" you must must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_2D:" &
                       //" the inputted file unit is already in use.")
    endif

   !If the data array has more than zero elements, then read in the data.
    if (size(fdata) .gt. 0) then

       !Initialize the data to zero.
        fdata = 0

       !Get the value for the "threading" flag.
        threading_flag = MPP_SINGLE
        if (present(threading)) then
            threading_flag = threading
        endif

       !Read in the data.
        if (threading_flag .eq. MPP_MULTI) then

           !For the multi-rank case, directly read in the data.
            call read_record(funit, &
                             field, &
                             size(fdata), &
                             fdata, &
                             tindex, &
                             start_in=start, &
                             axsiz_in=nread)
        elseif (threading_flag .eq. MPP_SINGLE) then

           !For the single-rank, first point to the I/O domain associated with
           !the inputted unstructured mpp domain.
            io_domain => null()
            io_domain => mpp_get_UG_io_domain(domain)

           !Get the pelist associated with the I/O domain.
            io_domain_npes = mpp_get_UG_domain_npes(io_domain)
            allocate(pelist(io_domain_npes))
            call mpp_get_UG_domain_pelist(io_domain, &
                                          pelist)
            io_domain => null()

           !Let only the root rank of the pelist read in the data.
            if (mpp_pe() .eq. pelist(1)) then
                call read_record(funit, &
                                 field, &
                                 size(fdata), &
                                 fdata, &
                                 tindex, &
                                 start_in=start, &
                                 axsiz_in=nread)
            endif

           !Send the data from the root rank to the rest of the ranks on the
           !pelist.
            if (mpp_pe() .eq. pelist(1)) then
                do p = 2,io_domain_npes
                    call mpp_send(fdata, &
                                  size(fdata), &
                                  pelist(p), &
                                  tag=COMM_TAG_1)
                enddo
                call mpp_sync_self()
            else
                call mpp_recv(fdata, &
                              size(fdata), &
                              pelist(1), &
                              block=.false., &
                              tag=COMM_TAG_1)
                call mpp_sync_self(check=EVENT_RECV)
            endif
            deallocate(pelist)
        else
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_read_r_2D:" &
                           //" threading should be MPP_SINGLE or MPP_MULTI")
        endif
    endif

   !Decided whether or not to compute a check-sum of the read-in data.  The
   !check-sum is calculated if the inputted field's checksum values are not
   !equal to the default checksum value for a field.
    compute_chksum = .false.
    if (any(field%checksum .ne. default_field%checksum)) then
        compute_chksum = .true.
    endif

   !If necessary, compute a check-sum of the read-in data.
    if (compute_chksum) then
#ifdef use_netCDF
        if (field%type .eq. NF_INT) then
            if (field%fill .eq. MPP_FILL_DOUBLE .or. field%fill .eq. &
                real(MPP_FILL_INT)) then
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=MPP_FILL_INT)
            else
                call mpp_error(NOTE, &
                               "mpp_io_unstructured_read_r_2D:" &
                               //" int field "//trim(field%name) &
                               //" found fill. Icebergs, or code using" &
                               //" defaults can safely ignore." &
                               //" If manually overriding compressed" &
                               //" restart fills, confirm this is what you" &
                               //" want.")
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=field%fill)
            endif
        else
            chk = mpp_chksum(fdata, &
                             mask_val=field%fill)
        endif
#endif
       !Print out the computed check-sum for the field.  This feature is
       !currently turned off.  Uncomment the following lines to turn it
       !back on.
!       if (mpp_pe() .eq. mpp_root_pe()) then
!           write(stdout(),'(A,Z16)') "mpp_read_compressed_2d chksum: " &
!                                     //trim(field%name)//" = ",chk
!           if (mod(chk,field%checksum(1)) .ne. 0) then
!               write(stdout(),'(A,Z16)') "File stored checksum: " &
!                                         //trim(field%name)//" = ", &
!                                         field%checksum(1)
!               call mpp_error(NOTE, &
!                              "mpp_io_unstructured_read_r_2D: " &
!                              //trim(field%name)//" failed!")
!           endif
!       endif
    endif

   !Stop the mpp timer.
    call mpp_clock_end(mpp_read_clock)

    return
end subroutine mpp_io_unstructured_read_r_2D

!------------------------------------------------------------------------------
!>Read in three-dimensional data for a field associated with an unstructured
!!mpp domain.
subroutine mpp_io_unstructured_read_r_3D(funit, &
                                         field, &
                                         domain, &
                                         fdata, &
                                         tindex, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/outputs
    integer(INT_KIND),intent(in)                       :: funit     !<A file unit returned by mpp_open.
    type(fieldtype),intent(in)                         :: field     !<A field whose data will be read in from the file.
    type(domainUG),intent(in)                          :: domain    !<An unstructured mpp domain.
    real,dimension(:,:,:),intent(inout)                :: fdata     !<The data that will be read in from the file.
    integer(INT_KIND),intent(in),optional              :: tindex    !<Time level index for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: start     !<Corner indices for a NetCDF file.
    integer(INT_KIND),dimension(:),intent(in),optional :: nread     !<Edge lengths for a NetCDF file.
    integer(INT_KIND),intent(in),optional              :: threading !<Flag telling whether one or multiple ranks will read the file.

   !Local variables
    integer(INT_KIND)                          :: threading_flag !<Flag telling whether one or multiple ranks will read the file.  This defaults to MPP_SINGLE.
    type(domainUG),pointer                     :: io_domain      !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist         !<A pelist.
    integer(INT_KIND)                          :: p              !<Loop variable.
    logical(INT_KIND)                          :: compute_chksum !<Flag telling whether or not a check-sum of the read-in data is calculated.
    integer(LONG_KIND)                         :: chk            !<Calculated check-sum for the read in data.

   !Start the mpp timer.
   !mpp_read_clock is a module variable.
    call mpp_clock_begin(mpp_read_clock)

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_3D:" &
                       //" you must must first call mpp_io_init.")
    endif

   !Make sure that another NetCDF file is not currently using the inputted
   !file unit.
    if (.not. mpp_file(funit)%valid) then
        call mpp_error(FATAL, &
                       "mpp_io_unstructured_read_r_3D:" &
                       //" the inputted file unit is already in use.")
    endif

   !If the data array has more than zero elements, then read in the data.
    if (size(fdata) .gt. 0) then

       !Initialize the data to zero.
        fdata = 0

       !Get the value for the "threading" flag.
        threading_flag = MPP_SINGLE
        if (present(threading)) then
            threading_flag = threading
        endif

       !Read in the data.
        if (threading_flag .eq. MPP_MULTI) then

           !For the multi-rank case, directly read in the data.
            call read_record(funit, &
                             field, &
                             size(fdata), &
                             fdata, &
                             tindex, &
                             start_in=start, &
                             axsiz_in=nread)
        elseif (threading_flag .eq. MPP_SINGLE) then

           !For the single-rank, first point to the I/O domain associated with
           !the inputted unstructured mpp domain.
            io_domain => null()
            io_domain => mpp_get_UG_io_domain(domain)

           !Get the pelist associated with the I/O domain.
            io_domain_npes = mpp_get_UG_domain_npes(io_domain)
            allocate(pelist(io_domain_npes))
            call mpp_get_UG_domain_pelist(io_domain, &
                                          pelist)
            io_domain => null()

           !Let only the root rank of the pelist read in the data.
            if (mpp_pe() .eq. pelist(1)) then
                call read_record(funit, &
                                 field, &
                                 size(fdata), &
                                 fdata, &
                                 tindex, &
                                 start_in=start, &
                                 axsiz_in=nread)
            endif

           !Send the data from the root rank to the rest of the ranks on the
           !pelist.
            if (mpp_pe() .eq. pelist(1)) then
                do p = 2,io_domain_npes
                    call mpp_send(fdata, &
                                  size(fdata), &
                                  pelist(p), &
                                  tag=COMM_TAG_1)
                enddo
                call mpp_sync_self()
            else
                call mpp_recv(fdata, &
                              size(fdata), &
                              pelist(1), &
                              block=.false., &
                              tag=COMM_TAG_1)
                call mpp_sync_self(check=EVENT_RECV)
            endif
            deallocate(pelist)
        else
            call mpp_error(FATAL, &
                           "mpp_io_unstructured_read_r_3D:" &
                           //" threading should be MPP_SINGLE or MPP_MULTI")
        endif
    endif

   !Decided whether or not to compute a check-sum of the read-in data.  The
   !check-sum is calculated if the inputted field's checksum values are not
   !equal to the default checksum value for a field.
    compute_chksum = .false.
    if (any(field%checksum .ne. default_field%checksum)) then
        compute_chksum = .true.
    endif

   !If necessary, compute a check-sum of the read-in data.
    if (compute_chksum) then
#ifdef use_netCDF
        if (field%type .eq. NF_INT) then
            if (field%fill .eq. MPP_FILL_DOUBLE .or. field%fill .eq. &
                real(MPP_FILL_INT)) then
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=MPP_FILL_INT)
            else
                call mpp_error(NOTE, &
                               "mpp_io_unstructured_read_r_3D:" &
                               //" int field "//trim(field%name) &
                               //" found fill. Icebergs, or code using" &
                               //" defaults can safely ignore." &
                               //" If manually overriding compressed" &
                               //" restart fills, confirm this is what you" &
                               //" want.")
                chk = mpp_chksum(ceiling(fdata), &
                                 mask_val=field%fill)
            endif
        else
            chk = mpp_chksum(fdata, &
                             mask_val=field%fill)
        endif
#endif
       !Print out the computed check-sum for the field.  This feature is
       !currently turned off.  Uncomment the following lines to turn it
       !back on.
!       if (mpp_pe() .eq. mpp_root_pe()) then
!           write(stdout(),'(A,Z16)') "mpp_read_compressed_2d chksum: " &
!                                     //trim(field%name)//" = ",chk
!           if (mod(chk,field%checksum(1)) .ne. 0) then
!               write(stdout(),'(A,Z16)') "File stored checksum: " &
!                                         //trim(field%name)//" = ", &
!                                         field%checksum(1)
!               call mpp_error(NOTE, &
!                              "mpp_io_unstructured_read_r_3D: " &
!                              //trim(field%name)//" failed!")
!           endif
!       endif
    endif

   !Stop the mpp timer.
    call mpp_clock_end(mpp_read_clock)

    return
end subroutine mpp_io_unstructured_read_r_3D

!------------------------------------------------------------------------------

!----------
