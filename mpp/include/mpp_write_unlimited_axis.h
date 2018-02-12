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
    subroutine MPP_WRITE_UNLIMITED_AXIS_1D_(unit,field,domain,data,nelems_io)
      integer,           intent(in)           :: unit
      type(fieldtype),   intent(inout)        :: field
      type(domain2D),    intent(inout)        :: domain
      MPP_TYPE_,         intent(inout)        :: data(:)
      integer,           intent(in)           :: nelems_io(:)  ! number of compressed elements from each
                                                               ! member of the io_domain. It MUST have the
                                                               ! same order as the io_domain pelist.
      integer, allocatable :: pelist(:)
      integer :: i,j,nelems,npes
      type(domain2d), pointer :: io_domain=>NULL()
      MPP_TYPE_, allocatable, dimension(:) :: rbuff

      call mpp_clock_begin(mpp_write_clock)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_WRITE_UNLIMITED_AXIS_1D_: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%valid )call mpp_error( FATAL, 'MPP_WRITE_UNLIMITED_AXIS_1D_: invalid unit number.' )

      io_domain=>mpp_get_io_domain(domain) 
      if (.not. ASSOCIATED(io_domain)) call mpp_error( FATAL, 'MPP_WRITE_UNLIMITED_AXIS_1D_: io_domain must be defined.' )
      npes = mpp_get_domain_npes(io_domain)
      allocate(pelist(npes))
      call mpp_get_pelist(io_domain,pelist)

      nelems = sum(nelems_io(:))

      if(mpp_file(unit)%write_on_this_pe) allocate(rbuff(nelems))

   !  Note that the gatherV implied here is asymmetric; only root needs to know the vector of recv size
      call mpp_gather(data,size(data),rbuff,nelems_io(:),pelist)

      if(mpp_file(unit)%write_on_this_pe) then
         field%size(1) = nelems  ! Correct the field size now that we have all the data
         call write_record(unit, field, nelems, rbuff)
         deallocate(rbuff)
      endif
      deallocate(pelist)

      call mpp_clock_end(mpp_write_clock)
      return
    end subroutine MPP_WRITE_UNLIMITED_AXIS_1D_
