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
! This subroutine pass data from unstructured domain2d to domain.
! First only implement for data at grid cell center.
SUBROUTINE mpp_pass_SG_to_UG_2D_(UG_domain, field_SG, field_UG)
  type(domainUG), intent(in) :: UG_domain
  MPP_TYPE_,          intent(inout) :: field_UG(:)
  MPP_TYPE_,             intent(in) :: field_SG(:,:)
  MPP_TYPE_ :: field3D_SG(size(field_SG,1),size(field_SG,2),1)
  MPP_TYPE_ :: field2D_UG(size(field_UG(:)), 1)
  pointer(ptr_SG, field3D_SG)
  pointer(ptr_UG, field2D_UG)
  ptr_SG = LOC(field_SG)
  ptr_UG = LOC(field_UG)
  call mpp_pass_SG_to_UG(UG_domain, field3D_SG, field2D_UG)

end SUBROUTINE mpp_pass_SG_to_UG_2D_


SUBROUTINE mpp_pass_SG_to_UG_3D_(UG_domain, field_SG, field_UG)
  type(domainUG), intent(in) :: UG_domain
  MPP_TYPE_,          intent(inout) :: field_UG(:,:)
  MPP_TYPE_,             intent(in) :: field_SG(:,:,:)
  MPP_TYPE_        :: buffer(mpp_domains_stack_size)
  character(len=8) :: text
  integer :: i, j, k, l, m, ke, from_pe, to_pe
  integer :: buffer_pos, pos, msgsize, ioff, joff

  !--- check if data is on data or computing domain
  if(size(field_SG,1) .EQ. UG_domain%SG_domain%x(1)%compute%size .AND.  &
     size(field_SG,2) .EQ. UG_domain%SG_domain%y(1)%compute%size) then
     ioff = 1 - UG_domain%SG_domain%x(1)%compute%begin
     joff = 1 - UG_domain%SG_domain%y(1)%compute%begin
  else if(size(field_SG,1) .EQ. UG_domain%SG_domain%x(1)%data%size .AND.  &
          size(field_SG,2) .EQ. UG_domain%SG_domain%y(1)%data%size) then
     ioff = 1 - UG_domain%SG_domain%x(1)%data%begin
     joff = 1 - UG_domain%SG_domain%y(1)%data%begin
  else
     call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D_: field_SG must match either compute domain or data domain.' )
  endif

  ke = size(field_SG,3)
  buffer_pos = 0
  !---pre-post receive
  do m = 1, UG_domain%SG2UG%nrecv
     msgsize = UG_domain%SG2UG%recv(m)%count*ke
     if( msgsize.GT.0 )then
        from_pe = UG_domain%SG2UG%recv(m)%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do 

  !---pack and send data
  do m = 1, UG_domain%SG2UG%nsend
     pos = buffer_pos
     msgsize = UG_domain%SG2UG%send(m)%count * ke
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_SG_to_UG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if

        do k = 1, ke
           do l = 1, UG_domain%SG2UG%send(m)%count
              i = UG_domain%SG2UG%send(m)%i(l)+ioff
              j = UG_domain%SG2UG%send(m)%j(l)+joff
              pos = pos+1
              buffer(pos) = field_SG(i,j,k)
           end do
        end do
        to_pe = UG_domain%SG2UG%send(m)%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  call mpp_sync_self(check=EVENT_RECV)
  !--- unpack the buffer
  buffer_pos = 0
  do m = 1, UG_domain%SG2UG%nrecv
     pos = buffer_pos
     do k = 1, ke
        do l = 1, UG_domain%SG2UG%recv(m)%count
           pos = pos+1
           i = UG_domain%SG2UG%recv(m)%i(l)
           field_UG(i,k) = buffer(pos)
        enddo
     enddo
     buffer_pos = pos
  enddo

  call mpp_sync_self( )

end SUBROUTINE mpp_pass_SG_to_UG_3D_


! This subroutine pass data from unstructured domain2d to domain.
! First only implement for data at grid cell center.
SUBROUTINE mpp_pass_UG_to_SG_2D_(UG_domain, field_UG, field_SG)
  type(domainUG), intent(in) :: UG_domain
  MPP_TYPE_,             intent(in) :: field_UG(:)
  MPP_TYPE_,          intent(inout) :: field_SG(:,:)
  MPP_TYPE_ :: field3D_SG(size(field_SG,1),size(field_SG,2),1)
  MPP_TYPE_ :: field2D_UG(size(field_UG(:)), 1)
  pointer(ptr_SG, field3D_SG)
  pointer(ptr_UG, field2D_UG)
  ptr_SG = LOC(field_SG)
  ptr_UG = LOC(field_UG)

  call mpp_pass_UG_to_SG(UG_domain, field2D_UG, field3D_SG)

end SUBROUTINE mpp_pass_UG_to_SG_2D_


SUBROUTINE mpp_pass_UG_to_SG_3D_(UG_domain, field_UG, field_SG)
  type(domainUG), intent(in) :: UG_domain
  MPP_TYPE_,             intent(in) :: field_UG(:,:)
  MPP_TYPE_,          intent(inout) :: field_SG(:,:,:)
  MPP_TYPE_        :: buffer(mpp_domains_stack_size)
  character(len=8) :: text
  integer :: i, j, k, l, m, ke, from_pe, to_pe
  integer :: buffer_pos, pos, msgsize, ioff, joff

  !--- check if data is on data or computing domain
  if(size(field_SG,1) .EQ. UG_domain%SG_domain%x(1)%compute%size .AND.  &
     size(field_SG,2) .EQ. UG_domain%SG_domain%y(1)%compute%size) then
     ioff = 1 - UG_domain%SG_domain%x(1)%compute%begin
     joff = 1 - UG_domain%SG_domain%y(1)%compute%begin
  else if(size(field_SG,1) .EQ. UG_domain%SG_domain%x(1)%data%size .AND.  &
          size(field_SG,2) .EQ. UG_domain%SG_domain%y(1)%data%size) then
     ioff = 1 - UG_domain%SG_domain%x(1)%data%begin
     joff = 1 - UG_domain%SG_domain%y(1)%data%begin
  else
     call mpp_error( FATAL, 'mpp_pass_UG_to_SG_3D_: field_SG must match either compute domain or data domain.' )
  endif

  ke = size(field_SG,3)
  buffer_pos = 0
  !---pre-post receive
  do m = 1, UG_domain%UG2SG%nrecv
     msgsize = UG_domain%UG2SG%recv(m)%count * ke
     if( msgsize.GT.0 )then
        from_pe = UG_domain%UG2SG%recv(m)%pe
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_UG_to_SG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if
        call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do 

  !---pack and send data
  do m = 1, UG_domain%UG2SG%nsend
     pos = buffer_pos
     msgsize = UG_domain%UG2SG%send(m)%count * ke
     if( msgsize.GT.0 )then
        mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, (buffer_pos+msgsize) )
        if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
           write( text,'(i8)' )mpp_domains_stack_hwm
           call mpp_error( FATAL, 'mpp_pass_UG_to_SG_3D: mpp_domains_stack overflow, '// &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
        end if

        do k = 1, ke
           do l = 1, UG_domain%UG2SG%send(m)%count
              pos = pos+1
              i = UG_domain%UG2SG%send(m)%i(l)
              buffer(pos) = field_UG(i,k)
           end do
        end do
        to_pe = UG_domain%UG2SG%send(m)%pe
        call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=to_pe, tag=COMM_TAG_1 )
        buffer_pos = buffer_pos + msgsize
     end if
  end do

  call mpp_sync_self(check=EVENT_RECV)
  !--- unpack the buffer
  buffer_pos = 0
  do m = 1, UG_domain%UG2SG%nrecv
     pos = buffer_pos
     do k = 1, ke
        do l = 1, UG_domain%UG2SG%recv(m)%count
           pos = pos+1
           i = UG_domain%UG2SG%recv(m)%i(l)+ioff
           j = UG_domain%UG2SG%recv(m)%j(l)+joff     
           field_SG(i,j,k) = buffer(pos)
        enddo
     enddo
     buffer_pos = pos
  enddo

  call mpp_sync_self( )

end SUBROUTINE mpp_pass_UG_to_SG_3D_



