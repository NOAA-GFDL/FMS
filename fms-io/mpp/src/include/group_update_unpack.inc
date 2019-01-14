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

if( group%k_loop_inside ) then
!$OMP parallel do default(none) shared(nunpack,group,nscalar,ptr,nvector,ksize,buffer_start_pos) &
!$OMP                          private(buffer_pos,pos,m,is, ie, js, je,rotation, &
!$OMP                                  ptr_field, ptr_fieldx, ptr_fieldy, n,k )
  do n = nunpack, 1, -1
     buffer_pos = group%unpack_buffer_pos(n) + buffer_start_pos
     pos = buffer_pos 
     is = group%unpack_is(n); ie = group%unpack_ie(n)
     js = group%unpack_js(n); je = group%unpack_je(n)
     if( group%unpack_type(n) == FIELD_S ) then
        do l=1,nscalar  ! loop over number of fields
           ptr_field = group%addrs_s(l)
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    field(i,j,k) = buffer(pos)
                 end do
              end do
           end do
        end do
     else if( group%unpack_type(n) == FIELD_X ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldx = group%addrs_x(l)
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    fieldx(i,j,k) = buffer(pos)
                 end do
              end do
           end do
        end do
     else if( group%unpack_type(n) == FIELD_Y ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldy = group%addrs_y(l)
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    fieldy(i,j,k) = buffer(pos)
                 end do
              end do
           end do
        end do
     endif
  enddo
else
!$OMP parallel do default(none) shared(nunpack,group,nscalar,ptr,nvector,ksize,buffer_start_pos) &
!$OMP                          private(buffer_pos,pos,m,is, ie, js, je,rotation, &
!$OMP                                  ptr_field, ptr_fieldx, ptr_fieldy,n,k)
  do nk = nunpack*ksize, 1, -1
     n = (nk-1)/ksize + 1
     k = mod((nk-1), ksize) + 1
     buffer_pos = group%unpack_buffer_pos(n) + buffer_start_pos
     pos = buffer_pos + (k-1)*group%unpack_size(n)
     is = group%unpack_is(n); ie = group%unpack_ie(n)
     js = group%unpack_js(n); je = group%unpack_je(n)
     if( group%unpack_type(n) == FIELD_S ) then
        do l=1,nscalar  ! loop over number of fields
           ptr_field = group%addrs_s(l)
           do j = js, je
              do i = is, ie
                 pos = pos + 1
                 field(i,j,k) = buffer(pos)
              end do
           end do
        end do
     else if( group%unpack_type(n) == FIELD_X ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldx = group%addrs_x(l)
           do j = js, je
              do i = is, ie
                 pos = pos + 1
                 fieldx(i,j,k) = buffer(pos)
              end do
           end do
        end do
     else if( group%unpack_type(n) == FIELD_Y ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldy = group%addrs_y(l)
           do j = js, je
              do i = is, ie
                 pos = pos + 1
                 fieldy(i,j,k) = buffer(pos)
              end do
           end do
        end do
     endif
  enddo
endif
