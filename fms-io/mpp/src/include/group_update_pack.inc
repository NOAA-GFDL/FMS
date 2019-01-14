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
!$OMP parallel do default(none) shared(npack,group,ptr,nvector,ksize,buffer_start_pos) &
!$OMP                          private(buffer_pos,pos,m,is,ie,js,je,rotation, &
!$OMP                                  ptr_field, ptr_fieldx, ptr_fieldy,n,k)
  do n = 1, npack
     buffer_pos = group%pack_buffer_pos(n) + buffer_start_pos
     pos  = buffer_pos
     is = group%pack_is(n); ie = group%pack_ie(n)
     js = group%pack_js(n); je = group%pack_je(n)
     rotation = group%pack_rotation(n)
     if( group%pack_type(n) == FIELD_S ) then
        select case( rotation )
        case(ZERO)
           do l=1, group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do k = 1, ksize
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              enddo
           enddo
        case( MINUS_NINETY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do k = 1, ksize
                 do i = is, ie
                    do j = je, js, -1
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           end do
        case( NINETY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do k = 1, ksize
                 do i = ie, is, -1
                    do j = js, je
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           end do
        case( ONE_HUNDRED_EIGHTY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do k = 1, ksize
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) = field(i,j,k)
                    end do
                 end do
              end do
           end do
        end select
     else if( group%pack_type(n) == FIELD_X  ) then
        select case( rotation )
        case(ZERO)
           do l=1, nvector  ! loop over number of fields
              ptr_fieldx = group%addrs_x(l)
              do k = 1, ksize
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = fieldx(i,j,k)
                    end do
                 end do
              end do
           end do
        case( MINUS_NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1, ksize
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1, ksize
                    do i = is, ie
                       do j = je, js, -1
                          pos = pos + 1
                          buffer(pos) = -fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           end if
        case( NINETY )
           do l=1, nvector  ! loop over number of fields
              ptr_fieldy = group%addrs_y(l)
              do k = 1, ksize
                 do i = ie, is, -1
                    do j = js, je
                       pos = pos + 1
                       buffer(pos) = fieldy(i,j,k)
                    end do
                 end do
              end do
           end do
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do k = 1, ksize
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) =  fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do k = 1, ksize
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) =  -fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           end if
        end select ! select case( rotation(n) )
     else if( group%pack_type(n) == FIELD_Y ) then
        select case( rotation )
        case(ZERO)
           do l=1, nvector  ! loop over number of fields
              ptr_fieldy = group%addrs_y(l)
              do k = 1, ksize
                 do j = js, je
                    do i = is, ie
                       pos = pos + 1
                       buffer(pos) = fieldy(i,j,k)
                    end do
                 end do
              end do
           end do
        case( MINUS_NINETY )
           do l=1,nvector  ! loop over number of fields
              ptr_fieldx = group%addrs_x(l)
              do k = 1, ksize
                 do i = is, ie
                    do j = je, js, -1
                       pos = pos + 1
                       buffer(pos) = fieldx(i,j,k)
                    end do
                 end do
              end do
           end do
        case( NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1, nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do k = 1, ksize
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do k = 1, ksize
                    do i = ie, is, -1
                       do j = js, je
                          pos = pos + 1
                          buffer(pos) = -fieldx(i,j,k)
                       end do
                    end do
                 end do
              end do
           end if
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1, ksize
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) =  fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do k = 1, ksize
                    do j = je, js, -1
                       do i = ie, is, -1
                          pos = pos + 1
                          buffer(pos) =  -fieldy(i,j,k)
                       end do
                    end do
                 end do
              end do
           end if
        end select ! select case( rotation(n) )
     endif
  enddo
else
!$OMP parallel do default(none) shared(npack,group,ptr,nvector,ksize,buffer_start_pos) &
!$OMP                          private(buffer_pos,pos,m,is,ie,js,je,rotation, &
!$OMP                                  ptr_field, ptr_fieldx, ptr_fieldy,n,k)
  do nk = 1, npack*ksize
     n = (nk-1)/ksize + 1
     k = mod((nk-1), ksize) + 1
     buffer_pos = group%pack_buffer_pos(n) + buffer_start_pos
     pos  = buffer_pos + (k-1)*group%pack_size(n)
     is = group%pack_is(n); ie = group%pack_ie(n)
     js = group%pack_js(n); je = group%pack_je(n)
     rotation = group%pack_rotation(n)
     if( group%pack_type(n) == FIELD_S ) then
        select case( rotation )
        case(ZERO)
           do l=1, group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    buffer(pos) = field(i,j,k)
                 end do
              end do
           enddo
        case( MINUS_NINETY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do i = is, ie
                 do j = je, js, -1
                    pos = pos + 1
                    buffer(pos) = field(i,j,k)
                 end do
              end do
           end do
        case( NINETY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do i = ie, is, -1
                 do j = js, je
                    pos = pos + 1
                    buffer(pos) = field(i,j,k)
                 end do
              end do
           end do
        case( ONE_HUNDRED_EIGHTY )
           do l=1,group%nscalar  ! loop over number of fields
              ptr_field = group%addrs_s(l)
              do j = je, js, -1
                 do i = ie, is, -1
                    pos = pos + 1
                    buffer(pos) = field(i,j,k)
                 end do
              end do
           end do
        end select
     else if( group%pack_type(n) == FIELD_X  ) then
        select case( rotation )
        case(ZERO)
           do l=1, nvector  ! loop over number of fields
              ptr_fieldx = group%addrs_x(l)
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    buffer(pos) = fieldx(i,j,k)
                 end do
              end do
           end do
        case( MINUS_NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do i = is, ie
                    do j = je, js, -1
                       pos = pos + 1
                       buffer(pos) = fieldy(i,j,k)
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do i = is, ie
                    do j = je, js, -1
                       pos = pos + 1
                       buffer(pos) = -fieldy(i,j,k)
                    end do
                 end do
              end do
           end if
        case( NINETY )
           do l=1, nvector  ! loop over number of fields
              ptr_fieldy = group%addrs_y(l)
              do i = ie, is, -1
                 do j = js, je
                    pos = pos + 1
                    buffer(pos) = fieldy(i,j,k)
                 end do
              end do
           end do
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) =  fieldx(i,j,k)
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) =  -fieldx(i,j,k)
                    end do
                 end do
              end do
           end if
        end select ! select case( rotation(n) )
     else if( group%pack_type(n) == FIELD_Y ) then
        select case( rotation )
        case(ZERO)
           do l=1, nvector  ! loop over number of fields
              ptr_fieldy = group%addrs_y(l)
              do j = js, je
                 do i = is, ie
                    pos = pos + 1
                    buffer(pos) = fieldy(i,j,k)
                 end do
              end do
           end do
        case( MINUS_NINETY )
           do l=1,nvector  ! loop over number of fields
              ptr_fieldx = group%addrs_x(l)
              do i = is, ie
                 do j = je, js, -1
                    pos = pos + 1
                    buffer(pos) = fieldx(i,j,k)
                 end do
              end do
           end do
        case( NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1, nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do i = ie, is, -1
                    do j = js, je
                       pos = pos + 1
                       buffer(pos) = fieldx(i,j,k)
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldx = group%addrs_x(l)
                 do i = ie, is, -1
                    do j = js, je
                       pos = pos + 1
                       buffer(pos) = -fieldx(i,j,k)
                    end do
                 end do
              end do
           end if
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) =  fieldy(i,j,k)
                    end do
                 end do
              end do
           else
              do l=1,nvector  ! loop over number of fields
                 ptr_fieldy = group%addrs_y(l)
                 do j = je, js, -1
                    do i = ie, is, -1
                       pos = pos + 1
                       buffer(pos) =  -fieldy(i,j,k)
                    end do
                 end do
              end do
           end if
        end select ! select case( rotation(n) )
     endif
  enddo
endif
