!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine GROUP_UPDATE_UNPACK_ (group, buffer, buffer_start_pos)
  type(mpp_group_update_type), intent(in) :: group
  MPP_TYPE_, intent(in) :: buffer(:)
  integer, intent(in) :: buffer_start_pos
  integer :: n, buffer_pos, pos, ltop, ni, nj, nk, nijk, l

  do n=group%nunpack,1,-1
     buffer_pos = group%unpack_buffer_pos(n) + buffer_start_pos
     ni = group%unpack_ie(n) - group%unpack_is(n) + 1
     nj = group%unpack_je(n) - group%unpack_js(n) + 1
     if( group%unpack_type(n) == FIELD_S ) then
       ltop = group%nscalar
       nk = group%ksize_s
     else
       ltop = group%nvector
       nk = group%ksize_v
     endif
     nijk = ni*nj*nk

     do l=1,ltop
        pos = buffer_pos + (l-1)*nijk + 1
        call unpack_inner_loop (buffer(pos:), group, group%unpack_type(n), n, l)
     end do
  enddo

contains

  subroutine unpack_inner_loop (buf, group, field_type, ni, l)
    MPP_TYPE_, intent(in) :: buf(:)
    type(mpp_group_update_type), intent(in), target :: group
    integer, intent(in) :: field_type !< FIELD_S, FIELD_X, or FIELD_Y
    integer, intent(in) :: ni
    integer, intent(in) :: l !< Index of the field within group

    integer, allocatable, dimension(:) :: lb, ub, n
    MPP_TYPE_, pointer :: field2d(:,:)
    MPP_TYPE_, pointer :: field3d(:,:,:)
    MPP_TYPE_, pointer :: field4d(:,:,:,:)

    integer :: ix, iy, iz1, iz2 !< Indices of the domain-decomposed dimensions
    integer :: nd !< Number of dimensions
    integer :: idx !< Index into the 1D buffer
    integer :: i1, i2, i3, i4, i(4)

    select case (field_type)
      case (FIELD_S)
        ix = group%ix_s
        iy = group%iy_s

        nd = size(group%shape_s)
        allocate (lb(nd), ub(nd), n(nd))

        lb = 1
        ub = group%shape_s

        select case(nd)
          case (2)
            call c_f_pointer(group%addrs_s(l), field2d, group%shape_s)
          case (3)
            call c_f_pointer(group%addrs_s(l), field3d, group%shape_s)
          case (4)
            call c_f_pointer(group%addrs_s(l), field4d, group%shape_s)
        end select
      case (FIELD_X)
        ix = group%ix_v
        iy = group%iy_v

        nd = size(group%shape_x)
        allocate (lb(nd), ub(nd), n(nd))

        lb = 1
        ub = group%shape_x

        select case(nd)
          case (2)
            call c_f_pointer(group%addrs_x(l), field2d, group%shape_x)
          case (3)
            call c_f_pointer(group%addrs_x(l), field3d, group%shape_x)
          case (4)
            call c_f_pointer(group%addrs_x(l), field4d, group%shape_x)
        end select
      case (FIELD_Y)
        ix = group%ix_v
        iy = group%iy_v

        nd = size(group%shape_y)
        allocate (lb(nd), ub(nd), n(nd))

        lb = 1
        ub = group%shape_y

        select case(nd)
          case (2)
            call c_f_pointer(group%addrs_y(l), field2d, group%shape_y)
          case (3)
            call c_f_pointer(group%addrs_y(l), field3d, group%shape_y)
          case (4)
            call c_f_pointer(group%addrs_y(l), field4d, group%shape_y)
        end select
    end select

    lb(ix) = group%unpack_is(ni)
    lb(iy) = group%unpack_js(ni)

    ub(ix) = group%unpack_ie(ni)
    ub(iy) = group%unpack_je(ni)

    n = ub - lb + 1

    select case (nd)
      case (2)
        do concurrent (i1=lb(1):ub(1), i2=lb(2):ub(2)) local(i, idx)
            i(1:2) = [i1, i2]
            idx = n(ix)*(i(iy)-lb(iy)) + (i(ix)-lb(ix)) + 1
            field2d(i1,i2) = buffer(idx)
        end do
      case (3)
        iz1 = 1 + 2 + 3 - ix - iy
        do concurrent (i1=lb(1):ub(1), i2=lb(2):ub(2), i3=lb(3):ub(3)) local(i, idx)
            i(1:3) = [i1, i2, i3]
            idx = n(ix)*( n(iy)*(i(iz1)-1) + (i(iy)-lb(iy))) + (i(ix)-lb(ix)) + 1
            field3d(i1,i2,i3) = buffer(idx)
        end do
      case (4)
        iz1 = trailz(sum(2**[1,2,3,4]) - 2**ix - 2**iy)
        iz2 = 1 + 2 + 3 + 4 - ix - iy - iz1
        do concurrent (i1=lb(1):ub(1), i2=lb(2):ub(2), i3=lb(3):ub(3), i4=lb(4):ub(4)) local(i, idx)
            i(1:4) = [i1, i2, i3, i4]
            idx = n(ix)*( n(iy)*( n(iz1)*(i(iz2)-1) + (i(iz1)-1)) + (i(iy)-lb(iy))) + (i(ix)-lb(ix)) + 1
            field4d(i1,i2,i3,i4) = buffer(idx)
        end do
    end select
  end subroutine unpack_inner_loop

end subroutine GROUP_UPDATE_UNPACK_
