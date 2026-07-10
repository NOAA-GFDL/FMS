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

subroutine GROUP_UPDATE_PACK_ (group, buffer, buffer_start_pos)
  type(mpp_group_update_type), intent(in) :: group
  MPP_TYPE_, intent(out) :: buffer(:)
  integer, intent(in) :: buffer_start_pos
  integer, parameter :: SWAP_IJ=1, REVERSE_I=2, REVERSE_J=4, MINUS=8 !< Flags which affect packing order
  integer :: n, buffer_pos, pos, ltop, ni, nj, nk, nijk, l

  do n=1,group%npack
     buffer_pos = group%pack_buffer_pos(n) + buffer_start_pos
     ni = group%pack_ie(n) - group%pack_is(n) + 1
     nj = group%pack_je(n) - group%pack_js(n) + 1
     if( group%pack_type(n) == FIELD_S ) then
       ltop = group%nscalar
       nk = group%ksize_s
     else
       ltop = group%nvector
       nk = group%ksize_v
     endif
     nijk = ni*nj*nk
     if( group%pack_type(n) == FIELD_S ) then
        select case( group%pack_rotation(n) )
        case(ZERO)
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_S, n, l, 0)
           enddo
        case( MINUS_NINETY )
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_S, n, l, SWAP_IJ + REVERSE_J)
           end do
        case( NINETY )
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_S, n, l, SWAP_IJ + REVERSE_I)
           end do
        case( ONE_HUNDRED_EIGHTY )
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_S, n, l, REVERSE_I + REVERSE_J)
           end do
        end select
     else if( group%pack_type(n) == FIELD_X  ) then
        select case( group%pack_rotation(n) )
        case(ZERO)
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, 0)
           end do
        case( MINUS_NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, SWAP_IJ + REVERSE_J)
              end do
           else
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, SWAP_IJ + REVERSE_J + MINUS)
              end do
           end if
        case( NINETY )
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, SWAP_IJ + REVERSE_I)
           end do
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, REVERSE_I + REVERSE_J)
              end do
           else
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, REVERSE_I + REVERSE_J + MINUS)
              end do
           end if
        end select
     else if( group%pack_type(n) == FIELD_Y ) then
        select case( group%pack_rotation(n) )
        case(ZERO)
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, 0)
           end do
        case( MINUS_NINETY )
           do l=1,ltop
              pos = buffer_pos + (l-1)*nijk + 1
              call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, SWAP_IJ + REVERSE_J)
           end do
        case( NINETY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, SWAP_IJ + REVERSE_I)
              end do
           else
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_X, n, l, SWAP_IJ + REVERSE_I + MINUS)
              end do
           end if
        case( ONE_HUNDRED_EIGHTY )
           if( BTEST(group%flags_v,SCALAR_BIT) ) then
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, REVERSE_I + REVERSE_J)
              end do
           else
              do l=1,ltop
                 pos = buffer_pos + (l-1)*nijk + 1
                 call pack_inner_loop (buffer(pos:), group, FIELD_Y, n, l, REVERSE_I + REVERSE_J + MINUS)
              end do
           end if
        end select
     endif
  enddo

contains

  subroutine pack_inner_loop (buf, group, field_type, ni, l, packing_flags)
    MPP_TYPE_, intent(out) :: buf(:)
    type(mpp_group_update_type), intent(in), target :: group
    integer, intent(in) :: field_type !< FIELD_S, FIELD_X, or FIELD_Y
    integer, intent(in) :: ni
    integer, intent(in) :: l !< Index of the field within group
    integer, intent(in) :: packing_flags !< Any combination of: SWAP_IJ, REVERSE_I, REVERSE_J, MINUS

    integer, allocatable, dimension(:) :: lb, ub, n, stride
    MPP_TYPE_, pointer :: field2d(:,:)
    MPP_TYPE_, pointer :: field3d(:,:,:)
    MPP_TYPE_, pointer :: field4d(:,:,:,:)

    integer :: ix, iy, iz1, iz2 !< Indices of the domain-decomposed dimensions
    integer :: nd !< Number of dimensions
    integer :: sgn !< Factor to multiply the filled values by
    integer :: i1, i2, i3, i4
    integer :: ij1, ij2
    integer :: is0, js0

    select case (field_type)
      case (FIELD_S)
        ix = group%ix_s
        iy = group%iy_s
        is0 = group%is_s
        js0 = group%js_s

        nd = size(group%shape_s)
        allocate (lb(nd), ub(nd), n(nd), stride(nd))

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
        is0 = group%is_x
        js0 = group%js_x

        nd = size(group%shape_x)
        allocate (lb(nd), ub(nd), n(nd), stride(nd))

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
        is0 = group%is_y
        js0 = group%js_y

        nd = size(group%shape_y)
        allocate (lb(nd), ub(nd), n(nd), stride(nd))

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

    lb(ix) = group%pack_is(ni) - is0 + 1
    lb(iy) = group%pack_js(ni) - js0 + 1

    ub(ix) = group%pack_ie(ni) - is0 + 1
    ub(iy) = group%pack_je(ni) - js0 + 1

    n = ub - lb + 1

    stride = 1

    if (iand(packing_flags, REVERSE_I).ne.0) then
      lb(ix) = group%pack_ie(ni) - is0 + 1
      ub(ix) = group%pack_is(ni) - is0 + 1
      stride(ix) = -1
    endif

    if (iand(packing_flags, REVERSE_J).ne.0) then
      lb(iy) = group%pack_je(ni) - js0 + 1
      ub(iy) = group%pack_js(ni) - js0 + 1
      stride(iy) = -1
    endif

    if (iand(packing_flags, SWAP_IJ).eq.0) then
      ij1 = ix
      ij2 = iy
    else
      ij1 = iy
      ij2 = ix
    endif

    sgn = 1
    if (iand(packing_flags, MINUS).ne.0) sgn = -1

    select case (nd)
      case (2)
        !$omp parallel loop
        do i1=lb(1),ub(1),stride(1)
          do i2=lb(2),ub(2),stride(2)
            block
              integer :: i(2), idx

              i = [i1, i2]
              idx = n(ij1)*stride(ij2)*(i(ij2)-lb(ij2)) + &
                    stride(ij1)*(i(ij1)-lb(ij1)) + 1
              buf(idx) = sgn * field2d(i1,i2)
            end block
          end do
        end do
      case (3)
        iz1 = 1 + 2 + 3 - ix - iy
        !$omp parallel loop
        do i1=lb(1),ub(1),stride(1)
          do i2=lb(2),ub(2),stride(2)
            do i3=lb(3),ub(3),stride(3)
              block
                integer :: i(3), idx

                i = [i1, i2, i3]
                idx = n(ij1)*( n(ij2)*(i(iz1)-1) + stride(ij2)*(i(ij2)-lb(ij2))) + &
                      stride(ij1)*(i(ij1)-lb(ij1)) + 1
                buf(idx) = sgn * field3d(i1,i2,i3)
              end block
            end do
          end do
        end do
      case (4)
        iz1 = trailz(sum(2**[1,2,3,4]) - 2**ix - 2**iy)
        iz2 = 1 + 2 + 3 + 4 - ix - iy - iz1
        !$omp parallel loop
        do i1=lb(1),ub(1),stride(1)
          do i2=lb(2),ub(2),stride(2)
            do i3=lb(3),ub(3),stride(3)
              do i4=lb(4),ub(4),stride(4)
                block
                  integer :: i(4), idx

                  i = [i1, i2, i3, i4]
                  idx = n(ij1)*( n(ij2)*( n(iz1)*(i(iz2)-1) + (i(iz1)-1)) + stride(ij2)*(i(ij2)-lb(ij2))) + &
                        stride(ij1)*(i(ij1)-lb(ij1)) + 1
                  buf(idx) = sgn * field4d(i1,i2,i3,i4)
                end block
              end do
            end do
          end do
        end do
    end select
  end subroutine pack_inner_loop
end subroutine GROUP_UPDATE_PACK_
