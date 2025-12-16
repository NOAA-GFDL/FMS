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
!$OMP parallel do shared(nunpack,group,nscalar,ptr,nvector,ksize,buffer_start_pos) &
!$OMP            private(buffer_pos,pos,m,is, ie, js, je,rotation, &
!$OMP                    ptr_field, ptr_fieldx, ptr_fieldy, n,k,ni,nj,idx) &
!$OMP                if (.not.use_device_ptr)
  do n = nunpack, 1, -1
     buffer_pos = group%unpack_buffer_pos(n) + buffer_start_pos
     pos = buffer_pos
     is = group%unpack_is(n); ie = group%unpack_ie(n); ni = ie-is+1
     js = group%unpack_js(n); je = group%unpack_je(n); nj = je-js+1
     if( group%unpack_type(n) == FIELD_S ) then
        do l=1,nscalar  ! loop over number of fields
           ptr_field = group%addrs_s(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(3) if(use_device_ptr) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+ksize*nj*ni)) &
           !$omp   map(from: field(is:ie,js:je,1:ksize))
#endif
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    idx = pos + (k-1)*nj*ni + (j-js)*ni + (i-is) + 1
                    field(i,j,k) = buffer(idx)
                 end do
              end do
           end do
           pos = pos + ksize*nj*ni
        end do
     else if( group%unpack_type(n) == FIELD_X ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldx = group%addrs_x(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(3) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+ksize*nj*ni)) &
           !$omp   map(from: fieldx(is:ie,js:je,1:ksize)) if(use_device_ptr)
#endif
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    idx = pos + (k-1)*nj*ni + (j-js)*ni + (i-is) + 1
                    fieldx(i,j,k) = buffer(idx)
                 end do
              end do
           end do
           pos = pos + ksize*nj*ni
        end do
     else if( group%unpack_type(n) == FIELD_Y ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldy = group%addrs_y(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(3) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+ksize*nj*ni)) &
           !$omp   map(from: fieldy(is:ie,js:je,1:ksize)) if(use_device_ptr)
#endif
           do k = 1, ksize
              do j = js, je
                 do i = is, ie
                    idx = pos + (k-1)*nj*ni + (j-js)*ni + (i-is) + 1
                    fieldy(i,j,k) = buffer(idx)
                 end do
              end do
           end do
           pos = pos + ksize*nj*ni
        end do
     endif
  enddo
else
#ifndef __NVCOMPILER_OPENMP_GPU
!$OMP parallel do default(none) shared(nunpack,group,nscalar,ptr,nvector,ksize,buffer_start_pos) &
!$OMP                          private(buffer_pos,pos,m,is, ie, js, je,rotation, &
!$OMP                                  ptr_field, ptr_fieldx, ptr_fieldy,n,k,ni,nj,idx)
#endif
  do nk = nunpack*ksize, 1, -1
     n = (nk-1)/ksize + 1
     k = mod((nk-1), ksize) + 1
     buffer_pos = group%unpack_buffer_pos(n) + buffer_start_pos
     pos = buffer_pos + (k-1)*group%unpack_size(n)
     is = group%unpack_is(n); ie = group%unpack_ie(n); ni = ie-is+1
     js = group%unpack_js(n); je = group%unpack_je(n); nj = je-js+1
     if( group%unpack_type(n) == FIELD_S ) then
        do l=1,nscalar  ! loop over number of fields
           ptr_field = group%addrs_s(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(2) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+nj*ni)) map(from: field(is:ie,js:je,k)) if(use_device_ptr)
#endif
           do j = js, je
              do i = is, ie
                 idx = pos + (j-js)*ni + (i-is) + 1
                 field(i,j,k) = buffer(idx)
              end do
           end do
           pos = pos + ni*nj
        end do
     else if( group%unpack_type(n) == FIELD_X ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldx = group%addrs_x(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(2) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+nj*ni)) map(from: fieldx(is:ie,js:je,k)) if(use_device_ptr)
#endif
           do j = js, je
              do i = is, ie
                 idx = pos + (j-js)*ni + (i-is) + 1
                 fieldx(i,j,k) = buffer(idx)
              end do
           end do
           pos = pos + ni*nj
        end do
     else if( group%unpack_type(n) == FIELD_Y ) then
        do l=1,nvector  ! loop over number of fields
           ptr_fieldy = group%addrs_y(l)
#ifdef __NVCOMPILER_OPENMP_GPU
           !$omp target teams distribute parallel do collapse(2) private(idx) &
           !$omp   map(to: buffer(pos+1:pos+nj*ni)) map(from: fieldy(is:ie,js:je,k)) if(use_device_ptr)
#endif
           do j = js, je
              do i = is, ie
               idx = pos + (j-js)*ni + (i-is) + 1
               fieldy(i,j,k) = buffer(idx)
              end do
           end do
           pos = pos + ni*nj
        end do
     endif
  enddo
endif
