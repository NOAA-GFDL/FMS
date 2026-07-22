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
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif

!> @author Andrew Brooks
!> @brief Test mpp_pelist_gather and mpp_pelist_scatter routines for generalized indices.
program test_mpp_pelist_gatscat_gen_ind
  use mpp_mod
  use platform_mod

  implicit none

  integer :: pe, npes, root
  integer :: storage_to_axis(3)
  integer :: perms(3,6)
  integer :: p

  call mpp_init(mpp_init_test_requests_allocated)
  call mpp_set_stack_size(3145746)

  pe   = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()

  perms = reshape([ &
        1, 2, 3, &
        1, 3, 2, &
        2, 1, 3, &
        2, 3, 1, &
        3, 1, 2, &
        3, 2, 1  ], [3, 6])

  if (pe == root) print *, '--- PELIST SCATTER/GATHER TESTS ---'

  do p = 1, 6

    storage_to_axis = perms(:,p)

    if (pe == root) then
      print *, '--------------------------------'
      print *, 'storage_to_axis =', storage_to_axis
    endif

    call test_scatter(npes, pe, root, storage_to_axis)
    call test_gather(npes, pe, root, storage_to_axis)

  enddo

  if (pe == root) print *, 'ALL PERMUTATIONS PASSED'

  call mpp_exit()

contains

!> @brief Map logical indices (i,j,k) to a unique scalar test value.
   pure real function val(i, j, k)
     integer, intent(in)  :: i, j, k

     val = 100.0*i + 10.0*j + k
   end function val

!> @brief Convert logical indices (i,j,k) into storage indices using storage_to_axis.
   subroutine permute(i, j, k, storage_to_axis, u, v, w)
     integer, intent(in)  :: i, j, k, storage_to_axis(3)
     integer, intent(out) :: u, v, w

     integer :: idx(3)

     idx = (/i, j, k/)
     u = idx(storage_to_axis(1))
     v = idx(storage_to_axis(2))
     w = idx(storage_to_axis(3))
   end subroutine permute

!> @brief Compute 1D domain decomposition in i-direction for each MPI rank.
   subroutine get_decomp(pe, npes, NI, NJ, is, ie, js, je)
     integer, intent(in)  :: pe, npes, NI, NJ
     integer, intent(out) :: is, ie, js, je

     integer :: chunk

     chunk = NI / npes
     is = pe*chunk + 1
     ie = (pe+1)*chunk
     js = 1
     je = NJ
   end subroutine get_decomp

!> @brief Construct pelist containing all MPI ranks participating in the operation.
  subroutine build_pelist(npes, pelist)
     integer, intent(in)  :: npes
     integer, intent(out) :: pelist(npes)

     integer :: i

     do i=0,npes-1
       pelist(i+1)=i
     enddo
   end subroutine build_pelist

!> @brief Allocate a 3D field in permuted layout, handling root and non-root cases.
   subroutine alloc_field(field, storage_to_axis, dims_logical, pe, root, is_global)
     real, allocatable, intent(out) :: field(:,:,:)
     integer, intent(in) :: storage_to_axis(3)
     integer, intent(in) :: dims_logical(3)
     integer, intent(in) :: pe, root
     logical, intent(in) :: is_global

     integer :: dims(3)

     if (is_global) then
       if (pe == root) then
         dims = dims_logical
         allocate(field(dims(storage_to_axis(1)), &
                        dims(storage_to_axis(2)), &
                        dims(storage_to_axis(3))))
       else
         allocate(field(1,1,1))
       endif
     else
       dims = dims_logical
       allocate(field(dims(storage_to_axis(1)), &
                      dims(storage_to_axis(2)), &
                      dims(storage_to_axis(3))))
     endif
   end subroutine alloc_field

!> @brief Populate a field with values from val(i,j,k) under the given permutation.
   subroutine fill_from_val(segment, is, ie, NJ, NK, storage_to_axis)
     real, intent(inout) :: segment(:,:,:)
     integer, intent(in) :: is, ie, NJ, NK, storage_to_axis(3)

     integer :: i, j, k
     integer :: u, v, w

     do i=is,ie
       do j=1,NJ
         do k=1,NK
           call permute(i-is+1, j, k, storage_to_axis, u, v, w)
           segment(u,v,w) = val(i,j,k)
         enddo
       enddo
     enddo
   end subroutine fill_from_val

!> @brief Verify field values against val(i,j,k) for local (scatter) or global (gather) domains.
   subroutine check_answer(field, is, ie, NI, NJ, NK, storage_to_axis, check_global)
     real, intent(in) :: field(:,:,:)
     integer, intent(in) :: is, ie, NI, NJ, NK
     integer, intent(in) :: storage_to_axis(3)
     logical, intent(in) :: check_global   ! .true. = gather, .false. = scatter

     integer :: i, j, k
     integer :: u, v, w
     integer :: istart, iend, iloc

     ! --- choose iteration domain ---
     if (check_global) then
       istart = 1
       iend   = NI
     else
       istart = is
       iend   = ie
     endif

     do i = istart, iend
       do j = 1, NJ
         do k = 1, NK

           if (check_global) then
             ! gather -> global indexing
             call permute(i, j, k, storage_to_axis, u, v, w)
           else
             ! scatter -> local indexing
             iloc = i - is + 1
             call permute(iloc, j, k, storage_to_axis, u, v, w)
           endif

           if (field(u,v,w) /= val(i,j,k)) then
             print *, 'FAIL'
             print *, 'i,j,k=', i,j,k
             print *, 'u,v,w=', u,v,w
             print *, 'got=', field(u,v,w)
             print *, 'expected=', val(i,j,k)
             call mpp_error(FATAL,'check_answer failed')
           endif

         enddo
       enddo
     enddo
   end subroutine check_answer

!> @brief Test mpp_scatter with pelist and permuted storage layouts.
   subroutine test_scatter(npes, pe, root, storage_to_axis)
     integer,intent(in) :: npes, pe, root
     integer,intent(in) :: storage_to_axis(3)

     integer :: pelist(npes)
     integer :: is, ie, js, je
     integer :: NI, NJ, NK
     integer :: axis_to_storage(3)
     real, allocatable :: global_perm(:,:,:)
     real, allocatable :: segment(:,:,:)

     NI = 8
     NJ = 5
     NK = 3

     call get_decomp(pe, npes, NI, NJ, is, ie, js, je)
     call build_pelist(npes, pelist)
     !call mpp_get_pelist(pelist)

     call alloc_field(global_perm, storage_to_axis, (/NI, NJ, NK/), pe, root, .true.)
     call alloc_field(segment, storage_to_axis, (/ie-is+1, NJ, NK/), pe, root, .false.)

     segment = -2.0

     if (pe == root) then
       call fill_from_val(global_perm, 1, NI, NJ, NK, storage_to_axis)
     endif

     ! Initialize axis_to_storage map
     axis_to_storage(storage_to_axis(1)) = 1
     axis_to_storage(storage_to_axis(2)) = 2
     axis_to_storage(storage_to_axis(3)) = 3

     ! --- scatter ---
     call mpp_sync()
     if (pe == root) then
       call mpp_scatter(is, ie, js, je, NK, pelist, segment, global_perm, axis_to_storage, .true.)
     else
       call mpp_scatter(is, ie, js, je, NK, pelist, segment, global_perm, axis_to_storage, .false.)
     endif
     call mpp_sync()

     call check_answer(segment, is, ie, NI, NJ, NK, storage_to_axis, .false.)

     if (pe == root) print *, 'SCATTER PASS'

     deallocate(segment)
     if (pe == root) deallocate(global_perm)
   end subroutine test_scatter

!> @brief Test mpp_gather with pelist and permuted storage layouts.
   subroutine test_gather(npes, pe, root, storage_to_axis)
     integer,intent(in) :: npes, pe, root
     integer,intent(in) :: storage_to_axis(3)

     integer :: pelist(npes)
     integer :: is, ie, js, je
     integer :: NI,NJ,NK
     integer :: axis_to_storage(3)
     real, allocatable :: segment(:,:,:)
     real, allocatable :: gather_data(:,:,:)

     NI = 8
     NJ = 5
     NK = 3

     call get_decomp(pe, npes, NI, NJ, is, ie, js, je)
     call build_pelist(npes, pelist)

     call alloc_field(gather_data, storage_to_axis, (/NI, NJ, NK/), pe, root, .true.)
     call alloc_field(segment, storage_to_axis, (/ie-is+1, NJ, NK/), pe, root, .false.)

     call fill_from_val(segment, is, ie, NJ, NK, storage_to_axis)

     ! Initialize axis_to_storage map
     axis_to_storage(storage_to_axis(1)) = 1
     axis_to_storage(storage_to_axis(2)) = 2
     axis_to_storage(storage_to_axis(3)) = 3

     ! --- GATHER ---
     call mpp_sync()
     if (pe == root) then
       call mpp_gather(is, ie, js, je, NK, pelist, segment, gather_data, axis_to_storage, .true.)
     else
       call mpp_gather(is, ie, js, je, NK, pelist, segment, gather_data, axis_to_storage, .false.)
     endif
     call mpp_sync()

     if (pe == root) then
       call check_answer(gather_data, is, ie, NI, NJ, NK, storage_to_axis, .true.)
       print *, 'GATHER PASS'
     endif

     deallocate(segment)
     if (pe == root) deallocate(gather_data)

   end subroutine test_gather

end program test_mpp_pelist_gatscat_gen_ind
