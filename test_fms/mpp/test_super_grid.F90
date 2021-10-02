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

!> @brief  This programs tests mpp_copy_domains and makes sure the domain was
!>         copied correctly tests mpp_create_super_grid_domain and makes sure the
!>         domain is correct
program test_super_grid
use   fms_mod,            only: fms_init, fms_end
use   mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, mpp_npes, FATAL
use   mpp_domains_mod,    only: domain2d, mpp_define_domains, mpp_copy_domain
use   mpp_domains_mod,    only: mpp_get_data_domain
use   mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_compute_domains
use   mpp_domains_mod,    only: mpp_get_global_domain, mpp_get_global_domains
use   mpp_domains_mod,    only: mpp_create_super_grid_domain

implicit none

type(domain2d)                     :: Domain            !< Domain
type(domain2d)                     :: Domain2           !< Domain
integer                            :: nlon              !< Number of points in x direction
integer                            :: nlat              !< Number of points in y direction
integer, dimension(2)              :: layout = (/1,6/)  !< Domain layout

call fms_init()

nlon = 360
nlat = 90

call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, xhalo=2, yhalo=2, name='test_supergrid')
call mpp_copy_domain(Domain, Domain2)
call compare_domains(Domain, Domain2, supergrid=.false.)

call mpp_create_super_grid_domain(Domain2)
call compare_domains(Domain, Domain2, supergrid=.true.)

call fms_end()

contains

subroutine compare_domains(Domain, Domain2, supergrid)
 type(domain2d), intent(inout)      :: Domain            !< Domain
 type(domain2d), intent(inout)      :: Domain2           !< Domain
 logical,        intent(in)         :: supergrid         !< Flag indicating if comparing supergrid

 integer, allocatable               :: pe_indices (:,:)  !< Indices for current pe
 integer, allocatable               :: pe_indices2(:,:)  !< Indices for current pe
 integer, allocatable               :: pes_indices(:,:)  !< Indices for all pes in the pelist
 integer, allocatable               :: pes_indices2(:,:) !< Indices for all pes in the pelist

 allocate(pes_indices (mpp_npes(),6))
 allocate(pes_indices2(mpp_npes(),6))
 allocate(pe_indices (1,6))
 allocate(pe_indices2(1,6))

 !< Check if the pe_indices are the same for the two domains
 !> Compute Domain:
 call mpp_get_compute_domain(Domain, xbegin=pe_indices (1,1), xend=pe_indices (1,2), &
                                   & ybegin=pe_indices (1,3), yend=pe_indices (1,4), &
                                   & xsize =pe_indices (1,5), ysize=pe_indices(1,6))
 if(supergrid) call get_expected_indices(pe_indices)

 call mpp_get_compute_domain(Domain2, xbegin=pe_indices2(1,1), xend=pe_indices2(1,2), &
                                   &  ybegin=pe_indices2(1,3), yend=pe_indices2(1,4), &
                                   &  xsize =pe_indices2(1,5), ysize=pe_indices2(1,6))
 call compare_pe_indices(pe_indices, pe_indices2, "compute domain")

 !> Data Domain:
 call mpp_get_data_domain(Domain, xbegin=pe_indices (1,1), xend=pe_indices (1,2), &
                                   & ybegin=pe_indices (1,3), yend=pe_indices (1,4), &
                                   & xsize =pe_indices (1,5), ysize=pe_indices(1,6))
 if(supergrid) call get_expected_indices(pe_indices)

 call mpp_get_data_domain(Domain2, xbegin=pe_indices2(1,1), xend=pe_indices2(1,2), &
                                   &  ybegin=pe_indices2(1,3), yend=pe_indices2(1,4), &
                                   &  xsize =pe_indices2(1,5), ysize=pe_indices2(1,6))
 call compare_pe_indices(pe_indices, pe_indices2, "data domain")

 !> Global Domain:
 call mpp_get_global_domain(Domain, xbegin=pe_indices (1,1), xend=pe_indices (1,2), &
                                   & ybegin=pe_indices (1,3), yend=pe_indices (1,4), &
                                   & xsize =pe_indices (1,5), ysize=pe_indices(1,6))
 if(supergrid) call get_expected_indices(pe_indices)

 call mpp_get_global_domain(Domain2, xbegin=pe_indices2(1,1), xend=pe_indices2(1,2), &
                                   &  ybegin=pe_indices2(1,3), yend=pe_indices2(1,4), &
                                   &  xsize =pe_indices2(1,5), ysize=pe_indices2(1,6))
 call compare_pe_indices(pe_indices, pe_indices2, "global domain")

 !> Compute Domains (get the indices for all of the ranks)
 call mpp_get_compute_domains(Domain,  xbegin=pes_indices (:,1), xend=pes_indices (:,2), &
                                     & ybegin=pes_indices (:,3), yend=pes_indices (:,4), &
                                     & xsize =pes_indices (:,5), ysize=pes_indices (:,6))
 if(supergrid) call get_expected_indices(pes_indices)

 call mpp_get_compute_domains(Domain2, xbegin=pes_indices2(:,1), xend=pes_indices2(:,2), &
                                     & ybegin=pes_indices2(:,3), yend=pes_indices2(:,4), &
                                     & xsize =pes_indices2(:,5), ysize=pes_indices2(:,6))
 call compare_pe_indices(pes_indices, pes_indices2, "compute domains")

 !> Global Domains (get the indices for all of the ranks)
 call mpp_get_global_domains(Domain,  xbegin=pes_indices (:,1), xend=pes_indices (:,2), &
                                     & ybegin=pes_indices (:,3), yend=pes_indices (:,4), &
                                     & xsize =pes_indices (:,5), ysize=pes_indices(:,6))
 if(supergrid) call get_expected_indices(pes_indices)

 call mpp_get_global_domains(Domain2, xbegin=pes_indices2(:,1), xend=pes_indices2(:,2), &
                                     & ybegin=pes_indices2(:,3), yend=pes_indices2(:,4), &
                                     & xsize =pes_indices2(:,5), ysize=pes_indices2(:,6))
 call compare_pe_indices(pes_indices, pes_indices2, "global domains")

end subroutine compare_domains

subroutine get_expected_indices(pe_indices)
  integer, intent(inout) :: pe_indices(:,:)  !< Indices from original domain

  integer :: i !< For do loop

  do i=1, size(pe_indices,1)
     pe_indices(i,1) = 2*pe_indices(i,1)-1 !< xbegin
     pe_indices(i,2) = 2*pe_indices(i,2)+1 !< xend
     pe_indices(i,3) = 2*pe_indices(i,3)-1 !< ybegin
     pe_indices(i,4) = 2*pe_indices(i,4)+1 !< yend
     pe_indices(i,5) = pe_indices(i,2) - pe_indices(i,1) +1 !< xsize
     pe_indices(i,6) = pe_indices(i,4) - pe_indices(i,3) +1 !< ysize
  enddo

end subroutine get_expected_indices

subroutine compare_pe_indices(pe_indices, pe_indices2, domain_type)
  integer, intent(in) :: pe_indices(:,:)  !< Indices from original domain
  integer, intent(in) :: pe_indices2(:,:) !< Indices from new domain
  character(len=*), intent(in) :: domain_type !< Type of domain "data", "compute", and "global"

  integer :: i,j !< For do loop

  do i = 1, size(pe_indices, 1)
     do j = 1, size(pe_indices,2)
        if (pe_indices(i, j) .ne. pe_indices2(i, j)) then
           call mpp_error(FATAL, "Compare pe_indices: the "//trim(domain_type)//" pe_indices are not the same ")
        endif
     enddo
  enddo
end subroutine compare_pe_indices

end program test_super_grid
