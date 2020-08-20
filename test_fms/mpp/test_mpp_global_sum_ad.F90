!***********************************************************************
!*                   Gnu Lesser General Public License
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
!************************************************************
!> @author Jessica Liptak
!> @email gfdl.climate.model.info@noaa.gov
!> @description Test the mpp_global_sum_ad interfaces with 32-bit and 64-bit
!! real and integer arrays. mpp_global_sum_ad computes the global adjoint sum of
!! a field.
program test_mpp_global_sum_ad

  use mpp_mod,         only : FATAL, MPP_DEBUG
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sum
  use mpp_mod,         only : mpp_init, stdout, stderr
  use mpp_mod,         only : mpp_get_current_pelist, mpp_broadcast
  use mpp_mod,         only : mpp_init_test_requests_allocated
  use mpp_domains_mod, only : BITWISE_EXACT_SUM
  use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN
  use mpp_domains_mod, only : domain2D
  use mpp_domains_mod, only : mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_sum
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : NORTH, EAST, CORNER, CENTER
  use mpp_domains_mod, only : mpp_global_sum_ad
  use mpp_io_mod,      only : mpp_io_init
  use platform_mod


  implicit none
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=2, stackmax=4000000
  integer :: layout(2)
  integer :: ierr
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2

  call mpp_init(test_level=mpp_init_test_requests_allocated)
  call mpp_domains_init(MPP_DEBUG)
  call mpp_io_init()
  call mpp_domains_set_stack_size(stackmax)

  pe = mpp_pe()
  npes = mpp_npes()

  call test_global_sum_ad_r4( 'Simple')
  call test_global_sum_ad_r4( 'Cyclic symmetry center')

  call test_global_sum_ad_r8( 'Simple')
  call test_global_sum_ad_r8( 'Cyclic symmetry center')

  call test_global_sum_ad_i4( 'Simple')
  call test_global_sum_ad_i4( 'Cyclic symmetry center')

  call test_global_sum_ad_i8( 'Simple')
  call test_global_sum_ad_i8( 'Cyclic symmetry center')


  call MPI_finalize(ierr)

contains

  !> test the 32-bit real global_sum_ad interfaces
  subroutine test_global_sum_ad_r4 (domain_type)
    character(len=*), intent(in) :: domain_type ! type of mpp domain to use
    ! local
    real(r4_kind) :: gsum_tl, gsum_ad
    real(r4_kind) :: gsum_tl_save, gsum_ad_save
    real(r4_kind) :: gsum_tl_bit, gsum_ad_bit
    real(r4_kind) :: gsum_tl_save_bit, gsum_ad_save_bit
    integer :: i,j,k, ishift, jshift, position
    integer :: isd, ied, jsd, jed

    type(domain2D) :: domain
    real(r4_kind), allocatable, dimension(:,:) :: x2, x2_ad, x2_ad_bit
    real(r4_kind), allocatable, dimension(:,:,:) :: x3, x3_ad, x3_ad_bit
    real(r4_kind), allocatable, dimension(:,:,:,:) :: x4, x4_ad, x4_ad_bit
    real(r4_kind), allocatable, dimension(:,:,:,:,:) :: x5, x5_ad, x5_ad_bit

    call generate_domain(domain, domain_type)

    call mpp_get_data_domain( domain, isd, ied, jsd, jed )

    position = CENTER

    ! test the 2D arrays
    allocate( x2(isd:ied,jsd:jed), x2_ad(isd:ied,jsd:jed), x2_ad_bit(isd:ied,jsd:jed) )

    x2=0.0
    do j = jsd, jed
      do i = isd, ied
        x2(i,j) = real(i+j, kind=r4_kind)
      enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x2, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x2, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x2_ad     = 0.
    x2_ad_bit = 0.
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x2_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x2_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum the original global sum and the adjoint global sum
    do j = jsd, jed
      do i = isd, ied
        gsum_ad_save     = gsum_ad_save + x2_ad(i,j)*x2(i,j)
        gsum_ad_save_bit = gsum_ad_save_bit + x2_ad_bit(i,j)*x2(i,j)
      enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif
    deallocate(x2, x2_ad, x2_ad_bit)

    ! test 3D arrays
    allocate( x3(isd:ied,jsd:jed,nz), x3_ad(isd:ied,jsd:jed,nz), x3_ad_bit(isd:ied,jsd:jed,nz) )

    x3=0.0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x3(i,j,k) = real(i+j+k, kind=r4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x3, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x3, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x3_ad     = 0.
    x3_ad_bit = 0.
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x3_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x3_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x3_ad(i,j,k)*x3(i,j,k)
           gsum_ad_save_bit = gsum_ad_save_bit + x3_ad_bit(i,j,k)*x3(i,j,k)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x3, x3_ad, x3_ad_bit)

    ! test 4D arrays
    allocate( x4(isd:ied,jsd:jed,nz,1), x4_ad(isd:ied,jsd:jed,nz,1), x4_ad_bit(isd:ied,jsd:jed,nz,1) )

    x4=0.0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x4(i,j,k,1) = real(i+j+k, kind=r4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x4, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x4, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x4_ad     = 0.
    x4_ad_bit = 0.
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x4_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x4_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x4_ad(i,j,k,1)*x4(i,j,k,1)
           gsum_ad_save_bit = gsum_ad_save_bit + x4_ad_bit(i,j,k,1)*x4(i,j,k,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x4, x4_ad, x4_ad_bit)

    ! test 5D arrays
    allocate( x5(isd:ied,jsd:jed,nz,1,1), x5_ad(isd:ied,jsd:jed,nz,1,1), &
              x5_ad_bit(isd:ied,jsd:jed,nz,1,1) )

    x5=0.0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x5(i,j,k,1,1) = real(i+j+k, kind=r4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x5, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x5, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit
    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x5_ad     = 0.
    x5_ad_bit = 0.
    ! adjoint sum of the global field
    call mpp_global_sum_ad( domain, x5_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x5_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x5_ad(i,j,k,1,1)*x5(i,j,k,1,1)
           gsum_ad_save_bit = gsum_ad_save_bit + x5_ad_bit(i,j,k,1,1)*x5(i,j,k,1,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r4,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x5, x5_ad, x5_ad_bit)

  end subroutine test_global_sum_ad_r4

  !> test 64-bit real global_sum_ad interfaces
  subroutine test_global_sum_ad_r8 (domain_type)
    character(len=*), intent(in) :: domain_type ! type of mpp domain to use
    ! local
    type(domain2D) :: domain
    real(r8_kind) :: gsum_tl, gsum_ad
    real(r8_kind) :: gsum_tl_save, gsum_ad_save
    real(r8_kind) :: gsum_tl_bit, gsum_ad_bit
    real(r8_kind) :: gsum_tl_save_bit, gsum_ad_save_bit
    integer :: i,j,k, ishift, jshift, position
    integer :: isd, ied, jsd, jed

    real(r8_kind), allocatable, dimension(:,:) :: x2, x2_ad, x2_ad_bit
    real(r8_kind), allocatable, dimension(:,:,:) :: x3, x3_ad, x3_ad_bit
    real(r8_kind), allocatable, dimension(:,:,:,:) :: x4, x4_ad, x4_ad_bit
    real(r8_kind), allocatable, dimension(:,:,:,:,:) :: x5, x5_ad, x5_ad_bit

    call generate_domain(domain, domain_type)
    call mpp_get_data_domain( domain, isd, ied, jsd, jed )
    position = CENTER

    ! test the 2D arrays
    allocate( x2(isd:ied,jsd:jed), x2_ad(isd:ied,jsd:jed), x2_ad_bit(isd:ied,jsd:jed) )

    x2=0.0
    do j = jsd, jed
      do i = isd, ied
        x2(i,j) = real(i+j, kind=r8_kind)
      enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0

    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x2, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x2, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x2_ad     = 0.
    x2_ad_bit = 0.
    ! adjoint sum of the global field
    call mpp_global_sum_ad( domain, x2_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x2_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum the the original global sum and the adjoint global sum
    do j = jsd, jed
      do i = isd, ied
        gsum_ad_save     = gsum_ad_save + x2_ad(i,j)*x2(i,j)
        gsum_ad_save_bit = gsum_ad_save_bit + x2_ad_bit(i,j)*x2(i,j)
      enddo
    enddo

    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif
    deallocate(x2, x2_ad, x2_ad_bit)

    ! test 3D arrays
    allocate( x3(isd:ied,jsd:jed,nz), x3_ad(isd:ied,jsd:jed,nz), x3_ad_bit(isd:ied,jsd:jed,nz) )

    x3=0.
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x3(i,j,k) = real(i+j+k, kind=r8_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x3, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x3, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x3_ad     = 0.
    x3_ad_bit = 0.
    ! adjoint sum of the global field
    call mpp_global_sum_ad( domain, x3_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x3_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum of the global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x3_ad(i,j,k)*x3(i,j,k)
           gsum_ad_save_bit = gsum_ad_save_bit + x3_ad_bit(i,j,k)*x3(i,j,k)
         enddo
       enddo
    enddo
    ! sum across all pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x3, x3_ad, x3_ad_bit)

    ! test 4D arrays
    allocate( x4(isd:ied,jsd:jed,nz,1), x4_ad(isd:ied,jsd:jed,nz,1), x4_ad_bit(isd:ied,jsd:jed,nz,1) )

    x4=0.
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x4(i,j,k,1) = real(i+j+k, kind=r8_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x4, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x4, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x4_ad     = 0.
    x4_ad_bit = 0.
    ! ajoint sum of the global field
    call mpp_global_sum_ad( domain, x4_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x4_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum of the adjoint global sum and the original global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x4_ad(i,j,k,1)*x4(i,j,k,1)
           gsum_ad_save_bit = gsum_ad_save_bit + x4_ad_bit(i,j,k,1)*x4(i,j,k,1)
         enddo
       enddo
    enddo
    ! sum across all pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x4, x4_ad)

    ! test 5D arrays
    allocate( x5(isd:ied,jsd:jed,nz,1,1), x5_ad(isd:ied,jsd:jed,nz,1,1), &
              x5_ad_bit(isd:ied,jsd:jed,nz,1,1) )

    x5=0.
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x5(i,j,k,1,1) = real(i+j+k, kind=r8_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0.0
    gsum_tl = 0.0
    gsum_tl_save = 0.0
    gsum_tl_bit = 0.0
    gsum_tl_save_bit = 0.0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x5, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x5, flags=BITWISE_EXACT_SUM  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x5_ad     = 0.
    x5_ad_bit = 0.
    ! global adjoint sum
    call mpp_global_sum_ad( domain, x5_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x5_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.
    ! sum of the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x5_ad(i,j,k,1,1)*x5(i,j,k,1,1)
           gsum_ad_save_bit = gsum_ad_save_bit + x5_ad_bit(i,j,k,1,1)*x5(i,j,k,1,1)
         enddo
       enddo
    enddo
    ! sum across all pes
    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (abs((gsum_ad_save-gsum_tl_save)/gsum_tl_save).lt.1E-7) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8"
      endif
      if (abs((gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit).lt.1E-7) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_r8,"//&
                "flags=BITWISE_EXACT_SUM"
      endif
    endif

    deallocate(x5, x5_ad, x5_ad_bit)

  end subroutine test_global_sum_ad_r8

  !> test the 32-bit integer global_sum_ad interfaces
  subroutine test_global_sum_ad_i4 (domain_type)
    character(len=*), intent(in) :: domain_type ! type of mpp domain to use
    ! local
    integer(i4_kind) :: gsum_tl, gsum_ad
    integer(i4_kind) :: gsum_tl_save, gsum_ad_save
    integer :: i,j,k, ishift, jshift, position
    integer :: isd, ied, jsd, jed

    type(domain2D) :: domain
    integer(i4_kind), allocatable, dimension(:,:) :: x2, x2_ad
    integer(i4_kind), allocatable, dimension(:,:,:) :: x3, x3_ad
    integer(i4_kind), allocatable, dimension(:,:,:,:) :: x4, x4_ad
    integer(i4_kind), allocatable, dimension(:,:,:,:,:) :: x5, x5_ad

    call generate_domain(domain, domain_type)

    call mpp_get_data_domain( domain, isd, ied, jsd, jed )

    position = CENTER

    ! test the 2D arrays
    allocate( x2(isd:ied,jsd:jed), x2_ad(isd:ied,jsd:jed))

    x2=0
    do j = jsd, jed
      do i = isd, ied
        x2(i,j) = int(i+j, kind=i4_kind)
      enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x2, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x2_ad     = 0

    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x2_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the original global sum and the adjoint global sum
    do j = jsd, jed
      do i = isd, ied
        gsum_ad_save     = gsum_ad_save + x2_ad(i,j)*x2(i,j)
      enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i4"
      endif
    endif
    deallocate(x2, x2_ad)

    ! test 3D arrays
    allocate( x3(isd:ied,jsd:jed,nz), x3_ad(isd:ied,jsd:jed,nz))

    x3 = 0
    do k = 1,nz
       do j = jsd,jed
         do i = isd, ied
           x3(i,j,k) = int(i+j+k, kind=i4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x3, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x3_ad     = 0
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x3_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x3_ad(i,j,k)*x3(i,j,k)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i4"
      endif
    endif

    deallocate(x3, x3_ad)

    ! test 4D arrays
    allocate( x4(isd:ied,jsd:jed,nz,1), x4_ad(isd:ied,jsd:jed,nz,1))

    x4=0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x4(i,j,k,1) = int(i+j+k, kind=i4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x4, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x4_ad     = 0
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x4_ad, gsum_ad, position = position )
    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save = gsum_ad_save + x4_ad(i,j,k,1)*x4(i,j,k,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i4"
      endif
    endif

    deallocate(x4, x4_ad)

    ! test 5D arrays
    allocate( x5(isd:ied,jsd:jed,nz,1,1), x5_ad(isd:ied,jsd:jed,nz,1,1))

    x5=0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x5(i,j,k,1,1) = int(i+j+k, kind=i4_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x5, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl
    gsum_ad      = gsum_tl

    x5_ad     = 0
    ! adjoint sum of the global field
    call mpp_global_sum_ad( domain, x5_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save = gsum_ad_save + x5_ad(i,j,k,1,1)*x5(i,j,k,1,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i4"
      endif
    endif

    deallocate(x5, x5_ad)

  end subroutine test_global_sum_ad_i4

  !> test the 64-bit integer global_sum_ad interfaces
  subroutine test_global_sum_ad_i8 (domain_type)
    character(len=*), intent(in) :: domain_type ! type of mpp domain to use
    ! local
    integer(i8_kind) :: gsum_tl, gsum_ad
    integer(i8_kind) :: gsum_tl_save, gsum_ad_save
    integer :: i,j,k, ishift, jshift, position
    integer :: isd, ied, jsd, jed

    type(domain2D) :: domain
    integer(i8_kind), allocatable, dimension(:,:) :: x2, x2_ad
    integer(i8_kind), allocatable, dimension(:,:,:) :: x3, x3_ad
    integer(i8_kind), allocatable, dimension(:,:,:,:) :: x4, x4_ad
    integer(i8_kind), allocatable, dimension(:,:,:,:,:) :: x5, x5_ad

    call generate_domain(domain, domain_type)

    call mpp_get_data_domain( domain, isd, ied, jsd, jed )

    position = CENTER

    ! test the 2D arrays
    allocate( x2(isd:ied,jsd:jed), x2_ad(isd:ied,jsd:jed))

    x2=0
    do j = jsd, jed
      do i = isd, ied
        x2(i,j) = int(i+j, kind=i8_kind)
      enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x2, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x2_ad     = 0

    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x2_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the original global sum and the adjoint global sum
    do j = jsd, jed
      do i = isd, ied
        gsum_ad_save     = gsum_ad_save + x2_ad(i,j)*x2(i,j)
      enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "2D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i8"
      endif
    endif
    deallocate(x2, x2_ad)

    ! test 3D arrays
    allocate( x3(isd:ied,jsd:jed,nz), x3_ad(isd:ied,jsd:jed,nz))

    x3=0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x3(i,j,k) = int(i+j+k, kind=i8_kind)
         enddo
       enddo
    enddo

    gsum_tl = 0
    gsum_ad = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x3, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x3_ad     = 0
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x3_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x3_ad(i,j,k)*x3(i,j,k)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "3D arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i8"
      endif
    endif

    deallocate(x3, x3_ad)

    ! test 4D arrays
    allocate( x4(isd:ied,jsd:jed,nz,1), x4_ad(isd:ied,jsd:jed,nz,1))

    x4=0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x4(i,j,k,1) = int(i+j+k, kind=i8_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x4, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl

    gsum_ad      = gsum_tl

    x4_ad     = 0
    ! adjoint sum of global field
    call mpp_global_sum_ad( domain, x4_ad, gsum_ad, position = position )
    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x4_ad(i,j,k,1)*x4(i,j,k,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "4d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i8"
      endif
    endif

    deallocate(x4, x4_ad)

    ! test 5D arrays
    allocate( x5(isd:ied,jsd:jed,nz,1,1), x5_ad(isd:ied,jsd:jed,nz,1,1))

    x5=0
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x5(i,j,k,1,1) = int(i+j+k, kind=i8_kind)
         enddo
       enddo
    enddo

    gsum_ad = 0
    gsum_tl = 0
    gsum_tl_save = 0
    ! global sum of the domain-decomposed array
    gsum_tl      = mpp_global_sum( domain, x5, position = position  )

    gsum_tl_save = gsum_tl*gsum_tl
    gsum_ad      = gsum_tl

    x5_ad     = 0
    ! adjoint sum of the global field
    call mpp_global_sum_ad( domain, x5_ad, gsum_ad, position = position )

    gsum_ad_save     = 0
    ! sum the the original global sum and the adjoint global sum
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save = gsum_ad_save + x5_ad(i,j,k,1,1)*x5(i,j,k,1,1)
         enddo
       enddo
    enddo
    ! sum across the pes
    call mpp_sum( gsum_ad_save )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
      if (gsum_ad_save .eq. gsum_tl_save) then
        print*, "5d arrays Passed Adjoint Dot Test: mpp_global_sum_ad_i8"
      endif
    endif

    deallocate(x5, x5_ad)

  end subroutine test_global_sum_ad_i8

  !> define the 2D test domain
  subroutine generate_domain(domain, domain_type)
    type(domain2D), intent(inout) :: domain ! 2D mpp domain
    character(len=*), intent(in) :: domain_type ! type of domain to generate
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )

    select case(trim(domain_type))
      case( 'Simple' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, &
             ehalo=ehalo, shalo=shalo, nhalo=nhalo, name=domain_type )
      case( 'Cyclic symmetry center')
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, &
             ehalo=ehalo, shalo=shalo, nhalo=nhalo, name=domain_type, symmetry = .true., &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
      case default
        call mpp_error( FATAL, 'test_mpp_global_sum_ad: no such test: '//trim(domain_type))
    end select
  end subroutine generate_domain

end program test_mpp_global_sum_ad
