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
!> @author Jessica Liptak
!> @brief Test mpp_update_domains_ad using different layouts and data precision
program test_mpp_update_domains_ad
  use mpp_mod, only : FATAL, WARNING, NOTE
  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error
  use mpp_mod, only : mpp_set_stack_size
  use mpp_mod, only : mpp_transmit, mpp_sum, mpp_sync
  use mpp_mod, only : mpp_init_test_requests_allocated
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN
  use mpp_domains_mod, only : CGRID_NE, MPP_DOMAIN_TIME
  use mpp_domains_mod, only : domain2D
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad, mpp_check_field
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_global_domain
  use mpp_io_mod, only : mpp_io_init
  use platform_mod, only : r4_kind, r8_kind

  implicit none

  integer :: ierr, id
  integer :: pe, npes
  integer :: nx=64, ny=64, nz=10, stackmax=10000000
  integer :: i, j, k, n
  integer :: layout(2)
  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  !> Initialize mpp and mpp IO modules
  call mpp_init(test_level=mpp_init_test_requests_allocated)
  call mpp_domains_init(MPP_DOMAIN_TIME)
  call mpp_io_init()
  call mpp_domains_set_stack_size(stackmax)
  pe = mpp_pe()
  npes = mpp_npes()
  !> run the tests
  if (pe == mpp_root_pe()) &
    print *, '--------------------> Calling test_halo_update_ad_r8(Simple) <-------------------'
  call test_halo_update_ad_r8('Simple')

  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling test_halo_update_ad_r4(Simple) <-------------------'
  call test_halo_update_ad_r4('Simple')

  call mpp_domains_exit()
  !> Finalize mpp
  call MPI_FINALIZE(ierr)
contains
  !> test calling mpp_halo_update_ad on a 3D 64-bit real data array
  subroutine test_halo_update_ad_r8( test_type )
    character(len=*), intent(in) :: test_type
    ! local
    type(domain2D) :: domain
    integer              :: shift, i, j, k
    logical              :: is_symmetry
    integer              :: is, ie, js, je, isd, ied, jsd, jed, pe
    real(kind=r8_kind), allocatable, dimension(:,:,:) :: x_ad, y_ad, x_fd, y_fd, x_save, y_save
    real(kind=r8_kind) :: ad_sum, fd_sum, sum_diff

    if(index(test_type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(test_type)
    case( 'Simple', 'Simple symmetry' )
      call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
      call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                              shalo=shalo, nhalo=nhalo, name=test_type, symmetry = is_symmetry )
    case default
      call mpp_error( FATAL, 'test_mpp_update_domains_ad_r8: '//test_type//' is not a valid test.')
    end select

!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain( domain, isd, ied, jsd, jed )

    shift=1
!---test 3d single fields----------------------------------------------------------
    allocate( x_fd(isd:ied,jsd:jed,nz) )
    allocate( x_ad(isd:ied,jsd:jed,nz) )
    allocate( x_save(isd:ied,jsd:jed,nz) )
    x_fd = 0.0; x_ad = 0.0; x_save = 0.0

    do k = 1,nz
      do j = js,je
        do i = is,ie
          x_fd(i,j,k) = i*j
        end do
      end do
    end do
    x_save = x_fd

    ! full update
    call mpp_update_domains( x_fd, domain )

    fd_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied
          fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    call mpp_update_domains_ad( x_ad, domain )

    ad_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied
          ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )
    call mpp_sync()
    pe = mpp_pe()
    sum_diff = 0.0
    sum_diff = abs(ad_sum-fd_sum)/fd_sum

    if( pe.EQ.mpp_root_pe() ) then
      if (sum_diff .lt. 1E-7) then
        call MPP_ERROR(NOTE, "Passed Adjoint Dot Test: mpp_update_domains_ad_r8(single 3D field)")
      else
        call MPP_ERROR(FATAL, "FAILED Adjoint Dot Test: mpp_update_domains_ad_r8(single 3D field)")
      endif
    endif

    deallocate (x_ad, x_fd, x_save)

    ! test 3d vector fields
    allocate( x_ad  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_fd  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_save(isd:ied+shift,jsd:jed  ,nz) )
    allocate( y_ad  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_fd  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_save(isd:ied  ,jsd:jed+shift,nz) )

    x_fd=0; y_fd=0
    do k = 1,nz
      do j = js,je
        do i = is,ie
          x_fd(i,j,k)=i*j
          y_fd(i,j,k)=i*j
        end do
      end do
    end do

    call mpp_update_domains( x_fd, y_fd, domain, gridtype=CGRID_NE)
    x_save=x_fd
    y_save=y_fd

    fd_sum = 0.
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
          fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
          fd_sum = fd_sum + y_fd(i,j,k)*y_fd(i,j,k)
        end do
      end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    y_ad = y_fd
    call mpp_update_domains_ad( x_ad, y_ad, domain, gridtype=CGRID_NE)

    ad_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
          ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
          ad_sum = ad_sum + y_ad(i,j,k)*y_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )
    call mpp_sync()
    pe = mpp_pe()
    sum_diff = 0.0
    sum_diff = abs(ad_sum-fd_sum)/fd_sum

    if ( pe.EQ.mpp_root_pe() ) then
      if (sum_diff .lt. 1E-7) then
        call MPP_ERROR(NOTE, "Passed Adjoint Dot Test: mpp_update_domains_ad_r8(vector 3D fields)")
      else
        call MPP_ERROR(FATAL,"FAILED Adjoint Dot Test: mpp_update_domains_ad_r8(vector 3D fields)")
      endif
    endif
    deallocate (x_ad, y_ad, x_fd, y_fd, x_save, y_save)

  end subroutine test_halo_update_ad_r8

  !> test calling mpp_halo_update_ad on a 3D 32-bit real data array
  subroutine test_halo_update_ad_r4( test_type )
    character(len=*), intent(in) :: test_type
    ! local
    type(domain2D) :: domain
    integer              :: shift, i, j, k
    logical              :: is_symmetry
    integer              :: is, ie, js, je, isd, ied, jsd, jed, pe
    real(kind=r4_kind), allocatable, dimension(:,:,:) :: x_ad, y_ad, x_fd, y_fd, x_save, y_save
    real(kind=r4_kind) :: ad_sum, fd_sum, sum_diff

    if(index(test_type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(test_type)
    case( 'Simple', 'Simple symmetry' )
      call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
      call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                              shalo=shalo, nhalo=nhalo, name=test_type, symmetry = is_symmetry )
    case default
      call mpp_error( FATAL, 'test_mpp_update_domains_ad_r4: '//test_type//' is not a valid test.')
    end select

    ! set up the x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain( domain, isd, ied, jsd, jed )

    shift=1
!---test 3d single fields----------------------------------------------------------
    allocate( x_fd(isd:ied,jsd:jed,nz) )
    allocate( x_ad(isd:ied,jsd:jed,nz) )
    allocate( x_save(isd:ied,jsd:jed,nz) )
    x_fd = 0.0; x_ad = 0.0; x_save = 0.0

    do k = 1,nz
      do j = js,je
        do i = is,ie
          x_fd(i,j,k) = i*j
        end do
      end do
    end do
    x_save = x_fd

    ! full update
    call mpp_update_domains( x_fd, domain )

    fd_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied
          fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    call mpp_update_domains_ad( x_ad, domain )

    ad_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied
          ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )
    call mpp_sync()
    pe = mpp_pe()
    sum_diff = 0.0
    sum_diff = abs(ad_sum-fd_sum)/fd_sum

    if( pe.EQ.mpp_root_pe() ) then
      !> @note: threshold differs from R8 test threshold because of expected
      !! increase in roundoff with reduced precision
      if (sum_diff .lt. 1E-6) then
        call MPP_ERROR(NOTE, "Passed Adjoint Dot Test: mpp_update_domains_ad_r4(single 3D field)")
      else
        call MPP_ERROR(NOTE, "FAILED Adjoint Dot Test: mpp_update_domains_ad_r4(single 3D field):")
      endif
    endif

    deallocate (x_ad, x_fd, x_save)

    ! test 3d vector fields
    allocate( x_ad  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_fd  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_save(isd:ied+shift,jsd:jed  ,nz) )
    allocate( y_ad  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_fd  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_save(isd:ied  ,jsd:jed+shift,nz) )

    x_fd=0; y_fd=0
    do k = 1,nz
      do j = js,je
        do i = is,ie
          x_fd(i,j,k)=i*j
          y_fd(i,j,k)=i*j
        end do
      end do
    end do

    call mpp_update_domains( x_fd, y_fd, domain, gridtype=CGRID_NE)
    x_save=x_fd
    y_save=y_fd

    fd_sum = 0.
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
          fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
          fd_sum = fd_sum + y_fd(i,j,k)*y_fd(i,j,k)
        end do
      end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    y_ad = y_fd
    call mpp_update_domains_ad( x_ad, y_ad, domain, gridtype=CGRID_NE)

    ad_sum = 0.0
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
          ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
          ad_sum = ad_sum + y_ad(i,j,k)*y_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )
    call mpp_sync()
    pe = mpp_pe()
    sum_diff = 0.0
    sum_diff = abs(ad_sum-fd_sum)/fd_sum

    if ( pe.EQ.mpp_root_pe() ) then
      !> @note: threshold differs from R8 test threshold because of expected increase
      !! in roundoff error with reduced precision
      if (sum_diff .lt. 1E-6) then
        call MPP_ERROR(NOTE, "Passed Adjoint Dot Test: mpp_update_domains_ad_r4(vector 3D fields)")
      else
        call MPP_ERROR(NOTE, "FAILED Adjoint Dot Test: mpp_update_domains_ad_r4(vector 3D fields)")
      endif
    endif
    deallocate (x_ad, y_ad, x_fd, y_fd, x_save, y_save)

  end subroutine test_halo_update_ad_r4
end program test_mpp_update_domains_ad
