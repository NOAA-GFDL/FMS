! Corner mosaic test with resource reallocation
! This test constructs a single-tile domain, then a three-tile mosaic, and goes back down to a single-tile domain.
! The number of PEs must be such that the single-tile domain uses all of the PEs (e.g., 64 PEs if the single-tile
! domain is 8x8). The layout of the corner mosaic is defined by ni and nj. The corner mosaic does not need to use all
! of the available PEs.

program test_corner_mosaic
  use mpp_mod, only: mpp_error, FATAL
  use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_declare_pelist, mpp_set_current_pelist
  use fms_mod, only: fms_init, fms_end
  use mpp_domains_mod, only: domain2d, mpp_define_mosaic, mpp_get_compute_domain, mpp_get_data_domain, &
                             mpp_get_tile_id, mpp_update_domains, mpp_domains_set_stack_size, CGRID_NE, &
                             DGRID_NE, mpp_define_domains, mpp_deallocate_domain
  use platform_mod
  use fms_string_utils_mod, only: string
  use random_numbers_mod
  implicit none

  type contact
    integer :: tile(2)
    integer :: is(2)
    integer :: ie(2)
    integer :: js(2)
    integer :: je(2)
  end type contact

  integer, parameter :: ntiles = 3
  integer, parameter :: num_contacts = 3
  integer, parameter :: halo = 2

  character(*), parameter :: name = "corner mosaic"
  integer :: global_indices(4,ntiles), layout(2,ntiles)
  integer, dimension(ntiles) :: pe_start, pe_end

  integer, parameter :: n=16
  integer, parameter :: layout_st(2) = [8, 8] ! Single-tile layout

  ! Corner mosaic layout:
  ! - Tile 1: 4x4
  ! - Tile 2: 6x4
  ! - Tile 3: 4x6
  integer, parameter :: ni(3) = [4, 6, 4] ! Corner mosaic layout (i axis)
  integer, parameter :: nj(3) = [4, 4, 6] ! Corner mosaic layout (j axis)

  real(r4_kind), parameter :: eps = 1e-3
  real(r4_kind), parameter :: halo_init_value = -123456._r4_kind
  integer, parameter :: stackmax=10000000

  integer :: isg, ieg, jsg, jeg
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: tid

  call fms_init
  call mpp_domains_set_stack_size(stackmax)

  call define_single_tile
  call define_corner_mosaic

  call mpp_set_current_pelist ! Restore original pelist
  call define_single_tile

  call fms_end

contains

  !!!!!!!!!
  !   --- !
  !  /|3| !
  ! ----- !
  ! |1|2| !
  ! ----- !
  !!!!!!!!!

  subroutine define_single_tile
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: tile_id
    integer :: npes

    npes = layout_st(1) * layout_st(2)
    if (mpp_npes().ne.npes) then
      call mpp_error(FATAL, "Cannot construct single-tile domain: Incorrect number of PEs")
    endif

    isg=1
    ieg=2*n
    jsg=1
    jeg=2*n

    call mpp_define_domains([isg,ieg,jsg,jeg], layout_st, domain, &
        symmetry=.true., whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    call single_tile_scalar_test(domain)
    call mpp_deallocate_domain(domain)
  end subroutine define_single_tile

  subroutine single_tile_scalar_test(domain)
    type(domain2d), intent(inout) :: domain
    real(r4_kind) :: global((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo))
    real(r4_kind) :: x(isd:ied, jsd:jed)
    integer :: i, j, h

    do concurrent (i=lbound(global,1):ubound(global,1), &
                   j=lbound(global,2):ubound(global,2))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n) h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n) h = 999000._r4_kind
      global(i,j) = i*1._r4_kind + j*.01_r4_kind + h
    enddo

    x(isd:ied, jsd:jed) = global(isd:ied, jsd:jed)
    call mpp_update_domains(x, domain)

    if (any(abs(x(isd:ied, jsd:jed) - global(isd:ied, jsd:jed)) .gt. eps)) then
      call mpp_error(FATAL, "Error: A-grid scalar update (single tile) did not produce expected result")
    else
      print "(A)", "A-grid scalar update (single tile) succeeded on PE " // string(mpp_pe())
    endif
  end subroutine single_tile_scalar_test

  subroutine define_corner_mosaic
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: tile_id
    integer, allocatable :: pelist(:)
    type(contact) :: contacts(num_contacts)
    integer :: npes, i

    npes = dot_product(ni, nj)

    if (mpp_npes().lt.npes) then
      call mpp_error(FATAL, "Cannot construct corner mosaic: Insufficient number of PEs")
    endif

    if (mpp_pe().gt.npes - 1) then
      return
    endif

    allocate (pelist(npes))
    do i=1,npes
      pelist(i) = i-1
    enddo
    call mpp_declare_pelist(pelist)
    call mpp_set_current_pelist(pelist)

    isg=1
    ieg=n
    jsg=1
    jeg=n

    contacts(:) = [ &
      contact( &
        tile = [1, 2], &
        is = [ieg, isg], &
        ie = [ieg, isg], &
        js = [jsg, jsg], &
        je = [jeg, jeg]), &
      contact( &
        tile = [2, 3], &
        is = [isg, isg], &
        ie = [ieg, ieg], &
        js = [jeg, jsg], &
        je = [jeg, jsg]), &
      contact( &
        tile = [1, 3], &
        is = [isg, isg], &
        ie = [ieg, isg], &
        js = [jeg, jeg], &
        je = [jeg, jsg]) &
    ]

    do i=1,ntiles
      global_indices(:, i) = [isg, ieg, jsg, jeg]
      layout(:, i) = [ni(i), nj(i)]

      pe_start(i) = dot_product(ni(1:i-1), nj(1:i-1))
      pe_end(i) = dot_product(ni(1:i), nj(1:i)) - 1
    enddo

    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contacts, &
        contacts(:)%tile(1), contacts(:)%tile(2), &
        contacts(:)%is(1), contacts(:)%ie(1), contacts(:)%js(1), contacts(:)%je(1), &
        contacts(:)%is(2), contacts(:)%ie(2), contacts(:)%js(2), contacts(:)%je(2), &
        pe_start, pe_end, symmetry=.true., whalo=halo, ehalo=halo, &
        shalo=halo, nhalo=halo, name=name)

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    tile_id = mpp_get_tile_id(domain)
    if (size(tile_id).ne.1) then
      call mpp_error(FATAL, "Too many tiles on current PE: Only one expected")
    endif
    tid = tile_id(1)

    call test_agrid_scalar_update(domain)
    call test_cgrid_vector_update(domain)
    call test_dgrid_vector_update(domain)

    call mpp_deallocate_domain(domain)
  end subroutine define_corner_mosaic

  subroutine agrid_scalar_init(global, x)
    real(r4_kind), intent(out) :: global((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo), ntiles)
    real(r4_kind), intent(out) :: x(isd:ied, jsd:jed)
    integer :: i, j, t, h

    do concurrent (i=lbound(global,1):ubound(global,1), &
                   j=lbound(global,2):ubound(global,2), &
                   t=lbound(global,3):ubound(global,3))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n) h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n) h = 999000._r4_kind
      global(i,j,t) = t*100._r4_kind + i*1._r4_kind + j*.01_r4_kind + h
    enddo

    x(isd:ied, jsd:jed) = global(isd:ied, jsd:jed, tid)

    ! 1-2 !
    global((n+1):(n+halo), 1:n, 1) = global(1:halo,       1:n, 2)
    global((1-halo):0,     1:n, 2) = global((n+1-halo):n, 1:n, 1)

    ! 2-3 !
    global(1:n, (n+1):(n+halo), 2) = global(1:n, 1:halo,       3)
    global(1:n, (1-halo):0,     3) = global(1:n, (n+1-halo):n, 2)

    ! 1-3 !
    do concurrent (i=1:n, j=1:halo)
      global(i,      n+j, 1) = global(j,     n-i+1,    3)
      global(j-halo, i,   3) = global(n-i+1, n-halo+j, 1)
    enddo
  end subroutine agrid_scalar_init

  !!!!!!!!!
  !   --- !
  !  /|3| !
  ! ----- !
  ! |1|2| !
  ! ----- !
  !!!!!!!!!

  subroutine cgrid_vector_init(gx, gy, x, y)
    real(r4_kind) :: gx((isg-halo):(ieg+halo+1), (jsg-halo):(jeg+halo), ntiles), &
                     gy((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo+1), ntiles)
    real(r4_kind) :: x(isd:ied+1, jsd:jed), &
                     y(isd:ied, jsd:jed+1)
    integer :: i, j, t, h, repeat

    do concurrent (i=lbound(gx,1):ubound(gx,1), &
                   j=lbound(gx,2):ubound(gx,2), &
                   t=lbound(gx,3):ubound(gx,3))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n+1) h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n)   h = 999000._r4_kind
      gx(i,j,t) = t*100._r4_kind + i*1._r4_kind + j*.01_r4_kind + h
    enddo

    do concurrent (i=lbound(gy,1):ubound(gy,1), &
                   j=lbound(gy,2):ubound(gy,2), &
                   t=lbound(gy,3):ubound(gy,3))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n)   h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n+1) h = 999000._r4_kind
      gy(i,j,t) = -1._r4_kind * (t*100._r4_kind + i*1._r4_kind + j*.01_r4_kind + h)
    enddo

    ! TODO: Copying the halos and the common rows/columns twice seems to make everything
    ! consistent. A single round of copying would probably suffice if the copying
    ! operations were more carefully ordered.
    do repeat=1,2
      ! 1,2 (x)
      gx((1-halo):1,       1:n, 2) = gx((n+1-halo):(n+1), 1:n, 1)
      gx((n+1):(n+1+halo), 1:n, 1) = gx(1:(1+halo),       1:n, 2)

      ! 1,2 (y)
      gy((1-halo):0,     1:(n+1), 2) = gy((n+1-halo):n, 1:(n+1), 1)
      gy((n+1):(n+halo), 1:(n+1), 1) = gy(1:halo,       1:(n+1), 2)

      ! 2,3 (x)
      gx(1:(n+1), (1-halo):0,     3) = gx(1:(n+1), (n+1-halo):n, 2)
      gx(1:(n+1), (n+1):(n+halo), 2) = gx(1:(n+1), 1:halo,       3)

      ! 2,3 (y)
      gy(1:n, (1-halo):1,       3) = gy(1:n, (n+1-halo):(n+1), 2)
      gy(1:n, (n+1):(n+1+halo), 2) = gy(1:n, 1:(1+halo),       3)

      ! 3(x), 1(y)
      gy(1:n,        (n+1):(n+1+halo), 1) = transpose(gx(1:(1+halo), n:1:-1,           3))
      gx((1-halo):1, 1:n,              3) = transpose(gy(n:1:-1,     (n+1-halo):(n+1), 1))

      ! 3(y), 1(x)
      gx(1:(n+1),    (n+1):(n+halo), 1) = -transpose(gy(1:halo,     (n+1):1:-1,   3))
      gy((1-halo):0, 1:(n+1),        3) = -transpose(gx((n+1):1:-1, (n+1-halo):n, 1))
    enddo

    x(isd:ied+1, jsd:jed)   = gx(isd:ied+1, jsd:jed,   tid)
    y(isd:ied,   jsd:jed+1) = gy(isd:ied,   jsd:jed+1, tid)
  end subroutine cgrid_vector_init

  subroutine dgrid_vector_init(gx, gy, x, y)
    real(r4_kind) :: gx((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo+1), ntiles), &
                     gy((isg-halo):(ieg+halo+1), (jsg-halo):(jeg+halo), ntiles)
    real(r4_kind) :: x(isd:ied, jsd:jed+1), &
                     y(isd:ied+1, jsd:jed)
    integer :: i, j, t, h, repeat

    do concurrent (i=lbound(gx,1):ubound(gx,1), &
                   j=lbound(gx,2):ubound(gx,2), &
                   t=lbound(gx,3):ubound(gx,3))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n) h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n+1)   h = 999000._r4_kind
      gx(i,j,t) = t*100._r4_kind + i*1._r4_kind + j*.01_r4_kind + h
    enddo

    do concurrent (i=lbound(gy,1):ubound(gy,1), &
                   j=lbound(gy,2):ubound(gy,2), &
                   t=lbound(gy,3):ubound(gy,3))
      h = 0._r4_kind
      if (i.lt.1 .or. i.gt.n+1)   h = 999000._r4_kind
      if (j.lt.1 .or. j.gt.n) h = 999000._r4_kind
      gy(i,j,t) = -1._r4_kind * (t*100._r4_kind + i*1._r4_kind + j*.01_r4_kind + h)
    enddo

    ! TODO: Copying the halos and the common rows/columns twice seems to make everything
    ! consistent. A single round of copying would probably suffice if the copying
    ! operations were more carefully ordered.
    do repeat=1,2
      ! 1,2 (x)
      gy((1-halo):1,       1:n, 2) = gy((n+1-halo):(n+1), 1:n, 1)
      gy((n+1):(n+1+halo), 1:n, 1) = gy(1:(1+halo),       1:n, 2)

      ! 1,2 (y)
      gx((1-halo):0,     1:(n+1), 2) = gx((n+1-halo):n, 1:(n+1), 1)
      gx((n+1):(n+halo), 1:(n+1), 1) = gx(1:halo,       1:(n+1), 2)

      ! 2,3 (x)
      gy(1:(n+1), (1-halo):0,     3) = gy(1:(n+1), (n+1-halo):n, 2)
      gy(1:(n+1), (n+1):(n+halo), 2) = gy(1:(n+1), 1:halo,       3)

      ! 2,3 (y)
      gx(1:n, (1-halo):1,       3) = gx(1:n, (n+1-halo):(n+1), 2)
      gx(1:n, (n+1):(n+1+halo), 2) = gx(1:n, 1:(1+halo),       3)

      ! 3(x), 1(y)
      gx(1:n,        (n+1):(n+1+halo), 1) = -transpose(gy(1:(1+halo), n:1:-1,           3))
      gy((1-halo):1, 1:n,              3) = -transpose(gx(n:1:-1,     (n+1-halo):(n+1), 1))

      ! 3(y), 1(x)
      gy(1:(n+1),    (n+1):(n+halo), 1) = transpose(gx(1:halo,     (n+1):1:-1,   3))
      gx((1-halo):0, 1:(n+1),        3) = transpose(gy((n+1):1:-1, (n+1-halo):n, 1))
    enddo

    x(isd:ied,   jsd:jed+1) = gx(isd:ied,   jsd:jed+1, tid)
    y(isd:ied+1, jsd:jed)   = gy(isd:ied+1, jsd:jed,   tid)
  end subroutine dgrid_vector_init

  subroutine test_agrid_scalar_update(domain)
    type(domain2d), intent(inout) :: domain
    real(r4_kind) :: global((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo), ntiles)
    real(r4_kind) :: x(isd:ied, jsd:jed)

    call agrid_scalar_init(global, x)
    call mpp_update_domains(x, domain)

    if (any(abs(x(isd:ied, jsd:jed) - global(isd:ied, jsd:jed, tid)) .gt. eps)) then
      call mpp_error(FATAL, "Error: A-grid scalar update did not produce expected result")
    else
      print "(A)", "A-grid scalar update succeeded on PE " // string(mpp_pe())
    endif
  end subroutine test_agrid_scalar_update

  subroutine test_cgrid_vector_update(domain)
    type(domain2d), intent(inout) :: domain
    real(r4_kind) :: gx((isg-halo):(ieg+halo+1), (jsg-halo):(jeg+halo), ntiles), &
                     gy((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo+1), ntiles)
    real(r4_kind) :: x(isd:ied+1, jsd:jed), &
                     y(isd:ied, jsd:jed+1)

    call cgrid_vector_init(gx, gy, x, y)
    call mpp_update_domains(x, y, domain, gridtype=CGRID_NE)

    if (any(abs(x(isd:ied+1, jsd:jed) - gx(isd:ied+1, jsd:jed, tid)) .gt. eps) .or. &
        any(abs(y(isd:ied, jsd:jed+1) - gy(isd:ied, jsd:jed+1, tid)) .gt. eps)) then
      call mpp_error(FATAL, "Error: C-grid vector update did not produce expected result")
    else
      print "(A)", "C-grid vector update succeeded on PE " // string(mpp_pe())
    endif
  end subroutine test_cgrid_vector_update

  subroutine test_dgrid_vector_update(domain)
    type(domain2d), intent(inout) :: domain
    real(r4_kind) :: gx((isg-halo):(ieg+halo), (jsg-halo):(jeg+halo+1), ntiles), &
                     gy((isg-halo):(ieg+halo+1), (jsg-halo):(jeg+halo), ntiles)
    real(r4_kind) :: x(isd:ied, jsd:jed+1), &
                     y(isd:ied+1, jsd:jed)

    call dgrid_vector_init(gx, gy, x, y)
    call mpp_update_domains(x, y, domain, gridtype=DGRID_NE)

    if (any(abs(x(isd:ied, jsd:jed+1) - gx(isd:ied, jsd:jed+1, tid)) .gt. eps) .or. &
        any(abs(y(isd:ied+1, jsd:jed) - gy(isd:ied+1, jsd:jed, tid)) .gt. eps)) then
      call mpp_error(FATAL, "Error: D-grid vector update did not produce expected result")
    else
      print "(A)", "D-grid vector update succeeded on PE " // string(mpp_pe())
    endif
  end subroutine test_dgrid_vector_update

  subroutine print_matrix(x)
    real(r4_kind), dimension(:,:), intent(in) :: x
    integer :: i,j
    integer :: iostat

    do j=size(x,2),1,-1
      do i=1,size(x,1)
        write(*,'(" ",F12.4)', advance="no", iostat=iostat) x(i,j)
      enddo
      write(*, "(A)") ""
    enddo
  end subroutine
end program test_corner_mosaic
