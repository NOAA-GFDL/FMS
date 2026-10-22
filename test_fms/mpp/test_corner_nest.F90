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

program test_nest_corner
  use platform_mod, only: r8_kind
  use test_domains_utility_mod, only: define_cubic_mosaic
  use mpp_mod
  use mpp_domains_mod
  use fms_mod
  use fms_string_utils_mod, only: string
  use answers_mod

  implicit none

  integer, parameter :: stackmax=10000000

  integer, parameter :: nc=20 !< 20x20 coarse grid
  integer, parameter :: nf=10 !< 10x10 fine grid. This corresponds to a 5x5 region on the coarse grid.
  integer, parameter :: nz=1
  integer, parameter :: halo=2

  integer, parameter :: ntiles_coarse=6
  integer, parameter :: ntiles_fine=3
  integer, parameter :: ntiles_all = ntiles_coarse+ntiles_fine

  integer, parameter :: layout_coarse(2) = [2,2] !< 2x2 PE layout for the coarse grid
  integer, parameter :: layout_fine(2) = [2,2] !< 2x2 PE layout for the fine grid

  type(domain2d)         :: domain !< Cubed sphere (coarse PEs) or single-tile fine-grid (fine PEs) domain
  type(domain2d)         :: domain_corner !< Domain representing the corner mosaic (only used on fine PEs)
  type(nest_domain_type) :: nest_domain !< Nest domain

  ! Nest parameters
  integer, parameter :: num_nests = 3 !< The corner nest system requires three nests
  integer, parameter :: refine_ratio = 2 !< 2:1 ratio of fine points to coarse points
  integer, parameter :: nest_level(3) = [1, 1, 1] !< All three nests are level 1 (directly on the coarse grid)
  integer, parameter :: tile_coarse(num_nests) = [1, 2, 3] !< Tile IDs of the coarse tiles
  integer, parameter :: tile_fine(num_nests)   = [7, 8, 9] !< Tile IDs of the corresponding fine tiles

  integer, parameter :: istart_coarse(num_nests) = [16,  1, 1]
  integer, parameter :: jstart_coarse(num_nests) = [16, 16, 1]
  integer, parameter :: icount_coarse(num_nests) = [ 5,  5, 5]
  integer, parameter :: jcount_coarse(num_nests) = [ 5,  5, 5]

  integer :: i
  integer, parameter :: npes(ntiles_coarse+ntiles_fine) = [(product(layout_coarse), i=1,ntiles_coarse), &
                                                           (product(layout_fine),   i=1,ntiles_fine)]

  integer, parameter :: npes_coarse = sum(npes(1:ntiles_coarse))
  integer, parameter :: npes_fine = sum(npes(ntiles_coarse+1:ntiles_all))
  integer, parameter :: npes_all = npes_coarse + npes_fine

  ! PE lists
  integer, target :: pelist_master(npes_all)
  integer, pointer :: pelist_coarse(:)
  integer, pointer :: pelist_fine(:)

  logical :: is_coarse, is_fine
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  call fms_init
  call mpp_domains_set_stack_size(stackmax)

  if (mpp_npes().ne.npes_all) call mpp_error(FATAL, "Test requires "//string(sum(npes))//" PEs: "// &
                                                     string(npes_coarse)//" for the coarse grid and "// &
                                                     string(npes_fine)// " for the fine grid")

  call mpp_get_current_pelist(pelist_master)
  pelist_coarse => pelist_master(1:npes_coarse)
  pelist_fine => pelist_master(npes_coarse+1:npes_all)

  call mpp_declare_pelist(pelist_coarse)
  call mpp_declare_pelist(pelist_fine)

  is_coarse = any(mpp_pe().eq.pelist_coarse)
  is_fine = any(mpp_pe().eq.pelist_fine)

  if (is_coarse.eqv.is_fine) then
    call mpp_error(FATAL, "PE "//string(mpp_pe())//" must be either coarse or fine (but not both)")
  endif

  call answers_mod_init(25)

  if (is_coarse) then
    call mpp_set_current_pelist(pelist_coarse)
    call define_cubic_mosaic_wrap(domain)
  endif

  if (is_fine) then
    call mpp_set_current_pelist(pelist_fine)
    call define_corner_mosaic(domain_corner)

    block
      integer :: pelist_this_tile(product(layout_fine))

      call mpp_get_tile_pelist(domain_corner, pelist_this_tile)
      call mpp_declare_pelist(pelist_this_tile)
      call mpp_set_current_pelist(pelist_this_tile)
    end block

    call mpp_define_domains([1,nf,1,nf], layout_fine, domain, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
                            symmetry=.true., name="Fine single tile", tile_id=get_tid(domain_corner))
  endif

  call mpp_set_current_pelist
  call mpp_define_nest_domains(nest_domain, domain, num_nests, nest_level, tile_fine, tile_coarse, &
       istart_coarse, icount_coarse, jstart_coarse, jcount_coarse, npes, &
       [refine_ratio, refine_ratio, refine_ratio], [refine_ratio, refine_ratio, refine_ratio], &
       name="corner_nest")

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(domain, isd, ied, jsd, jed)

  call update_cycle_scalar
  call update_cycle_vector_cgrid
  call update_cycle_vector_dgrid

  call answers_mod_end
  call fms_end

contains

  subroutine update_cycle_scalar
    real(r8_kind), allocatable, dimension(:,:,:) :: x
    integer :: i, j, k

    allocate(x(isd:ied, jsd:jed, nz))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Coarse to fine update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: wbuffer, ebuffer, sbuffer, nbuffer
      integer :: isw_c, iew_c, jsw_c, jew_c
      integer :: isw_f, iew_f, jsw_f, jew_f
      integer :: ise_c, iee_c, jse_c, jee_c
      integer :: ise_f, iee_f, jse_f, jee_f
      integer :: iss_c, ies_c, jss_c, jes_c
      integer :: iss_f, ies_f, jss_f, jes_f
      integer :: isn_c, ien_c, jsn_c, jen_c
      integer :: isn_f, ien_f, jsn_f, jen_f

      if (is_coarse) then
        x = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        call mpp_set_current_pelist(pelist_coarse)
        call mpp_update_domains(x, domain)
        call register_answer(x, "Initialized coarse-grid data (scalar)")

        allocate(wbuffer(1,1,1))
        allocate(ebuffer(1,1,1))
        allocate(sbuffer(1,1,1))
        allocate(nbuffer(1,1,1))
      endif

      if (is_fine) then
        x = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, WEST, 1)
        if (iew_c.ge.isw_c.and.jew_c.ge.jsw_c) then
          allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c, nz))
        else
          allocate(wbuffer(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, EAST, 1)
        if (iee_c.ge.ise_c.and.jee_c.ge.jse_c) then
          allocate(ebuffer(ise_c:iee_c, jse_c:jee_c, nz))
        else
          allocate(ebuffer(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, SOUTH, 1)
        if (ies_c.ge.iss_c.and.jes_c.ge.jss_c) then
          allocate(sbuffer(iss_c:ies_c, jss_c:jes_c, nz))
        else
          allocate(sbuffer(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, NORTH, 1)
        if (ien_c.ge.isn_c.and.jen_c.ge.jsn_c) then
          allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c, nz))
        else
          allocate(nbuffer(1,1,1))
        endif

        wbuffer=0._r8_kind
        ebuffer=0._r8_kind
        sbuffer=0._r8_kind
        nbuffer=0._r8_kind
      endif

      call mpp_set_current_pelist
      call mpp_update_nest_fine(x, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=1)

      if (is_fine) then
        if (iew_c.ge.isw_c.and.jew_c.ge.jsw_c) then
          call register_answer(wbuffer, "C2F West buffer (scalar)")
        endif

        if (iee_c.ge.ise_c.and.jee_c.ge.jse_c) then
          call register_answer(ebuffer, "C2F East buffer (scalar)")
        endif

        if (ies_c.ge.iss_c.and.jes_c.ge.jss_c) then
          call register_answer(sbuffer, "C2F South buffer (scalar)")
        endif

        if (ien_c.ge.isn_c.and.jen_c.ge.jsn_c) then
          call register_answer(nbuffer, "C2F North buffer (scalar)")
        endif
      endif
    end block

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Halo update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (is_fine) then
      associate (xc => x(isc:iec, jsc:jec, :))
        xc = xc + .000000001_r8_kind
      end associate

      call mpp_set_current_pelist(pelist_fine)
      call mpp_update_domains(x, domain_corner)
      call register_answer(x, "Fine data after domain update (scalar)")
    endif

    call mpp_set_current_pelist

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Fine to coarse update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: x_f2c
      integer :: is_c, ie_c, js_c, je_c
      integer :: is_f, ie_f, js_f, je_f
      integer :: if0, if1, jf0, jf1

      if (is_fine) then
        call mpp_get_f2c_index(nest_domain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, nest_level=1)
        allocate(x_f2c(is_c:ie_c, js_c:je_c, nz))

        do j=js_c,je_c
          jf0 = js_f + refine_ratio*(j-js_c)
          jf1 = js_f + refine_ratio*(j-js_c + 1) - 1

          do i=is_c, ie_c
            if0 = is_f + refine_ratio*(i-is_c)
            if1 = is_f + refine_ratio*(i-is_c + 1) - 1

            do k=1,nz
              x_f2c(i,j,k) = sum(x(if0:if1, jf0:jf1, k)) / refine_ratio**2
            enddo
          enddo
        enddo
      endif

      if (is_coarse) allocate(x_f2c(1,1,1))

      call mpp_update_nest_coarse(x_f2c, nest_domain, x, nest_level=1)

      if (is_coarse) then
        call register_answer(x, "F2C update result (scalar)")
      endif
    end block
  end subroutine update_cycle_scalar

  subroutine update_cycle_vector_cgrid
    real(r8_kind), allocatable, dimension(:,:,:) :: x, y
    integer :: i, j, k

    allocate(x(isd:ied+1, jsd:jed,   nz))
    allocate(y(isd:ied,   jsd:jed+1, nz))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Coarse to fine update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: wbufferx, ebufferx, sbufferx, nbufferx
      real(r8_kind), allocatable, dimension(:,:,:) :: wbuffery, ebuffery, sbuffery, nbuffery

      integer :: isw_cx, iew_cx, jsw_cx, jew_cx
      integer :: isw_fx, iew_fx, jsw_fx, jew_fx
      integer :: ise_cx, iee_cx, jse_cx, jee_cx
      integer :: ise_fx, iee_fx, jse_fx, jee_fx
      integer :: iss_cx, ies_cx, jss_cx, jes_cx
      integer :: iss_fx, ies_fx, jss_fx, jes_fx
      integer :: isn_cx, ien_cx, jsn_cx, jen_cx
      integer :: isn_fx, ien_fx, jsn_fx, jen_fx

      integer :: isw_cy, iew_cy, jsw_cy, jew_cy
      integer :: isw_fy, iew_fy, jsw_fy, jew_fy
      integer :: ise_cy, iee_cy, jse_cy, jee_cy
      integer :: ise_fy, iee_fy, jse_fy, jee_fy
      integer :: iss_cy, ies_cy, jss_cy, jes_cy
      integer :: iss_fy, ies_fy, jss_fy, jes_fy
      integer :: isn_cy, ien_cy, jsn_cy, jen_cy
      integer :: isn_fy, ien_fy, jsn_fy, jen_fy

      if (is_coarse) then
        x = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec+1
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        y = 0._r8_kind
        do j = jsc, jec+1
          do i = isc, iec
            y(i,j,:) = -(get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind)
          enddo
        enddo

        call mpp_set_current_pelist(pelist_coarse)
        call mpp_update_domains(x, y, domain, gridtype=CGRID_NE)
        call register_answer(x, "Initialized coarse-grid data (cgrid, x)")
        call register_answer(y, "Initialized coarse-grid data (cgrid, y)")

        allocate(wbufferx(1,1,1), wbuffery(1,1,1))
        allocate(ebufferx(1,1,1), ebuffery(1,1,1))
        allocate(sbufferx(1,1,1), sbuffery(1,1,1))
        allocate(nbufferx(1,1,1), nbuffery(1,1,1))
      endif

      if (is_fine) then
        x = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec+1
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        y = 0._r8_kind
        do j = jsc, jec+1
          do i = isc, iec
            y(i,j,:) = -(get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind)
          enddo
        enddo

        call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, 1, &
                               position=EAST)
        call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, 1, &
                               position=NORTH)
        if (iew_cx.ge.isw_cx.and.jew_cx.ge.jsw_cx) then
          allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx, nz))
          allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy, nz))
        else
          allocate(wbufferx(1,1,1))
          allocate(wbuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, 1, &
                               position=EAST)
        call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, 1, &
                               position=NORTH)
        if (iee_cx.ge.ise_cx.and.jee_cx.ge.jse_cx) then
          allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx, nz))
          allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy, nz))
        else
          allocate(ebufferx(1,1,1))
          allocate(ebuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, 1, &
                               position=EAST)
        call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, 1, &
                               position=NORTH)
        if (ies_cx.ge.iss_cx.and.jes_cx.ge.jss_cx) then
          allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx, nz))
          allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy, nz))
        else
          allocate(sbufferx(1,1,1))
          allocate(sbuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, 1, &
                               position=EAST)
        call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, 1, &
                               position=NORTH)
        if (ien_cx.ge.isn_cx.and.jen_cx.ge.jsn_cx) then
          allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx, nz))
          allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy, nz))
        else
          allocate(nbufferx(1,1,1))
          allocate(nbuffery(1,1,1))
        endif

        wbufferx=0._r8_kind
        ebufferx=0._r8_kind
        sbufferx=0._r8_kind
        nbufferx=0._r8_kind

        wbuffery=0._r8_kind
        ebuffery=0._r8_kind
        sbuffery=0._r8_kind
        nbuffery=0._r8_kind
      endif

      call mpp_set_current_pelist
      call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
                                nbufferx, nbuffery, nest_level=1, gridtype=CGRID_NE)

      if (is_fine) then
        if (iew_cx.ge.isw_cx.and.jew_cx.ge.jsw_cx) then
          call register_answer(wbufferx, "C2F West buffer (cgrid, x)")
          call register_answer(wbuffery, "C2F West buffer (cgrid, y)")
        endif

        if (iee_cx.ge.ise_cx.and.jee_cx.ge.jse_cx) then
          call register_answer(ebufferx, "C2F East buffer (cgrid, x)")
          call register_answer(ebuffery, "C2F East buffer (cgrid, y)")
        endif

        if (ies_cx.ge.iss_cx.and.jes_cx.ge.jss_cx) then
          call register_answer(sbufferx, "C2F South buffer (cgrid, x)")
          call register_answer(sbuffery, "C2F South buffer (cgrid, y)")
        endif

        if (ien_cx.ge.isn_cx.and.jen_cx.ge.jsn_cx) then
          call register_answer(nbufferx, "C2F North buffer (cgrid, x)")
          call register_answer(nbuffery, "C2F North buffer (cgrid, y)")
        endif
      endif
    end block

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Halo update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (is_fine) then
      associate (xc => x(isc:iec+1, jsc:jec, :))
        xc = xc + .000000001_r8_kind
      end associate

      associate (yc => y(isc:iec, jsc:jec+1, :))
        yc = yc - .000000001_r8_kind
      end associate

      call mpp_set_current_pelist(pelist_fine)
      call mpp_update_domains(x, y, domain_corner, gridtype=CGRID_NE)
      call register_answer(x, "Fine data after domain update (cgrid, x)")
      call register_answer(y, "Fine data after domain update (cgrid, y)")
    endif

    call mpp_set_current_pelist

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Fine to coarse update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: x_f2c, y_f2c
      integer :: is_cx, ie_cx, js_cx, je_cx
      integer :: is_fx, ie_fx, js_fx, je_fx
      integer :: is_cy, ie_cy, js_cy, je_cy
      integer :: is_fy, ie_fy, js_fy, je_fy
      integer :: if0, if1, jf0, jf1

      if (is_fine) then
        call mpp_get_f2c_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, nest_level=1, &
                               position=EAST)
        call mpp_get_f2c_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, nest_level=1, &
                               position=NORTH)
        allocate(x_f2c(is_cx:ie_cx, js_cx:je_cx, nz))
        allocate(y_f2c(is_cy:ie_cy, js_cy:je_cy, nz))

        do j=js_cx,je_cx
          jf0 = js_fx + refine_ratio*(j-js_cx)
          jf1 = js_fx + refine_ratio*(j-js_cx + 1) - 1

          do i=is_cx, ie_cx
            if0 = is_fx + refine_ratio*(i-is_cx)
            if1 = is_fx + refine_ratio*(i-is_cx + 1) - 1

            do k=1,nz
              x_f2c(i,j,k) = sum(x(if0:if1, jf0:jf1, k)) / refine_ratio**2
            enddo
          enddo
        enddo

        do j=js_cy,je_cy
          jf0 = js_fy + refine_ratio*(j-js_cy)
          jf1 = js_fy + refine_ratio*(j-js_cy + 1) - 1

          do i=is_cy, ie_cy
            if0 = is_fy + refine_ratio*(i-is_cy)
            if1 = is_fy + refine_ratio*(i-is_cy + 1) - 1

            do k=1,nz
              y_f2c(i,j,k) = sum(y(if0:if1, jf0:jf1, k)) / refine_ratio**2
            enddo
          enddo
        enddo
      endif

      if (is_coarse) allocate(x_f2c(1,1,1), y_f2c(1,1,1))

      call mpp_update_nest_coarse(x_f2c, y_f2c, nest_domain, x, y, nest_level=1, gridtype=CGRID_NE)

      if (is_coarse) then
        call register_answer(x, "F2C update result (cgrid, x)")
        call register_answer(y, "F2C update result (cgrid, y)")
      endif
    end block
  end subroutine update_cycle_vector_cgrid

  subroutine update_cycle_vector_dgrid
    real(r8_kind), allocatable, dimension(:,:,:) :: x, y
    integer :: i, j, k

    allocate(x(isd:ied,   jsd:jed+1, nz))
    allocate(y(isd:ied+1, jsd:jed,   nz))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Coarse to fine update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: wbufferx, ebufferx, sbufferx, nbufferx
      real(r8_kind), allocatable, dimension(:,:,:) :: wbuffery, ebuffery, sbuffery, nbuffery

      integer :: isw_cx, iew_cx, jsw_cx, jew_cx
      integer :: isw_fx, iew_fx, jsw_fx, jew_fx
      integer :: ise_cx, iee_cx, jse_cx, jee_cx
      integer :: ise_fx, iee_fx, jse_fx, jee_fx
      integer :: iss_cx, ies_cx, jss_cx, jes_cx
      integer :: iss_fx, ies_fx, jss_fx, jes_fx
      integer :: isn_cx, ien_cx, jsn_cx, jen_cx
      integer :: isn_fx, ien_fx, jsn_fx, jen_fx

      integer :: isw_cy, iew_cy, jsw_cy, jew_cy
      integer :: isw_fy, iew_fy, jsw_fy, jew_fy
      integer :: ise_cy, iee_cy, jse_cy, jee_cy
      integer :: ise_fy, iee_fy, jse_fy, jee_fy
      integer :: iss_cy, ies_cy, jss_cy, jes_cy
      integer :: iss_fy, ies_fy, jss_fy, jes_fy
      integer :: isn_cy, ien_cy, jsn_cy, jen_cy
      integer :: isn_fy, ien_fy, jsn_fy, jen_fy

      if (is_coarse) then
        x = 0._r8_kind
        do j = jsc, jec+1
          do i = isc, iec
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        y = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec+1
            y(i,j,:) = -(get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind)
          enddo
        enddo

        call mpp_set_current_pelist(pelist_coarse)
        call mpp_update_domains(x, y, domain, gridtype=DGRID_NE)
        call register_answer(x, "Initialized coarse-grid data (dgrid, x)")
        call register_answer(y, "Initialized coarse-grid data (dgrid, y)")

        allocate(wbufferx(1,1,1), wbuffery(1,1,1))
        allocate(ebufferx(1,1,1), ebuffery(1,1,1))
        allocate(sbufferx(1,1,1), sbuffery(1,1,1))
        allocate(nbufferx(1,1,1), nbuffery(1,1,1))
      endif

      if (is_fine) then
        x = 0._r8_kind
        do j = jsc, jec+1
          do i = isc, iec
            x(i,j,:) = get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind
          enddo
        enddo

        y = 0._r8_kind
        do j = jsc, jec
          do i = isc, iec+1
            y(i,j,:) = -(get_tid(domain) + i*.01_r8_kind + j*.00001_r8_kind)
          enddo
        enddo

        call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, 1, &
                               position=NORTH)
        call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, 1, &
                               position=EAST)
        if (iew_cx.ge.isw_cx.and.jew_cx.ge.jsw_cx) then
          allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx, nz))
          allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy, nz))
        else
          allocate(wbufferx(1,1,1))
          allocate(wbuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, 1, &
                               position=NORTH)
        call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, 1, &
                               position=EAST)
        if (iee_cx.ge.ise_cx.and.jee_cx.ge.jse_cx) then
          allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx, nz))
          allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy, nz))
        else
          allocate(ebufferx(1,1,1))
          allocate(ebuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, 1, &
                               position=NORTH)
        call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, 1, &
                               position=EAST)
        if (ies_cx.ge.iss_cx.and.jes_cx.ge.jss_cx) then
          allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx, nz))
          allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy, nz))
        else
          allocate(sbufferx(1,1,1))
          allocate(sbuffery(1,1,1))
        endif

        call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, 1, &
                               position=NORTH)
        call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, 1, &
                               position=EAST)
        if (ien_cx.ge.isn_cx.and.jen_cx.ge.jsn_cx) then
          allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx, nz))
          allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy, nz))
        else
          allocate(nbufferx(1,1,1))
          allocate(nbuffery(1,1,1))
        endif

        wbufferx=0._r8_kind
        ebufferx=0._r8_kind
        sbufferx=0._r8_kind
        nbufferx=0._r8_kind

        wbuffery=0._r8_kind
        ebuffery=0._r8_kind
        sbuffery=0._r8_kind
        nbuffery=0._r8_kind
      endif

      call mpp_set_current_pelist
      call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
                                nbufferx, nbuffery, nest_level=1, gridtype=DGRID_NE)

      if (is_fine) then
        if (iew_cx.ge.isw_cx.and.jew_cx.ge.jsw_cx) then
          call register_answer(wbufferx, "C2F West buffer (dgrid, x)")
          call register_answer(wbuffery, "C2F West buffer (dgrid, y)")
        endif

        if (iee_cx.ge.ise_cx.and.jee_cx.ge.jse_cx) then
          call register_answer(ebufferx, "C2F East buffer (dgrid, x)")
          call register_answer(ebuffery, "C2F East buffer (dgrid, y)")
        endif

        if (ies_cx.ge.iss_cx.and.jes_cx.ge.jss_cx) then
          call register_answer(sbufferx, "C2F South buffer (dgrid, x)")
          call register_answer(sbuffery, "C2F South buffer (dgrid, y)")
        endif

        if (ien_cx.ge.isn_cx.and.jen_cx.ge.jsn_cx) then
          call register_answer(nbufferx, "C2F North buffer (dgrid, x)")
          call register_answer(nbuffery, "C2F North buffer (dgrid, y)")
        endif
      endif
    end block

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Halo update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (is_fine) then
      associate (xc => x(isc:iec, jsc:jec+1, :))
        xc = xc + .000000001_r8_kind
      end associate

      associate (yc => y(isc:iec+1, jsc:jec, :))
        yc = yc - .000000001_r8_kind
      end associate

      call mpp_set_current_pelist(pelist_fine)
      call mpp_update_domains(x, y, domain_corner, gridtype=DGRID_NE)
      call register_answer(x, "Fine data after domain update (dgrid, x)")
      call register_answer(y, "Fine data after domain update (dgrid, y)")
    endif

    call mpp_set_current_pelist

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Fine to coarse update !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
      real(r8_kind), allocatable, dimension(:,:,:) :: x_f2c, y_f2c
      integer :: is_cx, ie_cx, js_cx, je_cx
      integer :: is_fx, ie_fx, js_fx, je_fx
      integer :: is_cy, ie_cy, js_cy, je_cy
      integer :: is_fy, ie_fy, js_fy, je_fy
      integer :: if0, if1, jf0, jf1

      if (is_fine) then
        call mpp_get_f2c_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, nest_level=1, &
                               position=NORTH)
        call mpp_get_f2c_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, nest_level=1, &
                               position=EAST)
        allocate(x_f2c(is_cx:ie_cx, js_cx:je_cx, nz))
        allocate(y_f2c(is_cy:ie_cy, js_cy:je_cy, nz))

        do j=js_cx,je_cx
          jf0 = js_fx + refine_ratio*(j-js_cx)
          jf1 = js_fx + refine_ratio*(j-js_cx + 1) - 1

          do i=is_cx, ie_cx
            if0 = is_fx + refine_ratio*(i-is_cx)
            if1 = is_fx + refine_ratio*(i-is_cx + 1) - 1

            do k=1,nz
              x_f2c(i,j,k) = sum(x(if0:if1, jf0:jf1, k)) / refine_ratio**2
            enddo
          enddo
        enddo

        do j=js_cy,je_cy
          jf0 = js_fy + refine_ratio*(j-js_cy)
          jf1 = js_fy + refine_ratio*(j-js_cy + 1) - 1

          do i=is_cy, ie_cy
            if0 = is_fy + refine_ratio*(i-is_cy)
            if1 = is_fy + refine_ratio*(i-is_cy + 1) - 1

            do k=1,nz
              y_f2c(i,j,k) = sum(y(if0:if1, jf0:jf1, k)) / refine_ratio**2
            enddo
          enddo
        enddo
      endif

      if (is_coarse) allocate(x_f2c(1,1,1), y_f2c(1,1,1))

      call mpp_update_nest_coarse(x_f2c, y_f2c, nest_domain, x, y, nest_level=1, gridtype=DGRID_NE)

      if (is_coarse) then
        call register_answer(x, "F2C update result (dgrid, x)")
        call register_answer(y, "F2C update result (dgrid, y)")
      endif
    end block
  end subroutine update_cycle_vector_dgrid

  subroutine define_cubic_mosaic_wrap(domain)
    type(domain2d), intent(inout) :: domain
    integer, parameter :: nc_arr(ntiles_coarse) = [(nc, i=1,ntiles_coarse)]
    integer :: global_indices(4, ntiles_coarse), layout(2, ntiles_coarse)
    integer :: pe_start(ntiles_coarse), pe_end(ntiles_coarse)
    integer :: i

    do i = 1, ntiles_coarse
      global_indices(:, i) = [1, nc, 1, nc]
      layout(:, i) = layout_coarse
      pe_start(i) = sum(npes(1:i-1))
      pe_end(i)   = sum(npes(1:i)) - 1
    enddo

    call define_cubic_mosaic("Cubic-Grid", domain, nc_arr, nc_arr, global_indices, layout, pe_start, pe_end, &
                             use_memsize=.false., halo=halo)
  end subroutine define_cubic_mosaic_wrap

  subroutine define_corner_mosaic(domain)
    type contact
      integer :: tile(2)
      integer :: is(2)
      integer :: ie(2)
      integer :: js(2)
      integer :: je(2)
    end type contact

    type(domain2d), intent(inout) :: domain
    integer, parameter :: isg=1, ieg=nf, jsg=1, jeg=nf
    integer, parameter :: num_contacts=3
    type(contact) :: contacts(num_contacts)
    integer :: global_indices(4,ntiles_fine), layout(2,ntiles_fine)
    integer :: pe_start(ntiles_fine), pe_end(ntiles_fine)
    integer :: i

    contacts(:) = [ &
      contact( &
        tile = [1, 2]+6, &
        is = [ieg, isg], &
        ie = [ieg, isg], &
        js = [jsg, jsg], &
        je = [jeg, jeg]), &
      contact( &
        tile = [2, 3]+6, &
        is = [isg, isg], &
        ie = [ieg, ieg], &
        js = [jeg, jsg], &
        je = [jeg, jsg]), &
      contact( &
        tile = [1, 3]+6, &
        is = [isg, isg], &
        ie = [ieg, isg], &
        js = [jeg, jeg], &
        je = [jeg, jsg]) &
    ]

    do i=1,ntiles_fine
      global_indices(:, i) = [isg, ieg, jsg, jeg]
      layout(:, i) = layout_fine

      pe_start(i) = sum(npes(1:ntiles_coarse+i-1))
      pe_end(i)   = sum(npes(1:ntiles_coarse+i)) - 1
    enddo

    call mpp_define_mosaic(global_indices, layout, domain, ntiles_fine, num_contacts, &
        contacts(:)%tile(1), contacts(:)%tile(2), &
        contacts(:)%is(1), contacts(:)%ie(1), contacts(:)%js(1), contacts(:)%je(1), &
        contacts(:)%is(2), contacts(:)%ie(2), contacts(:)%js(2), contacts(:)%je(2), &
        pe_start, pe_end, symmetry=.true., whalo=halo, ehalo=halo, &
        shalo=halo, nhalo=halo, name="Corner mosaic", tile_id=tile_fine)
  end subroutine define_corner_mosaic

  function get_tid(domain)
    type(domain2d), intent(in) :: domain
    integer :: get_tid
    integer :: tid(1)

    tid = mpp_get_tile_id(domain)
    get_tid = tid(1)
  end function get_tid

end program test_nest_corner
