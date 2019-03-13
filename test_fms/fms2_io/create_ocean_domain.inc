!> @brief Initialize an ocean style tripolar domain.
subroutine create_ocean_domain(nx, ny, npes, domain, layout, io_layout)

  integer, intent(in) :: nx !< Number of points in the x-direction.
  integer, intent(in) :: ny !< Number of points in the y-direction.
  integer, intent(in) :: npes !< Number of MPI ranks.
  type(domain2d), intent(inout) :: domain !< A tripolar domain.
  integer, dimension(2), intent(in) :: layout
  integer, dimension(2), intent(in) :: io_layout !< I/O layout.

  integer :: whalo
  integer :: ehalo
  integer :: shalo
  integer :: nhalo

  call mpp_domains_init()
! call mpp_define_layout((/1, nx, 1, ny/), npes, layout)
  whalo = 2
  ehalo = whalo
  shalo = whalo
  nhalo = whalo
  call mpp_define_domains((/1,nx,1,ny/), layout, domain, xflags=CYCLIC_GLOBAL_DOMAIN, &
                          yflags=FOLD_NORTH_EDGE, whalo=whalo, ehalo=ehalo, shalo=shalo, &
                          nhalo=nhalo, symmetry=.true., name="Tripolar Folded North Symmetry")
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_ocean_domain
