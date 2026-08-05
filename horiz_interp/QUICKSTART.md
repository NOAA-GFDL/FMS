# Introduction

`horiz_interp_mod` interpolates 2D and 3D data with `conservative order 1`,
`bilinear`, `bicubic`, and `spherical` (inverse-distance weighted) methods.
Not all methods support all input and output grids.  All longititude and 
latitude values must be in radians.

# Basic usage

## 2-step interpolation
1. Call `horiz_interp_init`.
2. Pass in grid to `horiz_interp_new` to compute and store weights in `Interp`.     
3. Pass `Interp` and `data` to `horiz_interp` to interpolate.
4. Free `Interp` with `horiz_interp_del` when done.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp_new, &
                             horiz_interp, horiz_interp_del, horiz_interp_type
use platform_mod, only:  r8_kind ! 64-bit floating point precision

implicit none

type(horiz_interp_type) :: Interp
real(r8_kind), allocatable :: lon_src(:, :), lat_src(:, :) ! source grid
real(r8_kind), allocatable :: lon_dst(:, :), lat_dst(:, :) ! destination grid
real(r8_kind), allocatable :: data_src(:,:), data_dst(:,:) ! data 

! ... allocate and define the source grid in lon_src and lat_src
! ... allocate and define the destination grid in lon_dst and lat_dst
! ... allocate data_src and compute data on source grid

! start program
call fms_init()

! ... set domains and MPI parallelization

! initialize horiz_interp_mod
call horiz_interp_init()

! Compute interpolation weights and mapping indices
call horiz_interp_new(Interp, lon_src, lat_src, lon_dst, lat_dst, interp_method="bilinear")

! Interpolate (can be reused for many data fields on the same grid)
call horiz_interp(Interp, data_src, data_dst)

! Release the interpolator when finished
call horiz_interp_del(Interp)

! finalize horiz_interp_mod
call horiz_interp_end()

! end program
call fms_end()
```

## "Solo" blackbox interpolation
1. Call `horiz_interp_init`.
2. Pass in the data and grid to `horiz_interp`.  The weights will not be saved.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp, horiz_interp_del
use platform_mod, only:  r8_kind ! 64-bit floating point precision

implicit none

real(r8_kind), allocatable :: lon_src(:, :), lat_src(:, :) ! source grid
real(r8_kind), allocatable :: lon_dst(:, :), lat_dst(:, :) ! destination grid
real(r8_kind), allocatable :: data_src(:,:), data_dst(:,:) ! data 

! ... allocate and define the source grid in lon_src and lat_src
! ... allocate and define the destination grid in lon_dst and lat_dst
! ... allocate data_src and compute data on source grid

! start program
call fms_init()

! ... set domains and MPI parallelization

! initialize horiz_interp_mod
call horiz_interp_init()

! Interpolate
call horiz_interp(data_src, lon_src, lat_src, lon_dst, lat_dst, data_dst, interp_method="conservative")

! Finalize horiz_interp_mod
call horiz_interp_end()

! end program
call fms_end()
```

## 1D destination grid example
Rectilinear grids such as a lat/lon grid where the longitude coordinates are 
identical along every line of latitude and the latitude coordinates are identical
along every line of longitude, can be represented as 1D arrays of longitude and 
latitude gridpoints and weight-generating algorithms can be simplified for efficiency.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp_new, &
                             horiz_interp, horiz_interp_del, horiz_interp_type
use platform_mod, only:  r8_kind ! 64-bit floating point precision
use constants_mod, only:  DEG_TO_RAD

implicit none

integer, parameter :: nlon_dst = 360 ! number of grid cells in x-direction
integer, parameter :: nlat_dst = 180 ! number of grid cells in y-direction

type(horiz_interp_type) :: Interp
real(r8_kind), allocatable :: lon_src(:, :), lat_src(:, :) ! source grid
real(r8_kind), allocatable :: lon_dst(:), lat_dst(:) ! destination grid
real(r8_kind), allocatable :: data_src(:,:), data_dst(:,:) ! data 

! ... allocate and define the source grid in lon_src and lat_src

! allocate and define the destination grid in lon_dst and lat_dst
lon_dst = [real(i, r8_kind)*DEG_TO_RAD, i=1, nlon_dst]
lat_dst = [dl*real(i, r8_kind)*DEG_TO_RAD, i=-nlat_dst/2, nlat_dst/2]

! ... allocate data_src and compute data on source grid

! start program
call fms_init()

! ... set domains and MPI parallelization

! initialize horiz_interp_mod
call horiz_interp_init()

! Compute interpolation weights and mapping indices
call horiz_interp_new(Interp, lon_src, lat_src, lon_dst, lat_dst, interp_method="bilinear")

! Interpolate (can be reused for many data fields on the same grid)
call horiz_interp(Interp, data_src, data_dst)

! Release the interpolator when finished
call horiz_interp_del(Interp)

! finalize horiz_interp_mod
call horiz_interp_end()

! end program
call fms_end()
```

# Interpolation methods
It is recommended to set the optional argument `interp_method` to one of the following:

* `"conserve"` — order 1 conservative interpolation.  Supports all combination of 1D and 2D source and destination grids.
* `"bilinear"` — bilinear interpolation.  Only supports 1D source with 2D destination, or 2D source and destination grids.
* `"bicubic"` — bicubic interpolation.  Only supports 1D source with 2D destination, or 1D source and destination grids.
* `"spherical"` — inverse-distance weighting over nearest neighbors.  Only supports 2D source and destination grids.

# Optional arguments
`Horiz_interp_new` supports the following optional arguments for each interpolation methods. 
All methods support `verbose` (sets the verbosity level, 0/1/2, default 0) and `interp_method` (selects the
method itself) are accepted by `horiz_interp_new` regardless of method and are not repeated below.

## `"conserve"`

* `mask_in` — mask for the source grid; excludes masked input cells from the interpolation.
  Values must be between 0 and 1.
* `mask_out` — returns the fractional area of each output grid cell covered by unmasked input cells.
* `is_latlon_in` / `is_latlon_out` — indicate whether the 2D input grid (`is_latlon_in`) and/or the
  destination grid (`is_latlon_out`) is a lat/lon grid.  If not provided, horiz_interp_new will automatically 
  check if the grids are rectilinear grids and trigger the "1D" algorithms for efficiency.
  
At the `horiz_interp` (interpolate) step, `mask_in` and `mask_out` are only used when `Interp%version=1`.

## `"bilinear"`

* `src_modulo` — `.true.` if the source grid is cyclic (periodic) in longitude.  Defaults to `.false.`
* `grid_at_center` — set to `.true.` if the input lon/lat coordinates (1D source grids only) are the
  grid cell center points. If `.false.` (default), the input coordinates are treated as cell
  boundaries and the centers are computed internally.

At the `horiz_interp` (interpolate) step:

* `missing_value` — a value in `data_in` to be treated as missing.
* `missing_permit` — the maximum number of surrounding missing values allowed before an output point
  is itself marked missing.
* `new_missing_handle` — set to `.true.` to turn on alternative missing-value handling.

## `"bicubic"`

* `src_modulo` — set to `.true.` if the source grid is cyclic (periodic) in longitude.
* `grid_at_center` — set to `.true.` if the input lon/lat coordinates are the grid cell center points.
  If `.false.` (default), the input coordinates are treated as cell boundaries and the centers are
  computed internally.

At the `horiz_interp` (interpolate) step:

* `missing_value` — a value in `data_in` to be treated as missing.
* `missing_permit` — the maximum number of surrounding missing values allowed before an output point
  is itself marked missing.

## `"spherical"`

* `num_nbrs` — number of nearest neighbors to use; defaults to 4.
* `max_dist` — maximum radius of influence [radians] around an output grid cell; defaults to 0.1 radians.
* `src_modulo` — set to `.true.` if the source grid is cyclic (periodic) in longitude.

At the `horiz_interp` (interpolate) step:

* `missing_value` — a value in `data_in` to be treated as missing.

## All methods, at the `horiz_interp` (interpolate) step

* `err_msg` — if present, interpolation errors (e.g. an uninitialized `Interp`) are returned in this
  string instead of triggering a fatal error.

## Notes

- `horiz_interp` also accepts 3D data arrays (`(i,j,k)`), interpolating each
  vertical level with the same horizontal weights.
- Both 32-bit (`r4`) and 64-bit (`r8`) real kinds are supported transparently
  through the generic interfaces.
- See `test_fms/horiz_interp/test_horiz_interp.F90` for complete working
  examples of every grid combination (1D/2D source × 1D/2D destination) and
  every interpolation method.
