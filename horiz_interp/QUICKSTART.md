# Introduction

`horiz_interp_mod` interpolates 2D and 3D data with `conservative order 1`,
`bilinear`, `bicubic`, and `spherical` (inverse-distance weighted) methods.
Not all methods support all input and output grids.  The units of source and destination
grids must be in radians.

# Basic usage

`Horiz_interp_new` and `horiz_interp` are Fortran generic interfaces that accept
varying lists of input arguments.  For example, `horiz_interp_new` will accept different combinations
of grids represented as 1D and 2D arrays (Rectilinear grids such as lat-lon grid can be defined
with 1D arrays).  `Horiz_interp` will accept arguments with and without `Interp` (an instance of `horiz_interp_type`):
when `horiz_interp` is called with `Interp` as the first argument followed by the input and output data, `horiz_interp`
will use the interpolation weights stored in `Interp` to interpolate the input data to the output data.  If `Interp`
is not the first argument, but instead the input and output grids are provided as the first arguments followed by
the input and output data, `horiz_interp` will take the "blackbox" approach:  weights will be computed
on-the-fly followed by data interpolation.  For the "blackbox" approach, weights are not stored and
are re-computed for every call.

## 2-step interpolation

1. Call `horiz_interp_init`.
2. Pass in grid to `horiz_interp_new` to compute and store weights in `Interp`.
3. Pass `Interp` and `data` to `horiz_interp` to interpolate.
4. Free `Interp` with `horiz_interp_del` when done.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp_new, horiz_interp_end, &
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

! Populate Interp with interpolation weights and mapping indices
! (Interp can be reused for data fields on the same grid)
call horiz_interp_new(Interp, lon_src, lat_src, lon_dst, lat_dst, interp_method="bilinear")

! Interpolate data into data_dst
call horiz_interp(Interp, data_src, data_dst)

! deallocate memory stored in Interp
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
use horiz_interp_mod, only: horiz_interp_init, horiz_interp, horiz_interp_del, horiz_interp_end
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

Rectilinear grids such as a lat-lon grids (the longitude coordinates are
identical along every line of latitude and the latitude coordinates are identical
along every line of longitude) can be represented as 1D array of longitude gridpoints and
1D array of latitude gridpoints, leading to simpler and more efficient weight-generating algorithms.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp_new, horiz_interp_end, &
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
integer :: i

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
! (Interp can be reused for data fields on the same grid)
call horiz_interp_new(Interp, lon_src, lat_src, lon_dst, lat_dst, interp_method="bilinear")

! Interpolate data into data_dst
call horiz_interp(Interp, data_src, data_dst)

! Deallocate memory in Interp
call horiz_interp_del(Interp)

! finalize horiz_interp_mod
call horiz_interp_end()

! end program
call fms_end()
```

## 3D data
If the data provided is a 3D array, e.g data(i,j,k), horiz_interp_mod will conduct spatial
interpolation for each horizontal slice along the k-axis.  Note, horiz_interp_mod does not
support vertical interpolation.

```fortran
use horiz_interp_mod, only: horiz_interp_init, horiz_interp_new, horiz_interp_end, &
                             horiz_interp, horiz_interp_del, horiz_interp_type
use platform_mod, only:  r8_kind ! 64-bit floating point precision

implicit none

type(horiz_interp_type) :: Interp
real(r8_kind), allocatable :: lon_src(:, :), lat_src(:, :) ! source grid
real(r8_kind), allocatable :: lon_dst(:,:), lat_dst(:,:) ! destination grid
real(r8_kind), allocatable :: data_src(:,:,:), data_dst(:,:,:) ! data

! ... allocate and define the source grid in lon_src and lat_src
! allocate and define the destination grid in lon_dst and lat_dst
! ... allocate data_src and compute data on source grid

! start program
call fms_init()

! ... set domains and MPI parallelization

! initialize horiz_interp_mod
call horiz_interp_init()

! Compute interpolation weights and mapping indices
! (Interp can be reused for data fields on the same grid)
call horiz_interp_new(Interp, lon_src, lat_src, lon_dst, lat_dst, interp_method="bilinear")

! Compute interpolation weights and mapping indices
! (Interp can be reused for data fields on the same grid)
call horiz_interp(Interp, data_src, data_dst)

! Deallocate memory in Interp
call horiz_interp_del(Interp)

! finalize horiz_interp_mod
call horiz_interp_end()

! end program
call fms_end()
```

# Interpolation methods
It is recommended to set the optional argument `interp_method` to one of the following:

* `"conservative"` — order 1 conservative interpolation.  Supports all combination of 1D and 2D source and destination grids.
* `"bilinear"` — bilinear interpolation.  Only supports 1D source with 2D destination, or 2D source and destination grids.
* `"bicubic"` — bicubic interpolation.  Only supports 1D source with 2D destination, or 1D source and destination grids.
* `"spherical"` — inverse-distance weighting over nearest neighbors.  Only supports 2D source and destination grids.


# Horiz_interp_new optional arguments

`Horiz_interp_new` supports the following optional arguments for each interpolation methods.
All methods support `verbose` (sets the verbosity level, 0/1/2, default 0) and `interp_method` (selects the
method itself) are accepted by `horiz_interp_new` regardless of method and are not repeated below.

## conservative

* `mask_in`:  mask for the source grid; excludes masked input cells from the interpolation.
  Values must be between 0 and 1.  If not provided, all source grid cells will be remapped to the destination grid.
* `mask_out`:   will be populated with destination grid cells area computed from the exchange grid.
* `is_latlon_in`/`is_latlon_out`:  .true. if the 2D input/destination grid is a lat/lon grid.  If not provided,
  `horiz_interp_new` will automatically check the grids and trigger the more efficient "1D" algorithms.


## bilinear

* `src_modulo`:  `.true.` if the source grid is cyclic (periodic) in longitude.  Defaults to `.false.` if not provided.
* `grid_at_center`:  set to `.true.` if the input lon/lat coordinates (1D source grids only) are the
  grid cell center points. If not provided, defaults to `.false.` and the input coordinates are treated as edge points
  to compute the centers.


## bicubic

* `src_modulo` — set to `.true.` if the source grid is cyclic (periodic) in longitude.
* `grid_at_center`:  set to `.true.` if the input lon/lat coordinates (1D source grids only) are the
  grid cell center points. If not provided, defaults to `.false.` and the input coordinates are treated as edge points
  to compute the centers internally.


## spherical

* `num_nbrs` — number of nearest neighbors to find.  Defaults to 4 if not provided.
* `max_dist` — maximum distance in radians to search for neighboring cells on the input grid for each output grid cell.  Defaults to 0.1 radians if not provided.
* `src_modulo` — set to `.true.` if the source grid is cyclic (periodic) in longitude.  Defaults to `.false.` if not provided.


## Horiz_interp optional arguments (when Interp is provided)

## conservative

* `mask_in`/ `mask_out`:  Used only when `Interp%version = 1` (when weights were generated from 1D representation of both input and output grids.)

## bilinear

* `missing_value`:  data on the source grid with values equal to the `missing_value` will be treated as missing.
* `missing_permit`:  maximum number of missing values permitted.  Defaults to 0 if not provided.  Is not used if `new_missing_handle` is `.true.`
* `new_missing_handle`:  set to `.true.` to turn on the new missing handle algorithm.  Defaults to `.false.`

## bicubic

* `missing_value`:  data on the source grid with values equal to the `missing_value` will be treated as missing.
* `missing_permit`:  maximum number of missing values permitted.  Defaults to 0 if not provided.  Is not used if `new_missing_handle` is `.true.`

## spherical

* `missing_value`:  data on the source grid with values equal to the `missing_value` will be treated as missing.


# Additional notes
- Horiz_interp_mod supports both 32-bit (`r4_kind`) and 64-bit (`r8_kind`) grids and data.
  If `Interp` was populated with weights computed from grids in 32-bit precision, data is expected to be in 32-bit precision.
  Likewise, if `Interp` was populated with weights computed from grids in 64-bit precision, data is expected to be in 64-bit precision.
