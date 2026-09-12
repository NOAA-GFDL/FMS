# QUICKSTART.md

# Introduction

`data_override_mod` reads forcing data described in `data_table` or
`data_table.yaml`, optionally interpolates it in space/time, and writes the
result into model arrays on ATM/OCN/LND/ICE grids.

Typical uses:
- Override a field with data from a file
- Override a field with a constant value (`fieldname_in_file: ""`)
- Run on-grid (`interp_method: none`) or off-grid (`bilinear`/`bicubic`)
- Restrict override to subregions
- Use multi-file time bridging (`prev_file_name`/`next_file_name`)
- Use external interpolation weights (for example from `fregrid`)

For full format details, see `README.MD` in this directory.

# Basic Usage

## 3-step workflow

1. Define and initialize model domain(s).
2. Call `data_override_init` once domains are available.
3. Call `data_override` (or `data_override_UG`) for each field and time.

```fortran
use fms_mod,          only: fms_init, fms_end
use mpp_domains_mod,  only: domain2d
use time_manager_mod, only: time_type, set_date
use data_override_mod, only: data_override_init, data_override
use platform_mod,     only: r8_kind

implicit none

type(domain2d) :: Ocn_domain
type(time_type) :: Time
real(r8_kind), allocatable :: sst(:, :)
logical :: used

call fms_init()

! ... define and initialize Ocn_domain and allocate sst ...
Time = set_date(2000, 1, 1, 0, 0, 0)

! Initialize data_override for the ocean domain.
call data_override_init(Ocean_domain_in=Ocn_domain)

! Apply override (example field names must exist in your data_table).
call data_override('OCN', 'sst_obs', sst, Time, override=used)

if (.not. used) then
  ! Field not found in data_table/data_table.yaml.
endif

call fms_end()
```

## Optional precision mode

If omitted, both r4 and r8 implementations are initialized.

```fortran
use platform_mod, only: r4_kind, r8_kind

call data_override_init(Ocean_domain_in=Ocn_domain, mode=r4_kind)
call data_override_init(Ocean_domain_in=Ocn_domain, mode=r8_kind)
```

# Interface Call Patterns

The generic interfaces dispatch by rank and real kind.

## Structured-grid interface: `data_override`

### 0D scalar

```fortran
real(r4_kind) :: co2_r4
real(r8_kind) :: co2_r8
logical :: used

! data_override will populate co2_r4 and co2_r8 with values from the data file
call data_override('OCN', 'co2_obs', co2_r4, Time, override=used)
call data_override('OCN', 'co2_obs', co2_r8, Time, override=used)
```

### 2D field

```fortran
real(r4_kind), allocatable :: fld2d_r4(:, :)
real(r8_kind), allocatable :: fld2d_r8(:, :)
logical :: used

call data_override('OCN', 'sst_obs', fld2d_r4, Time, override=used)
call data_override('OCN', 'sst_obs', fld2d_r8, Time, override=used)
```

### 3D field

```fortran
real(r4_kind), allocatable :: fld3d_r4(:, :, :)
real(r8_kind), allocatable :: fld3d_r8(:, :, :)
logical :: used

call data_override('LND', 'soil_temp_obs', fld3d_r4, Time, override=used)
call data_override('LND', 'soil_temp_obs', fld3d_r8, Time, override=used)
```

### Windowed compute-domain calls

When threading or splitting windows (for OpenMP parallelization), pass bounds relative to compute-domain
indices. This allows multiple threads to safely call `data_override` on different regions of the same array.

```fortran
call data_override('OCN', 'sst_obs', fld2d_r4(isw:iew, jsw:jew), Time, &
                   override=used, is_in=is_rel, ie_in=ie_rel, &
                   js_in=js_rel, je_in=je_rel)
```

## Unstructured-grid interface: `data_override_UG`

Initialize with `Land_domainUG_in` first.

```fortran
use mpp_domains_mod, only: domainUG

type(domainUG) :: Lnd_domain_ug
logical :: used

call data_override_init(Land_domainUG_in=Lnd_domain_ug)
```

### 1D UG field

```fortran
real(r4_kind), allocatable :: ug1d_r4(:)
real(r8_kind), allocatable :: ug1d_r8(:)

call data_override_UG('LND', 'sst_obs', ug1d_r4, Time, override=used)
call data_override_UG('LND', 'sst_obs', ug1d_r8, Time, override=used)
```

### 2D UG field

```fortran
real(r4_kind), allocatable :: ug2d_r4(:, :)
real(r8_kind), allocatable :: ug2d_r8(:, :)

call data_override_UG('LND', 'sst_obs', ug2d_r4, Time, override=used)
call data_override_UG('LND', 'sst_obs', ug2d_r8, Time, override=used)
```

# Minimal Data Table Examples

## YAML: file-based bilinear override

```yaml
data_table:
  - grid_name: OCN
    fieldname_in_model: sst_obs
    factor: 1.0
    override_file:
      - file_name: INPUT/sst_ice_clim.nc
        fieldname_in_file: SST
        interp_method: bilinear
```

## YAML: constant override

```yaml
data_table:
  - grid_name: OCN
    fieldname_in_model: co2_obs
    factor: 300.0
```

## Legacy ASCII

```text
"OCN", "sst_obs", "SST", "INPUT/sst_ice_clim.nc", "bilinear", 1.0
"OCN", "co2_obs", "",    "",                      "none",     300.0
```

# Initialization Notes

- `data_override_init` expects `INPUT/grid_spec.nc` to be available.
- Grid metadata can be read from classic `grid_spec.nc` variables or mosaic
  descriptions (`ocn_mosaic_file`/`gridfiles`).
- `grid_name` keys in YAML are uppercase (`ATM`, `OCN`, `LND`, `ICE`).
- Calls to `data_override` use 3-character grid ids like `'OCN'`.

# Optional Features

## Subregion override

Use YAML `subregion` with `inside_region` or `outside_region` to limit where
override applies.

## Multi-file time bridging

Use `multi_file` with `prev_file_name` and/or `next_file_name` when forcing
windows cross file boundaries.

## External weights

Use `external_weights` when you want interpolation weights supplied by file.
Currently, source support is `fregrid`.

# Cleanup

Unset any domains that should no longer be used:

```fortran
use data_override_mod, only: data_override_unset_domains

call data_override_unset_domains(unset_Ocean=.true.)
call data_override_unset_domains(unset_Land=.true., must_be_set=.false.)
```

# Troubleshooting

- `override=.false.` usually means the `(grid, fieldname_in_model)` pair was
  not found in the data table.
- A fatal init error often indicates grid/domain size mismatch against
  `grid_spec.nc` or mosaic tiles.
- For YAML usage, build with YAML support and use the `use_data_table_yaml`
  namelist setting as expected by your configuration.

# Additional Notes

- `data_override_mod` supports both `r4_kind` and `r8_kind`.
- Interpolated coordinates are internally handled in radians.
- This quickstart focuses on calling patterns; see `README.MD` for the full
  schema and detailed examples.
