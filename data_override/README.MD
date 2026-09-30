# Data_override

`data_override_mod` reads forcing data described by a legacy `data_table` or a
`data_table.yaml` file, interpolates it in space and time when requested, and
stores the result into model arrays on ATM, OCN, LND, or ICE grids. It also
handles the MPP domain decomposition. The public interfaces are
`data_override` for structured grids and `data_override_UG` for unstructured
land grids.

### Contents
- [Basic usage](README.MD#1-basic-usage)
- [YAML data table format](README.MD#2-yaml-data-table-format)
- [Building with YAML support](README.MD#3-building-with-yaml-support)
- [Converting a legacy data table](README.MD#4-converting-a-legacy-data-table)
- [Examples](README.MD#5-examples)
- [External weight file structure](README.MD#6-external-weight-file-structure)
- [Ensemble and nest support](README.MD#7-ensemble-and-nest-support)

#### 1. Basic usage

Define and initialize the component domain, call `data_override_init` after
the component model is initialized, then call `data_override` for each field
and model time. `data_override_init` requires `INPUT/grid_spec.nc`, typically
created by [fre-nctools](github.com/noaa-gfdl/frenc-tools), to describe the component grid. The dimensions in
that file must match the corresponding global MPP domain. Grid metadata may be
read from classic `grid_spec.nc` variables or from the mosaic descriptions
(`ocn_mosaic_file` and `gridfiles`).

The table file is normally `data_table` for the legacy format or
`data_table.yaml` for YAML. The selected table must be available in the run
directory; the YAML option also requires a build with YAML support and the
`use_data_table_yaml` namelist setting.

For example, the following YAML table defines the `sst_obs` field used below.
It reads `SST` from a NetCDF file, bilinearly interpolates it to the ocean
grid, and applies a factor of `1.0`:

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
logical :: override

call fms_init()

! Define and initialize Ocn_domain and allocate sst before this point.
Time = set_date(2000, 1, 1, 0, 0, 0)
call data_override_init(Ocean_domain_in=Ocn_domain)
call data_override('OCN', 'sst_obs', sst, Time, override=override)

if (.not. override) then
  ! The (OCN, sst_obs) entry was not found in the selected table.
endif

call fms_end()
```

The first argument is the three-character model grid ID (`'ATM'`, `'OCN'`,
`'LND'`, or `'ICE'`), and the second is the field name used by the model. The
third argument is the output value or array. For an array, the call replaces
the selected values with data from the table after spatial interpolation,
temporal interpolation to `Time`, and multiplication by `factor`. Outside a
configured subregion, existing array values are retained. The optional
`override` argument is set to `.true.` when the field is successfully read and
interpolated, and `.false.` when no matching table entry is found.

The generic `data_override` interface supports scalar, 2D, and 3D real values
in both `r4_kind` and `r8_kind`. A single scalar call illustrates the type
selection:

```fortran
real(r8_kind) :: co2
logical :: override

call data_override('OCN', 'co2_obs', co2, Time, override=override)
```

For threaded or multi-window calls, pass the window bounds relative to the
compute-domain indices. The array section and the four bounds must describe
the same window:

```fortran
call data_override('OCN', 'sst_obs', sst(isw:iew, jsw:jew), Time, &
                   override=override, is_in=is_rel, ie_in=ie_rel, &
                   js_in=js_rel, je_in=je_rel)
```

For an unstructured land grid, initialize with `Land_domainUG_in` and use
`data_override_UG`. It supports 1D and 2D real arrays:

```fortran
use mpp_domains_mod, only: domainUG
type(domainUG) :: Lnd_domain_ug
real(r8_kind), allocatable :: ug_field(:)

call data_override_init(Land_domainUG_in=Lnd_domain_ug)
call data_override_UG('LND', 'soil_moisture_obs', ug_field, Time, &
                      override=override)
```

If `mode` is omitted from `data_override_init`, both `r4_kind` and `r8_kind`
implementations are initialized. To initialize only one precision, pass
`mode=r4_kind` or `mode=r8_kind`.

```fortran
call data_override_init(Ocean_domain_in=Ocn_domain, mode=r8_kind)
```

When a domain is no longer needed, release it with
`data_override_unset_domains`. The optional `must_be_set` argument controls
whether attempting to unset an unset domain is fatal.

```fortran
use data_override_mod, only: data_override_unset_domains

call data_override_unset_domains(unset_Ocean=.true.)
call data_override_unset_domains(unset_Land=.true., must_be_set=.false.)
```

#### 2. YAML data table format

Each entry in `data_table.yaml` can have the following keys:

- **grid_name:** Target grid. The accepted values are `ATM`, `OCN`, `LND`,
  and `ICE`.
- **fieldname_in_model:** Field name as it appears in model code.
- **override_file:** Optional list describing data read from a NetCDF file.
  - **file_name:** File containing the variable. For `multi_file`, this is
    the middle file.
  - **fieldname_in_file:** Variable name in the NetCDF file.
  - **interp_method:** `bilinear`, `bicubic`, or `none`. Use `none` when the
    source data is already on the model grid. The LIMA format is not supported.
  - **multi_file:** Optional list with at most one entry. Use
    **prev_file_name** and **next_file_name** to provide the neighboring files
    needed for time interpolation across file boundaries.
  - **external_weights:** Optional list describing an external interpolation
    weight file.
    - **file_name:** External weight file.
    - **source:** Weight generator. The accepted value is `fregrid`.
- **factor:** Multiplier applied after the field is interpolated. It can also
  be used for unit conversion.
- **subregion:** Optional list limiting where the override is applied.
  - **type:** `inside_region` or `outside_region`.
  - **lon_start**, **lon_end:** Longitude bounds in the source grid units.
  - **lat_start**, **lat_end:** Latitude bounds in the source grid units.

#### 3. Building with YAML support

YAML tables require [libyaml](https://github.com/yaml/libyaml) to be installed
and linked with FMS, and FMS must be compiled with the `use_yaml` macro. With
Autotools:

```text
./configure --with-yaml
```

With CMake:

```text
cmake -DWITH_YAML=on
```

#### 4. Converting a legacy data table

The scripts
[`data_table_to_yaml.py`](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/main/data_table/data_table_to_yaml.py)
and
[`is_valid_data_table_yaml.py`](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/main/data_table/is_valid_data_table_yaml.py)
convert and validate YAML data tables.

#### 5. Examples

The following YAML entry bilinearly interpolates `sic` from a file to the ICE
grid, interpolates it to the model time, and multiplies it by `0.01`. The
result is written into the `icec` array by the corresponding call.

```yaml
data_table:
  - grid_name: ICE
    fieldname_in_model: sic_obs
    factor: 0.01
    override_file:
      - file_name: INPUT/hadisst_ice.data.nc
        fieldname_in_file: sic
        interp_method: bilinear
```

```fortran
call data_override('ICE', 'sic_obs', icec, Spec_Time)
```

A constant override omits `override_file`; `factor` is the value written to
the model field:

```yaml
data_table:
  - grid_name: ICE
    fieldname_in_model: sit_obs
    factor: 0.01
```

An on-grid field uses `interp_method: none`:

```yaml
data_table:
  - grid_name: OCN
    fieldname_in_model: runoff
    factor: 1.0
    override_file:
      - file_name: INPUT/runoff.daitren.clim.nc
        fieldname_in_file: runoff
        interp_method: none
```

The equivalent legacy entries are:

```text
"ICE", "sic_obs", "sic", "INPUT/hadisst_ice.data.nc", "bilinear", 0.01
"ICE", "sit_obs", "", "", "none", 0.01
"OCN", "runoff", "runoff", "INPUT/runoff.daitren.clim.nc", "none", 1.0
```

For a forcing period that crosses file boundaries, specify the previous and
next files. `file_name` is the middle file:

```yaml
data_table:
  - grid_name: ICE
    fieldname_in_model: sic_obs
    factor: 0.01
    override_file:
      - file_name: INPUT/hadisst_ice.data_yr1.nc
        fieldname_in_file: sic
        interp_method: bilinear
        multi_file:
          - prev_file_name: INPUT/hadisst_ice.data_yr0.nc
            next_file_name: INPUT/hadisst_ice.data_yr2.nc
```

External weights can be supplied for bilinear interpolation:

```yaml
data_table:
  - grid_name: ICE
    fieldname_in_model: sic_obs
    factor: 0.01
    override_file:
      - file_name: INPUT/hadisst_ice.data.nc
        fieldname_in_file: sic
        interp_method: bilinear
        external_weights:
          - file_name: INPUT/remap_file.nc
            source: fregrid
```

#### 6. External weight file structure

A bilinear weight file from `fregrid` has this structure:

```text
dimensions:
        nlon = 5 ;
        nlat = 6 ;
        three = 3 ;
        four = 4 ;
variables:
        int index(three, nlat, nlon) ;
        double weight(four, nlat, nlon) ;
```

- `nlon` and `nlat` must equal the global source-domain dimensions.
- `index(1,:,:)` is the source longitude index `i` nearest each model point.
- `index(2,:,:)` is the source latitude index `j` nearest each model point.
- `index(3,:,:)` is the source tile. It must be `1`; data_override does not
  support interpolation from cubed-sphere grids.
- The four corners are `(i,j)`, `(i,j+1)`, `(i+1,j)`, and `(i+1,j+1)`.
  Their weights are `weight(:,:,1)` through `weight(:,:,4)`, respectively.

#### 7. Ensemble and nest support

Each ensemble member can use its own YAML table. For two ensemble members, FMS
looks for `data_table_ens_01.yaml` and `data_table_ens_02.yaml`. Do not provide
both `data_table.yaml` and an ensemble-specific table: only one table choice is
allowed and the run will fail.

Nested domains can use `data_table_nest_01.yaml`, while the main grid uses
`data_table.yaml`. Unlike ensemble tables, a nest-specific table may coexist
with `data_table.yaml`.