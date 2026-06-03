The purpose of this document is to document the differences between the old (legacy) diag manager and the new (modern) diag manager.


## Contents
- [1. Diag Manager Rewrite Overview](README.md#1-overview-of-the-diag-manager-rewrite)
- [2. Diag Table Format](README.md#2-diag-table-format)
- [3. Scalar Axis](README.md#3-scalar-axis)
- [4. Average Time Variables](README.md#4-average-time-variables)
- [5. Subregional Files](README.md#5-subregional-files)
- [6. Global attributes](README.md#6-global-attributes)
- [7. Real attributes from diag_field_add_attribute calls](README.md#7-real-attributes-from-diag_field_add_attribute-calls)
- [8. History files data output "changes"](README.md#8-history-files-data-output-changes)

### 1. Overview of the Diag Manager Rewrite

The diag manager was completely rewritten to support YAML-formatted diagnostic tables and provide improved performance and maintainability. The rewrite maintains backward compatibility with the legacy ASCII-format diag tables through build-time configuration, but is only able to be used if FMS is built with libyaml support.

#### Module Organization and Naming Convention

The diag_manager is organized into functional modules, with new modules introduced in the rewrite using the `fms_diag_` prefix convention:

**Legacy Modules (original diag manager):**
- `diag_manager.F90` (top-level interface)
- `diag_axis.F90`
- `diag_data.F90`
- `diag_grid.F90`
- `diag_output.F90`
- `diag_table.F90`
- `diag_util.F90`
- `fms_diag_time_reduction.F90`
- `fms_diag_elem_weight_procs.F90`
- `fms_diag_fieldbuff_update.F90`

**New Modules (modern diag manager):**
- `diag_manager.F90` (top-level interface) Same routines as the legacy, but when use_modern_diag is true routines will make calls to fms_diag_object.F90 and friends to perform all operations.
- `fms_diag_object.F90` - Core diagnostic object implementation
- `fms_diag_field_object.F90` - Field-specific diagnostic data structures
- `fms_diag_file_object.F90` - File I/O management
- `fms_diag_axis_object.F90` - Axis and domain handling
- `fms_diag_bbox.F90` - Bounding box operations for subregional output
- `fms_diag_output_buffer.F90` - Output buffering and data management
- `fms_diag_reduction_methods.F90` - Reduction method implementations
- `fms_diag_input_buffer.F90` - Input buffer management
- `fms_diag_yaml.F90` - YAML diagnostic table parsing and handling
- `fms_diag_time_utils.F90` - Time utility functions
- `diag_data.F90` - stores all namelist parameters from diag_manager_nml

#### Enabling the Modern Diag Manager

The modern diag manager is optionally enabled via the `use_modern_diag` flag in the `diag_manager_nml` namelist, as seen below. FMS must be compiled with the `-Duse_yaml flag`.
By default, the legacy diag manager is used to maintain backward compatibility. When `use_modern_diag = .true.`, the modern implementation is invoked while maintaining the same public interface.

```
&diag_manager_nml
	use_modern_diag=.true.
/
```

### 2. Diag Table Format
The modern diag manager uses a YAML format instead of the legacy ascii table. A description of the YAML diag table can
be found [here](diag_yaml_format.md). A formal specification, in the form of a JSON schema, can be found in the
[gfdl_msd_schemas](https://github.com/NOAA-GFDL/gfdl_msd_schemas) repository on Github.

Options set in the file section will become the default for each variable in that file (ie. `kind: r8 below` will set all variables to use r8 kind),
unless otherwise specified. Ordering of the key-value pairs in the yaml file does not matter, as long as the appropriate sections have
matching indentation.

This barebones example creates a single netcdf file (per tile, if using a tiled domain):
```{yaml}
title: simple_diag_table
base_date: 1 1 1 0 0 0
diag_files:
- file_name: simple_diagnostics
  freq: 225 seconds
  time_units: seconds
  module: atm_mod
  kind: r8
  unlimdim: time
  varlist:
  - var_name: var1
    output_name: variable_one
	reduction: average
```

This is a more complex example utilizing subregional output and wildcard filenames:
```{yaml}
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: normal
  freq: 24 days
  time_units: hours
  unlimdim: records
  module: potato_mod
  kind: r8
  reduction: min
  varlist:
  - module: atm_mod 
    var_name: sst
    output_name: sst
    reduction: average
    kind: r4
    write_var: true
    attributes:
    - do_sst: .true.
  sub_region:
  - grid_type: latlon
    corner1: -80, 0
    corner2: -80, 75
    corner3: -60, 0
    corner4: -60, 75
- file_name: normal2
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: true
  module: atm_mod 
  reduction: none
  kind: r4
  varlist:
  - var_name: sstt
    output_name: sstt
    long_name: S S T
  - var_name: sstt2
    output_name: sstt2
    long_name: S S T
    write_var: false
  sub_region:
  - grid_type: index
    tile: 1
    corner1: 10, 15
    corner2: 20, 15
    corner3: 10, 25
    corner4: 20, 25
- file_name: normal3
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: false
- file_name: wild_card_name%4yr%2mo%2dy%2hr
  filename_time: end
  freq: 6 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours
  module: ocn_mod 
  reduction: average
  kind: r4
  varlist:
  - var_name: sst
    output_name: sst
  global_meta:
  - is_a_file: true
```

### 3. Scalar Axis
The old diag manager was adding a `scalar_axis` dimension of size 1 for scalar variables

```
dimensions:
    	scalar_axis = 1 ;
variables:
    	double p700(scalar_axis) ;
            	p700:_FillValue = 1.e+20 ;
            	p700:missing_value = 1.e+20 ;
```
The new diag manager will no longer have a dummy scalar axis dimension.

### 4. Average Time Variables
The old diag manager includes time bounds metadata in a non-standard convention (i.e. `average_T1`, `average_T2`, and `average_DT`)
1. `average_T1` is the start time for the averaging period (in the same time units as time)
2. `average_T2` is the end time for the averaging period
3. `average_DT` is the length of the averaging period, in days

```
	double average_T1(time) ;
		average_T1:_FillValue = 1.e+20 ;
		average_T1:missing_value = 1.e+20 ;
		average_T1:units = "days since 1979-01-01 00:00:00" ;
		average_T1:long_name = "Start time for average period" ;
	double average_T2(time) ;
		average_T2:_FillValue = 1.e+20 ;
		average_T2:missing_value = 1.e+20 ;
		average_T2:units = "days since 1979-01-01 00:00:00" ;
		average_T2:long_name = "End time for average period" ;
	double average_DT(time) ;
		average_DT:_FillValue = 1.e+20 ;
		average_DT:missing_value = 1.e+20 ;
		average_DT:units = "days" ;
		average_DT:long_name = "Length of average period" ;
```
These 3 variables are referenced as a variable attribute in each diagnostic. e.g.
```
		dis_liq:time_avg_info = "average_T1,average_T2,average_DT" ;
```

The new diag manager will not be using these non standard variables. Instead the time bounds information will be specified in CF standards. One time_bounds variable, with an extra nv dimension, that species both average_T1 and average_T2. The average_DT information can be obtained from the time_bnds variable.

```
	double nv(nv) ;
		nv:long_name = "vertex number" ;
	double time_bnds(time, nv) ;
		time_bnds:units = "days since 1979-01-01 00:00:00" ;
		time_bnds:long_name = "time axis boundaries" ;
```
This time_bounds variable is refernced as a variable attribute of time:
```
		time:bounds = "time_bnds" ;
```

### 5. Subregional Files

#### A. `is_subregional` global attribute:
Subregional files will have a global NetCDF attribute `is_subregional = True` set for non-global history files. This attribute will be used in PP tools.

#### B. Subregional dimension names:
In some cases, the old diag manager was adding `sub0X` to the dimension names where X is a number greater than 1. This was causing problems in PP tools that were expecting the dimension to have `sub01` in the name. The new diag manager will not have this problem.

#### C. Corner and center diagnostics:
In the old diag manager, if mixing variables that are corner variables, such as velocities={uo,vo,umo,vmo} and center variables, such as tracers={thetao,so,volcello} you sometimes ended up with a different number of variables per file. The extra files had duplicate data for the corner velocities because the two PEs shared the point at the edge. This happened with some grid/layouts/masks/subregion combinations and it caused problems with the combiner. The new diag manager will not have this problem.

### 6. Global attributes
#### A. Grid type and grid tile:
The old diag manager was adding the global attributes grid_type = "regular" and grid_tile = "N/A" for all files regardless of what the grid_type and the grid_title actually were. The new diag manager will no longer be doing this as they are not correct and don’t seem to be used.

#### B. Associated_files global attribute:
We were unable to reproduce the exact order of the associated_files global attribute, so users may see differences like

```
lake_area: 19790101.land_static.nc soil_area: 19790101.land_static.nc land_area: 19790101.land_static.nc <> land_area: 19790101.land_static.nc soil_area: 19790101.land_static.nc lake_area: 19790101.land_static.nc
```

### 7. Real attributes from diag_field_add_attribute calls
When real attributes were added to the file via a diag_field_add_attribute call, the old diag manager is always saving it as NF90_FLOAT regardless of the precision the data was [passed in](https://github.com/NOAA-GFDL/FMS/blob/ebb32649efa395ea14598f74c8d49e74d1408579/diag_manager/diag_manager.F90#L4532-L4543)

The new diag manager is going to write the attribute as it is passed in. This will cause differences when the model component was compiled with r8 as it will write the attribute as r8 instead of r4.

### 8. History files data output "changes"
When the model run time is less than then the output frequency (i.e if the module run time is 2 days and you are writing monthly diagnostics), the old diag manager was writing 9.96921e+36. The new diag manager is not going to write anything for this cases, so if you ncdump the output from the new diag manager, you will get:

```
 wa =
  _, _, _, _, _, _, ...
```

Similarly, when a variable was registered, but send_data was never called, the old diag manager was outputting the warning like

```
WARNING from PE     0: diag_manager_mod::closing_file: module/output_field soil/soil_fgw, skip one time level, maybe send_data never called
```

And writing out `9.96921e+36` for the variable. The new diag manager will also be outputting the warning, but it will not write out anything.


