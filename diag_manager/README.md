The purpose of this document is to document the differences between the old diag manager and the new (modern) diag manager. 

## Contents
- [1. Diag Table Format](README.md#1-diag-table-format)
- [2. Scalar Axis](README.md#2-scalar-axis)
- [3. Average Time Variables](README.md#3-average-time-variables)
- [4. Subregional Files](README.md#4-subregional-files)
- [5. Global attributes](README.md#5-global-attributes)
- [6. Real attributes from diag_field_add_attribute calls](README.md#6-real-attributes-from-diag_field_add_attribute-calls)
- [7. History files data output "changes"](README.md#7-history-files-data-output-changes)

### 1. Diag Table Format
The modern diag manager uses a YAML format instead of the legacy ascii table. A description of the YAML diag table can be found [here](diag_yaml_format.md).

### 2. Scalar Axis
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

### 3. Average Time Variables
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

### 4. Subregional Files

#### A. `is_subregional` global attribute:
Subregional files will have a global NetCDF attribute `is_subregional = True` set for non-global history files. This attribute will be used in PP tools. 

#### B. Subregional dimension names:
In some cases, the old diag manager was adding `sub0X` to the dimension names where X is a number greater than 1. This was causing problems in PP tools that were expecting the dimension to have `sub01` in the name. The new diag manager will not have this problem.

#### C. Corner and center diagnostics:
In the old diag manager, if mixing variables that are corner variables, such as velocities={uo,vo,umo,vmo} and center variables, such as tracers={thetao,so,volcello} you sometimes ended up with a different number of variables per file. The extra files had duplicate data for the corner velocities because the two PEs shared the point at the edge. This happened with some grid/layouts/masks/subregion combinations and it caused problems with the combiner. The new diag manager will not have this problem.

### 5. Global attributes
#### A. Grid type and grid tile:
The old diag manager was adding the global attributes grid_type = "regular" and grid_tile = "N/A" for all files regardless of what the grid_type and the grid_title actually were. The new diag manager will no longer be doing this as they are not correct and don’t seem to be used.

#### B. Associated_files global attribute:
We were unable to reproduce the exact order of the associated_files global attribute, so users may see differences like

```
lake_area: 19790101.land_static.nc soil_area: 19790101.land_static.nc land_area: 19790101.land_static.nc <> land_area: 19790101.land_static.nc soil_area: 19790101.land_static.nc lake_area: 19790101.land_static.nc
```

### 6. Real attributes from diag_field_add_attribute calls
When real attributes were added to the file via a diag_field_add_attribute call, the old diag manager is always saving it as NF90_FLOAT regardless of the precision the data was [passed in](https://github.com/NOAA-GFDL/FMS/blob/ebb32649efa395ea14598f74c8d49e74d1408579/diag_manager/diag_manager.F90#L4532-L4543)

The new diag manager is going to write the attribute as it is passed in. This will cause differences when the model component was compiled with r8 as it will write the attribute as r8 instead of r4.

### 7. History files data output "changes"
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


