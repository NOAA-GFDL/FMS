## Diag Table Yaml Format:

The purpose of this documents is to explain the diag_table yaml format. 

## Contents
- [1. Coverting from legacy ascii diag_table format](README.md#1-coverting-from-legacy-ascii-diag_table-format)
- [2. Diag table yaml sections](README.md#2-diag-table-yaml-sections)
- [2.1 Global Section](README.md#21-global-section)
- [2.2 File Section](README.md#22-file-section)
- [2.2.1 Flexible output timings](README.md#221-flexible-output-timings)
- [2.2.2 Coupled Model Diag Files](README.md#222-coupled-model-diag-files)
- [2.3 Variable Section](README.md#23-variable-section)
- [2.4 Variable Metadata Section](README.md#24-variable-metadata-section)
- [2.5 Global Meta Data Section](README.md#25-global-meta-data-section)
- [2.6 Sub_region Section](README.md#26-sub_region-section)
- [3. More examples](README.md#3-more-examples)

### 1. Coverting from legacy ascii diag_table format

To convert the legacy ascii diad_table format to this yaml format, the python script [**diag_table_to_yaml.py**](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/aafc3293d45df2fc173d3c7afd8b8b0adc18fde4/fms_yaml_tools/diag_table/diag_table_to_yaml.py#L23-L26) can be used. To confirm that your diag_table.yaml was created correctly, the python script [**is_valid_diag_table_yaml.py**](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/aafc3293d45df2fc173d3c7afd8b8b0adc18fde4/fms_yaml_tools/diag_table/is_valid_diag_table_yaml.py#L24-L27) can be used.

### 2. Diag table yaml sections
The diag_table.yaml is organized by file.  Each file has the required and optional key/value pairs for the file, an optional subsection defining any additional global metadata to add to the file, an optional subsection defining a subregion of the grid to output the data for and a required subsection for all of the variables in the file. Each variable has the required and optional key/value pairs for the variable and an optional subsection defining any additional variable attributes to add to the file. The hierarchical structure looks like this:

```yaml
title:
base_date:
diag_files:
- file1
  - #key/value pairs for file1
  varlist:
  - var1
    - #key/value pairs for var1
    attributes:
    - #atributes for var1
  global_metadata:
  - #global attributes for file1
  subregion:
  - #subregion for file1
```

### 2.1 Global Section
The diag_yaml requires “title” and the “baseDate”.
- The **title** is a string that labels the diag yaml.  The equivalent in the diag table would be the experiment.  It is recommended that each diag_yaml have a separate title label that is descriptive of the experiment that is using it.
- The **basedate** is an array of 6 integer indicating the base_date in the format [year month day hour minute second].

**Example:** 

In the YAML format:
```yaml
title: ESM4_piControl
base_date: 2022 5 26 12 3 1
```

In the legacy ascii format:
```
ESM4_piControl
2022 5 26 12 3 1
```

### 2.2 File Section
The files are listed under the diagFiles section as a dashed array. 

Below are the **required** keys needed to define each file.
- **file_name** is a string that defines the name of the file. Do not add ".nc" and "tileX" to the filename as this will handle by FMS. 
- **freq** is an integer that defines the frequency that data will be written. The acceptable values are:
  - =-1: output at the end of the run only 
  - =0: output every timestep 
  - \>0: output frequency
- **freq_units** is a string that defines the units of the frequency from above. The acceptable values are seconds, minutes, hours, days, months, years. 
- **time_units** is a string that defines units for time. The acceptable values are seconds, minutes, hours, days, months, years. 
- **unlimdim** is a string that defines the name of the unlimited dimension in the output netcdf file, usually “time”.
- **varlist** is a subsection that list all of the variable in the file

**Example:** The following creates a file with data written every 6 hours. 

In the YAML format:
```yaml
diag_files:
- file_name: atmos_6hours
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  varlist:
  - varinfo
```

In the legacy ascii format:
```
"atmos_6hours",      6,  "hours", 1, "hours", "time"
```

**NOTE:** The fourth column (file_format) has been deprecated. Netcdf files will always be written.

Below are some *optional* keys that may be added. 
- **write_file** is a logical that indicates if you want the file to be created (default is true). This is a new feature that is not supported by the legacy ascii data_table.
- **new_file_freq** is a integer that defines the frequency for closing the existing file
- **new_file_freq_units** is a string that defines the time units for creating a new file. Required if “new_file_freq” used. The acceptable values are seconds, minuts, hours, days, months, years. 
- **start_time** is an array of 6 integer indicating when to start the file for the first time. It is in the format [year month day hour minute second]. Requires “new_file_freq”
- **filename_time** is the time used to set the name of new files when using new_file_freq. The acceptable values are begin (which will use the begining of the file's time bounds), middle (which will use the middle of the file's time bounds), and end (which will use the end of the file's time bounds). The default is middle

**Example:** The following will create a new file every 6 hours starting at Jan 1 2020. Variable data will be written to the file every 6 hours.

In the YAML format:
```yaml
- file_name: ocn%4yr%2mo%2dy%2hr
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6
  new_file_freq_units: hours
  start_time: 2020 1 1 0 0 0
```

In the legacy ascii format:
```
"ocn%4yr%2mo%2dy%2hr",      6,  "hours", 1, "hours", "time", 6, "hours", "1901 1 1 0 0 0"
```

Because this is using the default `filename_time` (middle), this example will create the files:
```
ocn_2020_01_01_03.nc for time_bnds [0,6]
ocn_2020_01_01_09.nc for time_bnds [6,12]
ocn_2020_01_01_15.nc for time_bnds [12,18]
ocn_2020_01_01_21.nc for time_bnds [18,24]
```

**NOTE** If using the new_file_freq, there must be a way to distinguish each file, as it was done in the example above. 

- **file_duration** is an integer that defines how long the file should receive data after start time in “file_duration_units”.  This optional field can only  be used if the start_time field is present.  If this field is absent, then the file duration will be equal to the frequency for creating new files. The file_duration_units field must also be present if this field is present.
- **file_duration_units** is a string that defines the file duration units. The acceptable values are seconds, minutes, hours, days, months, years. 
- **global_meta** is a subsection that lists any additional global metadata to add to the file. This is a new feature that is not supported by the legacy ascii data_table.
- **sub_region** is a subsection that defines the four corners of a subregional section to capture.

### 2.2.1 Flexible output timings

In order to provide more flexibility in output timings, the new diag_table yaml format allows for different file frequencies for the same file by allowing the `freq`, `freq_units`, `new_file_freq`, `new_file_freq_units`, `file_duration`, `file_duration_units` keys to accept array of integers/strings. 

For example, 
``` yaml
- file_name: flexible_timing%4yr%2mo%2dy%2hr
  freq: 1 1 1
  freq_units: hours hours hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 3 1
  new_file_freq_units: hours hours hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 3 9
  file_duration_units: hours hours hours
  filename_time: begin
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
```
This will create a file every 6 hours for 12 hours
```
flexible_timing_0002_01_01_00.nc - using hourly averaged data from hour 0 to hour 6
flexible_timing_0002_01_01_06.nc - using hourly averaged data from hour 6 to hour 12
```

Then it will create a file every 3 hours for 3 hours
```
flexible_timing_0002_01_01_12.nc - using hourly averaged data from hour 12 to hour 15
```

Then it will create a file every 1 hour for 9 hours.
```
flexible_timing_0002_01_01_15.nc - using data from hour 15 to hour 16
flexible_timing_0002_01_01_16.nc - using data from hour 16 to hour 17
flexible_timing_0002_01_01_17.nc - using data from hour 17 to hour 18
flexible_timing_0002_01_01_18.nc - using data from hour 18 to hour 19
flexible_timing_0002_01_01_19.nc - using data from hour 19 to hour 20
flexible_timing_0002_01_01_20.nc - using data from hour 20 to hour 21
flexible_timing_0002_01_01_21.nc - using data from hour 21 to hour 22
flexible_timing_0002_01_01_22.nc - using data from hour 22 to hour 23
flexible_timing_0002_01_01_23.nc - using data from hour 23 to hour 24

```

### 2.2.2 Coupled Model Diag Files
In the *legacy ascii diag_table*, when running a coupled model (ATM + OCN) in a seperate PE list:
  - The ATM PEs ignored the files in the diag_table that contain "OCEAN" in the filename
  - The OCN PEs ignored the files in the diag_table that did not contain "OCEAN" in the filename

In the *yaml diag_table*:
  - The ATM PEs will ignore the files in the diag_table.yaml that contain the key/value pair `is_ocean: true`
  - The OCN PEs will ignore the files in the diag_table.yaml that do not contain the key/value pair `is_ocean: true`

### 2.3 Variable Section
The variables in each file are listed under the varlist section as a dashed array.

- **var_name:**  is a string that defines the variable name as it is defined in the register_diag_field call in the model
- **reduction:** is a string that describes the data reduction method to perform prior to writing data to disk. Acceptable values are average, diurnalXX (where XX is the number of diurnal samples), powXX (whre XX is the power level), min, max, none, rms, and sum.  
- **module:**  is a string that defines the module where the variable is registered in the model code
- **kind:** is a string that defines the type of variable  as it will be written out in the file. Acceptable values are r4, r8, i4, and i8

**Example:**

In the YAML format:
```yaml
  varlist:
  - module: moist
    var_name: precip
    reduction: average
    kind: r4
```

In the legacy ascii format:
```
"moist",     "precip",                         "precip",           "atmos_8xdaily",   "all", .true.,  "none", 2
```
**NOTE:** The fifth column (time_sampling) has be deprecated. The reduction_method (`.true.`) has been replaced with `average`. The output name was not included in the yaml because it is the same as the var_name. 

which corresponds to the following model code
```F90
id_precip = register_diag_field ( 'moist', 'precip', axes, Time)
```
where:
- `moist` corresonds to the module key in the diag_table.yaml
- `precip` corresponds to the var_name key in the diag_table.yaml
- `axes` are the ids of the axes the variable is a function of
- `Time` is the model time

Below are some *optional* keys that may be added. 
- **write_var:** is a logical that is set to false if the user doesn’t want the variable to be written to the file (default: true).
- **out_name:** is a string that defines the name of the variable that will be written to the file (default same as var_name)
- **long_name:** is a string defining the long_name attribute of the variable. It overwrites the long_name in the variable's register_diag_field call
- **attributes:** is a subsection with any additional metadata to add to the variable in the netcdf file. This is a new feature that is not supported by the legacy ascii data_table.
- **zbounds:** is a 2 member array of integers that define the bounds of the z axis (zmin, zmin), optional default is no limits. 

### 2.4 Variable Metadata Section
Any aditional variable attributes can be added for each varible can be listed under the attributes section as a dashed array. The key is attribute name and the value is the attribute value.

**Example:**

```yaml
    attributes:
    - attribute_name: attribute_value
      attribute_name: attribute_value
```

Although this was not supported by the legacy ascii data_table, with the legacy diag_manager, a call to `diag_field_add_attribute` could have been used to do the same thing.

```F90
call diag_field_add_attribute(diag_field_id, attribute_name, attribute_value)
```

### 2.5 Global Meta Data Section
Any aditional global attributes can be added for each file can be listed under the global_meta section as a dashed array.  The key is the attribute name and the value is the attribute value.

```yaml
  global_meta:
  - attribute_name: attribute_value
    attribute_name: attribute_value
```

### 2.6 Sub_region Section
The sub region can be listed under the sub_region section as a dashed array. The legacy ascii diag_table only allows regions to be defined using the latitude and longitude, and it only allowed rectangular sub regions. With the yaml diag_table, you can use indices to defined the sub_region and you can define **any** four corner shape. Each file can only have 1 sub_region defined. These are keys that can be used:
- **grid_type:** is a **required** string defining the method used to define the  fourth sub_region corners. The acceptable values are "latlon" if using latitude/longitude or "indices" if using the indices of the corners.
- **corner1:** is a **required** 2 member array of reals if using (grid_type="latlon") or integers if using (grid_type="indices") defining the x and y points of the first corner of a sub_grid.
- **corner2:** is a **required** 2 member array of reals if using (grid_type="latlon") or integers if using (grid_type="indices") defining the x and y points of the second corner of a sub_grid.
- **corner3:** is a **required** 2 member array of reals if using (grid_type="latlon") or integers if using (grid_type="indices") defining the x and y points of the third corner of a sub_grid.
- **corner4:** is a **required** 2 member array of reals if using (grid_type="latlon") or integers if using (grid_type="indices") defining the x and y points of the fourth corner of a sub_grid.
- **tile:** is an integer defining the tile number the sub_grid is on. It is **required** only if using (grid_type="indices").

**Exampe:**

```yaml
  sub_region:
  - grid_type: latlon
    corner1: -80, 0
    corner2: -80, 75
    corner3: -60, 0
    corner4: -60, 75
```

### 3. More examples
Bellow is a complete example of diag_table.yaml:
```yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: wild_card_name%4yr%2mo%2dy%2hr
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6
  new_file_freq_units: hours
  start_time: 2 1 1 0 0 0
  file_duration: 12
  file_duration_units: hours
  varlist:
  - module: test_diag_manager_mod
    var_name: sst
    reduction: average
    kind: r4
  global_meta:
  - is_a_file: true
- file_name: normal
  freq: 24
  freq_units: days
  time_units: hours
  unlimdim: records
  varlist:
  - module: test_diag_manager_mod
    var_name: sst
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
  freq_units: days
  time_units: hours
  unlimdim: records
  write_file: true
  varlist:
  - module: test_diag_manager_mod
    var_name: sstt
    reduction: average
    kind: r4
    long_name: S S T
  - module: test_diag_manager_mod
    var_name: sstt2
    reduction: average
    kind: r4
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
  freq_units: days
  time_units: hours
  unlimdim: records
  write_file: false
```
