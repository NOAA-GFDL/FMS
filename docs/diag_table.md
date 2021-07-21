@ingroup diag_manager
# The Diagnostic Manager Table

The _diagnostics table_ (`diag_table`) allows users to specify output fields and
sampling rates at run time.  The `diag_table` file consists of comma-separated
ASCII values.  The `diag_table` has three sections: **Global**, **File**, and
**Field** sections.  The **Global** section must be the first two lines of the
file, whereas the **File** and **Field** sections can be inter mixed to allow
the file to be organized as desired. Comments can be added to the `diag_table`
file by using the hash symbol (`#`) as the first character in the line.

All errors in the `diag_table` will throw a FMS `FATAL` error.  The [FRE NCtools
suite](https://github.com/NOAA-GFDL/FRE-NCtools) contains the utility,
`diag_table_chk` to check the `diag_table` for errors.  A brief usage statement
can be obtained by running `diag_table_chk --help`, and a more in-depth
description of the command and options can be viewed with the command `perldoc
diag_table_chk`.

## Global Section

The first two lines of the `diag_table` must contain a _title_ #and the _base
date_ of the experiment respectively.  The _title_ must be a Fortran
`CHARACTER` string.  The _base date_ is the reference time used for the #time
units, and must be greater than or equal to the model start time. The #_base
date_ consists of six space-separated integers in the following format.

```
year month day hour minute second
```

## File Section

File lines contain 6 required and 5 optional fields (optional fields are
surrounded with square brackets ([]).  File lines can be intermixed with the
field lines, but the file must be defined before any fields that are to be
written to the file.  File lines have the following format:<BR />

```
"file_name", output_freq, "output_freq_units", file_format, "time_axis_units", "time_axis_name" [, new_file_freq, "new_file_freq_units"[, "start_time"[, file_duration, "file_duration_units"]]]
```

with the following descriptions and Fortran types.

#### `CHARACTER(len=128) :: file_name`

Output file name without the trailing `.nc`.

A single file description can produce multiple files using special time string
suffix keywords.  This time string will append the time strings to the base file
name each time a new file is opened.  They syntax for the time string suffix
keywords are `%#tt` Where `#` is a mandatory single digit number specifying the
width of the field, and `tt` can be as follows:

`yr` &ndash; Years\
`mo` &ndash; Months\
`dy` &ndash; Days\
`hr` &ndash; Hours\
`mi` &ndash; Minutes\
`sc` &ndash; Seconds


Thus, a file name of `file2_yr_dy%1yr%3dy` will have a base file name of
`file2_yr_dy_1_001` if the file is created on year 1 day 1 of the model run.
**_NOTE:_** The time suffix keywords must be used if the optional fields
`new_file_freq` and `new_file_freq_units` are used, otherwise a FMS `FATAL`
error will occur.

#### `INTEGER :: output_freq`

How often to write fields to file.

`> 0` &ndash; Output frequency in `output_freq_units`.\
`= 0` &ndash; Output frequency every time set. (`output_freq_units` is ignored.)\
`=-1` &ndash; Output at end of run only. (`output_freq_units` is ignored.)


#### `CHARACTER(len=10) :: output_freq_units`

Time units for output.

The time units can be one of: `years`, `months`, `days`, `minutes`, `hours`, or
`seconds`.

#### `INTEGER :: file_format`

Output file format.

Currently only the _netCDF_ (value of `1`) file format is supported.

#### `CHARACTER(len=10) :: time_axis_units`

Time units for the output file time axis.

The time units can be one of: `years`, `months`, `days`, `minutes`, `hours`, or
`seconds`.

#### `CHARACTER(len=128) :: time_axis_name`

Axis name for the output file time axis.

The character sting must contain the case insensitive string `"time"`.


#### `INTEGER, OPTIONAL :: new_file_freq`

Frequency for closing the existing file, and creating a new file in
`new_file_freq_units`.

#### `CHARACTER(len=10), OPTIONAL :: new_file_freq_units`

Time units for creating a new file.  Can be either `years`, `months`, `days`,
`minutes`, `hours`, or `seconds`.  **_NOTE:_** If the `new_file_freq` field is
present, then this field must also be present.

#### `CHARACTER(len=25), OPTIONAL :: start_time`

Time to start the file for the first time.  The format of this string is the
same as the _global date_.  **_NOTE:_** The `new_file_freq` and
the `new_file_freq_units` fields must be present to use this field.

#### `INTEGER, OPTIONAL :: file_duration`

How long file should receive data after start time in `file_duration_units`.
This optional field can only be used if the `start_time` field is present.  If
this field is absent, then the file duration will be equal to the frequency for
creating new files.  **_NOTE:_** The `file_duration_units` field must
also be present if this field is present.

#### `CHARACTER(len=10), OPTIONAL :: file_duration_units`

File duration units.

Units an be one of `years`, `months`, `days`, `minutes`,
`hours`, or `seconds`.  **_NOTE:_** If the `file_duration` field is
present, then this field must also be present.

## Field Section

Field lines contain 8 fields.  Field lines can be intermixed with file lines.
Fields line can contain fields that are not written to any files.  The file name
for these fields is `null`.

Field lines have the following format:

```
"module_name", "field_name", "output_name", "file_name", "time_sampling", "reduction_method", "regional_section", packing
```

with the following descriptions.

#### `CHARACTER(len=128) :: module_name`

Module that contains the `field_name` variable.  (e.g. `atmos_mod`, `land_mod`)

#### `CHARACTER(len=128) :: field_name`

Module variable name that has data to be written to file.

#### `CHARACTER(len=128) :: output_name`

Name of the field as written in `file_name`.

#### `CHARACTER(len=128) :: file_name`

Name of the file where the field is to be written. **_Note:_** The file
`file_name` must be defined first.

#### `CHARACTER(len=50) :: time_sampling`

Currently not used.  Please use the string `"all"`.

#### `CHARACTER(len=50) :: reduction_method`

The data reduction method to perform prior to writing data to disk.  Valid
options are (redundant names are separated with commas):

##### `.TRUE.`, `"average"`, `"avg"`, `"mean"`

Average from the last time written to the current time.

##### `.FALSE.`, `"none"`

No reduction performed.  Write current time step value only.

##### `"rms"`
Calculate the root mean square from the last time written to the current time.

##### `"pow##"`
Calculate the mean of the power `##` from the last time written to the current time.

##### `"min"`

Minimum value from last write to current time.

##### `"max"`

Maximum value from last write to current time.

##### `"diurnal##"`

`##` diurnal averages

#### `CHARACTER(len=50) :: regional_section`

Bounds of the regional section to capture.  A value of `none` indicates a global
region.  The regional section has the following format:

```
lat_min, lat_max, lon_min, lon_max, vert_min, vert_max
```

Use `vert_min = -1` and `vert_max = -1` to get the entire vertical axis.
**_Note:_** Currently, the defined region _MUST_ be confined to a single
tile.

#### `INTEGER :: packing`

Fortran number `KIND` of the data written.  Valid values:

`= 1` &ndash; double precision\
`= 2` &ndash; float\
`= 4` &ndash; packed 16-bit integers\
`= 8` &ndash; packed 1-byte (not tested).

## Sample Diagnostic Manager Table

```
"diag manager test"
1999 1 1 0 0 0
# End the Global Section

# Output Files
"10_days",             10, "days", 1, "hours", "Time"
"file1_hr%hr3",         5, "days", 1, "hours", "Time", 15, "days"
"file2_yr_dy%yr1%dy3",  5, "days", 1, "hours", "Time", 10, "days", "1 1 7 0 0 0"
"file3_yr_dy%yr1%dy3",  5, "days", 1, "hours", "Time", 20, "days", "1 1 7 0 0 0", 5, "years"

# Output Variable
"ice_mod", "ice", "ice", "10_days", "all", .false., "none", 2

# File and Field
temp_local, 1, "days", 1, "hours", "Time"
"ocean_mod", "temp", "temp", "temp_local", "all", .FALSE., "5 259.5 -59.5 59.5 1 1", 2
```
