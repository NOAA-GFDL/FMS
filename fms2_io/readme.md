@ingroup fms2_io_mod

### FMS_io/mpp_io to FMS2_io conversion guide

Before introducing the FMS2_io module, subroutines and functions in fms_io and mpp_io modules were used to mainly handle read/writes to NetCDF files. However, there were duplicate routines present in both modules that led to redundancy, and the blackbox-like I/O processes restricted user flexibility. The FMS2_io module has thus been implemented for a cleaner set of I/O tools and to give users more control over the information being written/read to NetCDF files. This guide helps convert fms_io/mpp_io code to FMS2_io

### A. FMS2_io Fileobjs
FMS2_io provides three new derived types, which target the different I/O paradigms used in GFDL models.

**1. FmsNetcdfFile_t:** This type provides a thin wrapper over the netCDF4 library, but allows the user to assign a “pelist” to the file. If a pelist is assigned, only the first rank on the list directly interacts with the netCDF library, and performs broadcasts to relay the information read to the rest of the ranks on the list. When writing netcdf files, only the first rank in the pelist will perform the writes.

**2. FmsNetcdfDomainFile_t:** This type does everything that the FmsNetcdfFile_t type does and it adds support for “domain-decomposed” reads/writes. Here, "domain decomposed" refers to data that is on a user-defined mpp_domain and is decomposed in two dimensions, in which each MPI rank has its own section of the global data. This requires a [domain](https://github.com/NOAA-GFDL/FMS/blob/main/mpp/mpp_domains.F90#L379-L415) to be associated with the fileobj.

**3. FmsNetcdfUnstructuredDomainFile_t:** This type does everything that the FmsNetcdfFile_t type does and it adds support for “domain-decomposed” reads/writes on a user defined mpp_domains **unstructured** grid. This requires a [unstructured domain](https://github.com/NOAA-GFDL/FMS/blob/3329625ea48bc3a10a5726c9f251d6d47b33516d/mpp/mpp_domains.F90#L267-L284) to be associated with the fileobj.

*The FMS_io equivalent to these derived types is [restart_file_type](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms/fms_io.F90#L297-L323)*

### B. Writing Restarts

#### 1. Domain Decomposed Restarts

When writing domain decomposed restarts, FMS sets the filename for domain decomposed restarts in the following way:
- If the domain is a cubesphere (6 tiles) and the io_layout is (/1,1/), it will create restart files: "RESTART/filename.tile{tile_number}.res.nc"
- If the domain is a cubesphere (6 tiles) and the io_layout is (/x_layout,y_layout/), it write `x_layout*y_layout` restart files: "RESTART/filename.tile{tile_number}.res.nc.????"
- If the domain has only 1 tile and the tile number is 1, "tile{tile_number}" will not be added

Additionally, the user can append an appendix (i.e "nest{nest_number}" or ens{ensemble_number}") to the filename by calling `set_filename_appendix (appendix)`
- If the domain has more than 1 tile, the appendix will be added as "RESTART/filename.appendix.tile{tile_number}.res.nc"
- If the domain has 1 tile, the appendix will added as "RESTART/filename.res.appendix.nc"

**This is how a restart file can be written with FMS_io:**
```F90

use fms_io_mod,         only: restart_file_type, register_restart_field, save_restart
use mpp_domains_mod,    only: domain

integer                 :: id_restart     !< Id for restart variable
type(restart_file_type) :: fileobj        !< Fms_io fileobj
real, dimension(:,:,:)  :: variable_data  !< Variable data in the compute or data domain
type (domain2d)         :: domain         !< 2d mpp domain

id_restart = register_restart_field(fileobj, "filename", "variable_name", variable_data, domain=domain)
call save_restart(Atm_restart)
```
Metadata:
- FMS_io named the dimensions: "xaxis_1", "yaxis_1", zaxis_1" and "Time". With FMS2_io the user can name the axis whatever they like.
- FMS_io wrote the dimensions as variables as well. FMS2_io does not do this by default, the user can do this if they like.
- FMS_io wrote variable attribute: "longname = {same as variable} and "units = {"none"} to all variables by default. FMS2_io does not do this by default, the user can add real meta data if they like.

**This is how a restart file can be written with FMS2-io:**
```F90

use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_domains_mod,    only: domain2d, center

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data in the compute or data domain
type (domain2d)             :: domain         !< 2d mpp domain
character(len=8)            :: dim_names(4)   !< Array of dimension names

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

!> Create domain !<
!> Create an io_domain !<

if (open_file(fileobj, "filename", "overwrite", domain, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1), "x", position=center)
  call register_axis(fileobj, dim_names(2), "y", position=center)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)

  call register_restart_field(fileobj, 'variable_name', variable_data, dim_names)
  call write_restart(fileobj)
  call close_file(fileobj)
endif
```
- [open_file](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_domain_io.F90#L321)
  -  a logical function, outputs .true. if the file was opened successfully and .false. if it failed.
  -  Mangles the filename in the same manner as fms_io
  -  Opens the netcdf file to write (a `nf90_create` call)
  -  Set ups the pelist for each io_domain
  -  `is_restart` indicates that this is a restart file, so it adds ".res" to the filename and it allows user to use the `write_restart` and `register_restart_field` functionality
  -  **NOTE**: The filename needs to include the full path to the file, including the directory.
- [register_axis](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_domain_io.F90#L437)
  - writes the dimension metadata in the netcdf file (a `nf90_def_dim` call)
  - The "x" and "y" argument indicate that that dimension is domain decomposed in x/y. The only acceptable values are "x" and "y".
  - The `position=center` indicates the position of the axis (this is the default). The other acceptable values are position=east for "x" andposition=north for "y", in this cases the data is staggered, which may be on the volume face or corner instead of at the centroid.
  - The "unlimited" indicates that the dimension is unlimited (`nf90_unlimited`)
  - The integer "dimsize" indicates that this is not a domain decomposed dimension with a length equal to dimsize
- [register_restart_field](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/register_domain_restart_variable.inc)
  - Writes the variable metadata to the file (a `nf90_def_var` call)
  - Saves the data as pointers, which will be written to the netcdf file later
- [write_restart](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_domain_io.F90#L549)
  - Loops through the restart variables that were registered
  - Calculates and writes a global checksum for each variables
  - Writes the data to the file (a `nf90_put_var` call)
- [close_file](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_domain_io.F90#L424)
  - Cleans up the fileobj (deallocates any arrays)
  - Closes the netcdf file (a `nf90_close` call)

#### 2. Unstructured Domain Restarts
Restart files with domain decomposed variables in the unstructured domain can be written using the `FmsNetcdfUnstructuredDomainFile_t` fileobj.

```F90
use fms2_io_mod,        only: FmsNetcdfUnstructuredDomainFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_domains_mod,    only: domainug, center

type(FmsNetcdfUnstructuredDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
real, dimension(:,:)        :: variable_data  !< Variable data in the unstructured domain
type (domainug)             :: domain         !< 2d mpp domain
character(len=8)            :: dim_names(3)   !< Array of dimension names

dim_names(1) = "ucdim"
dim_names(2) = "zaxis_1"
dim_names(3) = "Time"

!> Create domain !<
!> Create an io_domain !<
!> Create an unstructured domain !<

if (open_file(fileobj, "filename", "overwrite", domain, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1))
  call register_axis(fileobj, dim_names(2), dimsize)
  call register_axis(fileobj, dim_names(3), unlimited)

  call register_restart_field(fileobj, 'variable_name', variable_data, dim_names)
  call write_restart(fileobj)
  call close_file(fileobj)
endif
```
- [open_file](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/fms_netcdf_unstructured_domain_io.F90#L69-L70)
  -  a logical function, outputs .true. if the file was opened successfully and .false. if it failed.
  -  Mangles the filename: (1) Adds ".tileXX" if you are on multiple tiles (2) Adds .XXXX to the end of the file if uncombined.
  -  Opens the netcdf file to write (a `nf90_create` call)
  -  Set ups the pelist for each io_domain
  -  `is_restart` indicates that this is a restart file, so it adds ".res" to the filename and it allows user to use the `write_restart` and `register_restart_field` functionality
  -  **NOTE**: The filename needs to include the full path to the file, including the directory.
- [register_axis](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_unstructured_domain_io.F90#L155)
  - writes the dimension metadata in the netcdf file (a `nf90_def_dim` call)
  - The lack of the third argument in the `register_axis` calls indicates that this dimension is in the unstructured domain
  - The "unlimited" indicates that the dimension is unlimited (`nf90_unlimited`)
  - The integer "dimsize" indicates that this is a normal dimension of length equal to dimsize
- [register_restart_field](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/register_unstructured_domain_restart_variable.inc)
  - Writes the variable metadata to the file (a `nf90_def_var` call)
  - Saves the data as pointers, which will be written to the netcdf file later
- [write_restart](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_unstructured_domain_io.F90#L193)
  - Loops through the restart variables that were registered
  - Calculates and writes a global checksum for each variables
  - Writes the data to the file (a `nf90_put_var` call)
- [close_file](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_unstructured_domain_io.F90#L146)
  - Cleans up the fileobj and deallocated any variables
  - Closes the netcdf file (a `nf90_close` call)

#### 3. Non-domain Decomposed Restarts
Restart files without domain decomposed variables can be written using the `FmsNetcdfFile_t` fileobj.

```F90

use fms2_io_mod,        only: FmsNetcdfFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_mod,            only: mpp_npes, mpp_get_current_pelist

type(FmsNetcdfFile_t)       :: fileobj        !< Fms2_io fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data
character(len=8)            :: dim_names(4)   !< Array of dimension names
integer, allocatable,       :: pes(:)         !< Array of the pes in the current pelist

!< Get the current pelist
allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

if (open_file(fileobj, "filename", "overwrite", pelist=pes, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1), 96)
  call register_axis(fileobj, dim_names(2), 96)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)

  call register_restart_field(fileobj, 'variable_name', variable_data, dim_names)
  call write_restart(fileobj)
  call close_file(fileobj)
endif
```

- [open_file](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1845)
  -  a logical function, outputs .true. if the file was opened successfully and .false. if it failed.
  -  Mangles the filename: ".res" is the only thing that is added to the filename
  -  Opens the netcdf file to write (a `nf90_create` call)
  -  `is_restart` indicates that this is a restart file, so it adds ".res" to the filename and it allows user to use the `write_restart` and `register_restart_field` functionality
  -  With the *optional* pelist argument, only the first rank interacts with the file (opens, writes, closes) and broadcasts the information to the rest of the ranks on the list.
  -  **NOTE**: The filename needs to include the full path to the file, including the directory.
- [register_axis](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/netcdf_io.F90#L727-L728)
  - writes the dimension metadata in the netcdf file (a `nf90_def_dim` call)
  - The "unlimited" indicates that the dimension is unlimited (`nf90_unlimited`)
  - The integer "dimsize" indicates that this is a normal dimension of length equal to dimsize
  - "x" and "y" cannot be added to the `register_axis` calls
- [register_restart_field](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms2_io.F90#L135-L140)
  - Writes the variable metadata to the file (a `nf90_def_var` call)
  - Saves the data as pointers, which will be written to the netcdf file later
- [write_restart](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1898)
  - Loops through the restart variables that were registered
  - Calculates and writes a global checksum for each variables
  - Writes the data to the file (a `nf90_put_var` call)
- [close_file](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1876)
  - Cleans up the fileobj
  - Closes the netcdf file (a `nf90_close` call)

#### 4. Other helpful information:
- By default, FMS2_io writes the file as `nf90_64bit_offset`, the user can change the netcdf file type by adding `nc_format="64bit", "classic", or "netcdf4"` to the `open_file` call (which will change the type for that file only) or by adding `netcdf_default_format="64bit", "classic", or "netcdf4"` to the fms2_io_nml (which will change the type for all files that don't have the nc_format in the open_file call)
- If the user wishes to not add ".res" to filename, the user can add `dont_add_res_to_filename=.true.` to the `open_file` call
- **Variable attributes** can be written by calling [register_variable_attribute](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/register_variable_attribute.inc). Scalar and 1d real and integers (32 and 64 bit) and string values are supported
```F90
call register_variable_attribute(fileobj, "varname", "attribute_name", value)
```
This interface can be used with any FMS2_io fileobj, but the open_file needs to be called before using it.

*A similar thing can be accomplished in fms_io with [mpp_write_meta](https://github.com/NOAA-GFDL/FMS/blob/main/mpp/include/mpp_io_write.inc#L22-L135)*

- **Global attributes** can be written by calling [register_global_attribute](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/register_global_attribute.inc). Scalar and 1d real and integers (32 and 64 bit) and scalar string values are supported
```F90
call register_global_attribute(fileobj, "global_attribute_name", value)
```
This interface can be used with any FMS2_io fileobj, but the open_file needs to be called before using it.

*A similar thing can be accomplished in fms_io with [mpp_write_meta](https://github.com/NOAA-GFDL/FMS/blob/main/mpp/include/mpp_io_write.inc#L22-L135)*

### C. Reading Restarts
The restarts can be read the same way as the writes. The only difference is that "read" is used in the `open_file` call and [read_restart](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/fms2_io.F90#L216-L219) is used instead of write restart.

Additionally, the `register_restart_field` has the optional argument, `is_optional`, where is .true. the code will not crash if the variable does not exist in the file when reading. *This argument is the same as `mandatory` in fms_io.

#### Other helpful information:
The following subroutines can be used with any of any of the FMS2_io fileobj.
- **Variable attributes** can be read by calling [get_variable_attribute](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/get_variable_attribute.inc). Scalar and 1d real and integers (32 and 64 bit) and string values are supported. To check if a variable attribute exists before reading, the logical function [variable_att_exists](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/netcdf_io.F90#L1046-L1047) can be used
```F90
if (variable_att_exists(fileobj, "varname", "attribute_name")) then
   call get_variable_attribute(fileobj, "varname", "attribute_name", value)
endif
```
*A similar thing can be accomplished in fms_io with [get_vart_att_value](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms/fms_io.F90#L466-L468)*
- **Global attributes** can be read by calling [get_global_attribute](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/get_global_attribute.inc). Scalar and 1d real and integers (32 and 64 bit) and scalar string values are supported. To check if a global attribute exist before reading, the logical function [global_att_exists](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1019) can be used.
```F90
if (global_att_exists(fileobj, "global_attribute_name)) then
   call get_global_attribute(fileobj, "global_attribute_name", value)
endif
```
*A similar thing can be accomplished in fms_io with [get_global_att_value](https://github.com/NOAA-GFDL/FMS/blob/main/fms/fms_io.F90#L461-L464)*

- **Reading dimension metadata**
  - [get_num_dimensions](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1079) can be used to get the total number of dimensions in a file
  - [get_dimension_names](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1107) can be used to get the dimension names in a file
  - [dimension_exists](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1161) is a logical function that can be used to check if a dimension exists
  - [get_unlimited_dimension_name](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1229) can be used to get the name of the unlimited dimension in a file
  - [is_dimension_unlimited](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1196) is a logical function that can be used to determine if a dimension is unlimited
  - [get_dimension_size](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1264) can be used to get the size of a dimension

 Sample usage:
 ```F90
 ndims_file = get_num_dimensions(fileobj)
 allocate(dim_names_file(ndims_file))
 call get_dimension_names(fileobj, dim_names_file)

 call get_unlimited_dimension_name(fileobj, dimension_name)
 if (is_dimension_unlimited(fileobj, "dimension_name")) call mpp_error(FATAL, "dimension_name is not an unlimited dimension")

 if (dimension_exists(fileobj, "dimension_name")) then
     call get_dimension_size(fileobj, "dimension_name", dimsize)
 endif
 ```
 *A similar thing can be accomplished with fms_io with [mpp_get_info](https://github.com/NOAA-GFDL/FMS/blob/main/mpp/include/mpp_io_util.inc#L41)*

- **Reading variable metadata**
  - [get_num_variables](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1295) can be used to get the number of variables in a file
  - [get_variable_names](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1323) can be used to get the variables names in a file
  - [variable_exists](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1377) is a logical function that can be used to check if a variable exists
  - [get_variable_num_dimensions](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1408) can be used to get the number of dimensions in a variable
  - [get_variable_size](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1500) can be used to get the size of the variable.
  - [get_variable_dimension_names](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L1439) can be used to get the names of the dimensions in a variable

Sample usage:
```F90
nvar = get_num_variables(fileobj)

allocate(var_names(nvar))
call get_variable_names(fileobj, var_names)

if (variable_exists(fileobj, "varname")) then
   ndims = get_variable_num_dimensions(fileobj, "varname")

   allocate(dim_sizes(ndims))
   call get_variable_size(fileobj, "varname", dim_sizes)

   allocate(dim_names(ndims))
   call get_variable_dimension_names(fileobj, "varname", dim_names)
endif

```
 *A similar thing can be accomplished with fms_io with [mpp_get_info](https://github.com/NOAA-GFDL/FMS/blob/main/mpp/include/mpp_io_util.inc#L41) and/or [field_size](https://github.com/NOAA-GFDL/FMS/blob/main/fms/fms_io.F90#L4965)*
### D. Reading/Writting Non-restarts
Reading and writing netcdf files that are not restarts can be done using `read_data` and `write_data` calls.

#### 1. Domain decomposed read/write:
```F90

use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_data
use mpp_domains_mod,    only: domain2d, center

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data in the compute or data domain
type (domain2d)             :: domain         !< 2d mpp domain
character(len=8)            :: dim_names(4)   !< Array of dimension names

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

if (open_file(fileobj, "filename", "overwrite", domain)) then
  call register_axis(fileobj, dim_names(1), "x", position=center)
  call register_axis(fileobj, dim_names(2), "y", position=center)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)

  call register_field(fileobj, 'variable_name', 'variable_type', dim_names)
  call write_data(fileobj, 'variable_name', variable_data)
  call close_file(fileobj)
endif
```
Difference from writing restarts:
- `open_file` does not have is_restart=.true., so the `write_restart` functionality cannot be used'
- [register_field](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/fms_netcdf_domain_io.F90#L527) is basically a wrapper for `nf90_def_var` which just adds the variable metadata to the file.
  - "variable_type" is a string which indicates the type you want the variable to be written as. The acceptable values are "int", "int64", "double", "float", and "char"
- [write_data](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/domain_write.inc)
  - The ranks send their data to the io_root pe. The io_root pe receives the data and writes it to the file.

Similarly for domain reads:
```F90
if (open_file(fileobj, "filename", "read", domain)) then
  call register_axis(fileobj, dim_names(1), "x", position=center)
  call register_axis(fileobj, dim_names(2), "y", position=center)

  call read_data(fileobj, 'variable_name', variable_data)
  call close_file(fileobj)
endif
```
- `register_axis` is required for the domain decomposed dimensions so the code knows which dimensions and therefore variables are domain decomposed.
- `register_field` is not required because the code can get the dimensions information from the file.
- [read_data](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/domain_read.inc)
  - The io_root pe reads the data and sends it to the other pes.

#### 2. Unstructured Domain non-restart read/writes

```F90
use fms2_io_mod,        only: FmsNetcdfUnstructuredDomainFile_t, register_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_data
use mpp_domains_mod,    only: domainug, center

type(FmsNetcdfUnstructuredDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
real, dimension(:,:)        :: variable_data  !< Variable data in the unstructured domain
type (domainug)             :: domain         !< 2d mpp domain
character(len=8)            :: dim_names(3)   !< Array of dimension names

dim_names(1) = "ucdim"
dim_names(2) = "zaxis_1"
dim_names(3) = "Time"

!> Create domain !<
!> Create an io_domain !<
!> Create an unstructured domain !<

if (open_file(fileobj, "filename", "overwrite", domain, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1))
  call register_axis(fileobj, dim_names(2), dimsize)
  call register_axis(fileobj, dim_names(3), unlimited)

  call register_field(fileobj, 'variable_name', 'variable_type', dim_names)
  call write_data(fileobj, 'variable_name', variable_data)
  call close_file(fileobj)
endif
```
Difference from writing restarts:
- `open_file` does not have is_restart=.true., so the `write_restart` functionality cannot be used'
- [register_field](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/fms_netcdf_unstructured_domain_io.F90#L178-L179) is basically a wrapper for `nf90_def_var` which just adds the variable metadata to the file.
  - "variable_type" is a string which indicates the type you want the variable to be written as. The acceptable values are "int", "int64", "double", "float", and "char"
- [write_data](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/unstructured_domain_write.inc)
  - The ranks send their data to the io_root pe. The io_root pe receives the data and writes it to the file.

Similarly for reads:
```F90
if (open_file(fileobj, "filename", "read", domain, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1))
  call write_data(fileobj, 'variable_name', variable_data)
  call close_file(fileobj)
endif
```
- `register_axis` is required for the dimensions that are on the unstructured domain, so the code knows which dimensions and therefore variables are in the unstructured domain
- `register_field` is not required because the code can get the dimensions information from the file.
- [read_data](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/unstructured_domain_read.inc)
  - The io_root pe reads the data and sends it to the other pes.

#### 3. Non-domain decomposed read/writes

Writing non-domain decomposed files is similar to writing domain decomposed files:

```F90
use fms2_io_mod,        only: FmsNetcdfFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_mod,            only: mpp_npes, mpp_get_current_pelist

type(FmsNetcdfFile_t)       :: fileobj        !< Fms2_io fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data
character(len=8)            :: dim_names(4)   !< Array of dimension names
integer, allocatable,       :: pes(:)         !< Array of the pes in the current pelist

!< Get the current pelist
allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

if (open_file(fileobj, "filename", "overwrite", pelist=pes, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1), 96)
  call register_axis(fileobj, dim_names(2), 96)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)

  call register_field(fileobj, 'variable_name', 'variable_type', dim_names)
  call write_data(fileobj)
  call close_file(fileobj)
endif
```
- `open_file` - the pelist argument is needed so that only the root pe writes the file
- `register_axis` is needed so that the code knows the size of each dimension
- `register_field` is needed so that the code knows what dimensions the variable is a function of
- [write_data](https://github.com/NOAA-GFDL/FMS/blob/main/fms2_io/include/netcdf_read_data.inc): the root pe in the pelist writes the data

To read non-domain decomposed restarts:

```F90
if (open_file(fileobj, "filename", "read", pelist=pes, is_restart=.true.)) then
  call read_data(fileobj, 'variable_name', variable_data)
  call close_file(fileobj)
endif
```
- Here the `register_axis` and `register_field` calls are not needed because the code can get the dimensions from the file

Additionally, the `corner` and `edge_lengths` optional arguments in the `read_data` calls can be used to read only a section of a file.
```F90
if (open_file(fileobj, "filename", "read", pelist=pes, is_restart=.true.)) then
  call read_data(fileobj, 'variable_name', variable_data, corner=(/10,2,2,3/), edge_lengths=(/2, 3, 4, 1/))
  call close_file(fileobj)
endif
```
Here the code is going to start reading at x=10, y=2, z=2, t=3 and read 2, 3, 4, 1 points in each direction. Note that `variable_data` must be the correct size or the code will fail.

Additionally, it is possible to do "compressed" reads and writes. Here "compressed" is when different ranks have different numbers of points for a specific
axis. For example:
```F90
if (open_file(fileobj, "filename", "overwrite", pelist=pes, is_restart=.true.)) then
  call register_axis(fileobj, "xaxis_1", mpp_pe()+1, is_compressed=.true.)
  call register_axis(fileobj, "yaxis_1", mpp_pe()+1, is_compressed=.true.)

  call register_field(fileobj, 'variable_name', 'variable_type', dim_names)
  call write_data(fileobj)
  call close_file(fileobj)
endif
```
- Because [is_compressed=.true.](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L747-L768), the root pe is going to gather the dimension size of each rank and add them to get the total length of the dimension. In the example above, rank 0 has a x/y axis of 1, rank 1 has a x/y of 2. The total dimension size of x/y is equal to 3 (the sum for all pes).

### E. Coupler Type Restarts
`coupler_type_register_restarts` permits interating a list of variables in the coupler_type and register each variable to the corresponding restart file and then read/write the file.

#### 1. Reading coupler_type restarts
In FMS_io, this was accomplished as:

```F90
type(coupler_2d_bc_type)              :: bc_type
type(restart_file_type), pointer  :: bc_rest_files(:)=> null() !< Array of fms2_io fileobjs

call coupler_type_register_restarts(bc_type, bc_rest_files, num_rest_files, domain)
do l = 1, num_rest_files
   call restore_state(bc_rest_files(l), directory='INPUT', &
                           nonfatal_missing_files=.true.)
enddo
```

In FMS2_io, this can be accomplished as:

```F90
type(coupler_2d_bc_type)              :: bc_type
type(FmsNetcdfDomainFile_t), pointer  :: bc_rest_files(:)=> null() !< Array of fms2_io fileobjs

call coupler_type_register_restarts(bc_type, bc_rest_files, num_rest_files, domain, to_read=.true.)

do l = 1, num_rest_files
   call read_restart(bc_rest_files(l))
   call close_file(bc_rest_files(l))
enddo
```
Note: `num_rest_files` is set inside `coupler_type_register_restarts`

Difference from fms_io
- `to_read` was added to the coupler_type_register_restarts subroutine
-  [coupler_type_register_restarts](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/coupler/coupler_types.F90#L3096) loops through the variables in the coupler_type, opens the file needed, register the axis, and registers the restart_variables.
- The file should be closed after writing it by calling `close_file`
- It is required that the `domain` has an io_domain or the code will fail

#### 2. Writting coupler type restarts
This is done the same way as the reads expect `write_restart` is used instead of `read_restart`.

### F. Boundary Conditions Restarts
Both FMS_io and FMS2_io have the functionality where one can read or write data where only some pes have a section of the data (i.e halos).

In FMS_io this was accomplished as:
```F90
type(restart_file_type) :: fileobj        !< Fms_io fileobj

id_restart = register_restart_field(fileobj, "filename", "variable_name", variable_data, indices, global_size, pelist, is_root_pe)

call save_restart_border(fileobj)
call free_restart_type(fileobj)

```

In FMS2_io this can be accomplished as:
```F90
type(FmsNetcdfFile_t)       :: fileobj        !< Fms2_io fileobj

if (open_file(fileobj, "filename", "overwrite", is_restart=.true., pelist=all_pelist)) then
    call register_restart_field(fileobj, "variable_name", variable_data, indices, global_size, pelist, is_root_pe)
    call write_restart_bc(fileobj)
    call close_file(fileobj)
endif
```

- `indices` are the starting and indices of the region for that pe (starting x, ending x, starting y, ending y)
- `global_size` is the size of the variable in (x,y,*z*)
- `pelist` is the list of pelist that have the data for the variable
- `is_root_pe` is a flag indicating if this is the root_pe from the pelist

Difference between FMS2_io and FMS_io:
  - In fms_io, `register_restart_field` is a integer function while in FMS2_io it is a subroutine
  - In fms_io, the first two arguments of `register_restart_field` were [fileobj and filename](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms/fms_io.F90#L2376-L2377) and in FMS2_io it is just the [fileobj](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/include/netcdf_add_restart_variable.inc#L350)
  - In FMS2_io, [restore_state_border](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms/fms_io.F90#L3668) was renamed to renamed to [read_restart_bc](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L2083)
  - In FMS2_io, [save_restart_border](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms/fms_io.F90#L3378) was renamed to [write_restart_bc](https://github.com/NOAA-GFDL/FMS/blob/b9fc6515c7e729909e59a0f9a1efc6eb1d3e44d1/fms2_io/netcdf_io.F90#L2114)
  - In FMS2_io, it is required for the file to be opened before registering the restarts.

### G. Ascii_io

#### 1. Writing ascii files
With fms_io, ascii writes can be accomplished with a combination of `mpp_open`, `mpp_close` from `mpp_io` and `write`.
```F90
integer :: unit_number
call mpp_open( unit_number, 'filename', nohdrs=.TRUE. )
if ( mpp_pe().EQ.mpp_root_pe() ) write (unit_number, *) 'FMS2_io is flawless and works as inteded'
call mpp_close(unit_number)
```

With FMS2_io, there is no real support for ascii writes. The user can use fortran's `open`, `write` and `close`.  For example,
```F90
integer :: unit_number

if ( mpp_pe().EQ.mpp_root_pe() ) then
   open(newunit = unit_number, file='filename', status='replace', form='formatted')
   write (unit_number, *) 'FMS2_io is flawless and works as inteded'
   close (unit_number)
endif
```
- The `if ( mpp_pe().EQ.mpp_root_pe() ) ` ensures that the file is only written once.
- The `newunit` ensures that an unused unit number is used.

#### 2. Reading ascii files
With fms_io:
- Opening an ascii file was accomplished with `mpp_open`, `open_namelist_file` or `open` calls.
- Reading an ascii file was accomplished with `read_data`, `read_distributed` or `read` calls.
- Closing an ascii file was sometimes done with `mpp_close`, `close_file` or `close` calls.

With FMS2_io, ascii reads can be accomplished using `ascii_reads`. In `ascii_reads`, the root_pe opens and reads the file into a a string buffer and broadcasts it to the other ranks. Each rank can then read from the buffer. For example:
```F90
character(len=:), dimension(:), allocatable :: restart_file !< Restart file saved as a string

call ascii_read('INPUT/coupler.res', restart_file)
read(restart_file(1), *) calendar_type
read(restart_file(2), *) date_init
read(restart_file(3), *) date
deallocate(restart_file)
```

### G. FMS2_io namelist
- **ncchksz:** Sets chunksize (in bytes) argument in netcdf file creation calls. The default is `64*1024`.
- **netcdf_default_format:** Sets the netcdf file type. The acceptable values are  "64bit", "classic", "netcdf4". This can be overwritten per file if you specify `nc_format` in the open_file call. The default is 64bit.
- **header_buffer_val:** Sets the netCDF header buffer size(in bytes). The default is 16384 bytes.
