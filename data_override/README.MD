@ingroup data_override_mod

# Data_override

### Contents
- [YAML Data Table format](README.MD#1-yaml-data-table-format)
- [How to use it?](README.MD#2-how-to-use-it)
- [Converting legacy data_table to data_table.yaml](README.MD#3-converting-legacy-data_table-to-data_tableyaml)
- [Examples](README.MD#4-examples)
- [External Weight File Structure](README.MD#5-external-weight-file-structure)

#### 1. YAML Data Table format:
Each entry in the data_table has the following key values:
- **grid_name:** Name of the grid to interpolate the data to. The acceptable values are "ICE", "OCN", "ATM", and "LND"
- **fieldname_in_model:** Name of the field as it is in the code to interpolate.
- **override_file:** Optional subsection with key/value pairs defining how to override from a netcdf file.
  - **file_name:** Name of the file where the variable is located, including the directory
  - **fieldname_in_file:** Name of the field as it is writen in the file
  - **interp_method:** Method used to interpolate the field. The acceptable values are "bilinear", "bicubic", and "none". "none" implies that the field in the file is already in the model grid. The LIMA format is no longer supported
  - **multi_file:** Optional subsection with key/value pairs to use multiple(3) input netcdf files instead of 1. Note that **file_name** must be the second file in the set when using multiple input netcdf files
    - **prev_file_name:** The name of the first file in the set
    - **next_file_name:** The name of the third file in the set
  - **external_weights:** Optional subsection with key/value pairs defining the external weights file to used for the interpolation.
    - **file_name:** Name of the file where the external weights are located, including the directory
    - **source:** Name of the source that generated the external weights. The only acceptable value is "fregrid"
- **factor:** Factor that will be multiplied after the data is interpolated
- **subregion:** Optional subsection with key/value pairs that define a subregion of the model grid to interpolate the data to.
  - **type:** The region type. The acceptable values are "inside_region" and "outside_region"
  - **lon_start:** The starting latitude in the same units as the grid data in the file
  - **lon_end:** The ending latitude in the same units as the grid data in the file
  - **lat_start:** The starting longitude in the same units as the grid data in the file
  - **lon_end:** The ending longitude in the same units as the grid data in the file

#### 2. How to use it?
In order to use the yaml data format, [libyaml](https://github.com/yaml/libyaml) needs to be installed and linked with FMS. Additionally, FMS must be compiled with -Duse_yaml macro. If using autotools, you can add `--with-yaml`, which will add the macro for you and check that libyaml is linked correctly.
```
./configure --with-yaml
```

If using cmake, you can use `-DWITH_YAML`.
```
cmake -DWITH_YAML=on;
```
#### 3. Converting legacy data_table to data_table.yaml
To convert legacy data_table to data_table.yaml, the python script, [`data_table_to_yaml.py`](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/main/data_table/data_table_to_yaml.py), can be used. Additionally, the python script, [`is_valid_data_table_yaml.py`](https://github.com/NOAA-GFDL/fms_yaml_tools/blob/main/data_table/is_valid_data_table_yaml.py), can be used to check if your data_table.yaml file is valid.

#### 4. Examples
**4.1** The following example is going to bilnearly interpolate `sic_obs` to the `ICE` grid using the `sic` data from the file: `INPUT/hadisst_ice.data.nc`, multiply by a factor of `0.01`, and interpolate to the model time.

In the **legacy format**, the data_table will look like:
```
"ICE", "sic_obs", "sic", "INPUT/hadisst_ice.data.nc", "bilinear", 0.01
```

In the **yaml format**, the data_table will look like
```
data_table:
 - grid_name           : ICE
   fieldname_in_model  : sic_obs
   override_file:
   - file_name         : INPUT/hadisst_ice.data.nc
     fieldname_in_file : sic
     interp_method     : bilinear
   factor              : 0.01
```
Which corresponds to the following model code:
```F90
call data_override('ICE', 'sic_obs', icec, Spec_Time)
```
where:
- `ICE` is the component domain for which the variable is being interpolated and corresponds to the grid_name in the data_table
- `sic_obs` corresponds to the fieldname_in_model in the data_table
- `icec` is the storage array that holds the interpolated data
- `Spec_Time` is the time to interpolate the data to.

Additionally, it is required to call data_override_init (in this case with the ICE domain). The grid_spec.nc file must also contain the coordinate information for the domain being used.

```F90
call data_override_init(Ice_domain_in=Ice_domain)
```

**4.2** The following example is going to simply multiply `sic_obs` by a factor of 0.01.

In the **legacy format**, the data_table will look like:
```
"ICE", "sit_obs", "", "INPUT/hadisst_ice.data.nc", "none",     2.0
```

In the **yaml format**, the data_table will look like:
``` yaml
data_table:
 - grid_name          : ICE
   fieldname_in_model : sit_obs
   factor             : 0.01
```

Which corresponds to the following model code:
```F90
call data_override('ICE', 'sit_obs', icec, Spec_Time)
```
where:
- `ICE` is the component domain for which the variable is being interpolated and corresponds to the grid_name in the data_table
- `sit_obs` corresponds to the fieldname_in_model in the data_table
- `icec` is the storage array that holds the interpolated data
- `Spec_Time` is the time to interpolate the data to.

Additionally, it is required to call data_override_init (in this case with the ICE domain). The grid_spec.nc file is still required to initialize data_override with the ICE domain.

```F90
call data_override_init(Ice_domain_in=Ice_domain)
```

**4.3** The following example is an `ongrid` case where it will simply read the `runoff` data in the `./INPUT/runoff.daitren.clim.nc` file, multiply it by a factor of `1.0` and interpolate to the model time. It assumes that "runoff" in the netcdf file is in the ocean grid.

In the **legacy format**, the data_table will look like:
```
"OCN", "runoff", "runoff", "./INPUT/runoff.daitren.clim.nc", "none" ,  1.0
```

In the **yaml format**, the data_table will look like:
``` yaml
data_table:
 - grid_name           : OCN
   fieldname_in_model  : runoff
   override_file:
   - file_name         : INPUT/runoff.daitren.clim.nc
     fieldname_in_file : runoff
     interp_method     : none
   factor              : 1.0
```

Which corresponds to the following model code:
```F90
call data_override('OCN', 'runoff', runoff_data, Spec_Time)
```
where:
- `OCN` is the component domain for which the variable is being interpolated and corresponds to the grid_name in the data_table
- `runoff` corresponds to the fieldname_in_model in the data_table
- `runoff_data` is the storage array that holds the interpolated data
- `Spec_Time` is the time to interpolate the data to.

Additionally, it is required to call data_override_init (in this case with the ocean domain). The grid_spec.nc file is still required to initialize data_override with the ocean domain and to determine if the data in the file is in the same grid as the ocean.

```F90
call data_override_init(Ocn_domain_in=Ocn_domain)
```

**4.4** The following example uses the multi-file capability
``` yaml
data_table:
 - grid_name           : ICE
   fieldname_in_model  : sic_obs
   override_file:
   - file_name         : INPUT/hadisst_ice.data_yr1.nc
     fieldname_in_file : sic
     interp_method     : bilinear
     multi_file:
     - next_file_name: INPUT/hadisst_ice.data_yr2.nc
       prev_file_name: INPUT/hadisst_ice.data_yr0.nc
   factor              : 0.01
```
Data override determines which file to use depending on the model time. This is to prevent having to combine the 3 yearly files into one, since the end of the previous file and the beginning of the next file are needed for yearly simulations.

**4.5** The following example uses the external weight file capability
``` yaml
data_table:
 - grid_name           : ICE
   fieldname_in_model  : sic_obs
   override_file:
   - file_name         : INPUT/hadisst_ice.data.nc
     fieldname_in_file : sic
     interp_method     : bilinear
     external_weights:
     - file_name: INPUT/remamp_file.nc
       source: fregrid
   factor              : 0.01
```

#### 5. External Weight File Structure

**5.1** Bilinear weight file example from fregrid

```
dimensions:
        nlon = 5 ;
        nlat = 6 ;
        three = 3 ;
        four = 4 ;
variables:
        int index(three, nlat, nlon) ;
        double weight(four, nlat, nlon) ;
```
- `nlon` and `nlat` must be equal to the size of the global domain.
- `index(1,:,:)` corresponds to the index (i) of the longitudes point in the data file, closest to each model lon, lat
- `index(2,:,:)` corresponds to the index (j) of the lattidude point in the data file, closest to each model lon, lat
- `index(3,:,:)` corresponds to the tile (it should be 1 since data_override does not support interpolation **from** cubesphere grids)
  - From there the four corners are (i,j), (i,j+1) (i+1) (i+1,j+1)
- The weights for the four corners
  - weight(:,:,1) -> (i,j)
  - weight(:,:,2) -> (i,j+1)
  - weight(:,:,3) -> (i+1,j)
  - weight(:,:,4) -> (i+1,j+1)
