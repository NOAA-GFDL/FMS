@ingroup data_override_mod

# Data_override

### Contents
- [YAML Data Table format](README.MD#1-yaml-data-table-format)
- [How to use it?](README.MD#2-how-to-use-it)
- [Converting legacy data_table to data_table.yaml](README.MD#3-converting-legacy-data_table-to-data_tableyaml)
- [Examples](README.MD#4-examples)

#### 1. YAML Data Table format:
Each entry in the data_table has the following key values:
- **gridname:** Name of the grid to interpolate the data to. The acceptable values are "ICE", "OCN", "ATM", and "LND"
- **fieldname_code:** Name of the field as it is in the code to interpolate.
- **fieldname_file:** Name of the field as it is writen in the file. **Required** only if overriding from a file
- **file_name:** Name of the file where the variable is located, including the directory. **Required** only if overriding from a file
- **interpol_method:** Method used to interpolate the field. The acceptable values are "bilinear", "bicubic", and "none". "none" implies that the field in the file is already in the model grid. The LIMA format is no longer supported. **Required** only if overriding from a file
- **factor:** Factor that will be multiplied after the data is interpolated

If it is desired to interpolate the data to a region of the model grid. The following **optional** arguments are available.
- **region_type:** The region type. The acceptable values are "inside_region" and "outside_region"
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
 - gridname          : ICE
   fieldname_code    : sic_obs
   fieldname_file    : sic
   file_name         : INPUT/hadisst_ice.data.nc
   interpol_method   : bilinear
   factor            : 0.01
```
Which corresponds to the following model code:
```F90
call data_override('ICE', 'sic_obs', icec, Spec_Time)
```
where:
- `ICE` corresponds to the gridname in the data_table
- `sic_obs` corresponds to the fieldname_code in the data_table
- `icec` is the variable to write the data to
- `Spec_Time` is the time to interpolate the data to.

Additionally, it is required to call data_override_init (in this case with the ICE domain). The grid_spec.nc file must also contain the coordinate information for the domain being used. 

```F90
call data_override_init(Ice_domain_in=Ice_domain)
```

**4.2** The following example is going to simply multiply `sic_obs` by a factor of 0.01.

In the **legacy format**, the data_table will look like:
```
"ICE", "sit_obs", "",    "INPUT/hadisst_ice.data.nc", "none",     2.0
```

In the **yaml format**, the data_table will look like:
```
data_table:
 - gridname          : ICE
   fieldname_code    : sit_obs
   factor            : 0.01
```

Which corresponds to the following model code:
```F90
call data_override('ICE', 'sit_obs', icec, Spec_Time)
```
where:
- `ICE` corresponds to the gridname in the data_table
- `sit_obs` corresponds to the fieldname_code in the data_table
- `icec` is the variable to write the data to
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
```
data_table:
 - gridname          : OCN
   fieldname_code    : runoff
   fieldname_file    : runoff
   file_name         : INPUT/runoff.daitren.clim.nc
   interpol_method   : none
   factor            : 1.0
```

Which corresponds to the following model code:
```F90
call data_override('OCN', 'runoff', runoff_data, Spec_Time)
```
where:
- `OCN` corresponds to the gridname in the data_table
- `runoff` corresponds to the fieldname_code in the data_table
- `runoff_data` is the variable to write the data to
- `Spec_Time` is the time to interpolate the data to.

Additionally, it is required to call data_override_init (in this case with the ocean domain). The grid_spec.nc file is still required to initialize data_override with the ocean domain and to determine if the data in the file is in the same grid as the ocean.

```F90
call data_override_init(Ocn_domain_in=Ocn_domain)
```
