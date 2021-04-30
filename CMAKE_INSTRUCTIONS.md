# Instructions for building FMS with cmake

## 1. Environment Variables

### For GNU compilers on Linux with the bash shell:

```
export FC=mpifort
export CC=mpicc
```

### For Intel Compilers on Linux with the bash shell:

```
export FC=mpiifort
export CC=mpiicc
```

### NetCDF is provided via the `nc-config` command:
```
export NetCDF_ROOT=`nc-config --prefix`
```


## 2. Build and install FMS with CMake
`<prefix>` is the full install directory for FMS provided by user

```
cd FMS
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<prefix> ..
make -j4
make install
```

### User configurable options:
By default, FMS is built without `OpenMP` and in `single precision (r4)`

The following build options are available:
```
-DOPENMP   "Build FMS with OpenMP support" DEFAULT: OFF
-D32BIT    "Build 32-bit (r4) FMS library" DEFAULT: ON
-D64BIT    "Build 64-bit (r8) FMS library" DEFAULT: OFF

-DINTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"     DEFAULT: ON
-DENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION" DEFAULT: ON
-DGFS_PHYS              "Enable compiler definition -DGFS_PHYS"              DEFAULT:OFF
-DLARGEFILE             "Enable compiler definition -Duse_LARGEFILE"         DEFAULT:OFF
```

## 3. Installation structure

When the above command finishes, the `<prefix>` will have an `include_r4` and a `lib` directory. The `lib ` directory will have these files:

```
libfms_r4.a
cmake/fms/fms-targets.cmake
cmake/fms/fms-targets-release.cmake
cmake/fms/fms-config.cmake
cmake/fms/fms-config-version.cmake
```

## 4. Using FMS in your application

FMS built with `cmake` provides transient targets with its package configuration.
To be able to look for FMS in your application, set the following environment variable:
```
export FMS_ROOT=<prefix>
```
where `<prefix>` is the full path to FMS installation.

To find FMS in your application `CMakeLists.txt`:

```
find_package(FMS)
```

To link your application with FMS library:
```
target_link_libraries(appName FMS::fms_r4)
```

If your application does not provide a means to locate the NetCDF installation via cmake, this may help:
```
-DCMAKE_MODULE_PATH=<FMS_SRC_DIR>/cmake
```
