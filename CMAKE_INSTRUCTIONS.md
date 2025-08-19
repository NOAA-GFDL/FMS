\page cmake Building with CMake
# Instructions for building FMS with cmake

## 1. Configuring the build

### Environment Variables

#### For GNU compilers on Linux with the bash shell:
```
export FC=mpifort
export CC=mpicc
```
#### For Intel Compilers on Linux with the bash shell:
```
export FC=mpiifort
export CC=mpiicc
```
MPI compiler wrappers may be named different on your system, so its best to ensure the compiler commands work prior to being set.

#### NetCDF is provided via the `nc-config` command:
```
export NetCDF_ROOT=`nc-config --prefix`
```

#### If building with yaml parser (-DWITH_YAML)
```
export LIBYAML_ROOT=<your libyaml install directory>
```

### Running CMake
It's best to create a build directory inside of the FMS folder to avoid building on top of the source code.
Once that is done, CMake can be ran to generate the necessary build files:
```
mkdir build && cd build
cmake [-DWITH_YAML=on] ..
```

### Setting a Build Type
CMake uses "build types" to set compiler flags. By default (.ie no -DCMAKE_BUILD_TYPE="Type" argument when cmake is run), it will use the "Release" build type which corresponds the the prod/opt flags used by the GFDL mkmf templates.

Similarly "Repro" and "Debug" also correspond to whatever options would be set by theircorresponding mkmf template files.

To reproduce the old behavior for the cmake build (2025.03 and prior), the "ReleaseUFS" and "DebugUFS" will use the same flags as before.

Please open an issue in the FMS repository if you have any issues or would like to request any changes to the current build type options.

### Setting custom flags
To override the default compiler flags, you can set them via CMake arguments:
```
cmake -DCMAKE_BUILD_TYPE="NoFlags" -DCMAKE_C_FLAGS="<c flags>" -DCMAKE_Fortran_FLAGS="<fortran flags>"
```
By using the "NoFlags" build type, only necessary compilation flags will be added besides what is set via the cmake arguments.

## 2. Build and install FMS with CMake
`<prefix>` is the full install directory for FMS provided by user

```
cd FMS
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<prefix> ..
make -j
make install
```

### User configurable options:
By default, FMS is built with `OpenMP` enabled and delivered in static library files.

FMS has mixed precision real support for most interfaces. By default, one library will be created,
`libfms`, that is compiled with r8 defaults but also contains overloaded r4 routines.

The 64BIT and 32BIT precision options will build distinct libraries when enabled with the given default
real size, libfms_r4 or libfms_r8. These option are provided for backwards compatibility, but are no longer
supported by our development team, since mixed precision can now be used with a single library.
Unit tests currently only work when no precision option is specified.

The following build options are available:
```
-DOPENMP      "Build FMS with OpenMP support"        DEFAULT: ON
-D32BIT       "Build 32-bit (r4) FMS library"        DEFAULT: OFF
-D64BIT       "Build 64-bit (r8) FMS library"        DEFAULT: OFF
-DFPIC        "Build with position independent code" DEFAULT: OFF
-DSHARED_LIBS "Build shared/dynamic libraries"       DEFAULT: OFF

-DCONSTANTS             "Build with <X> constants parameter definitions"     DEFAULT:GFDL  OPTIONS:GFS|GEOS|GFDL
-DINTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"     DEFAULT: ON
-DENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION" DEFAULT: ON
-DPORTABLE_KINDS        "Enable compiler definition -DPORTABLE_KINDS"        DEFAULT:OFF
-DGFS_PHYS              "Enable compiler definition -DGFS_PHYS"              DEFAULT:OFF
-DLARGEFILE             "Enable compiler definition -Duse_LARGEFILE"         DEFAULT:OFF
-DWITH_YAML             "Enable compiler definition -Duse_yaml"              DEFAULT:OFF
```

## 3. Installation structure

When the above command finishes, the `<prefix>` will have an `include` and a `lib` directory. The `lib ` directory will have these files:

```
libfms.a
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
