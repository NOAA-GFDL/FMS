# Unit Testing

FMS includes a suite of MPI unit tests for most subdirectories and modules that make up the library.

It consists of programs in the test_fms/ directory, with shell scripts to handle directory set up and input files.
test_lib.sh.in and tap-driver.sh provide additional helper functions used in the scripts and manage output.

The unit tests can be run either through the autotools build system, or through cmake via ctest.
The different build systems will run the exact same scripts and executables, but there are some differences in the commands used to start the testing.

## Autotools Unit Testing

### Running the Tests

1. Configure and build the code with autotools
```
mkdir build # create a build directory in FMS
cd build
autoreconf -if ../configure.ac
../configure <configure options>
make
```

2. Build and run suite
```
make check
```
This will compile any code not already compiled and then proceed to run the test scripts.

If a test fails, `make check` will stop running the rest of the unit tests once it is done with the current subdirectory.

To instead run all of the tests regardless of failures, the `-k` option can be given:
```
make check -k
```

### Debugging Output and Test Options for Autotools

Setting the environment variable TEST_VERBOSE will direct output to stdout as the test runs, while setting VERBOSE will only output on failure.
Logs are created for each test as well, with the name \<test script name\>.log in its corresponding test_fms/ directory.

To run an individual test:
```
make check -C test_fms/<test directory> TESTS=<test script name>
```

SKIP_TESTS can be set to in order to skip specific tests in a script. It uses the script name and test number, and takes ranges as well:
```
SKIP_TESTS="test_name.4 test_name.[1-3]"
```

Some options that effect the test suite can be set by passing options to the ./configure script that creates the makefiles
for the build system:

-    `--enable-code-coverage` allows for compilation with flags for coverage information.
     If enabled a coverage report can be generated with `make check-code-coverage`
-    `--enable-test-input=/path/to/input` turns on test scripts that require input netcdf files (interpolator, xgrid, data_override).
     This option is mainly used internally and in automated testing since we do not host the input data publicly.
-    `--with-yaml` compile with yaml input and enable its associated tests

## CMake Unit Testing

### Running the tests

1. Configure and build the code and tests with cmake
```
mkdir build
cd build
cmake <Any cmake options such as -DWITH_YAML=on> ..
make
```
2. Run the tests with ctest
```
ctest
```
Alternatively, the generated makefile can be used as well:
```
make test
```
### Debugging and Test Options for CMake

To rerun failed tests with debug output from a previous failed ctest command:
```
ctest --rerun-failed --verbose
```
All tests are labeled with the subdirectory they are in in the `test_fms` directory.
To run tests for a specific area of the code, such as the mpp modules:
```
ctest -L mpp
```
To run specific tests matching the given regex:
```
ctest -R <regex>
```

For more testing options with ctest, you can look to the [ctest documentation](https://cmake.org/cmake/help/latest/manual/ctest.1.html).
