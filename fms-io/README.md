# Fms I/O 2.0

# Overview
This fortran library aims to support all "distributed" I/O required by GFDL's
climate models.

# Requirements
This library requires a c and fortran MPI compiler.  By default, the compiler
flags in the provided makefiles are set up to use the Intel compiler suite.
These can be overridden by setting FC, CC, FFLAGS, and CFLAGS respectively.

# Source Code
The source code currently resides in
[this Gitlab repository](https://gitlab.gfdl.noaa.gov/Raymond.Menzel/fms-io).
To obtain the code, run

```
$ git clone https://gitlab.gfdl.noaa.gov/Raymond.Menzel/fms-io.git
$ cd fms-io
```

# Setting Up Your Environment
Scripts for both bash and tcsh are provided.  These scripts load the
necessary modules and alias make to point to the correct makefiles.

### Gaea
For bash:

```
$ source env/bash/gaea.bashrc
```

Using tcsh:

```
$ source env/tcsh/gaea.cshrc
```

### GFDL Linux Workstation
Using bash:

```
$ source env/bash/gfdl-ws.bashrc
```

Using tcsh:

```
$ source env/tcsh/gfdl-ws.cshrc
```

# Building
To build using default settings (assuming you have MPI and the intel compiler
suite installed), run

```
$ cd mpp
$ make install
$ cd -
$ make
```

# APIs
Documentation for the API was done with Doxygen. To view Doxygen-generated
HTML describing each API, please open the file API.html in a browser.

# Example Codes
A short, complete example program is included in the examples directory.
These examples demonstrate how to properly use each library function, and
can be built and run by simply running:

```
$ cd examples
$ make
$ make test
```
