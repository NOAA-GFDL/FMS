\page autotools Building with Autotools
# Autotools Build System Documentation

This document describes the autotools-based build system for FMS.

## Introduction to Autotools

Autoconf, automake, and libtool and the GNU/Linux standard build
tools. When building a package based on autotools, the user does
something like:

./configure && make check install

The configure step queries the system about many things, and contructs
makefiles.

The make step uses the generated Makefiles to build, test, and install
the software.

Standard environment variables and configure options can be used to
control many aspects of the build. Custom configure options can easily
be added to support additional needs.

### Caution Concerning Generated Files

Autotools creates many generated files, which should not be edited or
checked into the repository. Simply ignore them. Do not try to edit
generated Makefiles, do not move or rename any of the shell scripts
that autoreconf puts in place to let autotools work. These files have
been added to .gitignore and should never be added to the repo.

# How to Build FMS

Previously, everyone built FMS by checking out code from git. However,
with the new build system, only those who want to contribute to the
code base need check out the code from git.

## As an FMS Developer

All FMS developers will need a reasobably reacent version of tools
autoconf, automake, and libtool. These are available on package
management systems. (Ex. yum install automake autoconf libtool).

The process of building FMS from the repo is:

1. Clone repo and cd into repo directory.
2. Run autoreconf -i to build the developer build system.
3. Run ./configure to configure.
4. Run make to build.

## As an FMS User

Users start with a tarball, not the git repo. They do not have to have
any of the autotools installed. Thier build process is:

1. Unpack the tarball and cd into the directory.
2. Run ./configure --prefix=/my/installdir
3. Run make install

## Precious Flags for configure

Some environment variables are important to the autotools build
system, these are known as "precious" variables. One example is CC,
which should be set to the C compiler.

It's common to set some precious vars before the build. Commonly used
ones include:
* CC the C compiler
* FC the Fortran compiler
* CPPFLAGS C (and Fortran) pre-processor flags
* FCFLAGS Fortran compiler flags
* LDFLAGS Linker flags

## Standard Configure Options

The configure script has some standard options, including:
* --prefix allows user to specify install directory
* --disable-shared disables building of shared library
* --help prints message showing all options

## FMS Configure Options

Currently there are not FMS specific configure options, but probably
we will add some.

## Standard Make Targets

Some of the useful make targets include:
* make or make all - build code
* make install - build code (as needed) and install
* make check - build code (as needed) and run tests
* make clean - clean back build
* make distclean - clean configure output and build
* make dist - create a tarball for distribution
* make distcheck - create a tarball, unpack it, build and run tests, then then clean it.

