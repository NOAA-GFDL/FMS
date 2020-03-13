# This is an automake file for the fms2_io directory of the FMS package.

# Descend into include directory.
SUBDIRS = include

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_srcdir}/fms2_io/include

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfms2_io.la libfms_io_utils.la libnetcdf_io.la \
                     libfms_netcdf_domain_io.la libfms_netcdf_unstructured_domain_io.la \
                     libblackboxio.la

# Each convenience library depends on its source.
libfms2_io_la_SOURCES = fms2_io.F90
libfms_io_utils_la_SOURCES = fms_io_utils.F90
libnetcdf_io_la_SOURCES = netcdf_io.F90
libfms_netcdf_domain_io_la_SOURCES = fms_netcdf_domain_io.F90
libfms_netcdf_unstructured_domain_io_la_SOURCES = fms_netcdf_unstructured_domain_io.F90
libblackboxio_la_SOURCES = blackboxio.F90

# Each mod file depends on the .o file.
fms2_io_mod.mod: fms2_io.$(OBJEXT)
fms_io_utils_mod.mod: fms_io_utils.$(OBJEXT)
netcdf_io_mod.mod: netcdf_io.$(OBJEXT)
fms_netcdf_domain_io_mod.mod: fms_netcdf_domain_io.$(OBJEXT)
fms_netcdf_unstructured_domain_io_mod.mod: fms_netcdf_unstructured_domain_io.$(OBJEXT)
blackboxio.mod: blackboxio.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
fms2_io.$(OBJEXT): fms_io_utils_mod.mod netcdf_io_mod.mod fms_netcdf_domain_io_mod.mod \
                   fms_netcdf_unstructured_domain_io_mod.mod blackboxio.mod
netcdf_io.$(OBJEXT): fms_io_utils_mod.mod
fms_netcdf_domain_io.$(OBJEXT): fms_io_utils_mod.mod netcdf_io_mod.mod
fms_netcdf_unstructured_domain_io.$(OBJEXT): fms_io_utils_mod.mod netcdf_io_mod.mod
blackboxio.$(OBJEXT): fms_io_utils_mod.mod netcdf_io_mod.mod fms_netcdf_domain_io_mod.mod \
                      fms_netcdf_unstructured_domain_io_mod.mod

# Mod files are built and then installed as headers.
MODFILES = fms2_io_mod.mod fms_io_utils_mod.mod netcdf_io_mod.mod \
           fms_netcdf_domain_io_mod.mod fms_netcdf_unstructured_domain_io_mod.mod \
           blackboxio.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES)

CLEANFILES = *.mod
