# This is an automake file for the fms directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/constants
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/platform
AM_CPPFLAGS += -I${top_builddir}/memutils

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libfms_io.la libfms.la

# Each convenience library depends on its source.
libfms_la_SOURCES = fms.F90
libfms_io_la_SOURCES = fms_io.F90

# These inc files are in the fms directory.
FMS_INC_FILES = fms_io_unstructured_field_exist.inc \
fms_io_unstructured_get_file_name.inc \
fms_io_unstructured_register_restart_axis.inc \
fms_io_unstructured_setup_one_field.inc read_data_4d.inc \
fms_io_unstructured_file_unit.inc \
fms_io_unstructured_get_file_unit.inc \
fms_io_unstructured_register_restart_field.inc read_data_2d.inc \
write_data.inc fms_io_unstructured_get_field_size.inc \
fms_io_unstructured_read.inc fms_io_unstructured_save_restart.inc \
read_data_3d.inc

# Each mod file depends on the .o file.
fms_io_mod.mod: fms_io.$(OBJEXT)
fms_mod.mod: fms.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
fms.$(OBJEXT): fms_io_mod.mod

# Mod files are built and then installed as headers.
MODFILES = fms_io_mod.mod fms_mod.mod
BUILT_SOURCES = $(MODFILES)
include_HEADERS = $(MODFILES) $(FMS_INC_FILES)

CLEANFILES = *.mod
