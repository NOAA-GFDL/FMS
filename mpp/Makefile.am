# This is an automake file for the mpp directory of the MPP
# package.

# Ed Hartnett 2/22/19

# Descend into include directory.
SUBDIRS = include

# Include .h and .mod files.
AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_srcdir}/mpp/include

# Build these uninstalled convenience libraries.
noinst_LTLIBRARIES = libmpp_parameter.la libmpp_data.la libmpp.la \
libmpp_utilities.la libmpp_memutils.la libmpp_pset.la libmpp_efp.la \
libmpp_domains.la libmpp_io.la libmpp_c.la

# Each convenience library depends on its source.
libmpp_parameter_la_SOURCES = mpp_parameter.F90
libmpp_la_SOURCES = mpp.F90
libmpp_data_la_SOURCES = mpp_data.F90
libmpp_utilities_la_SOURCES = mpp_utilities.F90
libmpp_domains_la_SOURCES = mpp_domains.F90
libmpp_io_la_SOURCES = mpp_io.F90
libmpp_pset_la_SOURCES = mpp_pset.F90
libmpp_efp_la_SOURCES = mpp_efp.F90
libmpp_c_la_SOURCES = nsclock.c threadloc.c
libmpp_memutils_la_SOURCES = mpp_memutils.F90

# Each mod file depends on the .o file.
mpp_parameter_mod.mod: mpp_parameter.$(OBJEXT)
mpp_data_mod.mod: mpp_data.$(OBJEXT)
mpp_mod.mod: mpp.$(OBJEXT)
mpp_pset_mod.mod: mpp_pset.$(OBJEXT)
mpp_utilities_mod.mod: mpp_utilities.$(OBJEXT)
mpp_memutils_mod.mod: mpp_memutils.$(OBJEXT)
mpp_efp_mod.mod: mpp_efp.$(OBJEXT)
mpp_domains_mod.mod: mpp_domains.$(OBJEXT)
mpp_io_mod.mod: mpp_io.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
mpp_data.$(OBJEXT): mpp_parameter_mod.mod
mpp.$(OBJEXT): mpp_data_mod.mod
mpp_utilities.$(OBJEXT): mpp_mod.mod mpp_efp_mod.mod
mpp_memutils.$(OBJEXT): mpp_mod.mod
mpp_pset.$(OBJEXT): mpp_mod.mod
mpp_efp.$(OBJEXT): mpp_parameter_mod.mod mpp_mod.mod
mpp_domains.$(OBJEXT): mpp_data_mod.mod mpp_parameter_mod.mod \
  mpp_mod.mod mpp_memutils_mod.mod mpp_pset_mod.mod mpp_efp_mod.mod
mpp_io.$(OBJEXT): mpp_parameter_mod.mod mpp_mod.mod mpp_domains_mod.mod

# Mod files are built and then installed as headers.
MODFILES = mpp_parameter_mod.mod mpp_data_mod.mod mpp_mod.mod \
mpp_utilities_mod.mod mpp_memutils_mod.mod mpp_pset_mod.mod \
mpp_efp_mod.mod mpp_domains_mod.mod mpp_io_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
