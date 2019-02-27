# This is an automake file for the mpp directory of the MPP
# package.

# Ed Hartnett 2/22/19

# Turn off parallel builds in this directory.
.NOTPARALLEL:

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_srcdir}/mpp/include
#AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libmpp_parameter.la libmpp_data.la libmpp.la	\
libmpp_utilities.la libmpp_memutils.la libmpp_pset.la libmpp_efp.la	\
libmpp_domains.la libmpp_io.la libmpp_c.la

libmpp_parameter_la_SOURCES = mpp_parameter.F90
libmpp_la_SOURCES = mpp.F90
libmpp_data_la_SOURCES = mpp_data.F90
libmpp_utilities_la_SOURCES = mpp_utilities.F90
libmpp_domains_la_SOURCES = mpp_domains.F90
libmpp_io_la_SOURCES = mpp_io.F90
libmpp_pset_la_SOURCES = mpp_pset.F90
libmpp_efp_la_SOURCES = mpp_efp.F90
libmpp_c_la_SOURCES = nsclock.c affinity.c threadloc.c
libmpp_memutils_la_SOURCES = mpp_memutils.F90

SUBDIRS = include

# check_PROGRAMS = test_mpp
# TESTS = test_mpp
# test_mpp_SOURCES = test_mpp.F90

#test_mpp_domains.F90 test_mpp.F90 test_mpp_io.F90 test_mpp_pset.F90

mpp_data_mod.mod: mpp_data.$(OBJEXT)
mpp_domains_mod.mod: mpp_domains.$(OBJEXT)
mpp_efp_mod.mod: mpp_efp.$(OBJEXT)
mpp_io_mod.mod: mpp_io.$(OBJEXT)
mpp_memutils_mod.mod: mpp_memutils.$(OBJEXT)
mpp_mod.mod: mpp.$(OBJEXT)
mpp_parameter_mod.mod: mpp_parameter.$(OBJEXT)
mpp_pset_mod.mod: mpp_pset.$(OBJEXT)
mpp_utilities_mod.mod: mpp_utilities.$(OBJEXT)

MODFILES = mpp_parameter_mod.mod mpp_data_mod.mod mpp_mod.mod	\
mpp_utilities_mod.mod mpp_memutils_mod.mod mpp_pset_mod.mod	\
mpp_efp_mod.mod mpp_domains_mod.mod mpp_io_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

CLEANFILES = *.mod
