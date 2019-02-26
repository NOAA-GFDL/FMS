# This is an automake file for the mpp directory of the MPP
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -Iinclude

noinst_LTLIBRARIES = libmpp_parameter.la libmpp.la libmpp_data.la	\
libmpp_utilities.la libmpp_domains.la libmpp_io.la libmpp_c.la		\
libmpp_pset.la libmpp_efp.la

libmpp_parameter_la_SOURCES = mpp_parameter.F90
libmpp_la_SOURCES = mpp.F90
libmpp_data_la_SOURCES = mpp_data.F90
libmpp_utilities_la_SOURCES = mpp_utilities.F90
libmpp_domains_la_SOURCES = mpp_domains.F90
libmpp_io_la_SOURCES = mpp_io.F90
libmpp_pset_la_SOURCES = mpp_pset.F90
libmpp_efp_la_SOURCES = mpp_efp.F90
libmpp_c_la_SOURCES = nsclock.c affinity.c threadloc.c

SUBDIRS = include

# check_PROGRAMS = test_mpp
# TESTS = test_mpp
# test_mpp_SOURCES = test_mpp.F90

#test_mpp_domains.F90 test_mpp.F90 test_mpp_io.F90 test_mpp_pset.F90
