# This is an automake file for the mpp directory of the MPP
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -Iinclude

noinst_LTLIBRARIES = libmpp.la

libmpp_la_SOURCES = mpp_parameter.F90 nsclock.c affinity.c		\
threadloc.c mpp_data.F90 mpp.F90 mpp_memutils.F90 mpp_pset.F90		\
mpp_efp.F90 mpp_utilities.F90 mpp_domains.F90 mpp_io.F90		\
test_mpp_domains.F90 test_mpp.F90 test_mpp_io.F90 test_mpp_pset.F90
