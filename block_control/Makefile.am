# This is an automake file for the block_control directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp

noinst_LTLIBRARIES = libblock_control.la

libblock_control_la_SOURCES = block_control.F90

block_control_mod.mod: block_control.$(OBJEXT)
BUILT_SOURCES = block_control_mod.mod
include_HEADERS = block_control_mod.mod

CLEANFILES = *.mod
