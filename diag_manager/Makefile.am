# This is an automake file for the diag_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

AM_CPPFLAGS = -I${top_srcdir}/include
AM_CPPFLAGS += -I${top_builddir}/mpp
AM_CPPFLAGS += -I${top_builddir}/fms
AM_CPPFLAGS += -I${top_builddir}/time_manager
AM_CPPFLAGS += -I${top_builddir}/constants

noinst_LTLIBRARIES = libdiag_data.la libdiag_axis.la libdiag_grid.la	\
libdiag_output.la libdiag_util.la libdiag_manifest.la 		\
libdiag_table.la libdiag_manager.la

libdiag_axis_la_SOURCES = diag_axis.F90
libdiag_data_la_SOURCES = diag_data.F90
libdiag_grid_la_SOURCES = diag_grid.F90
libdiag_manager_la_SOURCES = diag_manager.F90
libdiag_manifest_la_SOURCES = diag_manifest.F90
libdiag_output_la_SOURCES = diag_output.F90
libdiag_table_la_SOURCES = diag_table.F90
libdiag_util_la_SOURCES = diag_util.F90

CLEANFILES = *.mod
