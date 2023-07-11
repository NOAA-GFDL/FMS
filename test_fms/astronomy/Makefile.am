#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the test_fms/astronomy directory of the
# FMS package.

# Caitlyn McAllister

# Find the fms and mpp mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = \
  test_set_orbital_parameters_r4 \
  test_set_orbital_parameters_r8 \
  test_get_orbital_parameters_r4 \
  test_get_orbital_parameters_r8 \
  test_diurnal_solar_r4 \
  test_diurnal_solar_r8 \
  test_dt_option_diurnal_r4 \
  test_dt_option_diurnal_r8 \
  test_diurnal_solar_cal_r4 \
  test_diurnal_solar_cal_r8 \
  test_daily_mean_r4 \
  test_daily_mean_r8 \
  test_daily_mean_cal_r4 \
  test_daily_mean_cal_r8 \
  test_annual_mean_solar_2d_r4 \
  test_annual_mean_solar_2d_r8 \
  test_annual_mean_solar_1d_r4 \
  test_annual_mean_solar_1d_r8 \
  test_annual_mean_solar_2level_r4 \
  test_annual_mean_solar_2level_r8

# This is the source code for the test.
test_set_orbital_parameters_r4_SOURCES = test_set_orbital_parameters.F90
test_set_orbital_parameters_r8_SOURCES = test_set_orbital_parameters.F90

test_get_orbital_parameters_r4_SOURCES = test_get_orbital_parameters.F90
test_get_orbital_parameters_r8_SOURCES = test_get_orbital_parameters.F90

test_diurnal_solar_r4_SOURCES = test_diurnal_solar.F90
test_diurnal_solar_r8_SOURCES = test_diurnal_solar.F90

test_dt_option_diurnal_r4_SOURCES = test_dt_option_diurnal.F90
test_dt_option_diurnal_r8_SOURCES = test_dt_option_diurnal.F90

test_diurnal_solar_cal_r4_SOURCES = test_diurnal_solar_cal.F90
test_diurnal_solar_cal_r8_SOURCES = test_diurnal_solar_cal.F90

test_daily_mean_r4_SOURCES = test_daily_mean.F90
test_daily_mean_r8_SOURCES = test_daily_mean.F90

test_daily_mean_cal_r4_SOURCES = test_daily_mean_cal.F90
test_daily_mean_cal_r8_SOURCES = test_daily_mean_cal.F90

test_annual_mean_solar_2d_r4_SOURCES = test_annual_mean_solar_2d.F90
test_annual_mean_solar_2d_r8_SOURCES = test_annual_mean_solar_2d.F90

test_annual_mean_solar_1d_r4_SOURCES = test_annual_mean_solar_1d.F90
test_annual_mean_solar_1d_r8_SOURCES = test_annual_mean_solar_1d.F90

test_annual_mean_solar_2level_r4_SOURCES = test_annual_mean_solar_2level.F90
test_annual_mean_solar_2level_r8_SOURCES = test_annual_mean_solar_2level.F90

# Set r4_kind and r8_kind
test_set_orbital_parameters_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_set_orbital_parameters_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_get_orbital_parameters_r4_CPPFLAGS = $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_get_orbital_parameters_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_diurnal_solar_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_diurnal_solar_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_dt_option_diurnal_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_dt_option_diurnal_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_diurnal_solar_cal_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_diurnal_solar_cal_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_daily_mean_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_daily_mean_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_daily_mean_cal_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_daily_mean_cal_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_annual_mean_solar_2d_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_annual_mean_solar_2d_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_annual_mean_solar_1d_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_annual_mean_solar_1d_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

test_annual_mean_solar_2level_r4_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r4_kind
test_annual_mean_solar_2level_r8_CPPFLAGS= $(AM_CPPFLAGS) -DTEST_AST_KIND_=r8_kind

# Run the test program.
TESTS = test_astronomy.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = test_astronomy.sh

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Clean up
CLEANFILES = input.nml *.out* *.dpi *.spi *.dyn *.spl
