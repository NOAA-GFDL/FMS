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

# This is an automake file for the test_fms/diag_manager directory of the FMS
# package.
# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I$(top_srcdir)/include  -I$(top_srcdir)/diag_manager/include -I$(MODDIR)

# Link to the FMS library.
LDADD = $(top_builddir)/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_diag_manager test_diag_manager_time \
 test_diag_yaml test_diag_ocean test_modern_diag test_diag_buffer \
 test_flexible_time test_diag_update_buffer test_reduction_methods check_time_none \
 check_time_min check_time_max check_time_sum check_time_avg test_diag_diurnal check_time_diurnal \
 check_time_pow check_time_rms check_subregional test_cell_measures test_var_masks \
 check_var_masks

# This is the source code for the test.
test_diag_manager_SOURCES = test_diag_manager.F90
test_diag_manager_time_SOURCES = test_diag_manager_time.F90
test_diag_update_buffer_SOURCES= test_diag_update_buffer.F90
test_diag_yaml_SOURCES = test_diag_yaml.F90
test_diag_ocean_SOURCES = test_diag_ocean.F90
test_modern_diag_SOURCES = test_modern_diag.F90
test_diag_buffer_SOURCES= test_diag_buffer.F90
test_flexible_time_SOURCES = test_flexible_time.F90
test_reduction_methods_SOURCES = testing_utils.F90 test_reduction_methods.F90
test_diag_diurnal_SOURCES = testing_utils.F90 test_diag_diurnal.F90
check_time_none_SOURCES = testing_utils.F90 check_time_none.F90
check_time_min_SOURCES = testing_utils.F90 check_time_min.F90
check_time_max_SOURCES = testing_utils.F90 check_time_max.F90
check_time_sum_SOURCES = testing_utils.F90 check_time_sum.F90
check_time_avg_SOURCES = testing_utils.F90 check_time_avg.F90
check_time_diurnal_SOURCES = testing_utils.F90 check_time_diurnal.F90
check_time_pow_SOURCES = testing_utils.F90 check_time_pow.F90
check_time_rms_SOURCES = testing_utils.F90 check_time_rms.F90
test_cell_measures_SOURCES = test_cell_measures.F90
check_subregional_SOURCES = check_subregional.F90
test_var_masks_SOURCES = test_var_masks.F90
check_var_masks_SOURCES = check_var_masks.F90

TEST_EXTENSIONS = .sh
SH_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(abs_top_srcdir)/test_fms/tap-driver.sh

# Run the test.
TESTS = test_diag_manager2.sh test_time_none.sh test_time_min.sh test_time_max.sh test_time_sum.sh \
        test_time_avg.sh test_time_pow.sh test_time_rms.sh test_time_diurnal.sh test_cell_measures.sh \
        test_subregional.sh test_var_masks.sh

testing_utils.mod: testing_utils.$(OBJEXT)

# Copy over other needed files to the srcdir
EXTRA_DIST = test_diag_manager2.sh check_crashes.sh test_time_none.sh test_time_min.sh test_time_max.sh \
             test_time_sum.sh test_time_avg.sh test_time_pow.sh test_time_rms.sh test_time_diurnal.sh \
             test_cell_measures.sh test_subregional.sh test_var_masks.sh

if USING_YAML
skipflag=""
else
skipflag="skip"
endif

TESTS_ENVIRONMENT = skipflag=${skipflag}

CLEANFILES = *.yaml input.nml *.nc *.out diag_table* *-files/* *.dpi *.spi *.dyn *.spl *.mod

