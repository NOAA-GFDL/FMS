#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL Flexible Modeling System (FMS).
#
# FMS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FMS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# @uramirez, Ed Hartnett, @underwoo

# Find the needed mod and inc files.
AM_CPPFLAGS = -I${top_srcdir}/include -I$(MODDIR)

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_mpp \
  test_super_grid \
  test_mpp_domains \
  test_redistribute_int \
  test_mpp_memuse \
  test_mpp_mem_dump \
  test_mpp_memutils_begin_end \
  test_mpp_print_memuse_stats_stderr \
  test_mpp_print_memuse_stats_file \
  test_mpp_memutils_begin_2x \
  test_mpp_memutils_end_before_begin \
  test_read_ascii_file \
  test_read_input_nml \
  test_stdout \
  test_stderr \
  test_mpp_get_ascii_lines \
  test_system_clock \
  test_mpp_broadcast \
  test_clock_init \
  test_domains_simple \
  test_mpp_npes \
  test_mpp_pe \
  test_mpp_root_pe \
  test_peset \
  test_mpp_update_domains \
  test_mpp_gatscat \
  test_mpp_sum \
  test_update_domains_performance \
  test_minmax \
  test_mpp_sendrecv \
  test_global_arrays \
  test_chksum_int \
  test_mpp_update_domains_ad \
  test_mpp_transmit \
  test_mpp_alltoall \
  test_mpp_global_field \
  test_mpp_global_field_ug \
  test_mpp_global_sum_ad

# These are the sources for the tests.
test_mpp_SOURCES = test_mpp.F90
test_mpp_domains_SOURCES = \
  compare_data_checksums.F90 \
  test_domains_utility_mod.F90 \
  test_mpp_domains.F90
test_mpp_memuse_SOURCES = test_mpp_memuse.F90
test_mpp_mem_dump_SOURCES = test_mpp_mem_dump.F90
test_mpp_memutils_begin_end_SOURCES = test_mpp_memutils_begin_end.F90
test_mpp_print_memuse_stats_stderr_SOURCES = test_mpp_print_memuse_stats_stderr.F90
test_mpp_print_memuse_stats_file_SOURCES = test_mpp_print_memuse_stats_file.F90
test_mpp_memutils_begin_2x_SOURCES = test_mpp_memutils_begin_2x.F90
test_mpp_memutils_end_before_begin_SOURCES = test_mpp_memutils_end_before_begin.F90
test_read_ascii_file_SOURCES = test_read_ascii_file.F90
test_read_input_nml_SOURCES = test_read_input_nml.F90
test_stdout_SOURCES = test_stdout.F90
test_stderr_SOURCES = test_stderr.F90
test_mpp_get_ascii_lines_SOURCES = test_mpp_get_ascii_lines.F90
test_system_clock_SOURCES = test_system_clock.F90
test_mpp_broadcast_SOURCES = test_mpp_broadcast.F90
test_clock_init_SOURCES = test_clock_init.F90
test_domains_simple_SOURCES = test_domains_simple.F90
test_mpp_npes_SOURCES = test_mpp_npes.F90
test_mpp_pe_SOURCES = test_mpp_pe.F90
test_mpp_root_pe_SOURCES = test_mpp_root_pe.F90
test_peset_SOURCES = test_peset.F90
test_mpp_update_domains_SOURCES = \
  fill_halo.F90 \
  compare_data_checksums.F90 \
  compare_data_checksums_int.F90 \
  test_mpp_update_domains_real.F90 \
  test_mpp_update_domains_int.F90 \
  test_mpp_update_domains_main.F90
test_mpp_gatscat_SOURCES = test_mpp_gatscat.F90
test_mpp_sendrecv_SOURCES = test_mpp_sendrecv.F90
test_mpp_sum_SOURCES = test_mpp_sum.F90
test_update_domains_performance_SOURCES = \
  compare_data_checksums.F90 \
  compare_data_checksums_int.F90 \
  test_update_domains_performance.F90
test_minmax_SOURCES = test_minmax.F90
test_mpp_update_domains_ad_SOURCES = test_mpp_update_domains_ad.F90
test_global_arrays_SOURCES = test_global_arrays.F90
test_chksum_int_SOURCES = test_chksum_int.F90
test_redistribute_int_SOURCES = test_redistribute_int.F90
test_mpp_transmit_SOURCES = test_mpp_transmit.F90
test_mpp_alltoall_SOURCES = test_mpp_alltoall.F90
test_mpp_global_field_SOURCES = \
  compare_data_checksums.F90 \
  compare_data_checksums_int.F90 \
  test_mpp_global_field.F90
test_mpp_global_field_ug_SOURCES = \
  compare_data_checksums.F90 \
  compare_data_checksums_int.F90 \
  test_mpp_global_field_ug.F90
test_mpp_global_sum_ad_SOURCES = test_mpp_global_sum_ad.F90
test_super_grid_SOURCES = test_super_grid.F90

# Run the test programs.
TESTS = test_mpp_domains2.sh \
  test_redistribute_int.sh \
  test_global_arrays.sh  \
  test_mpp2.sh \
  test_mpp_memuse \
  test_mpp_mem_dump \
  test_mpp_memutils_mod.sh \
  test_read_ascii_file.sh \
  test_read_input_nml2.sh \
  test_stdout.sh \
  test_stderr.sh \
  test_mpp_get_ascii_lines2.sh \
  test_system_clock.sh \
  test_mpp_broadcast.sh \
  test_clock_init.sh \
  test_mpp_npes.sh \
  test_mpp_pe.sh \
  test_mpp_root_pe.sh \
  test_peset.sh \
  test_mpp_update_domains.sh \
  test_mpp_sum.sh \
  test_mpp_gatscat.sh \
  test_update_domains_performance.sh \
  test_minmax.sh \
  test_mpp_sendrecv.sh \
  test_chksum_int.sh  \
  test_mpp_update_domains_ad.sh \
  test_mpp_transmit.sh \
  test_mpp_alltoall.sh \
  test_mpp_global_field.sh \
  test_mpp_global_field_ug.sh \
  test_mpp_global_sum_ad.sh \
  test_super_grid.sh

# These files will also be included in the distribution.
EXTRA_DIST = input_base.nml \
  test_mpp_domains2.sh \
  test_mpp2.sh \
  test_mpp_memutils_mod.sh \
  test_read_ascii_file.sh \
  test_read_input_nml2.sh \
  test_stdout.sh \
  test_stderr.sh \
  test_mpp_get_ascii_lines2.sh \
  base_ascii_5 \
  base_ascii_25 \
  base_ascii_0 \
  base_ascii_skip \
  base_ascii_long \
  test_system_clock.sh \
  test_mpp_broadcast.sh \
  test_clock_init.sh \
  test_mpp_npes.sh \
  test_mpp_pe.sh \
  test_mpp_root_pe.sh \
  test_peset.sh \
  test_mpp_update_domains.sh \
  test_mpp_sum.sh \
  test_mpp_gatscat.sh \
  test_update_domains_performance.sh \
  test_minmax.sh \
  test_mpp_sendrecv.sh \
  test_global_arrays.sh \
  test_chksum_int.sh \
  test_redistribute_int.sh \
  test_mpp_update_domains_ad.sh \
  test_mpp_transmit.sh \
  test_mpp_alltoall.sh \
  test_mpp_global_field.sh \
  test_mpp_global_field_ug.sh \
  test_mpp_global_sum_ad.sh \
  test_super_grid.sh

fill_halo.mod: fill_halo.$(OBJEXT)
compare_data_checksums.mod: compare_data_checksums.$(OBJEXT)
compare_data_checksums_int.mod: compare_data_checksums_int.$(OBJEXT)
test_mpp_update_domains_real.mod: test_mpp_update_domains_real.$(OBJEXT)
test_mpp_update_domains_int.mod: test_mpp_update_domains_int.$(OBJEXT)
test_domains_utility_mod.mod: test_domains_utility_mod.$(OBJEXT)

test_mpp_update_domains_real.$(OBJEXT): compare_data_checksums.mod fill_halo.mod test_mpp_update_domains
test_mpp_update_domains_int.$(OBJEXT): compare_data_checksums_int.mod fill_halo.mod

test_mpp_update_domains_main.$(OBJEXT): test_mpp_update_domains_real.mod test_mpp_update_domains_int.mod
test_update_domains_performance.$(OBJEXT): compare_data_checksums_int.mod compare_data_checksums.mod
test_mpp_global_field.$(OBJEXT): compare_data_checksums_int.mod compare_data_checksums.mod
test_mpp_global_field_ug.$(OBJEXT): compare_data_checksums_int.mod compare_data_checksums.mod
test_mpp_domains.$(OBJEXT): compare_data_checksums.mod test_domains_utility_mod.mod

# Clean up
CLEANFILES = input.nml input_alternative.nml input_blank.nml empty.nml *.out* *.tst* \
  ascii* test_numb* *.mod
