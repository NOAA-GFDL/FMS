#!/bin/sh

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

# Copyright (c) 2019-2020 Ed Hartnett, Seth Underwood

# Set common test settings.
. ../test-lib.sh

setup_test () {
  local tnum=$(printf "%2.2d" $my_test_count)
  # Clean up any remaining files from previous tests
  rm -f input.nml diag_test_${tnum}*.nc 2>&1 > /dev/null
  cat <<_EOF > input.nml
&test_diag_manager_nml
   layout = 1,1
   test_number = $my_test_count
   nlon = 10
   nlat = 10
   nlev = 10
   io_layout = 1,1
   numthreads = 1
   dt_step = 1
   months = 0
   days = 1
/

&diag_manager_nml
   max_field_attributes=3
   debug_diag_manager=.true.
   use_mpp_io = .false.
/

&ensemble_nml
   ensemble_size = 1
/
_EOF
}

# create and enter directory for in/output files
output_dir

my_test_count=1
cat <<_EOF > diag_table
test_diag_manager_01
1 3 1 0 0 0

#output files
 "diag_test_01",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1", "dat1", "diag_test_01",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too large in x and y direction (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=2
cat <<_EOF > diag_table
test_diag_manager_02
1 3 1 0 0 0

#output files
 "diag_test_02",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1", "dat1", "diag_test_02",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too large in x direction (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=3
cat <<_EOF > diag_table
test_diag_manager_03
1 3 1 0 0 0

#output files
 "diag_test_03",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1", "dat1", "diag_test_03",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too large in y direction (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=4
cat <<_EOF > diag_table
test_diag_manager_04
1 3 1 0 0 0

#output files
 "diag_test_04_file1",  1, "days", 1, "days", "time"
 "diag_test_04_file2",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_04_file1",  "all", .false., "none", 2
 "test_mod",              "dat2", "dat2", "diag_test_04_file2",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too small in x and y direction, checks for 2 time steps (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=5
cat <<_EOF > diag_table
test_diag_manager_05
1 3 1 0 0 0

#output files
 "diag_test_05_file1",  1, "days", 1, "days", "time"
 "diag_test_05_file2",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_05_file1",  "all", .false., "none", 2
 "test_mod",              "dat2", "dat2", "diag_test_05_file2",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too small in x directions, checks for 2 time steps (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=6
cat <<_EOF > diag_table
test_diag_manager_06
1 3 1 0 0 0

#output files
 "diag_test_06_file1",  1, "days", 1, "days", "time"
 "diag_test_06_file2",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_06_file1",  "all", .false., "none", 2
 "test_mod",              "dat2", "dat2", "diag_test_06_file2",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too small in y direction, checks for 2 time steps (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=7
cat <<_EOF > diag_table
test_diag_manager_07
1 3 1 0 0 0

#output files
 "diag_test_07",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1", "dat1", "diag_test_07",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too large in x and y, with halos, 2 time steps (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=8
cat <<_EOF > diag_table
test_diag_manager_08
1 3 1 0 0 0

#output files
 "diag_test_08_file1",  1, "days", 1, "days", "time"
 "diag_test_08_file2",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_08_file1",  "all", .false., "none", 2
 "test_mod",              "dat2", "dat2", "diag_test_08_file2",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too small in x and y, with halos, 2 time steps (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=9
cat <<_EOF > diag_table
test_diag_manager_09
1 3 1 0 0 0

#output files
 "diag_test_09",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "bk",   "bk",   "diag_test_09",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too small, 1D, static global data (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=10
cat <<_EOF > diag_table
test_diag_manager_10
1 3 1 0 0 0

#output files
 "diag_test_10",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "bk",   "bk",   "diag_test_10",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Data array is too large, 1D, static global data (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=11
cat <<_EOF > diag_table
test_diag_manager_11
1 3 1 0 0 0

#output files
 "diag_test_11",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_11",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Missing je_in as an input (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=12
cat <<_EOF > diag_table
test_diag_manager_12
1 3 1 0 0 0

#output files
 "diag_test_12",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_12",  "all", .false., "none", 2

# Test of the error check that duplicate field names do not appear in same file
 "test_mod",              "dat2", "dat2", "diag_test_12",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Catch duplicate field in diag_table (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=13
cat <<_EOF > diag_table
test_diag_manager_13
1 3 1 0 0 0

#output files
 "diag_test_13_file1",  1, "days",   1, "days", "time"
 "diag_test_13_file2",  1, "months", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_13_file1",  "all", .false., "none", 2

# Test of WARNING message that no data is written when run length is less than output interval
 "test_mod",              "dat2", "dat2", "diag_test_13_file2",  "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Output interval greater than runlength (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=14
cat <<_EOF > diag_table
test_diag_manager_14
1990 1 29 0 0 0

#output files
 "diag_test_14", 1, "months", 1, "days", "time"

#output variables
# Test of check for invalid date. (Jan 29 1990 + one month = Feb 29 1990)
 "test_mod",              "dat2", "dat2", "diag_test_14", "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Catch invalid date in register_diag_field call (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=15
cat <<_EOF > diag_table
test_diag_manager_15
1 3 1 0 0 0

#output files
 "diag_test_15",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1",           "dat1",           "diag_test_15", "all", .false., "none", 2
 "test_diag_manager_mod", "solar_constant", "solar_constant", "diag_test_15", "all", .false., "none", 2
_EOF
setup_test
test_expect_success "OpenMP thread test (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=16
cat <<_EOF > diag_table
test_diag_manager_16
1 3 1 0 0 0

#output files
 "diag_test_16", 1, "months", 1, "days", "time"

#output variables
# Test for output file name to be modified with appended string
 "test_diag_manager_mod", "dat2", "dat2", "diag_test_16", "all", .false., "none", 2
_EOF
setup_test
test_expect_success "Filename appendix added (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=17
cat <<_EOF > diag_table
test_diag_manager_17
1 3 1 0 0 0

#output files
 "diag_test_17", 1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat2", "dat2_rms", "diag_test_17", "all", "rms",  "none", 2
 "test_diag_manager_mod", "dat2", "dat2",     "diag_test_17", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Root-mean-square (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=18
cat <<_EOF > diag_table
test_diag_manager_18
1 3 1 0 0 0

#output files
 "diag_test_18_file1", 1, "days", 1, "days", "time"
 "diag_test_18_file2", 1, "days", 1, "days", "time"

#output variables
#"test_mod",              "dat2h_2", "dat2h_2", "diag_test_18_file1", "all", .true., "none", 2
"test_diag_manager_mod", "dat1",    "dat1",    "diag_test_18_file1", "all", .true., "none", 2
"test_diag_manager_mod", "dat2",    "dat2",    "diag_test_18_file1", "all", .true., "none", 2
"test_mod",              "dat2h",   "dat2h",   "diag_test_18_file1", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Added attributes, and cell_measures (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=19
cat <<_EOF > diag_table
test_diag_manager_19
1 3 1 0 0 0

#output files
 "diag_test_19_file1", 1, "days", 1, "days", "time"
 "diag_test_19_file2", 1, "days", 1, "days", "time"

#output variables
"test_diag_manager_mod", "dat1",    "dat1",    "diag_test_19_file1", "all", .true., "none", 2
"test_diag_manager_mod", "dat2",    "dat2",    "diag_test_19_file1", "all", .true., "none", 2
"test_mod",              "dat2h",   "dat2h",   "diag_test_19_file1", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Area and Volume same field (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=20
cat <<_EOF > diag_table
test_diag_manager_20
1 3 1 0 0 0

#output files
 "diag_test_20_file1", 1, "days", 1, "days", "time"
 "diag_test_20_file2", 1, "days", 1, "days", "time"

#output variables
"test_diag_manager_mod", "dat1",    "dat1",    "diag_test_20_file1", "all", .true., "none", 2
"test_diag_manager_mod", "dat2",    "dat2",    "diag_test_20_file1", "all", .true., "none", 2
"test_mod",              "dat2h",   "dat2h",   "diag_test_20_file1", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Get diag_field_id, ID found and not found (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=21
cat <<_EOF > diag_table
test_diag_manager_21
1 3 1 0 0 0

#output files
 "diag_test_21_file1", 1, "days", 1, "days", "time"

#output variables
"test_diag_manager_mod", "dat1",              "dat1",           "diag_test_21_file1", "all", .true., "none", 2
"test_diag_manager_mod", "solar_constant",    "solar_constant", "diag_test_21_file1", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Add axis attributes (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=22
cat <<_EOF > diag_table
test_diag_manager_22
1 3 1 0 0 0

#output files
 "diag_test_22_file1", 1, "days", 1, "days", "time"

#output variables
"test_diag_manager_mod", "dat1",              "dat1",           "diag_test_22_file1", "all", .true., "none", 2
"test_diag_manager_mod", "solar_constant",    "solar_constant", "diag_test_22_file1", "all", .true., "none", 2
_EOF
setup_test
test_expect_success "Get 'nv' axis id (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

my_test_count=23
cat <<_EOF > diag_table
test_diag_manager_23
1990 1 1 0 0 0

#output files
"unstructured_diag_test", 2, "days", 2, "days", "time",

#output variables
"UG_unit_test", "unstructured_real_scalar_field_data", "rsf_diag_1", "unstructured_diag_test", "all", .TRUE., "none", 1,
"UG_unit_test", "unstructured_real_1D_field_data", "unstructured_real_1D_field_data", "unstructured_diag_test", "all", .TRUE., "none", 1,
"UG_unit_test", "unstructured_real_2D_field_data", "unstructured_real_2D_field_data", "unstructured_diag_test", "all", .TRUE., "none", 1,
"UG_unit_test", "lon", "grid_xt", "unstructured_diag_test", "all", .TRUE., "none", 1,
"UG_unit_test", "lat", "grid_yt", "unstructured_diag_test", "all", .TRUE., "none", 1,
_EOF
setup_test
test_expect_success "Unstructured grid (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager
'

# test_diag_manager_time
cat <<_EOF > diag_table
test_diag_manager
2 1 1 0 0 0

#output files
"ocn_default%4yr%2mo%2dy%2hr",      6,  "hours", 1, "hours", "time", 6, "hours", "2 1 1 0 0 0"
"ocn_middle%4yr%2mo%2dy%2hr",      6,  "hours", 1, "hours", "time", 6, "hours", "2 1 1 0 0 0", 0, "", "middle"
"ocn_begin%4yr%2mo%2dy%2hr",      6,  "hours", 1, "hours", "time", 6, "hours", "2 1 1 0 0 0" 0, "", "begin"
"ocn_end%4yr%2mo%2dy%2hr",      6,  "hours", 1, "hours", "time", 6, "hours", "2 1 1 0 0 0" 0, "", "end"

#output variables
 "test_diag_manager_mod", "sst", "sst", "ocn_default%4yr%2mo%2dy%2hr",  "all", .true., "none", 2
 "test_diag_manager_mod", "sst", "sst", "ocn_middle%4yr%2mo%2dy%2hr",  "all", .true., "none", 2
 "test_diag_manager_mod", "sst", "sst", "ocn_begin%4yr%2mo%2dy%2hr",  "all", .true., "none", 2
 "test_diag_manager_mod", "sst", "sst", "ocn_end%4yr%2mo%2dy%2hr",  "all", .true., "none", 2
_EOF

rm -f input.nml && touch input.nml
test_expect_success "wildcard filenames (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager_time
'

rm -f input.nml diag_table

touch input.nml
cat <<_EOF > diag_table
test_diag_manager
2 1 1 0 0 0

#output files
"test_diurnal",         1, "hours",   1, "hours", "time"

#output variables
 "test_diag_manager_mod", "sst", "sst", "test_diurnal",  "all", "diurnal3", "none", 2
 "test_diag_manager_mod", "ice", "ice", "test_diurnal",  "all", "diurnal3", "none", 2
_EOF
test_expect_success "diurnal test (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager_time
'

test_done
