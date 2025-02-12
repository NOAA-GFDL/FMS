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
   do_diag_field_log=.true.
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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

my_test_count=`expr $my_test_count + 1`
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
my_test_count=`expr $my_test_count + 1`
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
 "test_diag_manager_mod", "sst", "sst", "test_diurnal",  "all", "diurnal4", "none", 2
 "test_diag_manager_mod", "ice", "ice", "test_diurnal",  "all", "diurnal4", "none", 2
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_success "diurnal test (test $my_test_count)" '
  mpirun -n 1 ../test_diag_manager_time
'
setup_test
my_test_count=`expr $my_test_count + 1`
test_expect_success "Test the diag update_buffer (test $my_test_count)" '
  mpirun -n 1 ../test_diag_update_buffer
'

## run tests that are ifdef'd out only if compiled with yaml
## otherwise just run the updated end to end to check for error
if [ -z "${skipflag}" ]; then

  cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: wild_card_name%4yr%2mo%2dy%2hr
  filename_time: end
  freq: 6 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours
  module: test_diag_manager_mod
  reduction: average
  kind: r4
  varlist:
  - var_name: sst
    output_name: sst
  global_meta:
  - is_a_file: true
- file_name: normal
  freq: 24 days
  time_units: hours
  unlimdim: records
# Here the module, kind and reduction are being overwritten with whats on the variable
  module: potato_mod
  kind: r8
  reduction: min
  varlist:
  - module: test_diag_manager_mod
    var_name: sst
    output_name: sst
    reduction: average
    kind: r4
    write_var: true
    attributes:
    - do_sst: .true.
  sub_region:
  - grid_type: latlon
    corner1: -80, 0
    corner2: -80, 75
    corner3: -60, 0
    corner4: -60, 75
- file_name: normal2
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: true
  module: test_diag_manager_mod
  reduction: none
  kind: r4
  varlist:
  - var_name: sstt
    output_name: sstt
    long_name: S S T
  - var_name: sstt2
    output_name: sstt2
    long_name: S S T
    write_var: false
  sub_region:
  - grid_type: index
    tile: 1
    corner1: 10, 15
    corner2: 20, 15
    corner3: 10, 25
    corner4: 20, 25
- file_name: normal3
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: false
_EOF

my_test_count=`expr $my_test_count + 1`
  test_expect_success "diag_yaml test with the simple diag table.yaml (test $my_test_count)" '
    mpirun -n 1 ../test_diag_yaml
  '

  cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: wild_card_name%4yr%2mo%2dy%2hr
  filename_time: end
  freq: 6 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours
  varlist:
  - module: test_diag_manager_mod
    var_name: sst
    output_name: sst
    reduction: average
    kind: r4
  global_meta:
  - is_a_file: true
- file_name: normal
  freq: 24 days
  time_units: hours
  unlimdim: records
  varlist:
  - module: test_diag_manager_mod
    var_name: sst
    output_name: sst
    reduction: average
    kind: r4
    write_var: true
    attributes:
    - do_sst: .true.
  sub_region:
  - grid_type: latlon
    corner1: -80, 0
    corner2: -80, 75
    corner3: -60, 0
    corner4: -60, 75
- file_name: normal2
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: true
  varlist:
  - module: test_diag_manager_mod
    var_name: sstt
    output_name: sstt
    reduction: none
    kind: r4
    long_name: S S T
  - module: test_diag_manager_mod
    var_name: sstt2
    output_name: sstt2
    reduction: none
    kind: r4
    long_name: S S T
    write_var: false
  sub_region:
  - grid_type: index
    tile: 1
    corner1: 10, 15
    corner2: 20, 15
    corner3: 10, 25
    corner4: 20, 25
- file_name: normal3
  freq: -1
  time_units: hours
  unlimdim: records
  write_file: false
_EOF
  cp diag_table.yaml diag_table.yaml_base

  my_test_count=`expr $my_test_count + 1`
  test_expect_success "diag_yaml test (test $my_test_count)" '
    mpirun -n 1 ../test_diag_yaml
  '
  . $top_srcdir/test_fms/diag_manager/check_crashes.sh
  my_test_count=`expr $my_test_count + 14`

  printf "&diag_manager_nml \n use_modern_diag = .true. \n/" | cat > input.nml
  cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: file1
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: test_diag_manager_mod
    var_name: sst1
    output_name: sst1
    reduction: average
    kind: r4
- file_name: file2
  freq: 6 hours
  time_units: hours
  unlimdim: time
  is_ocean: True
  varlist:
  - module: test_diag_manager_mod
    var_name: sst2
    output_name: sst2
    reduction: average
    kind: r4
- file_name: file3
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: test_diag_manager_mod
    var_name: sst3
    output_name: sst3
    reduction: average
    kind: r4
  - module: test_diag_manager_mod
    var_name: sst4
    output_name: sst4
    reduction: average
    kind: r4
_EOF

  my_test_count=`expr $my_test_count + 1`
  test_expect_success "Test the diag_ocean feature in diag_manager_init (test $my_test_count)" '
    mpirun -n 2 ../test_diag_ocean
  '


  printf "&diag_manager_nml \n use_modern_diag = .true. \n do_diag_field_log = .true. \n/" | cat > input.nml
  cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0

diag_files:
- file_name: static_file
  freq: -1
  time_units: hours
  unlimdim: time
  varlist:
  - module: atm_mod
    var_name: var7
    reduction: none
    kind: r4
  global_meta:
  - is_important: False
    has_important: True
- file_name: file1
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
  - module: ocn_mod
    var_name: var2
    output_name: potato
    reduction: average
    kind: r4
- file_name: file2
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: atm_mod
    var_name: var3
    reduction: average
    kind: r4
  - module: atm_mod
    var_name: var4
    output_name: i_on_a_sphere
    reduction: average
    kind: r8
  - module: atm_mod
    var_name: var6
    reduction: average
    kind: r8
  - module: atm_mod
    var_name: var4
    output_name: var4_bounded
    reduction: average
    kind: r8
    zbounds: 2.0 3.0
- file_name: file3
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: lnd_mod
    var_name: var5
    reduction: average
    kind: r4
  - module: atm_mod
    var_name: var7
    reduction: none
    kind: r4
- file_name: file4
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: lnd_mod
    var_name: var1
    reduction: average
    kind: r4
- file_name: file5
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: atm_mod
    var_name: var4
    reduction: average
    kind: r4
  sub_region:
  - grid_type: index
    tile: 1
    corner1: 10, 15
    corner2: 20, 15
    corner3: 10, 25
    corner4: 20, 25
- file_name: file6%4yr%2mo%2dy%2hr
  freq: 6 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
- file_name: file7
  freq: 6 hours
  time_units: hours
  unlimdim: time
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: none
    kind: r4
    attributes:
    - GFDL_name: var_var
- file_name: file8%4yr%2mo%2dy%2hr%2min
  freq: 1 hours,1 hours,1 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours, 3 hours, 1 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours, 3 hours, 9 hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
- file_name: file9%4yr%2mo%2dy%2hr%2min
  filename_time: begin
  freq: 1 hours,1 hours,1 hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 hours, 3 hours, 1 hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 hours, 3 hours, 9 hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
- file_name: file10_diurnal
  freq: 1 days
  time_units: hours
  unlimdim: time
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: diurnal12
    kind: r4
_EOF

  my_test_count=`expr $my_test_count + 1`
  test_expect_success "buffer functionality (test $my_test_count)" '
    mpirun -n 1 ../test_diag_buffer
  '

  my_test_count=`expr $my_test_count + 1`
  test_expect_success "Test the modern diag manager end to end (test $my_test_count)" '
    mpirun -n 6 ../test_modern_diag
  '

## print out a reference for the yaml output test, just uses the last diag table created
  cat <<_EOF > diag_out_ref.yaml
---
title: test_diag_manager
base_date: 2 1 1 0 0 0
diag_files:
- file_name: static_file
  freq: -1
  freq_units: days
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: atm_mod
    var_name: var7
    reduction: none
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: z
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - is_important: False
    has_important: True
- file_name: file1
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y x
  - module: ocn_mod
    var_name: var2
    reduction: average
    kind: r4
    output_name: potato
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time x y
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file2
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: atm_mod
    var_name: var3
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y3 x3
  - module: atm_mod
    var_name: var4
    reduction: average
    kind: r8
    output_name: i_on_a_sphere
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time z y3 x3
  - module: atm_mod
    var_name: var6
    reduction: average
    kind: r8
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time z
  - module: atm_mod
    var_name: var4
    reduction: average
    kind: r8
    output_name: var4_bounded
    long_name:
    units:
    zbounds: 2.00 3.00
    n_diurnal:
    pow_value:
    dimensions: time z_sub01 y3 x3
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file3
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: lnd_mod
    var_name: var5
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time grid_index
  - module: atm_mod
    var_name: var7
    reduction: none
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: z
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file4
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: lnd_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file5
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: atm_mod
    var_name: var4
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time z y3_sub01 x3_sub01
  sub_region:
  - grid_type: index
    tile: 1
    corner1: 10 15
    corner2: 20 15
    corner3: 10 25
    corner4: 20 25
  global_meta:
  - {}
- file_name: file6%4yr%2mo%2dy%2hr
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6
  new_file_freq_units: hours
  start_time: 2 1 1 0 0 0
  file_duration: 12
  file_duration_units: hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y x
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file7
  freq: 6
  freq_units: hours
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: none
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y x
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file8%4yr%2mo%2dy%2hr%2min
  freq: 1 1 1
  freq_units: hours hours hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 3 1
  new_file_freq_units: hours hours hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 3 9
  file_duration_units: hours hours hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y x
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file9%4yr%2mo%2dy%2hr%2min
  freq: 1 1 1
  freq_units: hours hours hours
  time_units: hours
  unlimdim: time
  new_file_freq: 6 3 1
  new_file_freq_units: hours hours hours
  start_time: 2 1 1 0 0 0
  file_duration: 12 3 9
  file_duration_units: hours hours hours
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal:
    pow_value:
    dimensions: time y x
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
- file_name: file10_diurnal
  freq: 1
  freq_units: days
  time_units: hours
  unlimdim: time
  new_file_freq:
  new_file_freq_units:
  start_time:
  file_duration:
  file_duration_units:
  varlist:
  - module: ocn_mod
    var_name: var1
    reduction: diurnal
    kind: r4
    output_name:
    long_name:
    units:
    zbounds:
    n_diurnal: 12
    pow_value:
    dimensions: time time_of_day_12 y x
  sub_region:
  - grid_type:
    tile:
    corner1:
    corner2:
    corner3:
    corner4:
  global_meta:
  - {}
...
_EOF

my_test_count=`expr $my_test_count + 1`
test_expect_success "check modern diag manager yaml output (test $my_test_count)" '
    mpirun -n 1 ../test_diag_out_yaml
'

printf "&diag_manager_nml \n use_modern_diag = .true. \n use_clock_average = .true. \n /" | cat > input.nml
cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0

diag_files:
- file_name: file1_clock
  freq: 1 days
  time_units: hours
  unlimdim: time
  varlist:
  - module: atm_mod
    var_name: var1
    reduction: average
    kind: r4
_EOF

my_test_count=`expr $my_test_count + 1`
  test_expect_success "Test the modern diag manager with use_clock_average = .true. (test $my_test_count)" '
    mpirun -n 1 ../test_flexible_time
  '

printf "&diag_manager_nml \n use_modern_diag = .true. \n use_clock_average = .false. \n /" | cat > input.nml
cat <<_EOF > diag_table.yaml
title: test_diag_manager
base_date: 2 1 1 0 0 0

diag_files:
- file_name: file1_forecast
  freq: 1 days
  time_units: hours
  unlimdim: time
  varlist:
  - module: atm_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name: var1_min
  - module: atm_mod
    var_name: var1
    reduction: average
    kind: r4
    output_name: var2_max
_EOF

my_test_count=`expr $my_test_count + 1`
  test_expect_success "Test the modern diag manager with use_clock_average = .false. (test $my_test_count)" '
    mpirun -n 1 ../test_flexible_time
  '
printf "&diag_manager_nml \n use_modern_diag = .false. \n use_clock_average = .true. \n /" | cat > input.nml
  test_expect_failure "Test if use_modern_diag = .false. and use_clock_average = .true. fails (test $my_test_count)" '
    mpirun -n 1 ../test_flexible_time
  '

else
  my_test_count=`expr $my_test_count + 1`
  test_expect_failure "test modern diag manager failure when compiled without -Duse-yaml flag (test $my_test_count)" '
    mpirun -n 6 ../test_modern_diag
  '
fi
test_done
