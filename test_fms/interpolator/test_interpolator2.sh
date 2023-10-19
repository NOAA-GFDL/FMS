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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/field_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# Tests to skip if input files not present
if test ! -z "$test_input_path" ; then
  rm -rf INPUT && mkdir INPUT
  cp $test_input_path/interpolator/INPUT/* INPUT
else
  SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).1"
fi

# Create files for test.
cat <<_EOF  > diag_table
test_diag_manager_01
1 3 1 0 0 0

#output files
 "diag_test_01",  1, "days", 1, "days", "time"

#output variables
 "test_diag_manager_mod", "dat1", "dat1", "diag_test_01",  "all", .false., "none", 2
_EOF

cat <<_EOF > input.nml
&interpolator_nml
/
_EOF

# Run test
test_expect_success "test interpolator" 'mpirun -n 1 ./test_interpolator'


#Run the daily interpolator tests when the file calendar is in units of days and calendar type is NOLEAP
cat <<EOF > test_interpolator.nml
&test_interpolator_nml
test_file_daily_noleap=.false.
test_file_daily_julian=.true.
test_file_yearly_noleap=.false.
test_file_yearly_julian=.false.
test_file_no_time=.false.
/
EOF
mkdir -p INPUT
test_expect_success "test_interpolator2 file data daily julian r4 unit tests" 'mpirun -n 1 ./test_interpolator2_r4'
test_expect_success "test_interpolator2 file data daily julian r8 unit tests" 'mpirun -n 1 ./test_interpolator2_r8'
rm -rf INPUT *.nc test_interpolator.nml


#Run the daily interpolator tests when the file calendar is in units of days and calendar type is JULIAN
cat <<EOF > test_interpolator.nml
&test_interpolator_nml
test_file_daily_noleap=.true.
test_file_daily_julian=.false.
test_file_yearly_noleap=.false.
test_file_yearly_julian=.false.
test_file_no_time=.false.
/
EOF
mkdir -p INPUT
test_expect_success "test_interpolator2 file data daily noleap r4 unit tests" 'mpirun -n 1 ./test_interpolator2_r4'
test_expect_success "test_interpolator2 file data daily noleap r8 unit tests" 'mpirun -n 1 ./test_interpolator2_r8'
rm -rf INPUT *.nc test_interpolator.nml


#Run the yearly interpolator tests when the file calendar is in units of years and calendar type is NOLEAP
cat <<EOF > test_interpolator.nml
&test_interpolator_nml
test_file_daily_noleap=.false.
test_file_daily_julian=.false.
test_file_yearly_noleap=.true.
test_file_yearly_julian=.false.
test_file_no_time=.false.
/
EOF
mkdir -p INPUT
test_expect_success "test_interpolator2 file data yearly noleap r4 unit tests" 'mpirun -n 1 ./test_interpolator2_r4'
test_expect_success "test_interpolator2 file data yearly noleap r8 unit tests" 'mpirun -n 1 ./test_interpolator2_r8'
rm -rf INPUT *.nc test_interpolator.nml


#Run the yearly interpolator tests when the file calendar is in units of years and calendar type is JULIAN
cat <<EOF > test_interpolator.nml
&test_interpolator_nml
test_file_daily_noleap=.false.
test_file_daily_julian=.false.
test_file_yearly_noleap=.false.
test_file_yearly_julian=.true.
test_file_no_time=.false.
/
EOF
mkdir -p INPUT
test_expect_success "test_interpolator2 file data yearly julian r4 unit tests" 'mpirun -n 1 ./test_interpolator2_r4'
test_expect_success "test_interpolator2 file data yearly julian r8 unit tests" 'mpirun -n 1 ./test_interpolator2_r8'
rm -rf INPUT *.nc test_interpolator.nml


#Run no_time_axis
cat <<EOF > test_interpolator.nml
&test_interpolator_nml
test_file_daily_noleap=.false.
test_file_daily_julian=.false.
test_file_yearly_noleap=.false.
test_file_yearly_julian=.false.
test_file_no_time=.true.
/
EOF
mkdir -p INPUT
test_expect_success "test_interpolator2 file data no time axis r4 unit tests" 'mpirun -n 1 ./test_interpolator2_r4'
test_expect_success "test_interpolator2 file data no time axis r8 unit tests" 'mpirun -n 1 ./test_interpolator2_r8'
rm -rf INPUT *.nc test_interpolator.nml


test_done
