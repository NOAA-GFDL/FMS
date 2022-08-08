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
test_expect_success "test interpolator" '
    mpirun -n 2 ./test_interpolator
'

rm -rf INPUT *.nc # remove any leftover io files to save space

test_done
