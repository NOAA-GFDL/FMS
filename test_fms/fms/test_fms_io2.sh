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
# execute tests in the test_fms/fms directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

setup_test () {
  rm -rf RESTART && mkdir -p RESTART

  cat <<_EOF > input.nml
&test_fms_io_nml
  io_layout = 1,1
/
_EOF
}

# Test the structured grid
setup_test
test_expect_success "test structured grid" '
  mpirun -n 6 ./test_fms_io
'

# Ensure the restart directory is empty
setup_test
test_expect_success "test unstructured grid" '
  mpirun -n 6 ./test_unstructured_fms_io
'

test_done
