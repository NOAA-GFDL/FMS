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
. ../test_common.sh

# Create the base input.nml file needed for the tests
cat <<_EOF > input.nml
&test_fms_io_nml
  io_layout = 1,1
/
_EOF

# Test the structured grid
rm -rf RESTART && mkdir RESTART
run_test test_fms_io 6

# Ensure the restart directory is empty
rm -rf RESTART && mkdir RESTART
run_test test_unstructured_fms_io 6
