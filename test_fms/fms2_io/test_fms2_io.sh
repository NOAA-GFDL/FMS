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
# execute tests in the test_fms/fms2_io directory.

# Authors: Raymond Menzel
# Jessica Liptak
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

# Warn if stack size not set
if test "`ulimit -s`" != "unlimited"; then
  say "Warning: system stack size not set to unlimited, test may fail"
fi

# make a dummy input.nml file for mpp_init to read
cat <<_EOF > input.nml
&dummy
  empty=true
/
_EOF

# run the tests
test_expect_success "FMS2 IO Test" '
  mpirun -n 6 ../test_fms2_io
'

test_done
