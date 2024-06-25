#!/bin/sh

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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ryan Mulhall 02/2021

# Set common test settings.
. ../test-lib.sh

output_dir

# ensure input.nml file present
cat <<_EOF  > input.nml
&test_stdlog_nml
  test_num = 1
/
_EOF
# Run test with one processor
test_expect_success "test stdlog and stdwarn" '
    mpirun -n 2 ../test_stdlog
'
sed -i 's/1/2/' input.nml
test_expect_failure "test stdlog and stdwarn with fatal output" '
    mpirun -n 2 ../test_stdlog
'
# move file so we don't overwrite
mv warnfile.*.out warnfile.000000.out.old
sed -i 's/2/3/' input.nml
test_expect_success "check stdwarn output" '
    mpirun -n 1 ../test_stdlog
'
test_done
