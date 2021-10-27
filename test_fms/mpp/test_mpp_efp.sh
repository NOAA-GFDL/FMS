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

# Ryan Mulhall 10/27/21 

# Set common test settings.
. ../test_common.sh

# makes directory to avoid other tests overwriting input.nml
mkdir test_mpp_efp_output && cd test_mpp_efp_output

cat <<_EOF > input.nml
&test_mpp_efp_nml
test_num = 1
/
_EOF

run_test ../test_mpp_efp 4
sed -i 's/test_num = 1/test_num = 2/' input.nml
run_test ../test_mpp_efp 4
sed -i 's/test_num = 2/test_num = 3/' input.nml
run_test ../test_mpp_efp 4

cd .. && rm -rf test_mpp_efp_output
