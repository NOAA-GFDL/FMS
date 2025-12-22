#!/bin/sh

#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# Tests to skip if input files not present
if test -z "$test_input_path" ; then
  SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).1"
else
  rm -rf INPUT && mkdir INPUT
  cp $test_input_path/exchange/INPUT/* INPUT
fi

# Copy file for test.
cat <<_EOF > input.nml
&xgrid_nml
  interp_method = 'second_order'
  make_exchange_reproduce=.true.
/

&xgrid_test_nml
  test_unstruct = .false.
  atm_layout=2,2
  lnd_layout= 2,2
  ice_layout=12,3
  nk_lnd = 14
  nk_ice = 6
  atm_field_name = 'depth'
  test_unstruct = .true.
/
_EOF

test_expect_success "Test exchange grid" '
  mpirun -n 12 ./test_xgrid
'

rm -rf INPUT *.nc # remove any leftover io files to save space

test_done
