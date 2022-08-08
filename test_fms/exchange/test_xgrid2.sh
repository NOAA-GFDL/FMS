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
