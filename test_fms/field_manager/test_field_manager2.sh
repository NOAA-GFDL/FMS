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

# Copy files for test.
cat <<_EOF > field_table
# Simplified field table to run the field table unit tests
 "TRACER", "ocean_mod", "biotic1"
           "diff_horiz", "linear", "slope=ok"
           "longname", "biotic one" /
 "TRACER", "ocean_mod", "age_ctl" /
 "TRACER", "atmos_mod","radon"
           "longname","radon-222"
           "units","VMR*1E21"
           "profile_type","fixed","surface_value=0.0E+00"
           "convection","all"/
 "TRACER", "land_mod", "sphum"
           "longname",     "specific humidity"
            "units",        "kg/kg" /
_EOF

cat <<_EOF > input.nml
&test_field_manager

/
_EOF

test_expect_success "field manager functional" '
  mpirun -n 2 ./test_field_manager
'

test_done
