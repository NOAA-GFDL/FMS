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

rm -rf input.nml field_table field_table.yaml
touch input.nml

# Copy files for test.
cat <<_EOF > field_table
# Simplified field table to run the field table unit tests
 "TRACER", "ocean_mod", "biotic1"
           "diff_horiz", "linear", "slope=ok"
           "longname", "biotic one" /
 "TRACER", "ocean_mod", "age_ctl" /
 "TRACER", "ocean_mod" "immadeup2"
           "longname", "im_made_up2_for_testing"
           "units", "atomic_units"
           "profile_type", "profile", "surface_value=1.e-12,bottom_value=1.e-9"/
 "TRACER", "atmos_mod","radon"
           "longname","radon-222"
           "units","VMR*1E21"
           "profile_type","fixed","surface_value=0.0E+00"
           "convection","all"/
 "TRACER", "atmos_mod" "immadeup"
           "longname", "im_made_up_for_testing"
           "units", "hbar"
           "profile_type", "profile", "surface_value=1.e-12,top_value=1.e-15"/
 "TRACER", "land_mod", "sphum"
           "longname",     "specific humidity"
            "units",        "kg/kg" /
_EOF

test_expect_success "tracer_manager r4 with the legacy field table" 'mpirun -n 2 ./test_tracer_manager_r4'
test_expect_success "tracer_manager r8 with the legacy field table" 'mpirun -n 2 ./test_tracer_manager_r8'

if [ $parser_skip ]; then
rm -rf field_table
cat <<_EOF > input.nml
&field_manager_nml
  use_field_table_yaml = .true.
/
_EOF

cat <<_EOF > field_table.yaml
field_table:
- field_type: tracer
  modlist:
  - model_type: atmos_mod
    varlist:
    - variable: radon
      longname: radon-222
      units: VMR*1E21
      profile_type:
      - value: fixed
        surface_value: 0.0e+00
      convection: all
  - model_type: atmos_mod
    varlist:
    - variable: immadeup
      longname: im_made_up_for_testing
      units: hbar
      profile_type:
      - value: profile
        surface_value: 1.02e-12
        top_value: 1.0e-15
  - model_type: ocean_mod
    varlist:
    - variable: biotic1
      diff_horiz:
      - value: linear
        slope: ok
      longname: biotic one
    - variable: age_ctl
  - model_type: ocean_mod
    varlist:
    - variable: immadeup2
      longname: im_made_up2_for_testing
      units: hbar
      profile_type:
      - value: profile
        surface_value: 1.0e-12
        bottom_value: 1.0e-9
  - model_type: land_mod
    varlist:
    - variable: sphum
      longname: specific humidity
      units: kg/kg
_EOF

test_expect_success "tracer_manager r4 with the yaml field table" 'mpirun -n 2 ./test_tracer_manager_r4'
test_expect_success "tracer_manager r8 with the yaml field table" 'mpirun -n 2 ./test_tracer_manager_r8'
fi

test_done
