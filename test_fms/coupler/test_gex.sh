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

# Set common test settings.
. ../test-lib.sh

output_dir


test_cmd="mpirun -n 1 ../test_gex"

# Create input.nml and field table (legacy field table)
prepare_legacy () {
  test_name=$1

  rm -f field_table.yaml

  cat <<EOF >field_table
"atm_to_lnd_ex", "coupler_mod", "dryoa"
"units=kg/m2/s" /
EOF

  cat <<EOF >input.nml
&field_manager_nml
  use_field_table_yaml = .false.
/

&test_gex_nml
  test_name = "$test_name"
/
EOF
}

# Create input.nml and field table (YAML field table)
prepare_yaml () {
  test_name=$1

  rm -f field_table

  cat <<EOF >field_table.yaml
field_table:
- field_type: atm_to_lnd_ex
  modlist:
  - model_type: coupler_mod
    varlist:
    - variable: dryoa
      units: kg/m2/s
EOF

  cat <<EOF >input.nml
&field_manager_nml
  use_field_table_yaml = .true.
/

&test_gex_nml
  test_name = "$test_name"
/
EOF
}

# Run the tests.

prepare_legacy default_test
test_expect_success "Test gex with atm_to_lnd tracer (legacy field_table)" "$test_cmd"

if [ -z "$parser_skip" ]; then
  prepare_yaml default_test
  test_expect_success "Test gex with atm_to_lnd tracer (YAML field_table)" "$test_cmd"
fi

prepare_legacy get_n_ex_invalid_model_src
test_expect_failure "Test gex_get_n_ex with invalid model_src" "$test_cmd"

prepare_legacy get_n_ex_invalid_model_rec
test_expect_failure "Test gex_get_n_ex with invalid model_rec" "$test_cmd"

prepare_legacy get_property_invalid_tracer
test_expect_failure "Test gex_get_property with invalid tracer ID" "$test_cmd"

prepare_legacy get_property_invalid_property
test_expect_failure "Test gex_get_property with invalid property ID" "$test_cmd"

test_done
