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
# execute tests in the test_fms/parser directory.

# Set common test settings.
. ../test_common.sh

touch input.nml

cat <<_EOF > data_table.yaml
data_table:
 - gridname          : "ICE"
   fieldname_code    : "sic_obs"
   fieldname_file    : "ice"
   file_name         : "INPUT/hadisst_ice.data.nc"
   interpol_method   : "bilinear"
   factor            : 0.01
 - gridname          : "WUT"
   fieldname_code    : "potato"
   fieldname_file    : "mullions"
   file_name         : "INPUT/potato.nc"
   interpol_method   : "bilinear"
   factor            : 1e-06
   region_type       : "inside_region"
   lat_start         : -89.1
   lat_end           : 89.8
   lon_start         : 3.4
   lon_end           : 154.4
   do_data_bug       : false
   use_data_bug      : True
_EOF

cat <<_EOF > diag_table.yaml
title: c384L49_esm5PIcontrol
baseDate: [1960 1 1 1 1 1 1]
diag_files:
-    fileName: "atmos_daily"
     freq: 24
     frequnit: hours
     timeunit: days
     unlimdim: time
     varlist:
     - varName: tdata
       reduction: False
       module: mullions
       mullions: 10
       fill_value: -999.9
     - varName: pdata
       outName: pressure
       reduction: False
       kind: double
       module: "moist"
-    fileName: atmos_8xdaily
     freq: 3
     frequnit: hours
     timeunit: days
     unlimdim: time
     varlist:
     - varName: tdata
       reduction: False
       module: "moist"
_EOF

run_test test_yaml_parser 1 $parser_skip
run_test parser_demo 1 $parser_skip
run_test parser_demo2 1 $parser_skip

printf "&check_crashes_nml \n missing_file = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n bad_conversion = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n missing_key = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_block_ids_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_num_blocks_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_nkeys_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_key_ids_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_key_name_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_key_value_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_value_from_key_bad_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_key_name_bad_key_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_key_value_bad_key_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

###
printf "&check_crashes_nml \n get_key_ids_bad_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_nkeys_bad_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_block_ids_bad_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_num_blocks_bad_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n get_value_from_key_bad_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 $parser_skip && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n wrong_buffer_size_key_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi

printf "&check_crashes_nml \n wrong_buffer_size_block_id = .true. \n/" | cat > input.nml
run_test check_crashes 1 && echo "It worked?"
if [ $? -eq 0 ]; then
  echo "The test should have failed"
  exit 3
fi
