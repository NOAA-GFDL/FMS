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
# execute tests in the test_fms/offloading directory.

# Set common test settings.
. ../test-lib.sh

cat > input.nml << _EOF
_EOF

test_expect_success "io offload with one model PE and one writer" '
  mpirun -n 2 ./test_io_offloading
'

cat > input.nml << _EOF
&offloading_io_nml
  max_files = 2
/
&test_io_offloading_nml
  mode = "normal"
/
_EOF

test_expect_success "io offload with 4 model pes and one writer" '
  rm -f offloaded_output*.nc && mpirun -n 5 ./test_io_offloading
'

cat > input.nml << _EOF
&offloading_io_nml
  max_files = 1
/
_EOF

test_expect_failure "io offload rejects more open files than max_files" '
  mpirun -n 2 ./test_io_offloading
'

# check error on bad input during register/attribute routines
for mode in bad-axis bad-dims bad-attribute bad-level; do
  cat > input.nml << _EOF
&test_io_offloading_nml
  mode = "$mode"
/
_EOF

  test_expect_failure "io offload rejects $mode" '
    mpirun -n 2 ./test_io_offloading
  '
done

cat > input.nml << _EOF
_EOF

test_expect_success "cubic io offload with 6 model pes and 6 writers" '
  mpirun -n 12 ./test_io_offloading_cubic
'

cat > input.nml << _EOF
&test_io_offloading_cubic_nml
  max_pending_writes = 1
/
_EOF

test_expect_success "cubic io offload with 12 model pes, 6 writers and one pending write" '
  mpirun -n 18 ./test_io_offloading_cubic
'

cat > input.nml << _EOF
&test_io_offloading_cubic_nml
  nwriters = 1
/
_EOF

test_expect_failure "cubic io offload rejects a writer count that is not a multiple of the tile count" '
  mpirun -n 7 ./test_io_offloading_cubic
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 1
  nwriters = 4
/
_EOF

test_expect_success "single tile io offload with 4 writers into one file" '
  rm -f offloaded_multi*.nc* && mpirun -n 8 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 1
  nwriters = 4
  writer_io_layout = 2, 2
/
_EOF

test_expect_success "single tile io offload with 4 writers into distributed files" '
  rm -f offloaded_multi*.nc* && mpirun -n 8 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 6
  nwriters = 12
/
_EOF

test_expect_success "cubic io offload with 2 writers per tile" '
  rm -f offloaded_multi*.nc* && mpirun -n 18 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 6
  nwriters = 12
  writer_layout = 1, 2
  writer_io_layout = 1, 2
/
_EOF

test_expect_success "cubic io offload with 2 writers per tile into distributed files" '
  rm -f offloaded_multi*.nc* && mpirun -n 18 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 6
  nwriters = 24
  writer_layout = 2, 2
/
_EOF

test_expect_success "cubic io offload with 4 writers per tile and 12 model pes" '
  rm -f offloaded_multi*.nc* && mpirun -n 36 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 6
  nwriters = 9
/
_EOF

test_expect_failure "cubic io offload rejects 9 writers for 6 tiles" '
  rm -f offloaded_multi*.nc* && mpirun -n 15 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
&test_io_offloading_multi_writer_nml
  ntiles = 6
  nwriters = 12
  writer_layout = 1, 2
  writer_io_layout = 2, 1
/
_EOF

test_expect_failure "io offload rejects a writer io_layout that does not divide the writer layout" '
  rm -f offloaded_multi*.nc* && mpirun -n 18 ./test_io_offloading_multi_writer
'

cat > input.nml << _EOF
_EOF

test_expect_success "test metadata transfer" '
  mpirun -n 4 ./test_metadata_transfer
'

test_done
