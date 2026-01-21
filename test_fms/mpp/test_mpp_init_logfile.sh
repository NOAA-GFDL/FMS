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
# execute tests in the test_fms/mpp directory.

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir
touch input.nml

fprefix="logfile.00000"
file0="${fprefix}""0.out"
file1="${fprefix}""1.out"
file2="${fprefix}""2.out"
file3="${fprefix}""3.out"
fcontent="test_log_line_88gd!ok5"

#Create two log files, each with some content.
echo ${fcontent} > ${file0}
echo ${fcontent} > ${file2}
#Make sure other possible log files are not in the system.
if test -f ${file1}; then
    rm  ${file1}
fi
if test -f ${file3}; then
    rm  ${file3}
fi

#Have mpp re-initialize the two log files created above:
test_expect_success "initialize mpp logfile" '
    mpirun -n 4 ../test_mpp_init_logfile
'

# Return sucess (0) only if the two "old" files have been replaced and the
# the two possible new ones are not present. Otherwise retun failure (1).
# Replacement is checked by the absence of the fcontent line.
if test $(grep  ${fcontent}  ${file0} | wc -l ) -ge 1 ||
   test $(grep  ${fcontent}  ${file1} | wc -l ) -ge 1 ||
   test -f ${file1} ||
   test -f ${file3} ;
then
  echo "ERROR: Test was unsuccessful."
  exit 1
else
    echo "Test has passed"
fi

test_done
