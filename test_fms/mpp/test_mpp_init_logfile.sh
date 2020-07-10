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

# Set common test settings.
. ../test_common.sh

skip_test="no"

echo "Calling test_mp_init_logfile"

# Get the number of available CPUs on the system
if [ $(command -v nproc) ]
then
    # Looks like a linux system
    nProc=$(nproc)
elif [ $(command -v sysctl) ]
then
    # Looks like a Mac OS X system
    nProc=$(sysctl -n hw.physicalcpu)
else
    nProc=-1
fi

# Do we need to oversubscribe
if [ ${nProc} -lt 0 ]
then
    # Couldn't get the number of CPUs, skip the test.
    skip_test="skip"
elif [ $nProc -lt 4 ]
then
    # Need to oversubscribe the MPI
    run_test test_mpp_init_logfile 4 $skip_test "true"
fi

fprefix="logfile.00000"
file0="${fprefix}""0.out"
file1="${fprefix}""1.out"
file2="${fprefix}""2.out"
file3="${fprefix}""3.out"

#Create two log files, each with some content.
#echo "test_log_line" > logfile.000000.out
#echo "test_log_line" > logfile.000002.out
echo "test_log_line" > ${file0}
echo "test_log_line" > ${file2}
#Make sure other possible log files are not in the system.
if test -f ${file1}; then
    rm  ${file1}
fi
if test -f ${file3}; then
    rm  ${file3}
fi

#Have mpp re-initialize the two log files created above:
run_test test_mpp_init_logfile 4 $skip_test

#Make sure the two "old" ones have been replaced and the
#two possible new ones are not present.
#TODO test number
if [ $(wc -l < ${file0}) -ge 1  ] || [ $(wc -l < ${file2}) -ge 1  ] ||
       [ -f ${file1} ] || [ -f ${file3} ]
then
  echo "ERROR: Test ? was unsuccessful."
  exit 1
else
    echo "Test ? has passed"
    exit 0
fi
