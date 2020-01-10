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

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

#bats $srcdir/test_fms/mpp/test_mpp_domains.bats

if [ "x$(uname -s)" = "xDarwin" ]
then
    skip_test=true
else
    skip_test=false
fi

#echo "1: Test update nest domain"
#sed "s/test_nest_domain = .false./test_nest_domain = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#run mpirun -n 2 ./test_mpp_domains

#echo "2:  Test Subset Update"
#sed "s/test_subset = .false./test_subset = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#run mpirun -n 2 ./test_mpp_domains

echo "3: Test Halosize Performance"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

#echo "4: Test Edge Update"
#sed "s/test_edge_update = .false./test_edge_update = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#mpirun -n 2 ./test_mpp_domains

#echo "5: Test Nonsym Edge"
#sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#mpirun -n 2 ./test_mpp_domains

echo "6: Test Performance"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
elif [ "x$TRAVIS" = "xtrue" ]
then
    echo "Fails on Travis"
else
    sed "s/test_performance = .false./test_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 6 ./test_mpp_domains
fi

echo "7: Test Global Sum"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_global_sum = .false./test_global_sum = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

echo "8: Test Cubic Grid Redistribute"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
elif [ "x$TRAVIS" = "xtrue" ]
then
    echo "Fails on Travis"
else
    sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 6 ./test_mpp_domains
fi

echo "9: Test Boundary"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
elif [ "x$TRAVIS" = "xtrue" ]
then
    echo "Fails on Travis in 32bit/Mixed mode"
else
    sed "s/test_boundary = .false./test_boundary = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

echo "10: Test Adjoint"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_adjoint = .false./test_adjoint = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

#echo "11: Test Unstruct"
#sed "s/test_unstruct = .false./test_unstruct = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#mpirun -n 2 ./test_mpp_domains

echo "12: Test Group"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_group = .false./test_group = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

echo "13: Test Interface"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_interface = .false./test_interface = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi

#echo "14: Test Check Parallel"
#echo "Does not work on Darwin or elsewhere"
#sed "s/check_parallel = .false./check_parallel = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#mpirun -n 6 ./test_mpp_domains

echo "15: Test Get Nbr"
if [ "$skip_test" = "true" ]
then
    echo "Does not work on Darwin"
else
    sed "s/test_get_nbr = .false./test_get_nbr = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
    mpirun -n 2 ./test_mpp_domains
fi
