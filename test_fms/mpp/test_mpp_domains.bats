setup () {
  if [ "x$(uname -s)" = "xDarwin" ]
  then
    skip_test=true
  else
    skip_test=false
  fi
}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1: Test update nest domain" {
    skip "To do"
    sed "s/test_nest_domain = .false./test_nest_domain = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "2:  Test Subset Update" {
    skip "To do"
    sed "s/test_subset = .false./test_subset = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "3: Test Halosize Performance" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "4: Test Edge Update" {
    skip "To do"
    sed "s/test_edge_update = .false./test_edge_update = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "5: Test Nonsym Edge" {
    skip "To do"
    sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "6: Test Performance" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    elif [ "x$TRAVIS" = "xtrue" ]
    then
      skip "Fails on Travis"
    fi
    sed "s/test_performance = .false./test_performance = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "7: Test Global Sum" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_global_sum = .false./test_global_sum = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "8: Test Cubic Grid Redistribute" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    elif [ "x$TRAVIS" = "xtrue" ]
    then
      skip "Fails on Travis"
    fi
    sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "9: Test Boundary" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    elif [ "x$TRAVIS" = "xtrue" ]
    then
      skip "Fails on Travis in 32bit/Mixed mode"
    fi
    sed "s/test_boundary = .false./test_boundary = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "10: Test Adjoint" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_adjoint = .false./test_adjoint = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "11: Test Unstruct" {
    skip "To do"
    sed "s/test_unstruct = .false./test_unstruct = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "12: Test Group" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_group = .false./test_group = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "13: Test Interface" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_interface = .false./test_interface = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "14: Test Check Parallel" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    elif [ "x$TRAVIS" = "xtrue" ]
    then
      skip "Fails on Travis"
    fi
    sed "s/check_parallel = .false./check_parallel = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "15: Test Get Nbr" {
    if [ "$skip_test" = "true" ]
    then
      skip "Does not work on Darwin"
    fi
    sed "s/test_get_nbr = .false./test_get_nbr = .true./" input.nml_base > input.nml
    run mpirun -n 2 ./test_mpp_domains
    [ "$status" -eq 0 ]
}
