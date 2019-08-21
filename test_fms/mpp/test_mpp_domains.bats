@test "1: Test update nest domain" {
    skip "To do"
    sed "s/test_nest_domain = .false./test_nest_domain = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "2:  Test Subset Update" {
    skip "To do"
    sed "s/test_subset = .false./test_subset = .true./" input.nml_base > input.nml
    run mpirun -n 26 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "3: Test Halosize Performance" {
    sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "4: Test Edge Update" {
    sed "s/test_edge_update = .false./test_edge_update = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "5: Test Nonsym Edge" {
    skip "To do"
    sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "6: Test Performance" {
    sed "s/test_performance = .false./test_performance = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "7: Test Global Sum" {
    sed "s/test_global_sum = .false./test_global_sum = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "8: Test Cubic Grid Redistribute" {
    sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "9: Test Boundary" {
    sed "s/test_boundary = .false./test_boundary = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "10: Test Adjoint" {
    sed "s/test_adjoint = .false./test_adjoint = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "11: Test Unstruct" {
    skip "To do"
    sed "s/test_unstruct = .false./test_unstruct = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "12: Test Group" {
    sed "s/test_group = .false./test_group = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "13: Test Interface" {
    sed "s/test_interface = .false./test_interface = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "14: Test Check Parallel" {
    sed "s/check_parallel = .false./check_parallel = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

@test "15: Test Get Nbr" {
    sed "s/test_get_nbr = .false./test_get_nbr = .true./" input.nml_base > input.nml
    run mpirun -n 6 ./test_mpp_domains
    [ "$status" -eq 0 ]
}

