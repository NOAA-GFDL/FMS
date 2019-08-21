@test "1" {
    skip "The input file is missing"
    run ./test_drifters
    [ "$status" -eq 0 ]
}
