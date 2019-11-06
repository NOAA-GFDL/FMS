setup () {
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
    skip_test=true
  elif [ $nProc -lt 4 ]
  then
    # Need to oversubscribe the MPI
    #
    # Is the oversubscribe option known?
    mpirun -oversubscribe -version 2>&1 > /dev/null
    if [ $? -eq 0 ]
    then
      # Looks like open MPI mpirun
      oversubscribe="-oversubscribe"
    fi
  fi
}

teardown () {
  # Echo the output.  Will only be done on a test failure.
  echo "$output"
}

@test "1" {
   if [ "$skip_test" = "true" ]
   then
     skip "Unable to determine number of available cpus"
   fi
   # Ensure an input.nml file exists
   touch input.nml
   run mpirun ${oversubscribe} -n 4 ./test_mpp
   [ "$status" -eq 0 ]
}
