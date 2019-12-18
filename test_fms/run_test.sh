# Function to set up and run tests. 
# Inputs: 
# {1} Name of test 
# {2} Number of processors 
# {3} Set to skip, if the test is meant to fail 
# {4} Set to true if you want mpi launcher to oversubscribe, IF possible
# For example: run_test test_time_manager 1

run_test()
{
    # If the tests is known to fail exit 
    if test "x${3}" == "xskip" ; then 
        exit 0
    fi     

    # If there is no mpi launcher just ./{job_script}
    if test "x$mpi_launcher" != "x" ; then  
        npes="-n ${2}"
    fi

    # Check if the oversubscribe flag is turned on 
    if test "x${4}" == "true" ; then 
    # Check if the your mpi launcher allows the oversubscribed option 
        if "x$oversubscribed" != "x"; then  
            $mpi_launcher $oversubscribed $npes ./${1}
            exit 1
        else #If you the mpi launcher doesn't allow the oversubscribed option, don't run the test
            exit 0
        fi
    fi 

    # Run test
    $mpi_launcher $npes ./${1}

}
