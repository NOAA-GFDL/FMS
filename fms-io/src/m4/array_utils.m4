include(`macros.m4')dnl
!> @brief Allocate arrays using an input array of sizes.
`subroutine allocate_array_'KIND`_'NUM_DIMS`d(buf, &'
                                              sizes)
    TYPE`(kind='KIND`),'dim_declare(NUM_DIMS) allocatable,intent(inout) :: buf !< Array that will be
                                                                               !! allocated.
    integer,dimension(NUM_DIMS),intent(in) :: sizes !< Array of dimension sizes.
    if (allocated(buf)) then
        deallocate(buf)
    endif
    allocate(buf(dim_sizes(1,NUM_DIMS,sizes)))
`end subroutine allocate_array_'KIND`_'NUM_DIMS`d'

!> @brief Put a section of an array into a larger array.
`subroutine put_array_section_'KIND`_'NUM_DIMS`d(section, &'
                                                 array, &
                                                 s, &
                                                 c)
    TYPE`(kind='KIND`),'dim_declare(NUM_DIMS) intent(in) :: section
    TYPE`(kind='KIND`),'dim_declare(NUM_DIMS) intent(inout) :: array
    integer,dimension(NUM_DIMS),intent(in) :: s !< Array of starting indices.
    integer,dimension(NUM_DIMS),intent(in) :: c !< Array of sizes.
    array(dim_slices(1,NUM_DIMS,s,c)) = section(dim_colons(1,NUM_DIMS))
`end subroutine put_array_section_'KIND`_'NUM_DIMS`d'


!> @brief Get a section of larger array.
`subroutine get_array_section_'KIND`_'NUM_DIMS`d(section, &'
                                                 array, &
                                                 s, &
                                                 c)
    TYPE`(kind='KIND`),'dim_declare(NUM_DIMS) intent(inout) :: section
    TYPE`(kind='KIND`),'dim_declare(NUM_DIMS) intent(in) :: array
    integer,dimension(NUM_DIMS),intent(in) :: s !< Array of starting indices.
    integer,dimension(NUM_DIMS),intent(in) :: c !< Array of sizes.
    section(dim_colons(1,NUM_DIMS)) = array(dim_slices(1,NUM_DIMS,s,c))
`end subroutine get_array_section_'KIND`_'NUM_DIMS`d'
