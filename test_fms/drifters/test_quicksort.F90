      program test_quicksort
        implicit none
        integer :: list(16) = (/6, 2, 3, 4, 1, 45, 3432, 3245, 32545, 66555, 32, 1,3, -43254, 324, 54/)
        print *,'before list=', list
        call qksrt_quicksort(size(list), list, 1, size(list))
        print *,'after  list=', list
      end program test_quicksort
