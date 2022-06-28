program test_output_yaml
#ifdef use_yaml

use fms_yaml_output_mod
use fms_string_utils_mod
implicit none

!> This program will print out the following delicious yaml:\
!! \verbatim
!! ---
!! name: time to eat
!! location: Bridgewater, NJ
!! order:
!! - Drink: Iced tea
!!   Food:
!!   - Main: pancake
!!     side: eggs
!!     sauce: hot
!!   - Appetizer: wings
!!     dip: ranch
!! - Drink: milk
!!   paper: coloring
!!   crayon: purple
!!   fork: plastic
!!   spoon: silver
!!   knife: none
!!   Food:
!!   - Main: cereal
!!     sauce: milk
!! - Drink: coffee
!!   fork: silver
!!   knife: steak
!!   Meal:
!!   - app: poppers
!!     sauce: tangy
!!   - main: steak
!!     side: mashed
!!     sauce: A1
!!   - dessert: cake
!!     topping: frosting
!! ...
!! \end verbatim
!! Great, now I have to create this long yaml for testing, lol.
type (fmsYamlOutKeys_type), allocatable :: k1 (:)
type (fmsYamlOutValues_type), allocatable :: v1 (:)
type (fmsYamlOutKeys_type), allocatable :: k2 (:)
type (fmsYamlOutValues_type), allocatable :: v2 (:)
type (fmsYamlOutKeys_type), allocatable :: k3 (:)
type (fmsYamlOutValues_type), allocatable :: v3 (:)
integer :: a1size = 1
integer :: a2 = 3
integer :: a3 
integer,allocatable :: a3each (:)
character(len=:), allocatable :: yaml_reference
character(len=string_len_parameter) :: tmpstr
integer :: i !< for looping
!> Set the number of "third level" elements and calculate a3
allocate (a3each(a2))
a3each(1) = 2
a3each(2) = 1
a3each(3) = 3
a3 = sum(a3each)
!> allocate all of the arrays
allocate(k1(a1size))
allocate(v1(a1size))
allocate(k2(a2))
allocate(v2(a2))
allocate(k3(a3))
allocate(v3(a3))

!> Copy the strings into the key/value pairings
call fms_f2c_string (k1(1)%key1,"name")
call fms_f2c_string (v1(1)%val1,"time to eat")
call fms_f2c_string (k1(1)%key2,"location")
call fms_f2c_string (v1(1)%val2,"Bridgewater, NJ")
call fms_f2c_string (k1(1)%level2key,"order")

call fms_f2c_string (k2(1)%key1,"Drink")
call fms_f2c_string (v2(1)%val1, "Iced tea")
call fms_f2c_string (k2(1)%level2key,"Food")
call fms_f2c_string (k3(1)%key1,"Main")
call fms_f2c_string (v3(1)%val1,"pancake")
call fms_f2c_string (k3(1)%key7,"side")
call fms_f2c_string (v3(1)%val7,"eggs")
call fms_f2c_string (k3(1)%key8,"sauce")
call fms_f2c_string (v3(1)%val8,"hot")
call fms_f2c_string (k3(2)%key1,"Appetizer")
call fms_f2c_string (v3(2)%val1,"wings")
call fms_f2c_string (k3(2)%key7,"dip")
call fms_f2c_string (v3(2)%val7,"ranch")

call fms_f2c_string (k2(2)%key1,"Drink")
call fms_f2c_string (v2(2)%val1, "Milk")
call fms_f2c_string (k2(2)%key2,"paper")
call fms_f2c_string (v2(2)%val2, "coloring")
call fms_f2c_string (k2(2)%key3,"crayon")
call fms_f2c_string (v2(2)%val3, "purple")
call fms_f2c_string (k2(2)%key4,"fork")
call fms_f2c_string (v2(2)%val4, "plastic")
call fms_f2c_string (k2(2)%key5,"spoon")
call fms_f2c_string (v2(2)%val5, "silver")
call fms_f2c_string (k2(2)%key12,"knife")
call fms_f2c_string (v2(2)%val12, "none")
call fms_f2c_string (k2(2)%level2key,"Food")
call fms_f2c_string (k3(3)%key1,"Main")
call fms_f2c_string (v3(3)%val1,"cereal")
call fms_f2c_string (k3(3)%key7,"sauce")
call fms_f2c_string (v3(3)%val7,"milk")

call fms_f2c_string (k2(3)%key1,"Drink")
call fms_f2c_string (v2(3)%val1, "coffee")
call fms_f2c_string (k2(3)%key2,"fork")
call fms_f2c_string (v2(3)%val2, "silver")
call fms_f2c_string (k2(3)%key13,"knife")
call fms_f2c_string (v2(3)%val13, "steak")
call fms_f2c_string (k2(3)%level2key,"Meal")
call fms_f2c_string (k3(4)%key1,"app")
call fms_f2c_string (v3(4)%val1,"poppers")
call fms_f2c_string (k3(4)%key7,"sauce")
call fms_f2c_string (v3(4)%val7,"tangy")
call fms_f2c_string (k3(5)%key4,"main")
call fms_f2c_string (v3(5)%val4,"steak")
call fms_f2c_string (k3(5)%key7,"side")
call fms_f2c_string (v3(5)%val7,"mashed")
call fms_f2c_string (k3(5)%key11,"sauce")
call fms_f2c_string (v3(5)%val11,"A1")
call fms_f2c_string (k3(6)%key10,"dessert")
call fms_f2c_string (v3(6)%val10,"cake")
call fms_f2c_string (k3(6)%key11,"topping")
call fms_f2c_string (v3(6)%val11,"frosting")
!> Write the yaml
call write_yaml_from_struct_3 (1, k1, v1, a2, k2, v2, 0, a3each, k3, v3)
#endif
end program test_output_yaml

