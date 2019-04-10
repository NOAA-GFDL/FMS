program test_mosaic

use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_ncontacts
use mosaic_mod, only : get_mosaic_grid_sizes, get_mosaic_contact

implicit none

integer              :: ntiles, ncontacts, n
integer, allocatable :: tile1(:), tile2(:), nx(:), ny(:)
integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:)
integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:)
character(len=128)   :: mosaic_file = "INPUT/mosaic.nc"

ntiles = get_mosaic_ntiles(mosaic_file)
ncontacts = get_mosaic_ncontacts(mosaic_file)
allocate(nx(ntiles), ny(ntiles))
allocate(tile1(ncontacts), tile2(ncontacts) )
allocate(istart1(ncontacts), iend1(ncontacts), jstart1(ncontacts), jend1(ncontacts) )
allocate(istart2(ncontacts), iend2(ncontacts), jstart2(ncontacts), jend2(ncontacts) )

call get_mosaic_grid_sizes(mosaic_file, nx, ny )
call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2)

! print out information

print '(a,i3,a,a)', "****** There is ", ntiles, " tiles in ", trim(mosaic_file)
do n = 1, ntiles
   print '(a,i3,a,i3,a,i3)', " tile = ", n, ", nx = ", nx(n), ", ny = ", ny(n)
end do

print '(a,i3,a,a)', "****** There is ", ncontacts, " contacts in ", trim(mosaic_file)
do n = 1, ncontacts
   print '(a,i3,a,i3,a,i3,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
           "contact=", n, ": tile1=", tile1(n), " tile2=", tile2(n),   &
           " is1=", istart1(n), " ie1=", iend1(n),                   &
           " js1=", jstart1(n), " je1=", jend1(n),                   &
           " is2=", istart2(n), " ie2=", iend2(n),                   &
           " js2=", jstart2(n), " je2=", jend2(n)
end do

deallocate(tile1, tile2, nx, ny)
deallocate(istart1, iend1, jstart1, jend1)
deallocate(istart2, iend2, jstart2, jend2)


end program test_mosaic
