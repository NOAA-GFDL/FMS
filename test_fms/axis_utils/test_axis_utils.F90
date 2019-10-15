!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

program test_axis_utils

use fms_mod,       only : fms_init, fms_end, file_exist, open_namelist_file, check_nml_error
use fms_mod,       only : close_file
use mpp_mod,       only : mpp_error, FATAL, stdout
use mpp_mod,       only : input_nml_file
use axis_utils_mod, only: interp_1d

implicit none



integer, parameter :: maxsize = 100

integer :: n_src = 0
integer :: n_dst = 0
real, dimension(MAXSIZE) :: grid_src = 0
real, dimension(MAXSIZE) :: grid_dst = 0
real, dimension(MAXSIZE) :: data_src = 0
real, dimension(MAXSIZE) :: out_linear_dst = 0
real, dimension(MAXSIZE) :: out_cubic_dst = 0
real :: diff

namelist / test_axis_utils_nml / n_src, n_dst, grid_src, grid_dst, data_src, out_linear_dst, out_cubic_dst

real, allocatable :: data_dst(:)
integer           :: unit, ierr, io

  call fms_init()

  !--- default option of data

 n_src = 31
  n_dst = 40
  grid_src(1:n_src) = (/ -63.6711465476916, -63.6711455476916, 166.564180735096, 401.25299580552, &
                         641.056493022762, 886.219516665347, 1137.35352761133, 1394.4936854079,   &
                         1657.17893448689, 1925.64572676068, 2200.13183483549, 2480.9124139255,   &
                         2768.35396680912, 3062.86513953019, 3675.47369643284, 4325.10564183322,  &
                         5020.19039479527, 5769.85432323481, 6584.25101514851, 7475.94655633703,  &
                         8462.01951335773, 9568.28246037887, 10178.3869413515, 10834.1425668942,  &
                         11543.5265942777, 12317.3907407535, 13170.4562394288, 14125.6466646843,  &
                         15225.8720618086, 16554.7859690842, 19697.1334102613   /)
  grid_dst(1:n_dst) = (/ 1002.9522552602, 1077.51144617887, 1163.37842788755, 1264.19848463606,  &
                         1382.57557953916, 1521.56713587855, 1684.76300370633, 1876.37817787584, &
                         2101.36166220498, 2365.52429149707, 2675.68881278444, 3039.86610206727, &
                         3467.4620678435, 3969.52058529847, 4553.81573511231, 5159.54844211827,  &
                         5765.28114912423, 6371.01385613019, 6976.74656313614, 7582.4792701421,  &
                         8188.21197714806, 8793.94468415402, 9399.67739115997, 10005.4100981659, &
                         10611.1428051719, 11216.8755121778, 11822.6082191838, 12428.3409261898, &
                         13034.0736331957, 13639.8063402017, 14245.5390472076, 14851.2717542136, &
                         15457.0044612196, 16062.7371682255, 16668.4698752315, 17274.2025822374, &
                         17879.9352892434, 18485.6679962493, 19091.4007032553, 19697.1334102613 /)
  data_src(1:n_src) = (/ 309.895999643929, 309.991081541887, 309.971074746584, 310.873654697145, &
                         311.946530606618, 312.862249229647, 314.821236806913, 315.001269608758, &
                         315.092410930288, 315.19010999336,  315.122964496815, 315.057882573487, &
                         314.998796850493, 314.984586411292, 315.782246062002, 318.142544345795, &
                         321.553905292867, 325.247730854554, 329.151282227113, 332.835673638378, &
                         336.810414210932, 341.64530983048,  344.155248759994, 346.650476976385, &
                         349.106430095269, 351.915323032738, 354.709396583792, 359.68904432446,  &
                         371.054289820675, 395.098187506342, 446.150726850039 /)
  out_linear_dst(1:n_src)  = (/ 313.772830731158,  314.354434665370,  314.839457748187, 314.910045389784, &
                         314.992925326750, 315.045359036116, 315.102449183793, 315.172180798747, &
                         315.147125910052, 315.084628301276, 315.017844853703, 314.985696136334, &
                         315.511400215769, 316.850602339827, 319.265015637373, 322.240565405352, &
                         325.225197414186, 328.129197765043, 330.773032206975, 333.265094095404, &
                         335.706729219375, 338.261085191662, 340.908425428092, 343.443630789374, &
                         345.801936337123, 347.975533834707, 350.119412036280, 352.278721832823, &
                         354.262698137434, 357.156236613905, 360.927523607037, 367.184696853722, &
                         375.236143688601, 386.195600996087, 396.945167256346, 406.786279175083, &
                         416.627391093823, 426.468503012560, 436.309614931300, 446.150726850039 /)
 out_cubic_dst(1:n_src) = (/ -313.942318474633, 314.503163655656, 314.913470956341, 315.055635714636, &
                         315.006673321389, 315.010907002311, 315.109298816167, 315.185728881530, &
                         315.155647446989, 315.081117180934, 315.017858886237, 314.979618544023, &
                         315.354523043308, 316.685513873638, 319.218943479033, 322.249746569912, &
                         325.225277554867, 328.168923274850, 330.835994100920, 333.255253114313, &
                         335.672089973660, 338.246525654852, 340.917577732601, 343.460878949881, &
                         345.833971678071, 347.972895520732, 350.130774046054, 352.280193191086, &
                         354.219167222638, 356.767905094900, 360.561704156599, 366.216858589885, &
                         374.647838146611, 385.613883532528, 397.240573301678, 408.134997959681, &
                         418.287826652131, 427.884458370414, 437.110292105921, 446.150726850039 /)

  !---reading namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, test_axis_utils_nml, iostat=io)
      ierr = check_nml_error(io,'test_axis_utils_nml')
#else
  if(file_exist('input.nml')) then
    unit =  open_namelist_file()
       ierr=1
    do while (ierr /= 0)
          read  (unit, nml=test_axis_utils_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'test_axis_utils_nml')  ! also initializes nml error codes
    enddo
 10    call close_file(unit)
  endif
#endif

  if(n_src >MAXSIZE) call mpp_error(FATAL, 'test_axis_utils: nml n_src is greater than MAXSIZE')
  if(n_dst >MAXSIZE) call mpp_error(FATAL, 'test_axis_utils: nml n_dst is greater than MAXSIZE')

  allocate(data_dst(n_dst) )


  !--- write out data
  unit = stdout()
  write(unit,*)' the source grid is ', grid_src(1:n_src)
  write(unit,*)' the destination grid is ', grid_dst(1:n_dst)
  write(unit,*)' the source data is ', data_src(1:n_src)

  !--- testing linear interpolation
  call interp_1d(grid_src(1:n_src), grid_dst(1:n_dst), data_src(1:n_src), data_dst, "linear")
  write(unit,*)' the destination data using linear interpolation is ', data_dst(1:n_dst)
  diff = sum(abs(data_dst - out_linear_dst(1:n_dst)))
  write(unit,*)' the total difference between the result and the expected result is ', diff
  if(diff > 1.0e-8) call mpp_error(FATAL, 'test_axis_utils: the result with linear interpolation is different')

  !--- testing cubic spline interpolation
  call interp_1d(grid_src(1:n_src), grid_dst(1:n_dst), data_src(1:n_src), data_dst, "cubic_spline")
  write(unit,*)' the destination data using cublic spline interpolation is ', data_dst(1:n_dst)
  diff = sum(abs(data_dst - out_cubic_dst(1:n_dst)))
  write(unit,*)' the total difference between the result and the expected result is ', diff
  if(diff > 1.0e-8) call mpp_error(FATAL, 'test_axis_utils: the result with cubic spline interpolation is different')

   call fms_end()
end program test_axis_utils
