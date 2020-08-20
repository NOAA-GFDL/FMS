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
!> @author Jessica Liptak
!> @brief Compare the checksums of 2D and 3D 32-bit or 64-bit integer arrays
module compare_data_checksums_int

use mpp_mod, only : mpp_root_pe, mpp_chksum, mpp_error, mpp_sync_self, mpp_pe
use mpp_mod, only : FATAL, NOTE
use platform_mod, only : i4_kind, i8_kind

implicit none
private


integer :: stdunit = 6

public :: compare_checksums_int

interface compare_checksums_int
  module procedure compare_checksums_2D_i4
  module procedure compare_checksums_3D_i4
  module procedure compare_checksums_2D_i8
  module procedure compare_checksums_3D_i8
end interface compare_checksums_int

contains

  !> compare checksums of 2D 32-bit integer arrays
  subroutine compare_checksums_2D_i4( a, b, chk_str )
    integer(kind=i4_kind), intent(in), dimension(:,:) :: a, b !< 2D arrays to compare
    character(len=*), intent(in) :: chk_str
    integer(kind=i8_kind) :: sum1, sum2
    integer :: i, j
    integer :: pe
    !> @note can't call mpp_sync here since there might be different number of tiles on each pe.
    call mpp_sync_self()
    pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
      call mpp_error(FATAL,'compare_checksums_2D_r4: sizes of a and b do not match')

    do j = 1, size(a,2)
      do i = 1, size(a,1)
        if(a(i,j) .ne. b(i,j)) then
          print*, "a =", a(i,j)
          print*, "b =", b(i,j)
          write(*,'(a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')"at the pe ", mpp_pe(), &
                ", at point (",i,", ", j, "),a=", a(i,j), ",b=", b(i,j)
          call mpp_error(FATAL, trim(chk_str)//': value mismatch at data point.')
        endif
      enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(chk_str)//': OK.' )
      !> @note in some cases, even though the checksum agree, the two arrays
      !! actually are different [e.g.,(1.1,-1.2) with (-1.1,1.2)].
      !! Thus, we need to check the values point-by-point.
    else
      call mpp_error( FATAL, trim(chk_str)//': checksums do not match.' )
    end if
  end subroutine compare_checksums_2D_i4

  !> Compare the checksums of 2 3D 32-bit real arrays
  subroutine compare_checksums_3D_i4( a, b, string )
     integer(kind=i4_kind), intent(in), dimension(:,:,:) :: a, b !< 3D 64-bit real arrays to compare
     character(len=*), intent(in) :: string
     integer(kind=i8_kind) :: sum1, sum2
     integer :: i, j, k
     integer :: pe
    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
     call mpp_sync_self()
     pe = mpp_pe()

     if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
       call mpp_error(FATAL,'compare_checkums_3d_r4: sizes of a and b do not match')

     do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                write(*,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)') trim(string)//" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': mismatch in checksums at data point.')
             endif
          enddo
       enddo
     enddo

     sum1 = mpp_chksum( a, (/pe/) )
     sum2 = mpp_chksum( b, (/pe/) )

     if( sum1.EQ.sum2 )then
       if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
     else
       write(stdunit, *)"sum1 =", sum1, mpp_pe()
       write(stdunit, *)"sum2 =", sum2, mpp_pe()
       write(stdunit,'(a,i3,a,i20,a,i20)')" at pe ", mpp_pe(), " sum(a)=", sum1, " sum(b)=", sum2
       call mpp_error( FATAL, trim(string)//': checksums do not match.' )
     end if
  end subroutine compare_checksums_3D_i4

  !> compare checksums of 2D 64-bit integer arrays
  subroutine compare_checksums_2D_i8( a, b, chk_str )
    integer(kind=i8_kind), intent(in), dimension(:,:) :: a, b !< 2D arrays to compare
    character(len=*), intent(in) :: chk_str
    integer(kind=i8_kind) :: sum1, sum2
    integer :: i, j
    integer :: pe

    !> @note can't call mpp_sync here since there might be different number of tiles on each pe.
    call mpp_sync_self()
    pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
      call mpp_error(FATAL,'compare_checksums_2d_r8: sizes of a and b do not match')

    do j = 1, size(a,2)
      do i = 1, size(a,1)
        if(a(i,j) .ne. b(i,j)) then
          print*, "a =", a(i,j)
          print*, "b =", b(i,j)
          write(*,'(a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')"at the pe ", mpp_pe(), &
                ", at point (",i,", ", j, "),a=", a(i,j), ",b=", b(i,j)
          call mpp_error(FATAL, trim(chk_str)//': value mismatch at data point.')
        endif
      enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(chk_str)//': OK.' )
      !> @note in some cases, even though the checksum agree, the two arrays
      !! actually are different [e.g.,(1.1,-1.2) with (-1.1,1.2)].
      !! Thus, we need to check the values point-by-point.
    else
      call mpp_error( FATAL, trim(chk_str)//': checksums do not match.' )
    end if
  end subroutine compare_checksums_2D_i8

  !> Compare the checksums of 2 3D 64-bit real arrays
  subroutine compare_checksums_3D_i8( a, b, string )
     integer(kind=i8_kind), intent(in), dimension(:,:,:) :: a, b !< 3D 64-bit real arrays to compare
     character(len=*), intent(in) :: string
     integer(kind=i8_kind) :: sum1, sum2
     integer :: i, j, k
     integer :: pe
    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
     call mpp_sync_self()
     pe = mpp_pe()

     if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
       call mpp_error(FATAL,'compare_checksums_3d_r8: size of a and b does not match')

     do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                write(*,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)') trim(string)//" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': mismatch in checksums at data point.')
             endif
          enddo
       enddo
     enddo

     sum1 = mpp_chksum( a, (/pe/) )
     sum2 = mpp_chksum( b, (/pe/) )

     if( sum1.EQ.sum2 )then
       if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
     else
       write(stdunit, *)"sum1 =", sum1, mpp_pe()
       write(stdunit, *)"sum2 =", sum2, mpp_pe()
       write(stdunit,'(a,i3,a,i20,a,i20)')" at pe ", mpp_pe(), " sum(a)=", sum1, " sum(b)=", sum2
       call mpp_error( FATAL, trim(string)//': checksums do not match.' )
     end if
  end subroutine compare_checksums_3D_i8

end module compare_data_checksums_int
