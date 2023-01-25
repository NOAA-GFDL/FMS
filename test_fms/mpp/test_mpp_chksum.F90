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
!> @author Ryan Mulhall
!> @email gfdl.climate.model.info@noaa.gov
!> @brief Test mpp_chksum routines
!> @description Tests mpp_chksum with 8 and 4 byte integer and real arrays with
!> single pe and distributed checksums
program test_mpp_chksum

  use fms

  implicit none

  integer :: test_num = 1
  integer :: nx = 96, ny = 96, npz = 63
  logical :: debug = .false.
  integer :: npes, root, pe, ierr

  namelist /test_mpp_chksum_nml/ test_num, nx, ny

  call mpp_init

  read(input_nml_file, test_mpp_chksum_nml, iostat=ierr)
  ierr = check_nml_error(ierr, 'test_mpp_chksum_nml')

  root = mpp_root_pe()
  pe = mpp_pe()
  npes = mpp_npes()

  if( test_num .eq. 1 ) then
    call test_chksum_simple
  else if( test_num .eq. 2 ) then
    call test_chksum_int
  else if( test_num .eq. 3) then
    call test_chksum_real
  else
    call mpp_error(FATAL, 'test_mpp_chksum: invalid test number given')
  endif

  call mpp_exit

  contains
  !> tests mixed precision checksums with ints
  subroutine test_chksum_int()
  integer                        :: out_unit
  integer(i8_kind), allocatable  :: data8(:), distData(:),temp(:)
  integer(i8_kind)               :: res4, res8, resDist
  integer(i4_kind), allocatable  :: data4(:)
  real, allocatable              :: rands(:)
  integer                        :: i, length

  out_unit = stdout()

  !> generate random arrays
  length = 1024
  allocate(rands(length), data8(length), data4(length), distData(length))
  call random_number(rands)
  do i = 1, length
    data8(i) = rands(i) * huge(data4(1))
    data4(i) = rands(i) * huge(data4(1))
    distData(i) = rands(i) * huge(distData(1))
  end do
  !>test mixed precision int checksums
  res4 = mpp_chksum(data4)
  res8 = mpp_chksum(data8)
  if(res4.NE.res8) then
    call mpp_error(FATAL, 'Test mpp_chksum_int: mixed precision checksums do not match')
  else
    call mpp_error(NOTE, 'Test mpp_chksum_int: mixed precision checksums match')
  endif
  !>test distributed int checksums
  call mpp_sync()
  call mpp_transmit( put_data=distData(1), plen=length, to_pe=ALL_PES, &
                     get_data=distData(1),glen=length, from_pe=root)
  call mpp_sync_self()
  allocate(temp(length/npes))
  temp = distData( pe*(length/npes)+1 : (pe+1)*(length/npes))!> distribute data for pelist
  resDist = mpp_chksum(distData(1:length), (/pe/))
  if(resDist.NE.mpp_chksum(temp)) then
    call mpp_error(FATAL, 'Test mpp_chksum_int: distributed checksums do not match')
  else
    call mpp_error(NOTE, 'Test mpp_chksum_int: distributed checksums match')
  endif
  deallocate(rands, data8, data4, distData, temp)

  end subroutine

  subroutine test_chksum_simple()
    integer :: n2, m, n = 1048576
    real, dimension(1024) :: a
    real, dimension(:), allocatable :: c

   if( modulo(n,npes).EQ.0 )then  !only set up for even division
     n2 = 1024
     a = 0.d0
     if( pe.EQ.root )call random_number(a(1:n2))

     call mpp_sync()
     call mpp_transmit( put_data=a(1), plen=n2, to_pe=ALL_PES, &
                        get_data=a(1), glen=n2, from_pe=root )
     call mpp_sync_self ()

     m= n2/npes

     allocate( c(m) )
     c = a(pe*m+1:pe*m+m)

     if( pe.EQ.root )then
        print *
        print *, '------------------ > Test mpp_chksum <------------------ '
        print *, 'This test shows that a whole array and a distributed array give identical checksums.'
     end if

     if ( mpp_chksum(a(1:n2),(/pe/)) .NE. mpp_chksum(c) ) then
       call mpp_error(FATAL, &
                     & 'Test mpp_chksum fails: a whole array and a distributed array did not give identical checksums')
     else
       print *, 'For pe=', pe, ' chksum(a(1:1024))=chksum(c(1:1024))='
     endif

   else
     call mpp_error(FATAL, 'Test mpp_chksum: cannot run this test since n cannot be evenly by npes')
   end if

  end subroutine test_chksum_simple

  !--- originally written by rusty
  !--- tests mixed precision real checksums
  subroutine test_chksum_real()
    real(kind=r4_kind), allocatable, dimension(:) :: u4_safe1d
    real(kind=r8_kind), allocatable, dimension(:) :: u8_safe1d
    real(kind=r4_kind), allocatable, dimension(:,:,:) :: u4, u4_safe
    real(kind=r8_kind), allocatable, dimension(:,:,:) :: u8, u8_safe
    integer :: nfile,fnum, num
    integer :: num_entry,ne,k,j,i,l
    integer :: isc, iec, jsc, jec, km
    integer, allocatable, dimension(:) :: isg, ieg, jsg, jeg
    integer :: size_bef, size_aft
    real(kind=r4_kind) :: temp4
    real(kind=r8_kind) :: temp8
    integer(kind=i4_kind)  :: mold4, cold4(1)
    integer(kind=i8_kind) :: mold8, cold8(1)
    integer(kind=i4_kind),  allocatable, dimension(:) :: sumi4
    integer(kind=i8_kind), allocatable, dimension(:) :: sumi8
    integer(kind=i4_kind),  allocatable, dimension(:,:,:) :: sumi43d
    integer(kind=i8_kind), allocatable, dimension(:,:,:) :: sumi83d
    character(len=4) :: fn
    type(domain2d) :: domain
    integer, dimension(2) :: layout
    integer(kind=i8_kind), allocatable :: check_root(:), check_pe(:)
    integer(kind=i8_kind)              :: check_num

    call mpp_set_current_pelist()

  !--- root pe reads in all data
    if (mpp_pe() == mpp_root_pe()) then
      allocate (u4_safe1d(nx*ny*npz))
      allocate (u8_safe1d(nx*ny*npz))
      allocate (u4_safe(nx,ny,npz))
      allocate (u8_safe(nx,ny,npz))
      allocate (isg(mpp_npes()))
      allocate (ieg(mpp_npes()))
      allocate (jsg(mpp_npes()))
      allocate (jeg(mpp_npes()))
      l = 0
      num = 0
    !--- set data incrementally
      do k = 1,npz
        do j = 1,ny
          do i = 1,nx
            l = l + 1
            temp4 = i + j * 1e-3_r4_kind + k * 1e-6_r4_kind
            temp8 = i + j * 1e-3_r8_kind + k * 1e-6_r8_kind
            u8_safe(i,j,k) = temp8
            u8_safe1D(l) = temp8
            u4_safe(i,j,k) = temp4
            u4_safe1D(l) = temp4
          enddo
        enddo
      enddo

    endif

    !--- split up data by pe
    isc = nx/npes * pe  + 1
    iec = nx/npes * (pe+1)
    jsc = ny/npes * pe + 1
    jec = ny/npes * (pe+1)
    if(iec == nx-1) iec = iec + 1
    if(jec == ny-1) jec = jec + 1
    if(debug) print *, 'pe', pe, 'bounds:', isc, iec, jsc, jec

    call mpp_sync()

    !--- all pes read in their data chunk
    !--- uses same value scheme as global
    allocate (u4(nx,ny,npz))
    allocate (u8(nx,ny,npz))
    do k = 1,npz
      do j = 1,ny
        do i = 1,nx
          if( (isc .le. i .and. i .ge. iec) .and. (jsc .le. j .and. j .ge. jec) ) then
            temp4 = i + j * 1e-3_r4_kind + k * 1e-6_r4_kind
            temp8 = i + j * 1e-3_r8_kind + k * 1e-6_r8_kind
            u4(i,j,k) = temp4
            u8(i,j,k) = temp8
          endif
        enddo
      enddo
    enddo

    !--- check output of r8 sums
    allocate(check_root(4))
    if (mpp_pe() == mpp_root_pe()) then
      !! unified checksums
      check_root(1) = mpp_chksum(u8_safe, (/0/) )
      check_root(2) = mpp_chksum(u8_safe1d, (/0/) )
      check_root(3) = mpp_chksum(u8_safe1d(nx*ny*npz:1:-1), (/0/))
      if( check_root(1) .ne. check_root(2) .or. check_root(2) .ne. check_root(3) ) &
          call mpp_error(FATAL,'test_mpp_chksum: r8 single pe checksums do not match')
      !! distributed
      check_root(4) = mpp_chksum(u8)
      call mpp_broadcast( check_root, 4, root)
      if(debug) then
        print *, 'KIND=8 chksums'
        print *, 'unified chksum is      :', check_root(1)
        print *, 'unified 1D chksum is   :', check_root(2)
        print *, 'unified 1D-R chksum is :', check_root(3)
        print *, 'distributed checksum is:', check_root(4)
      endif
    else
      check_num = mpp_chksum(u8)
      call mpp_broadcast( check_root, 4, root)
      if(check_num .ne. check_root(4)) &
        call mpp_error(FATAL, 'test_mpp_chksum: r8 distributed does not match root')
    endif

    call mpp_sync()
    check_root = 0

    !--- check output of r4 sums
    if (mpp_pe() == mpp_root_pe()) then
      !! unified checksums
      check_root(1) = mpp_chksum(u4_safe, (/0/) )
      check_root(2) = mpp_chksum(u4_safe1d, (/0/) )
      check_root(3) = mpp_chksum(u4_safe1d(nx*ny*npz:1:-1), (/0/))
      if( check_root(1) .ne. check_root(2) .or. check_root(2) .ne. check_root(3) ) &
          call mpp_error(FATAL,'test_mpp_chksum: r8 single pe checksums do not match')
      !! distributed
      check_root(4) = mpp_chksum(u4)
      call mpp_broadcast( check_root, 4, root)
      if(debug) then
        print *, 'KIND=4 chksums'
        print *, 'unified chksum is      :', check_root(1)
        print *, 'unified 1D chksum is   :', check_root(2)
        print *, 'unified 1D-R chksum is :', check_root(3)
        print *, 'distributed checksum is:', check_root(4)
      endif
    else
      check_num = mpp_chksum(u4)
      call mpp_broadcast( check_root, 4, root)
      if(debug) print *, pe, check_root, check_num
      if(check_root(4) .ne. check_num) call mpp_error(FATAL, 'test_mpp_chksum: r4 distributed does not match root')
    endif

    !! integer
    call mpp_sync()

    km = npz
    allocate (sumi4(nx * ny * npz))
    allocate (sumi8(nx * ny * npz))
    allocate (sumi43d(isc:iec,jsc:jec,1:km))
    allocate (sumi83d(isc:iec,jsc:jec,1:km))
    sumi4 = transfer(u4,cold4)
    size_bef = size(sumi8)
    sumi8 = transfer(u4,cold8)
    size_aft = size(sumi8)
    do k = 1,km
     do j = jsc,jec
      do i = isc,iec
       sumi43d(i,j,k) = transfer(u4(i,j,k),mold4)
       sumi83d(i,j,k) = transfer(u4(i,j,k),mold8)
      enddo
     enddo
    enddo

    deallocate(check_root)
    allocate(check_root(6))
    allocate(check_pe(6))

    if (mpp_pe() == mpp_root_pe()) then

      check_root(1) = mpp_chksum(int(sumi4,kind=i8_kind))
      call mpp_sync()
      check_root(2) = mpp_chksum(int(sumi4(nx*ny*npz:1:-1),kind=i8_kind))
      call mpp_sync()
      check_root(3) = mpp_chksum(sumi8)
      call mpp_sync()
      check_root(4) = mpp_chksum(int(sumi43d,kind=i8_kind))
      call mpp_sync()
      check_root(5) = mpp_chksum(int(sumi43d(iec:isc:-1,jec:jsc:-1,km:1:-1),kind=i8_kind))
      call mpp_sync()
      check_root(6) = mpp_chksum(sumi83d)
      call mpp_broadcast(check_root, 6, root)
      if(debug) then
        print *, 'INTEGER chksums'
        print *, 'distributed i4   array chksum:', check_root(1)
        print *, 'distributed i4-R array chksum:', check_root(2)
        print *, 'distributed i8   array chksum:', check_root(3)
        print *, 'distributed i4   value chksum:', check_root(4)
        print *, 'distributed i4-R value chksum:', check_root(5)
        print *, 'distributed i8   value chksum:', check_root(6)
      endif
    else
      check_pe(1) = mpp_chksum(int(sumi4,kind=i8_kind))
      call mpp_sync()
      check_pe(2) = mpp_chksum(int(sumi4(nx*ny*npz:1:-1),kind=i8_kind))
      call mpp_sync()
      check_pe(3) = mpp_chksum(sumi8)
      call mpp_sync()
      check_pe(4) = mpp_chksum(int(sumi43d,kind=i8_kind))
      call mpp_sync()
      check_pe(5) = mpp_chksum(int(sumi43d(iec:isc:-1,jec:jsc:-1,km:1:-1),kind=i8_kind))
      call mpp_sync()
      check_pe(6) = mpp_chksum(sumi83d)

      call mpp_broadcast(check_root, 6, root)
      if(debug) print *, 'pe:', pe, check_pe

      if( any(check_root .ne. check_pe) ) call mpp_error(FATAL, 'test_mpp_chksum: int checksums do not match')
    endif

  end subroutine
end program
