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
!> @brief Test mpp_chksum with mixed precision integers
!> @description Tests mpp_chksum with 8 and 4 byte integer arrays with
!> normal and distributed checksums
program test_mpp_chksum

  use fms

  implicit none

  integer :: test_num = 1
  !integer :: num_files = 36
  !character(len=32) :: dir = "./"
  character(len=32) :: safefile = "test_mpp_chksum_output.txt.out"
  integer :: nx = 96, ny = 96, npz = 32
!  integer :: nx = 32, ny = 32, npz = 16
  logical :: debug = .true. !< turns on debug output, writes to safefile
  integer :: npes, root, pe, ierr

  namelist /test_mpp_chksum_nml/ test_num !, nx, ny, npz, debug, safefile
  !namelist /base_nml/ num_files, dir, safefile, debug, nx, ny, npz
  !namelist /second_nml/ num_files, dir, safefile, debug, nx, ny, npz

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

  !call mpp_exit
  call MPI_FINALIZE

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
    integer, allocatable :: check_buff(:)
    integer, dimension(4) :: pelist

    call mpp_set_current_pelist()
    call mpp_domains_init()
    call mpp_domains_set_stack_size(3145746)

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

    !do nfile = 1,num_files
    !  fnum = 999 + nfile
    !  write(fn,'(i4)')  fnum
    !  if (debug) write(6,*) trim(dir)//'/fort.'//fn
    !  open(20,file=trim(dir)//'/fort.'//fn)
    !  read(20,*) isg(nfile),ieg(nfile),jsg(nfile),jeg(nfile),km
    !  num_entry = (ieg(nfile)-isg(nfile)+1)*(jeg(nfile)-jsg(nfile)+1)*km
    !  do ne = 1,num_entry
    !    l = l + 1
    !    read(20,*) i,j,k,temp8
    !    u8_safe(i,j,k) = temp8
    !    u8_safe1d(l) = temp8
    !  enddo
    !  rewind(20)
    !  read(20,*)
    !  do ne = 1,num_entry
    !    m = m + 1
    !    read(20,*) i,j,k,temp4
    !    u4_safe(i,j,k) = temp4
    !    u4_safe1d(m) = temp4
    !  enddo
    !  close(20)
    !enddo
      if (debug) then
        open (10,file=trim(safefile))
        do k = 1,npz
         do j = 1,ny
          do i = 1,nx
           write(10,'(3(i5,2x))') i,j,k
           write(10,*) u4_safe(i,j,k)
          enddo
         enddo
        enddo
        close(10)
      endif
    endif ! end root section

    call mpp_define_layout( (/ 1, nx, 1, ny /), mpp_npes(), layout)

    call mpp_define_domains( (/ 1, nx, 1, ny /), layout, domain)
    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)

    if (mpp_pe() == mpp_root_pe()) write(6,*) 'Root-pe here at sync'

    !call mpp_sync()

    !--- all pes read in their data chunk

    !fnum = 1000 + mpp_pe()
    !write(fn,'(i4)')  fnum
    !if (debug) write(6,*) trim(dir)//'/fort.'//fn
    !open(20,file=trim(dir)//'/fort.'//fn)
    !read(20,*) isc,iec,jsc,jec,km
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

    !num_entry = (iec-isc+1)*(jec-jsc+1)*km
    !do ne = 1,num_entry
    !  read(20,*) i,j,k,temp8
    !  u8(i,j,k) = temp8
    !enddo
    !rewind(20)
    !read(20,*)
    !do ne = 1,num_entry
    !  read(20,*) i,j,k,temp4
    !  u4(i,j,k) = temp4
    !enddo
    !close(20)
    !allocate(check_buff(0:npes-1))

    !! r8 sums
    if (mpp_pe() == mpp_root_pe()) then
      !check_buff(pe) = mpp_chksum(u8)
      !print *, mpp_chksum(u8)
      !call mpp_broadcast( check_buff, npes, root)
!!      print *, u8_safe
      if( debug) then
        print *, 'KIND=8 chksums'
        print *, 'unified chksum is      :', mpp_chksum(u8_safe,   (/0/))
        print *, 'unified 1D chksum is   :', mpp_chksum(u8_safe1d, (/0/))
        print *, 'unified 1D-R chksum is :', mpp_chksum(u8_safe1d(nx*ny*npz:1:-1), (/0/))
      endif

      !!- check root sums match then check pe sums

      !call mpp_transmit( mpp_chksum(u8), 1, pe, root_rbuf, 1, ALL_PES)
      !call mpp_broadcast()
    else
      !!- send pe sums

       ! print *, 'distributed checksum is:',mpp_chksum(u8)
      !call mpp_broadcast( check_buff, npes, root)
      !check_buff(pe) = mpp_chksum(u8)
      !open(10,file="/dev/null")
      !write(10,*) mpp_chksum(u8)
      !close(10)
      !call mpp_transmit( mpp_chksum(u8), 1, pe, pe_rbuf, 1, root)
    endif

    !!call mpp_sync()
    !! r4 sums
    if (mpp_pe() == mpp_root_pe()) then
      if(debug) then
        print *, 'KIND=4 chksums'
        print *, 'unified chksum is      :', mpp_chksum(u4_safe,   (/0/))
        print *, 'unified 1D chksum is   :', mpp_chksum(u4_safe1d, (/0/))
        print *, 'unified 1D-R chksum is :', mpp_chksum(u4_safe1d(nx*ny*npz:1:-1), (/0/))
       ! print *, 'distributed checksum is:', mpp_chksum(u4)
      endif
      !! - check sums match and pe sums
    else
      !!- send pe sums

      !open(10,file="/dev/null")
      !write(10,*) mpp_chksum(u4)
      !print *, pe, mpp_chksum(u4)
      !close(10)
    endif

    !return
    !call mpp_sync()
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

    if (mpp_pe() == mpp_root_pe()) then
      print *, 'INTEGER chksums'
      print *, 'distributed i4   array chksum:', mpp_chksum(int(sumi4,kind=i8_kind))
      call mpp_sync()
      print *, 'distributed i4-R array chksum:', mpp_chksum(int(sumi4(nx*ny*npz:1:-1),kind=i8_kind))
      call mpp_sync()
      print *, 'distributed i8   array chksum:', mpp_chksum(sumi8)
      call mpp_sync()
      print *, 'distributed i4   value chksum:', mpp_chksum(int(sumi43d,kind=i8_kind))
      call mpp_sync()
      print *, 'distributed i4-R value chksum:', mpp_chksum(int(sumi43d(iec:isc:-1,jec:jsc:-1,km:1:-1),kind=i8_kind))
      call mpp_sync()
      print *, 'distributed i8   value chksum:', mpp_chksum(sumi83d)
    else
      !call mpp_sync
      !open(10,file="/dev/null")
      !write(10,*) mpp_chksum(int(sumi4,kind=i8_kind))
      print *, mpp_chksum(int(sumi4,kind=i8_kind))
      call mpp_sync()
      !write(10,*) mpp_chksum(int(sumi4(num_entry:1:-1),kind=i8_kind))
      print *, mpp_chksum(int(sumi4(nx*ny*npz:1:-1),kind=i8_kind))
      call mpp_sync()
      !write(10,*) mpp_chksum(sumi8)
      print *, mpp_chksum(sumi8)
      call mpp_sync()
      !write(10,*) mpp_chksum(int(sumi43d,kind=i8_kind))
      print *, mpp_chksum(int(sumi43d,kind=i8_kind))
      call mpp_sync()
      !write(10,*) mpp_chksum(int(sumi43d(iec:isc:-1,jec:jsc:-1,km:1:-1),kind=i8_kind))
      print *, mpp_chksum(int(sumi43d(iec:isc:-1,jec:jsc:-1,km:1:-1),kind=i8_kind))
      call mpp_sync()
      !write(10,*) mpp_chksum(sumi83d)
      !call mpp_Sync()
      print *, mpp_chksum(sumi83d)
      !close(10)
    endif
    !!call mpp_sync()

    !!write(6,'(a,i4,2x,i,2x,i)') 'transfer function changes size of sumi8 array: ',mpp_pe(), size_bef, size_aft
  !  if (debug) then
  !    write(1000+mpp_pe(),'(i)') sumi8(:)
  !    write(2000+mpp_pe(),'(i)') sumi83d(:,:,:)
  !    write(3000+mpp_pe(),'(i)') sumi4(:)
  !    write(4000+mpp_pe(),'(i)') sumi43d(:,:,:)
  !  endif
    !!call mpp_domains_exit()

  end subroutine
end program
