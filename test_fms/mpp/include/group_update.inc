#define INDICES_MEM_ mem%lb(1):mem%ub(1), \
                     mem%lb(2):mem%ub(2), \
                     mem%lb(3):mem%ub(3)

#define INDICES_MEM_ISHIFT_ mem_ishift%lb(1):mem_ishift%ub(1), \
                            mem_ishift%lb(2):mem_ishift%ub(2), \
                            mem_ishift%lb(3):mem_ishift%ub(3)

#define INDICES_MEM_JSHIFT_ mem_jshift%lb(1):mem_jshift%ub(1), \
                            mem_jshift%lb(2):mem_jshift%ub(2), \
                            mem_jshift%lb(3):mem_jshift%ub(3)

#define INDICES_MEM_IJSHIFT_ mem_ijshift%lb(1):mem_ijshift%ub(1), \
                             mem_ijshift%lb(2):mem_ijshift%ub(2), \
                             mem_ijshift%lb(3):mem_ijshift%ub(3)

#define INDICES_COMPUTE_ compute%lb(1):compute%ub(1), \
                         compute%lb(2):compute%ub(2), \
                         compute%lb(3):compute%ub(3)

#define INDICES_COMPUTE_ISHIFT_ compute_ishift%lb(1):compute_ishift%ub(1), \
                                compute_ishift%lb(2):compute_ishift%ub(2), \
                                compute_ishift%lb(3):compute_ishift%ub(3)

#define INDICES_COMPUTE_JSHIFT_ compute_jshift%lb(1):compute_jshift%ub(1), \
                                compute_jshift%lb(2):compute_jshift%ub(2), \
                                compute_jshift%lb(3):compute_jshift%ub(3)

#define INDICES_COMPUTE_IJSHIFT_ compute_ijshift%lb(1):compute_ijshift%ub(1), \
                                 compute_ijshift%lb(2):compute_ijshift%ub(2), \
                                 compute_ijshift%lb(3):compute_ijshift%ub(3)

#define INDICES_COMPUTE_HALO_ compute_halo%lb(1):compute_halo%ub(1), \
                              compute_halo%lb(2):compute_halo%ub(2), \
                              compute_halo%lb(3):compute_halo%ub(3)

#define INDICES_DATA_ data%lb(1):data%ub(1), \
                      data%lb(2):data%ub(2), \
                      data%lb(3):data%ub(3)

#define INDICES_DATA_ISHIFT_ data_ishift%lb(1):data_ishift%ub(1), \
                             data_ishift%lb(2):data_ishift%ub(2), \
                             data_ishift%lb(3):data_ishift%ub(3)

#define INDICES_DATA_JSHIFT_ data_jshift%lb(1):data_jshift%ub(1), \
                             data_jshift%lb(2):data_jshift%ub(2), \
                             data_jshift%lb(3):data_jshift%ub(3)

#define INDICES_DATA_IJSHIFT_ data_ijshift%lb(1):data_ijshift%ub(1), \
                              data_ijshift%lb(2):data_ijshift%ub(2), \
                              data_ijshift%lb(3):data_ijshift%ub(3)

  subroutine TEST_GROUP_UPDATE_ (type, p)
    character(len=*), intent(in) :: type !< Name of the test to run
    integer, intent(in) :: p !< Permutation of array indices (ranges from 1 to 3!)

    real(FMS_TEST_KIND_), parameter :: zero = 0
    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: i, j, k, l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem

    integer, allocatable, dimension(:) :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:) :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:) :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:) :: layout2D, global_indices
    real(FMS_TEST_KIND_), allocatable, dimension(:,:,:,:) :: x1, y1, x2, y2
    real(FMS_TEST_KIND_), allocatable, dimension(:,:,:,:) :: a1, a2
    real(FMS_TEST_KIND_), allocatable, dimension(:,:,:) :: base
    integer :: id1, id2, id3
    logical :: folded_north
    logical :: cubic_grid
    character(len=3) :: text
    integer :: nx_save, ny_save
    type(mpp_group_update_type) :: group_update
    type(mpp_group_update_type), allocatable :: update_list(:)
    type(permutable_indices(3)) :: mem, mem_ishift, mem_jshift, mem_ijshift, &
                                   compute, compute_ishift, compute_jshift, compute_ijshift, compute_halo, &
                                   data, data_ishift, data_jshift, data_ijshift

    folded_north       = .false.
    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Folded-north' )
       ntiles = 1
       shift = 0
       num_contact = 2
       folded_north = .true.
       npes_per_tile = npes
       if(layout_tripolar(1)*layout_tripolar(2) == npes ) then
          layout = layout_tripolar
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       shift = 1
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from update_domains_performance ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_group_update: npes should be multiple of ntiles No test is done for '//trim(type))
          return
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(folded_north) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
       istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
            name = type, symmetry = .false.  )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )

    mem%lb = [ism, jsm, 1]
    mem%ub = [iem, jem, nz]
    call mem%permute(p)

    mem_ishift%lb = [ism, jsm, 1]
    mem_ishift%ub = [iem+shift, jem, nz]
    call mem_ishift%permute(p)

    mem_jshift%lb = [ism, jsm, 1]
    mem_jshift%ub = [iem, jem+shift, nz]
    call mem_jshift%permute(p)

    mem_ijshift%lb = [ism, jsm, 1]
    mem_ijshift%ub = [iem+shift, jem+shift, nz]
    call mem_ijshift%permute(p)

    compute%lb = [isc, jsc, 1]
    compute%ub = [iec, jec, nz]
    call compute%permute(p)

    compute_ishift%lb = [isc, jsc, 1]
    compute_ishift%ub = [iec+shift, jec, nz]
    call compute_ishift%permute(p)

    compute_jshift%lb = [isc, jsc, 1]
    compute_jshift%ub = [iec, jec+shift, nz]
    call compute_jshift%permute(p)

    compute_ijshift%lb = [isc, jsc, 1]
    compute_ijshift%ub = [iec+shift, jec+shift, nz]
    call compute_ijshift%permute(p)

    compute_halo%lb = [isc-1, jsc-1, 1]
    compute_halo%ub = [iec+1, jec+1, nz]
    call compute_halo%permute(p)

    data%lb = [isd, jsd, 1]
    data%ub = [ied, jed, nz]
    call data%permute(p)

    data_ishift%lb = [isd, jsd, 1]
    data_ishift%ub = [ied+shift, jed, nz]
    call data_ishift%permute(p)

    data_jshift%lb = [isd, jsd, 1]
    data_jshift%ub = [ied, jed+shift, nz]
    call data_jshift%permute(p)

    data_ijshift%lb = [isd, jsd, 1]
    data_ijshift%ub = [ied+shift, jed+shift, nz]
    call data_ijshift%permute(p)

    if(num_fields<1) then
       call mpp_error(FATAL, "test_mpp_domains: num_fields must be a positive integer")
    endif

    allocate(update_list(num_fields))

    id1 = mpp_clock_id( type//' group 2D', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' non-group 2D', flags=MPP_CLOCK_SYNC )
    id3 = mpp_clock_id( type//' non-block group 2D', flags=MPP_CLOCK_SYNC )

    allocate(a1(INDICES_MEM_, num_fields))
    allocate(x1(INDICES_MEM_ISHIFT_, num_fields))
    allocate(y1(INDICES_MEM_JSHIFT_, num_fields))
    allocate(a2(INDICES_MEM_, num_fields))
    allocate(x2(INDICES_MEM_ISHIFT_, num_fields))
    allocate(y2(INDICES_MEM_JSHIFT_, num_fields))
    allocate(base(INDICES_COMPUTE_IJSHIFT_))

    a1 = zero
    x1 = zero
    y1 = zero

    call arr_init(base)

    !--- Test for partial direction update
    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain, flags=WUPDATE+SUPDATE)
    end do

    do l = 1, num_fields
       a1(INDICES_COMPUTE_HALO_,l) = 999
       a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
    enddo

    a2 = a1
    call mpp_do_group_update(group_update, domain, zero)

    do l = 1, num_fields
       call mpp_update_domains( a2(:,:,:,l), domain, flags=WUPDATE+SUPDATE, complete=l==num_fields )
    enddo

    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l),type//' CENTER South West '//text)
    enddo

    call mpp_clear_group_update(group_update)

    !--- Test for DGRID update
    if(type == 'Cubic-Grid' ) then
       x1 = zero
       y1 = zero

       do l =1, num_fields
          call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=DGRID_NE)
       end do

       do l = 1, num_fields
          y1(INDICES_COMPUTE_ISHIFT_,l) = base(INDICES_COMPUTE_ISHIFT_) + l*1e3 + 1e6
          x1(INDICES_COMPUTE_JSHIFT_,l) = base(INDICES_COMPUTE_JSHIFT_) + l*1e3 + 2*1e6
       enddo
       x2 = x1
       y2 = y1
       call mpp_start_group_update(group_update, domain, zero)
       call mpp_complete_group_update(group_update, domain, zero)

       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=DGRID_NE, complete=l==num_fields )
       enddo

    !--- compare checksum
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(x1(INDICES_DATA_ISHIFT_,l),x2(INDICES_DATA_ISHIFT_,l),type//' DGRID X'//text)
          call compare_checksums(y1(INDICES_DATA_JSHIFT_,l),y2(INDICES_DATA_JSHIFT_,l),type//' DGRID Y'//text)
       enddo

       call mpp_clear_group_update(group_update)
    endif

    !--- Test for CGRID
    a1 = zero
    x1 = zero
    y1 = zero

    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain)
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE)
    end do

    do n = 1, num_iter
       a1 = zero
       x1 = zero
       y1 = zero

       do l = 1, num_fields
          a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
          x1(INDICES_COMPUTE_ISHIFT_,l) = base(INDICES_COMPUTE_ISHIFT_) + l*1e3 + 1e6
          y1(INDICES_COMPUTE_JSHIFT_,l) = base(INDICES_COMPUTE_JSHIFT_) + l*1e3 + 2*1e6
       enddo
       a2 = a1; x2 = x1; y2 = y1
       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update, domain, zero)
       call mpp_clock_end  (id1)

       call mpp_clock_begin(id2)
       do l = 1, num_fields
          call mpp_update_domains( a2(:,:,:,l), domain, complete=l==num_fields )
       enddo
       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=CGRID_NE, complete=l==num_fields )
       enddo
       call mpp_clock_end(id2)

       !--- compare checksum
       if( n == num_iter ) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l),type//' CENTER '//text)
             call compare_checksums(x1(INDICES_DATA_ISHIFT_,l),x2(INDICES_DATA_ISHIFT_,l),type//' CGRID X'//text)
             call compare_checksums(y1(INDICES_DATA_JSHIFT_,l),y2(INDICES_DATA_JSHIFT_,l),type//' CGRID Y'//text)
          enddo
       endif

       a1 = zero
       x1 = zero
       y1 = zero

       do l = 1, num_fields
          a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
          x1(INDICES_COMPUTE_ISHIFT_,l) = base(INDICES_COMPUTE_ISHIFT_) + l*1e3 + 1e6
          y1(INDICES_COMPUTE_JSHIFT_,l) = base(INDICES_COMPUTE_JSHIFT_) + l*1e3 + 2*1e6
       enddo
       call mpp_clock_begin(id3)
       call mpp_start_group_update(group_update, domain, zero)
       call mpp_complete_group_update(group_update, domain, zero)
       call mpp_clock_end  (id3)
       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l), type//' nonblock CENTER '//text)
          call compare_checksums(x1(INDICES_DATA_ISHIFT_,l),x2(INDICES_DATA_ISHIFT_,l), type//' nonblock CGRID X'//text)
          call compare_checksums(y1(INDICES_DATA_JSHIFT_,l),y2(INDICES_DATA_JSHIFT_,l), type//' nonblock CGRID Y'//text)
       enddo
       endif
    enddo

    call mpp_clear_group_update(group_update)

    !--- The following is to test overlapping start and complete
    if( num_fields > 1 ) then
       do l =1, num_fields
          call mpp_create_group_update(update_list(l), a1(:,:,:,l), domain)
          call mpp_create_group_update(update_list(l), x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE)
       end do

       do n = 1, num_iter

          a1 = zero
          x1 = zero
          y1 = zero

          do l = 1, num_fields
             a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
             x1(INDICES_COMPUTE_ISHIFT_,l) = base(INDICES_COMPUTE_ISHIFT_) + l*1e3 + 1e6
             y1(INDICES_COMPUTE_JSHIFT_,l) = base(INDICES_COMPUTE_JSHIFT_) + l*1e3 + 2*1e6
          enddo
          do l = 1, num_fields-1
             call mpp_start_group_update(update_list(l), domain, zero)
          enddo

          call mpp_complete_group_update(update_list(1), domain, zero)
          call mpp_start_group_update(update_list(num_fields), domain, zero)
          do l = 2, num_fields
             call mpp_complete_group_update(update_list(l), domain, zero)
          enddo
          !--- compare checksum
          if( n == num_iter ) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l), &
                                    type//' multiple nonblock CENTER '//text)
             call compare_checksums(x1(INDICES_DATA_ISHIFT_,l),x2(INDICES_DATA_ISHIFT_,l), &
                                    type//' multiple nonblock CGRID X'//text)
             call compare_checksums(y1(INDICES_DATA_JSHIFT_,l),y2(INDICES_DATA_JSHIFT_,l), &
                                    type//' multiple nonblock CGRID Y'//text)
          enddo
          endif
       enddo
    endif

    do l =1, num_fields
      call mpp_clear_group_update(update_list(l))
    enddo
    deallocate(update_list)

    !--- test scalar 4-D variable
    call mpp_create_group_update(group_update, a1(:,:,:,:), domain)

    a1 = zero
    x1 = zero
    y1 = zero

    do l = 1, num_fields
       a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
    enddo

    a2 = a1
    x2 = x1
    y2 = y1

    call mpp_clock_begin(id1)
    call mpp_do_group_update(group_update, domain, zero)
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id2)
    call mpp_update_domains( a2(:,:,:,:), domain )
    call mpp_clock_end(id2)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l),type//' 4D CENTER '//text)
    enddo

    a1 = 0
    do l = 1, num_fields
       a1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3
    enddo
    call mpp_clock_begin(id3)
    call mpp_start_group_update(group_update, domain, zero)
    call mpp_complete_group_update(group_update, domain, zero)
    call mpp_clock_end  (id3)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(INDICES_DATA_,l),a2(INDICES_DATA_,l), type//' nonblock 4D CENTER '//text)
    enddo



    !--- test for BGRID.
    deallocate(a1, x1, y1)
    deallocate(a2, x2, y2)
    call mpp_clear_group_update(group_update)

    allocate(a1(INDICES_MEM_IJSHIFT_, num_fields))
    allocate(x1(INDICES_MEM_IJSHIFT_, num_fields))
    allocate(y1(INDICES_MEM_IJSHIFT_, num_fields))
    allocate(a2(INDICES_MEM_IJSHIFT_, num_fields))
    allocate(x2(INDICES_MEM_IJSHIFT_, num_fields))
    allocate(y2(INDICES_MEM_IJSHIFT_, num_fields))

    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain, position=CORNER)
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=BGRID_NE)
    end do

    do n = 1, num_iter
       a1 = 0; x1 = 0; y1 = 0
       do l = 1, num_fields
          a1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3
          x1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3 + 1e6
          y1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3 + 2*1e6
       enddo

       a2 = a1
       x2 = x1
       y2 = y1

       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update, domain, zero)
       call mpp_clock_end  (id1)

       call mpp_clock_begin(id2)
       do l = 1, num_fields
          call mpp_update_domains( a2(:,:,:,l), domain, position=CORNER, complete=l==num_fields )
       enddo
       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=BGRID_NE, complete=l==num_fields )
       enddo
       call mpp_clock_end(id2)

       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(INDICES_DATA_IJSHIFT_,l),a2(INDICES_DATA_IJSHIFT_,l),type // ' CORNER ' // text)
          call compare_checksums(x1(INDICES_DATA_IJSHIFT_,l),x2(INDICES_DATA_IJSHIFT_,l),type // ' BGRID X' // text)
          call compare_checksums(y1(INDICES_DATA_IJSHIFT_,l),y2(INDICES_DATA_IJSHIFT_,l),type // ' BGRID Y' // text)
       enddo
       endif

       a1 = zero
       x1 = zero
       y1 = zero

       do l = 1, num_fields
          a1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3
          x1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3 + 1e6
          y1(INDICES_COMPUTE_IJSHIFT_,l) = base(INDICES_COMPUTE_IJSHIFT_) + l*1e3 + 2*1e6
       enddo
       call mpp_clock_begin(id3)
       call mpp_start_group_update(group_update, domain, zero)
       call mpp_complete_group_update(group_update, domain, zero)
       call mpp_clock_end  (id3)
       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(INDICES_DATA_IJSHIFT_,l),a2(INDICES_DATA_IJSHIFT_,l), &
                                 type//' nonblockCORNER '//text)
          call compare_checksums(x1(INDICES_DATA_IJSHIFT_,l),x2(INDICES_DATA_IJSHIFT_,l), &
                                 type//' nonblock BGRID X'//text)
          call compare_checksums(y1(INDICES_DATA_IJSHIFT_,l),y2(INDICES_DATA_IJSHIFT_,l), &
                                 type//' nonblock BGRID Y'//text)
       enddo
       endif
    enddo

    call mpp_clear_group_update(group_update)

    !-----------------------------------------------------------------------------
    !                   test for AGRID vector and scalar pair
    !-----------------------------------------------------------------------------
    deallocate(x1, y1)
    deallocate(x2, y2)

    allocate(x1(INDICES_MEM_, num_fields))
    allocate(y1(INDICES_MEM_, num_fields))
    allocate(x2(INDICES_MEM_, num_fields))
    allocate(y2(INDICES_MEM_, num_fields))

    x1 = zero
    y1 = zero

    do l = 1, num_fields
       x1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3 + 1e6
       y1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3 + 2*1e6
    enddo

    x2 = x1
    y2 = y1

    do l =1, num_fields
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=AGRID)
    end do

    do l = 1, num_fields
       call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=AGRID, complete=l==num_fields )
    enddo

    call mpp_start_group_update(group_update, domain, zero)
    call mpp_complete_group_update(group_update, domain, zero)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(INDICES_DATA_,l),x2(INDICES_DATA_,l),type//' AGRID X'//text)
       call compare_checksums(y1(INDICES_DATA_,l),y2(INDICES_DATA_,l),type//' AGRID Y'//text)
    enddo

    call mpp_clear_group_update(group_update)

    x1 = zero
    y1 = zero

    do l = 1, num_fields
       x1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3 + 1e6
       y1(INDICES_COMPUTE_,l) = base(INDICES_COMPUTE_) + l*1e3 + 2*1e6
    enddo

    x2 = x1
    y2 = y1

    do l =1, num_fields
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=AGRID, flags=SCALAR_PAIR)
    end do

    do l = 1, num_fields
       call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=AGRID, flags=SCALAR_PAIR, &
                              &  complete=l==num_fields)
    enddo

    call mpp_start_group_update(group_update, domain, zero)
    call mpp_complete_group_update(group_update, domain, zero)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(INDICES_DATA_,l),x2(INDICES_DATA_,l),type//' AGRID SCALAR_PAIR X'//text)
       call compare_checksums(y1(INDICES_DATA_,l),y2(INDICES_DATA_,l),type//' AGRID SCALAR_PAIR Y'//text)
    enddo

    call mpp_clear_group_update(group_update)

    deallocate(pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1)
    deallocate(istart2, iend2, jstart2, jend2)
    deallocate(layout2D, global_indices)

    deallocate(a1, x1, y1)
    deallocate(a2, x2, y2)
    deallocate(base)
    call mpp_deallocate_domain(domain)

end subroutine TEST_GROUP_UPDATE_

#undef INDICES_MEM_
#undef INDICES_MEM_ISHIFT_
#undef INDICES_MEM_JSHIFT_
#undef INDICES_MEM_IJSHIFT_
#undef INDICES_COMPUTE_
#undef INDICES_COMPUTE_ISHIFT_
#undef INDICES_COMPUTE_JSHIFT_
#undef INDICES_COMPUTE_IJSHIFT_
#undef INDICES_COMPUTE_HALO_
#undef INDICES_DATA_
#undef INDICES_DATA_ISHIFT_
#undef INDICES_DATA_JSHIFT_
#undef INDICES_DATA_IJSHIFT_
