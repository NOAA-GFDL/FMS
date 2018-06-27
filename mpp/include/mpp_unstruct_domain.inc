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

  !#####################################################################
  subroutine mpp_define_unstruct_domain(UG_domain, SG_domain, npts_tile, grid_nlev, ndivs, npes_io_group, grid_index, name)
     type(domainUG),   intent(inout) :: UG_domain
     type(domain2d), target,     intent(in) :: SG_domain 
     integer,                    intent(in) :: npts_tile(:) ! number of unstructured points on each tile
     integer,                    intent(in) :: grid_nlev(:) ! number of levels in each unstructured grid.
     integer,                    intent(in) :: ndivs       
     integer,                    intent(in) :: npes_io_group  ! number of processors in a io group. Only pe with same tile_id 
                                                              ! in the same group
     integer,                    intent(in) :: grid_index(:)
     character(len=*), optional, intent(in) :: name
     integer, dimension(size(npts_tile(:))) :: ndivs_tile, pe_start, pe_end
     integer, dimension(0:ndivs-1)          :: ibegin, iend, costs_list
     integer :: ntiles, ntotal_pts, ndivs_used, max_npts, cur_tile, cur_npts
     integer :: n, ts, te, p, pos, tile_id, ngroup, group_id, my_pos, i
     integer :: npes_in_group, is, ie, ntotal_costs, max_cost, cur_cost, costs_left
     integer :: npts_left, ndiv_left, cur_pos, ndiv, prev_cost, ioff
     real    :: avg_cost
     integer :: costs(size(npts_tile(:)))

     UG_domain%SG_domain => SG_domain
     ntiles = size(npts_tile(:))
     UG_domain%ntiles = ntiles

     !--- total number of points must be no less than ndivs
     if(sum(npts_tile)<ndivs) call mpp_error(FATAL, "mpp_define_unstruct_domain: total number of points is less than ndivs")
     !--- We are assuming nlev on each grid is at least one.
     do n = 1, size(grid_nlev(:))
        if(grid_nlev(n) < 1) call mpp_error(FATAL, "mpp_define_unstruct_domain: grid_nlev at some point is less than 1")
     enddo

     !-- costs for each tile.
     pos = 0
     do n = 1, ntiles
        costs(n) = 0
        do i = 1, npts_tile(n)
           pos = pos + 1
           costs(n) = costs(n) + grid_nlev(pos)
        enddo
     enddo
     ! compute number of divisions for each tile.
     ntotal_costs = sum(costs)
     !--- get the upper limit of ndivs for each tile.
     do n = 1, ntiles
        ndivs_tile(n) = CEILING(real(costs(n)*ndivs)/ntotal_costs)
     enddo

     ndivs_used = sum(ndivs_tile)
     do while (ndivs_used > ndivs)
        max_cost = 0
        cur_tile = 0
        do n = 1, ntiles
           if( ndivs_tile(n) > 1 ) then
              cur_cost = CEILING(real(costs(n))/(ndivs_tile(n)-1))
              if( max_cost == 0 .OR. cur_cost<max_cost) then
                 max_cost = cur_cost
                 cur_tile = n
              endif
           endif
        enddo
        ndivs_used = ndivs_used-1
        ndivs_tile(cur_tile) = ndivs_tile(cur_tile) - 1    
     enddo

     te = -1
     ioff = 0
     do n = 1, ntiles
        ts = te + 1
        te = ts + ndivs_tile(n) - 1
        costs_left = costs(n)
        ndiv_left = ndivs_tile(n)
        npts_left = npts_tile(n)
        cur_pos = 1
        do ndiv = 1, ndivs_tile(n)
           cur_cost = 0
           ibegin(ts+ndiv-1) = cur_pos
           avg_cost = real(costs_left)/ndiv_left
           do i = cur_pos, npts_tile(n)
              cur_cost = cur_cost + grid_nlev(i+ioff)
              costs_left = costs_left - grid_nlev(i+ioff)
              if(npts_left < ndiv_left ) then
                 call mpp_error(FATAL, "mpp_define_unstruct_domain: npts_left < ndiv_left")
              else if(npts_left == ndiv_left ) then
                 cur_pos = i + 1
                 exit
              else if(cur_cost .GE. avg_cost) then
                 prev_cost = cur_cost - grid_nlev(i+ioff)
                 if(i==cur_pos) then
                    cur_pos = i + 1
                    exit
                 else if( cur_cost - avg_cost .LE. avg_cost - prev_cost ) then
                    cur_pos = i + 1
                    exit
                 else
                    cur_pos = i
                    cur_cost = prev_cost
                    costs_left = costs_left + grid_nlev(i+ioff)
                    npts_left = npts_left+1
                    exit
                 endif
              endif
              npts_left = npts_left-1
           enddo
           iend(ts+ndiv-1) = cur_pos - 1
           costs_list(ts+ndiv-1) = cur_cost
           ndiv_left = ndiv_left-1
           npts_left = npts_left-1
        enddo
        pe_start(n) = ts
        pe_end(n) = te
        ioff = ioff+ npts_tile(n)
     enddo
     allocate(UG_domain%list(0:ndivs-1))
     do p = 0, ndivs-1
        UG_domain%list(p)%compute%begin = ibegin(p)
        UG_domain%list(p)%compute%end = iend(p)
        UG_domain%list(p)%compute%size = UG_domain%list(p)%compute%end - UG_domain%list(p)%compute%begin + 1
        UG_domain%list(p)%compute%max_size = 0
        UG_domain%list(p)%pos = p
        UG_domain%list(p)%pe = p + mpp_root_pe()
        pos = 0
        do n = 1, ntiles
           if( p .GE. pe_start(n) .AND. p .LE. pe_end(n) ) then
              UG_domain%list(p)%tile_id = n
              exit
           endif
           pos = pos + npts_tile(n)
        enddo
        is = UG_domain%list(p)%compute%begin+pos
        ie = UG_domain%list(p)%compute%end+pos
        UG_domain%list(p)%compute%begin_index = minval(grid_index(is:ie))
        UG_domain%list(p)%compute%end_index = maxval(grid_index(is:ie))
     enddo

     !--- write out domain decomposition from root pe
     if(mpp_pe() == mpp_root_pe() .and. present(name)) then
        write(stdout(),*) "unstruct domain name = ", trim(name)
        write(stdout(),*) UG_domain%list(:)%compute%size
     endif

     pos = mpp_pe() - mpp_root_pe() 
     UG_domain%pe = mpp_pe()
     UG_domain%pos = pos
     UG_domain%tile_id = UG_domain%list(pos)%tile_id
     p = pe_start(UG_domain%tile_id)
     UG_domain%tile_root_pe = UG_domain%list(p)%pe
     UG_domain%tile_npes = pe_end(UG_domain%tile_id) - pe_start(UG_domain%tile_id) + 1
     UG_domain%compute = UG_domain%list(pos)%compute
     UG_domain%compute%max_size = MAXVAL( UG_domain%list(:)%compute%size )
     UG_domain%global%begin = 1
     UG_domain%global%end   = npts_tile(UG_domain%tile_id)
     UG_domain%global%size  = UG_domain%global%end - UG_domain%global%begin + 1
     UG_domain%global%max_size = -1   ! currently this is not supposed to be used.
     pos = 0
     do n = 1, UG_domain%tile_id-1
        pos = pos + npts_tile(n)
     enddo
     UG_domain%global%begin_index = grid_index(pos+1)
     UG_domain%global%end_index = grid_index(pos+npts_tile(n))

     allocate(UG_domain%grid_index(UG_domain%compute%size))
     do n = 1, UG_domain%compute%size
        UG_domain%grid_index(n) = grid_index(pos+UG_domain%compute%begin+n-1)
     enddo

     !--- define io_domain
     allocate(UG_domain%io_domain)
     tile_id = UG_domain%tile_id
     UG_domain%io_domain%pe = UG_domain%pe
     !--- figure out number groups for current tile
     if(npes_io_group == 0) then
        ngroup = 1
     else
        ngroup = CEILING(real(ndivs_tile(tile_id))/ npes_io_group)
     endif

!----------
!ug support
     UG_domain%npes_io_group = npes_io_group
     UG_domain%io_layout = ngroup
!----------

     call mpp_compute_extent(1, ndivs_tile(tile_id), ngroup, ibegin(0:ngroup-1), iend(0:ngroup-1))
     my_pos = UG_domain%pe - UG_domain%tile_root_pe + 1
     do n = 0, ngroup-1
        if( my_pos .GE. ibegin(n) .AND. my_pos .LE. iend(n) ) then
           group_id = n
           exit
        endif
     enddo
     
     UG_domain%io_domain%tile_id            = group_id+1
     UG_domain%io_domain%compute            = UG_domain%compute
     UG_domain%io_domain%pe                 = UG_domain%pe
     UG_domain%io_domain%pos                = my_pos - ibegin(group_id) + 1
     UG_domain%io_domain%tile_root_pe       = ibegin(group_id) + UG_domain%tile_root_pe - 1
     pos = UG_domain%io_domain%tile_root_pe - mpp_root_pe() 
     UG_domain%io_domain%global%begin       = UG_domain%list(pos)%compute%begin
     UG_domain%io_domain%global%begin_index = UG_domain%list(pos)%compute%begin_index
     pos = iend(group_id) + UG_domain%tile_root_pe - mpp_root_pe() - 1
     UG_domain%io_domain%global%end         = UG_domain%list(pos)%compute%end
     UG_domain%io_domain%global%end_index   = UG_domain%list(pos)%compute%end_index
     UG_domain%io_domain%global%size        = UG_domain%io_domain%global%end -  UG_domain%io_domain%global%begin + 1
     
     npes_in_group = iend(group_id) - ibegin(group_id) + 1
     allocate(UG_domain%io_domain%list(0:npes_in_group-1))
     do n = 0, npes_in_group-1
        pos = UG_domain%io_domain%tile_root_pe - mpp_root_pe() + n
        UG_domain%io_domain%list(n)%compute = UG_domain%list(pos)%compute
        UG_domain%io_domain%list(n)%pos = n
        UG_domain%io_domain%list(n)%pe = UG_domain%list(pos)%pe
        UG_domain%io_domain%list(n)%tile_id = group_id+1
     enddo    
 
     call compute_overlap_SG2UG(UG_domain, SG_domain)
     call compute_overlap_UG2SG(UG_domain)

     return

  end subroutine mpp_define_unstruct_domain


  !####################################################################
  subroutine compute_overlap_SG2UG(UG_domain, SG_domain)
     type(domainUG),   intent(inout) :: UG_domain
     type(domain2d),             intent(in) :: SG_domain
     integer, dimension(0:size(SG_domain%list(:))-1) :: send_cnt, recv_cnt
     integer, dimension(0:size(SG_domain%list(:))-1) :: send_buffer_pos, recv_buffer_pos
     integer, dimension(:), allocatable              :: send_buffer, recv_buffer, index_list
     integer, dimension(:), allocatable              :: buffer_pos
     integer :: tile_id, nlist, nxg, begin_index, end_index, i, j
     integer :: m, n, list, l, isc, iec, jsc, jec, ibegin, iend, grid_index
     integer :: nrecv, nsend, send_pos, recv_pos, pos

     !--- figure out the recv index information.
     tile_id = UG_domain%tile_id
     nlist = size(SG_domain%list(:))
     nxg = SG_domain%x(1)%global%size
     begin_index = UG_domain%compute%begin_index
     end_index = UG_domain%compute%end_index
     pos = 0
     recv_cnt = 0
     allocate(index_list(UG_domain%compute%size))
     allocate(send_buffer(UG_domain%compute%size))
     index_list = -1
     do n = 0, nlist-1
        if(SG_domain%list(n)%tile_id(1) .NE. tile_id) cycle
        isc = SG_domain%list(n)%x(1)%compute%begin; iec = SG_domain%list(n)%x(1)%compute%end
        jsc = SG_domain%list(n)%y(1)%compute%begin; jec = SG_domain%list(n)%y(1)%compute%end
        ibegin = (jsc-1)*nxg + isc
        iend   = (jec-1)*nxg + iec
        if(ibegin > end_index .OR. iend < begin_index) cycle
        do l = 1, UG_domain%compute%size
           grid_index = UG_domain%grid_index(l)
           i = mod((grid_index-1), nxg) + 1
           j = (grid_index-1)/nxg + 1 
           if( i .GE. isc .AND. i .LE. iec .and. j .GE. jsc .AND. j .LE. jec ) then
              recv_cnt(n) = recv_cnt(n) + 1
              pos = pos + 1
              if(pos > UG_domain%compute%size) call mpp_error(FATAL, &
                  'compute_overlap_SG2UG: pos > UG_domain%compute%size')
              index_list(pos) = l
              send_buffer(pos) = grid_index
           endif
        enddo
     enddo

     !--- make sure sum(recv_cnt) == UG_domain%compute%size
     if( UG_domain%compute%size .NE. sum(recv_cnt) ) then
        print*,"pe=", mpp_pe(), UG_domain%compute%size, sum(recv_cnt)
        call mpp_error(FATAL, &
          "compute_overlap_SG2UG: UG_domain%compute%size .NE. sum(recv_cnt)")
     endif
     allocate(buffer_pos(0:nlist-1))
     pos = 0
     do list = 0,nlist-1
        buffer_pos(list) = pos
        pos = pos + recv_cnt(list)
     enddo

     nrecv = count( recv_cnt > 0 )
     UG_domain%SG2UG%nrecv = nrecv
     allocate(UG_domain%SG2UG%recv(nrecv))
     nrecv = 0
     pos = 0
     do list = 0,nlist-1
        m = mod( SG_domain%pos+nlist-list, nlist )
        if( recv_cnt(m) > 0 ) then
           nrecv = nrecv+1
           UG_domain%SG2UG%recv(nrecv)%count = recv_cnt(m)
           UG_domain%SG2UG%recv(nrecv)%pe = UG_domain%list(m)%pe
           allocate(UG_domain%SG2UG%recv(nrecv)%i(recv_cnt(m)))
           pos = buffer_pos(m)
           do l = 1, recv_cnt(m)
              pos = pos + 1
              UG_domain%SG2UG%recv(nrecv)%i(l) = index_list(pos)
           enddo
        endif
     enddo

     !--- figure out the send index information.
     send_cnt = recv_cnt
     recv_cnt = 0
     call mpp_alltoall(send_cnt,1,recv_cnt,1)
     !--- make sure sum(send_cnt) == UG_domain%compute%size
     if( UG_domain%compute%size .NE. sum(send_cnt) ) call mpp_error(FATAL, &
          "compute_overlap_SG2UG: UG_domain%compute%size .NE. sum(send_cnt)")
     allocate(recv_buffer(sum(recv_cnt)))
     send_buffer_pos = 0; recv_buffer_pos = 0
     send_pos = 0; recv_pos = 0
     do n = 0, nlist-1
        if(send_cnt(n) > 0) then
           send_buffer_pos(n) = send_pos
           send_pos = send_pos + send_cnt(n)
        endif
        if(recv_cnt(n) > 0) then
           recv_buffer_pos(n) = recv_pos 
           recv_pos = recv_pos + recv_cnt(n)
        endif
     enddo    

     call mpp_alltoall(send_buffer, send_cnt, send_buffer_pos, &
                       recv_buffer, recv_cnt, recv_buffer_pos)

     nsend = count( recv_cnt(:) > 0 )
     UG_domain%SG2UG%nsend = nsend
     allocate(UG_domain%SG2UG%send(nsend))
     nsend = 0
     isc = SG_domain%x(1)%compute%begin
     jsc = SG_domain%y(1)%compute%begin
     do list = 0,nlist-1
        m = mod( SG_domain%pos+list, nlist )
        if( recv_cnt(m) > 0 ) then
           nsend = nsend+1
           UG_domain%SG2UG%send(nsend)%count = recv_cnt(m)
           UG_domain%SG2UG%send(nsend)%pe = UG_domain%list(m)%pe
           allocate(UG_domain%SG2UG%send(nsend)%i(recv_cnt(m)))
           allocate(UG_domain%SG2UG%send(nsend)%j(recv_cnt(m)))
           pos = recv_buffer_pos(m)
           do l = 1, recv_cnt(m)
              grid_index = recv_buffer(pos+l)
              UG_domain%SG2UG%send(nsend)%i(l) = mod(grid_index-1,nxg) + 1
              UG_domain%SG2UG%send(nsend)%j(l) = (grid_index-1)/nxg + 1
           enddo
        endif
     enddo
     deallocate(send_buffer, recv_buffer, index_list, buffer_pos)

return

  end subroutine compute_overlap_SG2UG

  !####################################################################
  subroutine compute_overlap_UG2SG(UG_domain)
     type(domainUG),   intent(inout) :: UG_domain

     !--- UG2SG is the reverse of SG2UG
     UG_domain%UG2SG%nsend = UG_domain%SG2UG%nrecv
     UG_domain%UG2SG%send => UG_domain%SG2UG%recv
     UG_domain%UG2SG%nrecv = UG_domain%SG2UG%nsend
     UG_domain%UG2SG%recv => UG_domain%SG2UG%send

     return

  end subroutine compute_overlap_UG2SG

  !####################################################################
  subroutine mpp_get_UG_SG_domain(UG_domain,SG_domain)
     type(domainUG),   intent(inout) :: UG_domain
     type(domain2d),   pointer       :: SG_domain

     SG_domain => UG_domain%SG_domain

     return

  end subroutine mpp_get_UG_SG_domain

  !####################################################################
  function mpp_get_UG_io_domain(domain)
     type(domainUG), intent(in) :: domain
     type(domainUG), pointer    :: mpp_get_UG_io_domain

     if(ASSOCIATED(domain%io_domain)) then
        mpp_get_UG_io_domain => domain%io_domain
     else
        call mpp_error(FATAL, "mpp_get_UG_io_domain: io_domain is not defined, contact developer")
     endif

  end function mpp_get_UG_io_domain

  !#####################################################################
  subroutine mpp_get_UG_compute_domain( domain, begin, end, size)
    type(domainUG),  intent(in) :: domain
    integer, intent(out), optional :: begin, end, size

    if( PRESENT(begin)     )begin     = domain%compute%begin
    if( PRESENT(end)       )end       = domain%compute%end
    if( PRESENT(size)      )size      = domain%compute%size
    return
  end subroutine mpp_get_UG_compute_domain

  !#####################################################################
  subroutine mpp_get_UG_global_domain( domain, begin, end, size)
    type(domainUG),  intent(in) :: domain
    integer, intent(out), optional :: begin, end, size

    if( PRESENT(begin)     )begin     = domain%global%begin
    if( PRESENT(end)       )end       = domain%global%end
    if( PRESENT(size)      )size      = domain%global%size
    return
  end subroutine mpp_get_UG_global_domain

  !#####################################################################
  subroutine mpp_get_UG_compute_domains( domain, begin, end, size )
    type(domainUG),                   intent(in) :: domain
    integer, intent(out), optional, dimension(:) :: begin, end, size

    !we use shape instead of size for error checks because size is used as an argument
    if( PRESENT(begin) )then
       if( any(shape(begin).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'mpp_get_UG_compute_domains: begin array size does not match domain.' )
       begin(:) = domain%list(:)%compute%begin
    end if
    if( PRESENT(end) )then
       if( any(shape(end).NE.shape(domain%list)) ) &
            call mpp_error( FATAL, 'mpp_get_UG_compute_domains: end array size does not match domain.' )
            end(:) = domain%list(:)%compute%end
    end if
    if( PRESENT(size) )then
       if( any(shape(size).NE.shape(domain%list)) ) &
           call mpp_error( FATAL, 'mpp_get_UG_compute_domains: size array size does not match domain.' )
       size(:) = domain%list(:)%compute%size
    end if
    return
  end subroutine mpp_get_UG_compute_domains

  !#####################################################################
  subroutine mpp_get_UG_domains_index( domain, begin, end)
    type(domainUG),         intent(in) :: domain
    integer, intent(out), dimension(:) :: begin, end

    !we use shape instead of size for error checks because size is used as an argument
    if( any(shape(begin).NE.shape(domain%list)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_compute_domains: begin array size does not match domain.' )
    begin(:) = domain%list(:)%compute%begin_index
    if( any(shape(end).NE.shape(domain%list)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_compute_domains: end array size does not match domain.' )
         end(:) = domain%list(:)%compute%end_index
    return
  end subroutine mpp_get_UG_domains_index

  !#####################################################################
  function mpp_get_UG_domain_ntiles(domain)
    type(domainUG),  intent(in) :: domain
    integer :: mpp_get_UG_domain_ntiles

    mpp_get_UG_domain_ntiles = domain%ntiles
    return
  end function mpp_get_UG_domain_ntiles

  !#######################################################################
  subroutine mpp_get_ug_domain_tile_list(domain, tiles)
     type(domainUG), intent(in) :: domain
     integer,     intent(inout) :: tiles(:)
     integer                    :: i

     if( size(tiles(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'mpp_get_ug_domain_tile_list: tiles array size does not match domain.' )
     do i = 1, size(tiles(:))
        tiles(i) = domain%list(i-1)%tile_id
     end do

  end subroutine mpp_get_ug_domain_tile_list

  !#####################################################################
  function mpp_get_UG_domain_tile_id(domain)
    type(domainUG),  intent(in) :: domain
    integer :: mpp_get_UG_domain_tile_id

    mpp_get_UG_domain_tile_id = domain%tile_id
    return
  end function mpp_get_UG_domain_tile_id

  !####################################################################
  function mpp_get_UG_domain_npes(domain)
     type(domainUG), intent(in) :: domain
     integer :: mpp_get_UG_domain_npes

     mpp_get_UG_domain_npes = size(domain%list(:))
     return

  end function mpp_get_UG_domain_npes


  !####################################################################
  subroutine mpp_get_UG_domain_pelist( domain, pelist)
     type(domainUG), intent(in) :: domain
     integer,              intent(out) :: pelist(:)

     if( size(pelist(:)).NE.size(domain%list(:)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_domain_pelist: pelist array size does not match domain.' )

     pelist(:) = domain%list(:)%pe
     return

  end subroutine mpp_get_UG_domain_pelist

  !###################################################################
  subroutine mpp_get_UG_domain_tile_pe_inf( domain, root_pe, npes, pelist)
     type(domainUG),     intent(in) :: domain
     integer, optional, intent(out) :: root_pe, npes
     integer, optional, intent(out) :: pelist(:)

     if(present(root_pe)) root_pe = domain%tile_root_pe
     if(present(npes)) root_pe = domain%tile_npes

     if(present(pelist)) then
        if( size(pelist(:)).NE. domain%tile_npes ) &
           call mpp_error( FATAL, 'mpp_get_UG_domain_tile_pe_inf: pelist array size does not match domain.' )
        pelist(:) = domain%list(domain%pos:domain%pos+domain%tile_npes-1)%pe
     endif
     return

  end subroutine mpp_get_UG_domain_tile_pe_inf


  !####################################################################
  subroutine mpp_get_UG_domain_grid_index( domain, grid_index)
     type(domainUG), intent(in) :: domain
     integer,              intent(out) :: grid_index(:)

     if( size(grid_index(:)).NE.size(domain%grid_index(:)) ) &
         call mpp_error( FATAL, 'mpp_get_UG_domain_grid_index: grid_index array size does not match domain.' )

     grid_index(:) = domain%grid_index(:)
     return

  end subroutine mpp_get_UG_domain_grid_index

  !###################################################################
  subroutine mpp_define_null_UG_domain(domain)
     type(domainUG), intent(inout) :: domain

     domain%global%begin  = -1; domain%global%end  = -1; domain%global%size = 0
     domain%compute%begin = -1; domain%compute%end = -1; domain%compute%size = 0
     domain%pe = NULL_PE
     domain%ntiles = -1
     domain%pos = -1
     domain%tile_id = -1
     domain%tile_root_pe = -1

  end subroutine mpp_define_null_UG_domain

!##############################################################################
    subroutine mpp_broadcast_domain_ug( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domainUG), intent(inout) :: domain
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(7) :: msg, info         !pe and compute domain of each item in list
      integer                :: errunit

      errunit = stderr()
      if( .NOT.module_is_initialized ) &
                 call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN_ug: You must first call mpp_domains_init.' )

!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)

!am I part of this domain?
      native = ASSOCIATED(domain%list)

!set local list size
      if( native )then
          listsize = size(domain%list(:))
      else
          listsize = 0
      end if
      call mpp_max(listsize)

      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          domain%ntiles = -1
          domain%compute%begin =  1
          domain%compute%end   = -1
          domain%compute%begin_index =  1
          domain%compute%end_index   = -1
          domain%global %begin = -1
          domain%global %end   = -1
          domain%tile_id       = -1
          domain%tile_root_pe  = -1
      end if
!initialize values in info
      info(1) = domain%pe
      info(2) = domain%pos
      info(3) = domain%tile_id
      call mpp_get_UG_compute_domain( domain, info(4), info(5))
      info(6) = domain%compute%begin_index
      info(7) = domain%compute%end_index
!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( mpp_pe().EQ.pes(n) .AND. debug )write( errunit,* )'PE ', mpp_pe(), 'broadcasting msg ', msg
         call mpp_broadcast( msg, 7, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe            = msg(1)
             domain%list(listpos)%pos           = msg(2)
             domain%list(listpos)%tile_id       = msg(3)
             domain%list(listpos)%compute%begin = msg(4)
             domain%list(listpos)%compute%end   = msg(5)
             domain%list(listpos)%compute%begin_index = msg(6)
             domain%list(listpos)%compute%end_index   = msg(7)
             listpos = listpos + 1
             if( debug )write( errunit,* )'PE ', mpp_pe(), 'received domain from PE ', msg(1), 'ls,le=', msg(4:5)
         end if
      end do

    end subroutine mpp_broadcast_domain_ug

!------------------------------------------------------------------------------
function mpp_domain_UG_is_tile_root_pe(domain) result(is_root)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    logical(INT_KIND)         :: is_root

    if (domain%pe .eq. domain%tile_root_pe) then
        is_root = .true.
    else
        is_root = .false.
    endif

    return
end function mpp_domain_UG_is_tile_root_pe

!------------------------------------------------------------------------------
!HELP: There needs to be a subroutine to return the "io_layout" for
!      an unstructured domain, so I made one.  Someone should check
!      to see if this is correct.
function mpp_get_io_domain_UG_layout(domain) result(io_layout)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    integer(INT_KIND)         :: io_layout

    io_layout = domain%io_layout

    return
end function


!------------------------------------------------------------------
subroutine deallocate_unstruct_overlap_type(overlap)
  type(unstruct_overlap_type), intent(inout) :: overlap

  if(associated(overlap%i)) deallocate(overlap%i)
  if(associated(overlap%j)) deallocate(overlap%j)

end subroutine deallocate_unstruct_overlap_type

!------------------------------------------------------------------
subroutine deallocate_unstruct_pass_type(passobj)
  type(unstruct_pass_type), intent(inout) :: passobj
  integer :: n

  do n = 1, passobj%nsend
     call deallocate_unstruct_overlap_type(passobj%send(n))
  enddo
  do n = 1, passobj%nrecv
     call deallocate_unstruct_overlap_type(passobj%recv(n))
  enddo

  if(associated(passobj%send)) deallocate(passobj%send)
  if(associated(passobj%recv)) deallocate(passobj%recv)
 
end subroutine deallocate_unstruct_pass_type

!------------------------------------------------------------------
subroutine mpp_deallocate_domainUG(domain)

   !<Inputs/Outputs
    type(domainUG),intent(inout) :: domain

   !<Local variables
    integer(INT_KIND) :: i !<Loop variable.

    if (associated(domain%list)) then
        deallocate(domain%list)
        domain%list => null()
    endif

    if (associated(domain%io_domain)) then
        if (associated(domain%io_domain%list)) then
            deallocate(domain%io_domain%list)
            domain%io_domain%list => null()
        endif
        deallocate(domain%io_domain)
        domain%io_domain => null()
    endif

    call deallocate_unstruct_pass_type(domain%SG2UG)
    call deallocate_unstruct_pass_type(domain%UG2SG)

    if (associated(domain%grid_index)) then
        deallocate(domain%grid_index)
        domain%grid_index => null()
    endif

    if (associated(domain%SG_domain)) then
        domain%SG_domain => null()
    endif

    return
end subroutine mpp_deallocate_domainUG

  !###################################################################
  !> Overload the .eq. for UG
  function mpp_domainUG_eq( a, b )
    logical                    :: mpp_domainUG_eq
    type(domainUG), intent(in) :: a, b

    if (associated(a%SG_domain) .and. associated(b%SG_domain)) then
        if (a%SG_domain .ne. b%SG_domain) then
            mpp_domainUG_eq = .false.
            return
        endif
    elseif (associated(a%SG_domain) .and. .not. associated(b%SG_domain)) then
        mpp_domainUG_eq = .false.
        return
    elseif (.not. associated(a%SG_domain) .and. associated(b%SG_domain)) then
        mpp_domainUG_eq = .false.
        return
    endif

    mpp_domainUG_eq = (a%npes_io_group .EQ. b%npes_io_group) .AND. &
                      (a%pos .EQ. b%pos)                     .AND. &
                      (a%ntiles .EQ. b%ntiles)               .AND. &
                      (a%tile_id .EQ. b%tile_id)             .AND. &
                      (a%tile_npes .EQ. b%tile_npes)         .AND. &
                      (a%tile_root_pe .EQ. b%tile_root_pe) 

    if(.not. mpp_domainUG_eq) return

    mpp_domainUG_eq = ( a%compute%begin.EQ.b%compute%begin .AND. &
         a%compute%end  .EQ.b%compute%end   .AND. &
         a%global%begin .EQ.b%global%begin  .AND. &
         a%global%end   .EQ.b%global%end    .AND. &
         a%SG2UG%nsend  .EQ.b%SG2UG%nsend   .AND. &
         a%SG2UG%nrecv  .EQ.b%SG2UG%nrecv   .AND. &    
         a%UG2SG%nsend  .EQ.b%UG2SG%nsend   .AND. &    
         a%UG2SG%nrecv  .EQ.b%UG2SG%nrecv         &
    )
    
    return
  end function mpp_domainUG_eq

  !> Overload the .ne. for UG
  function mpp_domainUG_ne( a, b )
    logical                    :: mpp_domainUG_ne
    type(domainUG), intent(in) :: a, b

    mpp_domainUG_ne = .NOT. ( a.EQ.b )
    return
  end function mpp_domainUG_ne

#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef mpp_pass_SG_to_UG_2D_
#define mpp_pass_SG_to_UG_2D_ mpp_pass_SG_to_UG_r8_2d
#undef mpp_pass_SG_to_UG_3D_
#define mpp_pass_SG_to_UG_3D_ mpp_pass_SG_to_UG_r8_3d
#undef mpp_pass_UG_to_SG_2D_
#define mpp_pass_UG_to_SG_2D_ mpp_pass_UG_to_SG_r8_2d
#undef mpp_pass_UG_to_SG_3D_
#define mpp_pass_UG_to_SG_3D_ mpp_pass_UG_to_SG_r8_3d
#include <mpp_unstruct_pass_data.h>

#ifdef OVERLOAD_R4
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef mpp_pass_SG_to_UG_2D_
#define mpp_pass_SG_to_UG_2D_ mpp_pass_SG_to_UG_r4_2d
#undef mpp_pass_SG_to_UG_3D_
#define mpp_pass_SG_to_UG_3D_ mpp_pass_SG_to_UG_r4_3d
#undef mpp_pass_UG_to_SG_2D_
#define mpp_pass_UG_to_SG_2D_ mpp_pass_UG_to_SG_r4_2d
#undef mpp_pass_UG_to_SG_3D_
#define mpp_pass_UG_to_SG_3D_ mpp_pass_UG_to_SG_r4_3d
#include <mpp_unstruct_pass_data.h>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef mpp_pass_SG_to_UG_2D_
#define mpp_pass_SG_to_UG_2D_ mpp_pass_SG_to_UG_i4_2d
#undef mpp_pass_SG_to_UG_3D_
#define mpp_pass_SG_to_UG_3D_ mpp_pass_SG_to_UG_i4_3d
#undef mpp_pass_UG_to_SG_2D_
#define mpp_pass_UG_to_SG_2D_ mpp_pass_UG_to_SG_i4_2d
#undef mpp_pass_UG_to_SG_3D_
#define mpp_pass_UG_to_SG_3D_ mpp_pass_UG_to_SG_i4_3d
#include <mpp_unstruct_pass_data.h>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#undef mpp_pass_SG_to_UG_2D_
#define mpp_pass_SG_to_UG_2D_ mpp_pass_SG_to_UG_l4_2d
#undef mpp_pass_SG_to_UG_3D_
#define mpp_pass_SG_to_UG_3D_ mpp_pass_SG_to_UG_l4_3d
#undef mpp_pass_UG_to_SG_2D_
#define mpp_pass_UG_to_SG_2D_ mpp_pass_UG_to_SG_l4_2d
#undef mpp_pass_UG_to_SG_3D_
#define mpp_pass_UG_to_SG_3D_ mpp_pass_UG_to_SG_l4_3d
#include <mpp_unstruct_pass_data.h>

#undef MPP_GLOBAL_FIELD_UG_2D_
#define MPP_GLOBAL_FIELD_UG_2D_ mpp_global_field2D_ug_r8_2d
#undef MPP_GLOBAL_FIELD_UG_3D_
#define MPP_GLOBAL_FIELD_UG_3D_ mpp_global_field2D_ug_r8_3d
#undef MPP_GLOBAL_FIELD_UG_4D_
#define MPP_GLOBAL_FIELD_UG_4D_ mpp_global_field2D_ug_r8_4d
#undef MPP_GLOBAL_FIELD_UG_5D_
#define MPP_GLOBAL_FIELD_UG_5D_ mpp_global_field2D_ug_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field_ug.h>

#ifndef no_8byte_integers
#undef MPP_GLOBAL_FIELD_UG_2D_
#define MPP_GLOBAL_FIELD_UG_2D_ mpp_global_field2D_ug_i8_2d
#undef MPP_GLOBAL_FIELD_UG_3D_
#define MPP_GLOBAL_FIELD_UG_3D_ mpp_global_field2D_ug_i8_3d
#undef MPP_GLOBAL_FIELD_UG_4D_
#define MPP_GLOBAL_FIELD_UG_4D_ mpp_global_field2D_ug_i8_4d
#undef MPP_GLOBAL_FIELD_UG_5D_
#define MPP_GLOBAL_FIELD_UG_5D_ mpp_global_field2D_ug_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field_ug.h>
#endif

#ifdef OVERLOAD_R4
#undef MPP_GLOBAL_FIELD_UG_2D_
#define MPP_GLOBAL_FIELD_UG_2D_ mpp_global_field2D_ug_r4_2d
#undef MPP_GLOBAL_FIELD_UG_3D_
#define MPP_GLOBAL_FIELD_UG_3D_ mpp_global_field2D_ug_r4_3d
#undef MPP_GLOBAL_FIELD_UG_4D_
#define MPP_GLOBAL_FIELD_UG_4D_ mpp_global_field2D_ug_r4_4d
#undef MPP_GLOBAL_FIELD_UG_5D_
#define MPP_GLOBAL_FIELD_UG_5D_ mpp_global_field2D_ug_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field_ug.h>
#endif

#undef MPP_GLOBAL_FIELD_UG_2D_
#define MPP_GLOBAL_FIELD_UG_2D_ mpp_global_field2D_ug_i4_2d
#undef MPP_GLOBAL_FIELD_UG_3D_
#define MPP_GLOBAL_FIELD_UG_3D_ mpp_global_field2D_ug_i4_3d
#undef MPP_GLOBAL_FIELD_UG_4D_
#define MPP_GLOBAL_FIELD_UG_4D_ mpp_global_field2D_ug_i4_4d
#undef MPP_GLOBAL_FIELD_UG_5D_
#define MPP_GLOBAL_FIELD_UG_5D_ mpp_global_field2D_ug_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field_ug.h>


