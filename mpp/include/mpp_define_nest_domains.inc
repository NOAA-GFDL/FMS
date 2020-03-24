! -*-f90-*-


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

!#############################################################################
! Currently the contact will be limited to overlap contact.
subroutine mpp_define_nest_domains(nest_domain, domain, num_nest, nest_level, tile_fine, tile_coarse, &
                                  istart_coarse, icount_coarse, jstart_coarse, jcount_coarse, npes_nest_tile, &
                                  x_refine, y_refine, extra_halo, name)
  type(nest_domain_type),     intent(inout) :: nest_domain
  type(domain2D), target,     intent(in   ) :: domain
  integer,                    intent(in   ) :: num_nest
  integer,                    intent(in   ) :: nest_level(:)
  integer,                    intent(in   ) :: tile_fine(:), tile_coarse(:)
  integer,                    intent(in   ) :: istart_coarse(:), icount_coarse(:), jstart_coarse(:), jcount_coarse(:)
  integer,                    intent(in   ) :: npes_nest_tile(:)
  integer,                    intent(in   ) :: x_refine(:), y_refine(:)
  integer,          optional, intent(in   ) :: extra_halo
  character(len=*), optional, intent(in   ) :: name

  logical                                   :: concurrent
  integer                                   :: n, l, m, my_tile_coarse
  integer                                   :: nx_coarse, ny_coarse
  integer                                   :: nx_fine, ny_fine
  integer                                   :: npes, npes_level, prev_tile_coarse
  integer                                   :: extra_halo_local, npes_nest_top
  integer, dimension(:), allocatable        :: pes, pe_start_pos, pe_end_pos, pelist_level
  logical, dimension(:), allocatable        :: is_nest_fine, is_nest_coarse
  integer, dimension(num_nest)              :: istart_fine, iend_fine, jstart_fine, jend_fine
  integer, dimension(num_nest)              :: iend_coarse, jend_coarse
  integer                                   :: nnest, nlevels, ntiles_top, ntiles, pos
  logical                                   :: is_first

  if(PRESENT(name)) then
     if(len_trim(name) > NAME_LENGTH) then
        call mpp_error(FATAL, "mpp_domains_define.inc(mpp_define_nest_domain): "// &
             "the len_trim of optional argument name ="//trim(name)// &
             " is greater than NAME_LENGTH, change the argument name or increase NAME_LENGTH")
     endif
     nest_domain%name = name
  endif

  extra_halo_local = 0
  if(present(extra_halo)) then
     if(extra_halo .NE. 0) call mpp_error(FATAL, "mpp_define_nest_domains.inc: only support extra_halo=0, contact developer")
     extra_halo_local = extra_halo
  endif

  !---make sure dimension size is correct
  if(size(tile_fine(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(tile_fine) .NE. num_nest")
  if(size(tile_coarse(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(tile_coarse) .NE. num_nest")
  if(size(istart_coarse(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(istart_coarse) .NE. num_nest")
  if(size(icount_coarse(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(icount_coarse) .NE. num_nest")
  if(size(jstart_coarse(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(jstart_coarse) .NE. num_nest")
  if(size(jcount_coarse(:)) .NE. num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc: size(jcount_coarse) .NE. num_nest")

  do n = 1, num_nest
     if(istart_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: istart_coarse < 1")
     if(icount_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: iend_coarse < 1")
     if(jstart_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: jstart_coarse < 1")
     if(jcount_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: jend_coarse < 1")
     iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
     jend_coarse(n) = jstart_coarse(n) + jcount_coarse(n) - 1
     istart_fine(n) = 1 ; iend_fine(n) = icount_coarse(n)*x_refine(n)
     jstart_fine(n) = 1 ; jend_fine(n) = jcount_coarse(n)*y_refine(n)
  end do

  !--- make sure the nest level is monotonic and no jumping
  if(nest_level(1) .NE. 1)  call mpp_error(FATAL, "mpp_define_nest_domains.inc: nest_level(1) .NE. 1")
  do n = 2, num_nest
     if(nest_level(n) < nest_level(n-1)) call mpp_error(FATAL, "mpp_define_nest_domains.inc: nest_level is not monotone increasing")
     if(nest_level(n) > nest_level(n-1)+1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: nest_level(n) > nest_level(n-1)+1")
  enddo
  nlevels = nest_level(num_nest)

  !---make sure tile_fine and tile_nest are monotone increasing.
  do n = 2, num_nest
     if(tile_fine(n) < tile_fine(n-1)) call  mpp_error(FATAL, "mpp_define_nest_domains.inc: tile_fine is not monotone increasing")
     if(tile_coarse(n) < tile_coarse(n-1)) call  mpp_error(FATAL, "mpp_define_nest_domains.inc: tile_coarse is not monotone increasing")
  enddo

  allocate( pes(mpp_npes()) )
  call mpp_get_current_pelist(pes)

  nest_domain%num_nest = num_nest
  allocate(nest_domain%tile_fine(num_nest), nest_domain%tile_coarse(num_nest) )
  allocate(nest_domain%istart_fine(num_nest), nest_domain%iend_fine(num_nest) )
  allocate(nest_domain%jstart_fine(num_nest), nest_domain%jend_fine(num_nest) )
  allocate(nest_domain%istart_coarse(num_nest), nest_domain%iend_coarse(num_nest) )
  allocate(nest_domain%jstart_coarse(num_nest), nest_domain%jend_coarse(num_nest) )

  nest_domain%tile_fine = tile_fine(1:num_nest)
  nest_domain%tile_coarse = tile_coarse(1:num_nest)
  nest_domain%istart_fine = istart_fine(1:num_nest)
  nest_domain%iend_fine = iend_fine(1:num_nest)
  nest_domain%jstart_fine = jstart_fine(1:num_nest)
  nest_domain%jend_fine = jend_fine(1:num_nest)
  nest_domain%istart_coarse = istart_coarse(1:num_nest)
  nest_domain%iend_coarse = iend_coarse(1:num_nest)
  nest_domain%jstart_coarse = jstart_coarse(1:num_nest)
  nest_domain%jend_coarse = jend_coarse(1:num_nest)

  !--- make sure the tile_id of top level of grid is continuous and starting from 1
  if(mpp_pe()==mpp_root_pe()) then
     ntiles_top = domain%ntiles
     do n = 1, ntiles_top
        if(domain%tile_id_all(n) .NE. n) call mpp_error(FATAL, &
            "mpp_define_nest_domains.inc: top level grid tile_id should be 1, 2, ..,ntiles")
     enddo
  endif
  call mpp_broadcast(ntiles_top, mpp_root_pe())
  !--- make sure the nest grid tile_ids are continuous
  do n = 1, num_nest
     if(tile_fine(n) .NE. ntiles_top+n) then
        print*, "tile_fine, ntile_top, n=", tile_fine(n), ntiles_top, n, mpp_pe()
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: tile_id is not continuous")
     endif
  enddo

  allocate(pe_start_pos(ntiles_top+num_nest))
  allocate(pe_end_pos(ntiles_top+num_nest))

  do n = 2, ntiles_top
     if(npes_nest_tile(n) .NE. npes_nest_tile(n-1)) call mpp_error(FATAL, &
         "mpp_define_nest_domains.inc: all the tiles in top grid should use same number of MPI ranks")
  enddo

  npes_nest_top = npes_nest_tile(1)*ntiles_top

  !--- get the pe start and end pos for each tile
  do n = 1, ntiles_top
     pe_start_pos(n) = 1
     Pe_end_pos(n) = npes_nest_tile(1)*ntiles_top
  enddo
  ntiles = tile_fine(num_nest)
  if(ntiles .NE. ntiles_top + num_nest) call mpp_error(FATAL, "mpp_define_nest_domains.inc:ntiles .NE. ntiles_top + num_nest")
  do n = 1, num_nest
     pe_start_pos(ntiles_top+n) = pe_end_pos(ntiles_top+n-1) + 1
     pe_end_pos(ntiles_top+n)   = pe_end_pos(ntiles_top+n-1) + npes_nest_tile(tile_fine(n))
  enddo

  nest_domain%num_level = nlevels
  allocate(nest_domain%nest(nlevels))
  allocate(pelist_level(mpp_npes()))
  allocate(is_nest_fine(nlevels))
  allocate(is_nest_coarse(nlevels))

  !--- setup pelist for each level
  pos = 0
  is_nest_fine(:) = .false.
  is_nest_coarse(:) = .false.
  do l = 1, nlevels
     npes_level = 0
     pos = 0
     is_first = .true.
     prev_tile_coarse = 0
     !--- first get coarse processor
     do n = 1, num_nest
        if(nest_level(n) == l) then
           if(mpp_pe() .GE. pes(pe_start_pos(tile_fine(n))) .AND. mpp_pe() .LE. pes(pe_end_pos(tile_fine(n)))) then
              is_nest_fine(l) = .true.
           endif
           if(mpp_pe() .GE. pes(pe_start_pos(tile_coarse(n))) .AND. mpp_pe() .LE. pes(pe_end_pos(tile_coarse(n)))) then
              is_nest_coarse(l) = .true.
           endif
           if(pos==0 .OR. (l .NE. 1 .AND. prev_tile_coarse .NE. tile_coarse(n)) ) then
              do m = pe_start_pos(tile_coarse(n)), pe_end_pos(tile_coarse(n))
                 pos = pos+1
                 pelist_level(pos) = pes(m)
              enddo
              npes_level = npes_level + pe_end_pos(tile_coarse(n)) - pe_start_pos(tile_coarse(n)) + 1
           endif
           prev_tile_coarse = tile_coarse(n)
        endif
     enddo
     ! fine processor
     do n = 1, num_nest
        if(nest_level(n) == l) then
           do m = pe_start_pos(tile_fine(n)), pe_end_pos(tile_fine(n))
              pos = pos+1
              pelist_level(pos) = pes(m)
           enddo
           npes_level = npes_level + pe_end_pos(tile_fine(n)) - pe_start_pos(tile_fine(n)) + 1
        endif
     enddo

     allocate(nest_domain%nest(l)%pelist(npes_level))
     nest_domain%nest(l)%pelist(:) = pelist_level(1:npes_level)

     call mpp_declare_pelist(nest_domain%nest(l)%pelist)
     nest_domain%nest(l)%on_level = ANY(nest_domain%nest(l)%pelist(:)==mpp_pe())
     nest_domain%nest(l)%is_fine_pe = is_nest_fine(l)
     nest_domain%nest(l)%is_coarse_pe = is_nest_coarse(l)
     if(nest_domain%nest(l)%on_level .neqv. (is_nest_fine(l) .OR. is_nest_coarse(l))) then
        print*, "on_level=", nest_domain%nest(l)%on_level, is_nest_fine(l), is_nest_coarse(l), mpp_pe(),l
        call mpp_error(FATAL, "mpp_define_nest_domains.inc:on_level does not match is_nest_fine/is_nest_coarse")
     endif
     if(is_nest_fine(l) .and. is_nest_coarse(l)) then
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: is_nest_fine and is_nest_coarse can not both be true")
     endif
  enddo

  if(count(is_nest_fine)>1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: count(is_nest_fine)>1")
  if(count(is_nest_coarse)>1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: count(is_nest_coarse)>1")

  do l = 1, nlevels
     !--- setup for each level
     if(nest_domain%nest(l)%on_level) then
        call mpp_set_current_pelist(nest_domain%nest(l)%pelist)
        nnest = count(nest_level==l)
        nest_domain%nest(l)%num_nest = nnest
        allocate(nest_domain%nest(l)%tile_fine(nnest), nest_domain%nest(l)%tile_coarse(nnest) )
        allocate(nest_domain%nest(l)%istart_fine(nnest), nest_domain%nest(l)%iend_fine(nnest) )
        allocate(nest_domain%nest(l)%jstart_fine(nnest), nest_domain%nest(l)%jend_fine(nnest) )
        allocate(nest_domain%nest(l)%istart_coarse(nnest), nest_domain%nest(l)%iend_coarse(nnest) )
        allocate(nest_domain%nest(l)%jstart_coarse(nnest), nest_domain%nest(l)%jend_coarse(nnest) )
        my_tile_coarse = 0

        pos=0
        do n = 1, num_nest
           if(nest_level(n) ==l) then
              pos = pos+1
              nest_domain%nest(l)%tile_fine(pos) = tile_fine(n)
              nest_domain%nest(l)%tile_coarse(pos) = tile_coarse(n)
              nest_domain%nest(l)%istart_fine(pos) = istart_fine(n)
              nest_domain%nest(l)%iend_fine(pos) = iend_fine(n)
              nest_domain%nest(l)%jstart_fine(pos) = jstart_fine(n)
              nest_domain%nest(l)%jend_fine(pos) = jend_fine(n)
              nest_domain%nest(l)%istart_coarse(pos) = istart_coarse(n)
              nest_domain%nest(l)%iend_coarse(pos) = iend_coarse(n)
              nest_domain%nest(l)%jstart_coarse(pos) = jstart_coarse(n)
              nest_domain%nest(l)%jend_coarse(pos) = jend_coarse(n)
              if(l==1) then
                 my_tile_coarse = 1
              else if((mpp_pe() .GE. pes(pe_start_pos(tile_fine(n))) .AND. mpp_pe() .LE. pes(pe_end_pos(tile_fine(n)))) .OR.  &
                 (mpp_pe() .GE. pes(pe_start_pos(tile_coarse(n))) .AND. mpp_pe() .LE. pes(pe_end_pos(tile_coarse(n)))) ) then
                 my_tile_coarse = tile_coarse(n)
              endif
           endif
        enddo
        if(my_tile_coarse == 0) call mpp_error(FATAL, "mpp_define_nest_domains.inc: my_tile_coarse == 0")

        if(pos .NE. nest_domain%nest(l)%num_nest)  &
            call mpp_error(FATAL, "mpp_define_nest_domains.inc:pos .NE. nest_domain%nest(l)%num_nest")

        if(is_nest_fine(l)) then
           nest_domain%nest(l)%domain_fine=>domain
           allocate(nest_domain%nest(l)%domain_coarse)
        else if(is_nest_coarse(l)) then
           nest_domain%nest(l)%domain_coarse=>domain
           allocate(nest_domain%nest(l)%domain_fine)
        endif
!!!! DEBUG CODE ! has problems on coarse domain
!!$        print*, 'MPP_BROADCAST_DOMAIN: ', mpp_pe(), l, & !ASSOCIATED(nest_domain%nest(l)%domain_fine), &
!!$             nest_domain%nest(l)%domain_fine%tile_id(1), nest_domain%nest(l)%tile_fine, tile_fine
!!!! END DEBUG CODE
        call mpp_broadcast_domain(nest_domain%nest(l)%domain_fine, nest_domain%nest(l)%tile_fine)
        call mpp_broadcast_domain(nest_domain%nest(l)%domain_coarse, my_tile_coarse)
        call define_nest_level_type(nest_domain%nest(l), x_refine(l), y_refine(l), extra_halo_local)
     endif
  enddo

end subroutine mpp_define_nest_domains

subroutine define_nest_level_type(nest_domain, x_refine, y_refine, extra_halo)
  type(nest_level_type), intent(inout) :: nest_domain
  integer,               intent(in   ) :: extra_halo
  integer,               intent(in   ) :: x_refine, y_refine

  integer                                   :: n
  integer                                   :: nx_coarse, ny_coarse
  integer                                   :: nx_fine, ny_fine
  integer                                   :: npes, npes_fine, npes_coarse
  integer                                   :: extra_halo_local
  integer, allocatable                      :: pes_coarse(:)
  integer, allocatable                      :: pes_fine(:)
  integer, dimension(nest_domain%num_nest)  :: my_nest_id
  integer                                   :: my_num_nest, nnest

  npes = size(nest_domain%pelist(:))
  npes_coarse = size(nest_domain%domain_coarse%list(:))
  npes_fine = size(nest_domain%domain_fine%list(:))
  !--- pes_fine and pes_coarse should be subset of pelist
  allocate( pes_coarse(npes_coarse) )
  allocate( pes_fine  (npes_fine  ) )
  do n = 1, npes_coarse
     pes_coarse(n) = nest_domain%domain_coarse%list(n-1)%pe
     if( .NOT. ANY(nest_domain%pelist(:) == pes_coarse(n)) ) then
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: pelist_coarse is not subset of pelist")
     endif
  enddo
  do n = 1, npes_fine
     pes_fine(n) = nest_domain%domain_fine%list(n-1)%pe
     if( .NOT. ANY(nest_domain%pelist(:) == pes_fine(n)) ) then
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: pelist_fine is not subset of pelist")
     endif
  enddo

  allocate(nest_domain%pelist_fine(npes_fine))
  allocate(nest_domain%pelist_coarse(npes_coarse))
  nest_domain%pelist_fine = pes_fine
  nest_domain%pelist_coarse = pes_coarse
  if( nest_domain%is_fine_pe .neqv. ANY(pes_fine(:) == mpp_pe()) ) then
     call mpp_error(FATAL, "mpp_define_nest_domains.inc: nest_domain%is_fine_pe .neqv. ANY(pes_fine(:) == mpp_pe())")
  endif
  if( nest_domain%is_coarse_pe .neqv. ANY(pes_coarse(:) == mpp_pe()) ) then
     call mpp_error(FATAL, "mpp_define_nest_domains.inc: nest_domain%is_coarse_pe .neqv. ANY(pes_coarse(:) == mpp_pe())")
  endif

  !--- figure out the the corresponding nested region.
  !--- on coarse grid pe, it might overlap multiple fine regon.
  !--- on fine grid pe, it always only overlap at most 1 coarse region.
  my_num_nest= 0
  my_nest_id(:) = 0
  if( nest_domain%is_fine_pe ) then
     !--- figure out the nest number on current pe
     do n = 1, nest_domain%num_nest
        if(nest_domain%domain_fine%tile_id(1) == nest_domain%tile_fine(n)) then
           my_num_nest = my_num_nest + 1
           my_nest_id(my_num_nest) = n
           exit
        end if
     end do
     if(my_num_nest .NE. 1) then
        print*, "num_nest=", my_num_nest, nest_domain%domain_fine%tile_id(1), nest_domain%tile_fine(1)
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: my_num_nest .ne. 1 on fine pelist")
     endif
  else if( nest_domain%is_coarse_pe ) then
     my_num_nest = nest_domain%num_nest
     do n = 1, nest_domain%num_nest
        my_nest_id(n) = n
     enddo
  endif

  nest_domain%my_num_nest = my_num_nest
  if(my_num_nest>0) then
     allocate(nest_domain%my_nest_id(my_num_nest))
     nest_domain%my_nest_id(:) = my_nest_id(1:my_num_nest)
  endif

  !--- We are assuming the fine grid is fully overlapped with coarse grid.
  if( nest_domain%is_fine_pe ) then
     if( nest_domain%iend_fine(my_nest_id(1))-nest_domain%istart_fine(my_nest_id(1))+1 &
         .NE. nest_domain%domain_fine%x(1)%global%size .OR.  &
         nest_domain%jend_fine(my_nest_id(1))-nest_domain%jstart_fine(my_nest_id(1))+1 &
         .NE. nest_domain%domain_fine%y(1)%global%size ) then
        print*, "x size are", nest_domain%domain_fine%x(1)%global%size, &
                   nest_domain%istart_fine(my_nest_id(1)), nest_domain%iend_fine(my_nest_id(1))
        print*, "y size are", nest_domain%domain_fine%y(1)%global%size, &
                   nest_domain%jstart_fine(my_nest_id(1)), nest_domain%jend_fine(my_nest_id(1))
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: The fine global domain is not covered by coarse domain")
     endif
  endif

  ! only support concurrent run for fine and coarse domain, currently only check on coarse pe
  if(nest_domain%is_coarse_pe) then
!     if( npes_fine + npes_coarse .NE. npes ) then
!        print*, "On pe =", mpp_pe(), npes_fine, npes_coarse, npes
!        call mpp_error(FATAL, "mpp_domains_define.inc: size(pelist_coarse)+size(pelist_fine) .NE. size(pelist)")
!     endif
  endif

  !--- coarse grid and fine grid should be both symmetry or non-symmetry.
  if(nest_domain%domain_coarse%symmetry .neqv. nest_domain%domain_fine%symmetry) then
     print*,"symmetry is", nest_domain%domain_coarse%symmetry, nest_domain%domain_fine%symmetry, mpp_pe()
     call mpp_error(FATAL, "mpp_domains_define.inc: domain_coarse%symmetry .neqv. .NOT. domain_fine%symmetry")
  endif

  nest_domain%x_refine = x_refine
  nest_domain%y_refine = y_refine

  allocate( nest_domain%C2F_T, nest_domain%C2F_C, nest_domain%C2F_E, nest_domain%C2F_N )
  nest_domain%C2F_T%next => NULL()
  nest_domain%C2F_C%next => NULL()
  nest_domain%C2F_N%next => NULL()
  nest_domain%C2F_E%next => NULL()
  allocate( nest_domain%F2C_T, nest_domain%F2C_C, nest_domain%F2C_E, nest_domain%F2C_N )

  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_T, CENTER, "F2C T-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_E, EAST,   "F2C E-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_C, CORNER, "F2C C-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_N, NORTH,  "F2C N-cell")

  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_T, extra_halo, CENTER, "C2F T-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_E, extra_halo, EAST,   "C2F E-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_C, extra_halo, CORNER, "C2F C-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_N, extra_halo, NORTH,  "C2F N-cell")

  deallocate(pes_fine, pes_coarse)


end subroutine define_nest_level_type


!###############################################################################
subroutine compute_overlap_coarse_to_fine(nest_domain, overlap, extra_halo, position, name)
  type(nest_level_type),  intent(inout) :: nest_domain
  type(nestSpec),         intent(inout) :: overlap
  integer,                intent(in   ) :: extra_halo
  integer,                intent(in   ) :: position
  character(len=*),       intent(in   ) :: name

  type(domain2D), pointer :: domain_fine  =>NULL()
  type(domain2D), pointer :: domain_coarse=>NULL()
  type(overlap_type), allocatable :: overlapList(:)
  logical              :: is_first
  integer              :: tile_fine, tile_coarse
  integer              :: istart_fine, iend_fine, jstart_fine, jend_fine
  integer              :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
  integer              :: whalo, ehalo, shalo, nhalo
  integer              :: npes, npes_fine, npes_coarse, n, m
  integer              :: isg_fine, ieg_fine, jsg_fine, jeg_fine
  integer              :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
  integer              :: is_coarse, ie_coarse, js_coarse, je_coarse
  integer              :: is_coarse2, ie_coarse2, js_coarse2, je_coarse2
  integer              :: is_convert, ie_convert, js_convert, je_convert, rotate
  integer              :: is_convert2(2), ie_convert2(2), js_convert2(2), je_convert2(2), rotate2(2)
  integer              :: isc_fine, iec_fine, jsc_fine, jec_fine
  integer              :: isd_fine, ied_fine, jsd_fine, jed_fine
  integer              :: isc_east, iec_east, jsc_east, jec_east
  integer              :: isc_west, iec_west, jsc_west, jec_west
  integer              :: isc_south, iec_south, jsc_south, jec_south
  integer              :: isc_north, iec_north, jsc_north, jec_north
  integer              :: x_refine, y_refine, ishift, jshift
  integer              :: nsend, nrecv, dir, from_pe, l, nn
  integer              :: is, ie, js, je, msgsize, nconvert
  integer, allocatable :: msg1(:), msg2(:)
  integer, allocatable :: isl_coarse(:), iel_coarse(:), jsl_coarse(:), jel_coarse(:)
  integer, allocatable :: isl_fine(:), iel_fine(:), jsl_fine(:), jel_fine(:)
  integer, allocatable :: isgl_fine(:), iegl_fine(:), jsgl_fine(:), jegl_fine(:)
  integer              :: outunit


  outunit = stdout()
  domain_fine   => nest_domain%domain_fine
  domain_coarse => nest_domain%domain_coarse
  call mpp_get_domain_shift   (domain_coarse, ishift, jshift, position)
  npes          = mpp_npes()
  npes_fine     = size(nest_domain%pelist_fine(:))
  npes_coarse   = size(nest_domain%pelist_coarse(:))

  allocate(isl_coarse(npes_coarse), iel_coarse(npes_coarse))
  allocate(jsl_coarse(npes_coarse), jel_coarse(npes_coarse))
  allocate(isl_fine  (npes_fine  ), iel_fine  (npes_fine  ))
  allocate(jsl_fine  (npes_fine  ), jel_fine  (npes_fine  ))
  allocate(isgl_fine  (npes_fine  ), iegl_fine  (npes_fine  ))
  allocate(jsgl_fine  (npes_fine  ), jegl_fine  (npes_fine  ))

  call mpp_get_global_domain  (domain_fine,   xbegin=isg_fine,   xend=ieg_fine,   &
       ybegin=jsg_fine,   yend=jeg_fine)
  call mpp_get_compute_domain (domain_coarse, xbegin=isc_coarse, xend=iec_coarse, &
       ybegin=jsc_coarse, yend=jec_coarse)
  call mpp_get_compute_domain (domain_fine,   xbegin=isc_fine,   xend=iec_fine,   &
       ybegin=jsc_fine,   yend=jec_fine)
  call mpp_get_compute_domains(domain_coarse, xbegin=isl_coarse, xend=iel_coarse, &
       ybegin=jsl_coarse, yend=jel_coarse)
  call mpp_get_compute_domains(domain_fine,   xbegin=isl_fine,   xend=iel_fine,   &
       ybegin=jsl_fine,   yend=jel_fine)
  call mpp_get_global_domains(domain_fine,   xbegin=isgl_fine,   xend=iegl_fine,   &
       ybegin=jsgl_fine,   yend=jegl_fine)

  if( nest_domain%is_coarse_pe ) then
     allocate(overLaplist(npes_fine))
     overlap%xbegin = isc_coarse - domain_coarse%whalo
     overlap%xend   = iec_coarse + domain_coarse%ehalo + ishift
     overlap%ybegin = jsc_coarse - domain_coarse%shalo
     overlap%yend   = jec_coarse + domain_coarse%nhalo + jshift
  else
     allocate(overLaplist(npes_coarse))
     overlap%xbegin = isc_fine - domain_fine%whalo
     overlap%xend   = iec_fine + domain_fine%ehalo + ishift
     overlap%ybegin = jsc_fine - domain_fine%shalo
     overlap%yend   = jec_fine + domain_fine%nhalo + jshift
  endif

  overlap%extra_halo = extra_halo
  x_refine      = nest_domain%x_refine
  y_refine      = nest_domain%y_refine
  whalo         = domain_fine%whalo + extra_halo
  ehalo         = domain_fine%ehalo + extra_halo
  shalo         = domain_fine%shalo + extra_halo
  nhalo         = domain_fine%nhalo + extra_halo

  isd_fine = isc_fine - whalo
  ied_fine = iec_fine + ehalo
  jsd_fine = jsc_fine - shalo
  jed_fine = jec_fine + nhalo

  overlap%nsend = 0
  overlap%nrecv = 0
  call init_index_type(overlap%west)
  call init_index_type(overlap%east)
  call init_index_type(overlap%south)
  call init_index_type(overlap%north)
  nsend = 0
  nrecv = 0

  do nn = 1, nest_domain%num_nest

     tile_fine     = nest_domain%tile_fine(nn)
     tile_coarse   = nest_domain%tile_coarse(nn)
     istart_fine   = nest_domain%istart_fine(nn)
     iend_fine     = nest_domain%iend_fine(nn)
     jstart_fine   = nest_domain%jstart_fine(nn)
     jend_fine     = nest_domain%jend_fine(nn)
     istart_coarse = nest_domain%istart_coarse(nn)
     iend_coarse   = nest_domain%iend_coarse(nn)
     jstart_coarse = nest_domain%jstart_coarse(nn)
     jend_coarse   = nest_domain%jend_coarse(nn)

     !--- first compute the halo region and corresponding index in coarse grid.
     if( nest_domain%is_fine_pe .and.  domain_fine%tile_id(1) == tile_fine) then
        if( ieg_fine == iec_fine ) then   ! east halo
           is_coarse = iend_coarse
           ie_coarse = iend_coarse + ehalo
           js_coarse = jstart_coarse + ( jsc_fine - jsg_fine )/y_refine
           je_coarse = jstart_coarse + ( jec_fine - jsg_fine )/y_refine
           js_coarse = js_coarse - shalo
           je_coarse = je_coarse + nhalo

           overlap%east%is_me  = iec_fine + 1
           overlap%east%ie_me  = ied_fine
           overlap%east%js_me  = jsd_fine
           overlap%east%je_me  = jed_fine
           overlap%east%is_you = is_coarse
           overlap%east%ie_you = ie_coarse
           overlap%east%js_you = js_coarse
           overlap%east%je_you = je_coarse
        endif

        if( jsg_fine == jsc_fine ) then  ! south
           is_coarse = istart_coarse + ( isc_fine - isg_fine )/x_refine
           ie_coarse = istart_coarse + ( iec_fine - isg_fine )/x_refine
           is_coarse = is_coarse - whalo
           ie_coarse = ie_coarse + ehalo
           js_coarse = jstart_coarse - shalo
           je_coarse = jstart_coarse
           overlap%south%is_me  = isd_fine
           overlap%south%ie_me  = ied_fine
           overlap%south%js_me  = jsd_fine
           overlap%south%je_me  = jsc_fine-1
           overlap%south%is_you = is_coarse
           overlap%south%ie_you = ie_coarse
           overlap%south%js_you = js_coarse
           overlap%south%je_you = je_coarse
        endif

        if( isg_fine == isc_fine ) then ! west
           is_coarse = istart_coarse - whalo
           ie_coarse = istart_coarse
           js_coarse = jstart_coarse + ( jsc_fine - jsg_fine )/y_refine
           je_coarse = jstart_coarse + ( jec_fine - jsg_fine )/y_refine
           js_coarse = js_coarse - shalo
           je_coarse = je_coarse + nhalo
           overlap%west%is_me = isd_fine
           overlap%west%ie_me = isc_fine-1
           overlap%west%js_me = jsd_fine
           overlap%west%je_me = jed_fine
           overlap%west%is_you = is_coarse
           overlap%west%ie_you = ie_coarse
           overlap%west%js_you = js_coarse
           overlap%west%je_you = je_coarse
        endif

        if( jeg_fine == jec_fine ) then ! north
           is_coarse = istart_coarse + ( isc_fine - isg_fine )/x_refine
           ie_coarse = istart_coarse + ( iec_fine - isg_fine )/x_refine
           is_coarse = is_coarse - whalo
           ie_coarse = ie_coarse + ehalo
           js_coarse = jend_coarse
           je_coarse = jend_coarse + nhalo
           overlap%north%is_me = isd_fine
           overlap%north%ie_me = ied_fine
           overlap%north%js_me = jec_fine+1
           overlap%north%je_me = jed_fine
           overlap%north%is_you = is_coarse
           overlap%north%ie_you = ie_coarse
           overlap%north%js_you = js_coarse
           overlap%north%je_you = je_coarse
        endif

        !-------------------------------------------------------------------------
        !
        !                 Receiving
        !
        !-------------------------------------------------------------------------
        !--- loop through coarse pelist
        do n = 1, npes_coarse
           is_first = .true.
           do m = 1, 4
              select case (m)
              case (1)          !--- east halo receiving
                 dir = 1
                 is_coarse = overlap%east%is_you
                 ie_coarse = overlap%east%ie_you
                 js_coarse = overlap%east%js_you
                 je_coarse = overlap%east%je_you
              case (2)          !--- south halo receiving
                 dir = 3
                 is_coarse = overlap%south%is_you
                 ie_coarse = overlap%south%ie_you
                 js_coarse = overlap%south%js_you
                 je_coarse = overlap%south%je_you
              case (3)          !--- west halo receiving
                 dir = 5
                 is_coarse = overlap%west%is_you
                 ie_coarse = overlap%west%ie_you
                 js_coarse = overlap%west%js_you
                 je_coarse = overlap%west%je_you
              case (4)          !--- north halo receiving
                 dir = 7
                 is_coarse = overlap%north%is_you
                 ie_coarse = overlap%north%ie_you
                 js_coarse = overlap%north%js_you
                 je_coarse = overlap%north%je_you
              end select
              if( je_coarse .GE. js_coarse .AND. ie_coarse .GE. is_coarse ) then
                 ! convert coarse grid index to the nested grid coarse grid index.
                 nconvert = convert_index_to_nest(domain_coarse, 0, 0, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%list(n-1)%tile_id(1), &
                      isl_coarse(n), iel_coarse(n), jsl_coarse(n), jel_coarse(n), &
                      is_convert2, ie_convert2, js_convert2, je_convert2, rotate2)
                 do l = 1, nconvert
                    is_coarse2 = max( is_coarse, is_convert2(l) )
                    ie_coarse2 = min( ie_coarse, ie_convert2(l) )
                    js_coarse2 = max( js_coarse, js_convert2(l) )
                    je_coarse2 = min( je_coarse, je_convert2(l) )
                    if( ie_coarse2 .GE. is_coarse2 .AND. je_coarse2 .GE. js_coarse2 ) then
                       select case (m)
                       case (1)          !--- east halo
                          is_coarse2 = is_coarse2+ishift
                          ie_coarse2 = ie_coarse2+ishift
                          if(je_coarse2 == overlap%east%je_you) je_coarse2 = je_coarse2+jshift
                       case (2)          !--- south halo
                          if(ie_coarse2 == overlap%south%ie_you) ie_coarse2 = ie_coarse2+ishift
                       case (3)          !--- west halo
                          if(je_coarse2 == overlap%west%je_you) je_coarse2 = je_coarse2+jshift
                       case (4)          !--- north halo
                          if(ie_coarse2 == overlap%north%ie_you) ie_coarse2 = ie_coarse2+ishift
                          js_coarse2 = js_coarse2+jshift
                          je_coarse2 = je_coarse2+jshift
                       end select

                       if(is_first) then
                          nrecv = nrecv + 1
                          call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                          is_first = .false.
                       endif
                       rotate = -rotate2(l)
                       call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_coarse(n), &
                            is_coarse2, ie_coarse2, js_coarse2, je_coarse2 , dir,  rotate2(l))
                    endif
                 enddo
              endif
           enddo
        enddo

     endif
     !-----------------------------------------------------------------------
     !
     !                          Sending
     !
     !-----------------------------------------------------------------------

     if( nest_domain%is_coarse_pe ) then
        do n = 1, npes_fine
           if( domain_fine%list(n-1)%tile_id(1) .NE. tile_fine ) cycle
           is_first = .true.
           isg_fine = isgl_fine(n)
           ieg_fine = iegl_fine(n)
           jsg_fine = jsgl_fine(n)
           jeg_fine = jegl_fine(n)

           !--- to_pe's east
           if( ieg_fine == iel_fine(n) ) then
              dir = 1
              is_coarse = iend_coarse
              ie_coarse = iend_coarse + ehalo
              js_coarse = jstart_coarse + ( jsl_fine(n) - jsg_fine )/y_refine
              je_coarse = jstart_coarse + ( jel_fine(n) - jsg_fine )/y_refine
              js_coarse = js_coarse - shalo
              je_coarse = je_coarse + nhalo
              !--- convert the index to coarse grid index.
              nconvert = convert_index_to_coarse(domain_coarse, 0, 0, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%tile_id(1), is_coarse, ie_coarse, &
                   js_coarse, je_coarse, is_convert2, ie_convert2, &
                   js_convert2, je_convert2, rotate2)
              do l = 1, nconvert
                 is_coarse = max(isc_coarse, is_convert2(l))
                 ie_coarse = min(iec_coarse, ie_convert2(l))
                 js_coarse = max(jsc_coarse, js_convert2(l))
                 je_coarse = min(jec_coarse, je_convert2(l))
                 if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                    if(rotate2(l)==ZERO) then
                       is_coarse = is_coarse+ishift
                       ie_coarse = ie_coarse+ishift
                       if( je_coarse == je_convert2(l) ) je_coarse = je_coarse+jshift
                    else if(rotate2(l) == MINUS_NINETY) then
                       js_coarse = js_coarse+ishift
                       je_coarse = je_coarse+ishift
                       if(is_coarse==is_convert2(l)) is_coarse = is_coarse-jshift
                       is_coarse = is_coarse+jshift
                       ie_coarse = ie_coarse+jshift
                    else if(rotate2(l) == NINETY) then
                       if(ie_coarse==ie_convert2(l)) ie_coarse = ie_coarse+jshift
                    endif

                    if(is_first) then
                       nsend = nsend + 1
                       call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                       is_first = .false.
                    endif
                    rotate = -rotate2(l)
                    call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                         is_coarse, ie_coarse, js_coarse, je_coarse , dir, rotate)
                 endif
              enddo
           endif

           !--- to_pe's south
           if( jsg_fine == jsl_fine(n) ) then
              dir = 3
              is_coarse = istart_coarse + ( isl_fine(n) - isg_fine )/x_refine
              ie_coarse = istart_coarse + ( iel_fine(n) - isg_fine )/x_refine
              is_coarse = is_coarse - shalo
              ie_coarse = ie_coarse + nhalo
              js_coarse = jstart_coarse - shalo
              je_coarse = jstart_coarse
              !--- convert the index to coarse grid index.
              nconvert=convert_index_to_coarse(domain_coarse, 0, 0, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%tile_id(1), is_coarse, ie_coarse, &
                   js_coarse, je_coarse, is_convert2, ie_convert2, &
                   js_convert2, je_convert2, rotate2)
              do l = 1, nconvert
                 is_coarse = max(isc_coarse, is_convert2(l))
                 ie_coarse = min(iec_coarse, ie_convert2(l))
                 js_coarse = max(jsc_coarse, js_convert2(l))
                 je_coarse = min(jec_coarse, je_convert2(l))

                 if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                    if(rotate2(l)==ZERO .AND. ie_coarse==ie_convert2(l)) then
                       ie_coarse = ie_coarse+ishift
                    else if( rotate2(l) .NE. ZERO .AND. je_coarse == je_convert2(l) ) then
                       je_coarse = je_coarse+ishift
                    endif
                    if(is_first) then
                       nsend = nsend + 1
                       call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                       is_first = .false.
                    endif
                    rotate = -rotate2(l)
                    call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                         is_coarse, ie_coarse, js_coarse, je_coarse , dir, rotate)
                 endif
              enddo
           endif

           !--- to_pe's west
           if( isg_fine == isl_fine(n) ) then
              dir = 5
              is_coarse = istart_coarse - whalo
              ie_coarse = istart_coarse
              js_coarse = jstart_coarse + ( jsl_fine(n) - jsg_fine )/y_refine
              je_coarse = jstart_coarse + ( jel_fine(n) - jsg_fine )/y_refine
              js_coarse = js_coarse - shalo
              je_coarse = je_coarse + nhalo
              !--- convert the index to coarse grid index.
              nconvert=convert_index_to_coarse(domain_coarse, 0, 0, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%tile_id(1), is_coarse, ie_coarse, &
                   js_coarse, je_coarse, is_convert2, ie_convert2, &
                   js_convert2, je_convert2, rotate2)
              do l = 1, nconvert
                 is_coarse = max(isc_coarse, is_convert2(l))
                 ie_coarse = min(iec_coarse, ie_convert2(l))
                 js_coarse = max(jsc_coarse, js_convert2(l))
                 je_coarse = min(jec_coarse, je_convert2(l))
                 if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                    if(rotate2(l)==ZERO .and. je_coarse == je_convert2(l) ) then
                       je_coarse = je_coarse+jshift
                    else if(rotate2(l) .NE. ZERO .and. ie_coarse == ie_convert2(l) ) then
                       ie_coarse = ie_coarse+jshift
                    endif
                    if(is_first) then
                       nsend = nsend + 1
                       call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                       is_first = .false.
                    endif
                    rotate = -rotate2(l)
                    call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                         is_coarse, ie_coarse, js_coarse, je_coarse , dir, rotate)
                 endif
              enddo
           endif

           !--- to_pe's north
           if( jeg_fine == jel_fine(n) ) then
              dir = 7
              is_coarse = istart_coarse + ( isl_fine(n) - isg_fine )/x_refine
              ie_coarse = istart_coarse + ( iel_fine(n) - isg_fine )/x_refine
              is_coarse = is_coarse - shalo
              ie_coarse = ie_coarse + nhalo
              js_coarse = jend_coarse
              je_coarse = jend_coarse + nhalo
              !--- convert the index to coarse grid index.
              nconvert=convert_index_to_coarse(domain_coarse, 0, 0, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%tile_id(1), is_coarse, ie_coarse, &
                   js_coarse, je_coarse, is_convert2, ie_convert2, &
                   js_convert2, je_convert2, rotate2)
              do l = 1, nconvert
                 is_coarse = max(isc_coarse, is_convert2(l))
                 ie_coarse = min(iec_coarse, ie_convert2(l))
                 js_coarse = max(jsc_coarse, js_convert2(l))
                 je_coarse = min(jec_coarse, je_convert2(l))
                 if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                    if(rotate2(l)==ZERO) then
                       if(ie_coarse==ie_convert2(l)) ie_coarse = ie_coarse+ishift
                       js_coarse = js_coarse+jshift
                       je_coarse = je_coarse+jshift
                    else if(rotate2(l) == NINETY) then
                       if(js_coarse==js_convert2(l)) js_coarse = js_coarse-ishift
                       js_coarse = js_coarse+ishift
                       je_coarse = je_coarse+ishift
                       is_coarse = is_coarse+jshift
                       ie_coarse = ie_coarse+jshift
                    else if(rotate2(l) == MINUS_NINETY ) then
                       if(je_coarse==je_convert2(l)) je_coarse = je_coarse+ishift
                    endif
                    if(is_first) then
                       nsend = nsend + 1
                       call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                       is_first = .false.
                    endif
                    rotate = -rotate2(l)
                    call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                         is_coarse, ie_coarse, js_coarse, je_coarse , dir, rotate)
                 endif
              enddo
           endif
        enddo
     endif
  enddo

  !--- copy the overlapping into nest_domain data.
  overlap%nrecv = nrecv
  if( nrecv > 0 ) then
     allocate(overlap%recv(nrecv))
     do n = 1, nrecv
        call copy_nest_overlap( overlap%recv(n), overLaplist(n) )
!        call print_nest_overlap(overlap%recv(n), "C2F RECV")
        call deallocate_nest_overlap( overLaplist(n) )
     enddo
  endif

  overlap%nsend = nsend
  if( nsend > 0 ) then
     allocate(overlap%send(nsend))
     do n = 1, nsend
        call copy_nest_overlap( overlap%send(n), overLaplist(n) )
!       call print_nest_overlap(overlap%send(n), "C2F SEND")
        call deallocate_nest_overlap( overLaplist(n) )
     enddo
  endif
  if(allocated(overlaplist))deallocate(overLaplist)


  deallocate(isl_coarse, iel_coarse, jsl_coarse, jel_coarse)
  deallocate(isl_fine, iel_fine, jsl_fine, jel_fine)
  deallocate(isgl_fine, iegl_fine, jsgl_fine, jegl_fine)

  !--- add shift value accoring grid position
  if( nest_domain%is_fine_pe ) then
     if( ieg_fine == iec_fine ) then   ! east halo
        overlap%east%is_me  = overlap%east%is_me + ishift
        overlap%east%ie_me  = overlap%east%ie_me + ishift
        overlap%east%je_me  = overlap%east%je_me + jshift
        overlap%east%is_you = overlap%east%is_you + ishift
        overlap%east%ie_you = overlap%east%ie_you + ishift
        overlap%east%je_you = overlap%east%je_you + jshift
     endif

     if( jsg_fine == jsc_fine ) then  ! south
        overlap%south%ie_me  = overlap%south%ie_me + ishift
        overlap%south%ie_you  = overlap%south%ie_you + ishift
     endif

     if( isg_fine == isc_fine ) then ! west
        overlap%west%je_me  = overlap%west%je_me + jshift
        overlap%west%je_you  = overlap%west%je_you + jshift
     endif

     if( jeg_fine == jec_fine ) then ! north
        overlap%north%ie_me = overlap%north%ie_me + ishift
        overlap%north%js_me = overlap%north%js_me + jshift
        overlap%north%je_me = overlap%north%je_me + jshift
        overlap%north%ie_you = overlap%north%ie_you + ishift
        overlap%north%js_you = overlap%north%js_you + jshift
        overlap%north%je_you = overlap%north%je_you + jshift
     endif
  endif

  if(debug_message_passing) call debug_message_size(overlap, name)


end subroutine compute_overlap_coarse_to_fine

!###############################################################################
!-- This routine will compute the send and recv information between overlapped nesting
!-- region. The data is assumed on T-cell center.
subroutine compute_overlap_fine_to_coarse(nest_domain, overlap, position, name)
  type(nest_level_type),  intent(inout) :: nest_domain
  type(nestSpec),         intent(inout) :: overlap
  integer,                intent(in   ) :: position
  character(len=*),       intent(in   ) :: name

  !--- local variables

  type(domain2D), pointer :: domain_fine  =>NULL()
  type(domain2D), pointer :: domain_coarse=>NULL()
  type(overlap_type), allocatable :: overlapList(:)
  logical              :: is_first
  integer              :: tile_fine, tile_coarse
  integer              :: istart_fine, iend_fine, jstart_fine, jend_fine
  integer              :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
  integer              :: whalo, ehalo, shalo, nhalo
  integer              :: npes_fine, npes_coarse, n
  integer              :: isg_fine, ieg_fine, jsg_fine, jeg_fine
  integer              :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
  integer              :: is_coarse, ie_coarse, js_coarse, je_coarse
  integer              :: is_fine, ie_fine, js_fine, je_fine
  integer              :: isc_fine, iec_fine, jsc_fine, jec_fine
  integer              :: is_you, ie_you, js_you, je_you
  integer              :: x_refine, y_refine
  integer              :: nsend, nrecv, dir
  integer, allocatable :: isl_coarse(:), iel_coarse(:), jsl_coarse(:), jel_coarse(:)
  integer, allocatable :: isl_fine(:), iel_fine(:), jsl_fine(:), jel_fine(:)
  integer              :: is_convert, ie_convert, js_convert, je_convert
  integer              :: is_convert2(2), ie_convert2(2), js_convert2(2), je_convert2(2), rotate2(2)
  integer              :: rotate, is2, ie2, js2, je2, nconvert
  integer              :: xbegin_c, xend_c, ybegin_c, yend_c
  integer              :: ishift, jshift, l, is3, ie3, js3, je3, nn

  domain_fine   => nest_domain%domain_fine
  domain_coarse => nest_domain%domain_coarse
  npes_fine     = size(nest_domain%pelist_fine(:))
  npes_coarse   = size(nest_domain%pelist_coarse(:))

  allocate(isl_coarse(npes_coarse), iel_coarse(npes_coarse) )
  allocate(jsl_coarse(npes_coarse), jel_coarse(npes_coarse) )
  allocate(isl_fine(npes_fine), iel_fine(npes_fine) )
  allocate(jsl_fine(npes_fine), jel_fine(npes_fine) )
  call mpp_get_domain_shift   (domain_coarse, ishift, jshift, position)

  call mpp_get_compute_domain (domain_coarse, xbegin=isc_coarse, xend=iec_coarse, ybegin=jsc_coarse, yend=jec_coarse)
  call mpp_get_compute_domain (domain_fine,   xbegin=isc_fine,   xend=iec_fine,   ybegin=jsc_fine,   yend=jec_fine)
  call mpp_get_compute_domains(domain_coarse, xbegin=isl_coarse, xend=iel_coarse, ybegin=jsl_coarse, yend=jel_coarse)
  call mpp_get_compute_domains(domain_fine,   xbegin=isl_fine,   xend=iel_fine,   ybegin=jsl_fine,   yend=jel_fine)
  call mpp_get_global_domain (domain_fine,   xbegin=isg_fine,   xend=ieg_fine,   ybegin=jsg_fine,   yend=jeg_fine)
  overlap%center%is_you = 0; overlap%center%ie_you = -1
  overlap%center%js_you = 0; overlap%center%je_you = -1

  overlap%nsend = 0
  overlap%nrecv = 0
  call init_index_type(overlap%center)

  if( nest_domain%is_fine_pe ) then
     overlap%xbegin = 0; overlap%xend = -1
     overlap%ybegin = 0; overlap%yend = -1
  else
     overlap%xbegin = isc_coarse - domain_coarse%whalo
     overlap%xend   = iec_coarse + domain_coarse%ehalo + ishift
     overlap%ybegin = jsc_coarse - domain_coarse%shalo
     overlap%yend   = jec_coarse + domain_coarse%nhalo + jshift
     overlap%xsize_c  = overlap%xend - overlap%xbegin + 1
     overlap%ysize_c  = overlap%yend - overlap%ybegin + 1
     overlap%xbegin_f = 0
     overlap%xend_f   = -1
     overlap%ybegin_f = 0
     overlap%yend_f   = -1
     overlap%xbegin_c = 0
     overlap%xend_c   = -1
     overlap%ybegin_c = 0
     overlap%yend_c   = -1
  endif

  if(nest_domain%is_fine_pe) then
     nsend = 0
     allocate(overLaplist(npes_coarse))
     do nn = 1, nest_domain%num_nest
        tile_fine     = nest_domain%tile_fine(nn)
        tile_coarse   = nest_domain%tile_coarse(nn)
        istart_fine   = nest_domain%istart_fine(nn)
        iend_fine     = nest_domain%iend_fine(nn)
        jstart_fine   = nest_domain%jstart_fine(nn)
        jend_fine     = nest_domain%jend_fine(nn)
        istart_coarse = nest_domain%istart_coarse(nn)
        iend_coarse   = nest_domain%iend_coarse(nn)
        jstart_coarse = nest_domain%jstart_coarse(nn)
        jend_coarse   = nest_domain%jend_coarse(nn)
        x_refine      = nest_domain%x_refine
        y_refine      = nest_domain%y_refine

     !--- set up the data range for fine and coarse grid.
     !--- on coarse grid pelist, xbegin_f, ybegin_f, xend_f, yend_f is dummy value
     !--- on fine grid pelist, xbegin_c, xend_c, ybegin_c, yend_c are the coarse grid index that
     !--- the fine grid overlapped with.
     !--- One coarse grid box might overlap with multiple fine grid processor. We use
     !--- the west/south/southwest processor to store the coarse grid data.
        if(tile_fine .NE. domain_fine%tile_id(1)) cycle
        is_coarse = istart_coarse + (isc_fine-istart_fine)/x_refine
        ie_coarse = istart_coarse + (iec_fine-istart_fine)/x_refine
        if(mod(isc_fine-istart_fine, x_refine) .NE. 0 ) is_coarse = is_coarse + 1
        js_coarse = jstart_coarse + (jsc_fine-jstart_fine)/y_refine
        je_coarse = jstart_coarse + (jec_fine-jstart_fine)/y_refine
        if(mod(jsc_fine-jstart_fine, y_refine) .NE. 0 ) js_coarse = js_coarse + 1
        overlap%xbegin_c = is_coarse
        overlap%xend_c   = ie_coarse
        overlap%ybegin_c = js_coarse
        overlap%yend_c   = je_coarse
        overlap%xbegin_f = istart_fine + (overlap%xbegin_c-istart_coarse)*x_refine
        overlap%xend_f = istart_fine + (overlap%xend_c-istart_coarse+1)*x_refine - 1
        overlap%ybegin_f = jstart_fine + (overlap%ybegin_c-jstart_coarse)*y_refine
        overlap%yend_f = jstart_fine + (overlap%yend_c-jstart_coarse+1)*y_refine - 1
        xbegin_c = overlap%xbegin_c
        xend_c   = overlap%xend_c
        ybegin_c = overlap%ybegin_c
        yend_c   = overlap%yend_c
        !      if(iec_fine == ieg_fine) then
        overlap%xend_c   = overlap%xend_c + ishift
        overlap%xend_f   = overlap%xend_f + ishift
        !      endif
        !      if(jec_fine == jeg_fine) then
        overlap%yend_c   = overlap%yend_c + jshift
        overlap%yend_f   = overlap%yend_f + jshift
        !      endif

        overlap%xsize_c  = overlap%xend_c - overlap%xbegin_c + 1
        overlap%ysize_c  = overlap%yend_c - overlap%ybegin_c + 1

        !-----------------------------------------------------------------------------------------
        !
        !    Sending From fine to coarse.
        !    compute the send information from fine grid to coarse grid. This will only need to send
        !    the internal of fine grid to coarse grid.
        !-----------------------------------------------------------------------------------------
        do n = 1, npes_coarse
           nconvert = convert_index_to_nest(domain_coarse, ishift, jshift, tile_coarse, istart_coarse, iend_coarse, &
                jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%list(n-1)%tile_id(1), &
                isl_coarse(n), iel_coarse(n), jsl_coarse(n), jel_coarse(n), &
                is_convert2, ie_convert2, js_convert2, je_convert2, rotate2)
           is2 = xbegin_c; ie2 = xend_c
           js2 = ybegin_c; je2 = yend_c
           is3 = is2; js3 = js2
           do l = 1, nconvert
              if(rotate2(l) == NINETY .OR. rotate2(l) == MINUS_NINETY) then
                 ie3 = ie2 + jshift
                 je3 = je2 + ishift
              else
                 ie3 = ie2 + ishift
                 je3 = je2 + jshift
              endif
              is_coarse = max( is3, is_convert2(l) )
              ie_coarse = min( ie3,   ie_convert2(l) )
              js_coarse = max( js3, js_convert2(l) )
              je_coarse = min( je3,   je_convert2(l) )
              if(ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                 dir = 0
                 nsend = nsend + 1
                 call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                 call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_coarse(n), &
                      is_coarse, ie_coarse, js_coarse, je_coarse, dir, rotate2(l))
              endif
           enddo
        enddo
     enddo
     overlap%nsend = nsend
     if(nsend > 0) then
        allocate(overlap%send(nsend))
        do n = 1, nsend
           call copy_nest_overlap(overlap%send(n), overlaplist(n) )
!           call print_nest_overlap(overlap%send(n), "SEND")
           call deallocate_nest_overlap(overlaplist(n))
        enddo
     endif
     if(allocated(overlaplist))deallocate(overlaplist)
  endif
     !--------------------------------------------------------------------------------
     !   compute the recv information from fine grid to coarse grid. This will only need to send
     !   the internal of fine grid to coarse grid.
     !--------------------------------------------------------------------------------

  if( nest_domain%is_coarse_pe ) then
     nrecv = 0
     allocate(overLaplist(npes_fine))
     do nn = 1, nest_domain%num_nest
        tile_fine     = nest_domain%tile_fine(nn)
        tile_coarse   = nest_domain%tile_coarse(nn)
        istart_fine   = nest_domain%istart_fine(nn)
        iend_fine     = nest_domain%iend_fine(nn)
        jstart_fine   = nest_domain%jstart_fine(nn)
        jend_fine     = nest_domain%jend_fine(nn)
        istart_coarse = nest_domain%istart_coarse(nn)
        iend_coarse   = nest_domain%iend_coarse(nn)
        jstart_coarse = nest_domain%jstart_coarse(nn)
        jend_coarse   = nest_domain%jend_coarse(nn)
        x_refine      = nest_domain%x_refine
        y_refine      = nest_domain%y_refine

        dir = 0
        do n = 1, npes_fine
           if(tile_fine .NE. domain_fine%list(n-1)%tile_id(1)) cycle
           is_you = istart_coarse + (isl_fine(n)-istart_fine)/x_refine
           ie_you = istart_coarse + (iel_fine(n)-istart_fine)/x_refine
           if(mod(isl_fine(n)-istart_fine, x_refine) .NE. 0 ) is_you = is_you + 1
           js_you = jstart_coarse + (jsl_fine(n)-jstart_fine)/y_refine
           je_you = jstart_coarse + (jel_fine(n)-jstart_fine)/y_refine
           if(mod(jsl_fine(n)-jstart_fine, y_refine) .NE. 0 ) js_you = js_you + 1
           nconvert=convert_index_to_coarse(domain_coarse, ishift, jshift, tile_coarse, istart_coarse, iend_coarse, &
                      jstart_coarse, jend_coarse, domain_coarse%ntiles, domain_coarse%tile_id(1), is_you, ie_you, js_you, je_you, &
                      is_convert2, ie_convert2, js_convert2, je_convert2, rotate2)
           do l = 1, nconvert
              is2 = max(is_convert2(l), isc_coarse)
              ie2 = min(ie_convert2(l), iec_coarse+ishift)
              js2 = max(js_convert2(l), jsc_coarse)
              je2 = min(je_convert2(l), jec_coarse+jshift)

              if( ie2 .GE. is2 .AND. je2 .GE. js2 ) then
                 nrecv = nrecv + 1
                 call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                 call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_fine(n), &
                      is2, ie2, js2, je2, dir, rotate2(l))
              endif
           enddo
        enddo
     enddo
     overlap%nrecv = nrecv
        if(nrecv > 0) then
           allocate(overlap%recv(nrecv))
           do n = 1, nrecv
              call copy_nest_overlap(overlap%recv(n), overlaplist(n) )
!              call print_nest_overlap(overlap%recv(n), "RECV")
              call deallocate_nest_overlap( overLaplist(n) )
           enddo
        endif
        if(allocated(overlaplist))deallocate(overlaplist)

  endif

  if(debug_message_passing) call debug_message_size(overlap, name)

  deallocate(isl_coarse, iel_coarse, jsl_coarse, jel_coarse)
  deallocate(isl_fine, iel_fine, jsl_fine, jel_fine)

end subroutine compute_overlap_fine_to_coarse

function find_index(array, data, start_pos)
  integer, intent(in) :: array(:)
  integer, intent(in) :: data
  integer, intent(in) :: start_pos
  integer             :: find_index
  integer             :: i

  find_index = 0
  do i = start_pos, size(array)
     if(array(i) == data) then
        find_index = i
        exit
     endif
  enddo
  if(find_index == 0) then
     print*, "start_pos = ", start_pos, data, array
     call mpp_error(FATAL, "mpp_define_nest_domains.inc: can not find data in array")
  endif

end function find_index

subroutine  debug_message_size(overlap, name)
   type(nestSpec),   intent(in) :: overlap
   character(len=*), intent(in) :: name
   integer, allocatable :: msg1(:), msg2(:), msg3(:), pelist(:)
   integer :: m, n, l, npes, msgsize
   integer :: is, ie, js, je, outunit

      outunit = stdout()
      npes = mpp_npes()
      allocate(msg1(npes), msg2(npes), msg3(npes) )
      allocate(pelist(npes))
      call mpp_get_current_pelist(pelist)
      msg1 = 0
      msg2 = 0
      msg3 = 0
      l = 0
      do m = 1, overlap%nrecv
         msgsize = 0
         do n = 1, overlap%recv(m)%count
            is = overlap%recv(m)%is(n); ie = overlap%recv(m)%ie(n)
            js = overlap%recv(m)%js(n); je = overlap%recv(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         l = find_index(pelist, overlap%recv(m)%pe, l+1)
         msg2(l) = msgsize
      enddo
      l = 0
      do m = 1, overlap%nsend
         msgsize = 0
         do n = 1, overlap%send(m)%count
            is = overlap%send(m)%is(n); ie = overlap%send(m)%ie(n)
            js = overlap%send(m)%js(n); je = overlap%send(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         l = find_index(pelist, overlap%send(m)%pe, l+1)
         msg3(l) = msgsize
      enddo

      call mpp_alltoall(msg3, 1, msg1, 1)

      do m = 1, npes
         if(msg1(m) .NE. msg2(m)) then
            print*, "debug_message_size: My pe = ", mpp_pe(), ",name =", trim(name),", from pe=", &
                 pelist(m), ":send size = ", msg1(m), ", recv size = ", msg2(m)
            call mpp_error(FATAL, "debug_message_size: mismatch on send and recv size")
         endif
      enddo
      write(outunit,*)"NOTE from compute_overlap_fine_to_coarse: "// &
           "message sizes are matched between send and recv for "//trim(name)
      deallocate(msg1, msg2, msg3, pelist)

end subroutine  debug_message_size

!###############################################################################

subroutine init_index_type (indexData )
   type(index_type), intent(inout) :: indexData

     indexData%is_me  = 0
     indexData%ie_me  = -1
     indexData%js_me  = 0
     indexData%je_me  = -1
     indexData%is_you = 0
     indexData%ie_you = -1
     indexData%js_you = 0
     indexData%je_you = -1

end subroutine init_index_type

subroutine allocate_nest_overlap(overlap, count)
  type(overlap_type), intent(inout) :: overlap
  integer,            intent(in   ) :: count

  overlap%count = 0
  overlap%pe    = NULL_PE
  if( ASSOCIATED(overlap%is) ) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc: overlap is already been allocated")

  allocate(overlap%is          (count) )
  allocate(overlap%ie          (count) )
  allocate(overlap%js          (count) )
  allocate(overlap%je          (count) )
  allocate(overlap%dir         (count) )
  allocate(overlap%rotation    (count) )
  allocate(overlap%msgsize     (count) )

end subroutine allocate_nest_overlap

!##############################################################################
subroutine deallocate_nest_overlap(overlap)
  type(overlap_type), intent(inout) :: overlap

  overlap%count = 0
  overlap%pe    = NULL_PE
  deallocate(overlap%is)
  deallocate(overlap%ie)
  deallocate(overlap%js)
  deallocate(overlap%je)
  deallocate(overlap%dir)
  deallocate(overlap%rotation)
  deallocate(overlap%msgsize)

end subroutine deallocate_nest_overlap

!##############################################################################
subroutine insert_nest_overlap(overlap, pe, is, ie, js, je, dir, rotation)
  type(overlap_type), intent(inout) :: overlap
  integer,            intent(in   ) :: pe
  integer,            intent(in   ) :: is, ie, js, je
  integer,            intent(in   ) :: dir, rotation
  integer                           :: count

  if( overlap%count == 0 ) then
     overlap%pe = pe
  else
     if(overlap%pe .NE. pe) call mpp_error(FATAL,  &
          "mpp_define_nest_domains.inc: mismatch on pe")
  endif
  overlap%count = overlap%count+1
  count = overlap%count
  if(count > size(overlap%is(:))) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc: overlap%count > size(overlap%is), contact developer")
  overlap%is          (count) = is
  overlap%ie          (count) = ie
  overlap%js          (count) = js
  overlap%je          (count) = je
  overlap%dir         (count) = dir
  overlap%rotation    (count) = rotation
  overlap%msgsize     (count) = (ie-is+1)*(je-js+1)

end subroutine insert_nest_overlap

subroutine print_nest_overlap(overlap, msg)
  type(overlap_type), intent(in) :: overlap
  character(len=*),   intent(in) :: msg

  integer :: i
  write(1000+mpp_pe(),*) trim(msg), ",pe=",overlap%pe, overlap%count
  do i = 1, overlap%count
     write(1000+mpp_pe(),*) trim(msg), ",index=", overlap%is(i), overlap%ie(i),overlap%js(i),overlap%je(i)
     write(1000+mpp_pe(),*) trim(msg), ",rotation=", overlap%dir(i), overlap%rotation(i), overlap%msgsize(i)
  enddo
  call flush(1000+mpp_pe())

end subroutine print_nest_overlap

!#########################################################
subroutine copy_nest_overlap(overlap_out, overlap_in)
  type(overlap_type), intent(inout) :: overlap_out
  type(overlap_type), intent(in)    :: overlap_in

  if(overlap_in%count == 0) call mpp_error(FATAL, &
    "mpp_define_nest_domains.inc: overlap_in%count is 0")

  if(associated(overlap_out%is)) call mpp_error(FATAL, &
    "mpp_define_nest_domains.inc: overlap_out is already been allocated")

  call allocate_nest_overlap(overlap_out, overlap_in%count)
  overlap_out%count = overlap_in%count
  overlap_out%pe    = overlap_in%pe

  overlap_out%is(:) = overlap_in%is(1:overlap_in%count)
  overlap_out%ie(:)       = overlap_in%ie(1:overlap_in%count)
  overlap_out%js(:)       = overlap_in%js(1:overlap_in%count)
  overlap_out%je(:)       = overlap_in%je(1:overlap_in%count)
  overlap_out%is(:)       = overlap_in%is(1:overlap_in%count)
  overlap_out%dir(:)      = overlap_in%dir(1:overlap_in%count)
  overlap_out%rotation(:) = overlap_in%rotation(1:overlap_in%count)
  overlap_out%msgsize(:)  = overlap_in%msgsize(1:overlap_in%count)


end subroutine copy_nest_overlap


!#######################################################################
  ! this routine found the domain has the same halo size with the input
  ! whalo, ehalo,
function search_C2F_nest_overlap(nest_domain, nest_level, extra_halo, position)
    type(nest_domain_type), intent(inout) :: nest_domain
    integer,                intent(in)    :: extra_halo
    integer,                intent(in)    :: position, nest_level
    type(nestSpec),         pointer       :: search_C2F_nest_overlap
    type(nestSpec),        pointer        :: update_ref
    character(len=128)                    :: name

    if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(search_C2F_nest_overlap): nest_level should be between 1 and nest_domain%num_level")

    select case(position)
    case (CENTER)
       name = trim(nest_domain%name)//" T-cell"
       update_ref => nest_domain%nest(nest_level)%C2F_T
    case (CORNER)
       update_ref => nest_domain%nest(nest_level)%C2F_C
    case (NORTH)
       update_ref => nest_domain%nest(nest_level)%C2F_N
    case (EAST)
       update_ref => nest_domain%nest(nest_level)%C2F_E
    case default
       call mpp_error(FATAL,"mpp_define_nest_domains.inc(search_C2F_nest_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

    search_C2F_nest_overlap => update_ref

    do
       if(extra_halo == search_C2F_nest_overlap%extra_halo) then
            exit ! found domain
       endif
       !--- if not found, switch to next
       if(.NOT. ASSOCIATED(search_C2F_nest_overlap%next)) then
          allocate(search_C2F_nest_overlap%next)
          search_C2F_nest_overlap => search_C2F_nest_overlap%next
          call compute_overlap_coarse_to_fine(nest_domain%nest(nest_level), search_C2F_nest_overlap, extra_halo, position, name)
          exit
       else
          search_C2F_nest_overlap => search_C2F_nest_overlap%next
       end if

    end do

    update_ref => NULL()

  end function search_C2F_nest_overlap

!#######################################################################
  ! this routine found the domain has the same halo size with the input
  ! whalo, ehalo,
  function search_F2C_nest_overlap(nest_domain, nest_level, position)
    type(nest_domain_type), intent(inout) :: nest_domain
    integer,                intent(in)    :: position, nest_level
    type(nestSpec),         pointer       :: search_F2C_nest_overlap

    if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(search_F2C_nest_overlap): nest_level should be between 1 and nest_domain%num_level")

    select case(position)
    case (CENTER)
       search_F2C_nest_overlap => nest_domain%nest(nest_level)%F2C_T
    case (CORNER)
       search_F2C_nest_overlap => nest_domain%nest(nest_level)%F2C_C
    case (NORTH)
       search_F2C_nest_overlap => nest_domain%nest(nest_level)%F2C_N
    case (EAST)
       search_F2C_nest_overlap => nest_domain%nest(nest_level)%F2C_E
    case default
       call mpp_error(FATAL,"mpp_define_nest_domains.inc(search_F2C_nest_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

  end function search_F2C_nest_overlap

  !################################################################
  subroutine mpp_get_C2F_index(nest_domain, is_fine, ie_fine, js_fine, je_fine, &
                is_coarse, ie_coarse, js_coarse, je_coarse, dir, nest_level, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_fine, ie_fine, js_fine, je_fine
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer,                intent(in ) :: dir, nest_level
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

    if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_C2F_index): nest_level should be between 1 and nest_domain%num_level")

     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%nest(nest_level)%C2F_T
     case (EAST)
        update => nest_domain%nest(nest_level)%C2F_E
     case (CORNER)
        update => nest_domain%nest(nest_level)%C2F_C
     case (NORTH)
        update => nest_domain%nest(nest_level)%C2F_N
     case default
        call mpp_error(FATAL, "mpp_define_nest_domains.inc(mpp_get_C2F_index): invalid option argument position")
     end select

     select case(dir)
     case(WEST)
        is_fine = update%west%is_me
        ie_fine = update%west%ie_me
        js_fine = update%west%js_me
        je_fine = update%west%je_me
        is_coarse = update%west%is_you
        ie_coarse = update%west%ie_you
        js_coarse = update%west%js_you
        je_coarse = update%west%je_you
     case(EAST)
        is_fine = update%east%is_me
        ie_fine = update%east%ie_me
        js_fine = update%east%js_me
        je_fine = update%east%je_me
        is_coarse = update%east%is_you
        ie_coarse = update%east%ie_you
        js_coarse = update%east%js_you
        je_coarse = update%east%je_you
     case(SOUTH)
        is_fine = update%south%is_me
        ie_fine = update%south%ie_me
        js_fine = update%south%js_me
        je_fine = update%south%je_me
        is_coarse = update%south%is_you
        ie_coarse = update%south%ie_you
        js_coarse = update%south%js_you
        je_coarse = update%south%je_you
     case(NORTH)
        is_fine = update%north%is_me
        ie_fine = update%north%ie_me
        js_fine = update%north%js_me
        je_fine = update%north%je_me
        is_coarse = update%north%is_you
        ie_coarse = update%north%ie_you
        js_coarse = update%north%js_you
        je_coarse = update%north%je_you
     case default
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: invalid value for argument dir")
     end select


  end subroutine mpp_get_C2F_index

  !################################################################
  subroutine mpp_get_F2C_index_fine(nest_domain, is_coarse, ie_coarse, js_coarse, je_coarse, &
                is_fine, ie_fine, js_fine, je_fine, nest_level, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_fine, ie_fine, js_fine, je_fine
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer,                intent(in)  :: nest_level
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

    if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_F2C_index): nest_level should be between 1 and nest_domain%num_level")

     if(.not. nest_domain%nest(nest_level)%on_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_F2C_index_fine): nest_domain%nest(nest_level)%on_level is false")

     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%nest(nest_level)%F2C_T
     case (EAST)
        update => nest_domain%nest(nest_level)%F2C_E
     case (CORNER)
        update => nest_domain%nest(nest_level)%F2C_C
     case (NORTH)
        update => nest_domain%nest(nest_level)%F2C_N
     case default
        call mpp_error(FATAL, "mpp_define_nest_domains.inc(mpp_get_F2C_index): invalid option argument position")
     end select
     is_fine   = update%xbegin_f
     ie_fine   = update%xend_f
     js_fine   = update%ybegin_f
     je_fine   = update%yend_f
     is_coarse = update%xbegin_c
     ie_coarse = update%xend_c
     js_coarse = update%ybegin_c
     je_coarse = update%yend_c

  end subroutine mpp_get_F2C_index_fine

  !################################################################
  subroutine mpp_get_F2C_index_coarse(nest_domain, is_coarse, ie_coarse, js_coarse, je_coarse, nest_level, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer,                intent(in ) :: nest_level
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_F2C_index_coarse): nest_level should be between 1 and nest_domain%num_level")


     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%nest(nest_level)%F2C_T
     case (EAST)
        update => nest_domain%nest(nest_level)%F2C_E
     case (CORNER)
        update => nest_domain%nest(nest_level)%F2C_C
     case (NORTH)
        update => nest_domain%nest(nest_level)%F2C_N
     case default
        call mpp_error(FATAL, "mpp_define_nest_domains.inc(mpp_get_F2C_index_coarse): invalid option argument position")
     end select
     is_coarse = update%xbegin_c
     ie_coarse = update%xend_c
     js_coarse = update%ybegin_c
     je_coarse = update%yend_c

  end subroutine mpp_get_F2C_index_coarse

  subroutine get_coarse_index(rotate, is, ie, js, je, iadd, jadd, is_c, ie_c, js_c, je_c)
     integer, intent(in)  :: rotate, is, ie, js, je, iadd, jadd
     integer, intent(out) :: is_c, ie_c, js_c, je_c

     if(rotate == 0) then
        is_c = is; ie_c = ie
        js_c = js; je_c = je
     else
        is_c = js; ie_c = je
        js_c = is; je_c = ie
     endif
     is_c = is_c + iadd; ie_c = ie_c + iadd
     js_c = js_c + jadd; je_c = je_c + jadd

  end subroutine get_coarse_index

  !--- this routine will get number of nest.
  subroutine get_nnest(domain, num_nest, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                       x_refine, y_refine, nnest, t_coarse, ncross_coarse, rotate_coarse,                     &
                       is_coarse, ie_coarse, js_coarse, je_coarse, is_fine, ie_fine, js_fine, je_fine)
    type(domain2D), intent(in) :: domain
    integer, intent(in)  :: num_nest, istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
    integer, intent(in)  :: tile_coarse(:)
    integer, intent(in)  :: x_refine, y_refine
    integer, intent(out) :: nnest, is_coarse(:), ie_coarse(:), js_coarse(:), je_coarse(:)
    integer, intent(out) :: is_fine(:), ie_fine(:), js_fine(:), je_fine(:)
    integer, intent(out) :: t_coarse(:), ncross_coarse(:), rotate_coarse(:)
    integer :: is, ie, js, je, tile, isg, ieg, jsg, jeg
    integer :: ncross, rotate, i1, i2
    integer :: is_c, ie_c, js_c, je_c
    integer :: n, iadd, jadd


    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    nnest = 0
    do n = 1, num_nest
       is = istart_coarse(n); ie = iend_coarse(n)
       js = jstart_coarse(n); je = jend_coarse(n)
       tile = tile_coarse(n)
       iadd = 0 ; jadd = 0
       ncross = 0
       rotate = 0
       do while ( ie .GE. is .AND. je .GE. js)
          nnest = nnest+1
          t_coarse(nnest) = tile
          ncross_coarse(nnest) = ncross
          rotate_coarse(nnest) = rotate
          !--- rotate should be 0, 90 or -90.
          if(rotate .NE. 0 .AND. rotate .NE. 90 .AND. rotate .NE. -90) then
             call mpp_error(FATAL, "get_nnest: roate should be 0, 90 or -90")
          endif
          if( ieg .GE. ie .AND. jeg .GE. je) then
             is_coarse(nnest) = is; ie_coarse(nnest) = ie
             js_coarse(nnest) = js; je_coarse(nnest) = je
             call get_coarse_index(rotate, is_coarse(nnest), ie_coarse(nnest), js_coarse(nnest), je_coarse(nnest), &
                                   iadd, jadd, is_c, ie_c, js_c, je_c)
             is = ie + 1; js = je + 1
          else if( ieg .GE. ie ) then  ! jeg < je, will cross the north edge
             is_coarse(nnest) = is; ie_coarse(nnest) = ie
             js_coarse(nnest) = js; je_coarse(nnest) = jeg
             call get_coarse_index(rotate, is_coarse(nnest), ie_coarse(nnest), js_coarse(nnest), je_coarse(nnest), &
                                   iadd, jadd, is_c, ie_c, js_c, je_c)
             if(rotate ==0) then
                jadd = jadd + jeg
             else
                iadd = iadd + ieg
             endif
             js = 1; je = je-jeg
             ncross = ncross+1
             if(mod(tile,2) ==0) then ! tile 2 4 6
                tile = tile + 1
                if(tile>6) tile=tile-6
             else ! rotate 90 degree
                tile = tile + 2
                if(tile>6) tile=tile-6
                i1 = is; i2 = ie
                is = js; ie = je
                js = i1; je = i2
                rotate = rotate + 90
             endif


          else if( jeg .GE. je ) then ! ieg < ie, will cross the east edge
             is_coarse(nnest) = is; ie_coarse(nnest) = ieg
             js_coarse(nnest) = js; je_coarse(nnest) = je
             call get_coarse_index(rotate, is_coarse(nnest), ie_coarse(nnest), js_coarse(nnest), je_coarse(nnest), &
                                   iadd, jadd, is_c, ie_c, js_c, je_c)
             if(rotate ==0) then
                iadd = iadd + ieg
             else
                jadd = jadd + jeg
             endif
             is = 1; ie = ie-ieg
             ncross = ncross+1
             if(mod(tile,2) ==0) then ! rotate -90
                tile = tile + 2
                if(tile>6) tile=tile-6
                i1 = is; i2 = ie
                is = js; ie = je
                js = i1; je = i2
                rotate = rotate - 90
             else
                tile = tile + 1
                if(tile>6) tile=tile-6
             endif
          else
             call mpp_error(FATAL, "get_nnest: do not support cross the corner")
          endif

          !--- is_c:ie_c,js_c:je_c must be inside istart_coarse(n):iend_coarse(n), jstart_coarse(n):jend_coarse(n)
          if(is_c < istart_coarse(n)) call mpp_error(FATAL, "get_nnest: is_c < istart_coarse")
          if(ie_c > iend_coarse(n)) call mpp_error(FATAL, "get_nnest: ie_c > iend_coarse")
          if(js_c < jstart_coarse(n)) call mpp_error(FATAL, "get_nnest: js_c < jstart_coarse")
          if(je_c > jend_coarse(n)) call mpp_error(FATAL, "get_nnest: je_c > jend_coarse")
          is_fine(nnest) = (is_c - istart_coarse(n)) * x_refine + 1
          ie_fine(nnest) = (ie_c - istart_coarse(n)+1) * x_refine
          js_fine(nnest) = (js_c - jstart_coarse(n)) * y_refine + 1
          je_fine(nnest) = (je_c - jstart_coarse(n)+1) * y_refine

          !--- it should not cross the edge more than 3 times.
          if(ncross > 3) call mpp_error(FATAL, "get_nnest: nncross > 3")
       enddo
    enddo


  end subroutine get_nnest


!--- This routine will convert the global coarse grid index to nest grid index.
  function convert_index_to_nest(domain, ishift, jshift, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                       ntiles_coarse, tile_in, is_in, ie_in, js_in, je_in, is_out, ie_out, js_out, je_out, rotate_out)
    type(domain2D), intent(in) :: domain
    integer, intent(in)  :: ishift, jshift
    integer, intent(in)  :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
    integer, intent(in)  :: tile_coarse
    integer, intent(in)  :: ntiles_coarse, tile_in, is_in, ie_in, js_in, je_in
    integer, intent(out) :: is_out(:), ie_out(:), js_out(:), je_out(:), rotate_out(:)
    integer              :: convert_index_to_nest
    integer :: is, ie, js, je, tile, isg, ieg, jsg, jeg
    integer :: ncross, rotate, nout, diff, l, ntiles

    ntiles = ntiles_coarse
    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    is = istart_coarse; ie = iend_coarse
    js = jstart_coarse; je = jend_coarse
    tile = tile_coarse

    if(size(is_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_nest: size(is_out(:)) < 2")
    if(size(ie_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_nest: size(ie_out(:)) < 2")
    if(size(js_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_nest: size(js_out(:)) < 2")
    if(size(je_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_nest: size(je_out(:)) < 2")
    if(size(rotate_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_nest: size(rotate_out(:)) < 2")
    if( ie > ieg .AND. je > jeg) then
       call mpp_error(FATAL, "convert_index_to_nest: do not support cross the corner, contact developer")
    endif
    if( is > ieg .or. js > jeg) call mpp_error(FATAL,"convert_index_to_nest: is > ieg .or. js > jeg")


    nout = 0

    if(tile == tile_in) then
       nout = nout+1
       rotate_out(nout) = ZERO
       is_out(nout) = is_in; ie_out(nout) = ie_in + ishift
       js_out(nout) = js_in; je_out(nout) = je_in + jshift
    endif

    diff =  tile_in - tile
    if(diff < 0) diff = diff + ntiles
    ncross = -1
    if( ie > ieg ) then
       select case(diff)
       case (0)
          rotate = ZERO
          ncross = 4
       case (1)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = ZERO
             ncross = 1
          endif
       case (2)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = MINUS_NINETY
             ncross = 1
          endif
       case (3)
          rotate = MINUS_NINETY
          ncross = 2
       case (4)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = MINUS_NINETY
             ncross = 3
          endif
       case (5)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = ZERO
             ncross = 3
          endif
       case default
          call mpp_error(FATAL,"convert_index_to_nest: invalid value of diff")
       end select

       if(ncross > 0) then
          nout =nout+1
          rotate_out(nout) = rotate
          if(rotate_out(nout) == ZERO) then
             js_out(nout) = js_in
             je_out(nout) = je_in + jshift
             is_out(nout) = is_in+ncross*ieg
             ie_out(nout) = ie_in+ncross*ieg + ishift
          else if(rotate_out(nout) == MINUS_NINETY) then
             js_out(nout) = ieg-ie_in + 1
             je_out(nout) = ieg-is_in + 1 + ishift
             is_out(nout) = js_in+ncross*jeg
             ie_out(nout) = je_in+ncross*jeg + jshift
          endif
       endif
    else if( je > jeg ) then
       select case(diff)
       case (0)
          rotate = ZERO
          ncross = 4
       case (1)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = ZERO
             ncross = 1
          endif
       case (2)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = NINETY
             ncross = 1
          endif
       case (3)
          rotate = NINETY
          ncross = 2
       case (4)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = NINETY
             ncross = 3
          endif
       case (5)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = ZERO
             ncross = 3
          endif
       end select

       if(ncross > 0) then
          nout =nout+1
          rotate_out(nout) = rotate

          if(rotate_out(nout) == ZERO) then
             js_out(nout) = js_in
             je_out(nout) = je_in + jshift
             is_out(nout) = is_in+ncross*ieg
             ie_out(nout) = ie_in+ncross*ieg + ishift
          else if(rotate_out(nout) == NINETY) then
             is_out(nout) = ieg-je_in + 1
             ie_out(nout) = ieg-js_in+1 + jshift
             js_out(nout) = is_in+ncross*jeg
             je_out(nout) = ie_in+ncross*jeg + ishift
          endif
       endif
    endif

    convert_index_to_nest = nout

  end function convert_index_to_nest

  function convert_index_to_coarse(domain, ishift, jshift, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                             ntiles_coarse, tile_in, is_in, ie_in, js_in, je_in, is_out, ie_out, js_out, je_out, rotate_out)
    type(domain2D), intent(in) :: domain
    integer, intent(in)  :: ishift, jshift
    integer, intent(in)  :: istart_coarse, iend_coarse, jstart_coarse, jend_coarse
    integer, intent(in)  :: tile_coarse
    integer, intent(in)  :: ntiles_coarse, tile_in, is_in, ie_in, js_in, je_in
    integer, intent(out) :: is_out(:), ie_out(:), js_out(:), je_out(:), rotate_out(:)
    integer :: convert_index_to_coarse
    integer :: is, ie, js, je, isg, ieg, jsg, jeg
    integer :: ncross, rotate, ntiles, nout, diff, tile

    ntiles = ntiles_coarse
    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    is = istart_coarse; ie = iend_coarse
    js = jstart_coarse; je = jend_coarse
    tile = tile_coarse

    if(size(is_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_coarse: size(is_out(:)) < 2")
    if(size(ie_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_coarse: size(ie_out(:)) < 2")
    if(size(js_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_coarse: size(js_out(:)) < 2")
    if(size(je_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_coarse: size(je_out(:)) < 2")
    if(size(rotate_out(:)) < 2) call mpp_error(FATAL,"convert_index_to_coarse: size(rotate_out(:)) < 2")
    if( ie > ieg .AND. je > jeg) then
       call mpp_error(FATAL, "convert_index_to_coarse: do not support cross the corner, contact developer")
    endif
    if( is > ieg .or. js > jeg) call mpp_error(FATAL,"convert_index_to_coarse: is > ieg .or. js > jeg")

    nout = 0

    if(tile_coarse == tile_in) then
       nout = nout+1
       rotate_out(nout) = ZERO
       is_out(nout) = is_in; ie_out(nout) = ie_in + ishift
       js_out(nout) = js_in; je_out(nout) = je_in + jshift
    endif

    diff =  tile_in - tile
    if(diff < 0) diff = diff + ntiles
    ncross = -1
    if( ie > ieg ) then
       select case(diff)
       case (0)
          rotate = ZERO
          ncross = 4
       case (1)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = ZERO
             ncross = 1
          endif
       case (2)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = MINUS_NINETY
             ncross = 1
          endif
       case (3)
          rotate = MINUS_NINETY
          ncross = 2
       case (4)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = MINUS_NINETY
             ncross = 3
          endif
       case (5)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = ZERO
             ncross = 3
          endif
       case default
          call mpp_error(FATAL,"convert_index_to_coarse: invalid value of diff")
       end select

       if(ncross > 0) then
          nout =nout+1
          rotate_out(nout) = rotate
          if(rotate_out(nout) == ZERO) then
             js_out(nout) = js_in
             je_out(nout) = je_in + jshift
             is_out(nout) = is_in-ncross*ieg
             ie_out(nout) = ie_in-ncross*ieg + ishift
          else if(rotate_out(nout) == MINUS_NINETY) then
             is_out(nout) = ieg-je_in + 1
             ie_out(nout) = ieg-js_in + 1 + ishift
             js_out(nout) = is_in-ncross*jeg
             je_out(nout) = ie_in-ncross*jeg + jshift
          endif
       endif
    else if( je > jeg ) then
       select case(diff)
       case (0)
          rotate = ZERO
          ncross = 4
       case (1)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = ZERO
             ncross = 1
          endif
       case (2)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = NINETY
             ncross = 1
          endif
       case (3)
          rotate = NINETY
          ncross = 2
       case (4)
          if(mod(tile,2) ==0) then ! tile 2 4 6
             rotate = NINETY
             ncross = 3
          endif
       case (5)
          if(mod(tile,2) ==1) then ! tile 1 3 5
             rotate = ZERO
             ncross = 3
          endif
       end select

       if(ncross > 0) then
          nout =nout+1
          rotate_out(nout) = rotate

          if(rotate_out(nout) == ZERO) then
             js_out(nout) = js_in
             je_out(nout) = je_in + jshift
             is_out(nout) = is_in-ncross*ieg
             ie_out(nout) = ie_in-ncross*ieg + ishift
          else if(rotate_out(nout) == NINETY) then
             is_out(nout) = js_in - ncross*jeg
             ie_out(nout) = je_in - ncross*jeg + ishift
             js_out(nout) = jeg - ie_in + 1
             je_out(nout) = jeg - is_in + 1 + jshift
          endif
       endif
    endif

    convert_index_to_coarse = nout


  end function convert_index_to_coarse


  subroutine convert_index_back(domain, ishift, jshift, rotate, is_in, ie_in, js_in, je_in, is_out, ie_out, js_out, je_out)
    type(domain2D), intent(in) :: domain
    integer, intent(in)  :: ishift, jshift
    integer, intent(in)  :: is_in, ie_in, js_in, je_in, rotate
    integer, intent(out) :: is_out, ie_out, js_out, je_out
    integer :: isg, ieg, jsg, jeg
    integer :: ncross

    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    ncross = 0
    if( je_in > jeg+jshift .and. ie_in > ieg+ishift ) then
       call mpp_error(FATAL,"convert_index_back:  je_in > jeg .and. ie_in > ieg")
    else if (je_in > jeg+jshift) then
       ncross = je_in/jeg
       select case(rotate)
       case(0)
          is_out = is_in
          ie_out = ie_in
          js_out = js_in - ncross*jeg
          je_out = je_in - ncross*jeg
       case(90)
          is_out = js_in - ncross*jeg
          ie_out = je_in - ncross*jeg
          js_out = jeg - ie_in + 1
          je_out = jeg - is_in + 1
       case default
          call mpp_error(FATAL, "convert_index_back: rotate should be 0 or 90 when je_in>jeg")
       end select
    else if (ie_in > ieg+ishift) then
       ncross = ie_in/ieg
       select case(rotate)
       case(0)
          is_out = is_in - ncross*ieg
          ie_out = ie_in - ncross*ieg
          js_out = js_in
          je_out = je_in
       case(-90)
          js_out = is_in - ncross*ieg
          je_out = ie_in - ncross*ieg
          is_out = ieg - je_in + 1
          ie_out = ieg - js_in + 1
       case default
          call mpp_error(FATAL, "convert_index_back: rotate should be 0 or -90 when ie_in>ieg")
       end select
    else
       is_out = is_in
       ie_out = ie_in
       js_out = js_in
       je_out = je_in
    endif

  end subroutine convert_index_back



  function get_nest_vector_recv(nest_domain, update_x, update_y, ind_x, ind_y, start_pos, pelist)
    type(nest_level_type), intent(in) :: nest_domain
    type(nestSpec),    intent(in) :: update_x, update_y
    integer,          intent(out) :: ind_x(:), ind_y(:)
    integer,          intent(out) :: start_pos(:)
    integer,          intent(out) :: pelist(:)
    integer                       :: get_nest_vector_recv
    integer                       :: nlist, nrecv_x, nrecv_y, ntot, n
    integer                       :: ix, iy, rank_x, rank_y, cur_pos
    integer                       :: nrecv

    nlist = size(nest_domain%pelist)
    nrecv_x = update_x%nrecv
    nrecv_y = update_y%nrecv

    ntot = nrecv_x + nrecv_y

    n  = 1
    ix = 1
    iy = 1
    ind_x = -1
    ind_y = -1
    nrecv = 0
    cur_pos = 0
    do while (n<=ntot)
       if(ix <= nrecv_x ) then
          rank_x = update_x%recv(ix)%pe-mpp_pe()
          if(rank_x .LE. 0) rank_x = rank_x + nlist
       else
          rank_x = nlist+1
       endif
       if(iy <= nrecv_y ) then
          rank_y = update_y%recv(iy)%pe-mpp_pe()
          if(rank_y .LE. 0) rank_y = rank_y + nlist
       else
          rank_y = nlist+1
       endif
       nrecv = nrecv + 1
       start_pos(nrecv) = cur_pos
       if( rank_x == rank_y ) then
          n = n+2
          ind_x (nrecv) = ix
          ind_y (nrecv) = iy
          cur_pos = cur_pos + update_x%recv(ix)%totsize + update_y%recv(iy)%totsize
          pelist(nrecv) = update_x%recv(ix)%pe
          ix = ix + 1
          iy = iy + 1
       else if ( rank_x < rank_y ) then
          n = n+1
          ind_x (nrecv) = ix
          ind_y (nrecv) = -1
          cur_pos = cur_pos + update_x%recv(ix)%totsize
          pelist(nrecv) = update_x%recv(ix)%pe
          ix = ix + 1
       else if ( rank_y < rank_x ) then
          n = n+1
          ind_x (nrecv) = -1
          ind_y (nrecv) = iy
          cur_pos = cur_pos + update_y%recv(iy)%totsize
          pelist(nrecv) = update_y%recv(iy)%pe
          iy = iy+1
       endif
    end do

    get_nest_vector_recv = nrecv


  end function get_nest_vector_recv


  function get_nest_vector_send(nest_domain, update_x, update_y, ind_x, ind_y, start_pos, pelist)
    type(nest_level_type), intent(in) :: nest_domain
    type(nestSpec),    intent(in) :: update_x, update_y
    integer,          intent(out) :: ind_x(:), ind_y(:)
    integer,          intent(out) :: start_pos(:)
    integer,          intent(out) :: pelist(:)
    integer                       :: get_nest_vector_send
    integer                       :: nlist, nsend_x, nsend_y, ntot, n
    integer                       :: ix, iy, rank_x, rank_y, cur_pos
    integer                       :: nsend

    nlist = size(nest_domain%pelist_fine(:)) + size(nest_domain%pelist_coarse(:))
    nsend_x = update_x%nsend
    nsend_y = update_y%nsend

    ntot = nsend_x + nsend_y

    n  = 1
    ix = 1
    iy = 1
    ind_x = -1
    ind_y = -1
    nsend = 0
    cur_pos = 0
    do while (n<=ntot)
       if(ix <= nsend_x ) then
          rank_x = update_x%send(ix)%pe-mpp_pe()
          if(rank_x .LE. 0) rank_x = rank_x + nlist
       else
          rank_x = nlist+1
       endif
       if(iy <= nsend_y ) then
          rank_y = update_y%send(iy)%pe-mpp_pe()
          if(rank_y .LE. 0) rank_y = rank_y + nlist
       else
          rank_y = nlist+1
       endif
       nsend = nsend + 1
       start_pos(nsend) = cur_pos
       if( rank_x == rank_y ) then
          n = n+2
          ind_x (nsend) = ix
          ind_y (nsend) = iy
          cur_pos = cur_pos + update_x%send(ix)%totsize + update_y%send(iy)%totsize
          pelist(nsend) = update_x%send(ix)%pe
          ix = ix + 1
          iy = iy + 1
       else if ( rank_x < rank_y ) then
          n = n+1
          ind_x (nsend) = ix
          ind_y (nsend) = -1
          cur_pos = cur_pos + update_x%send(ix)%totsize
          pelist(nsend) = update_x%send(ix)%pe
          ix = ix + 1
       else if ( rank_y < rank_x ) then
          n = n+1
          ind_x (nsend) = -1
          ind_y (nsend) = iy
          cur_pos = cur_pos + update_y%send(iy)%totsize
          pelist(nsend) = update_y%send(iy)%pe
          iy = iy+1
       endif
    end do

    get_nest_vector_send = nsend


  end function get_nest_vector_send

  subroutine check_data_size_1d(module, str1, size1, str2, size2)
    character(len=*), intent(in) :: module, str1, str2
    integer,          intent(in) :: size1, size2


   if(size2 > 0 .AND. size1 .NE. size2 ) then
      print '(a, 3I5)', trim(module), mpp_pe(), size1, size2
      call mpp_error(FATAL, trim(module)//": mismatch between size of "//trim(str1)//" and "//trim(str2))
   endif

  end subroutine check_data_size_1d


  subroutine check_data_size_2d(module, str1, isize1, jsize1, str2, isize2, jsize2)
    character(len=*), intent(in) :: module, str1, str2
    integer,          intent(in) :: isize1, jsize1, isize2, jsize2


   if(isize2 > 0 .AND. jsize2 > 0 .AND. (isize1 .NE. isize2 .OR. jsize1 .NE. jsize2) ) then
      print '(a, 5I5)', trim(module), mpp_pe(), isize1, jsize1, isize2, jsize2
      call mpp_error(FATAL, trim(module)//": mismatch between size of "//trim(str1)//" and "//trim(str2))
   endif

  end subroutine check_data_size_2d

  function mpp_get_nest_coarse_domain(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     type(domain2d), pointer :: mpp_get_nest_coarse_domain

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_coarse_domain): nest_level should be between 1 and nest_domain%num_level")

     if(.not. nest_domain%nest(nest_level)%on_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_coarse_domain): nest_domain%nest(nest_level)%on_level is false")
     mpp_get_nest_coarse_domain => nest_domain%nest(nest_level)%domain_coarse

  end function mpp_get_nest_coarse_domain

  function mpp_get_nest_fine_domain(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     type(domain2d), pointer :: mpp_get_nest_fine_domain

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_fine_domain): nest_level should be between 1 and nest_domain%num_level")

     if(.not. nest_domain%nest(nest_level)%on_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_fine_domain): nest_domain%nest(nest_level)%on_level is false")
     mpp_get_nest_fine_domain => nest_domain%nest(nest_level)%domain_fine

  end function mpp_get_nest_fine_domain

  function mpp_get_nest_npes(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     integer                            :: mpp_get_nest_npes

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_npes): nest_level should be between 1 and nest_domain%num_level")

     mpp_get_nest_npes = size(nest_domain%nest(nest_level)%pelist(:))

  end function mpp_get_nest_npes

  subroutine mpp_get_nest_pelist(nest_domain, nest_level, pelist)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     integer,               intent(out) :: pelist(:)
     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_pelist): nest_level should be between 1 and nest_domain%num_level")

     if(size(pelist) .NE. size(nest_domain%nest(nest_level)%pelist)) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_pelist): size(pelist) .NE. size(nest_domain%nest(nest_level)%pelist)")

     pelist = nest_domain%nest(nest_level)%pelist

  end subroutine mpp_get_nest_pelist

  function mpp_get_nest_fine_npes(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     integer                            :: mpp_get_nest_fine_npes

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_fine_npes): nest_level should be between 1 and nest_domain%num_level")

     mpp_get_nest_fine_npes = size(nest_domain%nest(nest_level)%pelist_fine(:))

  end function mpp_get_nest_fine_npes

  subroutine mpp_get_nest_fine_pelist(nest_domain, nest_level, pelist)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     integer,               intent(out) :: pelist(:)
     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_fine_pelist): nest_level should be between 1 and nest_domain%num_level")

     if(size(pelist) .NE. size(nest_domain%nest(nest_level)%pelist_fine)) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_get_nest_fine_pelist): size(pelist) .NE. size(nest_domain%nest(nest_level)%pelist)")

     pelist = nest_domain%nest(nest_level)%pelist_fine

  end subroutine mpp_get_nest_fine_pelist



  function mpp_is_nest_fine(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     logical                            :: mpp_is_nest_fine

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_is_nest_fine): nest_level should be between 1 and nest_domain%num_level")

     if(.not. nest_domain%nest(nest_level)%on_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_is_nest_fine): nest_domain%nest(nest_level)%on_level is false")

     mpp_is_nest_fine = nest_domain%nest(nest_level)%is_fine_pe

  end function mpp_is_nest_fine

  function mpp_is_nest_coarse(nest_domain, nest_level)
     type(nest_domain_type), intent(in) :: nest_domain
     integer,                intent(in) :: nest_level
     logical                            :: mpp_is_nest_coarse

     if(nest_level < 1 .OR. nest_level > nest_domain%num_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_is_nest_coarse): nest_level should be between 1 and nest_domain%num_level")

     if(.not. nest_domain%nest(nest_level)%on_level) call mpp_error(FATAL, &
       "mpp_define_nest_domains.inc(mpp_is_nest_coarse): nest_domain%nest(nest_level)%on_level is false")

     mpp_is_nest_coarse = nest_domain%nest(nest_level)%is_coarse_pe

  end function mpp_is_nest_coarse
