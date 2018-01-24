! -*-f90-*-


!#############################################################################
! Currently the contact will be limited to overlap contact.
subroutine mpp_define_nest_domains(nest_domain, domain_fine, domain_coarse, num_nest, tile_fine, tile_coarse, &
                                  istart_coarse, icount_coarse, jstart_coarse, jcount_coarse,         &
                                  x_refine, y_refine, pelist, extra_halo, name) 
  type(nest_domain_type),     intent(inout) :: nest_domain
  type(domain2D), target,     intent(in   ) :: domain_fine, domain_coarse
  integer,                    intent(in   ) :: num_nest
  integer,                    intent(in   ) :: tile_fine(:), tile_coarse(:)
  integer,                    intent(in   ) :: istart_coarse(:), icount_coarse(:), jstart_coarse(:), jcount_coarse(:)
  integer,                    intent(in   ) :: x_refine, y_refine
  integer,          optional, intent(in   ) :: pelist(:)
  integer,          optional, intent(in   ) :: extra_halo
  character(len=*), optional, intent(in   ) :: name

  logical                                   :: concurrent
  integer                                   :: n
  integer                                   :: nx_coarse, ny_coarse
  integer                                   :: nx_fine, ny_fine
  integer                                   :: npes, npes_fine, npes_coarse
  integer                                   :: extra_halo_local
  integer, allocatable                      :: pes(:)
  integer, allocatable                      :: pes_coarse(:)
  integer, allocatable                      :: pes_fine(:)
  integer, dimension(num_nest)              :: istart_fine, iend_fine, jstart_fine, jend_fine
  integer, dimension(num_nest)              :: iend_coarse, jend_coarse
  integer, dimension(num_nest)              :: my_nest_id
  integer                                   :: my_num_nest

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

  do n = 1, num_nest
     if(istart_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: istart_coarse < 1")
     if(icount_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: iend_coarse < 1")
     if(jstart_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: jstart_coarse < 1")
     if(jcount_coarse(n) < 1) call mpp_error(FATAL, "mpp_define_nest_domains.inc: jend_coarse < 1")
     iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
     jend_coarse(n) = jstart_coarse(n) + jcount_coarse(n) - 1
     istart_fine(n) = 1 ; iend_fine(n) = icount_coarse(n)*x_refine
     jstart_fine(n) = 1 ; jend_fine(n) = jcount_coarse(n)*y_refine
  end do

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

  !--- check the pelist, Either domain_coarse%pelist = pelist or
  !--- domain_coarse%pelist + domain_fine%pelist = pelist
  if( PRESENT(pelist) )then
     allocate( pes(size(pelist(:))) )
     pes = pelist
  else
    allocate( pes(mpp_npes()) )
    call mpp_get_current_pelist(pes)
  end if

  npes = size(pes)
  npes_coarse = size(domain_coarse%list(:))
  npes_fine = size(domain_fine%list(:))
  !--- pes_fine and pes_coarse should be subset of pelist
  allocate( pes_coarse(npes_coarse) )
  allocate( pes_fine  (npes_fine  ) )
  do n = 1, npes_coarse
     pes_coarse(n) = domain_coarse%list(n-1)%pe
     if( .NOT. ANY(pes(:) == pes_coarse(n)) ) then
        call mpp_error(FATAL, "mpp_domains_define.inc: pelist_coarse is not subset of pelist")
     endif
  enddo
  do n = 1, npes_fine
     pes_fine(n) = domain_fine%list(n-1)%pe
     if( .NOT. ANY(pes(:) == pes_fine(n)) ) then
        call mpp_error(FATAL, "mpp_domains_define.inc: pelist_fine is not subset of pelist")
     endif
  enddo

  allocate(nest_domain%pelist_fine(npes_fine))
  allocate(nest_domain%pelist_coarse(npes_coarse))
  nest_domain%pelist_fine = pes_fine
  nest_domain%pelist_coarse = pes_coarse
  nest_domain%is_fine_pe = ANY(pes_fine(:) == mpp_pe())
  nest_domain%is_coarse_pe = ANY(pes_coarse(:) == mpp_pe())

  !--- figure out the current processor is fine pe or coarse pe, also the corresponding nested region.
  !--- on coarse grid pe, it might overlap multiple fine regon.
  !--- on fine grid pe, it always only overlap at most 1 coarse region.
  my_num_nest= 0 
  my_nest_id(:) = 0
  if( nest_domain%is_fine_pe ) then
     !--- figure out the nest number on current pe
     do n = 1, num_nest
        if(domain_fine%tile_id(1) == tile_fine(n)) then
           my_num_nest = my_num_nest + 1
           my_nest_id(my_num_nest) = n
           exit
        end if
     end do
     if(my_num_nest .NE. 1) then
        print*, "num_nest=", my_num_nest, domain_fine%tile_id(1), tile_fine(1)
        call mpp_error(FATAL, "mpp_define_nest_domains.inc: my_num_nest .ne. 1 on fine pelist")
     endif
  else if( nest_domain%is_coarse_pe ) then
     do n = 1, num_nest
        if(domain_coarse%tile_id(1) == tile_coarse(n)) then
           my_num_nest = my_num_nest + 1
           my_nest_id(my_num_nest) = n
        end if
     end do
  endif

  nest_domain%my_num_nest = my_num_nest
  if(my_num_nest>0) then
     allocate(nest_domain%my_nest_id(my_num_nest))
     nest_domain%my_nest_id(:) = my_nest_id(1:my_num_nest)
  endif
  !-- current only support my_num_nest == 1
  if(my_num_nest>1) call mpp_error(FATAL, "mpp_define_nest_domains.inc:  my_num_nest > 1, contact developer")

  !--- We are assuming the fine grid is fully overlapped with coarse grid.
  if( nest_domain%is_fine_pe ) then
     if( iend_fine(my_nest_id(1)) - istart_fine(my_nest_id(1)) + 1 .NE. domain_fine%x(1)%global%size .OR.  &
        jend_fine(my_nest_id(1)) - jstart_fine(my_nest_id(1)) + 1 .NE. domain_fine%y(1)%global%size ) then
print*, "x size are", domain_fine%x(1)%global%size, istart_fine(my_nest_id(1)), iend_fine(my_nest_id(1)) 
print*, "y size are", domain_fine%y(1)%global%size, jstart_fine(my_nest_id(1)), jend_fine(my_nest_id(1))
        call mpp_error(FATAL, ": The fine global domain is not covered by coarse domain")
     endif
  endif
  ! First computing the send and recv information from fine to coarse.
  if( npes == npes_coarse ) then
     concurrent = .false.
  else if( npes_fine + npes_coarse == npes ) then
     concurrent = .true.
  else
     call mpp_error(FATAL, "mpp_domains_define.inc: size(pelist_coarse) .NE. size(pelist) and "// &
              "size(pelist_coarse)+size(pelist_fine) .NE. size(pelist)")
  endif

  !--- coarse grid and fine grid should be both symmetry or non-symmetry.
  if(domain_coarse%symmetry .neqv. domain_fine%symmetry) then
     print*,"symmetry is", domain_coarse%symmetry, domain_fine%symmetry, mpp_pe()
     call mpp_error(FATAL, "mpp_domains_define.inc: domain_coarse%symmetry .neqv. .NOT. domain_fine%symmetry")
  endif

  nest_domain%x_refine = x_refine
  nest_domain%y_refine = y_refine
  nest_domain%domain_fine   => domain_fine
  nest_domain%domain_coarse => domain_coarse

  allocate( nest_domain%C2F_T, nest_domain%C2F_C, nest_domain%C2F_E, nest_domain%C2F_N )
  nest_domain%C2F_T%next => NULL()
  nest_domain%C2F_C%next => NULL()
  nest_domain%C2F_N%next => NULL()
  nest_domain%C2F_E%next => NULL()
  allocate( nest_domain%F2C_T, nest_domain%F2C_C, nest_domain%F2C_E, nest_domain%F2C_N )

  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_T, CENTER, trim(nest_domain%name)//" T-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_E, EAST,   trim(nest_domain%name)//" E-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_C, CORNER, trim(nest_domain%name)//" C-cell")
  call compute_overlap_fine_to_coarse(nest_domain, nest_domain%F2C_N, NORTH,  trim(nest_domain%name)//" N-cell")

  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_T, extra_halo_local, CENTER, trim(nest_domain%name)//" T-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_E, extra_halo_local, EAST,   trim(nest_domain%name)//" E-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_C, extra_halo_local, CORNER, trim(nest_domain%name)//" C-cell")
  call compute_overlap_coarse_to_fine(nest_domain, nest_domain%C2F_N, extra_halo_local, NORTH,  trim(nest_domain%name)//" N-cell")  

  deallocate(pes, pes_fine, pes_coarse)


end subroutine mpp_define_nest_domains

!###############################################################################
subroutine compute_overlap_coarse_to_fine(nest_domain, overlap, extra_halo, position, name) 
   type(nest_domain_type), intent(inout) :: nest_domain
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
   integer              :: isc_fine, iec_fine, jsc_fine, jec_fine
   integer              :: isd_fine, ied_fine, jsd_fine, jed_fine
   integer              :: isc_east, iec_east, jsc_east, jec_east
   integer              :: isc_west, iec_west, jsc_west, jec_west
   integer              :: isc_south, iec_south, jsc_south, jec_south
   integer              :: isc_north, iec_north, jsc_north, jec_north
   integer              :: x_refine, y_refine, ishift, jshift
   integer              :: nsend, nrecv, dir, from_pe, l
   integer              :: is, ie, js, je, msgsize
   integer, allocatable :: msg1(:), msg2(:)
   integer, allocatable :: isl_coarse(:), iel_coarse(:), jsl_coarse(:), jel_coarse(:)
   integer, allocatable :: isl_fine(:), iel_fine(:), jsl_fine(:), jel_fine(:)
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

   call mpp_get_global_domain  (domain_fine,   xbegin=isg_fine,   xend=ieg_fine,   &
                                ybegin=jsg_fine,   yend=jeg_fine, position=position)
   call mpp_get_compute_domain (domain_coarse, xbegin=isc_coarse, xend=iec_coarse, &
                                ybegin=jsc_coarse, yend=jec_coarse, position=position)
   call mpp_get_compute_domain (domain_fine,   xbegin=isc_fine,   xend=iec_fine,   &
                                ybegin=jsc_fine,   yend=jec_fine, position=position)
   call mpp_get_compute_domains(domain_coarse, xbegin=isl_coarse, xend=iel_coarse, &
                                ybegin=jsl_coarse, yend=jel_coarse, position=position)
   call mpp_get_compute_domains(domain_fine,   xbegin=isl_fine,   xend=iel_fine,   &
                                ybegin=jsl_fine,   yend=jel_fine, position=position)

   overlap%extra_halo = extra_halo


   if(nest_domain%my_num_nest>0) then
      n = nest_domain%my_nest_id(1)
      tile_fine     = nest_domain%tile_fine(n)
      tile_coarse   = nest_domain%tile_coarse(n)
      istart_fine   = nest_domain%istart_fine(n)
      iend_fine     = nest_domain%iend_fine(n)
      jstart_fine   = nest_domain%jstart_fine(n)
      jend_fine     = nest_domain%jend_fine(n)
      istart_coarse = nest_domain%istart_coarse(n)
      iend_coarse   = nest_domain%iend_coarse(n) + ishift
      jstart_coarse = nest_domain%jstart_coarse(n)
      jend_coarse   = nest_domain%jend_coarse(n) + jshift
   endif
   x_refine      = nest_domain%x_refine
   y_refine      = nest_domain%y_refine
   whalo         = domain_fine%whalo + extra_halo
   ehalo         = domain_fine%ehalo + extra_halo
   shalo         = domain_fine%shalo + extra_halo
   nhalo         = domain_fine%nhalo + extra_halo


   if( nest_domain%is_coarse_pe ) then
      overlap%xbegin = isc_coarse - domain_coarse%whalo
      overlap%xend   = iec_coarse + domain_coarse%ehalo
      overlap%ybegin = jsc_coarse - domain_coarse%shalo
      overlap%yend   = jec_coarse + domain_coarse%nhalo 
   else
      overlap%xbegin = isc_fine - domain_fine%whalo
      overlap%xend   = iec_fine + domain_fine%ehalo
      overlap%ybegin = jsc_fine - domain_fine%shalo
      overlap%yend   = jec_fine + domain_fine%nhalo
      !--- figure out the coarse grid index. The index does not contain halo and it is possibly cross the edge.


   endif

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

   !--- first compute the halo region and corresponding index in coarse grid.
   if( nest_domain%is_fine_pe ) then 
      if( ieg_fine == iec_fine .AND. domain_fine%tile_id(1) == tile_fine ) then   ! east halo
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

      if( jsg_fine == jsc_fine .AND. domain_fine%tile_id(1) == tile_fine) then  ! south
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

      if( isg_fine == isc_fine .AND. domain_fine%tile_id(1) == tile_fine) then ! west
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

      if( jeg_fine == jec_fine .AND. domain_fine%tile_id(1) == tile_fine) then ! north
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

      allocate(overLaplist(npes_coarse))

      !-------------------------------------------------------------------------
      !
      !                 Receiving  
      !
      !-------------------------------------------------------------------------
      !--- loop through coarse pelist
      nrecv = 0
      do n = 1, npes_coarse   
         if( domain_coarse%list(n-1)%tile_id(1) .NE. tile_coarse ) cycle
         is_first = .true.
         !--- east halo receiving
         is_coarse = overlap%east%is_you
         ie_coarse = overlap%east%ie_you
         js_coarse = overlap%east%js_you
         je_coarse = overlap%east%je_you
         if( je_coarse .GE. js_coarse .AND. ie_coarse .GE. is_coarse ) then
            dir = 1
            is_coarse = max( is_coarse, isl_coarse(n) )
            ie_coarse = min( ie_coarse, iel_coarse(n) )
            js_coarse = max( js_coarse, jsl_coarse(n) )
            je_coarse = min( je_coarse, jel_coarse(n) )
            if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
               if(is_first) then
                  nrecv = nrecv + 1
                  call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                  is_first = .false.
               endif
               call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_coarse(n), &
                    is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
            endif
         endif

         !--- south halo receiving
         is_coarse = overlap%south%is_you
         ie_coarse = overlap%south%ie_you
         js_coarse = overlap%south%js_you
         je_coarse = overlap%south%je_you
         if( je_coarse .GE. js_coarse .AND. ie_coarse .GE. is_coarse ) then
            dir = 3
            is_coarse = max( is_coarse, isl_coarse(n) )
            ie_coarse = min( ie_coarse, iel_coarse(n) )
            js_coarse = max( js_coarse, jsl_coarse(n) )
            je_coarse = min( je_coarse, jel_coarse(n) )

            if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
               if(is_first) then
                  nrecv = nrecv + 1
                  call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                  is_first = .false.
               endif
               call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_coarse(n), &
                    is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
            endif
         endif

         !--- west halo receiving
         is_coarse = overlap%west%is_you
         ie_coarse = overlap%west%ie_you
         js_coarse = overlap%west%js_you
         je_coarse = overlap%west%je_you
         if( je_coarse .GE. js_coarse .AND. ie_coarse .GE. is_coarse ) then
            dir = 5
            is_coarse = max( is_coarse, isl_coarse(n) )
            ie_coarse = min( ie_coarse, iel_coarse(n) )
            js_coarse = max( js_coarse, jsl_coarse(n) )
            je_coarse = min( je_coarse, jel_coarse(n) )

            if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
               if(is_first) then
                  nrecv = nrecv + 1
                  call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                  is_first = .false.
               endif
               call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_coarse(n), &
                    is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
            endif
         endif

         !--- north halo receiving
         is_coarse = overlap%north%is_you
         ie_coarse = overlap%north%ie_you
         js_coarse = overlap%north%js_you
         je_coarse = overlap%north%je_you
         if( je_coarse .GE. js_coarse .AND. ie_coarse .GE. is_coarse ) then
            dir = 7
            is_coarse = max( is_coarse, isl_coarse(n) )
            ie_coarse = min( ie_coarse, iel_coarse(n) )
            js_coarse = max( js_coarse, jsl_coarse(n) )
            je_coarse = min( je_coarse, jel_coarse(n) )

            if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
               if(is_first) then
                  nrecv = nrecv + 1
                  call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                  is_first = .false.
               endif
               call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_coarse(n), &
                    is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
            endif
         endif
      enddo

      !--- copy the overlapping into nest_domain data.
      overlap%nrecv = nrecv
      if( nrecv > 0 ) then
         allocate(overlap%recv(nrecv))
         do n = 1, nrecv
            call copy_nest_overlap( overlap%recv(n), overLaplist(n) )
            call deallocate_nest_overlap( overLaplist(n) )
         enddo
      endif
      if(allocated(overlaplist))deallocate(overlapList)
   endif
  !-----------------------------------------------------------------------
  !
  !                          Sending 
  !
  !-----------------------------------------------------------------------

   if( nest_domain%is_coarse_pe .and. nest_domain%my_num_nest>0 ) then 
      nsend = 0
      if(domain_coarse%tile_id(1) == tile_coarse) then
         isc_east = iend_coarse
         iec_east = iend_coarse + ehalo
         jsc_east = jstart_coarse - shalo
         jec_east = jend_coarse + nhalo
         isc_east = max(isc_coarse, isc_east)
         iec_east = min(iec_coarse, iec_east)
         jsc_east = max(jsc_coarse, jsc_east)
         jec_east = min(jec_coarse, jec_east)

         isc_south = istart_coarse - whalo
         iec_south = iend_coarse + ehalo
         jsc_south = jstart_coarse - shalo
         jec_south = jstart_coarse
         isc_south = max(isc_coarse, isc_south)
         iec_south = min(iec_coarse, iec_south)
         jsc_south = max(jsc_coarse, jsc_south)
         jec_south = min(jec_coarse, jec_south)

         isc_west = istart_coarse - whalo
         iec_west = istart_coarse
         jsc_west = jstart_coarse - shalo
         jec_west = jend_coarse + nhalo
         isc_west = max(isc_coarse, isc_west)
         iec_west = min(iec_coarse, iec_west)
         jsc_west = max(jsc_coarse, jsc_west)
         jec_west = min(jec_coarse, jec_west)

         isc_north = istart_coarse - whalo
         iec_north = iend_coarse + ehalo
         jsc_north = jend_coarse
         jec_north = jend_coarse + nhalo
         isc_north = max(isc_coarse, isc_north)
         iec_north = min(iec_coarse, iec_north)
         jsc_north = max(jsc_coarse, jsc_north)
         jec_north = min(jec_coarse, jec_north)
      else
         isc_west = 0; iec_west = -1; jsc_west = 0; jec_west = -1
         isc_east = 0; iec_east = -1; jsc_east = 0; jec_west = -1
         isc_south = 0; iec_south = -1; jsc_south = 0; jec_south = -1
         isc_north = 0; iec_north = -1; jsc_north = 0; jec_north = -1
      endif

      allocate(overLaplist(npes_fine))

      do n = 1, npes_fine
         if( domain_fine%list(n-1)%tile_id(1) .NE. tile_fine ) cycle
         is_first = .true.

         !--- to_pe's east
         if( ieg_fine == iel_fine(n) ) then
            dir = 1
            if( iec_east .GE. isc_east .AND. jec_east .GE. jsc_east ) then
               is_coarse = iend_coarse
               ie_coarse = iend_coarse + ehalo
               js_coarse = jstart_coarse + ( jsl_fine(n) - jsg_fine )/y_refine  
               je_coarse = jstart_coarse + ( jel_fine(n) - jsg_fine )/y_refine
               js_coarse = js_coarse - shalo
               je_coarse = je_coarse + nhalo
               is_coarse = max(isc_east, is_coarse)
               ie_coarse = min(iec_east, ie_coarse)
               js_coarse = max(jsc_east, js_coarse)
               je_coarse = min(jec_east, je_coarse)
               if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                  if(is_first) then
                     nsend = nsend + 1
                     call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                     is_first = .false.
                  endif
                  call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
               endif
            endif
         endif

         !--- to_pe's south
         if( jsg_fine == jsl_fine(n) ) then
            dir = 3
            if( iec_south .GE. isc_south .AND. jec_south .GE. jsc_south ) then
               is_coarse = istart_coarse + ( isl_fine(n) - isg_fine )/x_refine  
               ie_coarse = istart_coarse + ( iel_fine(n) - isg_fine )/x_refine
               is_coarse = is_coarse - shalo
               ie_coarse = ie_coarse + nhalo
               js_coarse = jstart_coarse - shalo
               je_coarse = jstart_coarse
               is_coarse = max(isc_south, is_coarse)
               ie_coarse = min(iec_south, ie_coarse)
               js_coarse = max(jsc_south, js_coarse)
               je_coarse = min(jec_south, je_coarse)
               if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                  if(is_first) then
                     nsend = nsend + 1
                     call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                     is_first = .false.
                  endif
                  call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
               endif
            endif
         endif

         !--- to_pe's west
         if( isg_fine == isl_fine(n) ) then
            dir = 5
            if( iec_west .GE. isc_west .AND. jec_west .GE. jsc_west ) then
               is_coarse = istart_coarse - whalo
               ie_coarse = istart_coarse
               js_coarse = jstart_coarse + ( jsl_fine(n) - jsg_fine )/y_refine  
               je_coarse = jstart_coarse + ( jel_fine(n) - jsg_fine )/y_refine
               js_coarse = js_coarse - shalo
               je_coarse = je_coarse + nhalo
               is_coarse = max(isc_west, is_coarse)
               ie_coarse = min(iec_west, ie_coarse)
               js_coarse = max(jsc_west, js_coarse)
               je_coarse = min(jec_west, je_coarse)
               if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                  if(is_first) then
                     nsend = nsend + 1
                     call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                     is_first = .false.
                  endif
                  call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
               endif
            endif
         endif

         !--- to_pe's north
         if( jeg_fine == jel_fine(n) ) then
            dir = 7
            if( iec_north .GE. isc_north .AND. jec_north .GE. jsc_north ) then
               is_coarse = istart_coarse + ( isl_fine(n) - isg_fine )/x_refine  
               ie_coarse = istart_coarse + ( iel_fine(n) - isg_fine )/x_refine
               is_coarse = is_coarse - shalo
               ie_coarse = ie_coarse + nhalo
               js_coarse = jend_coarse
               je_coarse = jend_coarse + nhalo
               is_coarse = max(isc_north, is_coarse)
               ie_coarse = min(iec_north, ie_coarse)
               js_coarse = max(jsc_north, js_coarse)
               je_coarse = min(jec_north, je_coarse)
               if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                  if(is_first) then
                     nsend = nsend + 1
                     call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
                     is_first = .false.
                  endif
                  call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_fine(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse , dir, ZERO)
               endif
            endif
         endif
      enddo

      !--- copy the overlapping into nest_domain data.
      overlap%nsend = nsend
      if( nsend > 0 ) then
         allocate(overlap%send(nsend))
         do n = 1, nsend
            call copy_nest_overlap( overlap%send(n), overLaplist(n) )
            call deallocate_nest_overlap( overLaplist(n) )
         enddo
      endif
      if(allocated(overlaplist))deallocate(overLaplist)
   endif


   deallocate(isl_coarse, iel_coarse, jsl_coarse, jel_coarse)
   deallocate(isl_fine, iel_fine, jsl_fine, jel_fine)

   if(debug_message_passing) then
       allocate(msg1(0:npes-1), msg2(0:npes-1) )
       msg1 = 0
       msg2 = 0
       do m = 1, overlap%nrecv
          msgsize = 0
          do n = 1, overlap%recv(m)%count
             is = overlap%recv(m)%is(n); ie = overlap%recv(m)%ie(n)
             js = overlap%recv(m)%js(n); je = overlap%recv(m)%je(n)
             msgsize = msgsize + (ie-is+1)*(je-js+1)
          end do
          from_pe = overlap%recv(m)%pe
          l = from_pe-mpp_root_pe()
          call mpp_recv( msg1(l), glen=1, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_1)
          msg2(l) = msgsize
       enddo

       do m = 1, overlap%nsend
          msgsize = 0
          do n = 1, overlap%send(m)%count
             is = overlap%send(m)%is(n); ie = overlap%send(m)%ie(n)
             js = overlap%send(m)%js(n); je = overlap%send(m)%je(n)
             msgsize = msgsize + (ie-is+1)*(je-js+1)
          end do
          call mpp_send( msgsize, plen=1, to_pe=overlap%send(m)%pe, tag=COMM_TAG_1)
       enddo
       call mpp_sync_self(check=EVENT_RECV)

       do m = 0, npes-1
          if(msg1(m) .NE. msg2(m)) then
             print*, "compute_overlap_coarse_to_fine: My pe = ", mpp_pe(), ",name =", trim(name),", from pe=", &
                     m+mpp_root_pe(), ":send size = ", msg1(m), ", recv size = ", msg2(m)
             call mpp_error(FATAL, "mpp_compute_overlap_coarse_to_fine: mismatch on send and recv size")
          endif
       enddo
       call mpp_sync_self()
       write(outunit,*)"NOTE from compute_overlap_coarse_to_fine: "// &
               "message sizes are matched between send and recv for "//trim(name)
       deallocate(msg1, msg2)
    endif


end subroutine compute_overlap_coarse_to_fine

!###############################################################################
!-- This routine will compute the send and recv information between overlapped nesting 
!-- region. The data is assumed on T-cell center. 
subroutine compute_overlap_fine_to_coarse(nest_domain, overlap, position, name)
   type(nest_domain_type), intent(inout) :: nest_domain
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
   integer              :: npes, npes_fine, npes_coarse, n, m
   integer              :: isg_fine, ieg_fine, jsg_fine, jeg_fine
   integer              :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
   integer              :: is_coarse, ie_coarse, js_coarse, je_coarse
   integer              :: is_fine, ie_fine, js_fine, je_fine
   integer              :: isc_fine, iec_fine, jsc_fine, jec_fine
   integer              :: is_you, ie_you, js_you, je_you
   integer              :: x_refine, y_refine, ishift, jshift
   integer              :: nsend, nrecv, dir, from_pe, l
   integer              :: is, ie, js, je, msgsize
   integer, allocatable :: msg1(:), msg2(:)
   integer, allocatable :: isl_coarse(:), iel_coarse(:), jsl_coarse(:), jel_coarse(:)
   integer, allocatable :: isl_fine(:), iel_fine(:), jsl_fine(:), jel_fine(:)
   integer              :: outunit

   outunit = stdout()
   domain_fine   => nest_domain%domain_fine
   domain_coarse => nest_domain%domain_coarse
   npes          = mpp_npes()
   npes_fine     = size(nest_domain%pelist_fine(:))
   npes_coarse   = size(nest_domain%pelist_coarse(:))

   allocate(isl_coarse(npes_coarse), iel_coarse(npes_coarse) )
   allocate(jsl_coarse(npes_coarse), jel_coarse(npes_coarse) )
   allocate(isl_fine(npes_fine), iel_fine(npes_fine) )
   allocate(jsl_fine(npes_fine), jel_fine(npes_fine) )

   call mpp_get_compute_domain (domain_coarse, xbegin=isc_coarse, xend=iec_coarse, ybegin=jsc_coarse, yend=jec_coarse)
   call mpp_get_compute_domain (domain_fine,   xbegin=isc_fine,   xend=iec_fine,   ybegin=jsc_fine,   yend=jec_fine)
   call mpp_get_compute_domains(domain_coarse, xbegin=isl_coarse, xend=iel_coarse, ybegin=jsl_coarse, yend=jel_coarse)
   call mpp_get_compute_domains(domain_fine,   xbegin=isl_fine,   xend=iel_fine,   ybegin=jsl_fine,   yend=jel_fine)
   call mpp_get_global_domain (domain_fine,   xbegin=isg_fine,   xend=ieg_fine,   ybegin=jsg_fine,   yend=jeg_fine)
   call mpp_get_domain_shift   (domain_coarse, ishift, jshift, position)
   overlap%center%is_you = 0; overlap%center%ie_you = -1
   overlap%center%js_you = 0; overlap%center%je_you = -1

   overlap%nsend = 0
   overlap%nrecv = 0
   call init_index_type(overlap%center)

   if(nest_domain%my_num_nest >0) then
      n = nest_domain%my_nest_id(1)
      tile_fine     = nest_domain%tile_fine(n)
      tile_coarse   = nest_domain%tile_coarse(n)
      istart_fine   = nest_domain%istart_fine(n)
      iend_fine     = nest_domain%iend_fine(n)
      jstart_fine   = nest_domain%jstart_fine(n)
      jend_fine     = nest_domain%jend_fine(n)
      istart_coarse = nest_domain%istart_coarse(n)
      iend_coarse   = nest_domain%iend_coarse(n)
      jstart_coarse = nest_domain%jstart_coarse(n)
      jend_coarse   = nest_domain%jend_coarse(n)
      x_refine      = nest_domain%x_refine
      y_refine      = nest_domain%y_refine
   endif

   !--- set up the data range for fine and coarse grid.
   !--- on coarse grid pelist, xbegin_f, ybegin_f, xend_f, yend_f is dummy value
   !--- on fine grid pelist, xbegin_c, xend_c, ybegin_c, yend_c are the coarse grid index that 
   !--- the fine grid overlapped with. 
   !--- One coarse grid box might overlap with multiple fine grid processor. We use 
   !--- the west/south/southwest processor to store the coarse grid data.

   if( nest_domain%is_fine_pe  ) then
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
      if(iec_fine == ieg_fine) then
         overlap%xend_c   = overlap%xend_c + ishift
         overlap%xend_f   = overlap%xend_f + ishift
      endif
      if(jec_fine == jeg_fine) then
         overlap%yend_c   = overlap%yend_c + jshift
         overlap%yend_f   = overlap%yend_f + jshift
      endif
      overlap%xsize_c  = overlap%xend_c - overlap%xbegin_c + 1
      overlap%ysize_c  = overlap%yend_c - overlap%ybegin_c + 1
   else
      overlap%xbegin_f = 0
      overlap%xend_f   = -1
      overlap%ybegin_f = 0
      overlap%yend_f   = -1
      overlap%xbegin_c = 0
      overlap%xend_c   = -1
      overlap%ybegin_c = 0
      overlap%yend_c   = -1
  
      if( nest_domain%my_num_nest > 0 ) then
        if(domain_coarse%tile_id(1) == tile_coarse) then
           is_coarse = max( istart_coarse, isc_coarse )
           ie_coarse = min( iend_coarse,   iec_coarse ) + ishift
           js_coarse = max( jstart_coarse, jsc_coarse )
           je_coarse = min( jend_coarse,   jec_coarse ) + jshift
           if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
              overlap%xbegin_c = is_coarse
              overlap%xend_c = ie_coarse
              overlap%ybegin_c = js_coarse
              overlap%yend_c = je_coarse
           endif
        endif
      endif
      overlap%xbegin = isc_coarse - domain_coarse%whalo
      overlap%xend   = iec_coarse + domain_coarse%ehalo + ishift
      overlap%ybegin = jsc_coarse - domain_coarse%shalo
      overlap%yend   = jec_coarse + domain_coarse%nhalo + jshift
      overlap%xsize_c  = overlap%xend - overlap%xbegin + 1
      overlap%ysize_c  = overlap%yend - overlap%ybegin + 1
   endif


   !-----------------------------------------------------------------------------------------
   !
   !    Sending From fine to coarse.  
   !    compute the send information from fine grid to coarse grid. This will only need to send
   !    the internal of fine grid to coarse grid.
   !-----------------------------------------------------------------------------------------  
   nsend = 0
   if( nest_domain%is_fine_pe ) then 
      allocate(overLaplist(npes_coarse))
      do n = 1, npes_coarse
         if(domain_coarse%list(n-1)%tile_id(1) == tile_coarse) then
            is_coarse = max( overlap%xbegin_c, isl_coarse(n) )
            ie_coarse = min( overlap%xend_c,   iel_coarse(n)+ishift )
            js_coarse = max( overlap%ybegin_c, jsl_coarse(n) )
            je_coarse = min( overlap%yend_c,   jel_coarse(n)+jshift )
            if(ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
               dir = 0
               nsend = nsend + 1
               call allocate_nest_overlap(overLaplist(nsend), MAXOVERLAP)
               call insert_nest_overlap(overLaplist(nsend), nest_domain%pelist_coarse(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse, dir, ZERO)
               write(mpp_pe()+1000,*) "send:", is_coarse, ie_coarse, js_coarse, je_coarse, nest_domain%pelist_coarse(n)
               call flush(mpp_pe()+1000)
            endif
         endif
      enddo
      overlap%nsend = nsend
      if(nsend > 0) then
         allocate(overlap%send(nsend))
         do n = 1, nsend
            call copy_nest_overlap(overlap%send(n), overlaplist(n) )      
            call deallocate_nest_overlap(overlaplist(n))  
         enddo
      endif
      if(allocated(overlaplist))deallocate(overlaplist)
   endif

   !--------------------------------------------------------------------------------
   !   compute the recv information from fine grid to coarse grid. This will only need to send
   !   the internal of fine grid to coarse grid.
   !--------------------------------------------------------------------------------
   
   if( nest_domain%is_coarse_pe .and. nest_domain%my_num_nest >0 ) then 
      nrecv = 0
      if(domain_coarse%tile_id(1) == tile_coarse) then
         if(overlap%xend_c .GE. overlap%xbegin_c .AND. overlap%yend_c .GE. overlap%ybegin_c ) then
            dir = 0
            allocate(overLaplist(npes_fine))
            do n = 1, npes_fine
               is_you = istart_coarse + (isl_fine(n)-istart_fine)/x_refine
               ie_you = istart_coarse + (iel_fine(n)-istart_fine)/x_refine
               if(mod(isl_fine(n)-istart_fine, x_refine) .NE. 0 ) is_you = is_you + 1
               js_you = jstart_coarse + (jsl_fine(n)-jstart_fine)/y_refine
               je_you = jstart_coarse + (jel_fine(n)-jstart_fine)/y_refine
               if(mod(jsl_fine(n)-jstart_fine, y_refine) .NE. 0 ) js_you = js_you + 1
               if( iel_fine(n) == ieg_fine ) ie_you = ie_you + ishift
               if( jel_fine(n) == jeg_fine ) je_you = je_you + jshift
               is_coarse = max(overlap%xbegin_c, is_you)
               ie_coarse = min(overlap%xend_c, ie_you)
               js_coarse = max(overlap%ybegin_c, js_you)
               je_coarse = min(overlap%yend_c, je_you)
               if( ie_coarse .GE. is_coarse .AND. je_coarse .GE. js_coarse ) then
                  nrecv = nrecv + 1
                  call allocate_nest_overlap(overLaplist(nrecv), MAXOVERLAP)
                  call insert_nest_overlap(overLaplist(nrecv), nest_domain%pelist_fine(n), &
                       is_coarse, ie_coarse, js_coarse, je_coarse, dir, ZERO)
              write(mpp_pe()+1000,*) "recv:", is_coarse, ie_coarse, js_coarse, je_coarse, nest_domain%pelist_fine(n)
               call flush(mpp_pe()+1000)
                  
               endif
            enddo
         endif
      endif
      overlap%nrecv = nrecv
      if(nrecv > 0) then
         allocate(overlap%recv(nrecv))
         do n = 1, nrecv
            call copy_nest_overlap(overlap%recv(n), overlaplist(n) )      
            call deallocate_nest_overlap( overLaplist(n) )
         enddo
      endif
      if(allocated(overlaplist))deallocate(overlaplist)

   endif

   deallocate(isl_coarse, iel_coarse, jsl_coarse, jel_coarse)
   deallocate(isl_fine, iel_fine, jsl_fine, jel_fine)

   if(debug_message_passing) then
      allocate(msg1(0:npes-1), msg2(0:npes-1) )
      msg1 = 0
      msg2 = 0
      do m = 1, overlap%nrecv
         msgsize = 0
         do n = 1, overlap%recv(m)%count
            is = overlap%recv(m)%is(n); ie = overlap%recv(m)%ie(n)
            js = overlap%recv(m)%js(n); je = overlap%recv(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         from_pe = overlap%recv(m)%pe
         l = from_pe-mpp_root_pe()
         call mpp_recv( msg1(l), glen=1, from_pe=from_pe, block=.FALSE., tag=COMM_TAG_2)
         msg2(l) = msgsize
      enddo

      do m = 1, overlap%nsend
         msgsize = 0
         do n = 1, overlap%send(m)%count
            is = overlap%send(m)%is(n); ie = overlap%send(m)%ie(n)
            js = overlap%send(m)%js(n); je = overlap%send(m)%je(n)
            msgsize = msgsize + (ie-is+1)*(je-js+1)
         end do
         call mpp_send( msgsize, plen=1, to_pe=overlap%send(m)%pe, tag=COMM_TAG_2)
      enddo
      call mpp_sync_self(check=EVENT_RECV)

      do m = 0, npes-1
         if(msg1(m) .NE. msg2(m)) then
            print*, "compute_overlap_fine_to_coarse: My pe = ", mpp_pe(), ",name =", trim(name),", from pe=", &
                 m+mpp_root_pe(), ":send size = ", msg1(m), ", recv size = ", msg2(m)
            call mpp_error(FATAL, "mpp_compute_overlap_coarse_to_fine: mismatch on send and recv size")
         endif
      enddo
      call mpp_sync_self()
      write(outunit,*)"NOTE from compute_overlap_fine_to_coarse: "// &
           "message sizes are matched between send and recv for "//trim(name)
      deallocate(msg1, msg2)
   endif

end subroutine compute_overlap_fine_to_coarse



!!$subroutine set_overlap_fine_to_coarse(nest_domain, position)
!!$   type(nest_domain_type), intent(inout) :: nest_domain
!!$   integer,                intent(in   ) :: position
!!$
!!$
!!$  call mpp_get_domain_shift(domain, ishift, jshift, position)
!!$  update_in => nest_domain%F2C_T
!!$  select case(position)
!!$  case (CORNER)
!!$     update_out => nest_domain%F2C_C
!!$  case (EAST)
!!$     update_out => nest_domain%F2C_E
!!$  case (NORTH)
!!$     update_out => nest_domain%F2C_N
!!$  case default
!!$     call mpp_error(FATAL, "mpp_domains_define.inc(set_overlap_fine_to_coarse): the position should be CORNER, EAST or NORTH")
!!$  end select
!!$
!!$  nsend = update_in%nsend
!!$  nrecv = update_in%nrecv
!!$  update_out%pe = update_in%pe
!!$  update_out%nsend = nsend
!!$  update_out%nrecv = nrecv
!!$
!!$  if( nsend > 0 ) then
!!$     allocate(update_out%send(nsend))
!!$     do n = 1, nsend
!!$        count = update_in%send(n)%count
!!$        call allocate_overlap_type(update_out%send(n), update_in%count, overlap_in%type)
!!$        do m = 1, count
!!$           update_out%send(n)%is      (count) = update_in%send(n)%is      (count)
!!$           update_out%send(n)%ie      (count) = update_in%send(n)%ie      (count) + ishift  
!!$           update_out%send(n)%js      (count) = update_in%send(n)%js      (count)
!!$           update_out%send(n)%je      (count) = update_in%send(n)%je      (count) + jshift    
!!$           update_out%send(n)%tileMe  (count) = update_in%send(n)%tileMe  (count)
!!$           update_out%send(n)%dir     (count) = update_in%send(n)%dir     (count)    
!!$           update_out%send(n)%rotation(count) = update_in%send(n)%rotation(count)
!!$        enddo
!!$     enddo
!!$  endif
!!$
!!$
!!$  if( nrecv > 0 ) then
!!$     allocate(update_out%recv(nrecv))
!!$     do n = 1, nrecv
!!$        count = update_in%recv(n)%count
!!$        call allocate_overlap_type(update_out%recv(n), update_in%count, overlap_in%type)
!!$        do m = 1, count
!!$           update_out%recv(n)%is      (count) = update_in%recv(n)%is      (count)
!!$           update_out%recv(n)%ie      (count) = update_in%recv(n)%ie      (count) + ishift  
!!$           update_out%recv(n)%js      (count) = update_in%recv(n)%js      (count)
!!$           update_out%recv(n)%je      (count) = update_in%recv(n)%je      (count) + jshift    
!!$           update_out%recv(n)%tileMe  (count) = update_in%recv(n)%tileMe  (count)
!!$           update_out%recv(n)%dir     (count) = update_in%recv(n)%dir     (count)    
!!$           update_out%recv(n)%rotation(count) = update_in%recv(n)%rotation(count)
!!$        enddo
!!$     enddo
!!$  endif
!!$
!!$end subroutine set_overlap_fine_to_coarse


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
function search_C2F_nest_overlap(nest_domain, extra_halo, position)
    type(nest_domain_type), intent(inout) :: nest_domain
    integer,                intent(in)    :: extra_halo
    integer,                intent(in)    :: position
    type(nestSpec),         pointer       :: search_C2F_nest_overlap
    type(nestSpec),        pointer        :: update_ref
    character(len=128)                    :: name

    select case(position)
    case (CENTER)
       name = trim(nest_domain%name)//" T-cell"
       update_ref => nest_domain%C2F_T
    case (CORNER)
       update_ref => nest_domain%C2F_C
    case (NORTH)
       update_ref => nest_domain%C2F_N
    case (EAST)
       update_ref => nest_domain%C2F_E
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
          call compute_overlap_coarse_to_fine(nest_domain, search_C2F_nest_overlap, extra_halo, position, name)
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
  function search_F2C_nest_overlap(nest_domain, position)
    type(nest_domain_type), intent(inout) :: nest_domain
    integer,                intent(in)    :: position
    type(nestSpec),         pointer       :: search_F2C_nest_overlap

    select case(position)
    case (CENTER)
       search_F2C_nest_overlap => nest_domain%F2C_T
    case (CORNER)
       search_F2C_nest_overlap => nest_domain%F2C_C
    case (NORTH)
       search_F2C_nest_overlap => nest_domain%F2C_N
    case (EAST)
       search_F2C_nest_overlap => nest_domain%F2C_E
    case default
       call mpp_error(FATAL,"mpp_define_nest_domains.inc(search_F2C_nest_overlap): position should be CENTER|CORNER|EAST|NORTH")
    end select

  end function search_F2C_nest_overlap

  !################################################################
  subroutine mpp_get_C2F_index(nest_domain, is_fine, ie_fine, js_fine, je_fine, &
                is_coarse, ie_coarse, js_coarse, je_coarse, dir, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_fine, ie_fine, js_fine, je_fine
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer,                intent(in ) :: dir
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%C2F_T
     case (EAST)
        update => nest_domain%C2F_E
     case (CORNER)
        update => nest_domain%C2F_C
     case (NORTH)
        update => nest_domain%C2F_N
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
                is_fine, ie_fine, js_fine, je_fine, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_fine, ie_fine, js_fine, je_fine
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%F2C_T
     case (EAST)
        update => nest_domain%F2C_E
     case (CORNER)
        update => nest_domain%F2C_C
     case (NORTH)
        update => nest_domain%F2C_N
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
  subroutine mpp_get_F2C_index_coarse(nest_domain, is_coarse, ie_coarse, js_coarse, je_coarse, position)

     type(nest_domain_type), intent(in ) :: nest_domain
     integer,                intent(out) :: is_coarse, ie_coarse, js_coarse, je_coarse
     integer, optional,      intent(in ) :: position

     integer                             :: update_position
     type(nestSpec), pointer             :: update => NULL()

     update_position = CENTER
     if(present(position)) update_position = position

     select case(update_position)
     case (CENTER)
        update => nest_domain%F2C_T
     case (EAST)
        update => nest_domain%F2C_E
     case (CORNER)
        update => nest_domain%F2C_C
     case (NORTH)
        update => nest_domain%F2C_N
     case default
        call mpp_error(FATAL, "mpp_define_nest_domains.inc(mpp_get_F2C_index): invalid option argument position")
     end select
     is_coarse = update%xbegin_c
     ie_coarse = update%xend_c
     js_coarse = update%ybegin_c
     je_coarse = update%yend_c

  end subroutine mpp_get_F2C_index_coarse



