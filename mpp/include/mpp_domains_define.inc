
! -*-f90-*-
! $Id: mpp_domains_define.inc,v 14.0.4.1 2007/04/25 13:26:14 z1l Exp $


  ! <SUBROUTINE NAME="mpp_define_layout2D" INTERFACE="mpp_define_layout">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="ndivs" TYPE="integer"></IN>
  !  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
  ! </SUBROUTINE>
  subroutine mpp_define_layout2D( global_indices, ndivs, layout )
    integer, intent(in) :: global_indices(:) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: ndivs !number of divisions to divide global domain
    integer, intent(out) :: layout(:)

    integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

    if(size(global_indices(:)) .NE. 4) call mpp_error(FATAL,"mpp_define_layout2D: size of global_indices should be 4")
    if(size(layout(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_layout2D: size of layout should be 2")

    isg = global_indices(1)
    ieg = global_indices(2)
    jsg = global_indices(3)
    jeg = global_indices(4)

    isz = ieg - isg + 1
    jsz = jeg - jsg + 1
    !first try to divide ndivs in the domain aspect ratio: if imperfect aspect, reduce idiv till it divides ndivs
    idiv = nint( sqrt(float(ndivs*isz)/jsz) )
    idiv = max(idiv,1) !for isz=1 line above can give 0
    do while( mod(ndivs,idiv).NE.0 )
       idiv = idiv - 1
    end do                 !will terminate at idiv=1 if not before
    jdiv = ndivs/idiv

    layout = (/ idiv, jdiv /)
    return
  end subroutine mpp_define_layout2D

  !############################################################################
  ! <SUBROUTINE NAME="mpp_define_mosaic_pelist">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
  ! </SUBROUTINE>
  ! NOTE: The following routine may need to revised to improve the capability.
  !  It is very hard to make it balance for all the situation.
  !  Hopefully some smart idea will come up someday.
  subroutine mpp_define_mosaic_pelist( sizes, pe_start, pe_end, pelist, costpertile)
    integer, dimension(:), intent(in)           :: sizes
    integer, dimension(:), intent(inout)        :: pe_start, pe_end
    integer, dimension(:), intent(in), optional :: pelist, costpertile
    integer, dimension(:), allocatable          :: costs, pes
    integer                                     :: ntiles, npes, totcosts, avgcost
    integer                                     :: ntiles_left, npes_left, pos, n, tile
    integer                                     :: cost_on_tile, cost_on_pe, npes_used, errunit

    ntiles = size(sizes(:))
    if(size(pe_start(:)) .NE. ntiles .OR. size(pe_end(:)) .NE. ntiles ) then
       call mpp_error(FATAL, "mpp_define_mosaic_pelist: size mismatch between pe_start/pe_end and sizes")
    end if

    allocate(costs(ntiles))   
    if(present(costpertile)) then
       if(size(costpertile(:)) .NE. ntiles ) then
          call mpp_error(FATAL, "mpp_define_mosaic_pelist: size mismatch between costpertile and sizes")
       end if
       costs = sizes*costpertile
    else
       costs = sizes
    end if

    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.mpp_pe()) )then
          errunit = stderr()
          write( errunit,* )'pe=', mpp_pe(), ' pelist=', pelist
          call mpp_error( FATAL, 'mpp_define_mosaic_pelist: pe must be in pelist.' )
       end if
       npes = size(pelist(:))
       allocate( pes(0:npes-1) )
       pes(:) = pelist(:)
    else
       npes = mpp_npes()
       allocate( pes(0:npes-1) )
       call mpp_get_current_pelist(pes)
    end if

    ntiles_left = ntiles
    npes_left = npes
    pos = pes(0)

    do while( ntiles_left > 0 )
       if( npes_left == 1 ) then ! all left tiles will on the last processor, imbalance possibly.
          do n = 1, ntiles
             if(costs(n) > 0) then
                pe_start(n) = pos
                pe_end(n) = pos
                costs(n)  = 0
             end if
          end do
          ntiles_left = 0
          npes_left = 0
       else
          totcosts = sum(costs)   
          avgcost  = CEILING(real(totcosts)/npes_left )
          tile = minval(maxloc(costs))
          cost_on_tile = costs(tile)  
          pe_start(tile) = pos
          ntiles_left = ntiles_left - 1
          costs(tile) = 0
          totcosts = totcosts - cost_on_tile
          if(cost_on_tile .GE. avgcost ) then
             npes_used = min(ceiling(real(cost_on_tile)/avgcost), npes_left)
             if( ntiles_left > 0 .AND. npes_used == npes_left ) npes_used = npes_used - 1
             pe_end(tile) = pos + npes_used - 1
             npes_left = npes_left - npes_used
             pos = pos + npes_used
          else
             !--- find other tiles to share the pe
             pe_end(tile) = pos
             cost_on_pe = cost_on_tile
             do while(ntiles_left>npes_left)  ! make sure all the pes are used.
                tile = minval(minloc(costs, costs> 0 ))
                cost_on_tile = costs(tile)
                cost_on_pe = cost_on_pe + cost_on_tile
                if(cost_on_pe > avgcost ) exit
                pe_start(tile) = pos
                pe_end(tile) = pos
                ntiles_left = ntiles_left - 1
                costs(tile) = 0
                totcosts = totcosts - cost_on_tile
             end do
             npes_left = npes_left - 1
             pos = pos + 1
          end if
       end if
    end do

    if(npes_left .NE. 0 ) call mpp_error(FATAL, "mpp_define_mosaic_pelist: the left npes should be zero")

  end subroutine mpp_define_mosaic_pelist

  !#####################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !              MPP_DEFINE_DOMAINS: define layout and decomposition            !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! <SUBROUTINE NAME="mpp_define_domains1D" INTERFACE="mpp_define_domains">>
  !   <IN NAME="global_indices" TYPE="integer" DIM="(2)"> </IN>
  !   <IN NAME="ndivs" TYPE="integer">  </IN>
  !   <INOUT NAME="domain" TYPE="type(domain1D)"> </INOUT>
  !   <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="flags" TYPE="integer">  </IN>
  !   <IN NAME="halo" TYPE="integer">  </IN>
  !   <IN NAME="extent" TYPE="integer" DIM="(0:)">  </IN>
  !   <IN NAME="maskmap" TYPE="logical" DIM="(0:)"> </IN>
  ! </SUBROUTINE>
  subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap, &
                                   memory_size, begin_halo, end_halo )
    !routine to divide global array indices among domains, and assign domains to PEs
    !domain is of type domain1D
    !ARGUMENTS:
    !      global_indices(2)=(isg,ieg) gives the extent of global domain
    !      ndivs is number of divisions of domain: even divisions unless extent is present.
    !      domain is the returned domain1D
    !      pelist (optional) list of PEs to which domains are to be assigned (default 0...npes-1)
    !                 size of pelist must correspond to number of mask=.TRUE. divisions
    !      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
    !      halo (optional) defines halo width (currently the same on both sides)
    !      extent (optional) array defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
    !      maskmap (optional) a division whose maskmap=.FALSE. is not assigned to any domain
    !  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
    !  as the input parameters permit
    integer, intent(in) :: global_indices(:) !(/ isg, ieg /)
    integer, intent(in) :: ndivs
    type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: flags, halo
    integer, intent(in), optional :: extent(0:)
    logical, intent(in), optional :: maskmap(0:)
    integer, intent(in), optional :: memory_size
    integer, intent(in), optional :: begin_halo, end_halo

    logical :: compute_domain_is_global, data_domain_is_global
    integer :: ndiv, n, isg, ieg, is, ie, i
    integer, allocatable :: pes(:)
    logical, allocatable :: mask(:)
    integer :: halosz, halobegin, haloend
    !used by symmetry algorithm
    integer :: imax, ndmax, ndmirror, errunit
    logical :: symmetrize, use_extent
    !statement functions
    logical :: even, odd
    even(n) = (mod(n,2).EQ.0)
    odd (n) = (mod(n,2).EQ.1)

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
    if(size(global_indices(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_domains1D: size of global_indices should be 2")
    !get global indices
    isg = global_indices(1)
    ieg = global_indices(2)
    if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
    !get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.mpp_pe()) )then
          errunit = stderr()
          write( errunit,* )'pe=', mpp_pe(), ' pelist=', pelist
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
       end if
       allocate( pes(0:size(pelist(:))-1) )
       pes(:) = pelist(:)
    else
       allocate( pes(0:mpp_npes()-1) )
       pes(:) = (/ (i,i=0,mpp_npes()-1) /)
    end if

    !get number of real domains: 1 mask domain per PE in pes
    allocate( mask(0:ndivs-1) )
    mask = .TRUE.                 !default mask
    if( PRESENT(maskmap) )then
       if( size(maskmap(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: maskmap array size must equal number of domain divisions.' )
       mask(:) = maskmap(:)
    end if
    if( count(mask).NE.size(pes(:)) ) &
         call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: number of TRUEs in maskmap array must match PE count.' )
    use_extent = .false.
    if( PRESENT(extent) )then
       if( size(extent(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
       use_extent = .true.
       if(ALL(extent ==0)) use_extent = .false.
    end if
    !get halosize
    halosz = 0
    if( PRESENT(halo) ) then
       halosz = halo
       !--- if halo is present, begin_halo and end_halo should not present
       if(present(begin_halo) .OR. present(end_halo) ) call mpp_error(FATAL, &
            "mpp_domains_define.inc: when halo is present, begin_halo and end_halo should not present")
    end if
    halobegin = halosz; haloend = halosz
    if(present(begin_halo)) halobegin = begin_halo
    if(present(end_halo))   haloend   = end_halo
    halosz = max(halobegin, haloend)
    !get flags
    compute_domain_is_global = .FALSE.
    data_domain_is_global    = .FALSE.
    domain%cyclic = .FALSE.
    if( PRESENT(flags) )then
       !NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
       compute_domain_is_global = ndivs.EQ.1
       !if compute domain is global, data domain must also be
       data_domain_is_global    = BTEST(flags,GLOBAL) .OR. compute_domain_is_global
       domain%cyclic  = BTEST(flags,CYCLIC) .AND. halosz.NE.0
    end if

    !set up links list
    allocate( domain%list(0:ndivs-1) )

    !set global domain
    domain%list(:)%global%begin     = isg
    domain%list(:)%global%end       = ieg
    domain%list(:)%global%size      = ieg-isg+1
    domain%list(:)%global%max_size  = ieg-isg+1
    domain%list(:)%global%is_global = .TRUE. !always

    !get compute domain
    if( compute_domain_is_global )then
       domain%list(:)%compute%begin = isg
       domain%list(:)%compute%end   = ieg
       domain%list(:)%compute%is_global = .TRUE.
       domain%list(:)%pe = pes(:)
       domain%pos = 0
    else
       domain%list(:)%compute%is_global = .FALSE.
       is = isg
       n = 0
       do ndiv=0,ndivs-1
          if( use_extent )then
             ie = is + extent(ndiv) - 1
             if( ndiv.EQ.ndivs-1 .AND. ie.NE.ieg ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extent array limits do not match global domain.' )
          else
             !modified for mirror-symmetry
             !original line
             !                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

             !problem of dividing nx points into n domains maintaining symmetry
             !i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
             !this will always work for nx even n even or odd
             !this will always work for nx odd, n odd
             !this will never  work for nx odd, n even: for this case we supersede the mirror calculation
             !                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
             !nx even n odd fails if n>nx/2
             symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                  (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

             !mirror domains are stored in the list and retrieved if required.
             if( ndiv.EQ.0 )then
                !initialize max points and max domains
                imax = ieg
                ndmax = ndivs
             end if
             !do bottom half of decomposition, going over the midpoint for odd ndivs
             if( ndiv.LT.(ndivs-1)/2+1 )then
                !domain is sized by dividing remaining points by remaining domains
                ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                ndmirror = (ndivs-1) - ndiv !mirror domain
                if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
                   !mirror extents, the max(,) is to eliminate overlaps
                   domain%list(ndmirror)%compute%begin = max( isg+ieg-ie, ie+1 )
                   domain%list(ndmirror)%compute%end   = max( isg+ieg-is, ie+1 )
                   imax = domain%list(ndmirror)%compute%begin - 1
                   ndmax = ndmax - 1
                end if
             else
                if( symmetrize )then
                   !do top half of decomposition by retrieving saved values
                   is = domain%list(ndiv)%compute%begin
                   ie = domain%list(ndiv)%compute%end
                else
                   ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                end if
             end if
          end if
          if( ie.LT.is )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents must be positive definite.' )
          domain%list(ndiv)%compute%begin = is
          domain%list(ndiv)%compute%end   = ie
          if( ndiv.GT.0 ) then
            if( is.NE.domain%list(ndiv-1)%compute%end+1 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
          endif
          if( ndiv.EQ.ndivs-1 .AND. domain%list(ndiv)%compute%end.NE.ieg ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
          if( mask(ndiv) )then
             domain%list(ndiv)%pe = pes(n)
             if( mpp_pe().EQ.pes(n) )domain%pos = ndiv
             n = n + 1
          else
             domain%list(ndiv)%pe = NULL_PE
          end if
          is = ie + 1
       end do
    end if

    domain%list(:)%compute%size  = domain%list(:)%compute%end - domain%list(:)%compute%begin + 1

    !get data domain
    !data domain is at least equal to compute domain
    domain%list(:)%data%begin = domain%list(:)%compute%begin
    domain%list(:)%data%end   = domain%list(:)%compute%end
    domain%list(:)%data%is_global = .FALSE.
    !apply global flags
    if( data_domain_is_global )then
       domain%list(:)%data%begin  = isg
       domain%list(:)%data%end    = ieg
       domain%list(:)%data%is_global = .TRUE.
    end if
    !apply margins
    domain%list(:)%data%begin = domain%list(:)%data%begin - halobegin
    domain%list(:)%data%end   = domain%list(:)%data%end   + haloend
    domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

    !--- define memory domain, if memory_size is not present or memory size is 0, memory domain size
    !--- will be the same as data domain size. if momory_size is present, memory_size should greater than
    !--- or equal to data size. The begin of memory domain will be always the same as data domain.
    domain%list(:)%memory%begin = domain%list(:)%data%begin
    domain%list(:)%memory%end   = domain%list(:)%data%end
    if( present(memory_size) ) then
       if(memory_size > 0) then
          if( domain%list(domain%pos)%data%size > memory_size ) call mpp_error(FATAL, &
               "mpp_domains_define.inc: data domain size is larger than memory domain size on this pe")
          domain%list(:)%memory%end   = domain%list(:)%memory%begin + memory_size - 1
       end if
    end if
    domain%list(:)%memory%size = domain%list(:)%memory%end - domain%list(:)%memory%begin + 1
    domain%list(:)%memory%is_global = domain%list(:)%data%is_global

    domain%compute = domain%list(domain%pos)%compute
    domain%data    = domain%list(domain%pos)%data
    domain%global  = domain%list(domain%pos)%global
    domain%memory  = domain%list(domain%pos)%memory
    domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
    domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
    domain%global%max_size  = domain%global%size
    domain%memory%max_size  = domain%memory%size

    !PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
    deallocate( pes, mask )
    return

  end subroutine mpp_define_domains1D

  ! <SUBROUTINE NAME="mpp_define_domains2D" INTERFACE="mpp_define_domains">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"> </IN>
  !  <IN NAME="layout" TYPE="integer" DIM="(2)"></IN>
  !  <INOUT NAME="domain" TYPE="type(domain2D)"></INOUT>
  !  <IN NAME="pelist" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="xflags, yflags" TYPE="integer"></IN>
  !  <IN NAME="xhalo, yhalo" TYPE="integer"></IN>
  !  <IN NAME="xextent, yextent" TYPE="integer" DIM="(0:)"></IN>
  !  <IN NAME="maskmap" TYPE="logical" DIM="(:,:)"></IN>
  !  <IN NAME="name" TYPE="character(len=*)"></IN>
  ! </SUBROUTINE>
  subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, xflags, yflags,    &
         xhalo, yhalo, xextent, yextent, maskmap, name, symmetry,  memory_size,               &
         whalo, ehalo, shalo, nhalo, is_mosaic, tile_count, complete )
    !define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
    integer,           intent(in) :: global_indices(:) !(/ isg, ieg, jsg, jeg /)
    integer,           intent(in) :: layout(:)
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
    integer, intent(in), optional :: xextent(0:), yextent(0:)
    logical, intent(in), optional :: maskmap(0:,0:)
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: symmetry
    logical, intent(in), optional :: is_mosaic ! indicate if calling mpp_define_domains from mpp_define_mosaic.
    integer, intent(in), optional :: memory_size(:)
    integer, intent(in), optional :: whalo, ehalo, shalo, nhalo ! halo size for West, East, South and North direction.
                                                                ! if whalo and ehalo is not present, 
                                                                ! will take the value of xhalo
                                                                ! if shalo and nhalo is not present, 
                                                                ! will take the value of yhalo
    integer, intent(in), optional :: tile_count                ! tile number on current pe, default value is 1
                                                                ! this is for the situation that multiple tiles on one processor.
    logical, intent(in), optional :: complete                   ! true indicate mpp_define_domain is completed for mosaic definition.    

    integer :: i, j, m, n, nn, xhalosz, yhalosz, memory_xsize, memory_ysize
    integer :: whalosz, ehalosz, shalosz, nhalosz
    integer :: ipos, jpos, pos, tile, ntile, msgsize, nlist
    integer :: ndivx, ndivy, isg, ieg, jsg, jeg, ishift, jshift, errunit, logunit
    logical :: from_mosaic, is_complete
    integer(LONG_KIND),save :: domain_cnt=0
    logical,    allocatable :: mask(:,:)
    integer,    allocatable :: pes(:), pesall(:), pearray(:,:), send_buf(:), recv_buf(:)
    character(len=8)        :: text

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
    if(size(global_indices(:)) .NE. 4) call mpp_error(FATAL,"mpp_define_domains2D: size of global_indices should be 4")
    if(size(layout(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_domains2D: size of layout should be 2")

    ndivx = layout(1); ndivy = layout(2)
    isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

    from_mosaic = .false.
    if(present(is_mosaic)) from_mosaic = is_mosaic
    is_complete = .true.
    if(present(complete)) is_complete = complete
    tile = 1
    if(present(tile_count)) tile = tile_count
  
    if( PRESENT(pelist) )then
       allocate( pes(0:size(pelist(:))-1) )
       pes = pelist
       if(from_mosaic) then
          allocate( pesall(0:mpp_npes()-1) )
          call mpp_get_current_pelist(pesall)   
       else
          allocate( pesall(0:size(pes(:))-1) )
          pesall = pes
       end if
    else
       allocate( pes(0:mpp_npes()-1) )
       allocate( pesall(0:mpp_npes()-1) )
       call mpp_get_current_pelist(pes)
       pesall = pes
    end if

    !--- when there is more than one tile on one processor, all the tile will limited on this processor
    if( tile > 1 .AND. size(pes(:)) > 1) call mpp_error(FATAL, &
         'MPP_DEFINE_DOMAINS2D: there are more than one tile on this pe, all the tile should be limited on this pe')

    !--- the position of current pe is changed due to mosaic, because  pes 
    !--- is only part of the pelist in mosaic (pesall). We assume the pe 
    !--- distribution are contious in mosaic.
    pos = mpp_pe() - mpp_root_pe()

    domain%symmetry = .FALSE.
    if(present(symmetry)) domain%symmetry = symmetry
    if(domain%symmetry) then
       ishift = 1; jshift = 1
    else
       ishift = 0; jshift = 0
    end if

    !Considering mosaic, the following will only be done on the pe in the pelist
    !when there is only one tile, all the current pe will be in the pelist.
    if( ANY(pes == mpp_pe()) ) then
       allocate( mask(0:ndivx-1,0:ndivy-1) )
       mask = .TRUE.
       if( PRESENT(maskmap) )then
          if( size(maskmap,1).NE.ndivx .OR. size(maskmap,2).NE.ndivy ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: maskmap array does not match layout.' )
          mask(:,:) = maskmap(:,:)
       end if
       !number of unmask domains in layout must equal number of PEs assigned
       n = count(mask)
       if( n.NE.size(pes(:)) )then
          write( text,'(i8)' )n
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for ' // &
               'this layout and maskmap. Use '//text//' PEs for this domain decomposition.' )
       end if

       !place on PE array; need flag to assign them to j first and then i
       allocate( pearray(0:ndivx-1,0:ndivy-1) )
       pearray(:,:) = NULL_PE
       ipos = NULL_PE; jpos = NULL_PE
       n = 0
       do j = 0,ndivy-1
          do i = 0,ndivx-1
             if( mask(i,j) )then
                pearray(i,j) = pes(n)
                if( pes(n).EQ.mpp_pe() )then
                   ipos = i
                   jpos = j
                end if
                n = n + 1
             end if
          end do
       end do
       if( ipos.EQ.NULL_PE .OR. jpos.EQ.NULL_PE ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: pelist must include this PE.' )
       if( debug ) then
         errunit = stderr()
         write( errunit, * )'pe, tile, ipos, jpos=', mpp_pe(), tile, ipos, jpos, ' pearray(:,jpos)=', &
                  pearray(:,jpos), ' pearray(ipos,:)=', pearray(ipos,:)
       endif
       !--- when tile is not equal to 1, the layout for that tile always ( 1, 1), so no need for pearray in domain
       if( tile == 1 ) then
          allocate( domain%pearray(0:ndivx-1,0:ndivy-1) )
          domain%pearray = pearray
       end if

       memory_xsize = 0; memory_ysize = 0
       if(present(memory_size)) then
          if(size(memory_size(:)) .NE. 2) call mpp_error(FATAL,"mpp_define_domains2D: size of memory_size should be 2")
          memory_xsize = memory_size(1)
          memory_ysize = memory_size(2)
       end if

       domain%pe  = mpp_pe()
       domain%pos  = pos
       if( .NOT. from_mosaic ) then
          allocate(domain%x(1), domain%y(1) )
          allocate(domain%tile_id(1))
          domain%tile_id        = 1
          domain%ntiles         = 1
          domain%max_ntile_pe   = 1
          domain%ncontacts      = 0
          domain%rotated_ninety = .FALSE.
       end if
       domain_cnt = domain_cnt + INT(1,KIND=LONG_KIND)
       domain%id = domain_cnt*DOMAIN_ID_BASE  ! Must be LONG_KIND arithmetic

       !do domain decomposition using 1D versions in X and Y, 
       call mpp_define_domains( global_indices(1:2), ndivx, domain%x(tile), &
            pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos), memory_xsize, whalo, ehalo )
       call mpp_define_domains( global_indices(3:4), ndivy, domain%y(tile), &
            pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:), memory_ysize, shalo, nhalo  )
       if( domain%x(tile)%list(ipos)%pe.NE.domain%y(tile)%list(jpos)%pe ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
       deallocate( mask, pearray )

       !set up fold
       domain%fold = 0
       if( PRESENT(xflags) )then
          if( BTEST(xflags,WEST) )domain%fold = domain%fold + FOLD_WEST_EDGE
          if( BTEST(xflags,EAST) )domain%fold = domain%fold + FOLD_EAST_EDGE
       end if
       if( PRESENT(yflags) )then
          if( BTEST(yflags,SOUTH) )domain%fold = domain%fold + FOLD_SOUTH_EDGE
          if( BTEST(yflags,NORTH) )domain%fold = domain%fold + FOLD_NORTH_EDGE
       end if

       !need to add code for EWS fold boundaries
       if( BTEST(domain%fold,WEST) .or. BTEST(domain%fold,EAST) .or. BTEST(domain%fold,SOUTH) )  &
            call mpp_error( FATAL, ' WEST, EAST and SOUTH fold are not implemented yet, contact author.' )

       if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,NORTH) )then
          if( domain%y(tile)%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%x(tile)%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in X must be even when there is a fold in Y.' )
          !check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x(tile)%list(i)%compute%size.NE.domain%x(tile)%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries ' // &
                                         'must line up (mirror-symmetric extents).' )
          end do
       end if
       if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
          if( domain%x(tile)%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%y(tile)%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in Y must be even when there is a fold in X.' )
          !check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y(tile)%list(i)%compute%size.NE.domain%y(tile)%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
       end if

       !set up domain%list
       if( mpp_pe().EQ.pes(0) .AND. PRESENT(name) )then
          logunit = stdlog()
          write( logunit, '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( logunit, '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
       end if
    end if   !     if( ANY(pes == mpp_pe()) ) 

    !--- set up domain%list.
    !--- set up 2-D domain decomposition for T, E, C, N and computing overlapping
    !--- when current tile is the last tile in the mosaic. 
    if(is_complete) then
       !--- set up domain%list
       call mpp_sync()
       ntile = size(domain%x(:))
       ! Because of mosaic, we need to replace mpp_transmit with mpp_send and mpp_recv.
       !--- pack the domain information and send to other pe.
       msgsize = ntile*19
       allocate(send_buf(msgsize))
       do n = 1, ntile
          nn = (n-1)*19
          call mpp_get_compute_domain( domain, send_buf(nn+1), send_buf(nn+2), send_buf(nn+3), send_buf(nn+4), tile_count = n )
          call mpp_get_data_domain   ( domain, send_buf(nn+5), send_buf(nn+6), send_buf(nn+7), send_buf(nn+8), tile_count = n )
          call mpp_get_global_domain ( domain, send_buf(nn+9), send_buf(nn+10), send_buf(nn+11), send_buf(nn+12), tile_count = n )
          call mpp_get_memory_domain ( domain, send_buf(nn+13), send_buf(nn+14), send_buf(nn+15), send_buf(nn+16) )
          send_buf(nn+17) = domain%tile_id(n)
          send_buf(nn+18) = domain%x(n)%pos
          send_buf(nn+19) = domain%y(n)%pos
       end do

       !the following are only used for storage of remote compute and data domain info
       nlist = size(pesall(:))
       if(.NOT. from_mosaic ) then
          allocate( domain%list(0:nlist-1) )        
          do i = 0, nlist-1
             allocate( domain%list(i)%x(1), domain%list(i)%y(1), domain%list(i)%tile_id(1) )
          end do
       end if

       do i = 0, nlist-1
          call mpp_send(send_buf(1), plen=msgsize, to_pe=pesall(mod(pos+nlist-i,nlist)) )
       end do

       !--- receive and unpack the domain information from other processors.
       do i = 0,nlist-1
          m = mod(pos+i,nlist) 
          domain%list(m)%pe = pesall(m)
          ntile = size(domain%list(m)%x(:))
          msgsize = 19*ntile
          allocate(recv_buf(msgsize))
          call mpp_recv(recv_buf(1), glen=msgsize, from_pe=pesall(m) ) 
          do n = 1, ntile
             nn = (n-1)*19
             domain%list(m)%x(n)%compute%begin = recv_buf(nn+1)
             domain%list(m)%x(n)%compute%end   = recv_buf(nn+2)
             domain%list(m)%y(n)%compute%begin = recv_buf(nn+3)
             domain%list(m)%y(n)%compute%end   = recv_buf(nn+4)
             domain%list(m)%x(n)%data%begin    = recv_buf(nn+5)
             domain%list(m)%x(n)%data%end      = recv_buf(nn+6)
             domain%list(m)%y(n)%data%begin    = recv_buf(nn+7)
             domain%list(m)%y(n)%data%end      = recv_buf(nn+8)
             domain%list(m)%x(n)%global%begin  = recv_buf(nn+9)
             domain%list(m)%x(n)%global%end    = recv_buf(nn+10)
             domain%list(m)%y(n)%global%begin  = recv_buf(nn+11)
             domain%list(m)%y(n)%global%end    = recv_buf(nn+12)
             domain%list(m)%x(n)%memory%begin  = recv_buf(nn+13)
             domain%list(m)%x(n)%memory%end    = recv_buf(nn+14)
             domain%list(m)%y(n)%memory%begin  = recv_buf(nn+15)
             domain%list(m)%y(n)%memory%end    = recv_buf(nn+16)
             domain%list(m)%tile_id(n)         = recv_buf(nn+17)
             domain%list(m)%x(n)%pos           = recv_buf(nn+18)
             domain%list(m)%y(n)%pos           = recv_buf(nn+19)
             domain%list(m)%pe                 = pesall(m)
          end do
          deallocate(recv_buf)
       end do

       call mpp_sync_self(pesall)
       deallocate(send_buf)

       do m = 0,nlist-1
          domain%list(m)%x(:)%compute%size = domain%list(m)%x(:)%compute%end - domain%list(m)%x(:)%compute%begin + 1
          domain%list(m)%y(:)%compute%size = domain%list(m)%y(:)%compute%end - domain%list(m)%y(:)%compute%begin + 1
          domain%list(m)%x(:)%data%size    = domain%list(m)%x(:)%data%end    - domain%list(m)%x(:)%data%begin + 1
          domain%list(m)%y(:)%data%size    = domain%list(m)%y(:)%data%end    - domain%list(m)%y(:)%data%begin + 1
          domain%list(m)%x(:)%global%size  = domain%list(m)%x(:)%global%end  - domain%list(m)%x(:)%global%begin + 1
          domain%list(m)%y(:)%global%size  = domain%list(m)%y(:)%global%end  - domain%list(m)%y(:)%global%begin + 1
          domain%list(m)%x(:)%memory%size  = domain%list(m)%x(:)%memory%end  - domain%list(m)%x(:)%memory%begin + 1
          domain%list(m)%y(:)%memory%size  = domain%list(m)%y(:)%memory%end  - domain%list(m)%y(:)%memory%begin + 1
       end do

       xhalosz = 0; yhalosz = 0
       if(present(xhalo)) xhalosz = xhalo
       if(present(yhalo)) yhalosz = yhalo
       whalosz = xhalosz; ehalosz = xhalosz
       shalosz = yhalosz; nhalosz = yhalosz
       if(present(whalo)) whalosz = whalo
       if(present(ehalo)) ehalosz = ehalo
       if(present(shalo)) shalosz = shalo
       if(present(nhalo)) nhalosz = nhalo
       domain%whalo = whalosz; domain%ehalo = ehalosz
       domain%shalo = shalosz; domain%nhalo = nhalosz

       allocate(domain%T, domain%E, domain%N, domain%C)
       call set_domain2d(domain, domain%T, CENTER, whalosz, ehalosz, shalosz, nhalosz, 0,           0)      
       call set_domain2d(domain, domain%C, CORNER, whalosz, ehalosz, shalosz, nhalosz, ishift, jshift)
       call set_domain2d(domain, domain%E, EAST,   whalosz, ehalosz, shalosz, nhalosz, ishift,      0)
       call set_domain2d(domain, domain%N, NORTH,  whalosz, ehalosz, shalosz, nhalosz,      0, jshift)
       call compute_overlaps(domain%T,      0,      0)
       call compute_overlaps(domain%C, ishift, jshift)
       call compute_overlaps(domain%E, ishift,      0)
       call compute_overlaps(domain%N,      0, jshift)
       !--- when ncontacts is nonzero, set_checkbound_point will be called in mpp_define
       if(domain%symmetry .AND. domain%ncontacts == 0 ) then
          allocate(domain%C%check, domain%E%check, domain%N%check )
          call set_checkbound( domain%C, domain%C%check )
          call set_checkbound( domain%E, domain%E%check )
          call set_checkbound( domain%N, domain%N%check )
          allocate(domain%C%bound, domain%E%bound, domain%N%bound)
          call set_boundary( domain%C, domain%C%bound )
          call set_boundary( domain%E, domain%E%bound )
          call set_boundary( domain%N, domain%N%bound )
       end if
    end if

    !print out decomposition, this didn't consider maskmap.
    if( mpp_pe() .EQ. pes(0) .AND. PRESENT(name) )then
       write(*,*) trim(name)//' domain decomposition'
       write (*,110) (domain%x(1)%list(i)%compute%size, i= 0, layout(1)-1)
       write (*,120) (domain%y(1)%list(i)%compute%size, i= 0, layout(2)-1)
110    format ('  X-AXIS = ',24i4,/,(11x,24i4))
120    format ('  Y-AXIS = ',24i4,/,(11x,24i4))
    endif

    deallocate( pes, pesall)

    return
  end subroutine mpp_define_domains2D

  !#####################################################################
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                    !
  !            set domain1d: set the 1-D domain decomposition          !
  !                                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_domain1d(domain_in, domain_out, shift, beginhalo, endhalo)
    type(domain1d), intent(in)    :: domain_in
    type(domain1d), intent(inout) :: domain_out
    integer,        intent(in)    :: shift
    integer,        intent(in)    :: beginhalo, endhalo

    integer :: n, ndivs
    integer :: bhalosz, ehalosz
  
    domain_out%cyclic  = domain_in%cyclic
    domain_out%pos     = domain_in%pos
    domain_out%shift   = shift
    domain_out%goffset = shift 
    domain_out%loffset = shift 

    ndivs = size(domain_in%list(:))
    allocate( domain_out%list(0:ndivs-1) )
    bhalosz = domain_in%compute%begin - domain_in%data%begin
    ehalosz = domain_in%data%end - domain_in%compute%end
    bhalosz = beginhalo + max(-sign(bhalosz, beginhalo), 0)
    ehalosz = endhalo   + max(-sign(ehalosz, endhalo  ), 0)

    do n = 0, ndivs-1
       domain_out%list(n)%pe                = domain_in%list(n)%pe
       domain_out%list(n)%global %begin     = domain_in%list(n)%global %begin 
       domain_out%list(n)%global %end       = domain_in%list(n)%global %end   + shift
       domain_out%list(n)%compute%begin     = domain_in%list(n)%compute%begin + min(beginhalo,0)
       domain_out%list(n)%compute%end       = domain_in%list(n)%compute%end   - min(endhalo,0)   + shift
       domain_out%list(n)%data   %begin     = domain_out%list(n)%compute%begin - bhalosz
       domain_out%list(n)%data   %end       = domain_out%list(n)%compute%end   + ehalosz
       domain_out%list(n)%memory %begin     = domain_in%list(n)%memory%begin
       domain_out%list(n)%memory %end       = domain_in%list(n)%memory%end   + shift
       domain_out%list(n)%global %size      = domain_out%list(n)%global %end - domain_out%list(n)%global %begin + 1
       domain_out%list(n)%compute%size      = domain_out%list(n)%compute%end - domain_out%list(n)%compute%begin + 1
       domain_out%list(n)%data   %size      = domain_out%list(n)%data   %end - domain_out%list(n)%data   %begin + 1
       domain_out%list(n)%memory %size      = domain_out%list(n)%memory %end - domain_out%list(n)%memory %begin + 1
       domain_out%list(n)%global %is_global = domain_in%list(n)%global %is_global
       domain_out%list(n)%compute%is_global = domain_in%list(n)%compute%is_global
       domain_out%list(n)%data   %is_global = domain_in%list(n)%data   %is_global
       domain_out%list(n)%memory %is_global = domain_in%list(n)%memory %is_global
   enddo

    domain_out%compute          = domain_out%list(domain_out%pos)%compute
    domain_out%data             = domain_out%list(domain_out%pos)%data
    domain_out%global           = domain_out%list(domain_out%pos)%global
    domain_out%memory           = domain_out%list(domain_out%pos)%memory
    domain_out%compute%max_size = MAXVAL( domain_out%list(:)%compute%size )
    domain_out%data%max_size    = MAXVAL( domain_out%list(:)%data%size )
    domain_out%global%max_size  = domain_out%global%size
    domain_out%memory%max_size  = domain_out%memory%size

    return

  end subroutine set_domain1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                    !
  !            set domain1d: set the 1-D domain decomposition          !
  !                                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_domain2d(domain_in, domain_out, position, whalo, ehalo, shalo, nhalo, ishift, jshift)
    type(domain2d), intent(in)    :: domain_in
    type(domain2d), intent(inout) :: domain_out
    integer,        intent(in)    :: position    
    integer,        intent(in)    :: whalo, ehalo, shalo, nhalo
    integer,        intent(in)    :: ishift, jshift
    integer                       :: ndivs, n, ntileMe, ntileNbr, tMe, tNbr, dir
    integer                       :: whalosz, ehalosz, shalosz, nhalosz
    type(domain1d), pointer       :: domx => NULL(), domy => NULL()

    whalosz = domain_in%whalo
    ehalosz = domain_in%ehalo
    shalosz = domain_in%shalo
    nhalosz = domain_in%nhalo
    domain_out%whalo = whalo
    domain_out%ehalo = ehalo
    domain_out%shalo = shalo
    domain_out%nhalo = nhalo
    if(whalo*ehalo < 0) call mpp_error(FATAL, &
        "mpp_domains_define.inc: whalo and ehalo can not have different sign")
    if(shalo*nhalo < 0) call mpp_error(FATAL, &
        "mpp_domains_define.inc: shalo and nhalo can not have different sign")
    if(whalo*shalo < 0) call mpp_error(FATAL, &
        "mpp_domains_define.inc: whalo and shalo can not have different sign")
    if( abs(whalo) > whalosz .OR. abs(ehalo) > ehalosz ) call mpp_error(FATAL, &
        "mpp_domains_define.inc: whalo and/or ehalo is larger than maximum halo size")
    if( abs(shalo) > shalosz .OR. abs(nhalo) > nhalosz ) call mpp_error(FATAL, &
        "mpp_domains_define.inc: shalo and/or nhalo is larger than maximum halo size")

    ntileMe = size(domain_in%x(:))
    ndivs = size(domain_in%list(:))
    !--- memory allocation
    allocate(domain_out%x(ntileMe), domain_out%y(ntileMe), domain_out%tile_id(ntileMe) )
    allocate(domain_out%rSpec(ntileMe) )
    allocate(domain_out%list(0:ndivs-1) )
    do n = 0, ndivs-1
       ntileNbr = size(domain_in%list(n)%x(:))
       allocate(domain_out%list(n)%x(ntileNbr),  domain_out%list(n)%y(ntileNbr) )
       allocate(domain_out%list(n)%tile_id(ntileNbr) )
       allocate(domain_out%list(n)%send(ntileNbr, ntileMe, 8), domain_out%list(n)%recv(ntileMe,ntileNbr, 8) )
    end do    

    do n = 1, ntileMe
       domain_out%rSpec(n)%count = 0 
       domain_out%rSpec(n)%total = 0
       call set_domain1d(domain_in%x(n), domain_out%x(n), ishift, whalo, ehalo)
       call set_domain1d(domain_in%y(n), domain_out%y(n), jshift, shalo, nhalo)
    end do

    domain_out%pos            = domain_in%pos
    domain_out%pe             = domain_in%pe
    domain_out%tile_id        = domain_in%tile_id
    domain_out%ntiles         = domain_in%ntiles
    domain_out%max_ntile_pe   = domain_in%max_ntile_pe
    domain_out%ncontacts      = domain_in%ncontacts
    domain_out%rotated_ninety = domain_in%rotated_ninety
    domain_out%position       = position
    domain_out%symmetry       = domain_in%symmetry
    domain_out%fold           = domain_in%fold

    whalosz = whalo + max(-sign(whalosz, whalo), 0)
    ehalosz = ehalo + max(-sign(ehalosz, ehalo), 0)
    shalosz = shalo + max(-sign(shalosz, shalo), 0)
    nhalosz = nhalo + max(-sign(nhalosz, nhalo), 0)

    do n = 0, ndivs-1
       ntileNbr = size(domain_in%list(n)%x(:))
       do dir = 1, 8
          do tMe = 1, ntileMe
             do tNbr = 1, ntileNbr
                domain_out%list(n)%send(tNbr, tMe, dir)%folded      = .FALSE.
                domain_out%list(n)%recv(tMe, tNbr, dir)%folded      = .FALSE.
                domain_out%list(n)%send(tNbr, tMe, dir)%rotation    = ZERO   
                domain_out%list(n)%recv(tMe, tNbr, dir)%rotation    = ZERO
                domain_out%list(n)%send(tNbr, tMe, dir)%overlap(:)  = .FALSE.
                domain_out%list(n)%recv(tMe, tNbr, dir)%overlap(:)  = .FALSE.
                domain_out%list(n)%send(tNbr, tMe, dir)%is(:)       = 0
                domain_out%list(n)%send(tNbr, tMe, dir)%ie(:)       = -1 
                domain_out%list(n)%send(tNbr, tMe, dir)%js(:)       = 0
                domain_out%list(n)%send(tNbr, tMe, dir)%je(:)       = -1 
                domain_out%list(n)%recv(tMe, tNbr, dir)%is(:)       = 0
                domain_out%list(n)%recv(tMe, tNbr, dir)%ie(:)       = -1 
                domain_out%list(n)%recv(tMe, tNbr, dir)%js(:)       = 0
                domain_out%list(n)%recv(tMe, tNbr, dir)%je(:)       = -1 
                domain_out%list(n)%send(tNbr, tMe, dir)%isMe        = 0
                domain_out%list(n)%send(tNbr, tMe, dir)%ieMe        = -1 
                domain_out%list(n)%send(tNbr, tMe, dir)%jsMe        = 0
                domain_out%list(n)%send(tNbr, tMe, dir)%jeMe        = -1 
                domain_out%list(n)%recv(tMe, tNbr, dir)%isMe        = 0
                domain_out%list(n)%recv(tMe, tNbr, dir)%ieMe        = -1 
                domain_out%list(n)%recv(tMe, tNbr, dir)%jsMe        = 0
                domain_out%list(n)%recv(tMe, tNbr, dir)%jeMe        = -1 
                domain_out%list(n)%send(tNbr, tMe, dir)%n           = 0
                domain_out%list(n)%recv(tMe, tNbr, dir)%n           = 0
                domain_out%list(n)%send(tNbr, tMe, dir)%is_refined  = .false.
                domain_out%list(n)%recv(tMe, tNbr, dir)%is_refined  = .false.
             end do
          end do
       end do
       domain_out%list(n)%overlap          = .FALSE.
       domain_out%list(n)%pe               = domain_in%list(n)%pe
       domain_out%list(n)%tile_id          = domain_in%list(n)%tile_id
       do tNbr = 1, ntileNbr
          domx => domain_out%list(n)%x(tNbr)
          domy => domain_out%list(n)%y(tNbr)
          domx%compute%begin  = domain_in%list(n)%x(tNbr)%compute%begin  + min(whalo,0)
          domx%compute%end    = domain_in%list(n)%x(tNbr)%compute%end    - min(ehalo,0)  + ishift
          domy%compute%begin  = domain_in%list(n)%y(tNbr)%compute%begin  + min(shalo,0)
          domy%compute%end    = domain_in%list(n)%y(tNbr)%compute%end    - min(nhalo,0)  + jshift
          domx%global%begin   = domain_in%list(n)%x(tNbr)%global%begin
          domx%global%end     = domain_in%list(n)%x(tNbr)%global%end     + ishift
          domy%global%begin   = domain_in%list(n)%y(tNbr)%global%begin   
          domy%global%end     = domain_in%list(n)%y(tNbr)%global%end     + jshift
          domx%memory%begin   = domain_in%list(n)%x(tNbr)%memory%begin   
          domx%memory%end     = domain_in%list(n)%x(tNbr)%memory%end     + ishift
          domy%memory%begin   = domain_in%list(n)%y(tNbr)%memory%begin   
          domy%memory%end     = domain_in%list(n)%y(tNbr)%memory%end     + jshift
          domx%data%begin     = domx%compute%begin - whalosz
          domx%data%end       = domx%compute%end   + ehalosz
          domy%data%begin     = domy%compute%begin - shalosz
          domy%data%end       = domy%compute%end   + nhalosz
          domx%compute%size   = domx%compute%end - domx%compute%begin + 1
          domy%compute%size   = domy%compute%end - domy%compute%begin + 1
          domx%data%size      = domx%data%end    - domx%data%begin    + 1
          domy%data%size      = domy%data%end    - domy%data%begin    + 1
          domx%global%size    = domx%global%end  - domx%global%begin + 1
          domy%global%size    = domy%global%end  - domy%global%begin + 1
          domx%memory%size    = domx%memory%end  - domx%memory%begin + 1
          domy%memory%size    = domy%memory%end  - domy%memory%begin + 1
       end do
   enddo

   return

  end subroutine set_domain2d

  !#####################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_define_mosaic: define mosaic domain                        !
!  NOTE: xflags and yflags is not in mpp_define_mosaic, because such relation !
!        are already defined in the mosaic relation.                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!??? do we need optional argument xextent and yextent
!??? how to specify pelist, we may use two dimensional variable pelist to represent.
!z1l: We assume the tilelist are in always limited to 1, 2, ... num_tile. If we want 
!     to remove this limitation, we need to add one more argument tilelist.
  subroutine mpp_define_mosaic( global_indices, layout, domain, num_tile, num_contact, tile1, tile2, &
                                istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, pe_start,      &
                                pe_end, pelist, whalo, ehalo, shalo, nhalo, xextent, yextent,                  &
                                maskmap, name, memory_size, symmetry, xflags, yflags )
    integer,          intent(in)           :: global_indices(:,:)  ! The size of first indice is 4, (/ isg, ieg, jsg, jeg /)
                                                                   ! The size of second indice is number of tiles in mosaic.
    integer,          intent(in)           :: layout(:,:)
    type(domain2D),   intent(inout)        :: domain
    integer,          intent(in)           :: num_tile             ! number of tiles in the mosaic
    integer,          intent(in)           :: num_contact          ! number of contact region between tiles.
    integer,          intent(in)           :: tile1(:), tile2(:)   ! tile number
    integer,          intent(in)           :: istart1(:), iend1(:) ! i-index in tile_1 of contact region
    integer,          intent(in)           :: jstart1(:), jend1(:) ! j-index in tile_1 of contact region
    integer,          intent(in)           :: istart2(:), iend2(:) ! i-index in tile_2 of contact region
    integer,          intent(in)           :: jstart2(:), jend2(:) ! j-index in tile_2 of contact region
    integer,          intent(in)           :: pe_start(:)          ! start pe of the pelist used in each tile
    integer,          intent(in)           :: pe_end(:)            ! end pe of the pelist used in each tile
    integer,          intent(in), optional :: pelist(:)            ! list of processors used in mosaic
    integer,          intent(in), optional :: whalo, ehalo, shalo, nhalo
    integer,          intent(in), optional :: xextent(:,:), yextent(:,:) 
    logical,          intent(in), optional :: maskmap(:,:,:)         
    character(len=*), intent(in), optional :: name
    integer,          intent(in), optional :: memory_size(2)
    logical,          intent(in), optional :: symmetry
    integer,          intent(in), optional :: xflags, yflags

    integer              :: n, m, ndivx, ndivy, list, nc, tile_count, nlist, nt, pos, n1, n2
    integer              :: whalosz, ehalosz, shalosz, nhalosz, xhalosz, yhalosz, t1, t2, tile
    integer              :: flags_x, flags_y
    logical, allocatable :: mask(:,:)
    integer, allocatable :: pes(:), xext(:), yext(:), pelist_tile(:), ntile_per_pe(:)
    logical              :: is_symmetry
    integer, allocatable :: align1(:), align2(:), is1(:), ie1(:), js1(:), je1(:), is2(:), ie2(:), js2(:), je2(:)
    integer, allocatable :: isgList(:), iegList(:), jsgList(:), jegList(:)
    real,    allocatable :: refine1(:), refine2(:)

    !--- the size of first indice of global_indices must be 4.
    if(size(global_indices, 1) .NE. 4) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of first dimension of global_indices is not 4')
    !--- the size of second indice of global_indices must be num_tile
    if(size(global_indices, 2) .NE. num_tile) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of second dimension of global_indices is not equal num_tile')    
    !--- the size of first indice of layout must be 2. The second dimension size of layout must equal num_tile.
    if(size(layout, 1) .NE. 2) call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of first dimension of layout is not 2')
    if(size(layout,2) .NE. num_tile)  call mpp_error(FATAL, &
         'mpp_domains_define.inc: The size of second dimension of layout is not equal num_tile')

    !--- setup pelist for the mosaic ---------------------
    nlist = mpp_npes()
    allocate(pes(0:nlist-1))
    if(present(pelist)) then
       if( nlist .NE. size(pelist(:))) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of pelist is not equal mpp_npes')
       pes = pelist
    else
       call mpp_get_current_pelist(pes)
    end if
    !--- pelist should be monotonic increasing by 1.
    do n = 1, nlist-1
       if(pes(n) - pes(n-1) .NE. 1) call mpp_error(FATAL, &
           'mpp_domains_define.inc: pelist is not monotonic increasing by 1')
    end do

    is_symmetry        = .FALSE.
    if(present(symmetry)) is_symmetry = symmetry

    if(size(pe_start(:)) .NE. num_tile .OR. size(pe_end(:)) .NE. num_tile ) call mpp_error(FATAL, &
         'mpp_domains_define.inc: size of pe_start and/or pe_end is not equal num_tile')
    !--- make sure pe_start and pe_end is in the pelist.   
    if( ANY( pe_start < pes(0) ) ) call mpp_error(FATAL, 'mpp_domains_define.inc: not all the pe_start are in the pelist')
    if( ANY( pe_end > pes(nlist-1)) ) call mpp_error(FATAL, 'mpp_domains_define.inc: not all the pe_end are in the pelist')

    !--- calculate number of tiles on each pe.
    allocate( ntile_per_pe(0:nlist-1) )
    ntile_per_pe = 0
    do n = 1, num_tile
       do m = pe_start(n) - mpp_root_pe(), pe_end(n) - mpp_root_pe()
          ntile_per_pe(m) = ntile_per_pe(m) + 1
       end do
    end do
    if(ANY(ntile_per_pe == 0)) call mpp_error(FATAL, &
         'mpp_domains_define.inc: At least one pe in pelist is not used by any tile in the mosaic')

    !--- check the size comformable of xextent and yextent
    if( PRESENT(xextent) ) then
       if(size(xextent,1) .GT. maxval(layout(1,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size mismatch between xextent and layout') 
       if(size(xextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of xextent is not eqaul num_tile')    
    end if
    if( PRESENT(yextent) ) then
       if(size(yextent,1) .GT. maxval(layout(2,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size mismatch between yextent and layout') 
       if(size(yextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: size of yextent is not eqaul num_tile')    
    end if

    !--- check the size comformable of maskmap
    !--- since the layout is different between tiles, so the actual size of maskmap for each tile is 
    !--- not diffrent. When define maskmap for multiple tiles, user can choose the maximum value
    !--- of layout of all tiles to the first and second dimension of maskmap. 
    if(present(maskmap)) then
       if(size(maskmap,1) .GT. maxval(layout(1,:)) .or. size(maskmap,2) .GT. maxval(layout(2,:))) &
            call mpp_error(FATAL, 'mpp_domains_define.inc: size mismatch between maskmap and layout')  
       if(size(maskmap,3) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define.inc: the third dimension of maskmap is not equal num_tile')  
    end if

    !--- define some mosaic information in domain type
    nt = ntile_per_pe(mpp_pe()-mpp_root_pe())
    allocate(domain%tile_id(nt), domain%x(nt), domain%y(nt) )
    allocate(domain%list(0:nlist-1))

    do n = 0, nlist-1
       nt = ntile_per_pe(mpp_root_pe()+n)
       allocate(domain%list(n)%x(nt), domain%list(n)%y(nt), domain%list(n)%tile_id(nt) )
    end do

    pe = mpp_pe()
    pos = 0
    do n = 1, num_tile
       if( pe .GE. pe_start(n) .AND. pe .LE. pe_end(n)) then
          pos = pos + 1
          domain%tile_id(pos) = n
       end if
    end do

    domain%rotated_ninety = .FALSE.
    domain%ntiles         = num_tile      
    domain%max_ntile_pe   = maxval(ntile_per_pe)
    domain%ncontacts      = num_contact
    
    !---call mpp_define_domain to define domain decomposition for each tile.
    tile_count = 0  ! tile number on current pe    
    do n = 1, num_tile
       allocate(mask(layout(1,n), layout(2,n)))
       allocate(pelist_tile(pe_start(n):pe_end(n)) )
       if( mpp_pe() .GE. pe_start(n) .AND. mpp_pe() .LE. pe_end(n) ) tile_count = tile_count + 1
       do m = pe_start(n), pe_end(n)
          pelist_tile(m) = m
       end do
       mask = .TRUE.
       if(present(maskmap))  mask = maskmap(1:layout(1,n), 1:layout(2,n), n)
       ndivx = layout(1,n); ndivy = layout(2,n)
       allocate(xext(ndivx), yext(ndivy))
       xext = 0; yext = 0
       if(present(xextent)) xext = xextent(1:ndivx,n)
       if(present(yextent)) yext = yextent(1:ndivy,n)
       ! when num_tile is one, we assume only folded_north and cyclic_x, cyclic_y boundary condition is the possible
       ! z1l: when we decide to support multiple-tile tripolar grid, we will redesign the following part.
       if(num_tile == 1) then
          flags_x = 0
          flags_y = 0
          if(PRESENT(xflags)) flags_x = xflags
          if(PRESENT(yflags)) flags_y = yflags
          do m = 1, num_contact
             if(istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
                if(istart2(m) .NE. iend2(m) ) call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, only cyclic condition is allowed for x-boundary")
                if(.NOT. BTEST(flags_x,CYCLIC))  flags_x = flags_x + CYCLIC_GLOBAL_DOMAIN
             else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
                if(jstart2(m) .NE. jend2(m) ) call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, only cyclic/folded-north condition is allowed for y-boundary")
                if( jstart1(m) == jstart2(m) ) then ! folded north
                   if(.NOT. BTEST(flags_y,NORTH) )  flags_y = flags_y + FOLD_NORTH_EDGE
                else 
                   if(.NOT. BTEST(flags_y,CYCLIC))  flags_y = flags_y + CYCLIC_GLOBAL_DOMAIN
                end if
             else 
               call mpp_error(FATAL,  &
                   "mpp_domains_define: for one tile mosaic, invalid boundary contact")
             end if
          end do
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pelist_tile, xflags = flags_x, &
                                  yflags = flags_y, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,           &
                                  xextent=xext, yextent=yext, maskmap=mask, name=name, symmetry=is_symmetry,      &
                                  memory_size = memory_size, is_mosaic = .true.)
       else
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pelist_tile,                   &
                                  whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, xextent=xext, yextent=yext, &
                                  maskmap=mask, name=name, symmetry=is_symmetry, memory_size = memory_size,       &
                                  is_mosaic = .true., tile_count = tile_count, complete = n==num_tile  )
       end if
       deallocate(mask, xext, yext, pelist_tile)
    end do

    deallocate(pes)

    if(num_contact == 0 .OR. num_tile == 1) return

    !--- loop through each contact region and find the contact for each tile ( including alignment )
    !--- we assume the tiles list is continuous and starting from 1.
    allocate(is1(num_contact), ie1(num_contact), js1(num_contact), je1(num_contact) )
    allocate(is2(num_contact), ie2(num_contact), js2(num_contact), je2(num_contact) )
    allocate(isgList(num_tile), iegList(num_tile), jsgList(num_tile), jegList(num_tile) )
    allocate(align1(num_contact), align2(num_contact), refine1(num_contact), refine2(num_contact))
    !--- get the global domain for each tile
    do n = 1, num_tile
       do list = 0, size(domain%list(:))-1
          do tile = 1, size(domain%list(list)%x(:))          
             if( domain%list(list)%tile_id(tile) == n ) then 
                isgList(n) = domain%list(list)%x(tile)%global%begin;  iegList(n) = domain%list(list)%x(tile)%global%end
                jsgList(n) = domain%list(list)%y(tile)%global%begin;  jegList(n) = domain%list(list)%y(tile)%global%end  
                exit
             end if
          end do
       end do
    end do

    !--- transfer the contact index to domain index.
    nc = 0
    do n = 1, num_contact
       t1 = tile1(n)
       t2 = tile2(n)
       is1(n) = istart1(n) + isgList(t1) - 1; ie1(n) = iend1(n) + isgList(t1) - 1
       js1(n) = jstart1(n) + jsgList(t1) - 1; je1(n) = jend1(n) + jsgList(t1) - 1
       is2(n) = istart2(n) + isgList(t2) - 1; ie2(n) = iend2(n) + isgList(t2) - 1
       js2(n) = jstart2(n) + jsgList(t2) - 1; je2(n) = jend2(n) + jsgList(t2) - 1
       call check_alignment( is1(n), ie1(n), js1(n), je1(n), isgList(t1), iegList(t1), jsgList(t1), jegList(t1), align1(n))
       call check_alignment( is2(n), ie2(n), js2(n), je2(n), isgList(t2), iegList(t2), jsgList(t2), jegList(t2), align2(n))
       if( (align1(n) == WEST .or. align1(n) == EAST ) .NEQV. (align2(n) == WEST .or. align2(n) == EAST ) )&
             domain%rotated_ninety=.true.   
    end do

    !--- calculate the refinement ratio between tiles
    do n = 1, num_contact
       n1 = max(abs(iend1(n) - istart1(n)), abs(jend1(n) - jstart1(n)) ) + 1
       n2 = max(abs(iend2(n) - istart2(n)), abs(jend2(n) - jstart2(n)) ) + 1
       refine1(n) = real(n2)/n1
       refine2(n) = real(n1)/n2
    end do    

    whalosz = 0; ehalosz = 0; shalosz = 0; nhalosz = 0
    if(present(whalo)) whalosz = whalo
    if(present(ehalo)) ehalosz = ehalo
    if(present(shalo)) shalosz = shalo
    if(present(nhalo)) nhalosz = nhalo
    xhalosz = max(whalosz, ehalosz)
    yhalosz = max(shalosz, nhalosz)

    domain%T%rotated_ninety = domain%rotated_ninety
    domain%E%rotated_ninety = domain%rotated_ninety
    domain%C%rotated_ninety = domain%rotated_ninety
    domain%N%rotated_ninety = domain%rotated_ninety
    !--- computing the overlap for the contact region with halo size xhalosz and yhalosz
    call define_contact_point( domain%T, num_contact, tile1, tile2, align1, align2, refine1, refine2, &
                               is1, ie1, js1, je1, is2, ie2, js2, je2, isgList, iegList, jsgList, jegList )

    call set_contact_point( domain%T, domain%C )
    call set_contact_point( domain%T, domain%E )
    call set_contact_point( domain%T, domain%N )

    !--- goffset setting is needed for exact global sum
    do m = 1, size(domain%tile_id(:))
       tile = domain%tile_id(m)
       do n = 1, num_contact
          if( tile1(n) == tile ) then
             if(align1(n) == EAST ) then
                domain%E%x(m)%goffset = 0
                domain%C%x(m)%goffset = 0
             end if
             if(align1(n) == NORTH ) then
                domain%N%y(m)%goffset = 0
                domain%C%y(m)%goffset = 0
             end if
          end if
          if( tile2(n) == tile ) then
             if(align2(n) == EAST ) then
                domain%E%x(m)%goffset = 0
                domain%C%x(m)%goffset = 0
             end if
             if(align2(n) == NORTH ) then
                domain%N%y(m)%goffset = 0
                domain%C%y(m)%goffset = 0
             end if
          end if
       end do
    end do

    !--- set the overlapping for boundary check if domain is symmetry
    if(domain%symmetry ) then
       allocate(domain%C%check, domain%E%check, domain%N%check )
       call set_checkbound( domain%C, domain%C%check )
       call set_checkbound( domain%E, domain%E%check )
       call set_checkbound( domain%N, domain%N%check )
       allocate(domain%C%bound, domain%E%bound, domain%N%bound)
       call set_boundary( domain%C, domain%C%bound )
       call set_boundary( domain%E, domain%E%bound )
       call set_boundary( domain%N, domain%N%bound )
    end if
    !--- release memory
    deallocate(align1, align2, is1, ie1, js1, je1, is2, ie2, js2, je2  )
    deallocate(isgList, iegList, jsgList, jegList )

    mosaic_defined = .true.
  end subroutine mpp_define_mosaic

!#####################################################################
  logical function mpp_mosaic_defined()
  ! Accessor function for value of mosaic_defined
    mpp_mosaic_defined = mosaic_defined
  end function mpp_mosaic_defined
!#####################################################################

    subroutine compute_overlaps( domain, ishift, jshift )
!computes remote domain overlaps
!assumes only one in each direction
!will calculate the overlapping for T,E,C,N-cell seperately.
      type(domain2D), intent(inout) :: domain
      integer, intent(in)           :: ishift, jshift
      integer :: i, j, m, nlist, tMe, tNbr
      integer :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff
      integer :: list, position, middle

      !--- since we restrict that if multiple tiles on one pe, all the tiles are limited to this pe. 
      !--- In this case, if ntiles on this pe is greater than 1, no overlapping between processor within each tile
      !--- In this case the overlapping exist only for tMe=1 and tNbr=1
      if(size(domain%x(:)) > 1) return

      position = domain%position
      !--- when there is only one tile, n will equal to np
      nlist = size(domain%list(:))

!send
      call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ioff, ysize=joff ) !cyclic offsets
      ioff = ioff - ishift
      joff = joff - jshift
      middle = (isg+ieg)/2+1
      tMe = 1; tNbr = 1

      do list = 0,nlist-1
         m = mod( domain%pos+list, nlist )
         if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
            !to_pe's eastern halo
            is = domain%list(m)%x(tNbr)%compute%end+1; ie = domain%list(m)%x(tNbr)%data%end
            js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end
            !--- to make sure the consistence between pes
            if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
                 .AND. ( jsc == je .or. jec == js ) ) then
               !--- do nothing, this point will come from other pe
            else 
               if( ie.GT.ieg )then
                  if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                     is = is-ioff; ie = ie-ioff
                  else if( BTEST(domain%fold,EAST) )then
                     i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                     j=js; js = jsg+jeg-je; je = jsg+jeg-j
                  end if
               end if
               !--- when the north face is folded, the east halo point at right side domain will be folded.
               !--- since folded-west is not implemented, we didn't consider it right now.
               !--- the position should be on CORNER or NORTH
               if( je == jeg .AND. BTEST(domain%fold,NORTH) .AND. (position == CORNER .OR. position == NORTH) & 
                    .AND. is .GE. middle .AND. domain%list(m)%x(tNbr)%data%end .LE. ieg ) then
                  call define_rect_overlap( domain%list(m)%send(tNbr,tMe,1), is, ie, js, je-1, isc, iec, jsc, jec)
                  is = domain%list(m)%x(tNbr)%compute%end+1; ie = domain%list(m)%x(tNbr)%data%end
                  js = je
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                  end select
                  call define_list_overlap( domain%list(m)%send(tNbr,tMe,1), is, ie, js, je, isc, iec, jsc, jec)
               else 
                  call define_rect_overlap( domain%list(m)%send(tNbr,tMe,1), is, ie, js, je, isc, iec, jsc, jec, &
                                            symmetry=domain%symmetry, direction=1)
               end if
            end if

!to_pe's SE halo
         is = domain%list(m)%x(tNbr)%compute%end+1; ie = domain%list(m)%x(tNbr)%data%end
         js = domain%list(m)%y(tNbr)%data%begin;    je = domain%list(m)%y(tNbr)%compute%begin-1
         if( ie.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,2), is, ie, js, je, isc, iec, jsc, jec)

!to_pe's southern halo
         is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end
         js = domain%list(m)%y(tNbr)%data%begin; je = domain%list(m)%y(tNbr)%compute%begin-1
         if( jsg.GT.js )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,3), is, ie, js, je, isc, iec, jsc, jec, &
                                            symmetry=domain%symmetry, direction=3)

!to_pe's SW halo
         is = domain%list(m)%x(tNbr)%data%begin; ie = domain%list(m)%x(tNbr)%compute%begin-1
         js = domain%list(m)%y(tNbr)%data%begin; je = domain%list(m)%y(tNbr)%compute%begin-1
         if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,4), is, ie, js, je, isc, iec, jsc, jec)

!to_pe's western halo
         is = domain%list(m)%x(tNbr)%data%begin;    ie = domain%list(m)%x(tNbr)%compute%begin-1
         js = domain%list(m)%y(tNbr)%compute%begin; je = domain%list(m)%y(tNbr)%compute%end
         if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            !--- when the north face is folded, some point at j=nj will be folded.
            !--- since folded-west is not implemented, we didn't consider it right now.
            !--- the position should be on CORNER or NORTH
            if( je == jeg .AND. BTEST(domain%fold,NORTH) .AND. (position == CORNER .OR. position == NORTH) &
                          .AND. ( (domain%x(tMe)%cyclic .AND. domain%list(m)%x(tNbr)%compute%begin == isg)            &
                           .OR. ( domain%list(m)%x(tNbr)%compute%begin-1 .GE. middle ) ) ) then
               call define_rect_overlap( domain%list(m)%send(tNbr,tMe,5), is, ie, js, je-1, isc, iec, jsc, jec)
               !--- consider at j = jeg for west edge.
               !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
               if(position == CORNER .AND. .NOT. domain%symmetry .AND. domain%list(m)%x(tNbr)%compute%begin == isg) then
                   call define_list_overlap( domain%list(m)%send(tNbr,tMe,5), ie, ie, je, je, isc, iec, jsc, jec)
               end if        
               is = domain%list(m)%x(tNbr)%data%begin;    ie = domain%list(m)%x(tNbr)%compute%begin-1
               je = domain%list(m)%y(tNbr)%compute%end;   js = je  
               if ( domain%list(m)%x(tNbr)%compute%begin == isg ) then
                  select case (position)
                  case(NORTH)
                     i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                  case(CORNER)
                     i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift  
                  end select
                  if(ie .GT. domain%x(tMe)%compute%end) call mpp_error( FATAL, &
                       'mpp_domains_util_mod: west edge ubound error send.' )     
               else 
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                  end select
               end if
               call define_list_overlap( domain%list(m)%send(tNbr,tMe,5), is, ie, js, je, isc, iec, jsc, jec)
            else
               call define_rect_overlap( domain%list(m)%send(tNbr,tMe,5), is, ie, js, je, isc, iec, jsc, jec, &
                                            symmetry=domain%symmetry, direction=5)
            end if            

!to_pe's NW halo
         is = domain%list(m)%x(tNbr)%data%begin;    ie = domain%list(m)%x(tNbr)%compute%begin-1
         js = domain%list(m)%y(tNbr)%compute%end+1; je = domain%list(m)%y(tNbr)%data%end
         if( isg.GT.is )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( je.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%send(tNbr,tMe,6)%folded = .TRUE.
                  select case ( position )
                  case(CENTER)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,6), is, ie, js, je, isc, iec, jsc, jec, &
                                      domain%list(m)%send(tNbr,tMe,6)%folded)
            !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%send(tNbr,tMe,6), is, is, js, je, isc, iec, jsc, jec)
         endif

!to_pe's northern halo
         is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end
         js = domain%list(m)%y(tNbr)%compute%end+1; je = domain%list(m)%y(tNbr)%data%end
         if( je.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%send(tNbr,tMe,7)%folded = .TRUE.
                  select case (position)
                  case(CENTER)
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
         !--- when domain symmetry and position is EAST or CORNER, the point when isc == ie, 
         !--- no need to send, because the data on that point will come from other pe. 
         !--- come from two pe ( there will be only one point on one pe. ).
         if( domain%symmetry .AND. (position == EAST .OR. position == CORNER ) &
              .AND. ( isc == ie .or. iec == is ) ) then
            !--- do nothing, this point will come from other pe
         else
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,7), is, ie, js, je, isc, iec, jsc, jec, &
                                      domain%list(m)%send(tNbr,tMe,7)%folded, symmetry=domain%symmetry, direction=7)
         end if
         !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
         if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
            is = is + ioff
            call define_list_overlap( domain%list(m)%send(tNbr,tMe,7), is, is, js, je, isc, iec, jsc, jec)
         endif

!to_pe's NE halo
         is = domain%list(m)%x(tNbr)%compute%end+1; ie = domain%list(m)%x(tNbr)%data%end
         js = domain%list(m)%y(tNbr)%compute%end+1; je = domain%list(m)%y(tNbr)%data%end
         if( ie.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if
         end if  
         if( je.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%send(tNbr,tMe,8)%folded = .TRUE.
                  select case(position)
                  case(CENTER)
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(tNbr,tMe,8), is, ie, js, je, isc, iec, jsc, jec, &
                                      domain%list(m)%send(tNbr,tMe,8)%folded)
         end if
      end do  ! end of send set up.

      ! begin setting up recv       
      do list = 0,nlist-1
         m = mod( domain%pos+nlist-list, nlist )
         if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe) ) then  ! only compute the overlapping within tile.
            call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
            !recv_e  
            isd = domain%x(tMe)%compute%end+1; ied = domain%x(tMe)%data%end
            jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end
            is=isc; ie=iec; js=jsc; je=jec
            if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
                 .AND. ( jsd == je .or. jed == js ) ) then
               ! --- do nothing, this point will come from other pe
            else
               if( ied.GT.ieg )then
                  if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                     is = is+ioff; ie = ie+ioff
                  else if( BTEST(domain%fold,EAST) )then
                     i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                     j=js; js = jsg+jeg-je; je = jsg+jeg-j
                  end if
               end if

               !--- when the north face is folded, the east halo point at right side domain will be folded.
               !--- since folded-west is not implemented, we didn't consider it right now.
               !--- the position should be on CORNER or NORTH
               if( jed == jeg .AND. BTEST(domain%fold,NORTH) .AND. (position == CORNER .OR. position == NORTH) &        
                    .AND. isd .GE. middle .AND. ied .LE. ieg ) then
                  call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,1), is, ie, js, je, isd, ied, jsd, jed-1)
                  is=isc; ie=iec; js=jsc; je=jec
                  jsd = jed
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                  end select
                  call define_list_overlap( domain%list(m)%recv(tMe,tNbr,1), is, ie, js, je, isd, ied, jsd, jed, .TRUE.)
               else 
                  call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,1), is, ie, js, je, isd, ied, jsd, jed, &
                                            symmetry=domain%symmetry, direction=1)
               end if
            end if
!recv_se     
         isd = domain%x(tMe)%compute%end+1; ied = domain%x(tMe)%data%end
         jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
         if( ied.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,2), is, ie, js, je, isd, ied, jsd, jed)

!recv_s      
         isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end
         jsd = domain%y(tMe)%data%begin;    jed = domain%y(tMe)%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
         call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,3), is, ie, js, je, isd, ied, jsd, jed, &
                                            symmetry=domain%symmetry, direction=3)

!recv_sw 
         isd = domain%x(tMe)%data%begin; ied = domain%x(tMe)%compute%begin-1
         jsd = domain%y(tMe)%data%begin; jed = domain%y(tMe)%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y(tMe)%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
         if( isd.LT.isg )then
             if( domain%x(tMe)%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,4), is, ie, js, je, isd, ied, jsd, jed)

!recv_w      
         isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
         jsd = domain%y(tMe)%compute%begin; jed = domain%y(tMe)%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         if( isd.LT.isg )then
             if( domain%x(tMe)%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            !--- when the north face is folded, some point at j=nj will be folded.
            !--- since folded-west is not implemented, we didn't consider it right now.
            !--- the position should be on CORNER or NORTH
            if( jed == jeg .AND. BTEST(domain%fold,NORTH) .AND. (position == CORNER .OR. position == NORTH) & 
                           .AND. ( (domain%x(tMe)%cyclic .AND. isd < isg ) .OR. ( ied .GE. middle ) ) )  then
               call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,5), is, ie, js, je, isd, ied, jsd, jed-1)
               !--- consider at j = jeg for west edge.
               !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
               if( position == CORNER .AND. .NOT. domain%symmetry .AND. isd < isg ) then
                   call define_list_overlap( domain%list(m)%recv(tMe,tNbr,5), is, ie, js, je, ied, ied, jed, jed)
               end if        
               is=isc; ie=iec; js=jsc; je=jec
               if(isd < isg) then
                  select case (position)
                  case(NORTH)
                     i=is; is = 2*isg-ie-1; ie = 2*isg-i-1     
                  case(CORNER)
                     ied = ied -1 + ishift
                     i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift 
                  end select
                  if(ie .GT. domain%x(tMe)%compute%end) call mpp_error( FATAL, &
                       'mpp_domains_util_mod: west edge ubound error recv.' )     
               else 
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                  end select
               end if
               call define_list_overlap( domain%list(m)%recv(tMe,tNbr,5), is, ie, js, je, isd, ied, jed, jed, .TRUE.)
            else
               call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,5), is, ie, js, je, isd, ied, jsd, jed, &
                                            symmetry=domain%symmetry, direction=5)
            end if            

!recv_nw     
         isd = domain%x(tMe)%data%begin;    ied = domain%x(tMe)%compute%begin-1
         jsd = domain%y(tMe)%compute%end+1; jed = domain%y(tMe)%data%end
         is=isc; ie=iec; js=jsc; je=jec
         if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(tMe,tNbr,6)%folded = .TRUE.
                  select case(position)
                  case(CENTER)
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
         if( isd.LT.isg )then
               if( domain%x(tMe)%cyclic .AND. is.GE.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if 
            call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,6), is, ie, js, je, isd, ied, jsd, jed, &
                                      domain%list(m)%recv(tMe,tNbr,6)%folded)
            !--- when north edge is folded, is will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%recv(tMe,tNbr,6), is, is, js, je, isd, ied, jsd, jed, &
                                         domain%list(m)%recv(tMe,tNbr,6)%folded )
         endif

!recv_n      
         isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end
         jsd = domain%y(tMe)%compute%end+1; jed = domain%y(tMe)%data%end
         is=isc; ie=iec; js=jsc; je=jec
         if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(tMe,tNbr,7)%folded = .TRUE.
                  select case(position)
                  case(CENTER)
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
         !--- when domain symmetry and position is EAST or CORNER, the point at i=isd will 
         !--- come from two pe ( there will be only one point on one pe. ).
         if( domain%symmetry .AND. (position == EAST .OR. position == CORNER ) &
              .AND. (isd == ie .or. ied == is ) ) then
            !--- do nothing, this point will come from other pe
         else
            call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,7), is, ie, js, je, isd, ied, jsd, jed, &
                                      domain%list(m)%recv(tMe,tNbr,7)%folded, symmetry=domain%symmetry, direction=7)
         end if
         !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
         if(is .LT. isg .AND. domain%x(tMe)%cyclic) then
            is = is + ioff
            call define_list_overlap( domain%list(m)%recv(tMe,tNbr,7), is, is, js, je, isd, ied, jsd, jed, &
                                      domain%list(m)%recv(tMe,tNbr,7)%folded)
         endif
!recv_ne     
         isd = domain%x(tMe)%compute%end+1; ied = domain%x(tMe)%data%end
         jsd = domain%y(tMe)%compute%end+1; jed = domain%y(tMe)%data%end
         is=isc; ie=iec; js=jsc; je=jec

         if( jed.GT.jeg )then
             if( domain%y(tMe)%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(tMe,tNbr,8)%folded = .TRUE.
                  select case(position)
                  case(CENTER)
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(EAST)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                     j=js; js = 2*jeg-je; je = 2*jeg-j
                  end select
             end if  
         end if  
         if( ied.GT.ieg )then
             if( domain%x(tMe)%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(tMe,tNbr,8), is, ie, js, je, isd, ied, jsd, jed, &
                       domain%list(m)%recv(tMe,tNbr,8)%folded)
         endif
      enddo ! end of recv do loop

      !--- Now calculate the overlapping for fold-edge. Currently we only consider about folded-north
      !--- for folded-north-edge, only need to consider to_pe's north(7) direction 
      !--- only position at NORTH and CORNER need to be considered
      if( BTEST(domain%fold,NORTH) .AND. ( position == NORTH .OR. position == CORNER) ) then
         if( domain%y(tMe)%data%begin .LE. jeg .AND. jeg .LE. domain%y(tMe)%data%end )then !fold is within domain
            !--- calculate the overlapping for sending
            if( domain%x(tMe)%pos .LT. (size(domain%x(tMe)%list(:))+1)/2 )then 
               do list = 0,nlist-1
                  m = mod( domain%pos+list, nlist )
                  if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe)) then  ! only compute the overlapping within tile.
                     call mpp_get_compute_domain( domain, isc, iec, jsc, jec )                     
                     js = domain%list(m)%y(tNbr)%compute%end;   je = js
                     if( je == jeg )then   ! fold is within domain.
                        is = domain%list(m)%x(tNbr)%compute%begin; ie = domain%list(m)%x(tNbr)%compute%end
                        select case (position)
                        case(NORTH)
                           is = max(is, middle)
                           i=is; is = isg+ieg-ie; ie = isg+ieg-i
                        case(CORNER)
                           is = max(is, middle)
                           i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                        end select
                        call define_list_overlap( domain%list(m)%send(tNbr,tMe,7), is, ie, js, je, isc, iec, jsc, jec)           
                     end if
                  end if
               end do
            end if

            !--- calculating overlapping for receving on north
            if( domain%x(tMe)%pos .GE. size(domain%x(tMe)%list(:))/2 )then 
               do list = 0,nlist-1
                  m = mod( domain%pos+list, nlist )
                  if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe)) then  ! only compute the overlapping within tile.
                     call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
                     jsd = domain%y(tMe)%compute%end;   jed = jsd
                     if( jed == jeg )then   ! fold is within domain.
                        isd = domain%x(tMe)%compute%begin; ied = domain%x(tMe)%compute%end
                        is=isc; ie=iec; js = jsc; je = jec
                        select case (position)
                        case(NORTH)
                           isd = max(isd, middle)
                           i=is; is = isg+ieg-ie; ie = isg+ieg-i
                        case(CORNER)
                           isd = max(isd, middle)
                           i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                        end select
                        call define_list_overlap( domain%list(m)%recv(tMe,tNbr,7), is, ie, js, je, isd, ied, jsd, jed, .TRUE.) 
                     end if
                  end if
               end do
            end if
         end if
      end if

      do m = 0,nlist-1
         if(domain%list(m)%tile_id(tNbr) == domain%tile_id(tMe)) then  
            if(ANY(domain%list(m)%send(:,:,:)%overlap(1))) domain%list(m)%overlap = .TRUE.
            if(ANY(domain%list(m)%send(:,:,:)%n > 0))      domain%list(m)%overlap = .TRUE.
            if(ANY(domain%list(m)%recv(:,:,:)%overlap(1))) domain%list(m)%overlap = .TRUE.
            if(ANY(domain%list(m)%recv(:,:,:)%n > 0))      domain%list(m)%overlap = .TRUE.
            if(ANY(domain%list(m)%recv(tMe,:,1)%overlap(1))) domain%x(tMe)%loffset = 0
            if(ANY(domain%list(m)%recv(tMe,:,7)%overlap(1))) domain%y(tMe)%loffset = 0
         end if
      end do
      if(domain%x(tMe)%cyclic) domain%x(tMe)%goffset = 0
      if(domain%y(tMe)%cyclic) domain%y(tMe)%goffset = 0

      domain%initialized = .true.

    end subroutine compute_overlaps

    !#####################################################################################
    ! The following two routines will be used in compute_overlaps
    subroutine define_rect_overlap(overlap, is1, ie1, js1, je1, is2, ie2, js2, je2, reverse, symmetry, direction)
      type(overlapSpec), intent(inout) :: overlap
      integer,           intent(in)    :: is1, ie1, js1, je1, is2, ie2, js2, je2
      logical, optional, intent(in)    :: reverse, symmetry
      integer, optional, intent(in)    :: direction

      logical :: is_reverse, is_symmetry
      integer :: is, ie, js, je, dir

      is_reverse = .FALSE.
      if(PRESENT(reverse)) is_reverse = reverse
      is_symmetry = .FALSE.
      if(PRESENT(symmetry)) is_symmetry = symmetry
      dir = 0
      if(PRESENT(direction)) dir = direction

      is = max(is1,is2); ie = min(ie1,ie2)
      js = max(js1,js2); je = min(je1,je2)

      !--- to avoid unnecessary ( duplicate overlap ) for symmetry domain
      if(is_symmetry .AND. (dir == 1 .OR. dir == 5)) then  ! x-direction
         if( ie .GE. is .AND. je .GT. js ) then
            overlap%overlap(1) = .TRUE.
            overlap%is(1) = is
            overlap%ie(1) = ie
            overlap%js(1) = js
            overlap%je(1) = je
            if(is_reverse) overlap%rotation = ONE_HUNDRED_EIGHTY
         end if
      else if(is_symmetry .AND. (dir == 3 .OR. dir == 7)) then ! y-direction
         if( ie .GT. is .AND. je .GE. js ) then
            overlap%overlap(1) = .TRUE.
            overlap%is(1) = is
            overlap%ie(1) = ie
            overlap%js(1) = js
            overlap%je(1) = je
            if(is_reverse) overlap%rotation = ONE_HUNDRED_EIGHTY
         end if
      else if(ie.GE.is .AND. je.GE.js )then

         overlap%overlap(1) = .TRUE.
         overlap%is(1) = is
         overlap%ie(1) = ie
         overlap%js(1) = js
         overlap%je(1) = je
         if(is_reverse) overlap%rotation = ONE_HUNDRED_EIGHTY
      endif

    end subroutine define_rect_overlap

    subroutine define_list_overlap(overlap, is1, ie1, js1, je1, is2, ie2, js2, je2, reverse)
      type(overlapSpec), intent(inout) :: overlap
      integer,         intent(in)   :: is1, ie1, js1, je1, is2, ie2, js2, je2
      logical, optional, intent(in)  :: reverse

      integer              :: is, ie, js, je, n, n_old, i, j
      integer, allocatable :: tmp(:)
      logical :: is_reverse

      is_reverse = .FALSE.
      if(PRESENT(reverse)) is_reverse = reverse

      is = max(is1,is2); ie = min(ie1,ie2)
      js = max(js1,js2); je = min(je1,je2)

      if(ie.GE.is .AND. je.GE.js )then
         n_old = overlap%n
         n = (ie-is+1)*(je-js+1) + n_old
         if( n_old > 0 ) then
            allocate(tmp(n_old))
            tmp = overlap%i
            deallocate(overlap%i)
            allocate( overlap%i(n) )
            overlap%i(1:n_old) = tmp
            tmp = overlap%j
            deallocate(overlap%j)
            allocate( overlap%j(n) )
            overlap%j(1:n_old) = tmp
            deallocate(tmp)
         else
            allocate( overlap%i(n) )
            allocate( overlap%j(n) )
         endif

         overlap%n = n
         n = n_old
         if(is_reverse) then  ! reverse when receiving
            do j = je, js, -1
               do i = ie, is, -1
                  n = n + 1
                  overlap%i(n) = i
                  overlap%j(n) = j
               end do
            end do
         else
            do j = js, je
               do i = is, ie
                  n = n + 1
                  overlap%i(n) = i
                  overlap%j(n) = j
               end do
            end do
         end if
      end if

      return
    end subroutine define_list_overlap


  !###################################################################################
  ! this routine setup the overlapping for mpp_update_domains for arbitrary halo update.
  ! should be the halo size defined in mpp_define_domains.
  ! xhalo_out, yhalo_out should not be exactly the same as xhalo_in, yhalo_in
  ! currently we didn't consider about tripolar grid situation, because in the folded north
  ! region, the overlapping is specified through list of points, not through rectangular.
  ! But will return back to solve this problem in the future.
  subroutine set_overlaps(domain_in, domain_out, whalo_out, ehalo_out, shalo_out, nhalo_out)
    type(domain2d),    intent(in) :: domain_in
    type(domain2d), intent(inout) :: domain_out
    integer,           intent(in) :: whalo_out, ehalo_out, shalo_out, nhalo_out
    integer                       :: nlist, list, m, isoff, ieoff, jsoff, jeoff, rotation
    integer                       :: whalo_in, ehalo_in, shalo_in, nhalo_in
    integer                       :: ntileMe, ntileNbr, tMe, tNbr, dir    
    type(overlapSpec), pointer    :: ptrIn  => NULL()
    type(overlapSpec), pointer    :: ptrOut => NULL()

    if( domain_in%fold .NE. 0) call mpp_error(FATAL, &
         "mpp_domains_define.inc: folded domain is not implemented for arbitrary halo update, contact developer")

    whalo_in = domain_in%whalo
    ehalo_in = domain_in%ehalo
    shalo_in = domain_in%shalo
    nhalo_in = domain_in%nhalo

    if( ASSOCIATED(domain_in%check) ) domain_out%check => domain_in%check

    if( .NOT. domain_in%initialized) call mpp_error(FATAL, &
         "mpp_domains_define.inc: domain_in is not defined yet")
    !--- first set up the domain
    call set_domain2d(domain_in, domain_out, domain_in%position, whalo_out, ehalo_out, shalo_out, nhalo_out, 0, 0)

    nlist = size(domain_in%list(:))
    isoff = whalo_in - abs(whalo_out)
    ieoff = ehalo_in - abs(ehalo_out)
    jsoff = shalo_in - abs(shalo_out)
    jeoff = nhalo_in - abs(nhalo_out)
    ntileMe = size(domain_in%x(:))
    !--- setting up overlap.
    do list = 0, nlist-1
       m = mod( domain_in%pos+list, nlist ) 
       domain_out%list(m)%overlap = domain_in%list(m)%overlap 
       if( .NOT. domain_in%list(m)%overlap ) cycle 
       ntileNbr = size(domain_in%list(m)%x(:))
       do dir = 1, 8  ! loop over 8 directions
          do tMe = 1, ntileMe
             do tNbr = 1, ntileNbr
                ptrIn  => domain_in %list(m)%send(tNbr,tMe,dir)
                ptrOut => domain_out%list(m)%send(tNbr,tMe,dir)
                if(ptrIn%overlap(1)) then
                   rotation = ptrIn%rotation
                   select case(dir)
                   case(1) ! to_pe's eastern halo
                      if(ehalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, 0, 0, rotation)
                      else if(ehalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, -ehalo_out, 0, 0, 0, rotation)
                      end if
                   case(2) ! to_pe's southeast halo
                      if(ehalo_out>0 .AND. shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, jsoff, 0, rotation)
                      else if(ehalo_out<0 .AND. shalo_out < 0) then ! three parts: southeast, south and east.
                         call set_single_overlap(ptrIn, ptrOut, -ehalo_out, 0, 0, shalo_out, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir-1), -ehalo_out, 0, jsoff, 0, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir+1), 0, -ieoff, 0, shalo_out, rotation)
                      end if
                   case(3) ! to_pe's southern halo
                      if(shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, jsoff, 0, rotation)
                      else if(shalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, 0, shalo_out, rotation)
                      end if
                   case(4) ! to_pe's southwest halo
                      if(whalo_out>0 .AND. shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, jsoff, 0, rotation)
                      else if(whalo_out<0 .AND. shalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, 0, shalo_out, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir-1), isoff, 0, 0, shalo_out, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir+1), 0, whalo_out, jsoff, 0, rotation)
                      end if
                   case(5) ! to_pe's western halo
                      if(whalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, 0, 0, rotation)
                      else if(whalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, 0, 0, rotation)
                      end if
                   case(6) ! to_pe's northwest halo
                      if(whalo_out>0 .AND. nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, 0, -jeoff, rotation)
                      else if(whalo_out<0 .AND. nhalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, -nhalo_out, 0, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir-1), 0, whalo_out, 0, -jeoff, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir+1), isoff, 0, -nhalo_out, 0, rotation)
                      end if
                   case(7) ! to_pe's northern halo
                      if(nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, 0, -jeoff, rotation)
                      else if(nhalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, -nhalo_out, 0, rotation)
                      end if
                   case(8) ! to_pe's northeast halo
                      if(ehalo_out>0 .AND. nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, 0, -jeoff, rotation)
                      else if(ehalo_out<0 .AND. nhalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, -ehalo_out, 0, -nhalo_out, 0, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,dir-1), 0, -ieoff, -nhalo_out, 0, rotation)
                         call set_single_overlap(ptrIn, domain_out%list(m)%send(tNbr,tMe,1), -ehalo_out, 0, 0, -jeoff, rotation)
                      end if
                   end select
                end if
             end do ! do tNbr = 1, ntileNbr
          end do ! end do tMe = 1, ntileMe
       end do ! do dir = 1, 8 
    end do ! end do list = 0, nlist-1

    !--------------------------------------------------
    !                     recving
    !---------------------------------------------------
    do list = 0, nlist-1
       m = mod( domain_in%pos+list, nlist )   
       if(.NOT. domain_in%list(m)%overlap ) cycle  
       ntileNbr = size(domain_in%list(m)%x(:))
       do dir = 1, 8  ! loop over 8 directions
          do tNbr = 1, ntileNbr
             do tMe = 1, ntileMe
                ptrIn  => domain_in %list(m)%recv(tMe,tNbr,dir)
                ptrOut => domain_out%list(m)%recv(tMe,tNbr,dir)   
                if(ptrIn%overlap(1)) then
                   select case(dir)
                   case(1) ! eastern halo
                      if(ehalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, 0, 0)
                      else if(ehalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, -ehalo_out, 0, 0, 0)
                      end if
                   case(2) ! southeast halo
                      if(ehalo_out>0 .AND. shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, jsoff, 0)
                      else if(ehalo_out<0 .AND. shalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, &
                              -ehalo_out, 0, 0, shalo_out)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir-1), -ehalo_out, 0, jsoff, 0)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir+1), 0, -ieoff, 0, shalo_out)
                      end if
                   case(3) ! southern halo
                      if(shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, jsoff, 0)
                      else if(shalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, 0, shalo_out)
                      end if
                   case(4) ! southwest halo
                      if(whalo_out>0 .AND. shalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, jsoff, 0)
                      else if(whalo_out<0 .AND. shalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, 0, shalo_out)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir-1), isoff, 0, 0, shalo_out)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir+1), 0, whalo_out, jsoff, 0)
                      end if
                   case(5) ! western halo
                      if(whalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, 0, 0)
                      else if(whalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, 0, 0)
                      end if
                   case(6) ! northwest halo
                      if(whalo_out>0 .AND. nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, isoff, 0, 0, -jeoff)
                      else if(whalo_out<0 .AND. nhalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, whalo_out, -nhalo_out, 0)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir-1), 0, whalo_out, 0, -jeoff)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir+1), isoff, 0, -nhalo_out, 0)
                      end if
                   case(7) ! northern halo
                      if(nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, 0, -jeoff)
                      else if(nhalo_out<0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, 0, -nhalo_out, 0)
                      end if
                   case(8) ! northeast halo
                      if(ehalo_out>0 .AND. nhalo_out > 0) then
                         call set_single_overlap(ptrIn, ptrOut, 0, -ieoff, 0, -jeoff)
                      else if(ehalo_out<0 .AND. nhalo_out < 0) then
                         call set_single_overlap(ptrIn, ptrOut, -ehalo_out, 0, -nhalo_out, 0)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,dir-1), 0, -ieoff, -nhalo_out, 0)
                         call set_single_overlap(ptrIn, domain_out%list(m)%recv(tMe,tNbr,1), -ehalo_out, 0, 0, -jeoff)
                      end if
                   end select
                end if
             end do ! end do tMe = 1, ntileMe
          end do ! do tNbr = 1, ntileNbr
       end do ! do dir = 1, 8 
    end do ! end do list = 0, nlist-1

    domain_out%initialized = .true.

    contains

    subroutine set_single_overlap(overlap_in, overlap_out, isoff, ieoff, jsoff, jeoff, rotation)
       type(overlapSpec),    intent(in) :: overlap_in
       type(overlapSpec), intent(inout) :: overlap_out
       integer,              intent(in) :: isoff, jsoff, ieoff, jeoff       
       integer, optional,    intent(in) :: rotation
       integer                          :: rotate
       integer                          :: level

       level = 1
       do level = 1, 3
          if( .NOT. overlap_out%overlap(level) ) exit
       end do
       if(level > 3) call mpp_error(FATAL,"mpp_domains_define.inc: level is greater than 3")

       rotate = ZERO
       overlap_out%overlap(level) = .true.
       if(present(rotation)) rotate = rotation
       overlap_out%rotation = overlap_in%rotation
       overlap_out%is_refined(level) = overlap_in%is_refined(1)

       select case(rotate)
       case(ZERO)
          overlap_out%is(level)      = overlap_in%is(1) + isoff
          overlap_out%ie(level)      = overlap_in%ie(1) + ieoff
          overlap_out%js(level)      = overlap_in%js(1) + jsoff
          overlap_out%je(level)      = overlap_in%je(1) + jeoff
       case(NINETY)
          overlap_out%is(level)      = overlap_in%is(1) - jeoff
          overlap_out%ie(level)      = overlap_in%ie(1) - jsoff
          overlap_out%js(level)      = overlap_in%js(1) + isoff
          overlap_out%je(level)      = overlap_in%je(1) + ieoff
       case(MINUS_NINETY)
          overlap_out%is(level)      = overlap_in%is(1) + jsoff
          overlap_out%ie(level)      = overlap_in%ie(1) + jeoff
          overlap_out%js(level)      = overlap_in%js(1) - ieoff
          overlap_out%je(level)      = overlap_in%je(1) - isoff
       case default
          call mpp_error(FATAL, "mpp_domains_define.inc: the value of rotation should be ZERO, NINETY or MINUS_NINETY")
       end select

    end subroutine set_single_overlap

  end subroutine set_overlaps

  !###################################################################################
  !--- compute the overlapping between tiles for the T-cell. 
  subroutine define_contact_point( domain, num_contact, tile1, tile2, align1, align2, refine1, refine2, istart1, iend1, &
                                   jstart1, jend1, istart2, iend2, jstart2, jend2,       &
                                   isgList, iegList, jsgList, jegList )
    type(domain2D),     intent(inout) :: domain
    integer,               intent(in) :: num_contact      ! number of contact regions
    integer, dimension(:), intent(in) :: tile1, tile2     ! tile number
    integer, dimension(:), intent(in) :: align1, align2   ! align direction of contact region
    real,    dimension(:), intent(in) :: refine1, refine2 ! refinement between tiles
    integer, dimension(:), intent(in) :: istart1, iend1   ! i-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: jstart1, jend1   ! j-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: istart2, iend2   ! i-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: jstart2, jend2   ! j-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: isgList, iegList ! i-global domain of each tile
    integer, dimension(:), intent(in) :: jsgList, jegList ! j-global domain of each tile

    type(contact_type), allocatable :: eCont(:), wCont(:), sCont(:), nCont(:)
    integer              :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer              :: isc1, iec1, jsc1, jec1, isc2, iec2, jsc2, jec2
    integer              :: isd1, ied1, jsd1, jed1, isd2, ied2, jsd2, jed2
    integer              :: is, ie, js, je, ioff, joff, isoff, ieoff, jsoff, jeoff
    integer              :: ntiles, max_contact
    integer              :: nlist, list, m, n, count, numS, numR
    integer              :: whalo, ehalo, shalo, nhalo
    integer              :: t1, t2
    integer              :: ntileMe, ntileNbr, tMe, tNbr, tileMe, dir
    integer              :: nxd, nyd, nxc, nyc
    integer, allocatable :: is1Send(:), ie1Send(:), js1Send(:), je1Send(:)
    integer, allocatable :: is2Send(:), ie2Send(:), js2Send(:), je2Send(:)
    integer, allocatable :: is2Recv(:), ie2Recv(:), js2Recv(:), je2Recv(:)
    integer, allocatable :: is1Recv(:), ie1Recv(:), js1Recv(:), je1Recv(:)
    integer, allocatable :: align1Recv(:), align2Recv(:), align1Send(:), align2Send(:)
    real,    allocatable :: refineRecv(:), refineSend(:)
    integer, allocatable :: rotateSend(:), rotateRecv(:), tileSend(:), tileRecv(:)

    if( domain%position .NE. CENTER )  call mpp_error(FATAL,  "mpp_domains_define.inc: " //&
       "routine define_contact_point can only be used to calculate overlapping for cell center.")

    ntiles      = domain%ntiles

    allocate(eCont(ntiles), sCont(ntiles), wCont(ntiles), nCont(ntiles) )
    eCont(:)%ncontact = 0; 

    do n = 1, ntiles
       eCont(n)%ncontact = 0; sCont(n)%ncontact = 0; wCont(n)%ncontact = 0; nCont(n)%ncontact = 0;
       allocate(eCont(n)%tile(num_contact), wCont(n)%tile(num_contact) )
       allocate(nCont(n)%tile(num_contact), sCont(n)%tile(num_contact) )
       allocate(eCont(n)%align1(num_contact), eCont(n)%align2(num_contact) )
       allocate(wCont(n)%align1(num_contact), wCont(n)%align2(num_contact) )
       allocate(sCont(n)%align1(num_contact), sCont(n)%align2(num_contact) )
       allocate(nCont(n)%align1(num_contact), nCont(n)%align2(num_contact) )
       allocate(eCont(n)%refine1(num_contact), eCont(n)%refine2(num_contact) )
       allocate(wCont(n)%refine1(num_contact), wCont(n)%refine2(num_contact) )
       allocate(sCont(n)%refine1(num_contact), sCont(n)%refine2(num_contact) )
       allocate(nCont(n)%refine1(num_contact), nCont(n)%refine2(num_contact) )
       allocate(eCont(n)%is1(num_contact), eCont(n)%ie1(num_contact), eCont(n)%js1(num_contact), eCont(n)%je1(num_contact))
       allocate(eCont(n)%is2(num_contact), eCont(n)%ie2(num_contact), eCont(n)%js2(num_contact), eCont(n)%je2(num_contact))
       allocate(wCont(n)%is1(num_contact), wCont(n)%ie1(num_contact), wCont(n)%js1(num_contact), wCont(n)%je1(num_contact))
       allocate(wCont(n)%is2(num_contact), wCont(n)%ie2(num_contact), wCont(n)%js2(num_contact), wCont(n)%je2(num_contact))
       allocate(sCont(n)%is1(num_contact), sCont(n)%ie1(num_contact), sCont(n)%js1(num_contact), sCont(n)%je1(num_contact))
       allocate(sCont(n)%is2(num_contact), sCont(n)%ie2(num_contact), sCont(n)%js2(num_contact), sCont(n)%je2(num_contact))
       allocate(nCont(n)%is1(num_contact), nCont(n)%ie1(num_contact), nCont(n)%js1(num_contact), nCont(n)%je1(num_contact))
       allocate(nCont(n)%is2(num_contact), nCont(n)%ie2(num_contact), nCont(n)%js2(num_contact), nCont(n)%je2(num_contact))
    end do

    !--- set up the east, south, west and north contact for each tile.
    do n = 1, num_contact
       t1 = tile1(n)
       t2 = tile2(n)
       select case(align1(n))
       case (EAST)
          call fill_contact( eCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
                             jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (WEST)
          call fill_contact( wCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
                             jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (SOUTH)
          call fill_contact( sCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
                             jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       case (NORTH)
          call fill_contact( nCont(t1), t2, istart1(n), iend1(n), jstart1(n), jend1(n), istart2(n), iend2(n), &
                             jstart2(n), jend2(n), align1(n), align2(n), refine1(n), refine2(n))
       end select
       select case(align2(n))
       case (EAST)
          call fill_contact( eCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
                             jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (WEST)
          call fill_contact( wCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
                             jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (SOUTH)
          call fill_contact( sCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
                             jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       case (NORTH)
          call fill_contact( nCont(t2), t1, istart2(n), iend2(n), jstart2(n), jend2(n), istart1(n), iend1(n), &
                             jstart1(n), jend1(n), align2(n), align1(n), refine2(n), refine1(n))
       end select
    end do

    !--- the tile number of current pe, halo size 
    whalo = domain%whalo
    ehalo = domain%ehalo
    shalo = domain%shalo
    nhalo = domain%nhalo

    !--- find if there is an extra point in x and y direction depending on position
    nlist = size(domain%list(:))

    max_contact = 4*num_contact ! should be enough

    !--- is2Send and is1Send will figure out the overlapping for sending from current pe.
    !--- is1Recv and iscREcv will figure out the overlapping for recving onto current pe.
    allocate( is1Send(max_contact), ie1Send(max_contact),  js1Send(max_contact),    je1Send(max_contact) )
    allocate( is2Send(max_contact), ie2Send(max_contact),  js2Send(max_contact),    je2Send(max_contact) )
    allocate( is2Recv(max_contact), ie2Recv(max_contact),  js2Recv(max_contact),    je2Recv(max_contact) )
    allocate( is1Recv(max_contact), ie1Recv(max_contact),  js1Recv(max_contact),    je1Recv(max_contact) )
    allocate( align1Recv(max_contact), align2Recv(max_contact), rotateSend(max_contact), rotateRecv(max_contact) )
    allocate( align1Send(max_contact), align2Send(max_contact), tileSend(max_contact), tileRecv(max_contact) )
    allocate( refineSend(max_contact), refineRecv(max_contact) )
    ntileMe = size(domain%x(:))
    refineSend = 1; refineRecv = 1

    !--------------------------------------------------------------------------------------------------
    !    loop over each tile on current domain to set up the overlapping for each tile
    !--------------------------------------------------------------------------------------------------

    do tMe = 1, ntileMe
       tileMe = domain%tile_id(tMe)
       rotateSend = ZERO; rotateRecv = ZERO

       !--- loop over all the contact region to figure out the index for overlapping region.
       count = 0
       do n = 1, eCont(tileMe)%ncontact  ! east contact
          count = count+1
          tileRecv(count)   = eCont(tileMe)%tile(n);    tileSend(count)   = eCont(tileMe)%tile(n)
          align1Recv(count) = eCont(tileMe)%align1(n);  align2Recv(count) = eCont(tileMe)%align2(n)
          align1Send(count) = eCont(tileMe)%align1(n);  align2Send(count) = eCont(tileMe)%align2(n)
          refineSend(count) = eCont(tileMe)%refine2(n); refineRecv(count) = eCont(tileMe)%refine1(n)
          is1Recv(count)    = eCont(tileMe)%is1(n) + 1; ie1Recv(count)    = is1Recv(count) + ehalo - 1  
          js1Recv(count)    = eCont(tileMe)%js1(n);     je1Recv(count)    = eCont(tileMe)%je1(n)
          select case(eCont(tileMe)%align2(n))
          case ( WEST )  ! w <-> e 
             is2Recv(count) = eCont(tileMe)%is2(n);     ie2Recv(count) = is2Recv(count) + ehalo - 1
             js2Recv(count) = eCont(tileMe)%js2(n);     je2Recv(count) = eCont(tileMe)%je2(n)
             ie1Send(count) = eCont(tileMe)%is1(n);     is1Send(count) = ie1Send(count) - whalo + 1  
             js1Send(count) = eCont(tileMe)%js1(n);     je1Send(count) = eCont(tileMe)%je1(n)
             ie2Send(count) = eCont(tileMe)%is2(n) - 1; is2Send(count) = ie2Send(count) - whalo + 1  
             js2Send(count) = eCont(tileMe)%js2(n);     je2Send(count) = eCont(tileMe)%je2(n)
          case ( SOUTH ) ! s <-> e
             rotateRecv(count) = NINETY;                rotateSend(count) = MINUS_NINETY
             js2Recv(count) = eCont(tileMe)%js2(n);     je2Recv(count) = js2Recv(count) + ehalo -1
             is2Recv(count) = eCont(tileMe)%is2(n);     ie2Recv(count) = eCont(tileMe)%ie2(n)
             ie1Send(count) = eCont(tileMe)%is1(n);     is1Send(count) = ie1Send(count) - shalo + 1  
             js1Send(count) = eCont(tileMe)%js1(n);     je1Send(count) = eCont(tileMe)%je1(n)
             is2Send(count) = eCont(tileMe)%is2(n);     ie2Send(count) = eCont(tileMe)%ie2(n)
             je2Send(count) = eCont(tileMe)%js2(n) - 1; js2Send(count) = je2Send(count) - shalo + 1  
          end select
       end do

       do n = 1, sCont(tileMe)%ncontact  ! south contact 
          count = count+1
          tileRecv(count)   = sCont(tileMe)%tile(n);    tileSend(count)   = sCont(tileMe)%tile(n)
          align1Recv(count) = sCont(tileMe)%align1(n);  align2Recv(count) = sCont(tileMe)%align2(n);
          align1Send(count) = sCont(tileMe)%align1(n);  align2Send(count) = sCont(tileMe)%align2(n);
          refineSend(count) = sCont(tileMe)%refine2(n); refineRecv(count) = sCont(tileMe)%refine1(n)
          is1Recv(count)    = sCont(tileMe)%is1(n);     ie1Recv(count)    = sCont(tileMe)%ie1(n)
          je1Recv(count)    = sCont(tileMe)%js1(n) - 1; js1Recv(count)    = je1Recv(count) - shalo + 1
          select case(sCont(tileMe)%align2(n))
          case ( NORTH ) ! n <-> s
             is2Recv(count) = sCont(tileMe)%is2(n);   ie2Recv(count) = sCont(tileMe)%ie2(n)
             je2Recv(count) = sCont(tileMe)%je2(n);   js2Recv(count) = je2Recv(count) - shalo + 1
             is1Send(count) = sCont(tileMe)%is1(n);   ie1Send(count) = sCont(tileMe)%ie1(n)
             js1Send(count) = sCont(tileMe)%js1(n);   je1Send(count) = js1Send(count) + nhalo -1
             is2Send(count) = sCont(tileMe)%is2(n);   ie2Send(count) = sCont(tileMe)%ie2(n)
             js2Send(count) = sCont(tileMe)%je2(n)+1; je2Send(count) = js2Send(count) + nhalo - 1
          case ( EAST )  ! e <-> s
             rotateRecv(count) = MINUS_NINETY;        rotateSend(count) = NINETY
             ie2Recv(count) = sCont(tileMe)%ie2(n);   is2Recv(count) = ie2Recv(count) - shalo + 1
             js2Recv(count) = sCont(tileMe)%js2(n);   je2Recv(count) = sCont(tileMe)%je2(n)
             is1Send(count) = sCont(tileMe)%is1(n);   ie1Send(count) = sCont(tileMe)%ie1(n)
             js1Send(count) = sCont(tileMe)%js1(n);   je1Send(count) = js1Send(count) + ehalo - 1
             is2Send(count) = sCont(tileMe)%ie2(n)+1; ie2Send(count) = is2Send(count) + ehalo - 1
             js2Send(count) = sCont(tileMe)%js2(n);   je2Send(count) = sCont(tileMe)%je2(n)
          end select
       end do

       do n = 1, wCont(tileMe)%ncontact  ! west contact
          count = count+1
          tileRecv(count)   = wCont(tileMe)%tile(n);    tileSend(count)   = wCont(tileMe)%tile(n)
          align1Recv(count) = wCont(tileMe)%align1(n);  align2Recv(count) = wCont(tileMe)%align2(n);
          align1Send(count) = wCont(tileMe)%align1(n);  align2Send(count) = wCont(tileMe)%align2(n);
          refineSend(count) = wCont(tileMe)%refine2(n); refineRecv(count) = wCont(tileMe)%refine1(n)
          ie1Recv(count)    = wCont(tileMe)%is1(n) - 1; is1Recv(count)    = ie1Recv(count) - whalo + 1
          js1Recv(count)    = wCont(tileMe)%js1(n);     je1Recv(count)    = wCont(tileMe)%je1(n)
          select case(wCont(tileMe)%align2(n))
          case ( EAST )  ! e <-> w
             ie2Recv(count) = wCont(tileMe)%ie2(n);   is2Recv(count) = ie2Recv(count) - whalo + 1
             js2Recv(count) = wCont(tileMe)%js2(n);   je2Recv(count) = wCont(tileMe)%je2(n) 
             is1Send(count) = wCont(tileMe)%is1(n);   ie1Send(count) = is1Send(count) + ehalo - 1 
             js1Send(count) = wCont(tileMe)%js1(n);   je1Send(count) = wCont(tileMe)%je1(n) 
             is2Send(count) = wCont(tileMe)%ie2(n)+1; ie2Send(count) = is2Send(count) + ehalo - 1 
             js2Send(count) = wCont(tileMe)%js2(n);   je2Send(count) = wCont(tileMe)%je2(n) 
          case ( NORTH ) ! n <-> w
             rotateRecv(count) = NINETY;              rotateSend(count) = MINUS_NINETY
             je2Recv(count) = wCont(tileMe)%je2(n);   js2Recv(count) = je2Recv(count) - whalo + 1
             is2Recv(count) = wCont(tileMe)%is2(n);   ie2Recv(count) = wCont(tileMe)%ie2(n)
             is1Send(count) = wCont(tileMe)%is1(n);   ie1Send(count) = is1Send(count) + nhalo - 1
             js1Send(count) = wCont(tileMe)%js1(n);   je1Send(count) = wCont(tileMe)%je1(n)
             js2Send(count) = wCont(tileMe)%je2(n)+1; je2Send(count) = js2Send(count) + nhalo - 1
             is2Send(count) = wCont(tileMe)%is2(n);   ie2Send(count) = wCont(tileMe)%ie2(n)
          end select
       end do

       do n = 1, nCont(tileMe)%ncontact  ! north contact 
          count = count+1
          tileRecv(count)   = nCont(tileMe)%tile(n);   tileSend(count)   = nCont(tileMe)%tile(n)
          align1Recv(count) = nCont(tileMe)%align1(n); align2Recv(count) = nCont(tileMe)%align2(n);
          align1Send(count) = nCont(tileMe)%align1(n); align2Send(count) = nCont(tileMe)%align2(n);
          refineSend(count) = nCont(tileMe)%refine2(n); refineRecv(count) = nCont(tileMe)%refine1(n)
          is1Recv(count)    = nCont(tileMe)%is1(n);    ie1Recv(count)    = nCont(tileMe)%ie1(n)
          js1Recv(count)    = nCont(tileMe)%je1(n)+1;  je1Recv(count)    = js1Recv(count) + nhalo - 1
          select case(nCont(tileMe)%align2(n))
          case ( SOUTH ) ! s <-> n
             is2Recv(count) = nCont(tileMe)%is2(n);   ie2Recv(count) = nCont(tileMe)%ie2(n)
             js2Recv(count) = nCont(tileMe)%js2(n);   je2Recv(count) = js2Recv(count) + nhalo - 1
             is1Send(count) = nCont(tileMe)%is1(n);   ie1Send(count) = nCont(tileMe)%ie1(n)
             je1Send(count) = nCont(tileMe)%je1(n);   js1Send(count) = je1Send(count) - shalo + 1
             is2Send(count) = nCont(tileMe)%is2(n);   ie2Send(count) = nCont(tileMe)%ie2(n)
             je2Send(count) = nCont(tileMe)%js2(n)-1; js2Send(count) = je2Send(count) - shalo + 1
          case ( WEST )  ! w <-> n
             rotateRecv(count) = MINUS_NINETY;        rotateSend(count) = NINETY
             is2Recv(count) = nCont(tileMe)%ie2(n);   ie2Recv(count) = is2Recv(count) + nhalo - 1
             js2Recv(count) = nCont(tileMe)%js2(n);   je2Recv(count) = nCont(tileMe)%je2(n)
             is1Send(count) = nCont(tileMe)%is1(n);   ie1Send(count) = nCont(tileMe)%ie1(n)
             je1Send(count) = nCont(tileMe)%je1(n);   js1Send(count) = je1Send(count) - whalo + 1
             ie2Send(count) = nCont(tileMe)%is2(n)-1; is2Send(count) = ie2Send(count) - whalo + 1
             js2Send(count) = nCont(tileMe)%js2(n);   je2Send(count) = nCont(tileMe)%je2(n)
          end select
       end do

       numS = count
       numR = count
       !--- figure out the index for corner overlapping,
       !--- fill_corner_contact will be updated to deal with the situation that there are multiple tiles on 
       !--- each side of six sides of cubic grid.
       if(.NOT. domain%rotated_ninety)  then
           call fill_corner_contact(eCont, sCont, wCont, nCont, isgList, iegList, jsgList, jegList, numR, numS, &
                                tileRecv, tileSend, is1Recv, ie1Recv, js1Recv, je1Recv, is2Recv, ie2Recv,   &
                                js2Recv, je2Recv, is1Send, ie1Send, js1Send, je1Send, is2Send, ie2Send,     &
                                js2Send, je2Send, align1Recv, align2Recv, align1Send, align2Send,           &
                                whalo, ehalo, shalo, nhalo, tileMe )
       end if

       isc = domain%x(tMe)%compute%begin; iec = domain%x(tMe)%compute%end
       jsc = domain%y(tMe)%compute%begin; jec = domain%y(tMe)%compute%end

       !--- compute the overlapping for send.
       do n = 1, numS
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             ntileNbr = size(domain%list(m)%x(:))
             do tNbr = 1, ntileNbr
                if( domain%list(m)%tile_id(tNbr) .NE. tileSend(n) ) cycle
                isc1 = max(isc, is1Send(n)); iec1 = min(iec, ie1Send(n))
                jsc1 = max(jsc, js1Send(n)); jec1 = min(jec, je1Send(n))  
                if( isc1 > iec1 .OR. jsc1 > jec1 ) cycle
                !--- loop over 8 direction to get the overlapping starting from east with clockwise.
                do dir = 1, 8
                   !--- get the to_pe's data domain.
                   select case ( dir )
                   case ( 1 )  ! eastern halo
                      if( align2Send(n) .NE. EAST ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%data%end
                      jsd = domain%list(m)%y(tNbr)%compute%begin; jed = domain%list(m)%y(tNbr)%compute%end
                   case ( 2 )  ! southeast halo
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%data%end
                      jsd = domain%list(m)%y(tNbr)%data%begin;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 3 )  ! southern halo
                      if( align2Send(n) .NE. SOUTH ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%begin; ied = domain%list(m)%x(tNbr)%compute%end
                      jsd = domain%list(m)%y(tNbr)%data%begin;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 4 )  ! southwest halo
                      isd = domain%list(m)%x(tNbr)%data%begin;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%data%begin;    jed = domain%list(m)%y(tNbr)%compute%begin-1
                   case ( 5 )  ! western halo
                      if( align2Send(n) .NE. WEST ) cycle
                      isd = domain%list(m)%x(tNbr)%data%begin;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%compute%begin; jed = domain%list(m)%y(tNbr)%compute%end
                   case ( 6 )  ! northwest halo
                      isd = domain%list(m)%x(tNbr)%data%begin;    ied = domain%list(m)%x(tNbr)%compute%begin-1
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%data%end
                   case ( 7 )  ! northern halo
                      if( align2Send(n) .NE. NORTH ) cycle
                      isd = domain%list(m)%x(tNbr)%compute%begin; ied = domain%list(m)%x(tNbr)%compute%end
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%data%end
                   case ( 8 )  ! northeast halo
                      isd = domain%list(m)%x(tNbr)%compute%end+1; ied = domain%list(m)%x(tNbr)%data%end
                      jsd = domain%list(m)%y(tNbr)%compute%end+1; jed = domain%list(m)%y(tNbr)%data%end
                   end select
                   isd = max(isd, is2Send(n)); ied = min(ied, ie2Send(n))
                   jsd = max(jsd, js2Send(n)); jed = min(jed, je2Send(n))  
                   if( isd > ied .OR. jsd > jed ) cycle
                   ioff = 0; joff = 0
                   nxd = ied - isd + 1
                   nyd = jed - jsd + 1
                   select case ( align2Send(n) )
                   case ( WEST, EAST )
                      joff = jsd - js2Send(n)
                      joff = get_refine_count( joff, refineSend(n), .true., dir )
                      nyd  = get_refine_count( nyd, refineSend(n), .false., dir  )
                   case ( SOUTH, NORTH )
                      ioff = isd - is2Send(n) 
                      ioff = get_refine_count( ioff, refineSend(n), .true., dir  )
                      nxd  = get_refine_count( nxd, refineSend(n), .false., dir  )
                   end select

                   !--- get the index in current pe.
                   select case ( rotateSend(n) )
                   case ( ZERO )
                      isc2 = is1Send(n) + ioff; iec2 = isc2 + nxd - 1
                      jsc2 = js1Send(n) + joff; jec2 = jsc2 + nyd - 1
                   case ( NINETY )                     ! N -> W or S -> E
                      iec2 = ie1Send(n) - joff; isc2 = iec2 - nyd + 1 
                      jsc2 = js1Send(n) + ioff; jec2 = jsc2 + nxd - 1
                   case ( MINUS_NINETY )               ! W -> N or E -> S
                      isc2 = is1Send(n) + joff; iec2 = isc2 + nyd - 1
                      jec2 = je1Send(n) - ioff; jsc2 = jec2 - nxd + 1 
                   end select
                   is = max(isc1,isc2); ie = min(iec1,iec2)
                   js = max(jsc1,jsc2); je = min(jec1,jec2)
                   if(ie.GE.is .AND. je.GE.js )then
                      if( domain%list(m)%send(tNbr,tMe,dir)%overlap(1) ) then
                          call mpp_error(FATAL, &
                           "mpp_domains_define.inc:  the overlapping for send should not already been defined")
                      end if
                      if(refineSend(n) .NE. 1) domain%list(m)%send(tNbr,tMe,dir)%is_refined(1) = .true.   
                      domain%list(m)%send(tNbr,tMe,dir)%overlap(1) = .TRUE.
                      domain%list(m)%send(tNbr,tMe,dir)%rotation = rotateSend(n)
                      domain%list(m)%send(tNbr,tMe,dir)%is(1) = is
                      domain%list(m)%send(tNbr,tMe,dir)%ie(1) = ie
                      domain%list(m)%send(tNbr,tMe,dir)%js(1) = js
                      domain%list(m)%send(tNbr,tMe,dir)%je(1) = je
                      if(refineSend(n) .NE. 1)  domain%list(m)%send(tNbr,tMe,dir)%is_refined(1) = .true.           
                   endif
                end do ! end do dir = 1, 8
             end do ! end do tNbr = 1, ntileNbr
             if( ANY(domain%list(m)%send(:,:,:)%overlap(1)) ) domain%list(m)%overlap = .true.
          end do ! end do list = 0, nlist-1
       end do ! end do n = 1, numS

       !--- compute the overlapping for recv.
       do n = 1, numR
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             ntileNbr = size(domain%list(m)%x(:))
             do tNbr = 1, ntileNbr
                if( domain%list(m)%tile_id(tNbr) .NE. tileRecv(n) ) cycle
                isc = domain%list(m)%x(tNbr)%compute%begin; iec = domain%list(m)%x(tNbr)%compute%end
                jsc = domain%list(m)%y(tNbr)%compute%begin; jec = domain%list(m)%y(tNbr)%compute%end
                isc = max(isc, is2Recv(n)); iec = min(iec, ie2Recv(n))
                jsc = max(jsc, js2Recv(n)); jec = min(jec, je2Recv(n))      
                if( isc > iec .OR. jsc > jec ) cycle
                !--- find the offset for this overlapping.
                ioff = 0; joff = 0
                nxc = iec - isc + 1; nyc = jec - jsc + 1
                select case ( align2Recv(n) )
                case ( WEST, EAST )
                   joff = jsc - js2Recv(n)
                   joff = get_refine_count( joff, 1/refineRecv(n), .true., dir  )
                   nyc  = get_refine_count( nyc, 1/refineRecv(n), .false., dir  )
                case ( NORTH, SOUTH )
                   ioff = isc - is2Recv(n) 
                   ioff = get_refine_count( ioff, 1/refineRecv(n), .true., dir  )
                   nxc  = get_refine_count( nxc, 1/refineRecv(n), .false., dir  )
                end select

                !--- get the index in current pe.
                select case ( rotateRecv(n) )
                case ( ZERO )
                   isd1 = is1Recv(n) + ioff; ied1 = isd1 + nxc - 1
                   jsd1 = js1Recv(n) + joff; jed1 = jsd1 + nyc - 1
                case ( NINETY )                      ! N -> W or S -> E
                   isd1 = is1Recv(n);        ied1 = isd1 + nyc - 1
                   jed1 = je1Recv(n) - ioff; jsd1 = jed1 - nxc + 1
                case ( MINUS_NINETY )                ! W -> N or E -> S
                   ied1 = ie1Recv(n) - joff; isd1 = ied1 - nyc + 1
                   jsd1 = js1Recv(n);        jed1 = jsd1 + nxc - 1
                end select

                !--- loop over 8 direction to get the overlapping starting from east with clockwise.
                do dir = 1, 8
                   select case ( dir )
                   case ( 1 )  ! eastern halo
                      if( align1Recv(n) .NE. EAST ) cycle
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%compute%begin; jed2 = domain%y(tMe)%compute%end
                   case ( 2 )  ! southeast halo
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 3 )  ! southern halo
                      if( align1Recv(n) .NE. SOUTH ) cycle
                      isd2 = domain%x(tMe)%compute%begin; ied2 = domain%x(tMe)%compute%end
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 4 )  ! southwest halo
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%data%begin;    jed2 = domain%y(tMe)%compute%begin-1
                   case ( 5 )  ! western halo
                      if( align1Recv(n) .NE. WEST ) cycle
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%compute%begin; jed2 = domain%y(tMe)%compute%end
                   case ( 6 )  ! northwest halo
                      isd2 = domain%x(tMe)%data%begin;    ied2 = domain%x(tMe)%compute%begin-1
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   case ( 7 )  ! northern halo
                      if( align1Recv(n) .NE. NORTH ) cycle
                      isd2 = domain%x(tMe)%compute%begin; ied2 = domain%x(tMe)%compute%end
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   case ( 8 )  ! northeast halo
                      isd2 = domain%x(tMe)%compute%end+1; ied2 = domain%x(tMe)%data%end
                      jsd2 = domain%y(tMe)%compute%end+1; jed2 = domain%y(tMe)%data%end
                   end select
                   is = max(isd1,isd2); ie = min(ied1,ied2)
                   js = max(jsd1,jsd2); je = min(jed1,jed2)
                   if(ie.GE.is .AND. je.GE.js )then
                      if( domain%list(m)%recv(tMe,tNbr,dir)%overlap(1) ) call mpp_error(FATAL, &
                           "mpp_domains_define.inc:  the overlapping for recv should not already been defined")
                      domain%list(m)%recv(tMe,tNbr,dir)%overlap(1) = .TRUE.
                      if(refineRecv(n) .NE. 1) domain%list(m)%recv(tMe,tNbr,dir)%is_refined(1) = .true.   
                      domain%list(m)%recv(tMe,tNbr,dir)%rotation = rotateRecv(n)
                      domain%list(m)%recv(tMe,tNbr,dir)%is(1) = is
                      domain%list(m)%recv(tMe,tNbr,dir)%ie(1) = ie
                      domain%list(m)%recv(tMe,tNbr,dir)%js(1) = js
                      domain%list(m)%recv(tMe,tNbr,dir)%je(1) = je
                      domain%list(m)%recv(tMe,tNbr,dir)%isMe  = is
                      domain%list(m)%recv(tMe,tNbr,dir)%ieMe  = ie
                      domain%list(m)%recv(tMe,tNbr,dir)%jsMe  = js
                      domain%list(m)%recv(tMe,tNbr,dir)%jeMe  = je
                      !--- get the index on the pe sending data when there is refinement
                      if( domain%list(m)%recv(tMe,tNbr,dir)%is_refined(1) ) then
                         select case ( align1Recv(n) )
                         case ( WEST, EAST )
                            isoff = 0;               ieoff = ie - is + 1
                            jsoff = js - js1Recv(n); jeoff = je - js1Recv(n) + 1
                            jsoff = get_refine_count( jsoff, refineRecv(n), .true., dir  )
                            jeoff = get_refine_count( jeoff, refineRecv(n), .false., dir  ) - 1
                         case ( NORTH, SOUTH )
                            isoff = is - is1Recv(n); ieoff = ie - is1Recv(n) + 1
                            jsoff = 0;               jeoff = je - js
                            isoff = get_refine_count( isoff, refineRecv(n), .true., dir  )
                            ieoff = get_refine_count( ieoff, refineRecv(n), .false., dir  ) - 1
                         end select

                         !--- get the index in neighbor pe.
                         select case ( rotateRecv(n) )
                         case ( ZERO )
                            isd2 = is2Recv(n) + isoff; ied2 = is2Recv(n) + ieoff
                            jsd2 = js2Recv(n) + jsoff; jed2 = js2Recv(n) + jeoff
                         case ( NINETY )                      ! N -> W or S -> E
                            isd2 = ie2Recv(n) - jeoff; ied2 = ie2Recv(n) - jsoff
                            jsd2 = js2Recv(n);         jed2 = js2Recv(n) + ieoff - isoff
                         case ( MINUS_NINETY )                ! W -> N or E -> S
                            isd2 = is2Recv(n);         ied2 = is2Recv(n) + jeoff - jsoff
                            jsd2 = je2Recv(n) - ieoff; jed2 = je2Recv(n) - isoff
                         end select
                         isd2= max(isd2,isc); ied2 = min(ied2, iec); jsd2= max(jsd2,jsc); jed2 = min(jed2, jec)
                         if( isd2 > ied2 .OR. jsd2 > jed2 ) call mpp_error(FATAL, &
                              "mpp_domains_define.inc:  the neighbor index is out of range")
                         domain%list(m)%recv(tMe,tNbr,dir)%is(1) = isd2
                         domain%list(m)%recv(tMe,tNbr,dir)%ie(1) = ied2
                         domain%list(m)%recv(tMe,tNbr,dir)%js(1) = jsd2
                         domain%list(m)%recv(tMe,tNbr,dir)%je(1) = jed2
                      end if
                   endif
                end do ! end do dir = 1, 8
             end do ! end do tNbr = 1, ntileNbr
             if( ANY(domain%list(m)%recv(:,:,:)%overlap(1)) ) domain%list(m)%overlap = .true.
          end do ! end do list = 0, nlist-1
       end do ! end do n = 1, numR
    end do   ! end do tMe = 1, ntileMe

    !--- release memory
    deallocate( is1Send, ie1Send, js1Send, je1Send, is2Send, ie2Send, js2Send, je2Send )
    deallocate( is2Recv, ie2Recv, js2Recv, je2Recv, is1Recv, ie1Recv, js1Recv, je1Recv )
    deallocate( align1Send, align2Send, align1Recv, align2Recv, refineSend, refineRecv)
    deallocate( rotateSend, rotateRecv, tileSend, tileRecv )
    do n = 1, ntiles
       deallocate(eCont(n)%tile, wCont(n)%tile, sCont(n)%tile, nCont(n)%tile )
       deallocate(eCont(n)%align1, wCont(n)%align1, sCont(n)%align1, nCont(n)%align1)
       deallocate(eCont(n)%align2, wCont(n)%align2, sCont(n)%align2, nCont(n)%align2)
       deallocate(eCont(n)%refine1, wCont(n)%refine1, sCont(n)%refine1, nCont(n)%refine1)
       deallocate(eCont(n)%refine2, wCont(n)%refine2, sCont(n)%refine2, nCont(n)%refine2)
       deallocate(eCont(n)%is1, eCont(n)%ie1, eCont(n)%js1, eCont(n)%je1 )
       deallocate(eCont(n)%is2, eCont(n)%ie2, eCont(n)%js2, eCont(n)%je2 )
       deallocate(wCont(n)%is1, wCont(n)%ie1, wCont(n)%js1, wCont(n)%je1 )
       deallocate(wCont(n)%is2, wCont(n)%ie2, wCont(n)%js2, wCont(n)%je2 )
       deallocate(sCont(n)%is1, sCont(n)%ie1, sCont(n)%js1, sCont(n)%je1 )
       deallocate(sCont(n)%is2, sCont(n)%ie2, sCont(n)%js2, sCont(n)%je2 )
       deallocate(nCont(n)%is1, nCont(n)%ie1, nCont(n)%js1, nCont(n)%je1 )
       deallocate(nCont(n)%is2, nCont(n)%ie2, nCont(n)%js2, nCont(n)%je2 )
    end do
    deallocate(eCont, wCont, sCont, nCont)

    !--- define the refined overlapping and group the overlapping 
    !--- between same tiles( in different direction or different pe ). 
    !--- set up the index position in the return buffer.
    call define_refine_overlap( domain )  

    domain%initialized = .true.

  end subroutine define_contact_point

  !--- currently assume no grid is divided over two pes.
  !--- we may remove this restriction in the future. 
  function get_refine_count(num , refine, use_floor, dir )
     integer, intent(in) :: num
     real,    intent(in) :: refine
     logical, intent(in) :: use_floor
     integer, intent(in) :: dir
     integer             :: get_refine_count
     real, parameter     :: EPSLN = 1.e-10
     real                :: rval

     rval = num*refine

     if(mod(dir,2) == 1) then ! WEST, EAST, SOUTH, NORTH
        get_refine_count = nint(rval)
        if(abs(rval-get_refine_count) > EPSLN) call mpp_error(FATAL, &
             "get_refine_count: The two adjacent tiles with refinement does not align well in E/W/S/N direction")
     else
        get_refine_count = floor(rval)
        if(use_floor) then
           if( rval - get_refine_count > 1 - EPSLN ) get_refine_count = get_refine_count + 1 ! Deal with truncation error
        else
           if(rval - get_refine_count  > EPSLN) get_refine_count = get_refine_count + 1 ! Deal with truncation error
        end if
     end if
     

     return

  end function get_refine_count


  !####################################################################################
  !--- The following will define the refined overlapping and group the overlapping 
  !--- between same tiles( in different direction or different pe ). 
  !--- set up the index position in the return buffer.
  subroutine define_refine_overlap(domain)
    type(domain2d), intent(inout) :: domain

    logical                        :: new_overlap
    integer,           parameter   :: MAXOVERLAP = 100
    integer, dimension(MAXOVERLAP) :: isMe, ieMe, jsMe, jeMe, isNbr, ieNbr, jsNbr, jeNbr
    integer, dimension(MAXOVERLAP) :: rotation, tileList, dirList 
    integer                        :: count, tMe, tNbr, ntileMe, ntileNbr, nlist
    integer                        :: list, m, n, dir, start, end, dirdiff
    type(overlapSpec), pointer     :: OverPtr => NULL()

    ntileMe = size(domain%x(:))
    nlist = size(domain%list(:))    

    do tMe = 1, ntileMe
       count = 0
       do dir = 1, 8  ! --- the overlapping will always take EAST, SOUTH, WEST and NORTH order.
          do list = nlist - 1, 0, -1
             m = mod( domain%pos+nlist-list, nlist )
             if( .NOT.domain%list(m)%overlap)cycle
             ntileNbr = size(domain%list(m)%x(:))
             do tNbr = ntileNbr, 1, -1
                overPtr => domain%list(m)%recv(tMe,tNbr,dir)
                if( overPtr%overlap(1) )then 
                   if(OverPtr%is_refined(1)) then
                      new_overlap = .true.
                      do n = 1, count
                         if( tileList(n) == domain%list(m)%tile_id(tNbr) ) then ! same tile
                            dirdiff = dir - dirList(n)   ! dir always no less than dirList(n). 
                            if( dir == 8 .AND. dirList(n) == 1) dirdiff = 1
                            if( dirdiff .LE. 1 ) then    ! same direction
                               if( mod(dirList(n),2) == 0 ) then
                                  if(mod(dir,2) == 0 ) call mpp_error(FATAL, &
                                       "mpp_get_mosaic_refine_overlap: corner overlap should be unique.")
                                  dirList(n) = dir
                               end if
                               if( dirList(n) == 1 .OR. dirList(n) == 5 ) then ! EAST, or WEST
                                  if( isMe(n) .NE. overPtr%isMe .OR. ieMe(n) .NE. overPtr%ieMe ) call mpp_error(FATAL, &
                                       "mpp_get_mosaic_refine_overlap: mismatch on my index in the overlapping 1" )
                                  jeMe(n) = max(jeMe(n),  overPtr%jeMe )
                                  jsMe(n) = min(jsMe(n),  overPtr%jsMe )
                               else if(dirList(n) == 3 .OR. dirList(n) == 7 ) then ! SOUTH or NORTH
                                  if( jsMe(n) .NE. overPtr%jsMe .OR. jeMe(n) .NE. overPtr%jeMe ) call mpp_error(FATAL, &
                                       "mpp_get_mosaic_refine_overlap: mismatch on my index in the overlapping 2" )
                                  ieMe(n) = max(ieMe(n),  overPtr%ieMe )
                                  isMe(n) = min(isMe(n),  overPtr%isMe )
                               end if
                               isNbr(n) = min(isNbr(n), overPtr%is(1))
                               ieNbr(n) = max(ieNbr(n), overPtr%ie(1))
                               jsNbr(n) = min(jsNbr(n), overPtr%js(1))
                               jeNbr(n) = max(jeNbr(n), overPtr%je(1))
                               OverPtr%index = n
                               new_overlap = .false.
                               exit
                            end if
                         end if
                      end do

                      if(new_overlap) then
                         count = count + 1
                         if(count > MAXOVERLAP) call mpp_error(FATAL,  &
                              "mpp_get_mosaic_refine_overlap: number of overlapping is larger than MAXOVERLAP")
                         tileList(count) = domain%list(m)%tile_id(tNbr)
                         dirList(count)  = dir
                         rotation(count) = overPtr%rotation
                         OverPtr%index = count
                         isMe(count) = overPtr%isMe;     ieMe(count) = overPtr%ieMe
                         jsMe(count) = overPtr%jsMe;     jeMe(count) = overPtr%jeMe
                         isNbr(count) = overPtr%is(1);     ieNbr(count) = overPtr%ie(1)
                         jsNbr(count) = overPtr%js(1);     jeNbr(count) = overPtr%je(1)
                      end if
                   end if
                end if
             end do
          end do
       end do
       if(count > 0) then
          allocate(domain%rSpec(tMe)%isMe(count),  domain%rSpec(tMe)%ieMe(count) )
          allocate(domain%rSpec(tMe)%jsMe(count),  domain%rSpec(tMe)%jeMe(count) )
          allocate(domain%rSpec(tMe)%isNbr(count), domain%rSpec(tMe)%ieNbr(count) )
          allocate(domain%rSpec(tMe)%jsNbr(count), domain%rSpec(tMe)%jeNbr(count) )
          allocate(domain%rSpec(tMe)%start(count), domain%rSpec(tMe)%end  (count) )
          allocate(domain%rSpec(tMe)%rotation(count) )
          domain%rSpec(tMe)%isMe = isMe(1:count); domain%rSpec(tMe)%isNbr = isNbr(1:count)
          domain%rSpec(tMe)%ieMe = ieMe(1:count); domain%rSpec(tMe)%ieNbr = ieNbr(1:count)
          domain%rSpec(tMe)%jsMe = jsMe(1:count); domain%rSpec(tMe)%jsNbr = jsNbr(1:count)
          domain%rSpec(tMe)%jeMe = jeMe(1:count); domain%rSpec(tMe)%jeNbr = jeNbr(1:count)
          domain%rSpec(tMe)%count = count
          domain%rSpec(tMe)%rotation = rotation(1:count)
          domain%rSpec(tMe)%total = sum( (domain%rSpec(tMe)%ieNbr-domain%rSpec(tMe)%isNbr+1)* &
                                         (domain%rSpec(tMe)%jeNbr-domain%rSpec(tMe)%jsNbr+1)  )
          end = 0
          do n = 1, count
             start = end + 1
             end = start + (ieNbr(n) - isNbr(n) + 1)*(jeNbr(n) - jsNbr(n) + 1) - 1
             domain%rSpec(tMe)%start(n) = start
             domain%rSpec(tMe)%end(n)   = end
          end do
       end if
    end do

 end subroutine define_refine_overlap

 !##############################################################################
 !--- always fill the contact according to index order.
  subroutine fill_contact(Contact, tile, is1, ie1, js1, je1, is2, ie2, js2, je2, align1, align2, refine1, refine2 )
     type(contact_type), intent(inout) :: Contact
     integer,            intent(in)    :: tile
     integer,            intent(in)    :: is1, ie1, js1, je1
     integer,            intent(in)    :: is2, ie2, js2, je2
     integer,            intent(in)    :: align1, align2
     real,               intent(in)    :: refine1, refine2 
     integer                           :: pos, n

     do pos = 1, Contact%ncontact
        select case(align1)
        case(WEST, EAST)  
           if( js1 < Contact%js1(pos) ) exit
        case(SOUTH, NORTH)  
           if( is1 < Contact%is1(pos) ) exit
        end select
     end do

     Contact%ncontact = Contact%ncontact + 1
     do n = Contact%ncontact, pos+1, -1  ! shift the data if needed.
        Contact%tile(n)   = Contact%tile(n-1)
        Contact%align1(n) = Contact%align1(n-1)
        Contact%align2(n) = Contact%align2(n-1)
        Contact%is1(n) = Contact%is1(n-1); Contact%ie1(n) = Contact%ie1(n-1)
        Contact%js1(n) = Contact%js1(n-1); Contact%je1(n) = Contact%je1(n-1)
        Contact%is2(n) = Contact%is2(n-1); Contact%ie2(n) = Contact%ie2(n-1)
        Contact%js2(n) = Contact%js2(n-1); Contact%je2(n) = Contact%je2(n-1)
     end do

     Contact%tile(pos)    = tile
     Contact%align1(pos)  = align1
     Contact%align2(pos)  = align2
     Contact%refine1(pos) = refine1
     Contact%refine2(pos) = refine2    
     Contact%is1(pos) = is1; Contact%ie1(pos) = ie1
     Contact%js1(pos) = js1; Contact%je1(pos) = je1
     Contact%is2(pos) = is2; Contact%ie2(pos) = ie2
     Contact%js2(pos) = js2; Contact%je2(pos) = je2

  end subroutine fill_contact

  !############################################################################
  ! this routine sets the overlapping between tiles for E,C,N-cell based on T-cell overlapping
  subroutine set_contact_point(domain_in, domain_out)
    type(domain2d),    intent(in) :: domain_in
    type(domain2d), intent(inout) :: domain_out
    integer                       :: ishift, jshift, nlist, list, m
    integer                       :: ntileMe, ntileNbr, tMe, tNbr, dir
    integer                       :: isoff1, ieoff1, isoff2, ieoff2, jsoff1, jeoff1, jsoff2, jeoff2
    type(overlapSpec),  pointer   :: ptrIn  => NULL()
    type(overlapSpec),  pointer   :: ptrOut => NULL()

    ishift = domain_out%x(1)%shift; jshift = domain_out%y(1)%shift
    nlist = size(domain_in%list(:))
    ntileMe = size(domain_in%x(:))
    !--- loop over the list of domains.
    do list = 0, nlist-1
       m = mod( domain_in%pos+list, nlist )
       if( .NOT. domain_in%list(m)%overlap ) cycle
       ntileNbr = size(domain_in%list(m)%x(:))
       !--- loop over 8 direction to get the overlapping starting from east with clockwise.
       do dir = 1, 8
          do tMe = 1, ntileMe
             do tNbr = 1, ntileNbr
                ptrIn  => domain_in %list(m)%send(tNbr,tMe,dir)
                ptrOut => domain_out%list(m)%send(tNbr,tMe,dir)          
                PtrOut%is_refined = PtrIn%is_refined
                ! only set overlapping between tiles for send ( ptrOut%overlap(1) is false )
                if(ptrIn%overlap(1) .AND. .not. ptrOut%overlap(1)) then  
                   ptrOut%overlap(1) = .TRUE.
                   ptrOut%rotation = ptrIn%rotation
                   select case ( dir )               
                   case ( 1 ) ! to_pe's eastern halo
                      select case(ptrIn%rotation)
                      case (ZERO)  !  W -> E
                         isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = jshift
                      case (NINETY) ! S -> E
                         isoff1 = 0;      ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
                      end select
                   case ( 2 ) ! to_pe's south-eastearn halo
                      select case(ptrIn%rotation)
                      case (ZERO)  
                         isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
                      case (NINETY) 
                         isoff1 = jshift; ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
                      case (MINUS_NINETY) 
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
                      end select
                   case ( 3 ) ! to_pe's southern halo
                      select case(ptrIn%rotation)
                      case (ZERO)  !  N -> S
                         isoff1 = 0;      ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
                      case (MiNUS_NINETY) ! E -> S
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = ishift
                      end select
                   case ( 4 ) ! to_pe's south-westearn halo
                      select case(ptrIn%rotation)
                      case (ZERO)  
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
                      case (NINETY) 
                         isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
                      case (MINUS_NINETY) 
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = ishift; jeoff1 = ishift
                      end select
                   case ( 5 ) ! to_pe's western halo
                      select case(ptrIn%rotation)
                      case (ZERO)  !  E -> W
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = jshift
                      case (NINETY) ! N -> W
                         isoff1 = 0;      ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
                      end select
                   case ( 6 ) ! to_pe's north-westearn halo
                      select case(ptrIn%rotation)
                      case (ZERO)  
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = jshift; jeoff1 = jshift
                      case (NINETY) 
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
                      case (MINUS_NINETY) 
                         isoff1 = jshift; ieoff1 = jshift; jsoff1 = ishift; jeoff1 = ishift
                      end select
                   case ( 7 ) ! to_pe's northern halo
                      select case(ptrIn%rotation)
                      case (ZERO)  !  S -> N
                         isoff1 = 0;      ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
                      case (MINUS_NINETY) ! W -> N
                         isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = ishift
                      end select
                   case ( 8 ) ! to_pe's north-eastearn halo
                      select case(ptrIn%rotation)
                      case (ZERO)  
                         isoff1 = ishift; ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
                      case (NINETY) 
                         isoff1 = 0;      ieoff1 = 0;      jsoff1 = ishift; jeoff1 = ishift
                      case (MINUS_NINETY) 
                         isoff1 = jshift; ieoff1 = jshift; jsoff1 = 0;      jeoff1 = 0
                      end select
                   end select
                   Ptrout%is(1) = Ptrin%is(1) + isoff1
                   Ptrout%ie(1) = Ptrin%ie(1) + ieoff1
                   Ptrout%js(1) = Ptrin%js(1) + jsoff1
                   Ptrout%je(1) = Ptrin%je(1) + jeoff1
                end if

                ptrIn  => domain_in %list(m)%recv(tMe,tNbr,dir)
                ptrOut => domain_out%list(m)%recv(tMe,tNbr,dir)
                PtrOut%is_refined = PtrIn%is_refined
                ! only set overlapping between tiles for recv ( ptrOut%overlap(1) is false )
                if(ptrIn%overlap(1)) then  
                   ptrOut%overlap(1) = .TRUE.
                   ptrOut%rotation = ptrIn%rotation
                   select case ( dir )               
                   case ( 1 ) ! E
                      isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = jshift
                   case ( 2 ) ! SE
                      isoff1 = ishift; ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
                   case ( 3 ) ! S
                      isoff1 = 0;      ieoff1 = ishift; jsoff1 = 0;      jeoff1 = 0
                   case ( 4 ) ! SW
                      isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = 0
                   case ( 5 ) ! W 
                      isoff1 = 0;      ieoff1 = 0;      jsoff1 = 0;      jeoff1 = jshift
                   case ( 6 ) ! NW
                     isoff1 = 0;       ieoff1 = 0;      jsoff1 = jshift; jeoff1 = jshift
                   case ( 7 ) ! N
                      isoff1 = 0;      ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
                   case ( 8 ) ! NE
                      isoff1 = ishift; ieoff1 = ishift; jsoff1 = jshift; jeoff1 = jshift
                   end select
                   Ptrout%is(1) = Ptrin%is(1) + isoff1
                   Ptrout%ie(1) = Ptrin%ie(1) + ieoff1
                   Ptrout%js(1) = Ptrin%js(1) + jsoff1
                   Ptrout%je(1) = Ptrin%je(1) + jeoff1
                   Ptrout%isMe  = Ptrin%isMe  + isoff1
                   Ptrout%ieMe  = Ptrin%ieMe  + ieoff1
                   Ptrout%jsMe  = Ptrin%jsMe  + jsoff1
                   Ptrout%jeMe  = Ptrin%jeMe  + jeoff1          
                   !--- figure out neighbor offset when there is refinemnet
                   if( PtrIn%is_refined(1) ) then
                      isoff2 = isoff1; ieoff2 = ieoff1; jsoff2 = jsoff1; jeoff2 = jeoff1
                      if( ptrOut%rotation == NINETY) then
                         select case ( dir )                                       
                         case ( 1, 5 ) !S -> E, N -> W
                            isoff2 = jsoff1; ieoff2 = jeoff1; jsoff2 = isoff1; jeoff2 = ieoff1
                         case ( 2 ) 
                            isoff2 = jshift; ieoff2 = jshift; jsoff2 = ishift; jeoff2 = ishift
                         case ( 4 ) 
                            isoff2 = jshift; ieoff2 = jshift; jsoff2 = 0;      jeoff2 = 0
                         case ( 6 ) 
                            isoff2 = 0;      ieoff2 = 0;      jsoff2 = 0;      jeoff2 = 0
                         case ( 8 ) 
                            isoff2 = 0;      ieoff2 = 0;      jsoff2 = ishift; jeoff2 = ishift
                         end select
                      else if( ptrOut%rotation == MINUS_NINETY) then
                         select case ( dir )                                       
                         case ( 3, 7 ) !E -> S, W -> N
                            isoff2 = jsoff1; ieoff2 = jeoff1; jsoff2 = isoff1; jeoff2 = ieoff1
                         case ( 2 )
                            isoff2 = 0;      ieoff2 = 0;      jsoff2 = 0;      jeoff2 = 0 
                         case ( 4 ) 
                            isoff2 = 0;      ieoff2 = 0;      jsoff2 = ishift; jeoff2 = ishift
                         case ( 6 ) 
                            isoff2 = jshift; ieoff2 = jshift; jsoff2 = ishift; jeoff2 = ishift
                         case ( 8 ) 
                            isoff2 = jshift; ieoff2 = jshift; jsoff2 = 0;      jeoff2 = 0
                         end select
                      end if
                      !--- calculate overlapping
                      Ptrout%is(1) = Ptrin%is(1) + isoff2
                      Ptrout%ie(1) = Ptrin%ie(1) + ieoff2
                      Ptrout%js(1) = Ptrin%js(1) + jsoff2
                      Ptrout%je(1) = Ptrin%je(1) + jeoff2
                   end if
                end if
             end do ! end do tNbr = 1, ntileNbr
          end do ! end do tMe = 1, ntileMe
       end do ! end do list = 0, nlist-1
       if( ANY(domain_out%list(m)%send(:,:,:)%overlap(1)) ) domain_out%list(m)%overlap = .true.
       if( ANY(domain_out%list(m)%recv(:,:,:)%overlap(1)) ) domain_out%list(m)%overlap = .true.
       do tMe = 1, ntileMe
       if( ANY(domain_out%list(m)%recv(tMe,:,1)%overlap(1)) ) domain_out%x(tMe)%loffset = 0
       if( ANY(domain_out%list(m)%recv(tMe,:,7)%overlap(1)) ) domain_out%y(tMe)%loffset = 0
    end do
    end do

    !--- define the refined overlapping and group the overlapping 
    !--- between same tiles( in different direction or different pe ). 
    !--- set up the index position in the return buffer.
    call define_refine_overlap( domain_out )  

    domain_out%initialized = .true.

  end subroutine set_contact_point

  !--- set up the overlapping for boundary check if the domain is symmetry. The check will be 
  !--- done on current pe for east boundary for E-cell, north boundary for N-cell, 
  !--- East and North boundary for C-cell
  subroutine set_checkbound( domain, check )
    type(domain2d),    intent(in) :: domain
    type(boundary), intent(inout) :: check
    integer                       :: nlist, list, m, ntileMe, ntileNbr, tMe, tNbr, count
    integer, parameter            :: MAXCOUNT = 100
    integer, dimension(MAXCOUNT)  :: dir, rotation, is, ie, js, je, tileMe

    if(domain%position == CENTER .OR. .NOT. domain%symmetry ) return

    nlist = size(domain%list(:))
    ntileMe = size(domain%x(:))
    allocate(check%send(0:nlist-1))
    allocate(check%recv(0:nlist-1))
    check%send(:)%count = 0
    check%recv(:)%count = 0
    !--- loop over the list of domains to find the boundary overlap for send
    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       ntileNbr = size(domain%list(m)%x(:))
       if( .NOT.  domain%list(m)%overlap ) cycle
       count = 0
       do tNbr = 1, ntileNbr
          do tMe = 1, ntileMe
             ! comparing east direction on currently pe
             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. domain%list(m)%send(tNbr,tMe,1)%overlap(1) ) then  
                if(.NOT. domain%list(m)%send(tNbr,tMe,1)%is_refined(1) ) then
                   count = count + 1
                   rotation(count) = domain%list(m)%send(tNbr,tMe,1)%rotation
                   tileMe(count) = tMe
                   dir(count) = 1
                   select case( rotation(count) ) 
                   case( ZERO ) ! W -> E
                      is(count) = domain%list(m)%send(tNbr,tMe,1)%is(1) - 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%send(tNbr,tMe,1)%js(1)
                      je(count) = domain%list(m)%send(tNbr,tMe,1)%je(1)
                   case( NINETY ) ! S -> E
                      is(count) = domain%list(m)%send(tNbr,tMe,1)%is(1)
                      ie(count) = domain%list(m)%send(tNbr,tMe,1)%ie(1)
                      js(count) = domain%list(m)%send(tNbr,tMe,1)%js(1) - 1
                      je(count) = js(count)
                   end select
                end if
             end if
             ! comparing north direction on currently pe
             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. domain%list(m)%send(tNbr,tMe,7)%overlap(1) ) then  
                !--- excluded folded edge since folded edge is updated on the boundary.
                if( (.NOT. domain%list(m)%send(tNbr,tMe,7)%is_refined(1) ) .AND. &
                     domain%list(m)%send(tNbr,tMe,7)%rotation .NE. ONE_HUNDRED_EIGHTY ) then
                   count=count+1
                   dir(count) = 4
                   rotation(count) = domain%list(m)%send(tNbr,tMe,7)%rotation
                   tileMe(count) = tMe
                   select case( rotation(count) ) 
                   case( ZERO ) ! S->N
                      is(count) = domain%list(m)%send(tNbr,tMe,7)%is(1)
                      ie(count) = domain%list(m)%send(tNbr,tMe,7)%ie(1)
                      js(count) = domain%list(m)%send(tNbr,tMe,7)%js(1) - 1
                      je(count) = js(count)
                   case( MINUS_NINETY ) ! W->N
                      is(count) = domain%list(m)%send(tNbr,tMe,7)%is(1) - 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%send(tNbr,tMe,7)%js(1)
                      je(count) = domain%list(m)%send(tNbr,tMe,7)%je(1)
                   end select
                end if
             end if
          end do ! end do tNbr = 1, ntileNbr
       end do ! end do tMe = 1, ntileMe
       if(count>0) then
          check%send(m)%count = count
          allocate(check%send(m)%is(count),  check%send(m)%ie(count) )
          allocate(check%send(m)%js(count),  check%send(m)%je(count) )
          allocate(check%send(m)%dir(count), check%send(m)%rotation(count) )
          allocate(check%send(m)%tileMe(count))
          check%send(m)%is(:)       = is(1:count)
          check%send(m)%ie(:)       = ie(1:count)
          check%send(m)%js(:)       = js(1:count)
          check%send(m)%je(:)       = je(1:count)
          check%send(m)%dir(:)      = dir(1:count)
          check%send(m)%tileMe(:)   = tileMe(1:count)
          check%send(m)%rotation(:) = rotation(1:count)
       end if
    end do ! end  do list = 0, nlist

    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       ntileNbr = size(domain%list(m)%x(:))
       if( .NOT.  domain%list(m)%overlap ) cycle
       count = 0
       do tMe = 1, ntileMe
          do tNbr = 1, ntileNbr
             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. domain%list(m)%recv(tMe,tNbr, 1)%overlap(1) ) then  
                if(.NOT. domain%list(m)%recv(tMe,tNbr, 1)%is_refined(1) ) then
                   count=count+1
                   dir(count) = 1
                   rotation(count) = domain%list(m)%recv(tMe,tNbr,1)%rotation
                   tileMe(count) = tMe
                   is(count) = domain%list(m)%recv(tMe,tNbr,1)%is(1) - 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%recv(tMe,tNbr,1)%js(1)
                   je(count) = domain%list(m)%recv(tMe,tNbr,1)%je(1)
                end if
             end if
             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. domain%list(m)%recv(tMe,tNbr,7)%overlap(1) ) then  
                !--- excluded folded edge since folded edge is updated on the boundary.
                if((.NOT.domain%list(m)%recv(tMe,tNbr,7)%is_refined(1) ) .AND. &
                     domain%list(m)%recv(tMe,tNbr,7)%rotation .NE. ONE_HUNDRED_EIGHTY ) then
                   count=count+1
                   dir(count) = 4
                   rotation(count) = domain%list(m)%recv(tMe,tNbr,7)%rotation
                   tileMe(count) = tMe
                   is(count) = domain%list(m)%recv(tMe,tNbr,7)%is(1)
                   ie(count) = domain%list(m)%recv(tMe,tNbr,7)%ie(1)
                   js(count) = domain%list(m)%recv(tMe,tNbr,7)%js(1) - 1
                   je(count) = js(count)
                end if
             end if
          end do ! end do tNbr = 1, ntileNbr
       end do ! end do tMe = 1, ntileMe
       if(count>0) then
          check%recv(m)%count = count
          allocate(check%recv(m)%is(count),     check%recv(m)%ie(count) )
          allocate(check%recv(m)%js(count),     check%recv(m)%je(count) )
          allocate(check%recv(m)%dir(count)  )
          allocate(check%recv(m)%tileMe(count), check%recv(m)%rotation(count) )
          check%recv(m)%is(:)       = is(1:count)
          check%recv(m)%ie(:)       = ie(1:count)
          check%recv(m)%js(:)       = js(1:count)
          check%recv(m)%je(:)       = je(1:count)
          check%recv(m)%dir(:)      = dir(1:count)
          check%recv(m)%tileMe(:)   = tileMe(1:count)
          check%recv(m)%rotation(:) = rotation(1:count)
       end if
    end do ! end  do list = 0, nlist

  end subroutine set_checkbound

  !#############################################################################
  !--- set up the overlapping for boundary if the domain is symmetry.
  subroutine set_boundary( domain, bound )
    type(domain2d),    intent(in) :: domain
    type(boundary), intent(inout) :: bound
    integer                       :: nlist, list, m, n, l, count, tMe, tNbr, ntileMe, ntileNbr, dr
    integer, parameter            :: MAXCOUNT = 100
    integer, dimension(MAXCOUNT)  :: dir, rotation, is, ie, js, je, isMe, ieMe, jsMe, jeMe, tileMe
    integer, dimension(size(domain%x(:)), 4)           :: nrecv
    integer, dimension(size(domain%x(:)), 4, MAXCOUNT) :: isl, iel, jsl, jel, islMe, ielMe, jslMe, jelMe


    if(domain%position == CENTER .OR. .NOT. domain%symmetry ) return

    nlist = size(domain%list(:))
    ntileMe = size(domain%x(:))
    allocate(bound%send(0:nlist-1))
    allocate(bound%recv(0:nlist-1))
    bound%send(:)%count = 0
    bound%recv(:)%count = 0
    !--- loop over the list of domains to find the boundary overlap for send

    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       ntileNbr = size(domain%list(m)%x(:))
       if( .NOT.  domain%list(m)%overlap ) cycle
       count = 0
       do tNbr = 1, ntileNbr
          do tMe = 1, ntileMe
             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. &
                domain%list(m)%send(tNbr,tMe,1)%overlap(1)) then                  ! east
                count=count+1
                dir(count) = 1
                rotation(count) = domain%list(m)%send(tNbr,tMe,1)%rotation
                tileMe(count) = tMe
                select case( rotation(count) ) 
                case( ZERO ) ! W -> E
                   is(count) = domain%list(m)%send(tNbr,tMe,1)%is(1) - 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%send(tNbr,tMe,1)%js(1)
                   je(count) = domain%list(m)%send(tNbr,tMe,1)%je(1)
                case( NINETY ) ! S -> E
                   is(count) = domain%list(m)%send(tNbr,tMe,1)%is(1)
                   ie(count) = domain%list(m)%send(tNbr,tMe,1)%ie(1)
                   js(count) = domain%list(m)%send(tNbr,tMe,1)%js(1) - 1
                   je(count) = js(count)
                end select
             end if
             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. &
                domain%list(m)%send(tNbr,tMe,3)%overlap(1)) then                 ! south
                count=count+1
                dir(count) = 2
                rotation(count) = domain%list(m)%send(tNbr,tMe,3)%rotation
                tileMe(count) = tMe
                select case( rotation(count) ) 
                case( ZERO ) ! N->S
                   is(count) = domain%list(m)%send(tNbr,tMe,3)%is(1)
                   ie(count) = domain%list(m)%send(tNbr,tMe,3)%ie(1)
                   js(count) = domain%list(m)%send(tNbr,tMe,3)%je(1) + 1
                   je(count) = js(count)
                case( MINUS_NINETY ) ! E->S
                   is(count) = domain%list(m)%send(tNbr,tMe,3)%ie(1) + 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%send(tNbr,tMe,3)%js(1)
                   je(count) = domain%list(m)%send(tNbr,tMe,3)%je(1)
                end select
             end if
             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. &
                domain%list(m)%send(tNbr,tMe,5)%overlap(1)) then                  ! west
                count=count+1
                dir(count) = 3
                rotation(count) = domain%list(m)%send(tNbr,tMe,5)%rotation
                tileMe(count) = tMe
                select case( rotation(count) ) 
                case( ZERO ) ! E->W
                   is(count) = domain%list(m)%send(tNbr,tMe,5)%ie(1) + 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%send(tNbr,tMe,5)%js(1)
                   je(count) = domain%list(m)%send(tNbr,tMe,5)%je(1)
                case( NINETY ) ! N->W
                   is(count) = domain%list(m)%send(tNbr,tMe,5)%is(1)
                   ie(count) = domain%list(m)%send(tNbr,tMe,5)%ie(1)
                   js(count) = domain%list(m)%send(tNbr,tMe,5)%je(1) + 1
                   je(count) = js(count)
                end select
             end if
             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. &
                domain%list(m)%send(tNbr,tMe,7)%overlap(1)) then                 ! south
                count=count+1
                dir(count) = 4
                rotation(count) = domain%list(m)%send(tNbr,tMe,7)%rotation
                tileMe(count) = tMe
                select case( rotation(count) ) 
                case( ZERO ) ! S->N
                   is(count) = domain%list(m)%send(tNbr,tMe,7)%is(1)
                   ie(count) = domain%list(m)%send(tNbr,tMe,7)%ie(1)
                   js(count) = domain%list(m)%send(tNbr,tMe,7)%js(1) - 1
                   je(count) = js(count)
                case( MINUS_NINETY ) ! W->N
                   is(count) = domain%list(m)%send(tNbr,tMe,7)%is(1) - 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%send(tNbr,tMe,7)%js(1)
                   je(count) = domain%list(m)%send(tNbr,tMe,7)%je(1)
                end select
             end if
          end do ! end do tNbr = 1, ntileNbr
       end do ! end do tMe = 1, ntileMe
       if(count>0) then
          bound%send(m)%count = count
          allocate(bound%send(m)%is(count),  bound%send(m)%ie(count) )
          allocate(bound%send(m)%js(count),  bound%send(m)%je(count) )
          allocate(bound%send(m)%dir(count), bound%send(m)%rotation(count) )
          allocate(bound%send(m)%tileMe(count))
          bound%send(m)%is(:)       = is(1:count)
          bound%send(m)%ie(:)       = ie(1:count)
          bound%send(m)%js(:)       = js(1:count)
          bound%send(m)%je(:)       = je(1:count)
          bound%send(m)%dir(:)      = dir(1:count)
          bound%send(m)%tileMe(:)   = tileMe(1:count)
          bound%send(m)%rotation(:) = rotation(1:count)
       end if
    end do ! end  do list = 0, nlist

    !--- loop over the list of domains to find the boundary overlap for recv

    nrecv = 0    

    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       ntileNbr = size(domain%list(m)%x(:))
       if( .NOT.  domain%list(m)%overlap ) cycle
       count = 0
       do tMe = 1, ntileMe
          do tNbr = 1, ntileNbr
             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. &
                domain%list(m)%recv(tMe,tNbr,1)%overlap(1)) then                  ! east
                count=count+1
                dir(count) = 1
                rotation(count) = domain%list(m)%recv(tMe,tNbr,1)%rotation
                isMe(count) = domain%list(m)%recv(tMe,tNbr,1)%isMe - 1
                ieMe(count) = isMe(count)
                jsMe(count) = domain%list(m)%recv(tMe,tNbr,1)%jsMe
                jeMe(count) = domain%list(m)%recv(tMe,tNbr,1)%jeMe
                tileMe(count) = tMe
                if( domain%list(m)%recv(tMe,tNbr,1)%is_refined(1)) then
                   select case( rotation(count) ) 
                   case( ZERO ) ! W -> E
                      is(count) = domain%list(m)%recv(tMe,tNbr,1)%is(1) - 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%recv(tMe,tNbr,1)%js(1)
                      je(count) = domain%list(m)%recv(tMe,tNbr,1)%je(1)
                   case( NINETY ) ! S -> E
                      is(count) = domain%list(m)%recv(tMe,tNbr,1)%is(1)
                      ie(count) = domain%list(m)%recv(tMe,tNbr,1)%ie(1)
                      js(count) = domain%list(m)%recv(tMe,tNbr,1)%js(1) - 1
                      je(count) = js(count)
                   end select
                else
                   is(count) = domain%list(m)%recv(tMe,tNbr,1)%is(1) - 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%recv(tMe,tNbr,1)%js(1)
                   je(count) = domain%list(m)%recv(tMe,tNbr,1)%je(1)
                end if
                nrecv(tMe, 1) = nrecv(tMe,1) + 1
                islMe(tMe,1,nrecv(tMe, 1)) = isMe(count)
                ielMe(tMe,1,nrecv(tMe, 1)) = ieMe(count)
                jslMe(tMe,1,nrecv(tMe, 1)) = jsMe(count)
                jelMe(tMe,1,nrecv(tMe, 1)) = jeMe(count)
                isl  (tMe,1,nrecv(tMe, 1)) = is  (count)
                iel  (tMe,1,nrecv(tMe, 1)) = ie  (count)
                jsl  (tMe,1,nrecv(tMe, 1)) = js  (count)
                jel  (tMe,1,nrecv(tMe, 1)) = je  (count)
             end if

             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. &
                domain%list(m)%recv(tMe,tNbr,3)%overlap(1)) then                 ! south
                count=count+1
                dir(count) = 2
                rotation(count) = domain%list(m)%recv(tMe,tNbr,3)%rotation
                isMe(count) = domain%list(m)%recv(tMe,tNbr,3)%isMe
                ieMe(count) = domain%list(m)%recv(tMe,tNbr,3)%ieMe
                jsMe(count) = domain%list(m)%recv(tMe,tNbr,3)%jeMe + 1
                jeMe(count) = jsMe(count)
                tileMe(count) = tMe
                if( domain%list(m)%recv(tMe,tNbr,3)%is_refined(1)) then
                   select case( rotation(count) ) 
                   case( ZERO ) ! N->S
                      is(count) = domain%list(m)%recv(tMe,tNbr,3)%is(1)
                      ie(count) = domain%list(m)%recv(tMe,tNbr,3)%ie(1)
                      js(count) = domain%list(m)%recv(tMe,tNbr,3)%je(1) + 1
                      je(count) = js(count)
                   case( MINUS_NINETY ) ! E->S
                      is(count) = domain%list(m)%recv(tMe,tNbr,3)%ie(1) + 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%recv(tMe,tNbr,3)%js(1)
                      je(count) = domain%list(m)%recv(tMe,tNbr,3)%je(1)
                   end select
                else
                   is(count) = domain%list(m)%recv(tMe,tNbr,3)%is(1)
                   ie(count) = domain%list(m)%recv(tMe,tNbr,3)%ie(1)
                   js(count) = domain%list(m)%recv(tMe,tNbr,3)%je(1) + 1
                   je(count) = js(count)
                end if
                nrecv(tMe, 2) = nrecv(tMe,2) + 1
                islMe(tMe,2,nrecv(tMe, 2)) = isMe(count)
                ielMe(tMe,2,nrecv(tMe, 2)) = ieMe(count)
                jslMe(tMe,2,nrecv(tMe, 2)) = jsMe(count)
                jelMe(tMe,2,nrecv(tMe, 2)) = jeMe(count)
                isl  (tMe,2,nrecv(tMe, 2)) = is  (count)
                iel  (tMe,2,nrecv(tMe, 2)) = ie  (count)
                jsl  (tMe,2,nrecv(tMe, 2)) = js  (count)
                jel  (tMe,2,nrecv(tMe, 2)) = je  (count)
             end if

             if( (domain%position == EAST .OR. domain%position == CORNER) .AND. &
                domain%list(m)%recv(tMe,tNbr,5)%overlap(1)) then                  ! west
                count=count+1
                dir(count) = 3
                rotation(count) = domain%list(m)%recv(tMe,tNbr,5)%rotation
                isMe(count) = domain%list(m)%recv(tMe,tNbr,5)%ieMe + 1 
                ieMe(count) = isMe(count)
                jsMe(count) = domain%list(m)%recv(tMe,tNbr,5)%jsMe
                jeMe(count) = domain%list(m)%recv(tMe,tNbr,5)%jeMe
                tileMe(count) = tMe
                if( domain%list(m)%recv(tMe,tNbr,5)%is_refined(1)) then
                   select case( rotation(count) ) 
                   case( ZERO ) ! E->W
                      is(count) = domain%list(m)%recv(tMe,tNbr,5)%ie(1) + 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%recv(tMe,tNbr,5)%js(1)
                      je(count) = domain%list(m)%recv(tMe,tNbr,5)%je(1)
                   case( NINETY ) ! S -> E
                      is(count) = domain%list(m)%recv(tMe,tNbr,5)%is(1)
                      ie(count) = domain%list(m)%recv(tMe,tNbr,5)%ie(1)
                      js(count) = domain%list(m)%recv(tMe,tNbr,5)%je(1) + 1
                      je(count) = js(count)
                   end select
                else
                   is(count) = domain%list(m)%recv(tMe,tNbr,5)%ie(1) + 1
                   ie(count) = is(count)
                   js(count) = domain%list(m)%recv(tMe,tNbr,5)%js(1)
                   je(count) = domain%list(m)%recv(tMe,tNbr,5)%je(1)
                end if
                nrecv(tMe, 3) = nrecv(tMe,3) + 1
                islMe(tMe,3,nrecv(tMe, 3)) = isMe(count)
                ielMe(tMe,3,nrecv(tMe, 3)) = ieMe(count)
                jslMe(tMe,3,nrecv(tMe, 3)) = jsMe(count)
                jelMe(tMe,3,nrecv(tMe, 3)) = jeMe(count)
                isl  (tMe,3,nrecv(tMe, 3)) = is  (count)
                iel  (tMe,3,nrecv(tMe, 3)) = ie  (count)
                jsl  (tMe,3,nrecv(tMe, 3)) = js  (count)
                jel  (tMe,3,nrecv(tMe, 3)) = je  (count)
             end if

             if( (domain%position == NORTH .OR. domain%position == CORNER) .AND. &
                domain%list(m)%recv(tMe,tNbr,7)%overlap(1)) then                 ! north
                count=count+1
                dir(count) = 4
                rotation(count) = domain%list(m)%recv(tMe,tNbr,7)%rotation
                isMe(count) = domain%list(m)%recv(tMe,tNbr,7)%isMe
                ieMe(count) = domain%list(m)%recv(tMe,tNbr,7)%ieMe
                jsMe(count) = domain%list(m)%recv(tMe,tNbr,7)%jsMe - 1
                jeMe(count) = jsMe(count)
                tileMe(count) = tMe
                if( domain%list(m)%recv(tMe,tNbr,7)%is_refined(1)) then
                   select case( rotation(count) ) 
                   case( ZERO ) ! S->N
                      is(count) = domain%list(m)%recv(tMe,tNbr,7)%is(1)
                      ie(count) = domain%list(m)%recv(tMe,tNbr,7)%ie(1)
                      js(count) = domain%list(m)%recv(tMe,tNbr,7)%js(1) - 1
                      je(count) = js(count)
                   case( MINUS_NINETY ) ! W->N
                      is(count) = domain%list(m)%recv(tMe,tNbr,7)%is(1) - 1
                      ie(count) = is(count)
                      js(count) = domain%list(m)%recv(tMe,tNbr,7)%js(1)
                      je(count) = domain%list(m)%recv(tMe,tNbr,7)%je(1)
                   end select
                else
                   is(count) = domain%list(m)%recv(tMe,tNbr,7)%is(1)
                   ie(count) = domain%list(m)%recv(tMe,tNbr,7)%ie(1)
                   js(count) = domain%list(m)%recv(tMe,tNbr,7)%js(1) - 1
                   je(count) = js(count)
                end if
                nrecv(tMe, 4) = nrecv(tMe,4) + 1
                islMe(tMe,4,nrecv(tMe, 4)) = isMe(count)
                ielMe(tMe,4,nrecv(tMe, 4)) = ieMe(count)
                jslMe(tMe,4,nrecv(tMe, 4)) = jsMe(count)
                jelMe(tMe,4,nrecv(tMe, 4)) = jeMe(count)
                isl  (tMe,4,nrecv(tMe, 4)) = is  (count)
                iel  (tMe,4,nrecv(tMe, 4)) = ie  (count)
                jsl  (tMe,4,nrecv(tMe, 4)) = js  (count)
                jel  (tMe,4,nrecv(tMe, 4)) = je  (count)
             end if
          end do ! end do tNbr = 1, ntileNbr
       end do ! end do tMe = 1, ntileMe
       if(count>0) then
          bound%recv(m)%count = count
          allocate(bound%recv(m)%isMe(count),   bound%recv(m)%ieMe(count) )
          allocate(bound%recv(m)%jsMe(count),   bound%recv(m)%jeMe(count) )
          allocate(bound%recv(m)%is(count),     bound%recv(m)%ie(count) )
          allocate(bound%recv(m)%js(count),     bound%recv(m)%je(count) )
          allocate(bound%recv(m)%dir(count),    bound%recv(m)%index(count)  )
          allocate(bound%recv(m)%tileMe(count), bound%recv(m)%rotation(count) )
          bound%recv(m)%isMe(:)     = isMe(1:count)
          bound%recv(m)%ieMe(:)     = ieMe(1:count)
          bound%recv(m)%jsMe(:)     = jsMe(1:count)
          bound%recv(m)%jeMe(:)     = jeMe(1:count)
          bound%recv(m)%is(:)       = is(1:count)
          bound%recv(m)%ie(:)       = ie(1:count)
          bound%recv(m)%js(:)       = js(1:count)
          bound%recv(m)%je(:)       = je(1:count)
          bound%recv(m)%dir(:)      = dir(1:count)
          bound%recv(m)%tileMe(:)   = tileMe(1:count)
          bound%recv(m)%rotation(:) = rotation(1:count)
       end if
    end do ! end  do list = 0, nlist
     
    !--- find the boundary index for each contact within the east boundary
    do list = 0, nlist-1
       m = mod( domain%pos+list, nlist )
       do n = 1, bound%recv(m)%count
          tMe = bound%recv(m)%tileMe(n)
          dr = bound%recv(m)%dir(n)
          bound%recv(m)%index(n) = 1
          do l = 1, nrecv(tMe,dr)
             if(dr == 1 .OR. dr == 3) then  ! EAST, WEST
                if( bound%recv(m)%jsMe(n) > jslMe(tMe, dr, l) ) then
                   bound%recv(m)%index(n) = bound%recv(m)%index(n)                  + &
                                            max(abs(jel(tMe, dr, l)-jsl(tMe, dr, l)), &
                                            abs(iel(tMe, dr, l)-isl(tMe, dr, l))) + 1 - domain%y(1)%shift
                end if
             else                             ! South, North
                if( bound%recv(m)%isMe(n) > islMe(tMe, dr, l) ) then
                   bound%recv(m)%index(n) = bound%recv(m)%index(n)                  + &
                                            max(abs(jel(tMe, dr, l)-jsl(tMe, dr, l)), &
                                            abs(iel(tMe, dr, l)-isl(tMe, dr, l))) + 1 - domain%x(1)%shift
                end if
             end if
          end do
       end do
    end do
   

  end subroutine set_boundary


  !#############################################################################

  subroutine fill_corner_contact(eCont, sCont, wCont, nCont, isg, ieg, jsg, jeg, numR, numS, tileRecv, tileSend, &
                                 is1Recv, ie1Recv, js1Recv, je1Recv, is2Recv, ie2Recv, js2Recv, je2Recv,         &
                                 is1Send, ie1Send, js1Send, je1Send, is2Send, ie2Send, js2Send, je2Send,         &
                                 align1Recv, align2Recv, align1Send, align2Send,                                 &
                                 whalo, ehalo, shalo, nhalo, tileMe)
    type(contact_type), dimension(:), intent(in) :: eCont, sCont, wCont, nCont
    integer, dimension(:),            intent(in) :: isg, ieg, jsg, jeg
    integer,                       intent(inout) :: numR, numS
    integer, dimension(:),         intent(inout) :: tileRecv, tileSend
    integer, dimension(:),         intent(inout) :: is1Recv, ie1Recv, js1Recv, je1Recv
    integer, dimension(:),         intent(inout) :: is2Recv, ie2Recv, js2Recv, je2Recv
    integer, dimension(:),         intent(inout) :: is1Send, ie1Send, js1Send, je1Send
    integer, dimension(:),         intent(inout) :: is2Send, ie2Send, js2Send, je2Send
    integer, dimension(:),         intent(inout) :: align1Recv, align2Recv, align1Send, align2Send
    integer,                          intent(in) :: tileMe, whalo, ehalo, shalo, nhalo
    integer                                      :: is1, ie1, js1, je1, is2, ie2, js2, je2
    integer                                      :: tn, tc, n, m
    logical                                      :: found_corner

    found_corner = .false.
    !--- southeast for recving
    if(eCont(tileMe)%ncontact > 0) then     
       if(eCont(tileMe)%js1(1) == jsg(tileMe) ) then
          tn = eCont(tileMe)%tile(1)
          if(econt(tileMe)%js2(1) > jsg(tn) ) then  ! the corner tile is tn.
             if( econt(tileMe)%js2(1) - jsg(tn) < shalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southeast tile for recv 1 is not tiled properly")
             found_corner = .true.; tc = tn
             is1 = eCont(tileMe)%ie1(1) + 1; je1 = eCont(tileMe)%js1(1) - 1
             is2 = eCont(tileMe)%is2(1);     je2 = eCont(tileMe)%js2(1) - 1
          else if(sCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
             if(sCont(tn)%is1(1) == isg(tn)) then ! corner is nc.
                found_corner = .true.; tc = sCont(tn)%tile(1)
                is1 = eCont(tileMe)%ie1(1) + 1; je1 = eCont(tileMe)%js1(1) - 1
                is2 = sCont(tn)%is2(1);         je2 = sCont(tn)%je2(1)
             end if
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found,
       n = sCont(tileMe)%ncontact
       if( n > 0) then
          if( sCont(tileMe)%ie1(n) == ieg(tileMe)) then
             tn = sCont(tileMe)%tile(n)
             if(scont(tileMe)%ie2(n) < ieg(tn) ) then  ! the corner tile is tn.
                if(ieg(tn) - scont(tileMe)%ie2(n) < ehalo ) call mpp_error(FATAL, &
                     "mpp_domains_define.inc: southeast tile for recv 2 is not tiled properly")
                found_corner = .true.; tc = tn
                is1 = sCont(tileMe)%ie1(n) + 1; je1 = sCont(tileMe)%js1(n) - 1
                is2 = sCont(tileMe)%ie2(n) + 1; je2 = sCont(tileMe)%je2(n)
             else if(eCont(tn)%ncontact >0) then ! the corner tile may be east tile of tn.
                m = eCont(tn)%ncontact
                if(eCont(tn)%je1(m) == jeg(tn)) then ! corner is nc.
                   found_corner = .true.; tc = eCont(tn)%tile(m)
                   is1 = sCont(tileMe)%ie1(n) + 1; je1 = sCont(tileMe)%js1(n) - 1
                   is2 = eCont(tn)%is2(m);         je2 = eCont(tn)%je2(m)
                end if
             end if
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tc; align1Recv(numR) = SOUTH_EAST;  align2Recv(numR) = NORTH_WEST
       is1Recv(numR) = is1;             ie1Recv(numR) = is1 + ehalo - 1
       js1Recv(numR) = je1 - shalo + 1; je1Recv(numR) = je1
       is2Recv(numR) = is2;             ie2Recv(numR) = is2 + ehalo - 1
       js2Recv(numR) = je2 - shalo + 1; je2Recv(numR) = je2  
    end if
 
    !--- southwest for recving
    found_corner = .false.
    if(wCont(tileMe)%ncontact > 0) then
       if(wCont(tileMe)%js1(1) == jsg(tileMe) ) then
          tn = wCont(tileMe)%tile(1)
          if(wcont(tileMe)%js2(1) > jsg(tn) ) then  ! the corner tile is tn.
             if( wcont(tileMe)%js2(1) - jsg(tn) < shalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southwest tile for recv 1 is not tiled properly")
             found_corner = .true.; tc = tn
             ie1 = wCont(tileMe)%is1(1) - 1; je1 = wCont(tileMe)%js1(1) - 1
             ie2 = wCont(tileMe)%is2(1);     je2 = wCont(tileMe)%js2(1) - 1
          else if(sCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
             n = sCont(tn)%ncontact
             if(sCont(tn)%ie1(n) == ieg(tn)) then ! corner is nc.
                found_corner = .true.;  tc = sCont(tn)%tile(n)
                ie1 = wCont(tileMe)%is1(1) - 1; je1 = wCont(tileMe)%js1(1) - 1
                ie2 = sCont(tn)%ie2(1);         je2 = sCont(tn)%je2(1)
             end if
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found,
       n = sCont(tileMe)%ncontact
       if( n > 0) then
          if( sCont(tileMe)%is1(1) == isg(tileMe)) then
             tn = sCont(tileMe)%tile(1)
             if(sCont(tileMe)%is2(1) > isg(tn) ) then  ! the corner tile is tn.
                if( scont(tileMe)%is2(1)-isg(tn) < whalo ) call mpp_error(FATAL, &
                     "mpp_domains_define.inc: southwest tile for recv 1 is not tiled properly")
                found_corner = .true.; tc = tn
                ie1 = sCont(tileMe)%is1(1) - 1; je1 = sCont(tileMe)%js1(1) - 1
                ie2 = sCont(tileMe)%is2(1) - 1; je2 = sCont(tileMe)%js2(1)
             else if(wCont(tn)%ncontact >0) then ! the corner tile may be west tile of tn.
                m = wCont(tn)%ncontact
                if(wCont(tn)%je1(m) == jeg(tn)) then ! corner is nc.
                   found_corner = .true.; tc = wCont(tn)%tile(m)
                   ie1 = sCont(tileMe)%is1(1) - 1; je1 = sCont(tileMe)%js1(1) - 1
                   ie2 = wCont(tn)%ie2(m);         je2 = wCont(tn)%je2(m)
                end if
             end if
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tc; align1Recv(numR) = SOUTH_WEST; align2Recv(numR) = NORTH_EAST
       is1Recv(numR) = ie1 - whalo + 1; ie1Recv(numR) = ie1
       js1Recv(numR) = je1 - shalo + 1; je1Recv(numR) = je1
       is2Recv(numR) = ie2 - whalo + 1; ie2Recv(numR) = ie2
       js2Recv(numR) = je2 - shalo + 1; je2Recv(numR) = je2 
    end if

    !--- northwest for recving
    found_corner = .false.
    n = wCont(tileMe)%ncontact
    if( n > 0) then
       if(wCont(tileMe)%je1(n) == jeg(tileMe) ) then
          tn = wCont(tileMe)%tile(n)
          if(wcont(tileMe)%je2(n) < jeg(tn) ) then  ! the corner tile is tn.
             if( jeg(tn) - wcont(tileMe)%je2(n) < nhalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: northwest tile for recv 1 is not tiled properly")
             found_corner = .true.; tc = tn
             ie1 = wCont(tileMe)%is1(n) - 1; js1 = wCont(tileMe)%je1(n) + 1
             ie2 = wCont(tileMe)%is2(n);     js2 = wCont(tileMe)%je2(n) + 1
          else if(nCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
             m = nCont(tn)%ncontact
             if(nCont(tn)%ie1(m) == ieg(tn)) then ! corner is nc.
                found_corner = .true.; tc = nCont(tn)%tile(m)
                ie1 = wCont(tileMe)%is1(n) - 1; js1 = wCont(tileMe)%je1(n) + 1
                ie2 = nCont(tn)%ie2(m);         js2 = nCont(tn)%js2(m)
             end if
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found,
       if( nCont(tileMe)%ncontact > 0) then
          if( nCont(tileMe)%is1(1) == isg(tileMe)) then
             tn = nCont(tileMe)%tile(1)
             if(nCont(tileMe)%is2(1) > isg(tn) ) then  ! the corner tile is tn.
                if( ncont(tileMe)%is2(1)-isg(tn) < whalo ) call mpp_error(FATAL, &
                     "mpp_domains_define.inc: northwest tile for recv 2 is not tiled properly")
                found_corner = .true.; tc = tn
                ie1 = nCont(tileMe)%is1(1) - 1; js1 = nCont(tileMe)%je1(1) + 1
                ie2 = nCont(tileMe)%is2(1) - 1; js2 = nCont(tileMe)%js2(1)
             else if(wCont(tn)%ncontact >0) then ! the corner tile may be west tile of tn.
                if(wCont(tn)%js1(1) == jsg(tn)) then ! corner is nc.
                   found_corner = .true.; tc = wCont(tn)%tile(1)
                   ie1 = nCont(tileMe)%is1(1) - 1; js1 = nCont(tileMe)%je1(1) + 1
                   ie2 = wCont(tn)%ie2(1);         js2 = wCont(tn)%js2(1)
                end if
             end if
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tc; align1Recv(numR) =NORTH_WEST;  align2Recv(numR) = SOUTH_EAST
       is1Recv(numR) = ie1 - whalo + 1; ie1Recv(numR) = ie1
       js1Recv(numR) = js1;             je1Recv(numR) = js1 + nhalo - 1
       is2Recv(numR) = ie2 - whalo + 1; ie2Recv(numR) = ie2
       js2Recv(numR) = js2;             je2Recv(numR) = js2 + nhalo - 1  
    end if

    !--- northeast for recving
    found_corner = .false.
    n = eCont(tileMe)%ncontact
    if( n > 0) then
       if(eCont(tileMe)%je1(n) == jeg(tileMe) ) then
          tn = eCont(tileMe)%tile(n)
          if(econt(tileMe)%je2(n) < jeg(tn) ) then  ! the corner tile is tn.
             if( jeg(tn) - econt(tileMe)%je2(n) < nhalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: northeast tile for recv 1 is not tiled properly")
             found_corner = .true.; tc = tn
             is1 = eCont(tileMe)%ie1(n) + 1; js1 = eCont(tileMe)%je1(n) + 1
             is2 = eCont(tileMe)%is2(1);     js2 = eCont(tileMe)%je2(1) + 1
          else if(nCont(tn)%ncontact >0) then ! the corner tile may be south tile of tn.
             if(nCont(tn)%is1(1) == isg(tn)) then ! corner is nc.
                found_corner = .true.; tc = nCont(tn)%tile(1)
                is1 = eCont(tileMe)%ie1(n) + 1; js1 = eCont(tileMe)%je1(n) + 1
                is2 = nCont(tn)%is2(1);         js2 = nCont(tn)%js2(1)
             end if
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found,
       n = nCont(tileMe)%ncontact
       if( n > 0) then
          if( nCont(tileMe)%ie1(n) == ieg(tileMe)) then
             tn = nCont(tileMe)%tile(n)
             if(nCont(tileMe)%ie2(n) < ieg(tn) ) then  ! the corner tile is tn.
                if(ieg(tn) - sCont(tileMe)%ie2(n) < ehalo ) call mpp_error(FATAL, &
                     "mpp_domains_define.inc: northeast tile for recv 2 is not tiled properly")
                found_corner = .true.; tc = tn
                is1 = sCont(tileMe)%ie1(n) + 1; js1 = sCont(tileMe)%je1(n) + 1
                is2 = sCont(tileMe)%ie2(n) + 1; js2 = sCont(tileMe)%js2(n)
             else if(eCont(tn)%ncontact >0) then ! the corner tile may be east tile of tn.
                if(eCont(tn)%js1(1) == jsg(tn)) then ! corner is nc.
                   found_corner = .true.; tc = eCont(tn)%tile(1)
                   is1 = sCont(tileMe)%ie1(n) + 1; js1 = sCont(tileMe)%je1(n) + 1
                   is2 = eCont(tn)%is2(m);         js2 = eCont(tn)%js2(m)
                end if
             end if
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tc; align1Recv(numR) =NORTH_EAST;  align2Recv(numR) = SOUTH_WEST
       is1Recv(numR) = is1; ie1Recv(numR) = is1 + ehalo - 1
       js1Recv(numR) = js1; je1Recv(numR) = js1 + nhalo - 1
       is2Recv(numR) = is2; ie2Recv(numR) = is2 + ehalo - 1
       js2Recv(numR) = js2; je2Recv(numR) = js2 + nhalo - 1  
    end if
    
    !--- to_pe's southeast for sending
    do n = 1, wCont(tileMe)%ncontact
       tn = wCont(tileMe)%tile(n)
       if(wCont(tileMe)%js2(n) == jsg(tn) ) then
          if(wcont(tileMe)%js1(n) > jsg(tileMe) ) then  ! send to  tile tn.
             if( wcont(tileMe)%js1(n) - jsg(tileMe) < shalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southeast tile for send 1 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
             is1Send(numS) = wCont(tileMe)%is1(n);     ie1Send(numS) = is1Send(numS) + ehalo - 1  
             je1Send(numS) = wCont(tileMe)%js1(n) - 1; js1Send(numS) = je1Send(numS) - shalo + 1
             is2Send(numS) = wCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
             je2Send(numS) = wCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
          end if
       end if
    end do
    do n = 1, nCont(tileMe)%ncontact
       tn = nCont(tileMe)%tile(n)
       if(nCont(tileMe)%ie2(n) == ieg(tn) ) then
          if(nCont(tileMe)%ie1(n) < ieg(tileMe) ) then  ! send to  tile tn.
             if( ieg(tileMe) - nCont(tileMe)%ie1(n) < ehalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southeast tile for send 2 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
             is1Send(numS) = nCont(tileMe)%ie1(n) + 1; ie1Send(numS) = is1Send(numS) + ehalo - 1  
             je1Send(numS) = nCont(tileMe)%je1(n)    ; js1Send(numS) = je1Send(numS) - shalo + 1
             is2Send(numS) = nCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
             je2Send(numS) = nCont(tileMe)%je2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1
          end if
       end if
    end do

    !--- found the corner overlap that is not specified through contact line.
    n = wCont(tileMe)%ncontact
    found_corner = .false.
    if( n > 0) then
       tn = wCont(tileMe)%tile(n)
       if( wCont(tileMe)%je1(n) == jeg(tileMe) .AND. wCont(tileMe)%je2(n) == jeg(tn) ) then
          m = nCont(tn)%ncontact
          if(m >0) then
             tc = nCont(tn)%tile(m)
             if( nCont(tn)%ie1(m) == ieg(tn) .AND. nCont(tn)%ie2(m) == ieg(tc) ) found_corner = .true.
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found, then starting from north contact
       if( nCont(tileMe)%ncontact > 0) then
          tn = nCont(tileMe)%tile(1)
          if( nCont(tileMe)%is1(1) == isg(tileMe) .AND. nCont(tileMe)%is2(1) == isg(tn) ) then
             if(wCont(tn)%ncontact >0) then 
                tc = wCont(tn)%tile(1)
                if( wCont(tn)%js1(1) == jsg(tn) .AND. wCont(tn)%js2(1) == jsg(tc) ) found_corner = .true.
             end if
          end if
       end if
    end if

    if(found_corner) then
       numS = numS+1; tileSend(numS) = tc
       align1Send(numS) = NORTH_WEST;  align2Send(numS) = SOUTH_EAST
       is1Send(numS) = isg(tileMe); ie1Send(numS) = is1Send(numS) + ehalo - 1  
       je1Send(numS) = jeg(tileMe); js1Send(numS) = je1Send(numS) - shalo + 1
       is2Send(numS) = ieg(tc) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
       je2Send(numS) = jsg(tc) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
    end if

    !--- to_pe's southwest for sending
    do n = 1, eCont(tileMe)%ncontact
       tn = eCont(tileMe)%tile(n)
       if(eCont(tileMe)%js2(n) == jsg(tn) ) then
          if(econt(tileMe)%js1(n) > jsg(tileMe) ) then  ! send to  tile tn.
             if( econt(tileMe)%js1(n) - jsg(tileMe) < shalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southwest tile for send 1 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
             ie1Send(numS) = eCont(tileMe)%ie1(n);     is1Send(numS) = ie1Send(numS) - whalo + 1  
             je1Send(numS) = eCont(tileMe)%js1(n) - 1; js1Send(numS) = je1Send(numS) - shalo + 1
             ie2Send(numS) = eCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
             je2Send(numS) = eCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
          end if
       end if
    end do
    do n = 1, nCont(tileMe)%ncontact
       tn = nCont(tileMe)%tile(n)
       if(nCont(tileMe)%is2(n) == isg(tn) ) then
          if(ncont(tileMe)%is1(n) > isg(tileMe) ) then  ! send to  tile tn.
             if( ncont(tileMe)%is1(n) - isg(tileMe) < whalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southwest tile for send 2 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
             ie1Send(numS) = nCont(tileMe)%is1(n) - 1; is1Send(numS) = ie1Send(numS) - whalo + 1  
             ie1Send(numS) = nCont(tileMe)%je1(n)    ; js1Send(numS) = je1Send(numS) - shalo + 1
             ie2Send(numS) = nCont(tileMe)%is2(n) - 1; is2Send(numS) = je2Send(numS) - whalo + 1
             je2Send(numS) = nCont(tileMe)%js2(n) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
          end if
       end if
    end do

    !--- found the corner overlap that is not specified through contact line.
    n = eCont(tileMe)%ncontact
    found_corner = .false.
    if( n > 0) then
       tn = eCont(tileMe)%tile(n)
       if( eCont(tileMe)%je1(n) == jeg(tileMe) .AND. eCont(tileMe)%je2(n) == jeg(tn) ) then
          if(nCont(tn)%ncontact >0) then
             tc = nCont(tn)%tile(1)
             if( nCont(tn)%is1(1) == isg(tn) .AND. nCont(tn)%is2(n) == isg(tc) ) found_corner = .true.
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found, then starting from north contact
       n = nCont(tileMe)%ncontact
       if( n > 0) then
          tn = nCont(tileMe)%tile(n)
          if( nCont(tileMe)%ie1(n) == ieg(tileMe) .AND. nCont(tileMe)%ie2(n) == ieg(tn) ) then
             if(eCont(tn)%ncontact >0) then
                tc = eCont(tn)%tile(1)
                if( eCont(tn)%js1(1) == jsg(tn) .AND. eCont(tn)%js2(n) == jsg(tc) ) found_corner = .true.
             end if
          end if
       end if
    end if

    if(found_corner) then
       numS = numS+1; tileSend(numS) = tc
       align1Send(numS) = NORTH_EAST;  align2Send(numS) = SOUTH_WEST
       ie1Send(numS) = ieg(tileMe); is1Send(numS) = ie1Send(numS) - whalo + 1  
       je1Send(numS) = jeg(tileMe); js1Send(numS) = je1Send(numS) - shalo + 1
       ie2Send(numS) = isg(tc) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
       je2Send(numS) = jsg(tc) - 1; js2Send(numS) = je2Send(numS) - shalo + 1  
    end if

    !--- to_pe's northwest for sending
    do n = 1, eCont(tileMe)%ncontact
       tn = eCont(tileMe)%tile(n)
       if(eCont(tileMe)%je2(n) == jeg(tn) ) then
          if(econt(tileMe)%je1(n) < jeg(tileMe) ) then  ! send to  tile tn.
             if( jeg(tileMe) - econt(tileMe)%je1(n) < nhalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: northwest tile for send 1 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
             ie1Send(numS) = eCont(tileMe)%ie1(n)    ; is1Send(numS) = ie1Send(numS) - whalo + 1  
             js1Send(numS) = eCont(tileMe)%je1(n) + 1; je1Send(numS) = js1Send(numS) + nhalo - 1
             ie2Send(numS) = eCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
             js2Send(numS) = eCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
          end if
       end if
    end do

    do n = 1, sCont(tileMe)%ncontact
       tn = sCont(tileMe)%tile(n)
       if(sCont(tileMe)%is2(n) == isg(tn) ) then
          if(scont(tileMe)%is1(n) > isg(tileMe) ) then  ! send to  tile tn.
             if( scont(tileMe)%is1(n) - isg(tileMe) < whalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southwest tile for send 2 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
             ie1Send(numS) = nCont(tileMe)%is1(n) - 1; is1Send(numS) = ie1Send(numS) - whalo + 1  
             js1Send(numS) = nCont(tileMe)%je1(n)    ; je1Send(numS) = js1Send(numS) + nhalo - 1
             ie2Send(numS) = nCont(tileMe)%is2(n) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
             js2Send(numS) = nCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
          end if
       end if
    end do

    !--- found the corner overlap that is not specified through contact line.
    n = eCont(tileMe)%ncontact
    found_corner = .false.
    if( n > 0) then
       tn = eCont(tileMe)%tile(1)
       if( eCont(tileMe)%js1(1) == jsg(tileMe) .AND. eCont(tileMe)%js2(1) == jsg(tn) ) then
          if(sCont(tn)%ncontact >0) then
             tc = sCont(tn)%tile(1)
             if( sCont(tn)%is1(1) == isg(tn) .AND. sCont(tn)%is2(1) == isg(tc) ) found_corner = .true.
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found, then starting from north contact
       n = sCont(tileMe)%ncontact
       found_corner = .false.
       if( n > 0) then
          tn = sCont(tileMe)%tile(n)
          if( sCont(tileMe)%ie1(n) == ieg(tileMe) .AND. sCont(tileMe)%ie2(n) == ieg(tn) ) then
             if(eCont(tn)%ncontact >0) then
                tc = eCont(tn)%tile(n)
                if( eCont(tn)%je1(n) == jeg(tn) .AND. eCont(tn)%je2(n) == jeg(tc) ) found_corner = .true.
             end if
          end if
       end if
    end if

    if(found_corner) then
       numS = numS+1; tileSend(numS) = tc
       align1Send(numS) = SOUTH_EAST;  align2Send(numS) = NORTH_WEST
       ie1Send(numS) = ieg(tileMe); is1Send(numS) = ie1Send(numS) - whalo + 1  
       js1Send(numS) = jsg(tileMe); je1Send(numS) = js1Send(numS) + nhalo - 1
       ie2Send(numS) = isg(tc) - 1; is2Send(numS) = ie2Send(numS) - whalo + 1
       js2Send(numS) = jeg(tc) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if

    !--- to_pe's northeast for sending
    do n = 1, wCont(tileMe)%ncontact
       tn = wCont(tileMe)%tile(n)
       if(wCont(tileMe)%je2(n) == jeg(tn) ) then
          if(wcont(tileMe)%je1(n) < jeg(tileMe) ) then  ! send to  tile tn.
             if( jeg(tileMe) - wcont(tileMe)%je1(n) < nhalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: northeast tile for send 1 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
             is1Send(numS) = wCont(tileMe)%is1(n)    ; ie1Send(numS) = is1Send(numS) + ehalo - 1  
             js1Send(numS) = wCont(tileMe)%je1(n) + 1; je1Send(numS) = js1Send(numS) + nhalo - 1
             is2Send(numS) = wCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
             js2Send(numS) = wCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
          end if
       end if
    end do

    do n = 1, sCont(tileMe)%ncontact
       tn = sCont(tileMe)%tile(n)
       if(sCont(tileMe)%ie2(n) == ieg(tn) ) then
          if(sCont(tileMe)%ie1(n) < ieg(tileMe) ) then  ! send to  tile tn.
             if( ieg(tileMe) - sCont(tileMe)%ie1(n) < ehalo ) call mpp_error(FATAL, &
                  "mpp_domains_define.inc: southeast tile for send 2 is not tiled properly")
             numS = numS+1; tileSend(numS) = tn
             align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
             is1Send(numS) = sCont(tileMe)%ie1(n) + 1; ie1Send(numS) = is1Send(numS) + ehalo - 1  
             js1Send(numS) = sCont(tileMe)%js1(n)    ; je1Send(numS) = js1Send(numS) + nhalo - 1
             is2Send(numS) = sCont(tileMe)%ie2(n) + 1; ie2Send(numS) = is1Send(numS) + ehalo - 1
             js2Send(numS) = sCont(tileMe)%je2(n) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
          end if
       end if
    end do

    !--- found the corner overlap that is not specified through contact line.
    n = wCont(tileMe)%ncontact
    found_corner = .false.
    if( n > 0) then
       tn = wCont(tileMe)%tile(1)
       if( wCont(tileMe)%js1(n) == jsg(tileMe) .AND. wCont(tileMe)%js2(n) == jsg(tn) ) then
          m = sCont(tn)%ncontact
          if(m >0) then
             tc = sCont(tn)%tile(m)
             if( sCont(tn)%ie1(m) == ieg(tn) .AND. sCont(tn)%ie2(m) == ieg(tc) ) found_corner = .true.
          end if
       end if
    end if
    if( .not. found_corner ) then  ! not found, then starting from north contact
       n = sCont(tileMe)%ncontact
       found_corner = .false.
       if( n > 0) then
          tn = sCont(tileMe)%tile(1)
          if( sCont(tileMe)%is1(1) == isg(tileMe) .AND. sCont(tileMe)%is2(1) == isg(tn) ) then
             m = wCont(tn)%ncontact
             if( m > 0 ) then
                tc = wCont(tn)%tile(m)
                if( wCont(tn)%je1(m) == jeg(tn) .AND. wCont(tn)%je2(m) == jeg(tc) ) found_corner = .true.
             end if
          end if
       end if
    end if
    if(found_corner) then
       numS = numS+1; tileSend(numS) = tc
       align1Send(numS) = SOUTH_WEST;  align2Send(numS) = NORTH_EAST
       is1Send(numS) = isg(tileMe); ie1Send(numS) = is1Send(numS) + ehalo - 1  
       js1Send(numS) = jsg(tileMe); je1Send(numS) = js1Send(numS) + nhalo - 1
       is2Send(numS) = ieg(tc) + 1; ie2Send(numS) = is2Send(numS) + ehalo - 1
       js2Send(numS) = jeg(tc) + 1; je2Send(numS) = js2Send(numS) + nhalo - 1  
    end if

  end subroutine fill_corner_contact

     !--- find the alignment direction, check if index is reversed, if reversed, exchange index.
  subroutine check_alignment( is, ie, js, je, isg, ieg, jsg, jeg, alignment )
    integer, intent(inout) :: is, ie, js, je, isg, ieg, jsg, jeg
    integer, intent(out)   :: alignment

    integer :: i, j

    if ( is == ie ) then      ! x-alignment
       if ( is == isg ) then
          alignment = WEST
       else if ( is == ieg ) then
          alignment = EAST
       else
          call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region is not on the x-boundary of the tile')
       end if
       if ( js > je ) then
          j = js; js = je; je = j
       end if
    else if ( js == je ) then ! y-alignment
       if ( js == jsg ) then
          alignment = SOUTH
       else if ( js == jeg ) then
          alignment = NORTH
       else
          call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region is not on the y-boundary of the tile')
       end if
       if ( is > ie ) then
          i = is; is = ie; ie = i
       end if
    else
       call mpp_error(FATAL, 'mpp_domains_define.inc: The contact region should be line contact' )
    end if

  end subroutine check_alignment
  !#####################################################################  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_MODIFY_DOMAIN: modify extent of domain                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_modify_domain1D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain1D)" > </IN>
!   <IN NAME="hbegin,hend" TYPE="integer,optional" > </IN>
!   <IN NAME="cbegin,cend" TYPE="integer,optional" > </IN>
!   <IN NAME="gbegin,gend" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain1D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain1D(domain_in,domain_out,cbegin,cend,gbegin,gend, hbegin, hend)
      ! </PUBLICROUTINE>
      type(domain1D), intent(in)    :: domain_in
      type(domain1D), intent(inout) :: domain_out
      integer, intent(in), optional :: hbegin, hend             ! halo size 
      integer, intent(in), optional :: cbegin, cend             ! extent of compute_domain
      integer, intent(in), optional :: gbegin, gend             ! extent of global domain
      integer :: ndivs, global_indices(2) !(/ isg, ieg /)
      integer                       :: flag
! get the global indices of the input domain
      global_indices(1) = domain_in%global%begin;  global_indices(2) = domain_in%global%end

      ! get the layout
      ndivs = size(domain_in%list(:))

! get the flag
      flag = 0
      if(domain_in%cyclic) flag = flag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%data%is_global) flag = flag + GLOBAL_DATA_DOMAIN

      call mpp_define_domains( global_indices, ndivs, domain_out, pelist = domain_in%list(:)%pe, &
           flags = flag, begin_halo = hbegin, end_halo = hend, extent = domain_in%list(:)%compute%size )
           
      if(present(cbegin)) domain_out%compute%begin = cbegin
      if(present(cend))   domain_out%compute%end = cend
      domain_out%compute%size = domain_out%compute%end - domain_out%compute%begin + 1
      if(present(gbegin)) domain_out%global%begin = gbegin
      if(present(gend))   domain_out%global%end = gend
      domain_out%global%size = domain_out%global%end - domain_out%global%begin + 1
      
    end subroutine mpp_modify_domain1D
! </SUBROUTINE>

  !#######################################################################
!----------------------------------------------------------------------------------
! <SUBROUTINE NAME="mpp_modify_domain2D" INTERFACE="mpp_modify_domain">
!   <IN NAME="domain_in" TYPE="type(domain2D)" > </IN>
!   <IN NAME="isc,iec" TYPE="integer,optional" > </IN>
!   <IN NAME="jsc,jec" TYPE="integer,optional" > </IN>
!   <IN NAME="isg,ieg" TYPE="integer,optional" > </IN>
!   <IN NAME="jsg,jeg" TYPE="integer,optional" > </IN>
!   <IN NAME="whalo,ehalo" TYPE="integer,optional" > </IN>
!   <IN NAME="shalo,nhalo" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain2D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain2D(domain_in, domain_out, isc, iec, jsc, jec, isg, ieg, jsg, jeg, whalo, ehalo, shalo, nhalo)
      ! </PUBLICROUTINE>
      type(domain2D), intent(in)    :: domain_in
      type(domain2D), intent(inout) :: domain_out
      integer, intent(in), optional :: isc, iec, jsc, jec
      integer, intent(in), optional :: isg, ieg, jsg, jeg
      integer, intent(in), optional :: whalo, ehalo, shalo, nhalo
      integer                       :: global_indices(4), layout(2)
      integer                       :: xflag, yflag
      
      if(present(whalo) .or. present(ehalo) .or. present(shalo) .or. present(nhalo) ) then
         ! get the global indices of the input domain
         global_indices(1) = domain_in%x(1)%global%begin;  global_indices(2) = domain_in%x(1)%global%end
         global_indices(3) = domain_in%y(1)%global%begin;  global_indices(4) = domain_in%y(1)%global%end

         ! get the layout
         layout(1) = size(domain_in%x(1)%list(:)); layout(2) = size(domain_in%y(1)%list(:))

         ! get the flag
         xflag = 0; yflag = 0
         if(domain_in%x(1)%cyclic) xflag = xflag + CYCLIC_GLOBAL_DOMAIN
         if(domain_in%x(1)%data%is_global) xflag = xflag + GLOBAL_DATA_DOMAIN
         if(domain_in%y(1)%cyclic) yflag = yflag + CYCLIC_GLOBAL_DOMAIN
         if(domain_in%y(1)%data%is_global) yflag = yflag + GLOBAL_DATA_DOMAIN

         call mpp_define_domains( global_indices, layout, domain_out, pelist = domain_in%list(:)%pe, &
                                  xflags = xflag, yflags = yflag,  whalo = whalo, ehalo = ehalo,     &
                                  shalo = shalo, nhalo = nhalo,                                      &
                                  xextent = domain_in%x(1)%list(:)%compute%size,                     &
                                  yextent = domain_in%y(1)%list(:)%compute%size,                     &
                                  symmetry=domain_in%symmetry,                                       &
                                  maskmap = domain_in%pearray .NE. NULL_PE )
                               
      else    
         call mpp_define_null_domain(domain_out)
         call mpp_modify_domain(domain_in%x(1), domain_out%x(1), isc, iec, isg, ieg)
         call mpp_modify_domain(domain_in%y(1), domain_out%y(1), jsc, jec, jsg, jeg)
      endif
         
    end subroutine mpp_modify_domain2D
! </SUBROUTINE>

  !#####################################################################


   subroutine mpp_define_null_domain1D(domain)
      type(domain1D), intent(inout) :: domain

      domain%global%begin  = -1; domain%global%end  = -1; domain%global%size = 0
      domain%data%begin    = -1; domain%data%end    = -1; domain%data%size = 0
      domain%compute%begin = -1; domain%compute%end = -1; domain%compute%size = 0
      domain%pe = NULL_PE

   end subroutine mpp_define_null_domain1D

  !#####################################################################


   subroutine mpp_define_null_domain2D(domain)
     type(domain2D), intent(inout) :: domain

      allocate(domain%x(1), domain%y(1))
      call mpp_define_null_domain(domain%x(1))
      call mpp_define_null_domain(domain%y(1))
      domain%pe = NULL_PE

   end subroutine mpp_define_null_domain2D

  !#####################################################################
