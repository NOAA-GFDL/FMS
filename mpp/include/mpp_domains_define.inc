
! -*-f90-*-
! $Id: mpp_domains_define.inc,v 13.0 2006/03/28 21:40:36 fms Exp $


  ! <SUBROUTINE NAME="mpp_define_layout2D" INTERFACE="mpp_define_layout">
  !  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
  !  <IN NAME="ndivs" TYPE="integer"></IN>
  !  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
  ! </SUBROUTINE>
  subroutine mpp_define_layout2D( global_indices, ndivs, layout )
    integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: ndivs !number of divisions to divide global domain
    integer, intent(out) :: layout(2)

    integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

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
  subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap )
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
    integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
    integer, intent(in) :: ndivs
    type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: flags, halo
    integer, intent(in), optional :: extent(0:)
    logical, intent(in), optional :: maskmap(0:)

    logical :: compute_domain_is_global, data_domain_is_global
    integer :: ndiv, n, isg, ieg, is, ie, i
    integer, allocatable :: pes(:)
    logical, allocatable :: mask(:)
    integer :: halosz
    !used by symmetry algorithm
    integer :: imax, ndmax, ndmirror
    logical :: symmetrize
    !statement functions
    logical :: even, odd
    even(n) = (mod(n,2).EQ.0)
    odd (n) = (mod(n,2).EQ.1)

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
    !get global indices
    isg = global_indices(1)
    ieg = global_indices(2)
    if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
    !get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
    if( PRESENT(pelist) )then
       if( .NOT.any(pelist.EQ.mpp_pe()) )then
          write( stderr(),* )'pe=', mpp_pe(), ' pelist=', pelist
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
    if( PRESENT(extent) )then
       if( size(extent(:)).NE.ndivs ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
    end if
    !get halosize
    halosz = 0
    if( PRESENT(halo) )halosz = halo

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
          if( PRESENT(extent) )then
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
    domain%list(:)%data%begin = domain%list(:)%data%begin - halosz
    domain%list(:)%data%end   = domain%list(:)%data%end   + halosz  
    domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

    domain%compute = domain%list(domain%pos)%compute
    domain%data = domain%list(domain%pos)%data
    domain%global = domain%list(domain%pos)%global
    domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
    domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
    domain%global%max_size  = domain%global%size

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
  subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, xflags, yflags, &
         xhalo, yhalo, xextent, yextent, maskmap, name, symmetry, tile_id, ntiles, topology_type )
    !define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
    integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
    integer, intent(in) :: layout(2)
    type(domain2D), intent(inout) :: domain
    integer, intent(in), optional :: pelist(0:)
    integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
    integer, intent(in), optional :: xextent(0:), yextent(0:)
    logical, intent(in), optional :: maskmap(0:,0:)
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: symmetry
    integer, intent(in), optional :: tile_id  ! id of current tile
    integer, intent(in), optional :: ntiles   ! total number of tiles for the mosaic.
    integer, intent(in), optional :: topology_type
    integer :: i, j, m, n
    integer :: ipos, jpos, pos, tid, nt
    integer :: ndivx, ndivy, isg, ieg, jsg, jeg, ishift, jshift
    integer(LONG_KIND),save :: domain_cnt=0

    logical, allocatable :: mask(:,:)
    integer, allocatable :: pes(:), pesall(:)
    character(len=8) :: text

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
    ndivx = layout(1); ndivy = layout(2)
    isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

    tid = 1    ! tile number, this will be needed when computing the contact region overlapping.
    if( PRESENT(tile_id)) tid = tile_id

    nt = 1
    if( PRESENT(ntiles) ) nt = ntiles

    if( PRESENT(pelist) )then
       allocate( pes(0:size(pelist(:))-1) )
       pes = pelist
       if(nt == 1) then
          allocate( pesall(0:size(pes(:))-1) )
          pesall = pes
       else
          allocate( pesall(0:mpp_npes()-1) )
          call mpp_get_current_pelist(pesall)   
       end if
    else
       allocate( pes(0:mpp_npes()-1) )
       allocate( pesall(0:mpp_npes()-1) )
       call mpp_get_current_pelist(pes)
       pesall = pes
       !          pes = (/ (i,i=0,mpp_npes()-1) /)
    end if

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
       allocate( domain%pearray(0:ndivx-1,0:ndivy-1) )
       domain%pearray(:,:) = NULL_PE
       ipos = NULL_PE; jpos = NULL_PE
       n = 0
       do j = 0,ndivy-1
          do i = 0,ndivx-1
             if( mask(i,j) )then
                domain%pearray(i,j) = pes(n)
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
       if( debug )write( stderr(), * )'pe, ipos, jpos=', mpp_pe(), ipos, jpos, ' pearray(:,jpos)=', &
                  domain%pearray(:,jpos), ' pearray(ipos,:)=', domain%pearray(ipos,:)

       !do domain decomposition using 1D versions in X and Y, 
       call mpp_define_domains( global_indices(1:2), ndivx, domain%x, &
            pack(domain%pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos) )
       call mpp_define_domains( global_indices(3:4), ndivy, domain%y, &
            pack(domain%pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:) )
       if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
            call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
       deallocate( mask )

       !--- set 1d domain for T,E,N,C-cell, 
       allocate(domain%T, domain%E, domain%N, domain%C)
       call set_domain1d(domain%x, domain%T%x, 0)
       call set_domain1d(domain%y, domain%T%y, 0)
       call set_domain1d(domain%x, domain%E%x, ishift)
       call set_domain1d(domain%y, domain%E%y, 0)
       call set_domain1d(domain%x, domain%N%x, 0)
       call set_domain1d(domain%y, domain%N%y, jshift)
       call set_domain1d(domain%x, domain%C%x, ishift)
       call set_domain1d(domain%y, domain%C%y, jshift)

       domain%pos = pos
       domain%pe  = mpp_pe()
       domain%tile_id = tid
       domain%ntiles  = nt
       domain%ncontacts = 0
       domain%topology_type = REGULAR
       if(present(topology_type) ) domain%topology_type = topology_type

       domain_cnt = domain_cnt + INT(1,KIND=LONG_KIND)
       domain%id = domain_cnt*DOMAIN_ID_BASE  ! Must be LONG_KIND arithmetic

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
          if( domain%y%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%x%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in X must be even when there is a fold in Y.' )
          !check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x%list(i)%compute%size.NE.domain%x%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries ' // &
                                         'must line up (mirror-symmetric extents).' )
          end do
       end if
       if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
          if( domain%x%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%y%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in Y must be even when there is a fold in X.' )
          !check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y%list(i)%compute%size.NE.domain%y%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
       end if

       !set up domain%list
       if( debug )write( stderr(),'(a,9i4)' )'pe, domain=', mpp_pe(), domain_info_buf(1:8)
       if( mpp_pe().EQ.pes(0) .AND. PRESENT(name) )then
          write( stdlog(), '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( stdlog(), '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
       end if

       call mpp_get_compute_domain( domain, domain_info_buf(1), domain_info_buf(2), domain_info_buf(3), domain_info_buf(4) )
       call mpp_get_data_domain   ( domain, domain_info_buf(5), domain_info_buf(6), domain_info_buf(7), domain_info_buf(8) )
       call mpp_get_global_domain ( domain, domain_info_buf(9), domain_info_buf(10), domain_info_buf(11), domain_info_buf(12) )
       domain_info_buf(13) = tid
    end if   !     if( ANY(pes == mpp_pe()) ) 

    call mpp_sync()
    n = size(pesall(:))

    !the following are only used for storage of remote compute and data domain info
    if(.NOT. associated(domain%list) )allocate( domain%list(0:n-1) ) 

    ! Because of mosaic, we need to replace mpp_transmit with mpp_send and mpp_recv.

    !--- send the data from the processor in pes ( pelist of current tile ) to all other pe.
    if(ANY(pes == mpp_pe())) then
       do i = 0,n-1
          call mpp_send(domain_info_buf(1), plen=13, to_pe=pesall(mod(pos+n-i,n)) )
       end do
    end if

    !--- receive data from the processor in pes.
    do i = 0,n-1
       m = mod(pos+i,n)
       if(ANY(pes == pesall(m)) ) then
          domain%list(m)%pe = pesall(m)
             call mpp_recv(domain_info_buf(14), glen=13, from_pe=pesall(m) ) 
             domain%list(m)%x%compute%begin = domain_info_buf(14)
             domain%list(m)%x%compute%end   = domain_info_buf(15)
             domain%list(m)%y%compute%begin = domain_info_buf(16)
             domain%list(m)%y%compute%end   = domain_info_buf(17)
             domain%list(m)%x%data%begin    = domain_info_buf(18)
             domain%list(m)%x%data%end      = domain_info_buf(19)
             domain%list(m)%y%data%begin    = domain_info_buf(20)
             domain%list(m)%y%data%end      = domain_info_buf(21) 
             domain%list(m)%x%global%begin  = domain_info_buf(22)
             domain%list(m)%x%global%end    = domain_info_buf(23)
             domain%list(m)%y%global%begin  = domain_info_buf(24)
             domain%list(m)%y%global%end    = domain_info_buf(25)  
             domain%list(m)%tile_id         = domain_info_buf(26) 
             if( mpp_pe() .EQ. pes(0) .AND. PRESENT(name) )write( stdlog(), '(2x,i3,x,4i5,3x,4i5)' ) &
                                                     pesall(m), domain_info_buf(14:)
             domain%list(m)%x%compute%size = domain%list(m)%x%compute%end - domain%list(m)%x%compute%begin + 1
             domain%list(m)%y%compute%size = domain%list(m)%y%compute%end - domain%list(m)%y%compute%begin + 1
             domain%list(m)%x%data%size    = domain%list(m)%x%data%end    - domain%list(m)%x%data%begin    + 1
             domain%list(m)%y%data%size    = domain%list(m)%y%data%end    - domain%list(m)%y%data%begin    + 1  
             domain%list(m)%x%global%size  = domain%list(m)%x%global%end  - domain%list(m)%x%global%begin  + 1
             domain%list(m)%y%global%size  = domain%list(m)%y%global%end  - domain%list(m)%y%global%begin  + 1  
       end if
    end do

    call mpp_sync_self( pesall )

    !--- set up 2-D domain decomposition for T, E, C, N and computing overlapping
    !--- when current tile is the last tile in the mosaic. 
    if(tid == nt) then
       call set_domain2d(domain, domain%T, CENTER, 0,      0     )
       call set_domain2d(domain, domain%E, EAST,   ishift, 0     )
       call set_domain2d(domain, domain%N, NORTH,  0,      jshift)
       call set_domain2d(domain, domain%C, CORNER, ishift, jshift)
       call compute_overlaps(domain%T, 0,      0)
       call compute_overlaps(domain%E, ishift, 0      )
       call compute_overlaps(domain%N, 0,      jshift )
       call compute_overlaps(domain%C, ishift, jshift )
    end if

    !PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)

    !print out decomposition, this didn't consider maskmap.
    if( mpp_pe() .EQ. pes(0) .AND. PRESENT(name) )then
       write(*,*) trim(name)//' domain decomposition'
       write (*,110) (domain%x%list(i)%compute%size, i= 0, layout(1)-1)
       write (*,120) (domain%y%list(i)%compute%size, i= 0, layout(2)-1)
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
  subroutine set_domain1d(domain_in, domain_out, shift)
    type(domain1d), intent(in)    :: domain_in
    type(domain1d), intent(inout) :: domain_out
    integer,        intent(in)    :: shift

    integer :: n, ndivs
  
    domain_out%cyclic = domain_in%cyclic
    domain_out%pos    = domain_in%pos

    ndivs = size(domain_in%list(:))
    allocate( domain_out%list(0:ndivs-1) )

    do n = 0, ndivs-1
       domain_out%list(n)%pe                = domain_in%list(n)%pe
       domain_out%list(n)%global %begin     = domain_in%list(n)%global %begin
       domain_out%list(n)%global %end       = domain_in%list(n)%global %end      + shift
       domain_out%list(n)%global %size      = domain_in%list(n)%global %size     + shift
       domain_out%list(n)%global %is_global = domain_in%list(n)%global %is_global
       domain_out%list(n)%compute%begin     = domain_in%list(n)%compute%begin
       domain_out%list(n)%compute%end       = domain_in%list(n)%compute%end      + shift
       domain_out%list(n)%compute%size      = domain_in%list(n)%compute%size     + shift
       domain_out%list(n)%compute%is_global = domain_in%list(n)%compute%is_global
       domain_out%list(n)%data   %begin     = domain_in%list(n)%data   %begin
       domain_out%list(n)%data   %end       = domain_in%list(n)%data   %end      + shift
       domain_out%list(n)%data   %size      = domain_in%list(n)%data   %size     + shift
       domain_out%list(n)%data   %is_global = domain_in%list(n)%data   %is_global
   enddo

    domain_out%compute          = domain_out%list(domain_out%pos)%compute
    domain_out%data             = domain_out%list(domain_out%pos)%data
    domain_out%global           = domain_out%list(domain_out%pos)%global
    domain_out%compute%max_size = MAXVAL( domain_out%list(:)%compute%size )
    domain_out%data%max_size    = MAXVAL( domain_out%list(:)%data%size )
    domain_out%global%max_size  = domain_out%global%size

    return

  end subroutine set_domain1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                    !
  !            set domain1d: set the 1-D domain decomposition          !
  !                                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_domain2d(domain_in, domain_out, position, ishift, jshift)
    type(domain2d), intent(in)    :: domain_in
    type(domain2d), intent(inout) :: domain_out
    integer,        intent(in)    :: position    
    integer,        intent(in)    :: ishift, jshift
    integer                       :: ndivs, n

    domain_out%pos         = domain_in%pos
    domain_out%pe          = domain_in%pe
    domain_out%tile_id     = domain_in%tile_id
    domain_out%ntiles      = domain_in%ntiles
    domain_out%topology_type = domain_in%topology_type

    ndivs = size(domain_in%list(:))
    allocate( domain_out%list(0:ndivs-1) )

    domain_out%position = position
    domain_out%symmetry = domain_in%symmetry
    domain_out%fold     = domain_in%fold

    do n = 0, ndivs-1
       domain_out%list(n)%pe               = domain_in%list(n)%pe
       domain_out%list(n)%x%compute%begin  = domain_in%list(n)%x%compute%begin
       domain_out%list(n)%x%compute%end    = domain_in%list(n)%x%compute%end   + ishift
       domain_out%list(n)%x%compute%size   = domain_in%list(n)%x%compute%size  + ishift
       domain_out%list(n)%y%compute%begin  = domain_in%list(n)%y%compute%begin
       domain_out%list(n)%y%compute%end    = domain_in%list(n)%y%compute%end   + jshift
       domain_out%list(n)%y%compute%size   = domain_in%list(n)%y%compute%size  + jshift
       domain_out%list(n)%x%data%begin     = domain_in%list(n)%x%data%begin
       domain_out%list(n)%x%data%end       = domain_in%list(n)%x%data%end      + ishift
       domain_out%list(n)%x%data%size      = domain_in%list(n)%x%data%size     + ishift
       domain_out%list(n)%y%data%begin     = domain_in%list(n)%y%data%begin
       domain_out%list(n)%y%data%end       = domain_in%list(n)%y%data%end      + jshift
       domain_out%list(n)%y%data%size      = domain_in%list(n)%y%data%size     + jshift
       domain_out%list(n)%x%global%begin   = domain_in%list(n)%x%global%begin
       domain_out%list(n)%x%global%end     = domain_in%list(n)%x%global%end      + ishift
       domain_out%list(n)%x%global%size    = domain_in%list(n)%x%global%size     + ishift
       domain_out%list(n)%y%global%begin   = domain_in%list(n)%y%global%begin
       domain_out%list(n)%y%global%end     = domain_in%list(n)%y%global%end      + jshift
       domain_out%list(n)%y%global%size    = domain_in%list(n)%y%global%size     + jshift
       domain_out%list(n)%tile_id          = domain_in%list(n)%tile_id
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
  subroutine mpp_define_mosaic( global_indices, layout, domain, num_tile, num_contact, tile1, tile2,       &
                                istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2, npes_tile, &
                                topology_type, pelist, xhalo, yhalo, xextent, yextent, maskmap, name         )
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
    integer,          intent(in)           :: npes_tile(:)         ! number of processor used in each tile
    character(len=*), intent(in), optional :: topology_type ! type of mosaic, value "cubic-grid" or "regular"(default)
    integer,          intent(in), optional :: pelist(:)            ! list of processors used in mosaic
    integer,          intent(in), optional :: xhalo, yhalo
    integer,          intent(in), optional :: xextent(:,:), yextent(:,:) 
    logical,          intent(in), optional :: maskmap(:,:,:)         
    character(len=*), intent(in), optional :: name

    integer              :: start, n, ndivx, ndivy, list, mType, nc, tileMe
    logical, allocatable :: mask(:,:)
    integer, allocatable :: pes(:)
    logical              :: is_symmetry
    integer, allocatable :: align1(:), align2(:), is1(:), ie1(:), js1(:), je1(:), is2(:), ie2(:), js2(:), je2(:)
    integer, allocatable :: isgList(:), iegList(:), jsgList(:), jegList(:), t1(:), t2(:)

    !--- the size of first indice of global_indices must be 4.
    if(size(global_indices, 1) .NE. 4) call mpp_error(FATAL, &
         'mpp_domains_define: The size of first dimension of global_indices is not 4')
    !--- the size of second indice of global_indices must be num_tile
    if(size(global_indices, 2) .NE. num_tile) call mpp_error(FATAL, &
         'mpp_domains_define: The size of second dimension of global_indices is not equal num_tile')    
    !--- the size of first indice of layout must be 2. The second dimension size of layout must equal num_tile.
    if(size(layout, 1) .NE. 2) call mpp_error(FATAL, &
         'mpp_domains_define: The size of first dimension of layout is not 2')
    if(size(layout,2) .NE. num_tile)  call mpp_error(FATAL, &
         'mpp_domains_define: The size of second dimension of layout is not equal num_tile')

    !--- setup pelist for the mosaic ---------------------
    allocate(pes(0:mpp_npes()-1))
    if(present(pelist)) then
       if(mpp_npes() .NE. size(pelist(:))) call mpp_error(FATAL, &
            'mpp_domains_define: size of pelist is not equal mpp_npes')
       pes = pelist
    else
       call mpp_get_current_pelist(pes)
    end if

    mType = REGULAR
    is_symmetry        = .FALSE.
    if( present(topology_type) ) then
       select case ( trim(topology_type) )
       case ( "regular" )
          ! do nothing
       case ( "cubic-grid" )
          is_symmetry = .TRUE.  ! cubic-grid is assume to be symmetry
          mType = CUBIC_GRID
       case default
          call mpp_error(FATAL, "mpp_domains_define: the value of optional arugment topology_type " // &
                                "should be regular of cubic-grid")
       end select
    end if

    !--- there should be 6 tiles in the mosaic.
    !--- there should be 12 contact lines in the mosaic
    if( mType == CUBIC_GRID ) then
       if( num_tile .NE. 6 ) call mpp_error(FATAL, 'mpp_domains_define: cubic-grid mosaic must have 6 tiles')
       if( num_contact .NE. 12) call mpp_error(FATAL, 'mpp_domains_define: cubic-grid mosaic must have 12 contact regions')
    end if

    if(size(npes_tile(:)) .NE. num_tile) call mpp_error(FATAL, &
         'mpp_domains_define: size of npeslist is not equal num_tile')      
    if(mpp_npes() .NE. sum(npes_tile)) call mpp_error(FATAL, &
         'mpp_domains_define: sum of npes_tile is not equal to mpp_npes')

    !--- check the size comformable of xextent and yextent
    if( PRESENT(xextent) ) then
       if(size(xextent,1) .GT. maxval(layout(1,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define: size mismatch between xextent and layout') 
       if(size(xextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define: size of xextent is not eqaul num_tile')    
    end if
    if( PRESENT(yextent) ) then
       if(size(yextent,1) .GT. maxval(layout(2,:)) ) call mpp_error(FATAL, &
            'mpp_domains_define: size mismatch between yextent and layout') 
       if(size(yextent,2) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define: size of yextent is not eqaul num_tile')    
    end if

    !--- check the size comformable of maskmap
    !--- since the layout is different between tiles, so the actual size of maskmap for each tile is 
    !--- not diffrent. When define maskmap for multiple tiles, user can choose the maximum value
    !--- of layout of all tiles to the first and second dimension of maskmap. 
    if(present(maskmap)) then
       if(size(maskmap,1) .GT. maxval(layout(1,:)) .or. size(maskmap,2) .GT. maxval(layout(2,:))) &
            call mpp_error(FATAL, 'mpp_domains_define: size mismatch between maskmap and layout')  
       if(size(maskmap,3) .NE. num_tile) call mpp_error(FATAL, &
            'mpp_domains_define: the third dimension of maskmap is not equal num_tile')  
    end if

    !---call mpp_define_domain to define domain decomposition for each tile.
    start = 0
    do n = 1, num_tile
       allocate(mask(layout(1,n), layout(2,n)))
       mask = .TRUE.
       if(present(maskmap))  mask = maskmap(1:layout(1,n), 1:layout(2,n), n)
       ndivx = layout(1,n); ndivy = layout(2,n)

       if(PRESENT(xextent) .AND. PRESENT(yextent)) then
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pes(start:start+npes_tile(n)-1), &
                                  xhalo=xhalo, yhalo=yhalo, xextent=xextent(1:ndivx,n),                             &
                                  yextent=yextent(1:ndivy,n), maskmap=mask, name=name, symmetry=is_symmetry,        &
                                  tile_id = n, ntiles = num_tile, topology_type = mType)
       else if (PRESENT(xextent)) then
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pes(start:start+npes_tile(n)-1), &
                                  xhalo=xhalo, yhalo=yhalo, xextent=xextent(1:ndivx,n), maskmap=mask, name=name,    &
                                  symmetry=is_symmetry, tile_id = n, ntiles = num_tile, topology_type = mType )
       else if (PRESENT(xextent)) then
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pes(start:start+npes_tile(n)-1), &
                                  xhalo=xhalo, yhalo=yhalo, yextent=yextent(1:ndivy,n), maskmap=mask, name=name,    &
                                  symmetry=is_symmetry, tile_id = n, ntiles = num_tile, topology_type = mType )
       else
          call mpp_define_domains(global_indices(:,n), layout(:,n), domain, pelist=pes(start:start+npes_tile(n)-1), &
                                  xhalo=xhalo, yhalo=yhalo, maskmap=mask, name=name, symmetry=is_symmetry,          &
                                  tile_id = n, ntiles = num_tile, topology_type = mType)
       end if
       start = start + npes_tile(n)
       deallocate(mask)
    end do

    deallocate(pes)

    domain%ncontacts = num_contact
    !--- loop through each contact region and find the contact for each tile ( including alignment )
    !--- we assume the tiles list is continuous and starting from 1.
    allocate(is1(2*num_contact), ie1(2*num_contact), js1(2*num_contact), je1(2*num_contact) )
    allocate(is2(2*num_contact), ie2(2*num_contact), js2(2*num_contact), je2(2*num_contact) )
    allocate(isgList(num_tile), iegList(num_tile), jsgList(num_tile), jegList(num_tile) )
    allocate(align1(2*num_contact), align2(2*num_contact), t1(2*num_contact), t2(2*num_contact) )
    !--- get the global domain for each tile
    do n = 1, num_tile
       do list = 0, size(domain%list(:))-1
          if( domain%list(list)%tile_id == n ) then 
             isgList(n) = domain%list(list)%x%global%begin;  iegList(n) = domain%list(list)%x%global%end
             jsgList(n) = domain%list(list)%y%global%begin;  jegList(n) = domain%list(list)%y%global%end  
             exit
          end if
       end do
    end do

    tileMe = domain%tile_id
    !--- transfer the contact index to domain index.
    nc = 0
    do n = 1, num_contact
       nc = nc + 1
       if(tile2(n) == tileMe) then   !--tileMe will be always the first
          t1(nc)  = tile2(n); t2(nc)  = tile1(n)
          is1(nc) = istart2(n) + isgList(tile2(n)) - 1; ie1(nc) = iend2(n) + isgList(tile2(n)) - 1
          js1(nc) = jstart2(n) + jsgList(tile2(n)) - 1; je1(nc) = jend2(n) + jsgList(tile2(n)) - 1
          is2(nc) = istart1(n) + isgList(tile1(n)) - 1; ie2(nc) = iend1(n) + isgList(tile1(n)) - 1
          js2(nc) = jstart1(n) + jsgList(tile1(n)) - 1; je2(nc) = jend1(n) + jsgList(tile1(n)) - 1
       else 
          t1(nc)  = tile1(n); t2(nc)  = tile2(n)
          is1(nc) = istart1(n) + isgList(tile1(n)) - 1; ie1(nc) = iend1(n) + isgList(tile1(n)) - 1
          js1(nc) = jstart1(n) + jsgList(tile1(n)) - 1; je1(nc) = jend1(n) + jsgList(tile1(n)) - 1
          is2(nc) = istart2(n) + isgList(tile2(n)) - 1; ie2(nc) = iend2(n) + isgList(tile2(n)) - 1
          js2(nc) = jstart2(n) + jsgList(tile2(n)) - 1; je2(nc) = jend2(n) + jsgList(tile2(n)) - 1
       end if
       call check_alignment( is1(nc), ie1(nc), js1(nc), je1(nc), isgList(tile1(n)), iegList(tile1(n)), &
                             jsgList(tile1(n)), jegList(tile1(n)), align1(nc) )
       call check_alignment( is2(nc), ie2(nc), js2(nc), je2(nc), isgList(tile2(n)), iegList(tile2(n)), &
                             jsgList(tile2(n)), jegList(tile2(n)), align2(nc) )

       !--- when the mosaic is not cubic-grid, we suppose there is no rotation between tiles.
       if( mType .NE. CUBIC_GRID ) then
          if ( (align1(nc) == EAST .OR. align1(nc) == WEST) .NEQV. (align2(nc) == EAST .OR. align2(nc) == WEST) ) then
              call mpp_error(FATAL,"mpp_domains_define_mod: No rotation between tiles is allowed " // &
                                   "for non cubic-grid, contact developer" )
          end if
       end if
       
       !--- when the contact is between one same tile, ont more contact is needed to be added
       if(tile1(n) == tile2(n)) then
          nc = nc + 1
          t1(nc)  = t2(nc-1);  t2(nc)  = t1(nc-1)
          is1(nc) = is2(nc-1); ie1(nc) = ie2(nc-1)
          js1(nc) = js2(nc-1); je1(nc) = je2(nc-1)
          is2(nc) = is1(nc-1); ie2(nc) = ie1(nc-1)
          js2(nc) = js1(nc-1); je2(nc) = je1(nc-1)
          align1(nc) = align2(nc-1); align2(nc) = align1(nc-1)
       end if
    end do

    !--- computing the overlap for the contact region. cyclic and folded-north boundary condition will
    !--- be part of contact region. 
    call define_contact_point( domain%T, nc, t1, t2, align1, align2, is1, ie1, js1, je1, &
                               is2, ie2, js2, je2, isgList, iegList, jsgList, jegList, mtype == CUBIC_GRID )

    call define_contact_point( domain%C, nc, t1, t2, align1, align2, is1, ie1, js1, je1, &
                               is2, ie2, js2, je2, isgList, iegList, jsgList, jegList, mtype == CUBIC_GRID )

    !--- for cubic grid, need to setup the communication for the corner to ensure consistency.
    if( mType == CUBIC_GRID ) call set_cubic_corner_point ( domain%C )

    call define_contact_point( domain%E, nc, t1, t2, align1, align2, is1, ie1, js1, je1, &
                               is2, ie2, js2, je2, isgList, iegList, jsgList, jegList, mtype == CUBIC_GRID )

    call define_contact_point( domain%N, nc, t1, t2, align1, align2, is1, ie1, js1, je1, &
                               is2, ie2, js2, je2, isgList, iegList, jsgList, jegList, mtype == CUBIC_GRID )

    !--- release memory
    deallocate(align1, align2, is1, ie1, js1, je1, is2, ie2, js2, je2  )
    deallocate(isgList, iegList, jsgList, jegList, t1, t2 )

  end subroutine mpp_define_mosaic

!#####################################################################

    subroutine compute_overlaps( domain, ishift, jshift )
!computes remote domain overlaps
!assumes only one in each direction
!will calculate the overlapping for T,E,C,N-cell seperately.
      type(domain2D), intent(inout) :: domain
      integer, intent(in)           :: ishift, jshift
      integer :: i, j, m, n
      integer :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff
      integer :: list, position, middle

      position = domain%position
      !--- when there is only one tile, n will equal to np
      n = size(domain%list(:))

      do i = 0, n-1
         do j = 1, 4
            domain%list(i)%send(:)%overlap(j)  = .FALSE.; domain%list(i)%recv(:)%overlap(j)  = .FALSE.
         end do
         domain%list(i)%send(:)%folded   = .FALSE.; domain%list(i)%recv(:)%folded   = .FALSE.
         domain%list(i)%send(:)%rotation = ZERO   ; domain%list(i)%recv(:)%rotation = ZERO
         domain%list(i)%send(:)%is       = 0      ; domain%list(i)%send(:)%ie       = -1 
         domain%list(i)%send(:)%js       = 0      ; domain%list(i)%send(:)%je       = -1 
         domain%list(i)%recv(:)%is       = 0      ; domain%list(i)%recv(:)%ie       = -1 
         domain%list(i)%recv(:)%js       = 0      ; domain%list(i)%recv(:)%je       = -1 
         domain%list(i)%send(:)%n        = 0      ; domain%list(i)%recv(:)%n        = 0
         domain%list(i)%overlap          = .FALSE.
      enddo

!send
      call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ioff, ysize=joff ) !cyclic offsets
      ioff = ioff - ishift
      joff = joff - jshift
      middle = (isg+ieg)/2+1

      do list = 0,n-1
         m = mod( domain%pos+list, n )
         if(domain%list(m)%tile_id == domain%tile_id) then  ! only compute the overlapping within tile.
!to_pe's eastern halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
            if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
                 .AND. ( jsc == je .or. jec == js ) ) then
               !--- do nothing, this point will come from other pe
            else 
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
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
                    .AND. is .GE. middle .AND. domain%list(m)%x%data%end .LE. ieg ) then
                  call define_rect_overlap( domain%list(m)%send(1), is, ie, js, je-1, isc, iec, jsc, jec)
                  is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
                  js = je
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                  end select
                  call define_list_overlap( domain%list(m)%send(1), is, ie, js, je, isc, iec, jsc, jec)
               else 
                  call define_rect_overlap( domain%list(m)%send(1), is, ie, js, je, isc, iec, jsc, jec)
               end if
            end if

!to_pe's SE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(2), is, ie, js, je, isc, iec, jsc, jec)

!to_pe's southern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(3), is, ie, js, je, isc, iec, jsc, jec)

!to_pe's SW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( jsg.GT.js )then
             if( domain%y%cyclic .AND. je.LT.jsc )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%send(4), is, ie, js, je, isc, iec, jsc, jec)

!to_pe's western halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
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
                          .AND. ( (domain%x%cyclic .AND. domain%list(m)%x%compute%begin == isg)            &
                           .OR. ( domain%list(m)%x%compute%begin-1 .GE. middle ) ) ) then
               call define_rect_overlap( domain%list(m)%send(5), is, ie, js, je-1, isc, iec, jsc, jec)
               !--- consider at j = jeg for west edge.
               !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
               if(position == CORNER .AND. .NOT. domain%symmetry .AND. domain%list(m)%x%compute%begin == isg) then
                   call define_list_overlap( domain%list(m)%send(5), ie, ie, je, je, isc, iec, jsc, jec)
               end if        
               is = domain%list(m)%x%data%begin;    ie = domain%list(m)%x%compute%begin-1
               je = domain%list(m)%y%compute%end;   js = je  
               if ( domain%list(m)%x%compute%begin == isg ) then
                  select case (position)
                  case(NORTH)
                     i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                  case(CORNER)
                     i=is; is = 2*isg-ie-2+2*ishift; ie = 2*isg-i-2+2*ishift  
                  end select
                  if(ie .GT. domain%x%compute%end) call mpp_error( FATAL, &
                       'mpp_domains_util_mod: west edge ubound error send.' )     
               else 
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                  end select
               end if
               call define_list_overlap( domain%list(m)%send(5), is, ie, js, je, isc, iec, jsc, jec)
            else
               call define_rect_overlap( domain%list(m)%send(5), is, ie, js, je, isc, iec, jsc, jec)
            end if            

!to_pe's NW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( isg.GT.is )then
             if( domain%x%cyclic .AND. ie.LT.isc )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
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
            call define_rect_overlap( domain%list(m)%send(6), is, ie, js, je, isc, iec, jsc, jec)
            !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%send(6), is, is, js, je, isc, iec, jsc, jec)
         endif

!to_pe's northern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
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
               call define_rect_overlap( domain%list(m)%send(7), is, ie, js, je, isc, iec, jsc, jec)
            end if
            !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%send(7), is, is, js, je, isc, iec, jsc, jec)
         endif

!to_pe's NE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( ie.GT.ieg )then
             if( domain%x%cyclic .AND. iec.LT.is )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if
         end if  
         if( je.GT.jeg )then
             if( domain%y%cyclic .AND. jec.LT.js )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
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
            call define_rect_overlap( domain%list(m)%send(8), is, ie, js, je, isc, iec, jsc, jec)
         end if
      end do  ! end of send set up.

      ! begin setting up recv       
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         if(domain%list(m)%tile_id == domain%tile_id) then  ! only compute the overlapping within tile.
         call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
!recv_e  
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
            if( domain%symmetry .AND. (position == NORTH .OR. position == CORNER ) &
                 .AND. ( jsd == je .or. jed == js ) ) then
               ! --- do nothing, this point will come from other pe
            else
         if( ied.GT.ieg )then
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
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
                  call define_rect_overlap( domain%list(m)%recv(1), is, ie, js, je, isd, ied, jsd, jed-1)
                  is=isc; ie=iec; js=jsc; je=jec
                  jsd = jed
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i  
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift      
                  end select
                  call define_list_overlap( domain%list(m)%recv(1), is, ie, js, je, isd, ied, jsd, jed, .TRUE.)
               else 
                  call define_rect_overlap( domain%list(m)%recv(1), is, ie, js, je, isd, ied, jsd, jed)
               end if
            end if
!recv_se     
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
         if( ied.GT.ieg )then
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(2), is, ie, js, je, isd, ied, jsd, jed)

!recv_s      
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(3), is, ie, js, je, isd, ied, jsd, jed)

!recv_sw 
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         if( jsd.LT.jsg )then
             if( domain%y%cyclic .AND. js.GT.jed )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
             end if  
         end if  
         if( isd.LT.isg )then
             if( domain%x%cyclic .AND. is.GT.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(4), is, ie, js, je, isd, ied, jsd, jed)

!recv_w      
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         if( isd.LT.isg )then
             if( domain%x%cyclic .AND. is.GT.ied )then !try cyclic offset
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
                           .AND. ( (domain%x%cyclic .AND. isd < isg ) .OR. ( ied .GE. middle ) ) )  then
               call define_rect_overlap( domain%list(m)%recv(5), is, ie, js, je, isd, ied, jsd, jed-1)
               !--- consider at j = jeg for west edge.
               !--- when the data is at corner and not symmetry, i = isg -1 will get from cyclic condition
               if( position == CORNER .AND. .NOT. domain%symmetry .AND. isd < isg ) then
                   call define_list_overlap( domain%list(m)%recv(5), is, ie, js, je, ied, ied, jed, jed)
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
                  if(ie .GT. domain%x%compute%end) call mpp_error( FATAL, &
                       'mpp_domains_util_mod: west edge ubound error recv.' )     
               else 
                  select case (position)
                  case(NORTH)
                     i=is; is = isg+ieg-ie; ie = isg+ieg-i
                  case(CORNER)
                     i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                  end select
               end if
               call define_list_overlap( domain%list(m)%recv(5), is, ie, js, je, isd, ied, jed, jed, .TRUE.)
            else
               call define_rect_overlap( domain%list(m)%recv(5), is, ie, js, je, isd, ied, jsd, jed)
            end if            

!recv_nw     
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(6)%folded = .TRUE.
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
               if( domain%x%cyclic .AND. is.GE.ied )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if 
            call define_rect_overlap( domain%list(m)%recv(6), is, ie, js, je, isd, ied, jsd, jed, domain%list(m)%recv(6)%folded)
            !--- when north edge is folded, is will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%recv(6), is, is, js, je, isd, ied, jsd, jed, domain%list(m)%recv(6)%folded )
         endif

!recv_n      
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(7)%folded = .TRUE.
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
               call define_rect_overlap( domain%list(m)%recv(7), is, ie, js, je, isd, ied, jsd, jed, domain%list(m)%recv(7)%folded)
            end if
            !--- when north edge is folded, ie will be less than isg when position is EAST and CORNER
            if(is .LT. isg .AND. domain%x%cyclic) then
               is = is + ioff
               call define_list_overlap( domain%list(m)%recv(7), is, is, js, je, isd, ied, jsd, jed, domain%list(m)%recv(7)%folded)
         endif
!recv_ne     
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec

         if( jed.GT.jeg )then
             if( domain%y%cyclic .AND. je.LT.jsd )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                  domain%list(m)%recv(8)%folded = .TRUE.
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
             if( domain%x%cyclic .AND. ie.LT.isd )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if  
         end if  
            call define_rect_overlap( domain%list(m)%recv(8), is, ie, js, je, isd, ied, jsd, jed, domain%list(m)%recv(8)%folded)
         endif
      enddo ! end of recv do loop

      !--- Now calculate the overlapping for fold-edge. Currently we only consider about folded-north
      !--- for folded-north-edge, only need to consider to_pe's north(7) direction 
      !--- only position at NORTH and CORNER need to be considered
      if( BTEST(domain%fold,NORTH) .AND. ( position == NORTH .OR. position == CORNER) ) then
         if( domain%y%data%begin .LE. jeg .AND. jeg .LE. domain%y%data%end )then !fold is within domain
            !--- calculate the overlapping for sending
            if( domain%x%pos .LT. (size(domain%x%list(:))+1)/2 )then 
               do list = 0,n-1
                  m = mod( domain%pos+list, n )
                  if(domain%list(m)%tile_id == domain%tile_id) then  ! only compute the overlapping within tile.
                     call mpp_get_compute_domain( domain, isc, iec, jsc, jec )                     
                     js = domain%list(m)%y%compute%end;   je = js
                     if( je == jeg )then   ! fold is within domain.
                        is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
                        select case (position)
                        case(NORTH)
                           is = max(is, middle)
                           i=is; is = isg+ieg-ie; ie = isg+ieg-i
                        case(CORNER)
                           is = max(is, middle)
                           i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                        end select
                        call define_list_overlap( domain%list(m)%send(7), is, ie, js, je, isc, iec, jsc, jec)           

                     end if
                  end if
               end do
            end if

            !--- calculating overlapping for receving on north
            if( domain%x%pos .GE. size(domain%x%list(:))/2 )then 
               do list = 0,n-1
                  m = mod( domain%pos+list, n )
                  if(domain%list(m)%tile_id == domain%tile_id) then  ! only compute the overlapping within tile.
                     call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
                     jsd = domain%y%compute%end;   jed = jsd
                     if( jed == jeg )then   ! fold is within domain.
                        isd = domain%x%compute%begin; ied = domain%x%compute%end
                        is=isc; ie=iec; js = jsc; je = jec
                        select case (position)
                        case(NORTH)
                           isd = max(isd, middle)
                           i=is; is = isg+ieg-ie; ie = isg+ieg-i
                        case(CORNER)
                           isd = max(isd, middle)
                           i=is; is = isg+ieg-ie-1+ishift; ie = isg+ieg-i-1+ishift
                        end select
                        call define_list_overlap( domain%list(m)%recv(7), is, ie, js, je, isd, ied, jsd, jed, .TRUE.) 
                     end if
                  end if
               end do
            end if
         end if
      end if

      do list = 0,n-1
         m = mod( domain%pos+list, n )
         if(domain%list(m)%tile_id == domain%tile_id) then  
            if(ANY(domain%list(m)%send(:)%overlap(1))) domain%list(m)%overlap(1) = .TRUE.
            if(ANY(domain%list(m)%send(:)%overlap(2))) domain%list(m)%overlap(1) = .TRUE.
            if(ANY(domain%list(m)%recv(:)%overlap(1))) domain%list(m)%overlap(1) = .TRUE.
            if(ANY(domain%list(m)%recv(:)%overlap(2))) domain%list(m)%overlap(1) = .TRUE.
         end if
      end do
      contains

      subroutine define_rect_overlap(overlap, is1, ie1, js1, je1, is2, ie2, js2, je2, reverse)
        type(overlapList), intent(inout) :: overlap
        integer,           intent(in)    :: is1, ie1, js1, je1, is2, ie2, js2, je2
        logical, optional, intent(in)    :: reverse
        logical :: is_reverse
        integer :: is, ie, js, je

        is_reverse = .FALSE.
        if(PRESENT(reverse)) is_reverse = reverse

        is = max(is1,is2); ie = min(ie1,ie2)
        js = max(js1,js2); je = min(je1,je2)

         if(ie.GE.is .AND. je.GE.js )then
            overlap%overlap(1) = .TRUE.
            overlap%is = is
            overlap%ie = ie
            overlap%js = js
            overlap%je = je
            if(is_reverse) overlap%rotation = ONE_HUNDRED_EIGHTY
         endif

      end subroutine define_rect_overlap

      subroutine define_list_overlap(overlap, is1, ie1, js1, je1, is2, ie2, js2, je2, reverse)
        type(overlapList), intent(inout) :: overlap
         integer,         intent(in)   :: is1, ie1, js1, je1, is2, ie2, js2, je2
        logical, optional, intent(in)    :: reverse

        integer              :: is, ie, js, je, n, n_old, i, j
        integer, allocatable :: tmp(:)
        logical :: is_reverse

        is_reverse = .FALSE.
        if(PRESENT(reverse)) is_reverse = reverse

         is = max(is1,is2); ie = min(ie1,ie2)
         js = max(js1,js2); je = min(je1,je2)

         if(ie.GE.is .AND. je.GE.js )then
           overlap%overlap(2) = .TRUE.
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

    end subroutine compute_overlaps


  !###################################################################################
  !--- here only defines overlapping between contact region. The corner point will be dealt later.
  subroutine define_contact_point( domain, num_contact, tile1, tile2, align1, align2, &
                                   istart1, iend1, jstart1, jend1, istart2, iend2,            &
                                   jstart2, jend2, isgList, iegList, jsgList, jegList, is_cubic )
    type(domain2D),     intent(inout) :: domain
    integer,               intent(in) :: num_contact
    integer, dimension(:), intent(in) :: tile1, tile2     ! tile number
    integer, dimension(:), intent(in) :: align1, align2   ! align direction of contact region
    integer, dimension(:), intent(in) :: istart1, iend1   ! i-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: jstart1, jend1   ! j-index in tile_1 of contact region
    integer, dimension(:), intent(in) :: istart2, iend2   ! i-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: jstart2, jend2   ! j-index in tile_2 of contact region
    integer, dimension(:), intent(in) :: isgList, iegList ! i-global domain of each tile
    integer, dimension(:), intent(in) :: jsgList, jegList ! j-global domain of each tile
    logical,               intent(in) :: is_cubic         ! indicate if the mosaic is cubic grid.

    integer              :: isgMe, iegMe, jsgMe, jegMe, isMe, ieMe, jsMe, jeMe
    integer              :: isgNbr, iegNbr, jsgNbr, jegNbr, isNbr, ieNbr, jsNbr, jeNbr
    integer              :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer              :: isc1, iec1, jsc1, jec1, isc2, iec2, jsc2, jec2
    integer              :: isd1, ied1, jsd1, jed1, isd2, ied2, jsd2, jed2
    integer              :: xhalo, yhalo, ishift, jshift, ioff, joff
    integer              :: nlist, list, n, m, l, position
    integer              :: tileMe, tileNbr, count
    integer              :: max_contact, numS, numR
    integer, allocatable :: iscSend(:), iecSend(:), jscSend(:), jecSend(:)
    integer, allocatable :: isdSend(:), iedSend(:), jsdSend(:), jedSend(:)
    integer, allocatable :: iscRecv(:), iecRecv(:), jscRecv(:), jecRecv(:)
    integer, allocatable :: isdRecv(:), iedRecv(:), jsdRecv(:), jedRecv(:)
    integer, allocatable :: alignMe(:), alignNbr(:), rotateSend(:), rotateRecv(:)
    integer, allocatable :: tileSend(:), tileRecv(:), alignSend(:)

    !--- the tile number, domain decomposition and halo size of current pe.
    tileMe = domain%tile_id
    isgMe = domain%x%global%begin;  iegMe = domain%x%global%end
    jsgMe = domain%y%global%begin;  jegMe = domain%y%global%end
    xhalo = domain%x%compute%begin - domain%x%data%begin
    yhalo = domain%y%compute%begin - domain%y%data%begin

    !--- find if there is an extra point in x and y direction depending on position
    ishift = 0; jshift = 0
    if( domain%symmetry ) then
       select case ( domain%position )
       case ( CORNER )
          ishift = 1; jshift = 1
       case ( EAST )
          ishift = 1
       case ( NORTH )
          jshift = 1
       end select
    end if

    nlist = size(domain%list(:))

    max_contact = 4*num_contact ! should be enough

    !--- isdSend and iscSend will figure out the overlapping for sending from current pe.
    !--- isdRecv and iscREcv will figure out the overlapping for recving onto current pe.
    allocate( iscSend(max_contact), iecSend(max_contact),  jscSend(max_contact),    jecSend(max_contact) )
    allocate( isdSend(max_contact), iedSend(max_contact),  jsdSend(max_contact),    jedSend(max_contact) )
    allocate( iscRecv(max_contact), iecRecv(max_contact),  jscRecv(max_contact),    jecRecv(max_contact) )
    allocate( isdRecv(max_contact), iedRecv(max_contact),  jsdRecv(max_contact),    jedRecv(max_contact) )
    allocate( alignMe(max_contact), alignNbr(max_contact), rotateSend(max_contact), rotateRecv(max_contact) )
    allocate( tileSend(max_contact), tileRecv(max_contact), alignSend(max_contact) )
    alignMe = CORNER;  alignNbr = CORNER
    rotateSend = ZERO; rotateRecv = ZERO

    !--- loop over all the contact region to figure out the index for overlapping region.
    count = 0
    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       count = count + 1
       tileNbr = tile2(n)
       isMe  = istart1(n); ieMe  = iend1(n)
       jsMe  = jstart1(n); jeMe  = jend1(n)
       isNbr = istart2(n); ieNbr = iend2(n)
       jsNbr = jstart2(n); jeNbr = jend2(n)
       alignMe(count) = align1(n); alignNbr(count) = align2(n); alignSend(count) = align2(n)
       isgNbr = isgList(tileNbr); iegNbr = iegList(tileNbr)
       jsgNbr = jsgList(tileNbr); jegNbr = jegList(tileNbr)
       
       !--- Based on the alignment, get the rotation for this contact.
       if (alignMe(count) == WEST .OR. alignMe(count) == EAST) then
          if (alignNbr(count) == SOUTH .OR. alignNbr(count) == NORTH) then
             rotateSend(count) = MINUS_NINETY; rotateRecv(count) = NINETY
          end if
       else if ( alignMe(count) == SOUTH .OR. alignMe(count) == NORTH ) then
          if (alignNbr(count) == WEST .OR. alignNbr(count) == EAST ) then
             rotateSend(count) = NINETY; rotateRecv(count) = MINUS_NINETY
          end if
       end if
       tileSend(count) = tileNbr; tileRecv(count) = tileNbr
       !--- define the overlapping region for send and receive.
        select case ( alignMe(count) )
        case(WEST)
           if( rotateSend(count) == ZERO ) then
              iscSend(count) = isMe;                      iecSend(count) = ieMe + xhalo - 1 + ishift
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + jshift
           else
              iscSend(count) = isMe;                      iecSend(count) = ieMe + yhalo - 1 + jshift
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + ishift
           end if
           isdRecv(count) = isMe - xhalo;              iedRecv(count) = ieMe - 1
           jsdRecv(count) = jsMe;                      jedRecv(count) = jeMe + jshift       
        case(EAST)
           if( rotateSend(count) == ZERO ) then
              iscSend(count) = isMe - xhalo + 1;          iecSend(count) = ieMe
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + jshift
           else
              iscSend(count) = isMe - yhalo + 1;          iecSend(count) = ieMe
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + ishift
           end if
           isdRecv(count) = isMe + 1;                  iedRecv(count) = ieMe + xhalo + ishift
           jsdRecv(count) = jsMe;                      jedRecv(count) = jeMe + jshift  
        case(SOUTH)
           if( rotateSend(count) == ZERO ) then
              iscSend(count) = isMe;                      iecSend(count) = ieMe + ishift
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + yhalo - 1 + jshift
           else
              iscSend(count) = isMe;                      iecSend(count) = ieMe + jshift
              jscSend(count) = jsMe;                      jecSend(count) = jeMe + xhalo - 1 + ishift
           end if
           isdRecv(count) = isMe;                      iedRecv(count) = ieMe + ishift
           jsdRecv(count) = jsMe - yhalo;              jedRecv(count) = jeMe - 1   
        case(NORTH)
           if( rotateSend(count) == ZERO ) then
              iscSend(count) = isMe;                      iecSend(count) = ieMe + ishift
              jscSend(count) = jsMe - yhalo + 1;          jecSend(count) = jeMe
           else
              iscSend(count) = isMe;                      iecSend(count) = ieMe + jshift
              jscSend(count) = jsMe - xhalo + 1;          jecSend(count) = jeMe
           end if
           isdRecv(count) = isMe;                      iedRecv(count) = ieMe + ishift
           jsdRecv(count) = jsMe + 1;                  jedRecv(count) = jeMe + yhalo + jshift
        end select

        select case ( alignNbr(count) )
        case(WEST)
           isdSend(count) = isNbr - xhalo;              iedSend(count) = ieNbr - 1
           jsdSend(count) = jsNbr;                      jedSend(count) = jeNbr + jshift
           if( rotateSend(count) == ZERO ) then
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + xhalo - 1 + ishift
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + jshift           
           else
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + yhalo - 1 + jshift
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + ishift    
           end if
        case(EAST)
           isdSend(count) = isNbr + 1;                  iedSend(count) = ieNbr + xhalo + ishift
           jsdSend(count) = jsNbr;                      jedSend(count) = jeNbr + jshift 
           if( rotateSend(count) == ZERO ) then  
              iscRecv(count) = isNbr - xhalo + 1;          iecRecv(count) = ieNbr
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + jshift
           else
              iscRecv(count) = isNbr - yhalo + 1;          iecRecv(count) = ieNbr
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + ishift
           end if
        case(SOUTH)
           isdSend(count) = isNbr;                      iedSend(count) = ieNbr + ishift
           jsdSend(count) = jsNbr - yhalo;              jedSend(count) = jeNbr - 1
           if( rotateSend(count) == ZERO ) then
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + ishift
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + yhalo - 1 + jshift 
           else
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + jshift
              jscRecv(count) = jsNbr;                      jecRecv(count) = jeNbr + xhalo - 1 + ishift 
           end if
        case(NORTH)
           isdSend(count) = isNbr;                      iedSend(count) = ieNbr + ishift
           jsdSend(count) = jsNbr + 1;                  jedSend(count) = jeNbr + yhalo + jshift
           if( rotateSend(count) == ZERO ) then
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + ishift
              jscRecv(count) = jsNbr - yhalo+ 1;           jecRecv(count) = jeNbr   
           else
              iscRecv(count) = isNbr;                      iecRecv(count) = ieNbr + jshift
              jscRecv(count) = jsNbr - xhalo+ 1;           jecRecv(count) = jeNbr  
           end if
        end select
    end do

    numS = count
    numR = count

    !--- figure out the index for corner overlapping
    !--- we assume no rotation for the contact region ( may remove this restricion in the future ) .
    !--- For cubic-grid, not overlapping on the corner tile.
    if( .NOT. is_cubic ) then 
        call compute_corner_overlap( xhalo, yhalo, ishift, jshift, tileMe, tile1, tile2, align1, align2, &
                                     istart1, iend1, jstart1 ,jend1, istart2, iend2, jstart2 ,jend2,     &
                                     isgList, iegList, jsgList, jegList,                                 &
                                     numS, numR, tileRecv, isdRecv, iedRecv, jsdRecv, jedRecv,           &
                                     iscRecv, iecRecv, jscRecv, jecRecv, tileSend, isdSend, iedSend,     &
                                     jsdSend, jedSend, iscSend, iecSend, jscSend, jecSend, alignSend )
    end if

    isc = domain%x%compute%begin; iec = domain%x%compute%end
    jsc = domain%y%compute%begin; jec = domain%y%compute%end
    !--- compute the overlapping for send.
    do n = 1, numS
       do list = 0, nlist-1
          m = mod( domain%pos+list, nlist )
          if( domain%list(m)%tile_id .NE. tileSend(n) ) cycle
          isc1 = isc; iec1 = iec; jsc1 = jsc; jec1 = jec
          if( rotateSend(n) .NE. ZERO .AND. domain%symmetry ) then
             select case ( domain%position )
             case ( EAST )   ! use N-cell
                iec1 = iec1 - 1; jec1 = jec1 + 1
             case ( NORTH )   ! use E-cell
                iec1 = iec1 + 1; jec1 = jec1 - 1
             end select
          end if
          isc1 = max(isc1, iscSend(n)); iec1 = min(iec1, iecSend(n))
          jsc1 = max(jsc1, jscSend(n)); jec1 = min(jec1, jecSend(n))  
          if( isc1 > iec1 .OR. jsc1 > jec1 ) cycle
          !--- loop over 8 direction to get the overlapping starting from east with clockwise.

          do l = 1, 8
             !--- get the to_pe's data domain.
             select case ( l )
             case ( 1 )  ! eastern halo
                if( alignNbr(n) .NE. EAST ) cycle
                isd = domain%list(m)%x%compute%end+1-ishift; ied = domain%list(m)%x%data%end
                jsd = domain%list(m)%y%compute%begin;        jed = domain%list(m)%y%compute%end
             case ( 2 )  ! southeast halo
                isd = domain%list(m)%x%compute%end+1; ied = domain%list(m)%x%data%end
                jsd = domain%list(m)%y%data%begin;    jed = domain%list(m)%y%compute%begin-1
                if(alignNbr(n) == EAST )  isd = isd - ishift
             case ( 3 )  ! southern halo
                if( alignNbr(n) .NE. SOUTH ) cycle
                isd = domain%list(m)%x%compute%begin; ied = domain%list(m)%x%compute%end
                jsd = domain%list(m)%y%data%begin;    jed = domain%list(m)%y%compute%begin-1
             case ( 4 )  ! southwest halo
                isd = domain%list(m)%x%data%begin;    ied = domain%list(m)%x%compute%begin-1
                jsd = domain%list(m)%y%data%begin;    jed = domain%list(m)%y%compute%begin-1
             case ( 5 )  ! western halo
                if( alignNbr(n) .NE. WEST ) cycle
                isd = domain%list(m)%x%data%begin;    ied = domain%list(m)%x%compute%begin-1
                jsd = domain%list(m)%y%compute%begin; jed = domain%list(m)%y%compute%end
             case ( 6 )  ! northwest halo
                isd = domain%list(m)%x%data%begin;    ied = domain%list(m)%x%compute%begin-1
                jsd = domain%list(m)%y%compute%end+1; jed = domain%list(m)%y%data%end
                if(alignNbr(n) == NORTH ) jsd = jsd - jshift
             case ( 7 )  ! northern halo
                if( alignNbr(n) .NE. NORTH ) cycle
                isd = domain%list(m)%x%compute%begin;        ied = domain%list(m)%x%compute%end
                jsd = domain%list(m)%y%compute%end+1-jshift; jed = domain%list(m)%y%data%end
             case ( 8 )  ! northeast halo
                isd = domain%list(m)%x%compute%end+1; ied = domain%list(m)%x%data%end
                jsd = domain%list(m)%y%compute%end+1; jed = domain%list(m)%y%data%end
                if(alignNbr(n) == NORTH ) jsd = jsd - jshift
                if(alignNbr(n) == EAST )  isd = isd - ishift
             end select
             isd = max(isd, isdSend(n)); ied = min(ied, iedSend(n))
             jsd = max(jsd, jsdSend(n)); jed = min(jed, jedSend(n))  
             if( isd > ied .OR. jsd > jed ) cycle
             ioff = 0; joff = 0
             select case ( alignNbr(n) )
             case ( WEST )
                joff = jsd - jsdSend(n)
             case ( EAST )
                joff = jsd - jsdSend(n)
!                if( l == 2 .OR. l == 8 ) ioff = ishift
             case ( SOUTH )
                ioff = isd - isdSend(n) 
             case ( NORTH )
                ioff = isd - isdSend(n) 
!                if( l == 6 .OR. l == 8 ) joff = jshift
             end select

             !--- get the index in current pe.
             select case ( rotateSend(n) )
             case ( ZERO )
                isc2 = iscSend(n) + ioff; iec2 = isc2 + ied - isd
                jsc2 = jscSend(n) + joff; jec2 = jsc2 + jed - jsd
             case ( NINETY )                     ! N -> W or S -> E
                iec2 = iecSend(n) - joff; isc2 = iec2 - ( jed - jsd ) 
                jsc2 = jscSend(n) + ioff; jec2 = jsc2 + ied - isd
             case ( MINUS_NINETY )               ! W -> N or E -> S
                isc2 = iscSend(n) + joff; iec2 = isc2 + jed - jsd
                jec2 = jecSend(n) - ioff; jsc2 = jec2 - ( ied - isd ) 
             end select
             position = domain%position
             if( rotateSend(n) .NE. ZERO ) then
                if ( position == EAST ) then
                   position = NORTH
                else if ( position == NORTH ) then
                   position = EAST 
                end if
             end if
             call define_overlap_region(domain%list(m)%send(l), isc1, iec1, jsc1, jec1, isc2, &
                  iec2, jsc2, jec2, alignMe(n), position, domain%symmetry, l )

             !--- define the overlap rotation for receiving.
             domain%list(m)%send(l)%rotation = rotateSend(n)
          end do
          if( ANY(domain%list(m)%send(:)%overlap(3)) ) domain%list(m)%overlap(3) = .true.
       end do
    end do

    !--- compute the overlapping for recv.
    do n = 1, numR
       do list = 0, nlist-1
          m = mod( domain%pos+list, nlist )
          if( domain%list(m)%tile_id .NE. tileRecv(n) ) cycle
          isc = domain%list(m)%x%compute%begin; iec = domain%list(m)%x%compute%end
          jsc = domain%list(m)%y%compute%begin; jec = domain%list(m)%y%compute%end
          if( rotateSend(n) .NE. ZERO .AND. domain%symmetry ) then
             select case ( domain%position )
             case ( EAST )   ! use N-cell
                iec = iec - 1; jec = jec + 1
             case ( NORTH )   ! use E-cell
                iec = iec + 1; jec = jec - 1
             end select
          end if
          isc = max(isc, iscRecv(n)); iec = min(iec, iecRecv(n))
          jsc = max(jsc, jscRecv(n)); jec = min(jec, jecRecv(n))      
          if( isc > iec .OR. jsc > jec ) cycle
          !--- find the offset for this overlapping.
          ioff = 0; joff = 0
          select case ( alignNbr(n) )
          case ( WEST, EAST )
             joff = jsc - jscRecv(n)
          case ( NORTH, SOUTH )
             ioff = isc - iscRecv(n) 
          end select

          !--- get the index in current pe.
          select case ( rotateRecv(n) )
          case ( ZERO )
             isd1 = isdRecv(n) + ioff; ied1 = isd1 + iec - isc
             jsd1 = jsdRecv(n) + joff; jed1 = jsd1 + jec - jsc
          case ( NINETY )                      ! N -> W or S -> E
             isd1 = isdRecv(n);        ied1 = isd1 + jec - jsc
             jed1 = jedRecv(n) - ioff; jsd1 = jed1 - ( iec - isc ) 
          case ( MINUS_NINETY )                ! W -> N or E -> S
             ied1 = iedRecv(n) - joff; isd1 = ied1 - ( jec - jsc ) 
             jsd1 = jsdRecv(n);        jed1 = jsd1 + iec - isc
          end select

          !--- loop over 8 direction to get the overlapping starting from east with clockwise.
          do l = 1, 8
             select case ( l )
             case ( 1 )  ! eastern halo
                if( alignMe(n) .NE. EAST ) cycle
                isd2 = domain%x%compute%end+1-ishift; ied2 = domain%x%data%end
                jsd2 = domain%y%compute%begin;        jed2 = domain%y%compute%end
             case ( 2 )  ! southeast halo
                isd2 = domain%x%compute%end+1; ied2 = domain%x%data%end
                jsd2 = domain%y%data%begin;    jed2 = domain%y%compute%begin-1
                if(alignMe(n) == EAST) isd2 = isd2 - ishift  
             case ( 3 )  ! southern halo
                if( alignMe(n) .NE. SOUTH ) cycle
                isd2 = domain%x%compute%begin; ied2 = domain%x%compute%end
                jsd2 = domain%y%data%begin;    jed2 = domain%y%compute%begin-1
             case ( 4 )  ! southwest halo
                isd2 = domain%x%data%begin;    ied2 = domain%x%compute%begin-1
                jsd2 = domain%y%data%begin;    jed2 = domain%y%compute%begin-1
             case ( 5 )  ! western halo
                if( alignMe(n) .NE. WEST ) cycle
                isd2 = domain%x%data%begin;    ied2 = domain%x%compute%begin-1
                jsd2 = domain%y%compute%begin; jed2 = domain%y%compute%end
             case ( 6 )  ! northwest halo
                isd2 = domain%x%data%begin;    ied2 = domain%x%compute%begin-1
                jsd2 = domain%y%compute%end+1; jed2 = domain%y%data%end
                if(alignMe(n) == NORTH) jsd2 = jsd2 -jshift
             case ( 7 )  ! northern halo
                if( alignMe(n) .NE. NORTH ) cycle
                isd2 = domain%x%compute%begin;        ied2 = domain%x%compute%end
                jsd2 = domain%y%compute%end+1-jshift; jed2 = domain%y%data%end
             case ( 8 )  ! northeast halo
                isd2 = domain%x%compute%end+1; ied2 = domain%x%data%end
                jsd2 = domain%y%compute%end+1; jed2 = domain%y%data%end
                if(alignMe(n) == NORTH ) jsd2 = jsd2 - jshift
                if(alignMe(n) == EAST )  isd2 = isd2 - ishift
             end select

             call define_overlap_region(domain%list(m)%recv(l), isd1, ied1, jsd1, jed1, isd2, &
                  ied2, jsd2, jed2, alignMe(n), domain%position, domain%symmetry, l )

             !--- define the overlap rotation for receiving.
             domain%list(m)%recv(l)%rotation = rotateRecv(n)
          end do
          if( ANY(domain%list(m)%recv(:)%overlap(3)) ) domain%list(m)%overlap(3) = .true.
       end do
    end do

    !--- release memory
    deallocate( iscSend, iecSend, jscSend, jecSend, isdSend, iedSend, jsdSend, jedSend )
    deallocate( iscRecv, iecRecv, jscRecv, jecRecv, isdRecv, iedRecv, jsdRecv, jedRecv )
    deallocate( alignMe, alignNbr, rotateSend, rotateRecv, tileSend, tileRecv, alignSend )

    contains

       subroutine define_overlap_region(overlap, is1, ie1, js1, je1, is2, ie2, js2, je2, alignment, &
                                        position, symmetry, direction)
         type(overlapList), intent(inout) :: overlap
         integer,            intent(in)   :: is1, ie1, js1, je1, is2, ie2, js2, je2
         integer,            intent(in)   :: alignment, position, direction
         logical,            intent(in)   :: symmetry

         integer :: is, ie, js, je
         is = max(is1,is2); ie = min(ie1,ie2)
         js = max(js1,js2); je = min(je1,je2)

         !--- remove unnecessary communication for symmetry domain.
         if( mod(direction,2) .NE. 0 ) then  ! east, south, west, north
            if( symmetry ) then
               select case ( alignment )
               case ( WEST, EAST )
                  if ( js == je .AND. ( position == CORNER .OR. position == NORTH ) ) return
               case ( SOUTH, NORTH )
                  if ( is == ie .AND. ( position == CORNER .OR. position == EAST ) ) return
               end select
            end if
         end if
         if(ie.GE.is .AND. je.GE.js )then
            overlap%overlap(3) = .TRUE.
            overlap%is = is
            overlap%ie = ie
            overlap%js = js
            overlap%je = je
         endif

       end subroutine define_overlap_region

     end subroutine define_contact_point

     !---------------------------------------------------------------------
     !--- this routine is used to calculate overlap between tiles for C,E,N-cell.
     
     subroutine set_contact_point ( domain_in, domain_out )
       type(domain2d),    intent(in) :: domain_in
       type(domain2d), intent(inout) :: domain_out

       integer :: nlist, ishift, jshift, m, l
       integer :: isoff, ieoff, jsoff, jeoff

       nlist = size(domain_in%list(:) )

       ishift = 0; jshift = 0
       if( domain_out%symmetry ) then
          select case ( domain_out%position )
          case ( CORNER )
             ishift = 1; jshift = 1
          case ( EAST )
             ishift = 1
          case ( NORTH )
             jshift = 1
          end select
       end if

       do m = 0, nlist-1
          do l = 1, 8      ! loop over 8 directions
             select case ( l )
             case ( 1 )                                                       !  east 
                isoff = 0;      ieoff = ishift; jsoff = 0;      jeoff = jshift
             case ( 2 )                                                       !  southeast 
                isoff = ishift; ieoff = ishift; jsoff = 0;      jeoff = 0
             case ( 3  )                                                      !  south 
                isoff = 0;      ieoff = ishift; jsoff = 0;      jeoff = 0
             case ( 4 )                                                       !  southwest 
                isoff = 0;      ieoff = 0;      jsoff = 0;      jeoff = 0
             case ( 5  )                                                      !  west 
                isoff = 0;      ieoff = 0;      jsoff = 0;      jeoff = jshift
             case ( 6 )                                                       !  northwest 
                isoff = 0;      ieoff = 0;      jsoff = jshift; jeoff = jshift
             case ( 7 )                                                       !  north
                isoff = 0;      ieoff = ishift; jsoff = 0;      jeoff = jshift
             case ( 8 )                                                       !  northeast 
                isoff = ishift; ieoff = ishift; jsoff = jshift; jeoff = jshift
             end select

             if( domain_in%list(m)%recv(l)%overlap(3) ) then
                domain_out%list(m)%recv(l)%overlap(3) = .TRUE.
                domain_out%list(m)%recv(l)%rotation   = domain_in%list(m)%recv(l)%rotation
                domain_out%list(m)%recv(l)%is         = domain_in%list(m)%recv(l)%is + isoff
                domain_out%list(m)%recv(l)%ie         = domain_in%list(m)%recv(l)%ie + ieoff
                domain_out%list(m)%recv(l)%js         = domain_in%list(m)%recv(l)%js + jsoff
                domain_out%list(m)%recv(l)%je         = domain_in%list(m)%recv(l)%je + jeoff
             end if

             if( domain_in%list(m)%send(l)%rotation .NE. ZERO ) then ! rotate 90 or -90 degree ( cubic-grid )
                ! offset will be different for some direction
                select case ( l )
                case ( 1 )                                                       !  east 
                   isoff = 0;      ieoff = jshift; jsoff = 0;      jeoff = ishift
                case ( 2 )                                                       !  southeast 
                   isoff = ishift; ieoff = ishift; jsoff = 0;      jeoff = 0
                case ( 3 )                                                       !  south
                   isoff = 0;      ieoff = 0;      jsoff = 0;      jeoff = ishift
                case ( 4 )                                                       !  southwest 
                   isoff = ishift; ieoff = ishift; jsoff = 0;      jeoff = 0
                case ( 5 )                                                      !  west no
                   isoff = 0;      ieoff = jshift; jsoff = 0;      jeoff = 0
                case ( 6 )                                                       !  northwest 
                   isoff = 0;      ieoff = 0;      jsoff = 0;      jeoff = 0
                case ( 7 )                                                       !  north
                   isoff = 0;      ieoff = jshift; jsoff = 0;      jeoff = ishift
                case ( 8 )                                                       !  northeast 
                   isoff = ishift; ieoff = ishift; jsoff = jshift; jeoff = jshift
                end select
             end if
             if( domain_in%list(m)%send(l)%overlap(3) ) then
                domain_out%list(m)%send(l)%overlap(3) = .TRUE.
                domain_out%list(m)%send(l)%rotation   = domain_in%list(m)%send(l)%rotation
                domain_out%list(m)%send(l)%is         = domain_in%list(m)%send(l)%is + isoff
                domain_out%list(m)%send(l)%ie         = domain_in%list(m)%send(l)%ie + ieoff
                domain_out%list(m)%send(l)%js         = domain_in%list(m)%send(l)%js + jsoff
                domain_out%list(m)%send(l)%je         = domain_in%list(m)%send(l)%je + jeoff
             end if


          end do
          domain_out%list(m)%overlap(3) = domain_in%list(m)%overlap(3)
       end do

     end subroutine set_contact_point

     !--- set up the overlapping of cubic grid corner point
     !--- for data on corner points of cubic, need to make sure three tiles have the same value.
     !   --- the 8 corner points will be decided by the following.
     !   --- 1,2,3: take value at 3, need to send southwest point from 3 to northeast point at 1.
     !   --- 2,3,4: take value at 4, need to send southwest point from 4 to northeast point at 2.
     !   --- 3,4,5: take value at 5, need to send southwest point from 5 to northeast point at 3.
     !   --- 4,5,6: take value at 6, need to send southwest point from 6 to northeast point at 4.
     !   --- 5,6,1: take value at 1, need to send southwest point from 1 to northeast point at 5.
     !   --- 6,1,2: take value at 2, need to send southwest point from 2 to northeast point at 6.
     !   --- 1,3,5: take value at 1, need to send northwest point from 1 to northwest point at 3 and 5.
     !   --- 2,4,6: take value at 2, need to send southeast point from 2 to southeast point at 4 and 6.

     subroutine set_cubic_corner_point ( domain )
       type(domain2d), intent(inout) :: domain

       integer :: isc1, iec1, jsc1, jec1, isc2, jsc2, iec2, jec2
       integer :: isgMe, iegMe, jsgMe, jegMe, isgNbr, iegNbr, jsgNbr, jegNbr
       integer :: tileMe, tileNbr, nlist, list, m

       isc1 = domain%x%compute%begin; iec1 = domain%x%compute%end
       jsc1 = domain%y%compute%begin; jec1 = domain%y%compute%end
       isgMe = domain%x%global%begin;  iegMe = domain%x%global%end
       jsgMe = domain%y%global%begin;  jegMe = domain%y%global%end
       tileMe = domain%tile_id
       nlist = size(domain%list(:))
       !--- set up send
       if( isc1 == isgMe .AND. jsc1 == jsgMe ) then  ! --- to_pe's northeast corner
          tileNbr = tileMe+4
          if( tileNbr > 6) tileNbr = tileNbr - 6
          !--- find the pe to send to.
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. tileNbr ) cycle
             iec2   = domain%list(m)%x%compute%end; jec2   = domain%list(m)%y%compute%end                         
             iegNbr = domain%list(m)%x%global%end;  jegNbr = domain%list(m)%y%global%end  
             if( iec2 == iegNbr .AND. jec2 == jegNbr ) then
                domain%list(m)%send(8)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%send(8)%i2 = isc1;  domain%list(m)%send(8)%j2 = jsc1
             end if
          end do
       end if
       !--- tile 1 to 3 and 5
       if( tileMe == 1 .AND. isc1 == isgMe .AND. jec1 == jegMe ) then !--- to_pe's northwest corner
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. 3 .AND. domain%list(m)%tile_id .NE. 5 ) cycle 
             isc2   = domain%list(m)%x%compute%begin; jec2   = domain%list(m)%y%compute%end                         
             isgNbr = domain%list(m)%x%global%begin;  jegNbr = domain%list(m)%y%global%end
             if( isc2 == isgNbr .AND. jec2 == jegNbr ) then
                domain%list(m)%send(6)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%send(6)%i2 = isc1;  domain%list(m)%send(6)%j2 = jec1
             end if
          end do
       end if
       !--- tile 2 to 4 and 6
       if( tileMe == 2 .AND. iec1 == iegMe .AND. jsc1 == jsgMe ) then !--- to_pe's southeast corner
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. 4 .AND. domain%list(m)%tile_id .NE. 6 ) cycle 
             iec2   = domain%list(m)%x%compute%end; jsc2   = domain%list(m)%y%compute%begin                        
             iegNbr = domain%list(m)%x%global%end;  jsgNbr = domain%list(m)%y%global%begin
             if( iec2 == iegNbr .AND. jsc2 == jsgNbr ) then
                domain%list(m)%send(2)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%send(2)%i2 = iec1;  domain%list(m)%send(2)%j2 = jsc1
             end if
          end do
       end if
       !--- set up recv
       if( iec1 == iegMe .AND. jec1 == jegMe ) then  ! --- northeast corner
          tileNbr = tileMe + 2
          if( tileNbr > 6) tileNbr = tileNbr - 6
          !--- find the pe to send to.
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. tileNbr ) cycle
             isc2   = domain%list(m)%x%compute%begin; jsc2   = domain%list(m)%y%compute%begin                        
             isgNbr = domain%list(m)%x%global%begin;  jsgNbr = domain%list(m)%y%global%begin  
             if( isc2 == isgNbr .AND. jsc2 == jsgNbr ) then
                domain%list(m)%recv(8)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%recv(8)%i2 = iec1;  domain%list(m)%recv(8)%j2 = jec1
             end if
          end do
       end if
       !--- tile  3 and 5 from 1
       if( (tileMe == 3 .or. tileMe == 5) .AND. isc1 == isgMe .AND. jec1 == jegMe ) then !--- northwest corner
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. 1 ) cycle 
             isc2   = domain%list(m)%x%compute%begin; jec2   = domain%list(m)%y%compute%end                         
             isgNbr = domain%list(m)%x%global%begin;  jegNbr = domain%list(m)%y%global%end
             if( isc2 == isgNbr .AND. jec2 == jegNbr ) then
                domain%list(m)%recv(6)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%recv(6)%i2 = isc1;  domain%list(m)%recv(6)%j2 = jec1
             end if
          end do
       end if
       !--- tile 4 and 6 from 2.
       if( (tileMe == 4 .OR. tileMe == 6) .AND. iec1 == iegMe .AND. jsc1 == jsgMe ) then !--- to_pe's southeast corner
          do list = 0, nlist-1
             m = mod( domain%pos+list, nlist )
             if( domain%list(m)%tile_id .NE. 2 ) cycle 
             iec2   = domain%list(m)%x%compute%end; jsc2   = domain%list(m)%y%compute%begin                        
             iegNbr = domain%list(m)%x%global%end;  jsgNbr = domain%list(m)%y%global%begin
             if( iec2 == iegNbr .AND. jsc2 == jsgNbr ) then
                domain%list(m)%recv(2)%overlap(4) = .true.
                domain%list(m)%overlap(4) = .true.
                domain%list(m)%recv(2)%i2 = iec1;  domain%list(m)%recv(2)%j2 = jsc1
             end if
          end do
       end if

     end subroutine set_cubic_corner_point

     !--- this routine is used to compute the overlap between corner tile.
     !--- the corner tile overlapped is not specified through contact region
     !--- 
     subroutine compute_corner_overlap( xhalo, yhalo, ishift, jshift, tileMe, tile1, tile2, align1, align2,   &
                                        istart1, iend1, jstart1 ,jend1, istart2, iend2, jstart2 ,jend2,       &
                                        isgList, iegList, jsgList, jegList, numS, numR, tileRecv,             &
                                        isdRecv, iedRecv, jsdRecv, jedRecv,                                   &
                                        iscRecv, iecRecv, jscRecv, jecRecv, tileSend, isdSend, iedSend,       &
                                        jsdSend, jedSend, iscSend, iecSend, jscSend, jecSend, alignSend )     
    integer,                  intent(in) :: xhalo, yhalo, ishift, jshift, tileMe
    integer, dimension(:),    intent(in) :: tile1, tile2     ! tile number
    integer, dimension(:),    intent(in) :: align1, align2   ! alignment
    integer, dimension(:),    intent(in) :: istart1, iend1   ! i-index in tile_1 of contact region
    integer, dimension(:),    intent(in) :: jstart1, jend1   ! j-index in tile_1 of contact region
    integer, dimension(:),    intent(in) :: istart2, iend2   ! i-index in tile_2 of contact region
    integer, dimension(:),    intent(in) :: jstart2, jend2   ! j-index in tile_2 of contact region
    integer, dimension(:),    intent(in) :: isgList, iegList ! global domain in i-direction of each tile
    integer, dimension(:),    intent(in) :: jsgList, jegList ! global domain in j-direction of each tile
    integer,               intent(inout) :: numS, numR
    integer, dimension(:), intent(inout) :: tileRecv         ! tile to recv from
    integer, dimension(:), intent(inout) :: isdRecv, iedRecv ! i-index on current tile: nbr --> current 
    integer, dimension(:), intent(inout) :: jsdRecv, jedRecv ! j-index on current tile: nbr --> current
    integer, dimension(:), intent(inout) :: iscRecv, iecRecv ! i-index on neighbor tile:  nbr --> current
    integer, dimension(:), intent(inout) :: jscRecv, jecRecv ! j-index on neighbor tile:  nbr --> current
    integer, dimension(:), intent(inout) :: tileSend         ! tile to send to
    integer, dimension(:), intent(inout) :: isdSend, iedSend ! i-index on neighbor tile: current --> nbr
    integer, dimension(:), intent(inout) :: jsdSend, jedSend ! j-index on neighbor tile: current --> nbr
    integer, dimension(:), intent(inout) :: iscSend, iecSend ! i-index on current tile:  current --> nbr
    integer, dimension(:), intent(inout) :: jscSend, jecSend ! j-index on current tile:  current --> nbr
    integer, dimension(:), intent(inout) :: alignSend        ! alignment of to-tile's 

    integer              :: n, m, l, ioff, joff
    integer              :: tileNbr, tileCor, num_contact
    logical              :: found_corner       
    integer              :: istart, iend, jstart, jend, alignNbr

    num_contact = size( istart1(:))
    !------------------------------------------------------------------------------
    !                            southeast corner
    !------------------------------------------------------------------------------
    found_corner = .FALSE.
    ioff = 0; joff = 0
    !--- define the srecv overlap for southeast corner.
    tileNbr = get_neighbor_tile( tileMe, EAST, tile1, tile2, align1, align2, jstart1, jstart2, jsgList, n )
    if( tileNbr > 0 ) then ! found nbr
       jstart = jstart2(n)
       if( jstart == jsgList(tileNbr) ) then
          tileCor = get_neighbor_tile( tileNbr, SOUTH, tile1, tile2, align1, align2, istart1, istart2, isgList, n )
          if(tileCor > 0) found_corner = .true.
       else
          joff = jegList(tileNbr) - jstart + 1
          tileCor = tileNbr
          found_corner = .true.
       end if
    else ! not found nbr for east contact, then found nbr for south contact
       tileNbr = get_neighbor_tile( tileMe, SOUTH, tile1, tile2, align1, align2, iend1, iend2, iegList, n)
       if( tileNbr > 0 ) then ! found nbr in south contact
          iend = iend2(n)
          if( iend == iegList(tileNbr) ) then
             tileCor = get_neighbor_tile( tileNbr, EAST, tile1, tile2, align1, align2, jend1, jend2, iegList, n )
             if(tileCor > 0) found_corner = .true.
          else
             ioff = iend - isgList(tileNbr) + 1
             tileCor = tileNbr
             found_corner = .true.
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tileCor
       isdRecv(numR) = iegList(tileMe) + 1 + ishift;        iedRecv(numR) = iegList(tileMe) + xhalo + ishift
       jsdRecv(numR) = jsgList(tileMe) - yhalo;             jedRecv(numR) = jsgList(tileMe) - 1
       iscRecv(numR) = isgList(tileCor) + ioff;             iecRecv(numR) = isgList(tileCor) + ioff + xhalo - 1
       jscRecv(numR) = jegList(tileCor) - joff - yhalo + 1; jecRecv(numR) = jegList(tileCor) - joff
    end if

    !--- define the send for to_pe's southeast corner.
    !--- loop over all the contact 
    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .false.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       jstart   = jstart1(n)
       jend     = jend1(n)
       if(alignNbr == EAST) then!--- to_pe's east contact
          if( jstart > jsgList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             ioff = 0; joff = jegList(tileMe) - jstart + 1
          end if
          if( jend == jegList(tileMe) ) then  ! may need to send data to north of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, NORTH, tile1, tile2, align1, align2, iend1, iend2, iegList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   iend = iend1(m)
                else
                   iend = iend2(m)
                end if
                if( iend == iegList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                end if
             end if
          end if

          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = SOUTH_EAST
             isdSend(numS) = iegList(tileCor) + 1 + ishift;      iedSend(numS) = iegList(tileCor) + xhalo + ishift
             jsdSend(numS) = jsgList(tileCor) - yhalo;           jedSend(numS) = jsgList(tileCor) - 1
             iscSend(numS) = isgList(tileMe) + ioff;             iecSend(numS) = isgList(tileMe) + ioff + xhalo - 1
             jscSend(numS) = jegList(tileMe) - joff - yhalo + 1; jecSend(numS) = jegList(tileMe) - joff
          end if
       end if
    end do

    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .FALSE.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       istart   = istart1(n)
       iend     = iend1(n)
       if(alignNbr == SOUTH) then !--- to_pe's SOUTH contact
          if( iend < iegList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             joff = 0; ioff =  iend - isgList(tileMe) + 1
          end if
          if( istart == isgList(tileMe) ) then  ! may need to send data to west of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, WEST, tile1, tile2, align1, align2, jstart1, jstart2, jsgList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   jstart = jstart1(m)
                else
                   jstart = jstart2(m)
                end if
                if( jstart == jsgList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                   do l = 1, numS
                      if( tileSend(l) == tileCor .AND. alignSend(l) == SOUTH_EAST ) then
                         found_corner = .false.
                         exit
                      end if
                   end do
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = SOUTH_EAST
             isdSend(numS) = iegList(tileCor) + 1 + ishift;      iedSend(numS) = iegList(tileCor) + xhalo + ishift
             jsdSend(numS) = jsgList(tileCor) - yhalo;           jedSend(numS) = jsgList(tileCor) - 1
             iscSend(numS) = isgList(tileMe) + ioff;             iecSend(numS) = isgList(tileMe) + ioff + xhalo - 1
             jscSend(numS) = jegList(tileMe) - joff - yhalo + 1; jecSend(numS) = jegList(tileMe) - joff
          end if
       end if
    end do

    !------------------------------------------------------------------------------
    !                            southwest corner
    !------------------------------------------------------------------------------
    found_corner = .FALSE.
    ioff = 0; joff = 0
    !--- define the recv overlap for southwest corner.
    tileNbr = get_neighbor_tile( tileMe, WEST, tile1, tile2, align1, align2, jstart1, jstart2, jsgList, n )
    if( tileNbr > 0 ) then ! found nbr
       jstart = jstart2(n)
       if( jstart == jsgList(tileNbr) ) then
          tileCor = get_neighbor_tile( tileNbr, SOUTH, tile1, tile2, align1, align2, iend1, iend2, iegList, n )
          if(tileCor > 0) found_corner = .true.
       else
          joff = jegList(tileNbr) - jstart + 1
          tileCor = tileNbr
          found_corner = .true.
       end if
    else ! not found nbr for west contact, then found nbr for south contact
       tileNbr = get_neighbor_tile( tileMe, SOUTH, tile1, tile2, align1, align2, istart1, istart2, isgList, n)
       if( tileNbr > 0 ) then ! found nbr in south contact
          istart = istart2(n)
          if( istart == isgList(tileNbr) ) then
             tileCor = get_neighbor_tile( tileNbr, WEST, tile1, tile2, align1, align2, jend1, jend2, iegList, n )
             if(tileCor > 0) found_corner = .true.
          else
             ioff = iegList(tileNbr) - istart + 1
             tileCor = tileNbr
             found_corner = .true.
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tileCor
       isdRecv(numR) = isgList(tileMe) -xhalo;              iedRecv(numR) = iegList(tileMe) - 1
       jsdRecv(numR) = jsgList(tileMe) - yhalo;             jedRecv(numR) = jsgList(tileMe) - 1
       iscRecv(numR) = iegList(tileCor) - ioff - xhalo + 1; iecRecv(numR) = iegList(tileCor) - ioff
       jscRecv(numR) = jegList(tileCor) - joff - yhalo + 1; jecRecv(numR) = jegList(tileCor) - joff
    end if

    !--- define the send for to_pe's southwest corner.
    !--- loop over all the contact 
    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .false.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       jstart   = jstart1(n)
       jend     = jend1(n)
       if(alignNbr == WEST) then!--- to_pe's west contact
          if( jstart > jsgList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             ioff = 0; joff = jegList(tileMe) - jstart + 1
          end if
          if( jend == jegList(tileMe) ) then  ! may need to send data to north of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, NORTH, tile1, tile2, align1, align2, istart1, istart2, isgList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   istart = istart1(m)
                else
                   istart = istart2(m)
                end if
                if( istart == isgList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = SOUTH_WEST
             isdSend(numS) = isgList(tileCor) - xhalo;           iedSend(numS) = isgList(tileCor) - 1
             jsdSend(numS) = jsgList(tileCor) - yhalo;           jedSend(numS) = jsgList(tileCor) - 1
             iscSend(numS) = iegList(tileMe) - ioff - xhalo + 1; iecSend(numS) = iegList(tileMe) - ioff
             jscSend(numS) = jegList(tileMe) - joff - yhalo + 1; jecSend(numS) = jegList(tileMe) - joff
          end if
       end if
    end do

    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .FALSE.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       istart   = istart1(n)
       iend     = iend1(n)
       if(alignNbr == SOUTH) then !--- to_pe's SOUTH contact
          if( istart > isgList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             joff = 0; ioff = iegList(tileMe) - istart + 1
          end if
          if( iend == iegList(tileMe) ) then  ! may need to send data to east of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, EAST, tile1, tile2, align1, align2, jstart1, jstart2, jsgList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   jstart = jstart1(m)
                else
                   jstart = jstart2(m)
                end if
                if( jstart == jsgList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                   do l = 1, numS 
                      if( tileSend(l) == tileCor .AND. alignSend(l) == SOUTH_WEST ) then
                         found_corner = .false.
                         exit
                      end if
                   end do
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = SOUTH_WEST
             isdSend(numS) = isgList(tileCor) - xhalo;           iedSend(numS) = isgList(tileCor) - 1
             jsdSend(numS) = jsgList(tileCor) - yhalo;           jedSend(numS) = jsgList(tileCor) - 1
             iscSend(numS) = iegList(tileMe) - ioff - xhalo + 1; iecSend(numS) = iegList(tileMe) - ioff
             jscSend(numS) = jegList(tileMe) - joff - yhalo + 1; jecSend(numS) = jegList(tileMe) - joff
          end if
       end if
    end do

    !------------------------------------------------------------------------------
    !                            northeast corner
    !------------------------------------------------------------------------------
    found_corner = .FALSE.
    ioff = 0; joff = 0
    !--- define the recv overlap for northeast corner.
    tileNbr = get_neighbor_tile( tileMe, EAST, tile1, tile2, align1, align2, jend1, jend2, jegList, n )
    if( tileNbr > 0 ) then ! found nbr
       jend = jend2(n)
       if( jend == jegList(tileNbr) ) then
          tileCor = get_neighbor_tile( tileNbr, NORTH, tile1, tile2, align1, align2, istart1, istart2, isgList, n )
          if(tileCor > 0) found_corner = .true.
       else
          joff = jend - jsgList(tileNbr) + 1
          tileCor = tileNbr
          found_corner = .true.
       end if
    else ! not found nbr for east contact, then found nbr for south contact
       tileNbr = get_neighbor_tile( tileMe, NORTH, tile1, tile2, align1, align2, iend1, iend2, iegList, n)
       if( tileNbr > 0 ) then ! found nbr in south contact
          iend = iend2(n)
          if( iend == iegList(tileNbr) ) then
             tileCor = get_neighbor_tile( tileNbr, EAST, tile1, tile2, align1, align2, jstart1, jstart2, isgList, n )
             if(tileCor > 0) found_corner = .true.
          else
             ioff = iend - isgList(tileNbr) + 1
             tileCor = tileNbr
             found_corner = .true.
          end if
       end if
    end if
    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tileCor
       isdRecv(numR) = iegList(tileMe) + 1 + ishift;        iedRecv(numR) = iegList(tileMe) + xhalo + ishift
       jsdRecv(numR) = jegList(tileMe) + 1 + jshift;        jedRecv(numR) = jegList(tileMe) + yhalo + jshift
       iscRecv(numR) = isgList(tileCor) + ioff;             iecRecv(numR) = isgList(tileCor) + ioff + xhalo - 1
       jscRecv(numR) = jsgList(tileCor) + joff;             jecRecv(numR) = jsgList(tileCor) + joff + xhalo - 1
    end if

    !--- define the send for to_pe's southeast corner.
    !--- loop over all the contact 
    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .false.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       jstart   = jstart1(n)
       jend     = jend1(n)
       if(alignNbr == EAST) then!--- to_pe's east contact
          if( jend < jegList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             ioff = 0; joff = jend - jsgList(tileMe) + 1
          end if
          if( jstart == jsgList(tileMe) ) then  ! may need to send data to north of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, SOUTH, tile1, tile2, align1, align2, iend1, iend2, iegList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   iend = iend1(m)
                else
                   iend = iend2(m)
                end if
                if( iend == iegList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                end if
             end if
          end if

          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = NORTH_EAST
             isdSend(numS) = iegList(tileCor) + 1 + ishift;      iedSend(numS) = iegList(tileCor) + xhalo + ishift
             jsdSend(numS) = jegList(tileCor) + 1 + jshift;      jedSend(numS) = jegList(tileCor) + yhalo + jshift
             iscSend(numS) = isgList(tileMe) + ioff;             iecSend(numS) = isgList(tileMe) + ioff + xhalo - 1
             jscSend(numS) = jsgList(tileMe) + joff;             jecSend(numS) = jsgList(tileMe) + joff + yhalo - 1
          end if
       end if
    end do

    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .FALSE.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       istart   = istart1(n)
       iend     = iend1(n)
       if(alignNbr == NORTH) then !--- to_pe's NORTH contact
          if( iend < iegList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             joff = 0; ioff =  iend - isgList(tileMe) + 1
          end if
          if( istart == isgList(tileMe) ) then  ! may need to send data to west of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, WEST, tile1, tile2, align1, align2, jend1, jend2, jegList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   jend = jend1(m)
                else
                   jend = jend2(m)
                end if
                if( jend == jegList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                   do l = 1, numS
                      if( tileSend(l) == tileCor .AND. alignSend(l) == NORTH_EAST ) then
                         found_corner = .false.
                         exit
                      end if
                   end do
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = NORTH_EAST
             isdSend(numS) = iegList(tileCor) + 1 + ishift;      iedSend(numS) = iegList(tileCor) + xhalo + ishift
             jsdSend(numS) = jegList(tileCor) + 1 + jshift;      jedSend(numS) = jegList(tileCor) + yhalo + jshift
             iscSend(numS) = isgList(tileMe) + ioff;             iecSend(numS) = isgList(tileMe) + ioff + xhalo - 1
             jscSend(numS) = jsgList(tileMe) + joff;             jecSend(numS) = jsgList(tileMe) + joff + yhalo - 1
          end if
       end if
    end do

    !------------------------------------------------------------------------------
    !                            northwest corner
    !------------------------------------------------------------------------------
    found_corner = .FALSE.
    ioff = 0; joff = 0
    !--- define the recv overlap for northwest corner.
    tileNbr = get_neighbor_tile( tileMe, WEST, tile1, tile2, align1, align2, jend1, jend2, jsgList, n )
    if( tileNbr > 0 ) then ! found nbr
       jend = jend2(n)
       if( jend == jegList(tileNbr) ) then
          tileCor = get_neighbor_tile( tileNbr, NORTH, tile1, tile2, align1, align2, iend1, iend2, iegList, n )
          if(tileCor > 0) found_corner = .true.
       else
          joff = jend - jsgList(tileNbr) + 1
          tileCor = tileNbr
          found_corner = .true.
       end if
    else ! not found nbr for west contact, then found nbr for south contact
       tileNbr = get_neighbor_tile( tileMe, NORTH, tile1, tile2, align1, align2, istart1, istart2, isgList, n)
       if( tileNbr > 0 ) then ! found nbr in north contact
          istart = istart2(n)
          if( istart == isgList(tileNbr) ) then
             tileCor = get_neighbor_tile( tileNbr, WEST, tile1, tile2, align1, align2, jstart1, jstart2, isgList, n )
             if(tileCor > 0) found_corner = .true.
          else
             ioff = iegList(tileNbr) - istart + 1
             tileCor = tileNbr
             found_corner = .true.
          end if
       end if
    end if

    if(found_corner) then
       numR = numR + 1
       tileRecv(numR) = tileCor
       isdRecv(numR) = isgList(tileMe) -xhalo;              iedRecv(numR) = iegList(tileMe) - 1
       jsdRecv(numR) = jegList(tileMe) + 1 + jshift;        jedRecv(numR) = jegList(tileMe) + yhalo + jshift
       iscRecv(numR) = iegList(tileCor) - ioff - xhalo + 1; iecRecv(numR) = iegList(tileCor) - ioff
       jscRecv(numR) = jsgList(tileCor) + joff;             jecRecv(numR) = jsgList(tileCor) + joff + xhalo - 1
    end if

    !--- define the send for to_pe's southeast corner.
    !--- loop over all the contact 
    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .false.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       jstart   = jstart1(n)
       jend     = jend1(n)
       if(alignNbr == WEST) then!--- to_pe's west contact
          if( jend < jegList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             ioff = 0; joff = jend - jsgList(tileMe) + 1
          end if
          if( jstart == jsgList(tileMe) ) then  ! may need to send data to south of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, SOUTH, tile1, tile2, align1, align2, istart1, istart2, isgList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   istart = istart1(m)
                else
                   istart = istart2(m)
                end if
                if( istart == isgList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = NORTH_WEST
             isdSend(numS) = isgList(tileCor) - xhalo;           iedSend(numS) = isgList(tileCor) - 1
             jsdSend(numS) = jegList(tileCor) + 1 + jshift;      jedSend(numS) = jegList(tileCor) + yhalo + jshift
             iscSend(numS) = iegList(tileMe) - ioff - xhalo + 1; iecSend(numS) = iegList(tileMe) - ioff
             jscSend(numS) = jsgList(tileMe) + joff;             jecSend(numS) = jsgList(tileMe) + joff + yhalo - 1
          end if

       end if
    end do

    do n = 1, num_contact
       if(tile1(n) .NE. tileMe .AND. tile2(n) .NE. tileMe) cycle
       found_corner = .FALSE.
       alignNbr = align2(n)
       tileNbr  = tile2(n)
       istart   = istart1(n)
       iend     = iend1(n)
       if(alignNbr == NORTH) then !--- to_pe's SOUTH contact
          if( istart > isgList(tileMe) ) then ! need to send data to this nbr
             found_corner = .true.
             tileCor = tileNbr
             joff = 0; ioff = iegList(tileMe) - istart + 1
          end if
          if( iend == iegList(tileMe) ) then  ! may need to send data to east of tileNbr.
             tileCor = get_neighbor_tile( tileNbr, EAST, tile1, tile2, align1, align2, jend1, jend2, jegList, m)
             if( tileCor > 0 ) then
                if( tile1(m) == tileCor ) then
                   jend = jend1(m)
                else
                   jend = jend2(m)
                end if
                if( jend == jegList(tileCor) ) then ! need to send data to tileCor
                   found_corner = .true.
                   ioff = 0; joff = 0
                   do l = 1, numS 
                      if( tileSend(l) == tileCor .AND. alignSend(l) == NORTH_WEST ) then
                         found_corner = .false.
                         exit
                      end if
                   end do
                end if
             end if
          end if
          if(found_corner) then
             numS = numS + 1
             tileSend(numS) = tileCor
             alignSend(numS) = NORTH_WEST
             isdSend(numS) = isgList(tileCor) - xhalo;           iedSend(numS) = isgList(tileCor) - 1
             jsdSend(numS) = jegList(tileCor) + 1 + jshift;      jedSend(numS) = jegList(tileCor) + yhalo + jshift
             iscSend(numS) = iegList(tileMe) - ioff - xhalo + 1; iecSend(numS) = iegList(tileMe) - ioff
             jscSend(numS) = jsgList(tileMe) + joff;             jecSend(numS) = jsgList(tileMe) + joff + yhalo - 1
          end if
       end if
    end do


  end subroutine compute_corner_overlap

  !#############################################################################
  !--- This function will find the neighbor tile
  function get_neighbor_tile( tileMe, alignment, tile1, tile2, align1, align2, iList1, iList2, iGlobal, index )
    integer,               intent(in)  :: tileMe, alignment 
    integer, dimension(:), intent(in)  :: tile1, tile2, align1, align2, iList1, iList2, iGlobal
    integer,               intent(out) :: index
    integer                            :: get_neighbor_tile
    integer                            :: n, num_contact

    get_neighbor_tile = 0
    index    = 0

    num_contact = size(tile1(:))

    do n = 1, num_contact
       if(tile1(n) == tileMe) then
          if(align1(n) == alignment .AND. iList1(n) == iGlobal(tileMe)) then  ! found 
             get_neighbor_tile = tile2(n)
             index = n
             exit
          end if
       else if(tile2(n) == tileMe) then
          if(align2(n) == alignment .AND. iList2(n) == iGlobal(tileMe)) then  ! found 
             get_neighbor_tile = tile1(n)
             index = n
             exit
          end if
       end if
    end do

    return

  end function get_neighbor_tile


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
          call mpp_error(FATAL, 'mpp_domains_define_mod: The contact region is not on the x-boundary of the tile')
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
          call mpp_error(FATAL, 'mpp_domains_define_mod: The contact region is not on the y-boundary of the tile')
       end if
       if ( is > ie ) then
          i = is; is = ie; ie = i
       end if
    else
       call mpp_error(FATAL, 'mpp_domains_define_mod: The contact region should be line contact' )
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
!   <IN NAME="halo" TYPE="integer,optional" > </IN>
!   <IN NAME="cbegin,cend" TYPE="integer,optional" > </IN>
!   <IN NAME="gbegin,gend" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain1D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain1D(domain_in,domain_out,cbegin,cend,gbegin,gend, halo)
      ! </PUBLICROUTINE>
      type(domain1D), intent(in)    :: domain_in
      type(domain1D), intent(inout) :: domain_out
      integer, intent(in), optional :: halo
      integer, intent(in), optional :: cbegin, cend             ! extent of compute_domain
      integer, intent(in), optional :: gbegin, gend             ! extent of global domain
      integer :: ndivs, global_indices(2) !(/ isg, ieg /)
      integer                       :: halosz, flag
! get the global indices of the input domain
      global_indices(1) = domain_in%global%begin;  global_indices(2) = domain_in%global%end

! get the halo
      halosz = 0
      if(present(halo)) halosz = halo

      ! get the layout
      ndivs = size(domain_in%list(:))

! get the flag
      flag = 0
      if(domain_in%cyclic) flag = flag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%data%is_global) flag = flag + GLOBAL_DATA_DOMAIN

      call mpp_define_domains( global_indices, ndivs, domain_out, pelist = domain_in%list(:)%pe, &
           flags = flag, halo = halosz, extent = domain_in%list(:)%compute%size )
           
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
!   <IN NAME="xhalo,yhalo" TYPE="integer,optional" > </IN>
!   <INOUT NAME="domain_out" TYPE="type(domain2D)" > </INOUT>

! <PUBLICROUTINE>
    subroutine mpp_modify_domain2D(domain_in, domain_out, isc, iec, jsc, jec, isg, ieg, jsg, jeg, xhalo, yhalo)
      ! </PUBLICROUTINE>
      type(domain2D), intent(in)    :: domain_in
      type(domain2D), intent(inout) :: domain_out
      integer, intent(in), optional :: isc, iec, jsc, jec
      integer, intent(in), optional :: isg, ieg, jsg, jeg
      integer, intent(in), optional :: xhalo, yhalo
      integer                       :: global_indices(4), layout(2)
      integer                       :: xhalosz, yhalosz, xflag, yflag
      
      if(present(xhalo) .or. present(yhalo)) then
! get the global indices of the input domain
         global_indices(1) = domain_in%x%global%begin;  global_indices(2) = domain_in%x%global%end
         global_indices(3) = domain_in%y%global%begin;  global_indices(4) = domain_in%y%global%end
         
! get the halo
      xhalosz = 0; yhalosz = 0
      if(present(xhalo)) xhalosz = xhalo
      if(present(yhalo)) yhalosz = yhalo
      
         ! get the layout
         layout(1) = size(domain_in%x%list(:)); layout(2) = size(domain_in%y%list(:))
         
! get the flag
      xflag = 0; yflag = 0
      if(domain_in%x%cyclic) xflag = xflag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%x%data%is_global) xflag = xflag + GLOBAL_DATA_DOMAIN
      if(domain_in%y%cyclic) yflag = yflag + CYCLIC_GLOBAL_DOMAIN
      if(domain_in%y%data%is_global) yflag = yflag + GLOBAL_DATA_DOMAIN
      
         call mpp_define_domains( global_indices, layout, domain_out, pelist = domain_in%list(:)%pe, &
                                  xflags = xflag, yflags = yflag,  xhalo = xhalosz, yhalo = yhalosz, &
                                  xextent = domain_in%x%list(:)%compute%size,                        &
                                  yextent = domain_in%y%list(:)%compute%size,                        &
                                  symmetry=domain_in%symmetry,                                       &
                                  maskmap = domain_in%pearray .NE. NULL_PE )
                               
      else    
  domain_out = NULL_DOMAIN2D
         call mpp_modify_domain(domain_in%x, domain_out%x, isc, iec, isg, ieg)
         call mpp_modify_domain(domain_in%y, domain_out%y, jsc, jec, jsg, jeg)
      endif
         
    end subroutine mpp_modify_domain2D
! </SUBROUTINE>

  !#####################################################################

