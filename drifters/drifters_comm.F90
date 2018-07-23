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
#include "fms_switches.h"

module drifters_comm_mod
#include <fms_platform.h>

#ifdef _SERIAL

#define _TYPE_DOMAIN2D integer
#define _NULL_PE 0

#else

  use mpp_mod,         only        : NULL_PE, FATAL, NOTE, mpp_error, mpp_pe, mpp_npes
  use mpp_mod,         only        : mpp_root_pe
  use mpp_mod,         only        : mpp_send, mpp_recv, mpp_sync_self
  use mpp_mod,         only        : COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
  use mpp_domains_mod, only        : domain2D
  use mpp_domains_mod, only        : mpp_get_neighbor_pe, mpp_define_domains, mpp_get_layout
  use mpp_domains_mod, only        : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only        : NORTH, SOUTH, EAST, WEST, CYCLIC_GLOBAL_DOMAIN
  use mpp_domains_mod, only        : NORTH_EAST, SOUTH_EAST, SOUTH_WEST, NORTH_WEST

#define _TYPE_DOMAIN2D type(domain2d)
#define _NULL_PE NULL_PE

#endif

  use drifters_core_mod, only: drifters_core_type, drifters_core_remove_and_add,  drifters_core_set_positions

  implicit none
  private

  public :: drifters_comm_type, drifters_comm_new, drifters_comm_del, drifters_comm_set_pe_neighbors
  public :: drifters_comm_set_domain, drifters_comm_update, drifters_comm_gather

  type drifters_comm_type
     ! compute domain
     real           :: xcmin, xcmax
     real           :: ycmin, ycmax
     ! data domain
     real           :: xdmin, xdmax
     real           :: ydmin, ydmax
     ! global valid min/max
     real           :: xgmin, xgmax
     real           :: ygmin, ygmax
     ! x/y period (can be be nearly infinite)
     logical        :: xperiodic, yperiodic
     ! neighbor domains
     integer        :: pe_N, pe_S, pe_E, pe_W, pe_NE, pe_SE, pe_SW, pe_NW
     ! starting/ending pe, set this to a value /= 0 if running concurrently
     integer        :: pe_beg, pe_end
  end type drifters_comm_type

contains

!===============================================================================
  subroutine drifters_comm_new(self)
    type(drifters_comm_type)   :: self
    
    self%xcmin = -huge(1.); self%xcmax = +huge(1.)
    self%ycmin = -huge(1.); self%ycmax = +huge(1.)

    self%xdmin = -huge(1.); self%xdmax = +huge(1.)
    self%ydmin = -huge(1.); self%ydmax = +huge(1.)

    self%xgmin = -huge(1.); self%xgmax = +huge(1.)
    self%ygmin = -huge(1.); self%ygmax = +huge(1.)

    self%xperiodic = .FALSE.; self%yperiodic = .FALSE.

    self%pe_N  = _NULL_PE
    self%pe_S  = _NULL_PE
    self%pe_E  = _NULL_PE
    self%pe_W  = _NULL_PE
    self%pe_NE = _NULL_PE
    self%pe_SE = _NULL_PE
    self%pe_SW = _NULL_PE
    self%pe_NW = _NULL_PE

    self%pe_beg =  0
    self%pe_end = -1
    
    
  end subroutine drifters_comm_new

!===============================================================================
  subroutine drifters_comm_del(self)
    type(drifters_comm_type)   :: self
    
    ! nothing to deallocate
    call drifters_comm_new(self)

  end subroutine drifters_comm_del

!===============================================================================
  subroutine drifters_comm_set_data_bounds(self, xmin, ymin, xmax, ymax)
    ! Set data domain bounds.
    type(drifters_comm_type)   :: self
    real, intent(in)           :: xmin, ymin, xmax, ymax
    
    self%xdmin = max(xmin, self%xdmin)
    self%xdmax = min(xmax, self%xdmax)
    self%ydmin = max(ymin, self%ydmin)
    self%ydmax = min(ymax, self%ydmax)

  end subroutine drifters_comm_set_data_bounds

!===============================================================================
  subroutine drifters_comm_set_comp_bounds(self, xmin, ymin, xmax, ymax)
    ! Set compute domain bounds.
    type(drifters_comm_type)   :: self
    real, intent(in)           :: xmin, ymin, xmax, ymax
    
    self%xcmin = max(xmin, self%xcmin)
    self%xcmax = min(xmax, self%xcmax)
    self%ycmin = max(ymin, self%ycmin)
    self%ycmax = min(ymax, self%ycmax)

  end subroutine drifters_comm_set_comp_bounds

!===============================================================================
  subroutine drifters_comm_set_pe_neighbors(self, domain)
    ! Set neighboring pe numbers.
    type(drifters_comm_type)   :: self
    _TYPE_DOMAIN2D, intent(inout) :: domain

#ifndef _SERIAL
! parallel code

    integer        :: pe_N, pe_S, pe_E, pe_W, pe_NE, pe_SE, pe_SW, pe_NW

    call mpp_get_neighbor_pe(domain, NORTH     , pe_N )
    call mpp_get_neighbor_pe(domain, NORTH_EAST, pe_NE)
    call mpp_get_neighbor_pe(domain, EAST      , pe_E )
    call mpp_get_neighbor_pe(domain, SOUTH_EAST, pe_SE)
    call mpp_get_neighbor_pe(domain, SOUTH     , pe_S )
    call mpp_get_neighbor_pe(domain, SOUTH_WEST, pe_SW)
    call mpp_get_neighbor_pe(domain, WEST      , pe_W )
    call mpp_get_neighbor_pe(domain, NORTH_WEST, pe_NW)

    if(pe_N  /= self%pe_N  .and. self%pe_N  == _NULL_PE) then
       self%pe_N  = pe_N 
    else if(pe_N  /= self%pe_N ) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: NORTH PE changed!.')
    endif 
    if(pe_NE /= self%pe_NE .and. self%pe_NE == _NULL_PE) then
       self%pe_NE = pe_NE
    else if(pe_NE /= self%pe_NE) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: NORTH-EAST PE changed!.')
    endif 
    if(pe_E  /= self%pe_E  .and. self%pe_E  == _NULL_PE) then
       self%pe_E  = pe_E 
    else if(pe_E  /= self%pe_E ) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: EAST PE changed!.')
    endif 
    if(pe_SE /= self%pe_SE .and. self%pe_SE == _NULL_PE) then
       self%pe_SE = pe_SE
    else if(pe_SE /= self%pe_SE) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: SOUTH-EAST PE changed!.')
    endif 
    if(pe_S  /= self%pe_S  .and. self%pe_S  == _NULL_PE) then
       self%pe_S  = pe_S 
    else if(pe_S  /= self%pe_S ) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: SOUTH PE changed!.')
    endif 
    if(pe_SW /= self%pe_SW .and. self%pe_SW == _NULL_PE) then
       self%pe_SW = pe_SW
    else if(pe_SW /= self%pe_SW) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: SOUTH-WEST PE changed!.')
    endif 
    if(pe_W  /= self%pe_W  .and. self%pe_W  == _NULL_PE) then
       self%pe_W  = pe_W 
    else if(pe_W  /= self%pe_W ) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: WEST PE changed!.')
    endif 
    if(pe_NW /= self%pe_NW .and. self%pe_NW == _NULL_PE) then
       self%pe_NW = pe_NW
    else if(pe_NW /= self%pe_NW) then
       call mpp_error( FATAL, 'drifters_comm_set_pe_neighbors: NORTH-WEST PE changed!.')
    endif 

#endif
! end of parallel code

  end subroutine drifters_comm_set_pe_neighbors

!===============================================================================
  subroutine drifters_comm_set_domain(self, domain, x, y, backoff_x, backoff_y)
    ! Set boundaries of domain and compute neighbors. This method can be called
    ! multiple times; the data domain will just be the intersection (overlap) of
    ! all domains (e.g domain_u, domain_v, etc). 
    type(drifters_comm_type)   :: self
    _TYPE_DOMAIN2D, intent(inout) :: domain
    real, intent(in)           :: x(:), y(:)           ! global axes
    integer, intent(in)        :: backoff_x, backoff_y ! >=0, data domain is reduced by "backoff_x,y" indices in x, resp. y
    
    ! compute/data domain start/end indices
    integer isc, iec, jsc, jec
    integer isd, ied, jsd, jed
    integer nx, ny, hx, hy, bckf_x, bckf_y, halox, haloy
    real dx, dy, xdmin, xdmax, ydmin, ydmax

#ifdef _SERIAL
    integer :: ibnds(1)

    ibnds = lbound(x); isc = ibnds(1)
    ibnds = ubound(x); iec = ibnds(1) - 1
    ibnds = lbound(y); jsc = ibnds(1)
    ibnds = ubound(y); jec = ibnds(1) - 1
#else
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
#endif

    self%xcmin = max(x(isc), self%xcmin)
    self%xcmax = min(x(iec), self%xcmax)
    self%ycmin = max(y(jsc), self%ycmin)
    self%ycmax = min(y(jec), self%ycmax)

    nx = iec - isc + 1
    ny = jec - jsc + 1

#ifdef _SERIAL
    isd = 1; ied = size(x); jsd = 1; jed = size(y)
#else
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
#endif

    hx = max(ied-iec, isc-isd)
    hy = max(jed-jec, jsc-jsd)
    bckf_x = min(backoff_x, hx)
    bckf_y = min(backoff_y, hy)

    halox = max(0, hx - bckf_x)
    haloy = max(0, hy - bckf_y)

    if(isd < 1) then
       dx = x(2) - x(1)
       xdmin = self%xcmin - dx*halox
    else
       xdmin = x(isd+bckf_x)
    endif

    if(ied > nx) then
       dx = x(nx) - x(nx-1)
       xdmax = self%xcmax + dx*halox
    else
       xdmax = x(ied-bckf_x)
    endif

    if(jsd < 1) then
       dy = y(2) - y(1)
       ydmin = self%ycmin - dy*haloy
    else
       ydmin = y(jsd+bckf_y)
    endif

    if(jed > ny) then
       dy = y(ny) - y(ny-1)
       ydmax = self%ycmax + dy*haloy
    else
       ydmax = y(jed-bckf_y)
    endif
    
    self%xdmin = max(xdmin, self%xdmin)
    self%ydmin = max(ydmin, self%ydmin)
    self%xdmax = min(xdmax, self%xdmax)
    self%ydmax = min(ydmax, self%ydmax)

    call drifters_comm_set_pe_neighbors(self, domain)

  end subroutine drifters_comm_set_domain

!===============================================================================
  subroutine drifters_comm_update(self, drfts, new_positions, &
       & comm, remove, max_add_remove)

    type(drifters_comm_type)   :: self
    type(drifters_core_type)   :: drfts
    real, intent(inout)           :: new_positions(:,:)
    integer, intent(in), optional :: comm ! MPI communicator
    logical, intent(in), optional :: remove(:) ! Set to True for particles that should be removed
    integer, intent(in), optional :: max_add_remove ! max no of particles to add/remove

#ifdef _SERIAL
! serial code
    
    drfts%positions(:, 1:drfts%np) = new_positions(:, 1:drfts%np)
    return

#else
! parallel code

    include 'mpif.h'


    integer nd, np, nar_est, ip, neigh_pe, irem, pe, npes, ntuples
    integer ntuples_tot, ndata, mycomm
#ifdef _USE_MPI
    integer ier 
#endif
    integer, allocatable :: iadd(:)
    integer, allocatable :: table_recv(:), table_send(:)
    real   , allocatable :: data_recv(:,:), data_send(:,:)
    integer, allocatable :: indices_to_remove(:)
    integer, allocatable :: ids_to_add(:)
    real   , allocatable :: positions_to_add(:,:)
    real                 :: x, y, xold, yold
    character(len=128) :: ermsg, notemsg
    logical            :: is_present
    integer            :: id, j, k, m, n, el
    logical            :: crossed_W, crossed_E, crossed_S, crossed_N
    logical            :: was_in_compute_domain, left_domain

    mycomm = MPI_COMM_WORLD
    if( present(comm) ) mycomm = comm

    nd = drfts%nd
    np = size(new_positions,2)
    if(np > 0 .and. nd < 2) call mpp_error( FATAL, &
         & 'drifters_comm_update: number of dimensions must be 2 or higher.' )

    nar_est = 100 
    if(present(max_add_remove)) nar_est = max(1, max_add_remove)

    pe   = mpp_pe()
    npes = mpp_npes()

    ! assume pe list is contiguous, self%pe_beg...self%pe_end
    allocate(iadd(self%pe_beg:self%pe_end))
    allocate(table_recv(self%pe_beg:self%pe_end))
    allocate(table_send(self%pe_beg:self%pe_end))
    allocate(data_recv(nar_est*(1+nd), self%pe_beg:self%pe_end))
    allocate(data_send(nar_est*(1+nd), self%pe_beg:self%pe_end))
    allocate(indices_to_remove(nar_est))

    table_send = 0
    table_recv = 0
    data_send  = 0
    data_recv  = 0

    iadd = 0
    irem = 0
    do ip = 1, np
       x = new_positions(1, ip)
       y = new_positions(2, ip)
       xold = drfts%positions(1, ip)
       yold = drfts%positions(2, ip)

       if(    xold<self%xcmin .or. xold>self%xcmax .or. &
            & yold<self%ycmin .or. yold>self%ycmax      ) then
          was_in_compute_domain = .FALSE.
       else
          was_in_compute_domain = .TRUE.
       endif

       ! check if drifters crossed compute domain boundary 
       
       crossed_W = .FALSE.
       crossed_E = .FALSE.
       crossed_S = .FALSE.
       crossed_N = .FALSE.
       if( was_in_compute_domain .and. &
            & (x<self%xcmin) .and. (xold>self%xcmin) ) crossed_W = .TRUE.
       if( was_in_compute_domain .and. &
            & (x>self%xcmax) .and. (xold<self%xcmax) ) crossed_E = .TRUE.
       if( was_in_compute_domain .and. &
            & (y<self%ycmin) .and. (yold>self%ycmin) ) crossed_S = .TRUE.
       if( was_in_compute_domain .and. &
            & (y>self%ycmax) .and. (yold<self%ycmax) ) crossed_N = .TRUE.

       neigh_pe = _NULL_PE
       if(crossed_N .and. .not. crossed_E .and. .not. crossed_W) neigh_pe = self%pe_N
       if(crossed_N .and.       crossed_E                      ) neigh_pe = self%pe_NE
       if(crossed_E .and. .not. crossed_N .and. .not. crossed_S) neigh_pe = self%pe_E
       if(crossed_S .and.       crossed_E                      ) neigh_pe = self%pe_SE
       if(crossed_S .and. .not. crossed_E .and. .not. crossed_W) neigh_pe = self%pe_S
       if(crossed_S .and.       crossed_W                      ) neigh_pe = self%pe_SW
       if(crossed_W .and. .not. crossed_S .and. .not. crossed_N) neigh_pe = self%pe_W
       if(crossed_N .and.       crossed_W                      ) neigh_pe = self%pe_NW

       if(neigh_pe /= _NULL_PE) then
          iadd(neigh_pe) = iadd(neigh_pe) + 1
          if(iadd(neigh_pe) > nar_est) then
             write(notemsg, '(a,i4,a,i4,a)') 'drifters_comm_update: exceeded nar_est (', &
                  & iadd(neigh_pe),'>',nar_est,').'
             call mpp_error( FATAL, notemsg)
          endif
          table_send(neigh_pe)  = table_send(neigh_pe) + 1
          k = ( iadd(neigh_pe)-1 )*(1+nd) + 1
          data_send(k       , neigh_pe) = drfts%ids(ip)
          data_send(k+1:k+nd, neigh_pe) = new_positions(:,ip)
       endif

       ! check if drifters left data domain

       left_domain = .FALSE.
       if(       (x<self%xdmin .and. (self%pe_W/=pe)) .or. &
            &    (x>self%xdmax .and. (self%pe_E/=pe)) .or. &
            &    (y<self%ydmin .and. (self%pe_S/=pe)) .or. &
            &    (y>self%ydmax .and. (self%pe_N/=pe)) ) then
          left_domain = .TRUE.
       endif

       ! remove if particle was tagged as such

       if(present(remove)) then
          if(remove(ip)) left_domain = .TRUE.
       endif

       if(left_domain) then
          irem = irem + 1
          if(irem > nar_est) then
             write(notemsg, '(a,i4,a,i4,a)') 'drifters_comm_update: exceeded nar_est (',&
                  & irem,'>',nar_est,').'
             call mpp_error( FATAL, notemsg)
          endif
          indices_to_remove(irem) = ip
       endif

    enddo


    ! update drifters' positions (remove whatever needs to be removed later)
    call drifters_core_set_positions(drfts, new_positions, ermsg)
    if(ermsg/='') call mpp_error( FATAL, ermsg)

    ! fill in table_recv from table_send. table_send contains the
    ! number of tuples that will be sent to another pe. table_recv
    ! will contain the number of tuples to be received. The indices 
    ! of table_send refer to the pe where the tuples should be sent to;
    ! the indices of table_recv refer to the pe number 
    ! (self%pe_beg..self%pe_end) from
    ! which the tuple should be received from.
    !
    ! table_send(to_pe) = ntuples; table_recv(from_pe) = ntuples

    ! the following is a transpose operation
    ! table_send(m)[pe] -> table_recv(pe)[m]
    do m = self%pe_beg, self%pe_end
#ifdef _USE_MPI
       call MPI_Scatter (table_send   , 1, MPI_INTEGER,  &
            &            table_recv(m), 1, MPI_INTEGER,  &
            &            m, mycomm, ier )
#else
       if(pe==m) then
          do k = self%pe_beg, self%pe_end
             call mpp_send(table_send(k), plen=1, to_pe=k, tag=COMM_TAG_1)
          enddo
       endif
       call mpp_recv(table_recv(m), glen=1, from_pe=m, tag=COMM_TAG_1)    
#endif
    enddo

    ! communicate new positions. data_send is an array of size n*(nd+1) times npes.
    ! Each column j of data_send contains the tuple (id, x, y, ..) to be sent to pe=j.
    ! Inversely, data_recv's column j contains tuples (id, x, y,..) received from pe=j.
    do m = self%pe_beg, self%pe_end
       ntuples = table_send(m)
       ndata   = ntuples*(nd+1)
       ! should be able to send ndata?
#ifdef _USE_MPI
       call MPI_Scatter (data_send     , nar_est*(1+nd), MPI_REAL8, &
            &            data_recv(1,m), nar_est*(1+nd), MPI_REAL8, &
            &            m, mycomm, ier )
#else
       if(pe==m) then
          do k = self%pe_beg, self%pe_end
             call mpp_send(data_send(1,k), plen=nar_est*(1+nd), to_pe=k, tag=COMM_TAG_2)
          enddo
       endif
       call mpp_recv(data_recv(1,m), glen=nar_est*(1+nd), from_pe=m, tag=COMM_TAG_2)           
#endif
    enddo

    ! total number of tuples will determine size of ids_to_add/positions_to_add
    ntuples_tot = 0
    do m = self%pe_beg, self%pe_end
       ntuples_tot = ntuples_tot + table_recv(m)
    enddo

    allocate(positions_to_add(nd, ntuples_tot))
    allocate(      ids_to_add(    ntuples_tot))

    ! fill positions_to_add and ids_to_add.
    k = 0
    do m = self%pe_beg, self%pe_end
       ! get ids/positions coming from all pes
       do n = 1, table_recv(m)
          ! iterate over all ids/positions coming from pe=m
          el = (n-1)*(nd+1) + 1
          id = int(data_recv(el, m))
          ! only add if id not already present in drfts
          ! this can happen if a drifter meanders about 
          ! the compute domain boundary
          is_present = .false.
          do j = 1, drfts%np
             if(id == drfts%ids(j)) then
                is_present = .true.
                write(notemsg, '(a,i4,a)') 'Drifter ', id, ' already advected (will not be added).'
                call mpp_error(NOTE, notemsg)
                exit
             endif
          enddo
          if(.not. is_present) then
             k = k + 1
             ids_to_add(k)         = id

             positions_to_add(1:nd, k) = data_recv(el+1:el+nd, m)

          endif
       enddo
    enddo
    
    ! remove and add
    if(irem > 0 .or. k > 0) then
       write(notemsg, '(i4,a,i4,a)') irem, ' drifter(s) will be removed, ', k,' will be added'
       call mpp_error(NOTE, notemsg)
!!$       if(k>0) print *,'positions to add ', positions_to_add(:,1:k)
!!$       if(irem>0) print *,'ids to remove: ', indices_to_remove(1:irem)
    endif
    call drifters_core_remove_and_add(drfts, indices_to_remove(1:irem), &
         & ids_to_add(1:k), positions_to_add(:,1:k), ermsg)
    if(ermsg/='') call mpp_error( FATAL, ermsg)

#ifndef _USE_MPI
    ! make sure unbuffered mpp_isend call returned before deallocating
    call mpp_sync_self()
#endif

    deallocate(ids_to_add)
    deallocate(positions_to_add)

    deallocate(iadd)
    deallocate(table_recv)
    deallocate(table_send)
    deallocate(data_recv)
    deallocate(data_send)
    deallocate(indices_to_remove)

#endif
! end of parallel code

  end subroutine drifters_comm_update

!===============================================================================
  subroutine drifters_comm_gather(self, drfts, dinp, &
       & lons, lats, do_save_lonlat, &
       & filename, &
       & root, mycomm)

    use drifters_input_mod, only : drifters_input_type, drifters_input_save

    type(drifters_comm_type)   :: self
    type(drifters_core_type)   :: drfts
    type(drifters_input_type)  :: dinp
    real, intent(in)           :: lons(:), lats(:)
    logical, intent(in)        :: do_save_lonlat
    character(len=*), intent(in)  :: filename
    integer, intent(in), optional :: root    ! root pe
    integer, intent(in), optional :: mycomm  ! MPI communicator

    character(len=128) :: ermesg

#ifdef _SERIAL
! serial code

       dinp%ids(1:drfts%np)          = drfts%ids(1:drfts%np)
       dinp%positions(:, 1:drfts%np) = drfts%positions(:, 1:drfts%np)

       if(do_save_lonlat) then

          call drifters_input_save(dinp, filename=filename, &
               & geolon=lons, geolat=lats, ermesg=ermesg)

       else

          call drifters_input_save(dinp, filename=filename, ermesg=ermesg)

       endif

#else
! parallel code


    integer :: npf, ip, comm, root_pe, pe, npes, nd, np, npdim, npmax, ier, nptot
    integer :: i, j, k, kk
    integer, allocatable ::  nps(:)
    real    :: x, y
    real, allocatable :: lons0(:), lats0(:), recvbuf(:,:)
    real    :: data(drfts%nd+3, drfts%np)
    include 'mpif.h'

    comm    = MPI_COMM_WORLD
    if(present(mycomm)) comm = mycomm

    root_pe = mpp_root_pe()
    if(present(root)) root_pe = root

    pe   = mpp_pe()
    npes = mpp_npes()

    nd = drfts%nd
    np = drfts%np
    npdim = drfts%npdim
    
    allocate(nps(self%pe_beg:self%pe_end))
    nps = 0
    
    ! npf= number of drifters in compute domain

    npf = 0
    do ip = 1, np
       x = drfts%positions(1, ip)
       y = drfts%positions(2, ip)
       if( x <= self%xcmax .and. x >= self%xcmin .and. &
        &  y <= self%ycmax .and. y >= self%ycmin) then
          npf = npf + 1
          data(1       , npf)   = real(drfts%ids(ip))
          data(1+1:1+nd, npf)   =      drfts%positions(:, ip)
          data(    2+nd, npf)   = lons(ip)
          data(    3+nd, npf)   = lats(ip)
       endif
    enddo

    ! gather number of drifters
#ifdef _USE_MPI
    call mpi_gather(npf, 1, MPI_INT, &
         &          nps, 1, MPI_INT, &
         &          root_pe, comm, ier)
    !!if(ier/=0) ermesg = 'drifters_write_restart: ERROR while gathering "npf"'
#else
    call mpp_send(npf, plen=1, to_pe=root_pe, tag=COMM_TAG_3)
    if(pe==root_pe) then
       do i = self%pe_beg, self%pe_end
          call mpp_recv(nps(i), glen=1, from_pe=i, tag=COMM_TAG_3)
       enddo
    endif
#endif
    
    ! Now we know the max number of drifters to expect from each PE, so allocate 
    ! recvbuf (first dim will be zero on all PEs except root).

    ! allocate recvbuf to receive all the data on root PE, strided by npmax*(nd+3)
    npmax = maxval(nps)
    allocate(recvbuf(npmax*(nd+3), self%pe_beg:self%pe_end))
    recvbuf = -1

    ! Each PE sends data to recvbuf on root_pe.
#ifdef _USE_MPI
    call mpi_gather(         data  ,     npf*(nd+3), MPI_REAL8, &
         &                  recvbuf,   npmax*(nd+3), MPI_REAL8, &
         &          root_pe, comm, ier)
    !!if(ier/=0) ermesg = 'drifters_write_restart: ERROR while gathering "data"'
#else
    if(npf > 0) call mpp_send(data(1,1), plen=npf*(nd+3), to_pe=root_pe, tag=COMM_TAG_4)
    if(pe==root_pe) then
       do i = self%pe_beg, self%pe_end
          if(nps(i) > 0) call mpp_recv(recvbuf(1, i), glen=nps(i)*(nd+3), from_pe=i, tag=COMM_TAG_4)
       enddo
    endif
#endif
    
    ! Set positions and ids
    if(pe == root_pe) then
       ! check dims 
       nptot = sum(nps) ! total number of drifters, across al PEs
       if(nptot /= size(dinp%ids)) then
          deallocate(dinp%ids      , stat=ier)
          deallocate(dinp%positions, stat=ier)
          allocate(dinp%ids(nptot))
          allocate(dinp%positions(nd, nptot))
          dinp%ids       = -1
          dinp%positions = -huge(1.)
       endif

       allocate(lons0(nptot), lats0(nptot))

       ! Set the new positions/ids in dinp, on PE=root. Also set 
       ! lons/lats, these arrays will hold garbage if x1, y1, etc. were 
       ! not passed as subroutine arguments, that's ok 'cause we won't
       ! save them.
       j = 1
       do i = self%pe_beg, self%pe_end
          do k = 1, nps(i)
             kk = (nd+3)*(k-1)
             dinp%ids(j)             = int(recvbuf(kk+1          , i))
             dinp%positions(1:nd, j) =     recvbuf(kk+1+1:kk+1+nd, i)
             lons0(j)                =     recvbuf(kk+2+nd, i)
             lats0(j)                =     recvbuf(kk+3+nd, i)
             j = j + 1
          enddo
       enddo

       if(do_save_lonlat) then

          call drifters_input_save(dinp, filename=filename, &
               & geolon=lons0, geolat=lats0, ermesg=ermesg)

       else

          call drifters_input_save(dinp, filename=filename, ermesg=ermesg)

       endif

       deallocate(lons0, lats0)

    endif

#ifndef _USE_MPI
    call mpp_sync_self()
#endif
    deallocate(nps    , stat=ier)
    deallocate(recvbuf, stat=ier)

#endif
! _end of parallel code

  end subroutine drifters_comm_gather


end module drifters_comm_mod

!===============================================================================
!===============================================================================
#ifdef _TEST_DRIFTERS_COMM
! set FMS=/home/ap/ia64/mpp/mpp_test/exec/
! set OPTS="-r8 -g"
! set OBJS="$FMS/mpp*.o $FMS/threadloc.o"
! set INCS="-I/usr/include -I/usr/local/include -I${FMS}"
! set LIBS="-L/usr/local/lib -lnetcdf -L/usr/lib -lmpi -lcprts"
! ifort $OPTS $INCS -D_TEST_DRIFTERS_COMM drifters_comm.F90 quicksort.F90 drifters_core.F90 $OBJS $LIBS
program main

  use drifters_core_mod
  use drifters_comm_mod
  use mpp_mod
  use mpp_domains_mod

  implicit none

  integer, parameter :: nd=2, npmax = 4
  integer :: nx, ny, halox, haloy, layout(2), i, j, npes, pe, it, nt
  _TYPE_DOMAIN2D      :: domain
  type(drifters_core_type)      :: drfts
  type(drifters_comm_type) :: drfts_com
  real, parameter                     :: xmin=0., xmax=1., ymin=0., ymax=1.
  real                                :: dx, dy, u0, v0, dt, Lx, Ly
  real, allocatable                   :: x(:), y(:)
  character(len=128)  :: ermsg = ''
  integer :: ids(npmax)
  real    :: positions(nd, npmax), velocity(nd, npmax)
  integer :: io_status
!!$  integer :: stackmax=4000000

  namelist /drifters_comm_nml/ nx, ny, halox, haloy, u0, v0, dt, nt 
  

  call mpp_init !(MPP_DEBUG)
  !call mpp_set_stack_size(3145746)

  ! default input values
  nx = 11
  ny = 21
  halox = 2
  haloy = 2
  u0    = 1.0
  v0    = 0.0
  dt    = 0.1
  nt    = 10

  ! read input
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, drifters_comm_nml, iostat=io_status)
#else
  open(unit=1, file='input.nml', form='formatted')
  read(1, drifters_comm_nml)
  close(unit=1)
  if(mpp_pe()==0) write(*,drifters_comm_nml)
#endif

  ! create global domain
  Lx = xmax - xmin
  Ly = ymax - ymin
  dx = Lx/real(nx-1)
  dy = Ly/real(ny-1)
  allocate(x(nx), y(ny))
  x = xmin + (/ ( i*dx, i=0, nx-1) /)
  y = ymin + (/ ( j*dy, j=0, ny-1) /)

  ! decompose domain
  call mpp_domains_init ! (MPP_DEBUG)
!!$  call mpp_domains_set_stack_size(stackmax)

  npes = mpp_npes()
  call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
  if(mpp_pe()==0) print *,'LAYOUT: ', layout
  call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halox, yhalo=haloy,&
       & xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)

  ! set up drifters' communicator
  call drifters_comm_new(drfts_com)
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)
  ! repeated calls just for the fun of it
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)
  
  ! create drifters
  call drifters_core_new(drfts, nd, npmax, ermsg)
  if(ermsg /= '') print *,ermsg

  pe = mpp_pe()
  ids = (/ (i+100*pe, i=1,npmax) /)
  call drifters_core_set_ids(drfts, ids, ermsg)
  if(ermsg /= '') print *,ermsg

  ! position particles
  if(pe == 0) then
     positions(:, 1) = (/ (drfts_com%xcmin + drfts_com%xcmax)/2., &
          &               (drfts_com%ycmin + drfts_com%ycmax)/2. /)
     !positions(:, 2) = (/0.,0.01/)
     call drifters_core_set_positions(drfts, positions(:, 1:1), ermsg)
     if(ermsg /= '') print *,ermsg
  endif

  ! push drifters
  velocity(:,1) = (/u0, v0/)
  do it = 1, nt
     positions(:,1:drfts%np) = xmin + &
          & modulo(drfts%positions(:,1:drfts%np) + dt*velocity(:,1:drfts%np)-xmin, xmax-xmin)
     ! this will redistribute the drifters and update the positions
     call drifters_comm_update(drfts_com, drfts, positions(:,1:drfts%np))
     
     if(drfts%np > 0) then 
        do i=1,drfts%np
           print '(a,i6,a,i3,a,i3,a, i3, a,2f10.6)', 'PE: ',pe, ' it=', it, ' np=', drfts%np, ' ip=', i, &
                & ' x,y=', drfts%positions(1,i), drfts%positions(2,i)
        enddo
     endif
!!$     call drifters_print(drfts, ermsg)
!!$     if(ermsg /= '') print *,ermsg
  enddo

  ! clean up
  call drifters_core_del(drfts, ermsg)
  if(ermsg /= '') print *,ermsg
  call drifters_comm_del(drfts_com)
  call mpp_domains_exit
  call mpp_exit

end program main
#endif
! _TEST_DRIFTERS_COMM
