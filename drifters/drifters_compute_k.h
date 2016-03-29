! -*-f90-*-

subroutine drifters_compute_k_XXX(self, positions, u, v, &
#if _DIMS >= 3
     & w, &
#endif
     & k, ermesg)

  use cloud_interpolator_mod
  type(drifters_type) :: self
  real, intent(in)    :: positions(:,:)
#if _DIMS == 2
  real, intent(in)    :: u(:,:)
  real, intent(in)    :: v(:,:)
#endif
#if _DIMS == 3
  real, intent(in)    :: u(:,:,:)
  real, intent(in)    :: v(:,:,:)
  real, intent(in)    :: w(:,:,:)
#endif
  real, intent(out)   :: k(:,:)
  character(len=*), intent(out) :: ermesg

  integer, parameter :: nd = _DIMS ! number of dims
  integer i, ip, np, ij(nd), ier, nsizes_u(nd), nsizes_v(nd)
#if _DIMS >= 3
  integer nsizes_w(nd)
#endif
  real fvals(2**nd), ts(nd)
  real pos(nd, self%core%np)

  ermesg = ''

  nsizes_u(1) = size(u, 1)
  nsizes_u(2) = size(u, 2)

  nsizes_v(1) = size(v, 1)
  nsizes_v(2) = size(v, 2)

#if _DIMS >= 3
  nsizes_u(3) = size(u, 3)
  nsizes_v(3) = size(v, 3)
  nsizes_w(1) = size(w, 1)
  nsizes_w(2) = size(w, 2)
  nsizes_w(3) = size(w, 3)
#endif

  np = self%core%np

  ! correct for periodicity
  if(self%comm%xperiodic) then
     do ip = 1, np
        pos(1,ip) = self%comm%xgmin + modulo(positions(1,ip)-self%comm%xgmin, self%comm%xgmax-self%comm%xgmin)
     enddo
  else
     pos(1,:) = positions(1,1:np)
  endif
  if(self%comm%yperiodic) then
     do ip = 1, np
        pos(2,ip) = self%comm%ygmin + modulo(positions(2,ip)-self%comm%ygmin, self%comm%ygmax-self%comm%ygmin)
     enddo
  else
     pos(2,:) = positions(2,1:np)
  endif

#if _DIMS >= 3
  pos(3,:) = positions(3,1:self%core%np)
#endif

  do ip = 1, np

     ! iterate over particles

     k(:, ip) = huge(1.)

     ! u-component...
     call cld_ntrp_locate_cell(self%xu, pos(1,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(1,ip), ' axis min/max=', minval(self%xu), maxval(self%xu)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xu(i))/(self%xu(i+1)-self%xu(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yu, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(2,ip), ' axis min/max=', minval(self%yu), maxval(self%yu)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yu(i))/(self%yu(i+1)-self%yu(i))
     ij(2) = i

#if _DIMS >= 3
     call cld_ntrp_locate_cell(self%zu, pos(3,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(3,ip), ' axis min/max=', minval(self%zu), maxval(self%zu)
     endif
#endif
     i = max(1, i)
     ts(3) = (pos(3,ip) - self%zu(i))/(self%zu(i+1)-self%zu(i))
     ij(3) = i
#endif

     call cld_ntrp_get_cell_values(nsizes_u, _FLATTEN(u), ij, fvals, ier)
     call cld_ntrp_linear_cell_interp(fvals, ts, k(1, ip), ier)
     k(1, ip) = self%dt * k(1, ip)

     ! v-component...
     call cld_ntrp_locate_cell(self%xv, pos(1,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(1,ip), ' axis min/max=', minval(self%xv), maxval(self%xv)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xv(i))/(self%xv(i+1)-self%xv(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yv, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(2,ip), ' axis min/max=', minval(self%yv), maxval(self%yv)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yv(i))/(self%yv(i+1)-self%yv(i))
     ij(2) = i

#if _DIMS >= 3
     call cld_ntrp_locate_cell(self%zv, pos(3,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(3,ip), ' axis min/max=', minval(self%zv), maxval(self%zv)
     endif
#endif
     i = max(1, i)
     ts(3) = (pos(3,ip) - self%zv(i))/(self%zv(i+1)-self%zv(i))
     ij(3) = i
#endif

     call cld_ntrp_get_cell_values(nsizes_v, _FLATTEN(v), ij, fvals, ier)
     call cld_ntrp_linear_cell_interp(fvals, ts, k(2, ip), ier)
     k(2, ip) = self%dt * k(2, ip)


#if _DIMS >= 3
     ! w-component...
     call cld_ntrp_locate_cell(self%xw, pos(1,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(1,ip), ' axis min/max=', minval(self%xw), maxval(self%xw)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xw(i))/(self%xw(i+1)-self%xw(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yw, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(2,ip), ' axis min/max=', minval(self%yw), maxval(self%yw)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yw(i))/(self%yw(i+1)-self%yw(i))
     ij(2) = i

     call cld_ntrp_locate_cell(self%zw, pos(3,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        print *,'***PE: ', _MPP_PE,' i=', i, 'pos=', pos(3,ip), ' axis min/max=', minval(self%zw), maxval(self%zw)
     endif
#endif
     i = max(1, i)
     ts(3) = (pos(3,ip) - self%zw(i))/(self%zw(i+1)-self%zw(i))
     ij(3) = i

     call cld_ntrp_get_cell_values(nsizes_w, _FLATTEN(w), ij, fvals, ier)
     call cld_ntrp_linear_cell_interp(fvals, ts, k(3, ip), ier)
     k(3, ip) = self%dt * k(3, ip)
#endif

  enddo

end subroutine drifters_compute_k_XXX

