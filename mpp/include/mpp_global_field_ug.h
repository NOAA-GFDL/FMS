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
    subroutine MPP_GLOBAL_FIELD_UG_2D_( domain, local, global, flags, default_data)
      type(domainUG), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:)
      MPP_TYPE_, intent(out) :: global(:)
      integer, intent(in), optional :: flags
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),1)
      MPP_TYPE_ :: global3D(size(global,1),1)
      pointer( lptr,  local3D )
      pointer( gptr, global3D )
      lptr = LOC( local)
      gptr = LOC(global)
      call mpp_global_field_UG( domain, local3D, global3D, flags, default_data )

    end subroutine MPP_GLOBAL_FIELD_UG_2D_

    subroutine MPP_GLOBAL_FIELD_UG_3D_( domain, local, global, flags, default_data)
!get a global field from a local field
!local field may be on compute OR data domain
      type(domainUG), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(domain%compute%begin:,:)
      MPP_TYPE_, intent(out) :: global(domain%global%begin:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_, intent(in), optional :: default_data

     integer :: l, k, m, n, nd, nwords, lpos, rpos, ioff, joff, from_pe, tile_id
      integer :: ke, lsc, lec, ls, le, nword_me
      integer :: ipos, jpos
      logical :: xonly, yonly, root_only, global_on_this_pe
      MPP_TYPE_ :: clocal (domain%compute%size*size(local,2))
      MPP_TYPE_ :: cremote(domain%compute%max_size*size(local,2))
      integer :: stackuse
      character(len=8) :: text

      pointer( ptr_local,  clocal  )
      pointer( ptr_remote, cremote )

      stackuse = size(clocal(:))+size(cremote(:))
      if( stackuse.GT.mpp_domains_stack_size )then
          write( text, '(i8)' )stackuse
          call mpp_error( FATAL, &
               'MPP_DO_GLOBAL_FIELD_UG user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
      end if
      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, stackuse )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal(:))+1))

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_UG: must first call mpp_domains_init.' )

      root_only = .FALSE.
      if( PRESENT(flags) ) root_only = BTEST(flags, ROOT_GLOBAL)

      global_on_this_pe =  .NOT. root_only .OR. domain%pe == domain%tile_root_pe
      ipos = 0; jpos = 0
      if(global_on_this_pe ) then
         if(size(local,2).NE.size(global,2) ) call mpp_error( FATAL, &
              'MPP_GLOBAL_FIELD_UG: mismatch of third dimension size of global and local')
         if( size(global,1).NE.domain%global%size)then
            call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_UG: incoming arrays do not match domain.' )
         endif
      endif

      if( size(local,1).NE.domain%compute%size )then
         call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_UG_: incoming field array must match compute domain ' )
      end if

      ke  = size(local,2)
      lsc = domain%compute%begin; lec = domain%compute%end

      nword_me = (lec-lsc+1)*ke

! make contiguous array from compute domain
      m = 0
      if(global_on_this_pe) then
         if(PRESENT(default_data)) then
            global = default_data
         else
            global = 0
         endif

         do k = 1, ke
            do l = lsc, lec
               m = m + 1
               clocal(m) = local(l,k)
               global(l,k) = clocal(m) !always fill local domain directly
            end do
         end do
      else
         do k = 1, ke
            do l = lsc, lec
               m = m + 1
               clocal(m) = local(l,k)
            end do
         end do
      endif

      !fill off-domains (note loops begin at an offset of 1)
      tile_id = domain%tile_id
      nd = size(domain%list(:))
      if(root_only) then
         if(domain%pe .NE. domain%tile_root_pe) then
            call mpp_send( clocal(1), plen=nword_me, to_pe=domain%tile_root_pe, tag=COMM_TAG_1 )
         else
            do n = 1,nd-1
               rpos = mod(domain%pos+n,nd)
               if( domain%list(rpos)%tile_id .NE. tile_id ) cycle
               nwords = domain%list(rpos)%compute%size * ke
               call mpp_recv(cremote(1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_1 )
               m = 0
               ls = domain%list(rpos)%compute%begin; le = domain%list(rpos)%compute%end

               do k = 1,ke
                  do l = ls, le
                     m = m + 1
                     global(l,k) = cremote(m)
                  end do
               end do
            end do
         endif
      else
         do n = 1,nd-1
            lpos = mod(domain%pos+nd-n,nd)
            if( domain%list(lpos)%tile_id.NE. tile_id ) cycle ! global field only within tile
            call mpp_send( clocal(1), plen=nword_me, to_pe=domain%list(lpos)%pe, tag=COMM_TAG_2 )
         end do
         do n = 1,nd-1
            rpos = mod(domain%pos+n,nd)
            if( domain%list(rpos)%tile_id .NE. tile_id ) cycle ! global field only within tile
            nwords = domain%list(rpos)%compute%size * ke
            call mpp_recv( cremote(1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_2 )
            m = 0
            ls = domain%list(rpos)%compute%begin; le = domain%list(rpos)%compute%end

            do k = 1,ke
               do l = ls, le
                  m = m + 1
                  global(l,k) = cremote(m)
               end do
            end do
         enddo
      endif

      call mpp_sync_self()

    end subroutine MPP_GLOBAL_FIELD_UG_3D_


   subroutine MPP_GLOBAL_FIELD_UG_4D_( domain, local, global, flags, default_data )
      type(domainUG), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),size( local,2)*size( local,3))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2)*size(global,3))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field_UG( domain, local3D, global3D, flags, default_data )
    end subroutine MPP_GLOBAL_FIELD_UG_4D_

    subroutine MPP_GLOBAL_FIELD_UG_5D_( domain, local, global, flags, default_data )
      type(domainUG), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),size( local,2)*size( local,3)*size( local,4))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2)*size(global,3)*size(global,4))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field_UG( domain, local3D, global3D, flags, default_data )
    end subroutine MPP_GLOBAL_FIELD_UG_5D_



