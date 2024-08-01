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
!> @ingroup time_interp
!> @addtogroup time_interp_external2_mod
!> @{

    !> @brief 2D interpolation for @ref time_interp_external
    subroutine TIME_INTERP_EXTERNAL_2D_(index, time, data_in, interp, verbose,horz_interp, mask_out, &
               is_in, ie_in, js_in, je_in, window_id)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real(FMS_TI_KIND_), dimension(:,:), intent(inout) :: data_in
      integer, intent(in), optional :: interp
      logical, intent(in), optional :: verbose
      type(horiz_interp_type),intent(in), optional :: horz_interp
      logical, dimension(:,:), intent(out), optional :: mask_out !< set to true where output data is valid
      integer,                  intent(in), optional :: is_in, ie_in, js_in, je_in
      integer,                  intent(in), optional :: window_id

      real(FMS_TI_KIND_), dimension(size(data_in,1), size(data_in,2), 1) :: data_out
      logical, dimension(size(data_in,1), size(data_in,2), 1) :: mask3d

      data_out(:,:,1) = data_in(:,:) ! fill initial values for the portions of array that are not touched by 3d routine
      call time_interp_external(index, time, data_out, interp, verbose, horz_interp, mask3d, &
                                   is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      data_in(:,:) = data_out(:,:,1)
      if (PRESENT(mask_out)) mask_out(:,:) = mask3d(:,:,1)

      return
    end subroutine TIME_INTERP_EXTERNAL_2D_



    !> 3D interpolation for @ref time_interp_external
    !! Provide data from external file interpolated to current model time.
    !! Data may be local to current processor or global, depending on
    !! "init_external_field" flags.
    subroutine TIME_INTERP_EXTERNAL_3D_(index, time, time_data, interp,verbose,horz_interp, mask_out, is_in, ie_in, &
                                      &  js_in, je_in, window_id)

        integer,                    intent(in)           :: index !< index of external field from previous call
                                                                  !! to init_external_field
        type(time_type),            intent(in)           :: time !< target time for data
        real(FMS_TI_KIND_), dimension(:,:,:),  intent(inout)           :: time_data !< global or local data array
        integer,                    intent(in), optional :: interp
        logical,                    intent(in), optional :: verbose !< flag for debugging
        type(horiz_interp_type),    intent(in), optional :: horz_interp
        logical, dimension(:,:,:), intent(out), optional :: mask_out !< set to true where output data is valid
        integer,                    intent(in), optional :: is_in, ie_in, js_in, je_in
        integer,                    intent(in), optional :: window_id

        integer :: nx, ny, nz, interp_method, t1, t2
        integer :: i1, i2, isc, iec, jsc, jec, mod_time
        integer :: yy, mm, dd, hh, min, ss
        character(len=256) :: err_msg
        character(len=FMS_PATH_LEN) :: filename

        integer :: isw, iew, jsw, jew, nxw, nyw
            ! these are boundaries of the updated portion of the "data" argument
            ! they are calculated using sizes of the "data" and isc,iec,jsc,jsc
            ! fileds from respective input field, to center the updated portion
            ! in the output array

        real(FMS_TI_KIND_) :: w1,w2
        logical :: verb
        character(len=16) :: message1, message2
        integer, parameter :: kindl = FMS_TI_KIND_

        nx = size(time_data,1)
        ny = size(time_data,2)
        nz = size(time_data,3)

        interp_method = LINEAR_TIME_INTERP
        if (PRESENT(interp)) interp_method = interp
        verb=.false.
        if (PRESENT(verbose)) verb=verbose
        if (debug_this_module) verb = .true.

        if (index < 1.or.index > num_fields) &
            call mpp_error(FATAL, &
                    & 'invalid index in call to time_interp_ext -- field was not initialized or failed to initialize')

        isc=loaded_fields(index)%isc;iec=loaded_fields(index)%iec
        jsc=loaded_fields(index)%jsc;jec=loaded_fields(index)%jec

        if( loaded_fields(index)%numwindows == 1 ) then
            nxw = iec-isc+1
            nyw = jec-jsc+1
        else
            if(.not. present(is_in) .or. .not. present(ie_in) .or. .not. present(js_in) .or. .not. present(je_in))then
                call mpp_error(FATAL, 'time_interp_external: is_in, ie_in, js_in and je_in must be present ' // &
                                    'when numwindows > 1, field='//trim(loaded_fields(index)%name))
            endif
            nxw = ie_in - is_in + 1
            nyw = je_in - js_in + 1
            isc = isc + is_in - 1
            iec = isc + ie_in - is_in
            jsc = jsc + js_in - 1
            jec = jsc + je_in - js_in
        endif

        isw = (nx-nxw)/2+1; iew = isw+nxw-1
        jsw = (ny-nyw)/2+1; jew = jsw+nyw-1

        if (nx < nxw .or. ny < nyw .or. nz < loaded_fields(index)%siz(3)) then
            write(message1,'(i6,2i5)') nx,ny,nz
            call mpp_error(FATAL,'field '//trim(loaded_fields(index)%name)// &
                        ' Array size mismatch in time_interp_external. Array "data" is too small. shape(data)=' &
                        //message1)
        endif
        if(PRESENT(mask_out)) then
            if (size(mask_out,1) /= nx .or. size(mask_out,2) /= ny .or. size(mask_out,3) /= nz) then
            write(message1,'(i6,2i5)') nx,ny,nz
            write(message2,'(i6,2i5)') size(mask_out,1),size(mask_out,2),size(mask_out,3)
            call mpp_error(FATAL,'field '//trim(loaded_fields(index)%name)// &
                                 ' array size mismatch in time_interp_external.'// &
                                 ' Shape of array "mask_out" does not match that of array "data".'// &
                                 ' shape(data)='//message1//' shape(mask_out)='//message2)
            endif
        endif

        ! only one record in the file => time-independent field
        if (loaded_fields(index)%siz(4) == 1) then
            call load_record(loaded_fields(index),1,horz_interp, is_in, ie_in ,js_in, je_in,window_id)
            i1 = find_buf_index(1,loaded_fields(index)%ibuf)
            if( loaded_fields(index)%region_type == NO_REGION ) then
                where(loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1))
                time_data(isw:iew,jsw:jew,:) = real( loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i1),&
                                                FMS_TI_KIND_)
                elsewhere
    !               time_data(isw:iew,jsw:jew,:) = time_interp_missing !field(index)%missing? Balaji
                    time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index)%missing, FMS_TI_KIND_)
                end where
            else
                where(loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1))
                time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i1),&
                                               FMS_TI_KIND_)
                end where
            endif
            if(PRESENT(mask_out)) mask_out(isw:iew,jsw:jew,:) = loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1)
        ! otherwise do interpolation
        else
            ! using time_interp_modulo
            if(loaded_fields(index)%have_modulo_times) then
                call time_interp(time,loaded_fields(index)%modulo_time_beg, loaded_fields(index)%modulo_time_end, &
                                 loaded_fields(index)%time(:), w2, t1, t2, &
                                 loaded_fields(index)%correct_leap_year_inconsistency, err_msg=err_msg)
                if(err_msg .NE. '') then
                    filename = trim(loaded_fields(index)%fileobj%path)
                    call mpp_error(FATAL,"time_interp_external 1: "//trim(err_msg)//&
                             ",file="//trim(filename)//",field="//trim(loaded_fields(index)%name) )
                endif
            ! using time_interp_list
            else
                if(loaded_fields(index)%modulo_time) then
                    mod_time=1
                else
                    mod_time=0
                endif
                call time_interp(time,loaded_fields(index)%time(:),w2,t1,t2,modtime=mod_time, err_msg=err_msg)
                if(err_msg .NE. '') then
                    filename = trim(loaded_fields(index)%fileobj%path)
                    call mpp_error(FATAL,"time_interp_external 2: "//trim(err_msg)//&
                                   ",file="//trim(filename)//",field="//trim(loaded_fields(index)%name) )
                endif
            endif
            w1 = 1.0_kindl -w2
            if (verb) then
                call get_date(time,yy,mm,dd,hh,min,ss)
                write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                    'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',min,':',ss
                write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
            endif

            call load_record(loaded_fields(index),t1,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
            call load_record(loaded_fields(index),t2,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
            i1 = find_buf_index(t1,loaded_fields(index)%ibuf)
            i2 = find_buf_index(t2,loaded_fields(index)%ibuf)
            if(i1<0.or.i2<0) &
                call mpp_error(FATAL,'time_interp_external : records were not loaded correctly in memory')

            if (verb) then
                write(outunit,*) 'ibuf= ',loaded_fields(index)%ibuf
                write(outunit,*) 'i1,i2= ',i1, i2
            endif

            if( loaded_fields(index)%region_type == NO_REGION ) then
                where(loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1) .and. &
                      loaded_fields(index)%mask(isc:iec,jsc:jec,:,i2))
                time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i1), kindl)&
                      * w1 + real(loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i2), kindl) * w2
                elsewhere
    !               time_data(isw:iew,jsw:jew,:) = time_interp_missing !loaded_fields(index)%missing? Balaji
                time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index)%missing, kindl)
                end where
            else
                where(loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1) .and. &
                      loaded_fields(index)%mask(isc:iec,jsc:jec,:,i2))
                time_data(isw:iew,jsw:jew,:) = real( loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i1), kindl)&
                      * w1 + real(loaded_fields(index)%domain_data(isc:iec,jsc:jec,:,i2), kindl) * w2
                end where
            endif
            if(PRESENT(mask_out)) &
                mask_out(isw:iew,jsw:jew,:) = &
                                            loaded_fields(index)%mask(isc:iec,jsc:jec,:,i1).and.&
                                            loaded_fields(index)%mask(isc:iec,jsc:jec,:,i2)
        endif

    end subroutine TIME_INTERP_EXTERNAL_3D_
!</SUBROUTINE> NAME="time_interp_external"

    !> @brief Scalar interpolation for @ref time_interp_external
    subroutine TIME_INTERP_EXTERNAL_0D_(index, time, time_data, verbose)

      integer, intent(in) :: index
      type(time_type), intent(in) :: time
      real(FMS_TI_KIND_), intent(inout) :: time_data
      logical, intent(in), optional :: verbose

      integer :: t1, t2
      integer :: i1, i2, mod_time
      integer :: yy, mm, dd, hh, min, ss
      character(len=256) :: err_msg
      character(len=FMS_PATH_LEN) :: filename

      real(FMS_TI_KIND_) :: w1,w2
      logical :: verb
      integer, parameter :: kindl = FMS_TI_KIND_

      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.

      if (index < 1.or.index > num_fields) &
           call mpp_error(FATAL, &
                     & 'invalid index in call to time_interp_ext -- field was not initialized or failed to initialize')

      if (loaded_fields(index)%siz(4) == 1) then
         ! only one record in the file => time-independent loaded_fields
         call load_record_0d(loaded_fields(index),1)
         i1 = find_buf_index(1,loaded_fields(index)%ibuf)
         time_data = real(loaded_fields(index)%domain_data(1,1,1,i1), FMS_TI_KIND_)
      else
        if(loaded_fields(index)%have_modulo_times) then
          call time_interp(time,loaded_fields(index)%modulo_time_beg, loaded_fields(index)%modulo_time_end, &
                           loaded_fields(index)%time(:), w2, t1, t2, &
                           loaded_fields(index)%correct_leap_year_inconsistency, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(loaded_fields(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 3:"//trim(err_msg)//&
                    ",file="//trim(filename)//",field="//trim(loaded_fields(index)%name) )
          endif
        else
          if(loaded_fields(index)%modulo_time) then
            mod_time=1
          else
            mod_time=0
          endif
          call time_interp(time,loaded_fields(index)%time(:),w2,t1,t2,modtime=mod_time, err_msg=err_msg)
          if(err_msg .NE. '') then
             filename = trim(loaded_fields(index)%fileobj%path)
             call mpp_error(FATAL,"time_interp_external 4:"//trim(err_msg)// &
                    ",file="//trim(filename)//",field="//trim(loaded_fields(index)%name) )
          endif
        endif
         w1 = 1.0_kindl-w2
         if (verb) then
            call get_date(time,yy,mm,dd,hh,min,ss)
            write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
                 'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',min,':',ss
            write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
         endif
         call load_record_0d(loaded_fields(index),t1)
         call load_record_0d(loaded_fields(index),t2)
         i1 = find_buf_index(t1,loaded_fields(index)%ibuf)
         i2 = find_buf_index(t2,loaded_fields(index)%ibuf)

         if(i1<0.or.i2<0) &
              call mpp_error(FATAL,'time_interp_external : records were not loaded correctly in memory')
         time_data = real(loaded_fields(index)%domain_data(1,1,1,i1), FMS_TI_KIND_)*w1 &
              & + real(loaded_fields(index)%domain_data(1,1,1,i2), FMS_TI_KIND_)*w2
         if (verb) then
            write(outunit,*) 'ibuf= ',loaded_fields(index)%ibuf
            write(outunit,*) 'i1,i2= ',i1, i2
         endif
      endif

    end subroutine TIME_INTERP_EXTERNAL_0D_

! ============================================================================
