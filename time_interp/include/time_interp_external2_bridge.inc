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

    !> @brief 2D interpolation for @ref time_interp_external_bridge
    subroutine TIME_INTERP_EXTERNAL_BRIDGE_2D_(index1, index2, time, data_in, interp, verbose,horz_interp, mask_out, &
               is_in, ie_in, js_in, je_in, window_id)

      integer, intent(in) :: index1 !< index of first external field
      integer, intent(in) :: index2 !< index of second external field
      type(time_type), intent(in) :: time !< target time for data
      real(FMS_TI_KIND_), dimension(:,:), intent(inout) :: data_in !< global or local data array
      integer, intent(in), optional :: interp !< hardcoded to linear
      logical, intent(in), optional :: verbose !< flag for debugging
      type(horiz_interp_type),intent(in), optional :: horz_interp !< horizontal interpolation type
      logical, dimension(:,:), intent(out), optional :: mask_out !< set to true where output data is valid
      integer,                  intent(in), optional :: is_in, ie_in, js_in, je_in !< horizontal indices for load_record
      integer,                  intent(in), optional :: window_id !< harcoded to 1 in load_record

      real(FMS_TI_KIND_), dimension(size(data_in,1), size(data_in,2), 1) :: data_out !< 3d version of data_in
      logical, dimension(size(data_in,1), size(data_in,2), 1) :: mask3d !< 3d version of mask_out

      data_out(:,:,1) = data_in(:,:) ! fill initial values for the portions of array that are not touched by 3d routine
      call time_interp_external_bridge(index1, index2, time, data_out, interp, verbose, horz_interp, mask3d, &
                                   is_in=is_in,ie_in=ie_in,js_in=js_in,je_in=je_in,window_id=window_id)
      data_in(:,:) = data_out(:,:,1)
      if (PRESENT(mask_out)) mask_out(:,:) = mask3d(:,:,1)

      return
    end subroutine TIME_INTERP_EXTERNAL_BRIDGE_2D_

    !> 3D interpolation for @ref time_interp_external
    !! Provide data from external file interpolated to current model time.
    !! Data may be local to current processor or global, depending on
    !! "init_external_field" flags.
    subroutine TIME_INTERP_EXTERNAL_BRIDGE_3D_(index1, index2, time, time_data, interp,verbose,horz_interp, mask_out, &
                                             &  is_in, ie_in, js_in, je_in, window_id)

      integer,                    intent(in) :: index1 !< index of first external field
      integer,                    intent(in) :: index2 !< index of second external field
      type(time_type),            intent(in) :: time   !< target time for data
      real(FMS_TI_KIND_), dimension(:,:,:),  intent(inout) :: time_data !< global or local data array
      integer,                    intent(in), optional :: interp !< hardcoded to linear
      logical,                    intent(in), optional :: verbose !< flag for debugging
      type(horiz_interp_type),    intent(in), optional :: horz_interp !< horizontal interpolation type
      logical, dimension(:,:,:), intent(out), optional :: mask_out !< set to true where output data is valid
      integer,                    intent(in), optional :: is_in, ie_in !<  x horizontal indices for load_record
      integer,                    intent(in), optional :: js_in, je_in !< y horizontal indices for load_record
      integer,                    intent(in), optional :: window_id !< harcoded to 1 in load_record

      type(time_type) :: time1 !< time type associated with index1 of external field
      type(time_type) :: time2 !< time type associated with index2 of external field
      integer :: dims1(4) !< dimensions XYZT of index1
      integer :: dims2(4) !< dimensions XYZT of index2
      integer :: nx, ny, nz !< size in X,Y,Z of array
      integer :: interp_method !< hardcoded to linear in time
      integer :: t1, t2 !< temporary to store time index
      integer :: i1, i2 !< temporary to store time index
      integer :: isc, iec, jsc, jec !< start/end arrays in X,Y
      integer :: isc1, iec1, jsc1, jec1 !< start/end arrays in X,Y for field1
      integer :: isc2, iec2, jsc2, jec2 !< start/end arrays in X,Y for field2
      integer :: yy, mm, dd, hh, minu, ss !< year, month, day, hour, minute, second for date

      integer :: isw, iew, jsw, jew, nxw, nyw !< these are boundaries of the updated portion of the "data" argument
                                              !! they are calculated using sizes of the "data" and isc,iec,jsc,jsc
                                              !! fileds from respective input field, to center the updated portion
                                              !! in the output array

      real(FMS_TI_KIND_) :: w1 !< interp weight for index1
      real(FMS_TI_KIND_) :: w2 !< interp weight for index2
      logical :: verb !< verbose
      character(len=16) :: message1 !< temp string
      character(len=16) :: message2 !< temp string
      integer, parameter :: kindl = FMS_TI_KIND_

      nx = size(time_data,1)
      ny = size(time_data,2)
      nz = size(time_data,3)

      interp_method = LINEAR_TIME_INTERP
      if (PRESENT(interp)) interp_method = interp
      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.

      if (index1 < 1 .or. index1 > num_fields) &
           call mpp_error(FATAL,'invalid index1 in call to time_interp_external_bridge:' // &
                                '-- field was not initialized or failed to initialize')
      if (index2 < 1 .or. index2 > num_fields) &
           call mpp_error(FATAL,'invalid index2 in call to time_interp_external_bridge:' // &
                                ' -- field was not initialized or failed to initialize')

      isc1=loaded_fields(index1)%isc;iec1=loaded_fields(index1)%iec
      jsc1=loaded_fields(index1)%jsc;jec1=loaded_fields(index1)%jec
      isc2=loaded_fields(index2)%isc;iec2=loaded_fields(index2)%iec
      jsc2=loaded_fields(index2)%jsc;jec2=loaded_fields(index2)%jec

      if ((isc1 /= isc2) .or. (iec1 /= iec2) .or. (jsc1 /= jsc2) .or. (jec1 /= jec2)) then
         call mpp_error(FATAL, 'time_interp_external_bridge: dimensions must be the same in both index1 and index2 ')
      endif

      isc=isc1 ; iec=iec1 ; jsc=jsc1 ; jec=jec1

      if (trim(loaded_fields(index2)%name) /= trim(loaded_fields(index1)%name)) then
         call mpp_error(FATAL, 'time_interp_external_bridge: cannot bridge between different fields.' // &
                        'field1='//trim(loaded_fields(index1)%name)//' field2='//trim(loaded_fields(index2)%name))
      endif

      if ((loaded_fields(index1)%numwindows == 1) .and. (loaded_fields(index2)%numwindows == 1)) then
         nxw = iec-isc+1
         nyw = jec-jsc+1
      else
         if(.not. present(is_in) .or. .not. present(ie_in) .or. .not. present(js_in) .or. .not. present(je_in))then
             call mpp_error(FATAL, 'time_interp_external: is_in, ie_in, js_in and je_in must be present ' // &
                                    'when numwindows > 1, field='//trim(loaded_fields(index1)%name))
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

      if (nx < nxw .or. ny < nyw .or. nz < loaded_fields(index1)%siz(3)) then
         write(message1,'(i6,2i5)') nx,ny,nz
         call mpp_error(FATAL,'field '//trim(loaded_fields(index1)%name)// &
                              ' Array size mismatch in time_interp_external_bridge.'// &
                              ' Array "data" is too small. shape(time_data)='//message1)
      endif
      if (nx < nxw .or. ny < nyw .or. nz < loaded_fields(index2)%siz(3)) then
         write(message1,'(i6,2i5)') nx,ny,nz
         call mpp_error(FATAL,'field '//trim(loaded_fields(index2)%name)// &
                              ' Array size mismatch in time_interp_external_bridge.'// &
                              ' Array "data" is too small. shape(time_data)='//message1)
      endif

      if(PRESENT(mask_out)) then
        if (size(mask_out,1) /= nx .or. size(mask_out,2) /= ny .or. size(mask_out,3) /= nz) then
          write(message1,'(i6,2i5)') nx,ny,nz
          write(message2,'(i6,2i5)') size(mask_out,1),size(mask_out,2),size(mask_out,3)
          call mpp_error(FATAL,'field '//trim(loaded_fields(index1)%name)// &
                               ' array size mismatch in time_interp_external_bridge.'// &
                               ' Shape of array "mask_out" does not match that of array "data".'// &
                               ' shape(time_data)='//message1//' shape(mask_out)='//message2)
        endif
      endif

      if ((loaded_fields(index1)%have_modulo_times) .or. (loaded_fields(index2)%have_modulo_times)) then
          call mpp_error(FATAL, 'time_interp_external_bridge: field '//trim(loaded_fields(index1)%name)// &
                                ' array cannot have modulo time')
      endif
      if ((loaded_fields(index1)%modulo_time) .or. (loaded_fields(index2)%modulo_time)) then
          call mpp_error(FATAL, 'time_interp_external_bridge: field '//trim(loaded_fields(index1)%name)// &
                                ' array cannot have modulo time')
      endif

      dims1 = get_external_field_size(index1)
      dims2 = get_external_field_size(index2)

      t1 = dims1(4)
      t2 = 1

      time1 = loaded_fields(index1)%time(t1)
      time2 = loaded_fields(index2)%time(t2)
      w2= (time - time1) // (time2 - time1)
      w1 = 1.0_kindl-w2

      if (verb) then
         call get_date(time,yy,mm,dd,hh,minu,ss)
         write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
              'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',minu,':',ss
         write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
      endif

      call load_record(loaded_fields(index1),t1,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
      call load_record(loaded_fields(index2),t2,horz_interp, is_in, ie_in ,js_in, je_in, window_id)
      i1 = find_buf_index(t1,loaded_fields(index1)%ibuf)
      i2 = find_buf_index(t2,loaded_fields(index2)%ibuf)
      if(i1<0.or.i2<0) &
           call mpp_error(FATAL,'time_interp_external_bridge : records were not loaded correctly in memory')

      if (verb) then
         write(outunit,*) 'ibuf= ',loaded_fields(index2)%ibuf
         write(outunit,*) 'i1,i2= ',i1, i2
      endif

      if((loaded_fields(index1)%region_type == NO_REGION) .and. (loaded_fields(index2)%region_type == NO_REGION) ) then
         where(loaded_fields(index1)%mask(isc:iec,jsc:jec,:,i1).and.loaded_fields(index2)%mask(isc:iec,jsc:jec,:,i2))
            time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index1)%domain_data(isc:iec,jsc:jec,:,i1), kindl)*w1 + &
                 real(loaded_fields(index2)%domain_data(isc:iec,jsc:jec,:,i2), kindl)*w2
         elsewhere
            time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index1)%missing, kindl)
         end where
      else
         where(loaded_fields(index1)%mask(isc:iec,jsc:jec,:,i1).and.loaded_fields(index2)%mask(isc:iec,jsc:jec,:,i2))
            time_data(isw:iew,jsw:jew,:) = real(loaded_fields(index1)%domain_data(isc:iec,jsc:jec,:,i1), kindl)*w1 + &
                 real(loaded_fields(index2)%domain_data(isc:iec,jsc:jec,:,i2), kindl)*w2
         end where
      endif

      if(PRESENT(mask_out)) &
           mask_out(isw:iew,jsw:jew,:) = &
                                     loaded_fields(index1)%mask(isc:iec,jsc:jec,:,i1).and.&
                                     loaded_fields(index2)%mask(isc:iec,jsc:jec,:,i2)

    end subroutine TIME_INTERP_EXTERNAL_BRIDGE_3D_

    subroutine TIME_INTERP_EXTERNAL_BRIDGE_0D_(index1, index2, time, time_data, verbose)

      integer, intent(in) :: index1 !< index of first external field
      integer, intent(in) :: index2 !< index of second external field
      type(time_type), intent(in) :: time !< target time for data
      real(FMS_TI_KIND_), intent(inout) :: time_data !< global or local data array
      logical, intent(in), optional :: verbose !< flag for debugging

      integer :: t1, t2 !< temporary to store time index
      integer :: i1, i2 !< temporary to store time index
      integer :: yy, mm, dd, hh, minu, ss !< year, month, day, hour, minute, second for date
      type(time_type) :: time1 !< time type associated with index1 of external field
      type(time_type) :: time2 !< time type associated with index2 of external field
      integer :: dims1(4) !< dimensions XYZT of index1
      integer :: dims2(4) !< dimensions XYZT of index2

      real(FMS_TI_KIND_) :: w1 !< interp weight for index1
      real(FMS_TI_KIND_) :: w2 !< interp weight for index2
      logical :: verb !< verbose
      integer, parameter :: kindl = FMS_TI_KIND_

      verb=.false.
      if (PRESENT(verbose)) verb=verbose
      if (debug_this_module) verb = .true.

      if (index1 < 0 .or. index1 > num_fields) &
           call mpp_error(FATAL,'invalid index1 in call to time_interp_ext' // &
                                ' -- field was not initialized or failed to initialize')
      if (index2 < 0 .or. index2 > num_fields) &
           call mpp_error(FATAL,'invalid index2 in call to time_interp_ext' // &
                                ' -- field was not initialized or failed to initialize')

      if ((loaded_fields(index1)%have_modulo_times) .or. (loaded_fields(index2)%have_modulo_times)) then
          call mpp_error(FATAL, 'time_interp_external_bridge: field '//trim(loaded_fields(index1)%name)// &
                                ' array cannot have modulo time')
      endif
      if ((loaded_fields(index1)%modulo_time) .or. (loaded_fields(index2)%modulo_time)) then
          call mpp_error(FATAL, 'time_interp_external_bridge: field '//trim(loaded_fields(index1)%name)// &
                                ' array cannot have modulo time')
      endif

      dims1 = get_external_field_size(index1)
      dims2 = get_external_field_size(index2)

      t1 = dims1(4)
      t2 = 1

      time1 = loaded_fields(index1)%time(t1)
      time2 = loaded_fields(index2)%time(t2)
      w2= (time - time1) // (time2 - time1)
      w1 = 1.0_kindl-w2

      if (verb) then
         call get_date(time,yy,mm,dd,hh,minu,ss)
         write(outunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') &
              'target time yyyy/mm/dd hh:mm:ss= ',yy,'/',mm,'/',dd,hh,':',minu,':',ss
         write(outunit,*) 't1, t2, w1, w2= ', t1, t2, w1, w2
      endif
      call load_record_0d(loaded_fields(index2),t1)
      call load_record_0d(loaded_fields(index2),t2)
      i1 = find_buf_index(t1,loaded_fields(index2)%ibuf)
      i2 = find_buf_index(t2,loaded_fields(index2)%ibuf)

      if(i1<0.or.i2<0) &
           call mpp_error(FATAL,'time_interp_external : records were not loaded correctly in memory')
      time_data = real(loaded_fields(index2)%domain_data(1,1,1,i1), FMS_TI_KIND_)*w1 &
                + real(loaded_fields(index2)%domain_data(1,1,1,i2), FMS_TI_KIND_)*w2
      if (verb) then
         write(outunit,*) 'ibuf= ',loaded_fields(index2)%ibuf
         write(outunit,*) 'i1,i2= ',i1, i2
      endif

    end subroutine TIME_INTERP_EXTERNAL_BRIDGE_0D_


! ============================================================================
