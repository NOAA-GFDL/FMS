    subroutine MPP_WRITE_( unit, field, data, tstamp )
      integer, intent(in) :: unit
      type(fieldtype), intent(in) :: field
      MPP_TYPE_, intent(in) :: data MPP_RANK_
      real, intent(in), optional :: tstamp

      MPP_WRITE_RECORD_
      return
    end subroutine MPP_WRITE_
