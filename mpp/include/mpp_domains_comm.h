    subroutine MPP_ASSOC_CAF_FIELD_3D_(is,ie,js,je,ke,f,cafptr)
      use mpp_datatype_mod, only: CAFPNTR_TYPE_3D_
      integer, intent(in) :: is,ie,js,je,ke
      MPP_TYPE_,target, intent(inout) :: f(is:ie,js:je,ke)
      type(CAFPNTR_TYPE_3D_),intent(inout) :: cafptr

      cafptr%pfield =>f
    end subroutine MPP_ASSOC_CAF_FIELD_3D_
