      TYPE CAFPNTR_TYPE_3D_
        MPP_TYPE_,dimension(:,:,:), pointer :: pfield
      END TYPE CAFPNTR_TYPE_3D_
      type(CAFPNTR_TYPE_3D_), dimension(MAX_DOMAIN_FIELDS) :: CAFPNTR_3D_[0:*]

      TYPE CAFPNTR_TYPE_1D_
        MPP_TYPE_,dimension(:), pointer :: pfield
      END TYPE CAFPNTR_TYPE_1D_
      type(CAFPNTR_TYPE_1D_) :: CAFPNTR_1D_[0:*]

