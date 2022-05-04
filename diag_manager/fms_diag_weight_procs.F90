MODULE fms_diag_weight_procs_mod
   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION weigh_field ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL :: weigh_field
      END FUNCTION weigh_field
   END INTERFACE

   ABSTRACT INTERFACE
      FUNCTION weigh_field_3d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:,:,:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d
      END FUNCTION weigh_field_3d
   END INTERFACE


   ABSTRACT INTERFACE
      FUNCTION weigh_field_1d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d
      END FUNCTION weigh_field_1d
   END INTERFACE

   TYPE :: FmsWeightProcCfg_t
      INTEGER :: pow_value
      procedure (weigh_field), pointer, nopass:: fwf_0d_ptr   !! A pointer to the field weighing function.
      procedure (weigh_field_1d),  pointer, nopass:: fwf_1d_ptr !! A pointer to the 3D field weight function.
      procedure (weigh_field_3d), pointer, nopass :: fwf_3d_ptr !! A pointer to the 3D field weight function.
   CONTAINS
      procedure :: initialize => initialize_weight_proc_cfg
   END TYPE FmsWeightProcCfg_t


CONTAINS

   PURE REAL FUNCTION weigh_field_0d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weigh_field_0d_p1 = field_val * weight
   END FUNCTION weigh_field_0d_p1

   PURE REAL FUNCTION weigh_field_0d_p2 ( field_val, weight,  pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      REAL :: fTw
      fTw =  field_val * weight
      weigh_field_0d_p2 = fTw * fTw
   END FUNCTION weigh_field_0d_p2

   PURE REAL FUNCTION weigh_field_0d_pp ( field_val, weight, pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weigh_field_0d_pp = (field_val * weight) ** pow_value
   END FUNCTION weigh_field_0d_pp

   PURE FUNCTION weigh_field_1d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_p1
      ALLOCATE(weigh_field_1d_p1(size(field_val)))
      DO i = 1, size(field_val)
         weigh_field_1d_p1(i) = field_val(i) * weight
      END DO
   END FUNCTION weigh_field_1d_p1

   PURE FUNCTION weigh_field_1d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_p2
      !!TODO:verify the allocate
      !!TODO: name change: weigh_field (weight_the_feild) -- >> weight_the_field
      ALLOCATE(weigh_field_1d_p2, mold=field_val)
      DO i = 1, size(field_val)
         weigh_field_1d_p2(i) = field_val(i) * field_val(i) * weight * weight
      END DO
   END FUNCTION weigh_field_1d_p2

   PURE FUNCTION weigh_field_1d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weigh_field_1d_pp
      ALLOCATE(weigh_field_1d_pp(size(field_val)))
      DO i = 1, size(field_val)
         weigh_field_1d_pp(i) = (field_val(i) * weight) ** pow_value
      END DO
   END FUNCTION weigh_field_1d_pp

   PURE FUNCTION weigh_field_3d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p1
      ALLOCATE(weigh_field_3d_p1, mold=field_val)
      weigh_field_3d_p1 = field_val * weight
   END FUNCTION weigh_field_3d_p1

   PURE FUNCTION weigh_field_3d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_p2
      ALLOCATE(weigh_field_3d_p2, mold=field_val)
      weigh_field_3d_p2 = field_val * field_val * weight * weight
   END FUNCTION weigh_field_3d_p2

   PURE FUNCTION weigh_field_3d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weigh_field_3d_pp
      ALLOCATE(weigh_field_3d_pp, mold=field_val)
      weigh_field_3d_pp = (field_val * weight) ** pow_value
   END FUNCTION weigh_field_3d_pp



   !> @brief The constructor for the WEIGHT_PROC_CFG type.
   !! @return Returns a newly allocated list iterator.
   subroutine initialize_weight_proc_cfg ( this, pow_value )
      class (FmsWeightProcCfg_t), intent(inout) :: this
      integer, intent(in) :: pow_value

      this%pow_value = pow_value

      if ( pow_value == 1) then
         this%fwf_0d_ptr => weigh_field_0d_p1
         this%fwf_1D_ptr => weigh_field_1d_p1
         this%fwf_3D_ptr => weigh_field_3d_p1
      else if ( pow_value == 2 ) then
         this%fwf_0d_ptr => weigh_field_0d_p2
         this%fwf_1d_ptr => weigh_field_1d_p2
         this%fwf_3D_ptr => weigh_field_3d_p2
      else
         this%fwf_0d_ptr => weigh_field_0d_pp
         this%fwf_1D_ptr => weigh_field_1d_pp
         this%fwf_3D_ptr => weigh_field_3d_pp
      end if

   end subroutine initialize_weight_proc_cfg


  END MODULE fms_diag_weight_procs_mod
