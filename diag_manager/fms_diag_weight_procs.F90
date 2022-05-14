MODULE fms_diag_weight_procs_mod
   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION weight_the_field ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL :: weight_the_field
      END FUNCTION weight_the_field
   END INTERFACE

   ABSTRACT INTERFACE
      FUNCTION weight_the_field_3d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:,:,:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d
      END FUNCTION weight_the_field_3d
   END INTERFACE

   ABSTRACT INTERFACE
      FUNCTION weight_the_field_1d ( field_val, weight, pow_value )
         REAL, INTENT(in) :: field_val(:)
         REAL, INTENT(in) :: weight
         INTEGER, INTENT(in) :: pow_value
         REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d
      END FUNCTION weight_the_field_1d
   END INTERFACE

   TYPE :: FmsWeightProcCfg_t
      INTEGER :: pow_value
      procedure (weight_the_field), pointer, nopass::fwf_0d_ptr   !! A pointer to the field weighting function.
      procedure (weight_the_field_1d), pointer, nopass::fwf_1d_ptr !! A pointer to the 3D field weighting function.
      procedure (weight_the_field_3d), pointer, nopass ::fwf_3d_ptr !! A pointer to the 3D field weighting function.
   CONTAINS
      procedure :: initialize => initialize_weight_proc_cfg
   END TYPE FmsWeightProcCfg_t

CONTAINS

   PURE REAL FUNCTION weight_the_field_0d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weight_the_field_0d_p1 = field_val * weight
   END FUNCTION weight_the_field_0d_p1

   PURE REAL FUNCTION weight_the_field_0d_p2 ( field_val, weight,  pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      REAL :: fTw
      fTw =  field_val * weight
      weight_the_field_0d_p2 = fTw * fTw
   END FUNCTION weight_the_field_0d_p2

   PURE REAL FUNCTION weight_the_field_0d_pp ( field_val, weight, pow_value  )
      REAL, INTENT(in) :: weight
      REAL, INTENT(in) :: field_val
      INTEGER, INTENT(in) :: pow_value
      weight_the_field_0d_pp = (field_val * weight) ** pow_value
   END FUNCTION weight_the_field_0d_pp

   PURE FUNCTION weight_the_field_1d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_p1
      ALLOCATE(weight_the_field_1d_p1(size(field_val)))
      DO i = 1, size(field_val)
         weight_the_field_1d_p1(i) = field_val(i) * weight
      END DO
   END FUNCTION weight_the_field_1d_p1

   PURE FUNCTION weight_the_field_1d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_p2
      !!TODO:verify the allocate
      ALLOCATE(weight_the_field_1d_p2, mold=field_val)
      DO i = 1, size(field_val)
         weight_the_field_1d_p2(i) = field_val(i) * field_val(i) * weight * weight
      END DO
   END FUNCTION weight_the_field_1d_p2

   PURE FUNCTION weight_the_field_1d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      INTEGER :: i
      REAL, DIMENSION(:), ALLOCATABLE :: weight_the_field_1d_pp
      ALLOCATE(weight_the_field_1d_pp(size(field_val)))
      DO i = 1, size(field_val)
         weight_the_field_1d_pp(i) = (field_val(i) * weight) ** pow_value
      END DO
   END FUNCTION weight_the_field_1d_pp

   PURE FUNCTION weight_the_field_3d_p1 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_p1
      ALLOCATE(weight_the_field_3d_p1, mold=field_val)
      weight_the_field_3d_p1 = field_val * weight
   END FUNCTION weight_the_field_3d_p1

   PURE FUNCTION weight_the_field_3d_p2 ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_p2
      ALLOCATE(weight_the_field_3d_p2, mold=field_val)
      weight_the_field_3d_p2 = field_val * field_val * weight * weight
   END FUNCTION weight_the_field_3d_p2

   PURE FUNCTION weight_the_field_3d_pp ( field_val, weight, pow_value )
      REAL, INTENT(in) :: field_val(:,:,:)
      REAL, INTENT(in) :: weight
      INTEGER, INTENT(in) :: pow_value
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: weight_the_field_3d_pp
      ALLOCATE(weight_the_field_3d_pp, mold=field_val)
      weight_the_field_3d_pp = (field_val * weight) ** pow_value
   END FUNCTION weight_the_field_3d_pp

   !> @brief The constructor for the WEIGHT_PROC_CFG type.
   !! @return Returns a newly allocated list iterator.
   subroutine initialize_weight_proc_cfg ( this, pow_value )
      class (FmsWeightProcCfg_t), intent(inout) :: this
      integer, intent(in) :: pow_value

      this%pow_value = pow_value

      if ( pow_value == 1) then
         this%fwf_0d_ptr => weight_the_field_0d_p1
         this%fwf_1D_ptr => weight_the_field_1d_p1
         this%fwf_3D_ptr => weight_the_field_3d_p1
      else if ( pow_value == 2 ) then
         this%fwf_0d_ptr => weight_the_field_0d_p2
         this%fwf_1d_ptr => weight_the_field_1d_p2
         this%fwf_3D_ptr => weight_the_field_3d_p2
      else
         this%fwf_0d_ptr => weight_the_field_0d_pp
         this%fwf_1D_ptr => weight_the_field_1d_pp
         this%fwf_3D_ptr => weight_the_field_3d_pp
      end if

   end subroutine initialize_weight_proc_cfg


  END MODULE fms_diag_weight_procs_mod
