
! -*-f90-*-
! $Id: mpp_domains_reduce.inc,v 13.0 2006/03/28 21:40:40 fms Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_REDUCE: get global max/min of field                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#ifndef no_4byte_reals
#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i8_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i8_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i8_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ maxval
#undef REDUCE_LOC_
#define REDUCE_LOC_ maxloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#undef MPP_GLOBAL_REDUCE_2D_
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i4_2d
#undef MPP_GLOBAL_REDUCE_3D_
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i4_3d
#undef MPP_GLOBAL_REDUCE_4D_
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i4_4d
#undef MPP_GLOBAL_REDUCE_5D_
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#undef REDUCE_VAL_
#define REDUCE_VAL_ minval
#undef REDUCE_LOC_
#define REDUCE_LOC_ minloc
#undef MPP_REDUCE_
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                   MPP_GLOBAL_SUM: global sum of field                       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#ifndef no_4byte_reals
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>
#endif

#ifndef no_8byte_integers
#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>
#endif

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_2d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_3d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_4d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#undef MPP_GLOBAL_SUM_
#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_5d
#undef MPP_EXTRA_INDICES_
#define MPP_EXTRA_INDICES_ ,:,:,:
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_FIELD: get global field from domain field           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_field.h>

#ifndef no_8byte_integers
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l8_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l8_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l8_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l8_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_global_field.h>
#endif

#ifndef no_4byte_reals
#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_field.h>
#endif

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field.h>

#undef MPP_GLOBAL_FIELD_2D_
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l4_2d
#undef MPP_GLOBAL_FIELD_3D_
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l4_3d
#undef MPP_GLOBAL_FIELD_4D_
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l4_4d
#undef MPP_GLOBAL_FIELD_5D_
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l4_5d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_global_field.h>

!****************************************************
#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_r8_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_r8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_c8_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_c8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#ifndef no_8byte_integers
#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_i8_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_i8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_l8_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_l8_3d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#endif

#ifndef no_4byte_reals
#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_r4_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_r4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_c4_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_c4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#endif

#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_i4_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_i4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#undef MPP_DO_GLOBAL_FIELD_3Dnew_
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_l4_3d
#undef MPP_DO_GLOBAL_FIELD_3Dold_
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_l4_3d
#undef MPP_TYPE_
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)
#undef CAFPNTR_TYPE_3D_
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#include <mpp_do_global_field_caf.h>
#else
#include <mpp_do_global_field_new.h>
#endif

