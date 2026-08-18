#ifndef __NVCOMPILER

#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_r8
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_r8
#undef ARR_INIT_
#define ARR_INIT_ arr_init_r8
#include <mpp_pack.fh>

#ifdef OVERLOAD_C8
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_c8
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_c8
#undef ARR_INIT_
#define ARR_INIT_ arr_init_c8
#include <mpp_pack.fh>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_i8
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_i8
#undef ARR_INIT_
#define ARR_INIT_ arr_init_i8
#include <mpp_pack.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(l8_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_l8
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_l8
#undef ARR_INIT_
#define ARR_INIT_ arr_init_l8
#include <mpp_pack.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_r4
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_r4
#undef ARR_INIT_
#define ARR_INIT_ arr_init_r4
#include <mpp_pack.fh>

#ifdef OVERLOAD_C4
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_c4
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_c4
#undef ARR_INIT_
#define ARR_INIT_ arr_init_c4
#include <mpp_pack.fh>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_i4
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_i4
#undef ARR_INIT_
#define ARR_INIT_ arr_init_i4
#include <mpp_pack.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(l4_kind)
#undef ARR2VEC_
#define ARR2VEC_ arr2vec_l4
#undef VEC2ARR_
#define VEC2ARR_ vec2arr_l4
#undef ARR_INIT_
#define ARR_INIT_ arr_init_l4
#include <mpp_pack.fh>

#else

#undef MPP_TYPE_
#define MPP_TYPE_ real(r8_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_r8_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_r8_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_r8_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_r8_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_r8_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_r8_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_r8_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_r8_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_r8_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_r8_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_r8_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_r8_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_r8_5d
#include <mpp_pack_compat.fh>

#ifdef OVERLOAD_C8
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c8_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_c8_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_c8_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_c8_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_c8_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_c8_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_c8_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_c8_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_c8_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_c8_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_c8_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_c8_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_c8_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_c8_5d
#include <mpp_pack_compat.fh>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(i8_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_i8_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_i8_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_i8_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_i8_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_i8_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_i8_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_i8_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_i8_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_i8_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_i8_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_i8_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_i8_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_i8_5d
#include <mpp_pack_compat.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(l8_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_l8_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_l8_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_l8_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_l8_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_l8_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_l8_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_l8_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_l8_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_l8_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_l8_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_l8_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_l8_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_l8_5d
#include <mpp_pack_compat.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ real(r4_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_r4_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_r4_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_r4_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_r4_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_r4_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_r4_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_r4_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_r4_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_r4_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_r4_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_r4_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_r4_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_r4_5d
#include <mpp_pack_compat.fh>

#ifdef OVERLOAD_C4
#undef MPP_TYPE_
#define MPP_TYPE_ complex(c4_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_c4_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_c4_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_c4_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_c4_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_c4_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_c4_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_c4_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_c4_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_c4_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_c4_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_c4_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_c4_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_c4_5d
#include <mpp_pack_compat.fh>
#endif

#undef MPP_TYPE_
#define MPP_TYPE_ integer(i4_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_i4_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_i4_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_i4_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_i4_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_i4_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_i4_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_i4_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_i4_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_i4_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_i4_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_i4_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_i4_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_i4_5d
#include <mpp_pack_compat.fh>

#undef MPP_TYPE_
#define MPP_TYPE_ logical(l4_kind)
#undef ARR2VEC_2D
#define ARR2VEC_2D arr2vec_l4_2d
#undef ARR2VEC_3D
#define ARR2VEC_3D arr2vec_l4_3d
#undef ARR2VEC_4D
#define ARR2VEC_4D arr2vec_l4_4d
#undef ARR2VEC_5D
#define ARR2VEC_5D arr2vec_l4_5d
#undef VEC2ARR_2D
#define VEC2ARR_2D vec2arr_l4_2d
#undef VEC2ARR_3D
#define VEC2ARR_3D vec2arr_l4_3d
#undef VEC2ARR_4D
#define VEC2ARR_4D vec2arr_l4_4d
#undef VEC2ARR_5D
#define VEC2ARR_5D vec2arr_l4_5d
#undef ARR_INIT_1D
#define ARR_INIT_1D arr_init_l4_1d
#undef ARR_INIT_2D
#define ARR_INIT_2D arr_init_l4_2d
#undef ARR_INIT_3D
#define ARR_INIT_3D arr_init_l4_3d
#undef ARR_INIT_4D
#define ARR_INIT_4D arr_init_l4_4d
#undef ARR_INIT_5D
#define ARR_INIT_5D arr_init_l4_5d
#include <mpp_pack_compat.fh>

#endif
