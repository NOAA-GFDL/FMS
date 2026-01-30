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
