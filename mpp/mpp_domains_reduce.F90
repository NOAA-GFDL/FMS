module mpp_domains_reduce_mod
#include <fms_platform.h>

use mpp_mod,           only : mpp_error, FATAL, WARNING, mpp_npes, mpp_max, mpp_min, mpp_sum
use mpp_mod,           only : mpp_broadcast, mpp_transmit, mpp_sync_self, mpp_sync, mpp_send, mpp_recv 
use mpp_mod,           only : mpp_pe,mpp_root_pe,stdout,mpp_chksum
use mpp_parameter_mod, only : NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, XUPDATE, YUPDATE
use mpp_datatype_mod,  only : domain1D, domain2D, DomainCommunicator2D
use mpp_data_mod,      only : module_is_initialized=>mpp_domains_is_initialized, pe, mpp_domains_stack
use mpp_data_mod,      only : ptr_domains_stack, mpp_domains_stack_size, mpp_domains_stack_hwm, debug_gsm
use mpp_domains_comm_mod, only : mpp_global_field_init_comm, mpp_global_field_free_comm

implicit none
private

  public :: mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum

  character(len=128), public :: version= &
       '$Id: mpp_domains_reduce.F90,v 12.0 2005/04/14 17:58:12 fms Exp $'
  character(len=128), public :: tagname= &
       '$Name: lima $'


! <INTERFACE NAME="mpp_global_field">
!  <OVERVIEW>
!    Fill in a global array from domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_field</TT> is used to get an entire
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_field</TT>, and each will have a complete global field
!    at the end. Please note that a global array of rank 3 or higher could
!    occupy a lot of memory.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_field( domain, local, global, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="local" TYPE="MPP_TYPE_">
!    <TT>local</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="global" TYPE="MPP_TYPE_">
!    <TT>global</TT> is dimensioned on the corresponding global domain.
!  </OUT>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT> can be given the value <TT>XONLY</TT> or
!    <TT>YONLY</TT>, to specify a globalization on one axis only.
!  </IN>
! </INTERFACE>
  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
     module procedure mpp_global_field2D_c4_2d
     module procedure mpp_global_field2D_c4_3d
     module procedure mpp_global_field2D_c4_4d
     module procedure mpp_global_field2D_c4_5d
#endif
     module procedure mpp_global_field2D_i4_2d
     module procedure mpp_global_field2D_i4_3d
     module procedure mpp_global_field2D_i4_4d
     module procedure mpp_global_field2D_i4_5d
     module procedure mpp_global_field2D_l4_2d
     module procedure mpp_global_field2D_l4_3d
     module procedure mpp_global_field2D_l4_4d
     module procedure mpp_global_field2D_l4_5d
  end interface

  interface mpp_do_global_field
     module procedure mpp_do_global_field2Dnew_r8_3d
     module procedure mpp_do_global_field2Dold_r8_3d
     module procedure mpp_do_global_field2Dnew_c8_3d
     module procedure mpp_do_global_field2Dold_c8_3d
#ifndef no_8byte_integers
     module procedure mpp_do_global_field2Dnew_i8_3d
     module procedure mpp_do_global_field2Dold_i8_3d
     module procedure mpp_do_global_field2Dnew_l8_3d
     module procedure mpp_do_global_field2Dold_l8_3d
#endif
#ifndef no_4byte_reals
     module procedure mpp_do_global_field2Dnew_r4_3d
     module procedure mpp_do_global_field2Dold_r4_3d
     module procedure mpp_do_global_field2Dnew_c4_3d
     module procedure mpp_do_global_field2Dold_c4_3d
#endif
     module procedure mpp_do_global_field2Dnew_i4_3d
     module procedure mpp_do_global_field2Dold_i4_3d
     module procedure mpp_do_global_field2Dnew_l4_3d
     module procedure mpp_do_global_field2Dold_l4_3d
  end interface

! <INTERFACE NAME="mpp_global_max">
!  <OVERVIEW>
!    Global max/min of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_max</TT> is used to get the maximum value of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT> or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5. The dimension of <TT>locus</TT> must equal the rank of
!    <TT>field</TT>.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_max</TT>, and each will have the result upon exit.
!    
!    The function <TT>mpp_global_min</TT>, with an identical syntax. is
!    also available.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_global_max( domain, field, locus )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">  
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="locus" TYPE="integer" DIM="(:)">
!    <TT>locus</TT>, if present, can be used to retrieve the location of
!    the maximum (as in the <TT>MAXLOC</TT> intrinsic of f90).
!  </OUT>
! </INTERFACE>

  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
#endif
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
#endif
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

! <INTERFACE NAME="mpp_global_sum">
!  <OVERVIEW>
!    Global sum of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_sum</TT> is used to get the sum of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT>, <TT>complex</TT>, or <TT>real</TT>; of 4-byte or
!    8-byte kind; of rank up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_sum( domain, field, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT>, if present, must have the value
!    <TT>BITWISE_EXACT_SUM</TT>. This produces a sum that is guaranteed to
!    produce the identical result irrespective of how the domain is
!    decomposed. This method does the sum first along the ranks beyond 2,
!    and then calls <LINK
!    SRC="#mpp_global_field"><TT>mpp_global_field</TT></LINK> to produce a
!    global 2D array which is then summed. The default method, which is
!    considerably faster, does a local sum followed by <LINK
!    SRC="mpp.html#mpp_sum"><TT>mpp_sum</TT></LINK> across the domain
!    decomposition.
!  </IN>
!  <NOTE>
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_sum</TT>, and each will have the result upon exit.
!  </NOTE>
! </INTERFACE>

  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
#endif
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_REDUCE: get global max/min of field                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#ifndef no_4byte_reals
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#ifndef no_8byte_integers
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                   MPP_GLOBAL_SUM: global sum of field                       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#ifndef no_4byte_reals
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>
#endif

#ifndef no_8byte_integers
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>
#endif

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_FIELD: get global field from domain field           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_field.h>

#ifndef no_8byte_integers
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l8_5d
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_global_field.h>
#endif

#ifndef no_4byte_reals
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_field.h>
#endif

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l4_5d
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_global_field.h>

!****************************************************
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_r8_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_r8_3d
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_r8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_c8_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_c8_3d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_c8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#ifndef no_8byte_integers
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_i8_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_i8_3d
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_i8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_l8_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_l8_3d
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_l8_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#endif

#ifndef no_4byte_reals
#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_r4_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_r4_3d
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_r4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_c4_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_c4_3d
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_c4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#endif

#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_i4_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_i4_3d
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)                                    
#define CAFPNTR_TYPE_3D_ cafptr_i4_3d_type
#include <mpp_do_global_field_caf.h>                                    
#else                                    
#include <mpp_do_global_field_new.h>                                    
#endif

#define MPP_DO_GLOBAL_FIELD_3Dnew_ mpp_do_global_field2Dnew_l4_3d
#define MPP_DO_GLOBAL_FIELD_3Dold_ mpp_do_global_field2Dold_l4_3d
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_do_global_field_old.h>
#if defined(use_GSM)
#include <mpp_do_global_field_gsm.h>
#elif defined(use_CAF)
#define CAFPNTR_TYPE_3D_ cafptr_l4_3d_type
#include <mpp_do_global_field_caf.h>
#else
#include <mpp_do_global_field_new.h>
#endif

end module mpp_domains_reduce_mod
