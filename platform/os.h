#ifdef __sgi
#ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!parallel machine types
#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(sgi_mipspro)
#define SGICRAY_MPP
#endif

!most compilers support Cray pointers
!if you find a compiler that doesn't, #undef this inside a suitable #ifdef
#define use_CRI_pointers

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#define SHORT_KIND 2
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#define SHORT_KIND 2
#endif

#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints at any price
#define no_8byte_integers
#define LONG_KIND 4
#endif
