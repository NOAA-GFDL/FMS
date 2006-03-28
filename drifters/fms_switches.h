#if !defined(use_libMPI)
#define _SERIAL
!DEC$ MESSAGE:'Compiling in serial mode'
#else
#undef _SERIAL
!DEC$ MESSAGE:'Compiling in MPI mode (with or without MPP) '
#endif
