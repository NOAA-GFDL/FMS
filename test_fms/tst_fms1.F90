! This is a sample test program for the FMS library.

      program tst_fms1
      use constants_mod
      implicit none

      print *, ''
      print *,'*** Testing FMS library...'

      ! Insert test code here.
      ! constanst_init is a no-op routine.  This is just to ensure the library
      ! can be linked into an executable
      call constants_init()

      print *,'*** SUCCESS!'
      end
