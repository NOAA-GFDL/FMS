!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!! @brief unit test for the mpp_root_pe() function
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests some of the procedures in the sat_vap_pressure_mod module.

program test_sat_vap_pressure

use            fms_mod, only: fms_init, fms_end
use            mpp_mod, only: mpp_error, FATAL
use       platform_mod, only: r4_kind, r8_kind
use      constants_mod, only: RDGAS, RVGAS, TFREEZE

use sat_vapor_pres_mod, only: sat_vapor_pres_init, &
                              compute_qs, compute_mrs, &
                              lookup_es, lookup_es2, lookup_es3, &
                              TCMIN, TCMAX


implicit none

integer, parameter :: N=4
real(r8_kind), dimension(N) :: TABLE, DTABLE, TABLE2, DTABLE2, TABLE3, DTABLE3
real(r8_kind), dimension(N) :: D2TABLE, D2TABLE2, D2TABLE3

call fms_init()

call sat_vapor_pres_init()
call compute_tables()
call test_compute_qs_k_0d()
call test_mrs_k_0d()
call test_lookup_es()
call test_lookup_es2()
call test_lookup_es3()

call fms_end()

contains
  !-----------------------------------------------------------------------
  subroutine test_compute_qs_k_0d()

    implicit none

    real(kind=r4_kind) :: temp4, press4, answer4, qsat4
    real(kind=r8_kind) :: temp8, press8, answer8, qsat8
    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)

    !test 1:   press is 0.  Therefore answer should be eps=EPSILO=RDGAS/RVGAS
    temp4 = 270.0_r4_kind ; press4 = 0.0_r4_kind ; answer4=real(EPSILO,r4_kind)
    temp8 = 270.0_r8_kind ; press8 = 0.0_r8_kind ; answer8=EPSILO
    ! test r4
    call compute_qs(temp4, press4, qsat4)
    if( qsat4 .ne. answer4 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r4")
    ! test r8
    call compute_qs(temp8, press8, qsat8)
    if( qsat8 .ne. answer8 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r8")


  end subroutine test_compute_qs_k_0d
  !-----------------------------------------------------------------------
  subroutine test_mrs_k_0d()

    implicit none
    real(kind=r4_kind) :: temp4, press4, answer4, mrsat4
    real(kind=r8_kind) :: temp8, press8, answer8, mrsat8
    real(kind=r8_kind), parameter :: EPSILO=real(RDGAS,r8_kind)/real(RVGAS, r8_kind)

    !press is 0.  Therefore answer should be eps=EPSILO=RDGAS/RVGAS
    temp4 = 270.0_r4_kind ; press4 = 0.0_r4_kind ; answer4=real(EPSILO,r4_kind)
    temp8 = 270.0_r8_kind ; press8 = 0.0_r8_kind ; answer8=EPSILO
    ! test r4
    call compute_mrs(temp4, press4, mrsat4)
    if( mrsat4 .ne. answer4 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r4")
    ! test r8
    call compute_mrs(temp8, press8, mrsat8)
    if( mrsat8 .ne. answer8 )call mpp_error(FATAL,"ERROR: test_compute_qs_0d fails r8")

  end subroutine test_mrs_k_0d
  !-----------------------------------------------------------------------
  subroutine test_lookup_es

    implicit none
    real(kind=r4_kind) :: temp4, esat4, answer4
    real(kind=r8_kind) :: temp8, esat8, answer8

    real(kind=r4_kind), dimension(1) :: temp4_1d, esat4_1d, answer4_1d
    real(kind=r8_kind), dimension(1) :: temp8_1d, esat8_1d, answer8_1d

    real(kind=r4_kind), dimension(1,1) :: temp4_2d, esat4_2d, answer4_2d
    real(kind=r8_kind), dimension(1,1) :: temp8_2d, esat8_2d, answer8_2d

    real(kind=r4_kind), dimension(1,1,1) :: temp4_3d, esat4_3d, answer4_3d
    real(kind=r8_kind), dimension(1,1,1) :: temp8_3d, esat8_3d, answer8_3d


    !-----0d test-------!
    temp4 = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8 = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8 = TABLE(1)
    answer4 = real(answer8, r4_kind)
    ! test r4
    call lookup_es(temp4,esat4)
    if(esat4 .ne. answer4) then
       write(*,*) 'Expected ', answer4, ' but got ', esat4
       call mpp_error(FATAL,'ERROR: test_lookup_es_0d fails r4')
    end if
    ! test r8
    call lookup_es(temp8, esat8)
    if(esat8 .ne. answer8) then
       write(*,*) 'Expected ', answer8, ' but got ', esat8
       call mpp_error(FATAL,'ERROR:  test_lookup_es_0d fails r8')
    end if


    !-----1d test-------!
    temp4_1d(1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_1d(1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_1d = TABLE(1)
    answer4_1d = real(answer8_1d, r4_kind)
    ! test r4
    call lookup_es(temp4_1d,esat4_1d)
    if(esat4_1d(1) .ne. answer4_1d(1) ) then
       write(*,*) 'Expected ', answer4_1d(1), ' but got ', esat4_1d(1)
       call mpp_error(FATAL,'ERROR: test_lookup_es_1d fails r4')
    end if
    ! test r8
    call lookup_es(temp8_1d, esat8_1d)
    if(esat8_1d(1) .ne. answer8_1d(1) ) then
       write(*,*) 'Expected ', answer8_1d(1), ' but got ', esat8_1d(1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es_1d fails r8')
    end if


    !-----2d test-------!
    temp4_2d(1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_2d(1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_2d(1,1) = TABLE(1)
    answer4_2d(1,1) = real(answer8_2d(1,1), r4_kind)
    ! test r4
    call lookup_es(temp4_2d,esat4_2d)
    if(esat4_2d(1,1) .ne. answer4_2d(1,1)) then
       write(*,*) 'Expected ', answer4_2d(1,1), ' but got ', esat4_2d(1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es_2d fails r4')
    end if
    ! test r8
    call lookup_es(temp8_2d, esat8_2d)
    if(esat8_2d(1,1) .ne. answer8_2d(1,1)) then
       write(*,*) 'Expected ', answer8_2d(1,1), ' but got ', esat8_2d(1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es_2d fails r8')
    end if


    !-----3d test-------!
    temp4_3d(1,1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_3d(1,1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_3d(1,1,1) = TABLE(1)
    answer4_3d(1,1,1) = real(answer8_3d(1,1,1), r4_kind)
    ! test r4
    call lookup_es(temp4_3d,esat4_3d)
    if(esat4_3d(1,1,1) .ne. answer4_3d(1,1,1)) then
       write(*,*) 'Expected ', answer4_3d(1,1,1), ' but got ', esat4_3d(1,1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es_3d fails r4')
    end if
    ! test r8
    call lookup_es(temp8_3d, esat8_3d)
    if(esat8_3d(1,1,1) .ne. answer8_3d(1,1,1)) then
       write(*,*) 'Expected ', answer8_3d(1,1,1), ' but got ', esat8_3d(1,1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es_3d fails r8')
    end if

  end subroutine test_lookup_es
  !-----------------------------------------------------------------------
  subroutine test_lookup_es2

    implicit none
    real(kind=r4_kind) :: temp4, esat4, answer4
    real(kind=r8_kind) :: temp8, esat8, answer8

    real(kind=r4_kind), dimension(1) :: temp4_1d, esat4_1d, answer4_1d
    real(kind=r8_kind), dimension(1) :: temp8_1d, esat8_1d, answer8_1d

    real(kind=r4_kind), dimension(1,1) :: temp4_2d, esat4_2d, answer4_2d
    real(kind=r8_kind), dimension(1,1) :: temp8_2d, esat8_2d, answer8_2d

    real(kind=r4_kind), dimension(1,1,1) :: temp4_3d, esat4_3d, answer4_3d
    real(kind=r8_kind), dimension(1,1,1) :: temp8_3d, esat8_3d, answer8_3d


    !-----0d test-------!
    temp4 = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8 = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8 = TABLE2(1)
    answer4 = real(answer8, r4_kind)
    ! test r4
    call lookup_es2(temp4,esat4)
    if(esat4 .ne. answer4) then
       write(*,*) 'Expected ', answer4, ' but got ', esat4
       call mpp_error(FATAL,'ERROR: test_lookup_es2_0d fails r4')
    end if
    ! test r8
    call lookup_es2(temp8, esat8)
    if(esat8 .ne. answer8) then
       write(*,*) 'Expected ', answer8, ' but got ', esat8
       call mpp_error(FATAL,'ERROR:  test_lookup_es2_0d fails r8')
    end if

        !-----1d test-------!
    temp4_1d(1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_1d(1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_1d = TABLE2(1)
    answer4_1d = real(answer8_1d, r4_kind)
    ! test r4
    call lookup_es2(temp4_1d,esat4_1d)
    if(esat4_1d(1) .ne. answer4_1d(1) ) then
       write(*,*) 'Expected ', answer4_1d(1), ' but got ', esat4_1d(1)
       call mpp_error(FATAL,'ERROR: test_lookup_es2_1d fails r4')
    end if
    ! test r8
    call lookup_es2(temp8_1d, esat8_1d)
    if(esat8_1d(1) .ne. answer8_1d(1) ) then
       write(*,*) 'Expected ', answer8_1d(1), ' but got ', esat8_1d(1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es2_1d fails r8')
    end if

    !-----2d test-------!
    temp4_2d(1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_2d(1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_2d(1,1) = TABLE2(1)
    answer4_2d(1,1) = real(answer8_2d(1,1), r4_kind)
    ! test r4
    call lookup_es2(temp4_2d,esat4_2d)
    if(esat4_2d(1,1) .ne. answer4_2d(1,1)) then
       write(*,*) 'Expected ', answer4_2d(1,1), ' but got ', esat4_2d(1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es2_2d fails r4')
    end if
    ! test r8
    call lookup_es2(temp8_2d, esat8_2d)
    if(esat8_2d(1,1) .ne. answer8_2d(1,1)) then
       write(*,*) 'Expected ', answer8_2d(1,1), ' but got ', esat8_2d(1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es2_2d fails r8')
    end if


    !-----3d test-------!
    temp4_3d(1,1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_3d(1,1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_3d(1,1,1) = TABLE2(1)
    answer4_3d(1,1,1) = real(answer8_3d(1,1,1), r4_kind)
    ! test r4
    call lookup_es2(temp4_3d,esat4_3d)
    if(esat4_3d(1,1,1) .ne. answer4_3d(1,1,1)) then
       write(*,*) 'Expected ', answer4_3d(1,1,1), ' but got ', esat4_3d(1,1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es2_3d fails r4')
    end if
    ! test r8
    call lookup_es2(temp8_3d, esat8_3d)
    if(esat8_3d(1,1,1) .ne. answer8_3d(1,1,1)) then
       write(*,*) 'Expected ', answer8_3d(1,1,1), ' but got ', esat8_3d(1,1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es2_3d fails r8')
    end if

  end subroutine test_lookup_es2
  !-----------------------------------------------------------------------
  subroutine test_lookup_es3

    implicit none
    real(kind=r4_kind) :: temp4, esat4, answer4
    real(kind=r8_kind) :: temp8, esat8, answer8

    real(kind=r4_kind), dimension(1) :: temp4_1d, esat4_1d, answer4_1d
    real(kind=r8_kind), dimension(1) :: temp8_1d, esat8_1d, answer8_1d

    real(kind=r4_kind), dimension(1,1) :: temp4_2d, esat4_2d, answer4_2d
    real(kind=r8_kind), dimension(1,1) :: temp8_2d, esat8_2d, answer8_2d

    real(kind=r4_kind), dimension(1,1,1) :: temp4_3d, esat4_3d, answer4_3d
    real(kind=r8_kind), dimension(1,1,1) :: temp8_3d, esat8_3d, answer8_3d


    !-----0d test-------!
    temp4 = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8 = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8 = TABLE3(1)
    answer4 = real(answer8, r4_kind)
    ! test r4
    call lookup_es3(temp4,esat4)
    if(esat4 .ne. answer4) then
       write(*,*) 'Expected ', answer4, ' but got ', esat4
       call mpp_error(FATAL,'ERROR: test_lookup_es3_0d fails r4')
    end if
    ! test r8
    call lookup_es3(temp8, esat8)
    if(esat8 .ne. answer8) then
       write(*,*) 'Expected ', answer8, ' but got ', esat8
       call mpp_error(FATAL,'ERROR:  test_lookup_es3_0d fails r8')
    end if


    !-----1d test-------!
    temp4_1d(1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_1d(1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_1d = TABLE3(1)
    answer4_1d = real(answer8_1d, r4_kind)
    ! test r4
    call lookup_es3(temp4_1d,esat4_1d)
    if(esat4_1d(1) .ne. answer4_1d(1) ) then
       write(*,*) 'Expected ', answer4_1d(1), ' but got ', esat4_1d(1)
       call mpp_error(FATAL,'ERROR: test_lookup_es3_1d fails r4')
    end if
    ! test r8
    call lookup_es3(temp8_1d, esat8_1d)
    if(esat8_1d(1) .ne. answer8_1d(1) ) then
       write(*,*) 'Expected ', answer8_1d(1), ' but got ', esat8_1d(1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es3_1d fails r8')
    end if


    !-----2d test-------!
    temp4_2d(1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_2d(1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_2d(1,1) = TABLE3(1)
    answer4_2d(1,1) = real(answer8_2d(1,1), r4_kind)
    ! test r4
    call lookup_es3(temp4_2d,esat4_2d)
    if(esat4_2d(1,1) .ne. answer4_2d(1,1)) then
       write(*,*) 'Expected ', answer4_2d(1,1), ' but got ', esat4_2d(1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es3_2d fails r4')
    end if
    ! test r8
    call lookup_es3(temp8_2d, esat8_2d)
    if(esat8_2d(1,1) .ne. answer8_2d(1,1)) then
       write(*,*) 'Expected ', answer8_2d(1,1), ' but got ', esat8_2d(1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es3_2d fails r8')
    end if


    !-----3d test-------!
    temp4_3d(1,1,1) = real(TCMIN,r4_kind) + real(TFREEZE,r4_kind) !tminl corresponding to TABLE(1)
    temp8_3d(1,1,1) = real(TCMIN,r8_kind) + real(TFREEZE,r8_kind) !tminl corresponding to TABLE(1)

    !get answers.  The TABLE is computed with r8_kind precision
    answer8_3d(1,1,1) = TABLE3(1)
    answer4_3d(1,1,1) = real(answer8_3d(1,1,1), r4_kind)
    ! test r4
    call lookup_es3(temp4_3d,esat4_3d)
    if(esat4_3d(1,1,1) .ne. answer4_3d(1,1,1)) then
       write(*,*) 'Expected ', answer4_3d(1,1,1), ' but got ', esat4_3d(1,1,1)
       call mpp_error(FATAL,'ERROR: test_lookup_es3_3d fails r4')
    end if
    ! test r8
    call lookup_es3(temp8_3d, esat8_3d)
    if(esat8_3d(1,1,1) .ne. answer8_3d(1,1,1)) then
       write(*,*) 'Expected ', answer8_3d(1,1,1), ' but got ', esat8_3d(1,1,1)
       call mpp_error(FATAL,'ERROR:  test_lookup_es3_3d fails r8')
    end if

  end subroutine test_lookup_es3
  !-----------------------------------------------------------------------
  !subroutine test_lookup_des

  !-----------------------------------------------------------------------
  !subroutine test_lookup_des2
  !subroutine test_lookup_des3
  !subroutine test_lookup_es_des
  !subroutine test_lookup_es2_des2
  !subroutine test_lookup_es3_des3



  !-----------------------------------------------------------------------
  subroutine compute_tables

    ! increment used to generate derivative table
    real(kind=r8_kind), dimension(3) :: tem, es
    real(kind=r8_kind) :: dtres, tminl, dtinvl, tepsl, hdtinv, tinrc, tfact
    integer :: i, t

    integer, parameter ::esres=10


    t = (TCMAX-TCMIN)*esres
    dtres  = (real(TCMAX,r8_kind)-real(TCMIN,r8_kind))/real(t,r8_kind)
    tminl  = real(tcmin,r8_kind)+real(TFREEZE,r8_kind)  ! minimum valid temp in table
    dtinvl = 1.0_r8_kind/dtres
    tepsl  = 0.5_r8_kind*dtres
    tinrc  = 0.1_r8_kind*dtres
    tfact  = 5.0_r8_kind*dtinvl
    hdtinv = 0.5_r8_kind*dtinvl

    do i = 1, N
       tem(1) = tminl + dtres*real(i-1,r8_kind)
       tem(2) = tem(1)-tinrc
       tem(3) = tem(1)+tinrc
       es = compute_es_k (tem, real(TFREEZE,r8_kind))
       TABLE(i)  = es(1)
       DTABLE(i) = (es(3)-es(2))*tfact
    enddo
    do i = 2, N-1
       D2TABLE(i) = 0.25_r8_kind*dtinvl*(DTABLE(i+1)-DTABLE(i-1))
    enddo
    ! one-sided derivatives at boundaries
    D2TABLE(1) = 0.50_r8_kind*dtinvl*(DTABLE(2)-DTABLE(1))

    do i = 1, N
        tem(1) = tminl + dtres*real(i-1,r8_kind)
        tem(2) = tem(1)-tinrc
        tem(3) = tem(1)+tinrc
        !   pass in flag to force all values to be wrt liquid
        es = compute_es_liq_k (tem, real(TFREEZE,r8_kind))
        TABLE2(i) = es(1)
        DTABLE2(i) = (es(3)-es(2))*tfact
     enddo
     do i = 2, N-1
        D2TABLE2(i) = 0.25_r8_kind*dtinvl*(DTABLE2(i+1)-DTABLE2(i-1))
     enddo
     ! one-sided derivatives at boundaries
     D2TABLE2(1) = 0.50_r8_kind*dtinvl*(DTABLE2(2)-DTABLE2(1))

    do i = 1, N
       tem(1) = tminl + dtres*real(i-1,r8_kind)
       tem(2) = tem(1)-tinrc
       tem(3) = tem(1)+tinrc
       !   pass in flag to force all values to be wrt liquid
       es = compute_es_liq_ice_k (tem, real(TFREEZE,r8_kind))
       TABLE3(i) = es(1)
       DTABLE3(i) = (es(3)-es(2))*tfact
    enddo
    do i = 2, N-1
       D2TABLE3(i) = 0.25_r8_kind*dtinvl*(DTABLE3(i+1)-DTABLE3(i-1))
    enddo
    ! one-sided derivatives at boundaries
    D2TABLE3(1) = 0.50_r8_kind*dtinvl*(DTABLE3(2)-DTABLE3(1))

  end subroutine compute_tables
  !-----------------------------------------------------------------------
  function compute_es_k(tem, TFREEZE) result (es)
    real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
    real(kind=r8_kind) :: es(size(tem,1))

    real(kind=r8_kind)  :: x, esice, esh2o, TBASW, TBASI
    integer :: i

    real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind
    real(kind=r8_kind), parameter :: ESBASI = 610.71_r8_kind

    real(r8_kind), parameter :: one=1.0_r8_kind
    real(r8_kind), parameter :: ten=10.0_r8_kind

    TBASW = TFREEZE+100.0_r8_kind
    TBASI = TFREEZE
    do i = 1, size(tem)

       !  compute es over ice

       if (tem(i) < TBASI) then
          x = -9.09718_r8_kind*(TBASI/tem(i)-one)   &
               -3.56654_r8_kind*log10(TBASI/tem(i)) &
               +0.876793_r8_kind*(one-tem(i)/TBASI) + log10(ESBASI)
          esice =ten**(x)
       else
          esice = 0.0_r8_kind
       endif

       !  compute es over water greater than -20 c.
       !  values over 100 c may not be valid
       !  see smithsonian meteorological tables page 350.

       if (tem(i) > -20.0_r8_kind+TBASI) then
          x = -7.90298_r8_kind*(TBASW/tem(i)-one)   &
               +5.02808_r8_kind*log10(TBASW/tem(i)) &
               -1.3816d-07*(ten**((one-tem(i)/TBASW)*11.344_r8_kind)-one) &
               +8.1328d-03*(ten**((TBASW/tem(i)-one)*(-3.49149_r8_kind))-one) &
               +log10(ESBASW)
          esh2o = ten**(x)
       else
          esh2o = 0.0_r8_kind
       endif

       !  derive blended es over ice and supercooled water between -20c and 0c

       if (tem(i) <= -20.0_r8_kind+TBASI) then
          es(i) = esice
       else if (tem(i) >= TBASI) then
          es(i) = esh2o
       else
          es(i) = 0.05_r8_kind*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.0_r8_kind)*esh2o)
       endif

    enddo

  end function compute_es_k
!-----------------------------------------------------------------------
 function compute_es_liq_k(tem, TFREEZE) result (es)
 real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
 real(kind=r8_kind) :: es(size(tem,1))

 real(kind=r8_kind) :: x, esh2o, TBASW
 integer :: i

 real(kind=r8_kind), parameter :: one=1.0_r8_kind
 real(kind=r8_kind), parameter :: ten=10.0_r8_kind
 real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind


   TBASW = TFREEZE+100.0_r8_kind

   do i = 1, size(tem)
!  compute es over water for all temps.
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.
         x = -7.90298_r8_kind*(TBASW/tem(i)-one) &
             +5.02808_r8_kind*log10(TBASW/tem(i)) &
             -real(1.3816d-07,r8_kind)*(ten**((one-tem(i)/TBASW)*11.344_r8_kind)-one)    &
             +real(8.1328d-03,r8_kind)*(ten**((TBASW/tem(i)-one)*(-3.49149_r8_kind))-one)&
             +log10(ESBASW)
         esh2o = ten**(x)
         es(i) = esh2o

   enddo

 end function compute_es_liq_k
 !-----------------------------------------------------------------------
 function compute_es_liq_ice_k(tem, TFREEZE) result (es)

 real(kind=r8_kind), intent(in) :: tem(:), TFREEZE
 real(kind=r8_kind) :: es(size(tem,1))

 real(kind=r8_kind)    :: x, TBASW, TBASI
 integer :: i

 real(kind=r8_kind), parameter :: ESBASW = 101324.60_r8_kind
 real(kind=r8_kind), parameter :: ESBASI = 610.71_r8_kind
 real(kind=r8_kind), parameter :: one=  1.0_r8_kind
 real(kind=r8_kind), parameter :: ten= 10.0_r8_kind

   TBASW = TFREEZE+100.0_r8_kind
   TBASI = TFREEZE

   do i = 1, size(tem)

     if (tem(i) < TBASI) then

!  compute es over ice

         x = -9.09718_r8_kind*(TBASI/tem(i)-one)   &
             -3.56654_r8_kind*log10(TBASI/tem(i)) &
             +0.876793_r8_kind*(one-tem(i)/TBASI) + log10(ESBASI)
         es(i) =ten**(x)
     else

!  compute es over water
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

          x = -7.90298_r8_kind*(TBASW/tem(i)-one) &
              +5.02808_r8_kind*log10(TBASW/tem(i)) &
              -real(1.3816d-07,r8_kind)*(ten**((one-tem(i)/TBASW)*11.344_r8_kind)-one)      &
              +real(8.1328d-03,r8_kind)*(ten**((TBASW/tem(i)-one)*(-3.49149_r8_kind))-one) &
             +log10(ESBASW)
         es(i) = ten**(x)
     endif
   enddo

 end function compute_es_liq_ice_k
 !-----------------------------------------------------------------------
end program test_sat_vap_pressure
