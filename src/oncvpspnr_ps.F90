!
! Non-relativistic pseudopotential stage: oncvpspnr.x with the
! all-electron atom read from the file written by oncvpspnr_ae.x
!
 program oncvpspnr_ps
 use oncvpsp_main_m, only: oncvpsp_main
 implicit none
 call oncvpsp_main(0,2)
 end program oncvpspnr_ps
