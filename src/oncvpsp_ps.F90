!
! Scalar-relativistic pseudopotential stage: oncvpsp.x with the
! all-electron atom read from the file written by oncvpsp_ae.x
!
 program oncvpsp_ps
 use oncvpsp_main_m, only: oncvpsp_main
 implicit none
 call oncvpsp_main(1,2)
 end program oncvpsp_ps
