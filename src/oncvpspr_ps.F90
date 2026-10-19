!
! Fully-relativistic pseudopotential stage: reads the all-electron data
! file written by oncvpspr_ae.x (see oncvpsp_main.F90 and ae_file_m)
!
 program oncvpspr_ps
 implicit none
 call oncvpsp_main(2,2)
 end program oncvpspr_ps
