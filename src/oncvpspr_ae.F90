!
! Fully-relativistic all-electron stage: writes the all-electron data
! file read by oncvpspr_ps.x (see oncvpsp_main_m.F90 and ae_file_m)
!
 program oncvpspr_ae
 use oncvpsp_main_m, only: oncvpsp_main
 implicit none
 call oncvpsp_main(2,1)
 end program oncvpspr_ae
