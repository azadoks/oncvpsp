!
! Fully-relativistic all-electron stage: writes the all-electron data
! file read by oncvpspr_ps.x (see oncvpsp_main.F90 and ae_file_m)
!
 program oncvpspr_ae
 implicit none
 call oncvpsp_main(2,1)
 end program oncvpspr_ae
