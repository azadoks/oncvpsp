!
! Non-relativistic all-electron stage: writes the all-electron data
! file read by oncvpspnr_ps.x (see oncvpsp_main.F90 and ae_file_m)
!
 program oncvpspnr_ae
 implicit none
 call oncvpsp_main(0,1)
 end program oncvpspnr_ae
