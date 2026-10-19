!
! Scalar-relativistic all-electron stage: writes the all-electron data
! file read by oncvpsp_ps.x (see oncvpsp_main.F90 and ae_file_m)
!
 program oncvpsp_ae
 implicit none
 call oncvpsp_main(1,1)
 end program oncvpsp_ae
