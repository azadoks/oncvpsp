!
! Fully-relativistic all-electron stage: writes the all-electron data
! file read by oncvpspr_ps.x (see oncvpsp_r.F90 and ae_file_m)
!
#define ONCVPSP_AE_ONLY
#define ONCVPSP_EXE 'oncvpspr_ae.x'
#include "oncvpsp_r.F90"
