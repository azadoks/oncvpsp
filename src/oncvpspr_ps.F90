!
! Fully-relativistic pseudopotential stage: reads the all-electron data
! file written by oncvpspr_ae.x (see oncvpsp_r.F90 and ae_file_m)
!
#define ONCVPSP_PS_ONLY
#define ONCVPSP_EXE 'oncvpspr_ps.x'
#include "oncvpsp_r.F90"
