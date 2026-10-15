!
! Non-relativistic all-electron stage: writes the all-electron data
! file read by oncvpspnr_ps.x (see oncvpsp.F90 and ae_file_m)
!
#define ONCVPSP_NONREL
#define ONCVPSP_AE_ONLY
#define ONCVPSP_EXE 'oncvpspnr_ae.x'
#include "oncvpsp.F90"
