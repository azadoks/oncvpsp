!
! Non-relativistic pseudopotential stage: oncvpspnr.x with the
! all-electron atom read from the file written by oncvpspnr_ae.x
!
#define ONCVPSP_NONREL
#define ONCVPSP_PS_ONLY
#define ONCVPSP_EXE 'oncvpspnr_ps.x'
#include "oncvpsp.F90"
