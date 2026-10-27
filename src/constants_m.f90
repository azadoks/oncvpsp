!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!> Kind and physical and numerical constants shared by all routines

 module constants_m
 use, intrinsic :: iso_fortran_env, only: real64
 implicit none
 private
 public :: dp, pi, mxprj, mxstate, mxcnf, mxl, mxlin, mxkap
 public :: nrhomod, nlev, nerrlv, nerrdg, dqout, well_onset_default
 public :: diag_match
 public :: ae_tolres

!> kind of all real variables
 integer, parameter :: dp=real64

!> pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!> maximum number of projectors for each l
 integer, parameter :: mxprj=5

!> maximum number of states (core plus valence) of a configuration
 integer, parameter :: mxstate=30

!> maximum number of configurations: the reference and up to mxcnf-1 tests
 integer, parameter :: mxcnf=5

!> maximum l of the pseudopotentials and projectors (check_data): the arrays
!> of the projectors have mxl+1 channels (l+1), those of the semi-local
!> potentials mxl+2, the last for the local potential (lloc=4)
 integer, parameter :: mxl=3

!> largest lmax the input readers accept, so that they cannot overrun before
!> check_data limits lmax to mxl: the per-l input arrays (rc, ep, ncon, ...)
!> have mxlin+1 slots, rc(mxl+2) holding the radius of the local potential
 integer, parameter :: mxlin=5

!> maximum number of kappa per l (nkap): j=l+1/2, and j=l-1/2 when fully
!> relativistic
 integer, parameter :: mxkap=2

!> number of columns of the model core charge rhomod: the charge and its
!> first nrhomod-1 radial derivatives
 integer, parameter :: nrhomod=5

!> number of levels of each configuration: reference, maximally ionized and
!> excited
 integer, parameter :: nlev=3

!> number of energy-error levels of the convergence tables and profiles, from
!> 1e-2 Ha down by half decades
 integer, parameter :: nerrlv=7

!> number of the levels of cvgplt, the nerrlv of the tables and those of
!> the DIAGNOSTICS block only, down to 1e-7 Ha
 integer, parameter :: nerrdg=11

!> spacing (bohr**-1) of the q mesh of the "Energy error per electron"
!> tables, as their interpolation assumes it (eresid moves the mesh points
!> onto its own q mesh, so the actual spacing differs slightly)
 real(dp), parameter :: dqout=0.5000001_dp

!> default of well_onset, the onset of the well of the well states as a
!> multiple of rc (wellstate): rc amesh**5 on the default mesh, where it was
!> before it became an input
 real(dp), parameter :: well_onset_default=1.006_dp**5

!> radius at which run_diag matches the pseudo to the all-electron
!> functions, as a multiple of max(rc(l), rc(lloc)): the mesh point 5 past
!> irc on the default mesh, where it was before it became continuous; a
!> multiple independent of the mesh, as amesh**5 would put a first-order
!> mesh error in the diagnostic (doc/smoothness_plan.md, X-1)
 real(dp), parameter :: diag_match=1.006_dp**5

!> convergence threshold of the all-electron SCF (sratom, relatom) for the
!> residual of the screening potential (Ha) and for the change of each
!> eigenvalue (relative to max(|E|,1 Ha)): 3 times the residual's
!> round-off floor (3-4e-12 for Nd, Bi), so that the energies are those of
!> the fixed point, not of the path to it. At 1e-8 the terms of the total
!> energy depended on that path to 3e-7 Ha (or at the residual's floor
!> where a finer mesh puts it above this: ae_converged)
 real(dp), parameter :: ae_tolres=1.0e-11_dp

 end module constants_m
