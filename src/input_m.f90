!
! Copyright (c) 2026 by Austin Zadoks
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
!> The input of a calculation, as read (read_input_text, read_input_toml):
!> the values of the input file, or their defaults
 module input_m

 use constants_m, only: dp, mxstate, mxcnf, mxlin, mxkap
 implicit none
 private
 public :: input_t, config_t

!> A configuration of the atom: its states and their occupations
 type :: config_t
!> Number of valence states
   integer :: nv=0
!> principal quantum numbers of the states, the core states first (as
!> many as the reference's nc)
   integer :: na(mxstate)=0
!> angular momenta of the states, ordered as na
   integer :: la(mxstate)=0
!> occupancies of the states, ordered as na
   real(dp) :: fa(mxstate)=0.0d0
 end type config_t

!> The values of the input file
 type :: input_t
! [Atom and reference configuration]
!> Atomic symbol
   character(len=2) :: atsym
!> Atomic number
   real(dp) :: zz
!> Number of core states
   integer :: nc
!> Number of valence states
   integer :: nv
!> Exchange-correlation functional code
   integer :: iexc
!> Pseudopotential output file format flag
   character(len=4) :: psfile
!> Principal quantum number array
   integer :: na(mxstate)
!> Orbital angular momentum quantum number array
   integer :: la(mxstate)
!> Occupation number array
   real(dp) :: fa(mxstate)
! [Pseudopotential and optimization]
!> Maximum angular momentum
   integer :: lmax
!> Whether the well states take ep as given (wellstate; optional on the
!> lmax line, 0 or 1)
   logical :: ep_as_given
!> Onset of the well of the well states as a multiple of rc (optional on
!> the lmax line, 0 or less for the default)
   real(dp) :: well_onset
!> The tail of the well states (optional on the lmax line): 0 "well",
!> Hamann's well; 1 "free", the free tail (free_tail_m); 2 "confined",
!> a wall that rises without bound, with the state's n as given
   integer :: well_tail
!> The regularization of the pseudo-wave optimization (wave_regularization,
!> optional on the lmax line after well_tail, 0 for none)
   real(dp) :: wave_regularization
!> The normalization of the well states (optional on the lmax line, after
!> wave_regularization): 0 that of well_tail (all space for "well", inside
!> rc for "free" and "confined"), 1 all space, 2 inside rc
   integer :: well_norm
!> Core radii for pseudopotentials
   real(dp) :: rc(mxlin+1)
!> Pseudopotential energies, the input value for both j=l+1/2 and l-1/2
!> (only the fully-relativistic code uses the second)
   real(dp) :: ep(mxlin+1, mxkap)
!> Number of matching constraints for each pseudopotential
   integer :: ncon(mxlin+1)
!> Number of basis functions for each pseudopotential
   integer :: nbas(mxlin+1)
!> Residual-energy cutoff wave vector for each pseudopotential
   real(dp) :: qcut(mxlin+1)
! [Local potential]
!> Angular momentum used for local potential
   integer :: lloc
!> Type of polynomial continuation of the lloc=4 local potential
   integer :: lpopt
!> Local potential offset at origin
   real(dp) :: dvloc0
! [Vanderbilt-Kleinman-Bylander projectors]
!> Number of projectors for each angular momentum
   integer :: nproj(mxlin+1)
!> Energy added to ep for the second projector (Ha)
   real(dp) :: debl(mxlin+1)
!> n of the first well state of each l with well_tail 2 (optional on the
!> projector line, 0 for the default)
   integer :: well_n(mxlin+1)
!> Weight of the first pseudo-wave of each l with no occupied state in the
!> cutoff of the pseudo-waves (optional on the projector line after well_n,
!> negative for the default 2(2l+1))
   real(dp) :: ecut_weight(mxlin+1)
! [Model core charge]
!> Model core charge flag
   integer :: icmod
!> Scaling factor for core charge
   real(dp) :: fcfact
!> Scaling factor range minimum for optimization
   real(dp) :: fcfact_min
!> Scaling factor range maximum for optimization
   real(dp) :: fcfact_max
!> Scaling factor step size for optimization
   real(dp) :: fcfact_step
!> Core charge width factor
   real(dp) :: rcfact
!> Share of the cutoff the model core must fit (icmod 6, 7)
   real(dp) :: slack
!> The cutoff of the pseudo-waves (icmod 6, 7; optional after slack): 0
!> "max", the largest over the ground-state waves (those of occupied
!> states and each empty channel's first), 1 "weighted", their weighted
!> mean, 2 "max_all", the largest over all waves
   integer :: budget_cutoff
!> Core charge width factor minimum for optimization
   real(dp) :: rcfact_min
!> Core charge width factor maximum for optimization
   real(dp) :: rcfact_max
!> Core charge width factor step size for optimization
   real(dp) :: rcfact_step
! [Log derivative analysis]
!> Lower bound of energy range for log derivative analysis
   real(dp) :: epsh1
!> Upper bound of energy range for log derivative analysis
   real(dp) :: epsh2
!> Energy step size for log derivative analysis
   real(dp) :: depsh
!> Radius for log derivative analysis
   real(dp) :: rxpsh
! [Output grid]
!> Maximum radius for output grid
   real(dp) :: rlmax
!> Grid spacing for output grid
   real(dp) :: drl
!> Most points of the log mesh of the UPF file (0 for the linear output grid)
   integer :: upf_mesh
! [Test configurations]
!> Number of test configurations
   integer :: ncnf
!> The test configurations, cnf(2:ncnf+1) as read (their valence states);
!> oncvpsp_run puts the reference configuration in cnf(1) and the core
!> states in each
   type(config_t) :: cnf(mxcnf)
! [Log mesh]
!> Ratio of successive points of the log radial mesh
!> (0 or less for the default)
   real(dp) :: amesh
!> First point of the log radial mesh (0 or less for the default)
   real(dp) :: rr1
!> Outer radius of the log radial mesh (0 or less for the default)
   real(dp) :: rmax
! [Atom, optional]
!> Speed of light (0 or less for the default)
   real(dp) :: c_light
 end type input_t

 end module input_m
