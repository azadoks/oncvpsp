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

 use constants_m, only: dp, mxstate, mxcnf, mxlin, mxkap, well_onset_default
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

!> The values of the input file. Each component is defaulted to the value
!> an input that leaves it out gets, so that a reader (whose inp is
!> intent(out)) sets only what it reads
 type :: input_t
! [Atom and reference configuration]
!> Atomic symbol
   character(len=2) :: atsym=' '
!> Atomic number
   real(dp) :: zz=0.0d0
!> Number of core states
   integer :: nc=0
!> Number of valence states
   integer :: nv=0
!> Exchange-correlation functional code
   integer :: iexc=0
!> Pseudopotential output file format flag
   character(len=4) :: psfile=' '
!> Principal quantum number array
   integer :: na(mxstate)=0
!> Orbital angular momentum quantum number array
   integer :: la(mxstate)=0
!> Occupation number array
   real(dp) :: fa(mxstate)=0.0d0
! [Pseudopotential and optimization]
!> Maximum angular momentum
   integer :: lmax=0
!> Whether the well states take ep as given (wellstate; optional on the
!> lmax line, 0 or 1)
   logical :: ep_as_given=.false.
!> Onset of the well of the well states as a multiple of rc (optional on
!> the lmax line, 0 or less for the default)
   real(dp) :: well_onset=well_onset_default
!> The tail of the well states (optional on the lmax line): 0 "well",
!> Hamann's well; 1 "free", the free tail (free_tail_m); 2 "confined",
!> a wall that rises without bound, with the state's n as given
   integer :: well_tail=0
!> The regularization of the pseudo-wave optimization (wave_regularization,
!> optional on the lmax line after well_tail, 0 for none)
   real(dp) :: wave_regularization=0.0d0
!> The normalization of the well states (optional on the lmax line, after
!> wave_regularization): 0 that of well_tail (all space for "well", inside
!> rc for "free" and "confined"), 1 all space, 2 inside rc
   integer :: well_norm=0
!> Core radii for pseudopotentials
   real(dp) :: rc(mxlin+1)=0.0d0
!> Pseudopotential energies, the input value for both j=l+1/2 and l-1/2
!> (only the fully-relativistic code uses the second)
   real(dp) :: ep(mxlin+1, mxkap)=0.0d0
!> Number of matching constraints for each pseudopotential
   integer :: ncon(mxlin+1)=0
!> Number of basis functions for each pseudopotential
   integer :: nbas(mxlin+1)=0
!> Residual-energy cutoff wave vector for each pseudopotential
   real(dp) :: qcut(mxlin+1)=0.0d0
! [Local potential]
!> Angular momentum used for local potential
   integer :: lloc=0
!> Type of polynomial continuation of the lloc=4 local potential
   integer :: lpopt=0
!> Local potential offset at origin
   real(dp) :: dvloc0=0.0d0
! [Vanderbilt-Kleinman-Bylander projectors]
!> Number of projectors for each angular momentum
   integer :: nproj(mxlin+1)=0
!> Energy added to ep for the second projector (Ha)
   real(dp) :: debl(mxlin+1)=0.0d0
!> n of the first well state of each l with well_tail 2 (optional on the
!> projector line, 0 for the default)
   integer :: well_n(mxlin+1)=0
!> Weight of the first pseudo-wave of each l with no occupied state in the
!> cutoff of the pseudo-waves (optional on the projector line after well_n,
!> negative for the default 2(2l+1))
   real(dp) :: ecut_weight(mxlin+1)=-1.0d0
! [Model core charge]
!> Model core charge flag
   integer :: icmod=0
!> Scaling factor for core charge
   real(dp) :: fcfact=0.0d0
!> Scaling factor range minimum for optimization
   real(dp) :: fcfact_min=1.5d0
!> Scaling factor range maximum for optimization
   real(dp) :: fcfact_max=6.0d0
!> Scaling factor step size for optimization
   real(dp) :: fcfact_step=0.5d0
!> Core charge width factor
   real(dp) :: rcfact=0.0d0
!> Share of the cutoff the model core must fit (icmod 6, 7)
   real(dp) :: slack=1.0d0
!> The cutoff of the pseudo-waves (icmod 6, 7; optional after slack): 0
!> "max", the largest over the ground-state waves (those of occupied
!> states and each empty channel's first), 1 "weighted", their weighted
!> mean, 2 "max_all", the largest over all waves
   integer :: budget_cutoff=0
!> Core charge width factor minimum for optimization
   real(dp) :: rcfact_min=1.0d0
!> Core charge width factor maximum for optimization
   real(dp) :: rcfact_max=1.9d0
!> Core charge width factor step size for optimization
   real(dp) :: rcfact_step=0.1d0
! [Log derivative analysis]
!> Lower bound of energy range for log derivative analysis
   real(dp) :: epsh1=0.0d0
!> Upper bound of energy range for log derivative analysis
   real(dp) :: epsh2=0.0d0
!> Energy step size for log derivative analysis
   real(dp) :: depsh=0.0d0
!> Radius for log derivative analysis
   real(dp) :: rxpsh=-1.0d0
! [Output grid]
!> Maximum radius for output grid
   real(dp) :: rlmax=0.0d0
!> Grid spacing for output grid
   real(dp) :: drl=0.0d0
!> Most points of the log mesh of the UPF file (0 for the linear output grid)
   integer :: upf_mesh=0
! [Test configurations]
!> Number of test configurations
   integer :: ncnf=0
!> The test configurations, cnf(2:ncnf+1) as read (their valence states).
!> cnf(1) is not set: the configurations are numbered from the reference
!> (na, la, fa), 1, and oncvpsp_run makes the list of all of them, the
!> core states in each
   type(config_t) :: cnf(mxcnf)
! [Log mesh]
!> Ratio of successive points of the log radial mesh
!> (0 or less for the default)
   real(dp) :: amesh=-1.0d0
!> First point of the log radial mesh (0 or less for the default)
   real(dp) :: rr1=-1.0d0
!> Outer radius of the log radial mesh (0 or less for the default)
   real(dp) :: rmax=-1.0d0
! [Atom, optional]
!> Speed of light (0 or less for the default)
   real(dp) :: c_light=-1.0d0
 end type input_t

 end module input_m
