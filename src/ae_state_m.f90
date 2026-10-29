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
!> The first stage of a run (oncvpsp_run): the all-electron atom of the
!> reference configuration on its log mesh, solved (aeatom) or read from
!> the all-electron data file (ae_file_m), and the all-electron halves of
!> the test configurations, solved as they are wanted (run_config_ae)
!> or read. Everything the later stages take from the all-electron atom
!> is in ae_state_t, and they do not change it, so that one state serves
!> any number of pseudopotentials of the same atom and mesh.

 module ae_state_m
 use constants_m, only: dp, mxstate, mxkap, mxcnf
 use adams_m, only: adams_t
 use config_ae_m, only: ae_config_t, ae_cache_t
 use failure_m, only: warnings_t
 implicit none
 private
 public :: ae_state_t, solve_ae

!> the all-electron atom of a run
 type :: ae_state_t
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
   integer :: irel=0
!> number of kappas per l: 1, or 2 for irel=2
   integer :: nkap=0
!> speed of light (Hartree atomic units)
   real(dp) :: c_light=0.0d0
!> size of the log radial mesh
   integer :: mmax=0
!> log of the mesh ratio
   real(dp) :: al=0.0d0
!> log radial mesh (bohr)
   real(dp), allocatable :: rr(:)
!> the mesh as the output echoes it (make_log_mesh)
   real(dp) :: log_mesh(4)=0.0d0
!> order and coefficients of the Adams-Moulton steps
   type(adams_t) :: adams
!> aeatom's iteration count
   integer :: it=0
!> total energy (Ha)
   real(dp) :: etot=0.0d0
!> its terms (ae_energies)
   real(dp) :: eterm(4)=0.0d0
!> residual of the screening potential (ae_residual)
   real(dp) :: res=0.0d0
!> eigenvalues of the core and valence states (Ha), per kappa; for
!> irel=2, the l=0 states' second column equals the first
   real(dp) :: ea(mxstate,mxkap)=0.0d0
!> core charge density
   real(dp), allocatable :: rhoc(:)
!> valence charge density
   real(dp), allocatable :: rho(:)
!> all-electron potential (Ha)
   real(dp), allocatable :: vfull(:)
!> the all-electron halves of the configurations, the reference first
!> (run_config_ae; read with the file)
   type(ae_config_t) :: aec(mxcnf)
!> the all-electron atoms solved for the configurations
   type(ae_cache_t) :: cache
!> the warnings of solve_ae
   type(warnings_t) :: warn
 end type ae_state_t

 contains

 subroutine solve_ae(irel,nkap,istage,inp,cnf,ae_filename,ae,ierr,iout,ifail)

!! the all-electron atom of the input inp: its mesh, and the reference
!! configuration solved by aeatom or, with istage=2, read with the test
!! configurations from ae_filename; the description of a Libxc functional
!! to iout before it. An aeatom that does not converge is its ierr, which
!! the caller reports
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input values (input_m), checked (check_data)
!! @param cnf  the reference and test configurations (input_m)
!! @param ae_filename  the all-electron data file of istage 1 and 2
!! @param ae  the all-electron atom
!! @param ierr  0 on success; the non-zero values differ by routine:
!!              those of aeatom (failure_m)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use input_m, only: input_t, config_t
 use adams_m, only: adams_init
 use relativity_m, only: c_light_of
 use main_setup_m, only: make_log_mesh
 use exc_libxc_m, only: exc_libxc_info
 use ae_file_m, only: read_ae_file
 use config_ae_m, only: aeatom_cached
 use failure_m, only: stage_ae

!Input variables
 integer, intent(in) :: irel,nkap,istage
 type(input_t), intent(in) :: inp
 character(len=*), intent(in) :: ae_filename
 integer, intent(in) :: iout

!Input/output variables
 type(config_t), intent(inout) :: cnf(mxcnf)

!Output variables
 type(ae_state_t), intent(out) :: ae
 integer, intent(out) :: ierr,ifail

!Local variables
 integer :: ii

 ierr=0
 ifail=0
 ae%warn%stage=stage_ae
 ae%irel=irel
 ae%nkap=nkap
 ae%c_light=c_light_of(inp%c_light)
 call make_log_mesh(irel,inp%zz,inp%amesh,inp%rr1,inp%rmax,inp%adams_order, &
&                   ae%al,ae%mmax,ae%rr,ae%log_mesh)
 ae%adams=adams_init(inp%adams_order)
 allocate(ae%rhoc(ae%mmax),ae%rho(ae%mmax),ae%vfull(ae%mmax))

 if(inp%iexc<0) then
   call exc_libxc_info(inp%iexc,iout,ifail)
   if(ifail/=0) return
 end if
 if(istage==2) then
   call read_ae_file(ae_filename,irel,ae%c_light,nkap,inp%atsym,inp%zz, &
&                    inp%nc,inp%nv,inp%iexc,inp%na,inp%la,inp%fa,ae%mmax, &
&                    ae%rr,inp%ncnf,cnf,ae%it,ae%etot,ae%eterm,ae%res,ae%ea, &
&                    ae%rhoc,ae%rho,ae%vfull,ae%aec,iout,ifail)
   if(ifail/=0) return
 else
   call aeatom_cached(ae%cache,irel,ae%c_light,nkap,inp%na,inp%la,ae%ea, &
&                     inp%fa,inp%nc,inp%nc+inp%nv,ae%it,ae%rhoc,ae%rho,ae%rr, &
&                     ae%vfull,inp%zz,ae%mmax,ae%adams,inp%iexc,ae%etot, &
&                     ae%eterm,ae%res,ierr,ae%warn,iout,ifail)
   if(ifail/=0) return
   if(ierr/=0) return
 end if

!the s states have one kappa
 if(irel==2) then
   do ii=1,inp%nc+inp%nv
     if(inp%la(ii)==0) ae%ea(ii,2)=ae%ea(ii,1)
   end do
 end if
 return
 end subroutine solve_ae

 end module ae_state_m
