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
! The confined atom, oncvpsp_confine.x (doc/smoothness_plan.md,
! diagnostics item 7): the all-electron atom and the pseudopotential in a
! Wigner-Seitz sphere, at the frozen potentials of the all-electron data
! file of an _ae.x (ae_file_m) and the pseudopotential data file of a run
! with --ps-file (ps_file_m), for any relativistic treatment.
!
! For each l and j and each radius rs of the sphere, the energies of the
! band edges of the Wigner-Seitz rule: where r phi'/phi = 0 at rs (the
! bottom of a band), -l-1 (its centre) and phi = 0 (its top), the roots of
! the phases atan(r phi'/phi) - pi (nodes inside rs) of phsft_node_phase,
! which fall monotonically with the energy, so a scan brackets every root
! whatever its step and Brent's method takes it to round-off. The
! pseudopotential's phase counts instead the poles of its r phi'/phi at rs
! below the energy, the states of phi(rs) = 0 in the sphere: its wave
! function can gain a node inside the projectors (Pb p, I s with a large
! debl), but r phi'/phi falls between its poles for any Hermitian
! potential. They are the eigenvalues of ws_eigenvalues in the basis of
! that condition, each taken to the root of u(rs) beside it. The
! pseudopotential's phases are shifted by the multiple of pi by which
! they differ from the all-electron ones at most energies, so an
! all-electron root and the pseudopotential root of the same phase are
! the same band's edge; n is the all-electron state's nodes inside rs plus
! l+1. For each pseudopotential edge, the cutoff that the state needs in
! the sphere, from the variational energies of ws_basis_m: a cutoff for
! every channel at the energies of the valence bands, with no tail beyond
! rs to invent. The output, on stdout:
!   ws_r l ikap rs(bohr)
!   ws l ikap edge n e_ae e_ps            (Ha; edge 1 bottom, 2 centre,
!                                          3 top)
!   ws_ecut l ikap edge n tol ecut        (Ha: the cutoff at which the
!                                          energy error falls to tol)
!   ws_floor l ikap edge n ecut err       (the error at the largest cutoff)
!   ws_conv l ikap edge n ecut err        (--table: every basis size)
!   ws_missing l ikap edge n e_ae         (no pseudopotential root of that
!                                          phase in the scan)
!   ws_core l ikap edge n e_ae            (the same, the band of a core
!                                          state, which the
!                                          pseudopotential leaves out)
!   ws_paired_off l ikap edge n e_ps      (a pseudopotential edge of
!                                          another phase, nearer the
!                                          all-electron one than its
!                                          partner)
!   ws_extra l ikap edge e_ps             (a pseudopotential edge with no
!                                          all-electron partner: a state
!                                          the all-electron atom has not,
!                                          such as a ghost)
! for l=0..3 (or those of --channels) and each kappa of l, ikap 1 for
! j=l+1/2 and 2 for j=l-1/2.
!
! Then the atom compressed: the self-consistent all-electron atom (aeatom)
! and pseudo-atom (psatom) of the reference configuration in a smooth wall
! 1/2 k x**2 exp(-dwall/x), x = r - r0 beyond r0 (by default the largest
! core radius), for a few strengths k, each from the free atom:
!   wall_r0 r0(bohr)
!   wall k etot_ae etot_ps de_ae de_ps    (Ha, k in Ha/bohr**2; de the
!                                          change from the free atom, k=0)
!   wall_eig k n l ikap e_ae e_ps         (the valence eigenvalues, Ha)
!   wall_failed k atom                    (atom ae or ps: no
!                                          self-consistent solution)
! and last a summary of both:
!   confine_summary max_de ecut4 ecut5 n_missing n_extra n_paired_off
!                   n_wall_failed max_dde ecut4_top ecut5_top
! max_de the largest |e_ps - e_ae| of the edges below 0 (Ha), ecut4 and
! ecut5 the largest cutoffs of those edges at 1e-4 and 1e-5 Ha (Infinity
! if one does not reach it), the counts of the ws_missing, ws_extra,
! ws_paired_off and wall_failed lines, max_dde the largest
! |de_ps - de_ae| of the walls (Ha), and ecut4_top and ecut5_top those of
! the top edges alone, whose errors have no tail from the potential's
! slope at rs (ws_basis_m); NaN and -1 for the parts left out.

 program oncvpsp_confine
 use constants_m, only: dp, mxstate, mxkap
 use ps_file_m, only: ps_file_t, ps_file_free
 use cli_m, only: cli_next_arg, cli_count_list, cli_radii, cli_channels, &
&                 cli_energies
 use ae_ps_pair_m, only: load_pair
 use confine_atom_m, only: confine_atom, kdef
 use confine_summary_m, only: confine_summary_t
 use build_info_m, only: write_build_info
 implicit none
 integer, parameter :: iout=6

!> the most radii of a run
 integer, parameter :: mxq=100
!> the most strengths of a run
 integer, parameter :: mxk=20

 character(len=1024) :: ae_filename,ps_filename,arg
 integer :: iarg,irel,mmax,nq,ios,nk,nc,nv,iexc
 integer :: na(mxstate),la(mxstate)
 real(dp) :: zz,ecutmax,ee(3),rq(mxq),r0,kwall(mxk)
 real(dp) :: fa(mxstate),ea(mxstate,mxkap)
!> the speed of light of the files
 real(dp) :: c_light
 real(dp), allocatable :: rr(:),vfull(:)
 logical :: table,sphere
!> the angular momenta of the output
 logical :: chan(0:3)
 type(ps_file_t) :: pf
 type(confine_summary_t) :: sm
 integer :: ifail

 ae_filename=''
 ps_filename=''
 ee(:)=0.0d0
 r0=0.0d0
 table=.false.
 sphere=.true.
 nk=size(kdef)
 kwall(1:nk)=kdef(:)
 ecutmax=200.0d0
 chan(:)=.true.
 nq=0
 iarg=1
 do while(iarg<=command_argument_count())
   call get_command_argument(iarg,arg)
   select case(arg)
   case('-h','--help')
     call write_help()
     stop
   case('-v','--version')
     call write_build_info(6,'oncvpsp_confine.x','any relativistic treatment')
     stop
   case('-a','--ae-file')
     call cli_next_arg('oncvpsp_confine',iarg,ae_filename)
   case('-p','--ps-file')
     call cli_next_arg('oncvpsp_confine',iarg,ps_filename)
   case('-r','--radii')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     call cli_radii('oncvpsp_confine',arg,mxq,nq,rq)
   case('-l','--channels')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     call cli_channels('oncvpsp_confine',arg,chan)
   case('-e','--energies')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     call cli_energies('oncvpsp_confine',arg,ee)
   case('-c','--ecut-max')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     read(arg,*,iostat=ios) ecutmax
     if(ios/=0 .or. ecutmax<=0.0d0) then
       write(iout,'(/a)') 'oncvpsp_confine: ERROR --ecut-max must be positive'
       error stop 1
     end if
   case('-t','--table')
     table=.true.
   case('-s','--no-sphere')
     sphere=.false.
   case('-w','--walls')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     if(trim(arg)=='none') then
       nk=0
     else
       nk=cli_count_list(arg)
       if(nk>mxk) then
         write(iout,'(/a,i0,a)') 'oncvpsp_confine: ERROR more than ',mxk, &
&             ' walls'
         error stop 1
       end if
       read(arg,*,iostat=ios) kwall(1:nk)
       if(ios/=0 .or. any(kwall(1:nk)<=0.0d0)) then
         write(iout,'(/2a)') &
&              'oncvpsp_confine: ERROR --walls must be positive ', &
&              'numbers separated by commas, or none'
         error stop 1
       end if
     end if
   case('-R','--wall-radius')
     call cli_next_arg('oncvpsp_confine',iarg,arg)
     read(arg,*,iostat=ios) r0
     if(ios/=0 .or. r0<=0.0d0) then
       write(iout,'(/a)') &
&           'oncvpsp_confine: ERROR --wall-radius must be positive'
       error stop 1
     end if
   case default
     write(iout,'(/2a)') 'oncvpsp_confine: ERROR unknown argument ',trim(arg)
     error stop 1
   end select
   iarg=iarg+1
 end do
 if(len_trim(ae_filename)==0 .or. len_trim(ps_filename)==0) then
   write(iout,'(/a)') 'oncvpsp_confine: ERROR requires --ae-file and --ps-file'
   error stop 1
 end if

 call load_pair('oncvpsp_confine',trim(ae_filename),trim(ps_filename),irel, &
&               c_light,zz,mmax,rr,vfull,pf,nc=nc,nv=nv,iexc=iexc,na=na,la=la, &
&               fa=fa,ea=ea,iout=iout,ifail=ifail)
 if(ifail/=0) error stop 1
 call confine_atom(pf,irel,c_light,zz,mmax,rr,vfull,nc,nv,iexc,na,la,fa,ea, &
&                  chan,sphere,nq,rq(1:nq),ee,ecutmax,table,nk,kwall(1:nk), &
&                  r0,sm,iout,ifail)
 if(ifail/=0) error stop 1

 call ps_file_free(pf)
 deallocate(rr,vfull)
 stop

 contains

 subroutine write_help()
!! the usage: what the program does and its options
 write(iout,'(a)') 'Usage: oncvpsp_confine.x -a AE_FILE -p PS_FILE [options]'
 write(iout,'(a)') ''
 write(iout,'(a)') &
&      'A pseudopotential and its all-electron atom in a Wigner-Seitz'
 write(iout,'(a)') 'sphere of radius rs, at their frozen potentials: for each l'
 write(iout,'(a)') '(and j), the band edges of the Wigner-Seitz rule, where'
 write(iout,'(a)') 'r phi''/phi = 0 at rs (bottom), -l-1 (centre) and phi = 0'
 write(iout,'(a)') &
&      '(top), both atoms'' in the energy window; and the cutoff the'
 write(iout,'(a)') 'pseudopotential''s state needs in the sphere, from its'
 write(iout,'(a)') 'variational energy in the spherical Bessel functions of the'
 write(iout,'(a)') 'same boundary condition. The output, on stdout:'
 write(iout,'(a)') '  ws_r l ikap rs(bohr)'
 write(iout,'(a)') '  ws l ikap edge n e_ae(Ha) e_ps(Ha)'
 write(iout,'(a)') '  ws_ecut l ikap edge n tol(Ha) ecut(Ha)'
 write(iout,'(a)') '  ws_floor l ikap edge n ecut(Ha) err(Ha)'
 write(iout,'(a)') '  ws_conv l ikap edge n ecut(Ha) err(Ha)       (-t)'
 write(iout,'(a)') '  ws_missing l ikap edge n e_ae(Ha)'
 write(iout,'(a)') 'edge 1 bottom, 2 centre, 3 top; n the all-electron state''s'
 write(iout,'(a)') 'nodes inside rs plus l+1; ikap 1 for j = l+1/2 (and every l'
 write(iout,'(a)') &
&      'if not fully relativistic), 2 for j = l-1/2. ws_ecut for the'
 write(iout,'(a)') 'errors 1e-2 to 1e-6 Ha that are reached; ws_missing for an'
 write(iout,'(a)') 'all-electron edge with no pseudopotential partner.'
 write(iout,'(a)') '  ws_core l ikap edge n e_ae(Ha)'
 write(iout,'(a)') '  ws_paired_off l ikap edge n e_ps(Ha)'
 write(iout,'(a)') '  ws_extra l ikap edge e_ps(Ha)'
 write(iout,'(a)') 'ws_core for one that is a core state''s band, ws_paired_off'
 write(iout,'(a)') 'for a pseudopotential edge of another phase nearer it than'
 write(iout,'(a)') 'its partner, ws_extra for a pseudopotential edge with no'
 write(iout,'(a)') 'all-electron partner (a ghost). Then the walls:'
 write(iout,'(a)') '  wall_r0 r0(bohr)'
 write(iout,'(a)') '  wall k etot_ae etot_ps de_ae de_ps     (Ha, Ha/bohr**2)'
 write(iout,'(a)') '  wall_eig k n l ikap e_ae(Ha) e_ps(Ha)'
 write(iout,'(a)') '  wall_failed k ae|ps'
 write(iout,'(a)') 'and last:'
 write(iout,'(a)') '  confine_summary max_de ecut4 ecut5 n_missing n_extra'
 write(iout,'(a)') '                  n_paired_off n_wall_failed max_dde'
 write(iout,'(a)') '                  ecut4_top ecut5_top'
 write(iout,'(a)') 'the largest |e_ps - e_ae| and the largest cutoffs at 1e-4'
 write(iout,'(a)') 'and 1e-5 Ha of the edges below 0 (Infinity if one does not'
 write(iout,'(a)') 'reach it), the counts of those lines, the largest'
 write(iout,'(a)') '|de_ps - de_ae| of the walls, and the cutoffs of the top'
 write(iout,'(a)') &
&      'edges alone, free of the tail from the potential''s slope at'
 write(iout,'(a)') 'rs; NaN and -1 for a part left out.'
 write(iout,'(a)') ''
 write(iout,'(a)') 'Options:'
 write(iout,'(a)') '  -h, --help'
 write(iout,'(a)') '      Show this help and exit'
 write(iout,'(a)') '  -v, --version'
 write(iout,'(a)') '      Show the build and exit'
 write(iout,'(a)') '  -a, --ae-file FILE'
 write(iout,'(a)') '      The all-electron data file of an _ae.x (required)'
 write(iout,'(a)') '  -p, --ps-file FILE'
 write(iout,'(a)') '      The pseudopotential data file of a run of oncvpsp.x,'
 write(iout,'(a)') &
&      '      oncvpspnr.x, oncvpspr.x or their _ps.x with --ps-file,'
 write(iout,'(a)') '      for the same atom, treatment and mesh (required)'
 write(iout,'(a)') '  -r, --radii R1,R2,...'
 write(iout,'(a)') &
&      '      The radii of the sphere (bohr), each beyond every core'
 write(iout,'(a)') '      radius. Default 1.25, 1.5 and 2 times the largest'
 write(iout,'(a)') '  -l, --channels L1,L2,...'
 write(iout,'(a)') '      The angular momenta, from 0 to 3 (both j of each when'
 write(iout,'(a)') '      fully relativistic). Default all four'
 write(iout,'(a)') '  -e, --energies EMIN,EMAX,DE'
 write(iout,'(a)') '      The window of the edges and the step of the scan of'
 write(iout,'(a)') &
&      '      the phases (Ha). Default the input''s epsh1,epsh2, and'
 write(iout,'(a)') '      0.1'
 write(iout,'(a)') '  -c, --ecut-max ECUT'
 write(iout,'(a)') '      The largest cutoff of the basis (Ha). Default 200'
 write(iout,'(a)') '  -t, --table'
 write(iout,'(a)') '      Also print the error at every basis size'
 write(iout,'(a)') '  -s, --no-sphere'
 write(iout,'(a)') '      Leave out the Wigner-Seitz sphere'
 write(iout,'(a)') '  -w, --walls K1,K2,... | none'
 write(iout,'(a)') '      The strengths k of the wall (Ha/bohr**2), or none.'
 write(iout,'(a)') '      Default 0.01, 0.03 and 0.1'
 write(iout,'(a)') '  -R, --wall-radius R0'
 write(iout,'(a)') '      The onset of the wall (bohr). Default the largest'
 write(iout,'(a)') '      core radius'
 write(iout,'(a)') ''
 write(iout,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 end subroutine write_help

 end program oncvpsp_confine
