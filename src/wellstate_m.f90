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

 module wellstate_m

 implicit none
 private
 public :: wellstate

 contains

 subroutine wellstate(nnin,ll,kap,rc,ep,ep_as_given,onset,confwall,nowell,rr, &
&                     vfull,vwell,uu,up,zz,mmax,adams,mch,irel,c_light, &
&                     confined, &
&                     eptgt,wmesh,iout,ifail)

!! creates quantum well to confine positive-energy state, and calculates
!! the resulting all-electron wave function. The well rises from
!! ron = onset rc as cwell x**3/(1+x**3) exp(-0.01/x), x = (r-ron)/rwell.
!! onset is the input well_onset, at least 1: the well is flat to every
!! order at ron, so it leaves the function and its derivatives at rc those
!! of the atom. Its default 1.006**5 is rc amesh**5 on the default mesh,
!! where the onset was before, and does not depend on the mesh. The well
!! is that of before times a factor that makes it C-infinity at the onset
!! and is within 5 % of 1 from x = 0.2 on, so the states it confines are
!! those of before (doc/smoothness_plan.md, R-12). A bound state of this n and
!! l localized inside 0.75 rc replaces ep by its own energy, one inside
!! 1.5 rc replaces it by 0.05, and a negative ep becomes 0.25 (the
!! behaviour of before); with ep_as_given the input ep is always used instead
!! (D-1): a bound state at or above ep is an input error, and one below ep
!! inside 1.5 rc is raised to ep by a wall as close as that takes.
!! The state is solved on the log mesh continued beyond rr(mmax) until its
!! WKB decay beyond the wall at ep reaches machine precision (well_mesh_m),
!! so the mesh's end is not a box for a wide well (numerics_audit C-19)
!!
!! @param nnin  principal quantum number of well state (input); that of the well
!!              state found, which may have more nodes (output)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param rc  core radius (bohr)
!! @param ep  target energy for well state; the well state's energy on return
!! @param ep_as_given  use ep as given (D-1), instead of the replacements above
!! @param onset  onset of the well as a multiple of rc (well_onset, >= 1)
!! @param confwall  the confined well (well_tail = "confined"): a wall
!!                  x**3 exp(-0.01/x) Ha that rises without bound, and
!!                  the state of n = nnin only. Its eigenvalues are then
!!                  nondegenerate, each a smooth function of the width,
!!                  and any ep above the atom's own state of that n has a
!!                  width; Hamann's well holds only the states below its
!!                  asymptote, ep + 0.5 Ha, and moved to more nodes where
!!                  the wall would have been too steep, which made the
!!                  state's tail, and the kernels, jump
!! @param nowell  no well: only the energy is wanted (the free tail of
!!                free_tail_m, for which main solves the state at eptgt), so
!!                return after the replacements above, with uu, up zero and
!!                vwell vfull; a bound state at or above ep is then no error
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param vwell  potential with binding well, for the ikap=2 calculation when
!!               irel=2
!! @param uu  all-electron well-bound wave function
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mch  matching mesh point for inward-outward integrations
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param confined  whether uu is confined by the well; not when it is the
!!                  localized bound state (vwell is then vfull)
!! @param eptgt  the energy the well's eigenvalue (ep on return) was found
!!               for: the input ep or its replacement above
!! @param wmesh  the confined state on the log mesh continued as far as it needs
!!               (well_mesh_m), with the well's potential; mx = 0 when uu is not
!!               confined. The state is solved on it, so it does not depend on
!!               rmax; uu, up and vwell are its first mmax points
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_convergence, fail_input, fail_not_found, failed, &
&                    msg_wellstate_confined_ep, msg_wellstate_conv, &
&                    msg_wellstate_jump
 use constants_m, only: dp, eig_floor
 use well_mesh_m, only: well_mesh_t, wm_build, wm_short
 use boundstate_m, only: boundstate
 use root_m, only: root_tol
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: ll,kap,mmax,irel
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: rr(mmax),vfull(mmax)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rc
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: onset
 logical, intent(in) :: confwall
 logical, intent(in) :: nowell
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: ep
 integer, intent(inout) :: nnin !(nnin is actually in/out)
 integer, intent(inout) :: mch

!Output variables
 real(dp), intent(out) :: uu(mmax,2),up(mmax,2),vwell(mmax)
 logical, intent(out) :: confined
 real(dp), intent(out) :: eptgt
 type(well_mesh_t), intent(out) :: wmesh
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: cwell,et,xx,rwell,rwmax,rwmin,rwscale,ron,umax
 real(dp) :: fw,fwmin,fwmax,tolw

! the well width solves E(rwell) = ep to round-off (doc/smoothness_plan.md,
! Q-9): converged at a few ulps of ep, or when the bracket of rwell has
! closed to a few ulps (E then within eps of ep, else an error)
 real(dp),parameter :: eps=1.0d-10
 integer :: ierr,iumax,itrwell,nn,nnloop,iside,nnhi
!the continued mesh, vfull and the well on it, and the state
 real(dp), allocatable :: rx(:),vx(:),vw(:),ux(:,:),upx(:,:)
 integer :: mx,mneed,mchx
 logical :: convg,bound,done

 ifail=0
 ron=rc*onset
 confined=.true.
 eptgt=ep

!check for bound state with specified ll and nnin
 call find_bound(ifail)
 if(ifail/=0) return

 call target_energy(done,ifail)
 if(ifail/=0) return
 if(done) return

!the energy the well's eigenvalue is found for
 eptgt=ep

 if(nowell) then
   vwell(:)=vfull(:)
   uu(:,:)=0.0d0
   up(:,:)=0.0d0
   write(iout,'(/a,i3,a,i3/a,f8.4/a)') &
&           '   Wellstate for l =',ll,'  n =',nnin, &
&           '     eigenvalue = ',ep, &
&           '     free tail, no well'
   return
 end if

!scale factor to decrease well width until trial energy exceeds target
 rwscale=0.8d0

!asymptotic potential value to give target bound state a "binding" energy
!of 0.5 Ha; the confined wall is 1 Ha at x = 1 and has no asymptote
 cwell=ep+0.5d0
 if(confwall) cwell=1.0d0
 nnhi=10
 if(confwall) nnhi=0

!loop which will increment number of nodes before well wall gets too steep
!(the confined well keeps the n it is given)
 do nnloop=0,nnhi
   nn=nnin+nnloop

   rwell=8.0d0*rc
   rwmax=0.0d0
   rwmin=0.0d0
   fwmin=0.0d0
   fwmax=0.0d0
   iside=0
   tolw=root_tol(ep,eig_floor)

   convg=.false.

   do itrwell=1,100

!create well potential, on the log mesh continued until a state at ep
!decays in it to machine precision (well_mesh_m)
     call well_potential(ifail)
     if(ifail/=0) return

!find bound state in well
     if(allocated(ux)) deallocate(ux,upx)
     allocate(ux(mx,2),upx(mx,2))
     ux(:,:)=0.0d0
     upx(:,:)=0.0d0
     et=ep
     call boundstate(irel,c_light,nn,ll,kap,ierr,et,rx,zz,vw,ux,upx,mx,adams, &
&                    mchx,iout,ifail)
     if(ifail/=0) return
! a trial well whose state did not converge (ierr/=0: a box state, ierr=10,
! can stand in for it) only steers the bracketing; the state taken must
! have converged (checked below)

     fw=et-ep
     if(abs(fw)<=tolw) then
       ep=et
       convg=.true.
       exit
     end if

!Illinois regula falsi after proper rwell has been bracketed: E falls as
!the well widens, fwmin>0 at rwmin and fwmax<0 at rwmax; the end kept
!twice in a row has its value halved
     if(rwmin>0.0d0 .and. rwmax>0.0d0) then
       if(fw>0.0d0) then
         rwmin=rwell
         fwmin=fw
         if(iside==1) fwmax=0.5d0*fwmax
         iside=1
       else
         rwmax=rwell
         fwmax=fw
         if(iside==-1) fwmin=0.5d0*fwmin
         iside=-1
       end if
       if(rwmax-rwmin<=root_tol(rwmax,0.0d0)) then
         if(abs(fw)>eps) then
           write(iout,'(/a,i0,a,i0,a,es12.4)') &
&                'wellstate: ERROR E(rwell) jumps across ep, n=',nn, &
&                ' l=',ll,' E-ep=',fw
           ifail=failed(fail_not_found,msg_wellstate_jump)
           return
         end if
         ep=et
         convg=.true.
         exit
       end if
       rwell=(rwmin*fwmax - rwmax*fwmin)/(fwmax - fwmin)
       cycle
     end if

     if(et>ep) then
       rwmin=rwell
       fwmin=fw
       if(rwmax>0.0d0) then
         rwell=0.5d0*(rwmax+rwmin)
       else
         rwell=rwell/rwscale
       end if
       cycle
     end if

     if(et<ep) then
       rwmax=rwell
       fwmax=fw
       if(rwmin>0.0d0) then
         rwell=0.5d0*(rwmax+rwmin)
       else
         rwell=rwell*rwscale
! a wall closer than this is too steep for a state of this n, unless the
! state is bound without the well: it then only needs raising; the
! confined well's energy rises without bound as it narrows
         if(confwall) then
           if(rwell<1.0d-2*rc) exit
         else if(.not. bound .and. rwell<max(1.0d0,0.5d0*rc)) then
           exit
         end if
       end if
       cycle
     end if

   end do !itrwell

   if(convg) then
     if(ierr/=0) then
       write(iout,'(a,i3,a,i3)') &
&            'wellstate: ERROR lschfb convergence error, n=',nn,' l=',ll
       ifail=failed(fail_convergence,msg_wellstate_conv)
       return
     end if
     exit
   end if
 end do !nnloop

 if(.not. convg) then
   if(irel==2) then
     write(iout,'(a,a,i3,a,i3,a,i3,a,f8.4)') &
&     'ERROR wellstate_r: well potential iteration failed', &
&     ' to converge, n=',nn,' l=',ll,' kap=', kap,' ep=',ep
   else
     write(iout,'(a,a,i3,a,i3,a,f8.4)') &
&     'ERROR wellstate: well potential iteration failed', &
&     ' to converge, n=',nn,' l=',ll,' ep=',ep
   end if
   if(confwall) then
     write(iout,'(/a)') &
&          'wellstate: ERROR with the confined well, ep must be '// &
&         'above the atom''s own state of this n: raise ep or lower well_n'
   end if
   ifail=failed(fail_input,msg_wellstate_confined_ep)
   return
 end if

!the state and the well on the continued mesh, and their first mmax points
 wmesh%mx=mx
 wmesh%rr=rx
 wmesh%vv=vw
 wmesh%uu=ux
 wmesh%up=upx
 vwell(:)=vw(1:mmax)
 uu(:,:)=ux(1:mmax,:)
 up(:,:)=upx(1:mmax,:)
 mch=min(mchx,mmax)

 call report()

!carry over node count
 nnin=nn

 return

 contains

 subroutine find_bound(ifail)

!! the bound state of nnin and ll in vfull, if any (ierr 0, energy et), and
!! whether it is localized inside 1.5 rc (bound)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ii

 ifail=0
 uu(:,:)=0.0d0 ; up(:,:)=0.0d0
 et=-0.1d0
 call boundstate(irel,c_light,nnin,ll,kap,ierr,et,rr,zz,vfull,uu,up,mmax, &
&                adams,mch, &
&                iout,ifail)
 if(ifail/=0) return

!if bound state is found, check its localization
 bound=.false.
 if(ierr==0) then
   umax=0.0d0
   iumax=1
   do ii=mmax,1,-1
     if(abs(uu(ii,1))>=umax) then
       umax=abs(uu(ii,1))
     else
       iumax=ii+1
       exit
     end if
   end do
   bound=rr(iumax)<1.5d0*rc
 end if
 return
 end subroutine find_bound

 subroutine target_energy(done,ifail)

!! the energy the well is built for: the input ep (ep_as_given), or its
!! replacements above; done when the localized bound state replaces it,
!! with no well (uu, up its own, vwell vfull)

 use failure_m, only: fail_construction, failed, msg_wellstate_bound_below

!Output variables
 logical, intent(out) :: done
 integer, intent(out) :: ifail

 ifail=0
 done=.false.
 if(ep_as_given) then
!a bound state at or above ep cannot be raised to ep by the well: an input
!error; one below ep and localized inside 1.5 rc only needs raising, by a
!wall as close as that takes
   if(ierr==0 .and. et>=ep .and. .not. nowell) then
     write(iout,'(/a,i0,a,i0,a,f12.6/a,f12.6)') &
&          'wellstate: ERROR a bound state exists for n=',nnin,' l=',ll, &
&          ' at energy',et,'wellstate: ERROR above the requested ep=',ep
     ifail=failed(fail_construction,msg_wellstate_bound_below)
     return
   end if
 else
   bound=.false.
   if(ierr==0) then
!if bound state is localized compared to rc, use it and its energy
     if(rr(iumax)<0.75d0*rc) then
       ep=et
       write(iout,'(/a,i2,a,i2/a,f12.8/a)') &
&            'WARNING wellstate: localized bound state found for n=', &
&             nnin,' l=', ll,'WARNING this state with energy',ep, &
&            'WARNING will be used for the first projector'

!make sure vwell exists for ikap=2 calculation in this case
       vwell(:)=vfull(:)
       confined=.false.
       eptgt=ep
       done=.true.

       return

!if moderately delocalized use very low-energy scattering state
     else if(rr(iumax)<1.5d0*rc) then

       ep=0.05d0
       write(iout,'(/2a,i2,a,i2/a,f12.8/a)') &
&            'WARNING wellstate: moderately localized bound state found', &
&            ' for n=', &
&            nnin,' l= ',ll,'WARNING scattering state with energy',ep, &
&            'WARNING will be used for the first projector'
     end if
   end if

!if specified ep is retained and is negative
   if(ep<0.0d0) then
     ep=0.25d0
     write(iout,'(/a,i2,a,i2/a,f12.8/a)') &
&           'WARNING wellstate: negative energy specified for n=', &
&           nnin,' l= ',ll, &
&           'WARNING scattering state with energy',ep, &
&           'WARNING will be used for the first projector'
   end if
 end if

 return
 end subroutine target_energy

 subroutine well_potential(ifail)

!! the well of width rwell on the log mesh continued until a state at ep
!! decays in it to machine precision (well_mesh_m): rx, vx and vw on mx
!! points

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ii

 ifail=0
 mx=mmax
 do
   call wm_build(rr,vfull,mmax,mx,rx,vx)
   if(allocated(vw)) deallocate(vw)
   allocate(vw(mx))
   vw(:)=vx(:)

!   the well starts past rc, where the matching data of the projectors are
!   taken, and rises smoothly
   do ii=1,mx
     if(rx(ii)<=ron) cycle
     xx=(rx(ii)-ron)/rwell
     if(confwall) then
       vw(ii)=vw(ii)+cwell*xx**3*exp(-0.01d0/xx)
     else
       vw(ii)=vw(ii)+cwell*xx**3/(1.0d0+xx**3)*exp(-0.01d0/xx)
     end if
   end do
   call wm_short(ll,ep,rx,vw,mx,mneed,iout,ifail)
   if(ifail/=0) return
   if(mneed==0) exit
   mx=mneed
 end do
 return
 end subroutine well_potential

 subroutine report()

!! the well state's energy and well on output

 if(confwall .and. irel==2) then
   write(iout,'(/a,i3,a,i3,a,i3/a,f8.4/a,f8.4)') &
&           '   Wellstate for l =',ll,'  n =',nn,' kap=',kap, &
&           '     eigenvalue = ',et, &
&           '     confined, wall 1 Ha at r =',ron+rwell
 else if(confwall) then
   write(iout,'(/a,i3,a,i3/a,f8.4/a,f8.4)') &
&           '   Wellstate for l =',ll,'  n =',nn, &
&           '     eigenvalue = ',et, &
&           '     confined, wall 1 Ha at r =',ron+rwell
 else if(cwell>0.0d0) then
   if(irel==2) then
     write(iout,'(/a,i3,a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',nn,' kap=',kap, &
&             '     eigenvalue = ',et, &
&             '     asymptotic potential = ',cwell, &
&             '     half-point radius =',ron+rwell
   else
     write(iout,'(/a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',nn, &
&             '     eigenvalue = ',et, &
&             '     asymptotic potential = ',cwell, &
&             '     half-point radius =',ron+rwell
   end if
 end if

 if(mx>mmax) then
   write(iout,'(a,f8.2,a)') '     mesh continued to r =',rx(mx),' bohr'
 end if
 return
 end subroutine report

 end subroutine wellstate

 end module wellstate_m
