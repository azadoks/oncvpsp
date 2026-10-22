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
 subroutine wellstate(nnin,ll,kap,rc,ep,epgiven,onset,nowell,rr,vfull,vwell, &
&                     uu,up,zz,mmax,mch,irel,confined,eptgt)

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
!! behaviour of before); with epgiven the input ep is always used instead
!! (D-1): a bound state at or above ep is an input error, and one below ep
!! inside 1.5 rc is raised to ep by a wall as close as that takes
!!
!! @param nnin  principal quantum number of well state (input); that of the well
!!              state found, which may have more nodes (output)
!! @param ll  angular momentum
!! @param kap  l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param rc  core radius (bohr)
!! @param ep  target energy for well state; the well state's energy on return
!! @param epgiven  use ep as given (D-1), instead of the replacements above
!! @param onset  onset of the well as a multiple of rc (well_onset, >= 1)
!! @param nowell  no well: only the energy is wanted (the free tail of
!!                free_tail_m, for which main solves the state at eptgt), so
!!                return after the replacements above, with uu, up zero and
!!                vwell vfull; a bound state at or above ep is then no error
!! @param rr  log radial mesh
!! @param vfull  all-electron potential
!! @param vwell  potential with binding well, for the ikap=2 calculation when
!!               irel=2
!! @param uu  all-electron well-bound wave function
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param zz  atomic number
!! @param mmax  size of log radial mesh
!! @param mch  matching mesh point for inward-outward integrations
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param confined  whether uu is confined by the well; not when it is the
!!                  localized bound state (vwell is then vfull)
!! @param eptgt  the energy the well's eigenvalue (ep on return) was found
!!               for: the input ep or its replacement above

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vfull(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: ll,kap,mmax,irel
 real(dp), intent(in) :: rc
 logical, intent(in) :: epgiven
 real(dp), intent(in) :: onset
 logical, intent(in) :: nowell

!Input/output variables
 real(dp), intent(inout) :: ep
 integer, intent(inout) :: nnin !(nnin is actually in/out)
 integer, intent(inout) :: mch

!Output variables
 real(dp), intent(out) :: uu(mmax,2),up(mmax,2),vwell(mmax)
 logical, intent(out) :: confined
 real(dp), intent(out) :: eptgt

!Local variables
 real(dp) :: cwell,et,xx,rwell,rwmax,rwmin,rwscale,ron,umax
 real(dp) :: fw,fwmin,fwmax,tolw

! the well width solves E(rwell) = ep to round-off (doc/smoothness_plan.md,
! Q-9): converged at a few ulps of ep, or when the bracket of rwell has
! closed to a few ulps (E then within eps of ep, else an error)
 real(dp),parameter :: eps=1.0d-10
 integer :: ii,iter,ierr,iumax,l1,itrwell,nn,nnloop,iside
 logical :: convg,bound

 l1=ll+1
 ron=rc*onset
 confined=.true.
 eptgt=ep

!check for bound state with specified ll and nnin

 uu(:,:)=0.0d0 ; up(:,:)=0.0d0
 et=-0.1d0
 call boundstate(irel,nnin,ll,kap,ierr,et,rr,zz,vfull,uu,up,mmax,mch)

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

 if(epgiven) then
!a bound state at or above ep cannot be raised to ep by the well: an input
!error; one below ep and localized inside 1.5 rc only needs raising, by a
!wall as close as that takes
   if(ierr==0 .and. et>=ep .and. .not. nowell) then
     write(6,'(/a,i2,a,i2,a,f12.8/a,f12.8)') &
&          'wellstate: ERROR a bound state exists for n=',nnin,' l=',ll, &
&          ' at energy',et,'wellstate: ERROR above the requested ep=',ep
     error stop 1
   end if
 else
   bound=.false.
   if(ierr==0) then
!if bound state is localized compared to rc, use it and its energy
     if(rr(iumax)<0.75d0*rc) then
       ep=et
       write(6,'(/a,i2,a,i2/a,f12.8/a)') &
&            'WARNING wellstate: localized bound state found for n=', &
&             nnin,' l=', ll,'WARNING this state with energy',ep, &
&            'WARNING will be used for the first projector'

!make sure vwell exists for ikap=2 calculation in this case
       vwell(:)=vfull(:)
       confined=.false.
       eptgt=ep

       return

!if moderately delocalized use very low-energy scattering state
     else if(rr(iumax)<1.5d0*rc) then

       ep=0.05d0
       write(6,'(/2a,i2,a,i2/a,f12.8/a)') &
&            'WARNING wellstate: moderately localized bound state found', &
&            ' for n=', &
&            nnin,' l= ',ll,'WARNING scattering state with energy',ep, &
&            'WARNING will be used for the first projector'
     end if
   end if

!if specified ep is retained and is negative
   if(ep<0.0d0) then
     ep=0.25d0
     write(6,'(/a,i2,a,i2/a,f12.8/a)') &
&           'WARNING wellstate: negative energy specified for n=', &
&           nnin,' l= ',ll, &
&           'WARNING scattering state with energy',ep, &
&           'WARNING will be used for the first projector'
   end if
 end if

!the energy the well's eigenvalue is found for
 eptgt=ep

 if(nowell) then
   vwell(:)=vfull(:)
   uu(:,:)=0.0d0
   up(:,:)=0.0d0
   write(6,'(/a,i3,a,i3/a,f8.4/a)') &
&           '   Wellstate for l =',ll,'  n =',nnin, &
&           '     eigenvalue = ',ep, &
&           '     free tail, no well'
   return
 end if

!scale factor to decrease well width until trial energy exceeds target
 rwscale=0.8d0

!asymptotic potential value to give target bound state a "binding" energy
!of 0.5 Ha
 cwell=ep+0.5d0

!loop which will increment number of nodes before well wall gets too steep
 do nnloop=0,10
   nn=nnin+nnloop

   rwell=8.0d0*rc
   rwmax=0.0d0
   rwmin=0.0d0
   fwmin=0.0d0
   fwmax=0.0d0
   iside=0
   tolw=4.0d0*epsilon(1.0d0)*max(abs(ep),0.2d0)

   convg=.false.

   do itrwell=1,100


!create well potential
     do ii=1,mmax
       vwell(ii)=vfull(ii)
     end do

!   the well starts past rc, where the matching data of the projectors are
!   taken, and rises smoothly
     do ii=1,mmax
       if(rr(ii)<=ron) cycle
       xx=(rr(ii)-ron)/rwell
       vwell(ii)=vwell(ii)+cwell*xx**3/(1.0d0+xx**3)*exp(-0.01d0/xx)
     end do

!find bound state in well
     et=ep
     call boundstate(irel,nn,ll,kap,ierr,et,rr,zz,vwell,uu,up,mmax,mch)
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
       if(rwmax-rwmin<=4.0d0*epsilon(1.0d0)*rwmax) then
         if(abs(fw)>eps) then
           write(6,'(a,i3,a,i3,a,1p,e10.2)') &
&                'wellstate: ERROR E(rwell) jumps across ep, n=',nn, &
&                ' l=',ll,' E-ep=',fw
           error stop 1
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
! state is bound without the well: it then only needs raising
         if(.not. bound .and. rwell<max(1.0d0,0.5d0*rc)) exit
       end if
       cycle
     end if

   end do !itrwell

   if(convg) then
     if(ierr/=0) then
       write(6,'(a,i3,a,i3)') &
&            'wellstate: ERROR lschfb convergence error, n=',nn,' l=',ll
       error stop 1
     end if
     exit
   end if
 end do !nnloop

 if(.not. convg) then
   if(irel==2) then
     write(6,'(a,a,i3,a,i3,a,i3,a,f8.4)') &
&     'ERROR wellstate_r: well potential iteration failed', &
&     ' to converge, n=',nn,' l=',ll,' kap=', kap,' ep=',ep
   else
     write(6,'(a,a,i3,a,i3,a,f8.4)') &
&     'ERROR wellstate: well potential iteration failed', &
&     ' to converge, n=',nn,' l=',ll,' ep=',ep
   end if
   error stop 1
 end if

 if(cwell>0.0d0) then
   if(irel==2) then
     write(6,'(/a,i3,a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',nn,' kap=',kap, &
&             '     eigenvalue = ',et, &
&             '     asymptotic potential = ',cwell, &
&             '     half-point radius =',ron+rwell
   else
     write(6,'(/a,i3,a,i3/a,f8.4,/a,f8.4/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',nn, &
&             '     eigenvalue = ',et, &
&             '     asymptotic potential = ',cwell, &
&             '     half-point radius =',ron+rwell
   end if
 end if

!carry over node count

 nnin=nn

 return
 end subroutine wellstate
