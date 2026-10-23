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
 subroutine lschvkbbe(nn,ll,nvkb,ierr,ee,uld,emin,emax, &
&   rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mch)

!! integrates radial schroedinger equation for pseudopotential with
!! Vanderbilt-Kleinman-Bylander non-local projectors finding energy at
!! which desired log derivative uld is matched at point mch, the root
!! with nn-ll-1 eigenvalues below it on (0, rr(mch)] (outwf_index; a node
!! count would also count nodes that enter at r = 0 with the energy)
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param nvkb  = number of VKB projectors to be used
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param uld  log derivative up/uu to be matched at mch
!! @param emin  externally generaated estimate of lower bound for ee
!! @param emax  externally generaated estimate of upper bound for ee
!! @param rr  log radial mesh
!! @param vloc  local part of psp on the mesh
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param evkb  coefficients of VKB projectors
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 implicit none

!Input variables
 real(dp), intent(in) :: uld
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb)
 integer, intent(in) :: nn,ll,nvkb,mmax
 integer, intent(in) :: mch
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Input/output variables
 real(dp), intent(inout) :: emin,emax
 real(dp), intent(inout) :: ee  !in/out, really - needs starting guess

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr


!Local variables
 real(dp) :: cn
 real(dp) :: de,deprev
 real(dp) :: eps,ro,sc
 real(dp) :: sls,sn,uout,upin,upout
 real(dp) :: amesh,al,als
 real(dp) :: xkap
 integer :: ii,it,nin,nint,neig
 type(outwf_t) :: ow

 real(dp), allocatable :: upp(:),cf(:)
 allocate(upp(mmax),cf(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=1.0d-12
 deprev=huge(1.0d0)
 ierr = 60

 sls=ll*(ll+1)


 if(ee>emax) ee=1.25d0*emax
 if(ee<emin) ee=0.75d0*emin
 if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2

! the grids of the outward integration (outwf_m), once for this potential
 call outwf_setup(ll,rr,mmax,vloc,vlocpw,nvkb,vkbpw,ow)

! return point for bound state convergence
 do nint=1,60

! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ee)*rr(ii)**2
   end do

   nin=mch

! outward integration
   call outwf_solve(ow,ee,evkb,mch,mch,uu,up,neig)

   uout=uu(mch)
   upout=up(mch)

   if(neig-nn+ll+1==0) then

     upin=uld*uout

! perform normalization sum

     ro=rr(1)/sqrt(amesh)
     sn=ro**(2*ll+3)/real(2*ll+3,dp)

     do ii=1,nin-3
       sn=sn+al*rr(ii)*uu(ii)**2
     end do

     sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0

! normalize u

     cn=1.0d0/sqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

     do ii=1,nin
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
       up(ii)=0.0d0
     end do

! perturbation theory for energy shift

     de=0.5d0*uout*(upout-upin)/(al*rr(mch))


! convergence test and possible exit

     if(abs(de)<max(abs(ee),0.2d0)*eps) then
       if(abs(de)<=8.0d0*epsilon(1.0d0)*max(abs(ee),0.2d0) .or. &
&         abs(de)>=0.5d0*abs(deprev)) then
         ierr = 0
         exit
       end if
     end if
     deprev=de

     if(de>0.0d0) then
       emin=ee
     else
       emax=ee
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)

   else if(neig-nn+ll+1<0) then
! too few eigenvalues below ee
     emin=ee
     ee=0.5d0*(emin+emax)

   else
! too many eigenvalues below ee
     emax=ee
     ee=0.5d0*(emin+emax)
   end if
 end do

 call outwf_free(ow)
 deallocate(upp,cf)
 return

 end subroutine lschvkbbe
