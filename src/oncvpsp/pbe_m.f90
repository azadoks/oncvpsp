!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! calculates PBE exchange-correlation potential and energy density
!
! The core of this routine was obtained from John Perdew a long time
! ago and does not follow the coding style of the rest of ONCVPSP
!
module pbe_m
    use constants_m, only: dp, pi, fourpi, ifourpi, third
    private
    public :: excggc
contains

subroutine excggc(rho,vxc,exc,r,mmax)
!
! Error (from Perdew''s group) corrected 12/27/00 in expression for FSS
! in routine exchpbe. (See NOTES file.)
!
    implicit none
    integer, parameter :: dp=kind(1.0d0)
!
    real*8 rho(*),vxc(*),exc(*),r(*)
    real*8 amesh,al
    real*8 fk,sk,g,ec
    integer i, mmax
!
!      common/GAS/fk,sk,g,ec,ecrs,eczet
!
    real(dp), allocatable :: d(:),dpn(:),dppn(:),dpr(:)
    real(dp), allocatable :: dppr(:),dlap(:)
    real(dp) :: pi,pi4,pi4i
    real(dp) :: c11,c12,c13,c14,c15,c21,c22,c23,c24,c25
    real(dp) :: thrd
    real(dp) :: conf,conrs
    real(dp) :: s,u,v,ex,vx
    real(dp) :: rs,zet,vcup,vcdn
    real(dp) :: t,uu,vv,ww,h,dvcup,dvcdn
!
    al = 0.01d0 * dlog(r(101) / r(1))
    amesh = dexp(al)
!
    thrd = 1.0d0 / 3.0d0
    pi=4.0d0*datan(1.0d0)
    pi4=4.0d0 * pi

    allocate(d(mmax),dpn(mmax),dppn(mmax),dpr(mmax))
    allocate(dppr(mmax),dlap(mmax))
!
    conf = (3.d0*pi**2)**thrd
    conrs = (3.d0/(4.d0*pi))**thrd
!
    c11 =   2.0d0 / 24.0d0
    c12 = -16.0d0 / 24.0d0
    c13 =   0.0d0 / 24.0d0
    c14 =  16.0d0 / 24.0d0
    c15 =  -2.0d0 / 24.0d0
!
    c21 =   -1.0d0 / 12.0d0
    c22 =   16.0d0 / 12.0d0
    c23 =  -30.0d0 / 12.0d0
    c24 =   16.0d0 / 12.0d0
    c25 =   -1.0d0 / 12.0d0
!
! d is properly scaled charge density
!
    pi4i = 1.0d0 / pi4
    do i = 1, mmax
        d(i) = pi4i * rho(i)
    end do
!
! n derivatives of d
!
    do i = 3, mmax - 2
        dpn(i) =  c11*d(i-2) + c12*d(i-1) + c14*d(i+1) + c15*d(i+2)
        dppn(i) = c21*d(i-2) + c22*d(i-1) + c23*d(i)   + c24*d(i+1)&
        &+ c25*d(i+2)
    end do
!
! r derivatives of d
!
    do i = 3, mmax - 2
        dpr(i) = dpn(i) / (al * r(i))
        dppr(i) = (dppn(i) - al * dpn(i)) / (al * r(i))**2
        dlap(i) = (dppn(i) + al * dpn(i)) / (al * r(i))**2
    end do
!
! set up input for Perdew''s subroutines and call them
!
    do i = 3, mmax - 2
!
!     SUBROUTINE EXCH(D,S,U,V,EX,VX)
!  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  INPUT D : DENSITY
!  INPUT S:  ABS(GRAD D)/(2*KF*D)
!  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!
        if(d(i) .gt. 1.0d-18) then
!
            fk = conf * d(i) ** thrd
!
            s = dabs(dpr(i)) / (2.0d0 * fk * d(i))
!
            u = dabs(dpr(i)) * dppr(i) / (d(i)**2 * (2.0d0 * fk)**3)
!
            v = dlap(i) / (d(i) * (2.0d0 * fk)**2)
!
!          call exchpbe(d(i),s,u,v,ex,vx)
            call EXCHPBE(d(i),s,u,v,1,1,ex,vx)
!
!     SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!
!
            rs = conrs / d(i)**thrd
!
            zet = 0.0d0
            g = 1.0d0
!
!          call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,zlfc)
!
!     SUBROUTINE CORGGA(RS,ZET,T,UU,VV,WW,H,DVCUP,DVCDN)
!  GGA91 CORRELATION
!  INPUT RS: SEITZ RADIUS
!  INPUT ZET: RELATIVE SPIN POLARIZATION
!  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
!  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
!  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
!  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
!  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
! note g = 1 for unpolarized case
!
            sk = dsqrt(4.0d0 * fk / pi)
!
            t = dabs(dpr(i)) / (d(i) * 2.0d0 * sk)
!
            uu = dabs(dpr(i)) * dppr(i) / (d(i)**2 * (2.0d0 * sk)**3)
!
            vv = dlap(i) / (d(i) * (2.0d0 * sk)**2)
!
            ww = 0.0d0
!
!          call corpbe(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn)
            call CORPBE(rs,zet,t,uu,vv,ww,1,1,ec,vcup,vcdn,&
            &h,dvcup,dvcdn)
!
        else
            ex = 0.0d0
            vx = 0.0d0
            ec = 0.0d0
            vcup = 0.0d0
            vcdn = 0.0d0
            h = 0.0d0
            dvcup = 0.0d0
            dvcdn = 0.0d0
        end if
!
        vxc(i) = vx + 0.5d0 * (vcup + vcdn + dvcup + dvcdn)
!
        exc(i) = ex + ec + h
    end do
!
! assume end values are approximately constant
!
    vxc(1) = vxc(3)
    vxc(2) = vxc(3)
    exc(1) = exc(3)
    exc(2) = exc(3)
    vxc(mmax - 1) = vxc(mmax - 2)
    vxc(mmax) = vxc(mmax - 2)
    exc(mmax - 1) = exc(mmax - 2)
    exc(mmax) = exc(mmax - 2)
!
    deallocate(d,dpn,dppn,dpr)
    deallocate(dppr,dlap)
!
    return
end subroutine excggc
!
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke''s modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  input lgga:  (=0=>don''t put in gradient corrections, just LDA)
!  input lpot:  (=0=>don''t get potential and don''t need U and V)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!   	e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!	e_x[PBE]=e_x[unif]*FxPBE(s)
!	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
    parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
    parameter(pi=3.14159265358979323846264338327950d0)
    parameter(ax=-0.738558766382022405884230032680836d0)
    parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
    exunif = AX*rho**THRD
    if(lgga.eq.0)then
        ex=exunif
        vx=ex*thrd4
        return
    endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
    S2 = S*S
    P0=1.d0+ul*S2
    FxPBE = 1d0+uk-uk/P0
    EX = exunif*FxPBE
    if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
    Fs=2.d0*uk*ul/(P0*P0)
    Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24)
    VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
    RETURN
end subroutine EXCHPBE
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,H,DVCUP,DVCDN)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
!       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!       :  UU,VV,WW, only needed for PBE potential
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    IMPLICIT REAL*8 (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
    parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
    parameter(sixthm=thrdm/2.d0)
    parameter(thrd4=4.d0*thrd)
    parameter(GAM=0.5198420997897463295344212145565d0)
    parameter(fzz=8.d0/(9.d0*GAM))
    parameter(gamma=0.03109069086965489503494086371273d0)
    parameter(bet=0.06672455060314922d0,delt=bet/gamma)
    parameter(eta=1.d-12)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
    rtrs=dsqrt(rs)
    CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,&
    &0.49294D0,rtrs,EU,EURS)
    CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,&
    &0.62517D0,rtRS,EP,EPRS)
    CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,&
    &0.49671D0,rtRS,ALFM,ALFRSM)
    ALFC = -ALFM
    Z4 = ZET**4
    F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
    EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
    ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
    FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
    ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU&
    &-(1.D0-Z4)*ALFM/FZZ)
    COMM = EC -RS*ECRS/3.D0-ZET*ECZET
    VCUP = COMM + ECZET
    VCDN = COMM - ECZET
    if(lgga.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
    G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
    G3 = G**3
    PON=-EC/(G3*gamma)
    B = DELT/(DEXP(PON)-1.D0)
    B2 = B*B
    T2 = T*T
    T4 = T2*T2
    RS2 = RS*RS
    RS3 = RS2*RS
    Q4 = 1.D0+B*T2
    Q5 = 1.D0+B*T2+B2*T4
    H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
    if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
    G4 = G3*G
    T6 = T4*T2
    RSTHRD = RS/3.D0
    GZ=(((1.d0+zet)**2+eta)**sixthm-&
    &((1.d0-zet)**2+eta)**sixthm)/3.d0
    FAC = DELT/B+1.D0
    BG = -3.D0*B2*EC*FAC/(BET*G4)
    BEC = B2*FAC/(BET*G3)
    Q8 = Q5*Q5+DELT*Q4*Q5*T2
    Q9 = 1.D0+2.D0*B*T2
    hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
    hRS = -RSTHRD*hB*BEC*ECRS
    FACT0 = 2.D0*DELT-6.D0*B
    FACT1 = Q5*Q9+Q4*Q9*Q9
    hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
    hRST = RSTHRD*T2*hBT*BEC*ECRS
    hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
    hT = 2.d0*BET*G3*Q9/Q8
    hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
    FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
    FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
    hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
    COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
    PREF = HZ-GZ*T2*HT/G
    FACT5 = GZ*(2.D0*HT+T*HTT)/G
    COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
    DVCUP = COMM + PREF
    DVCDN = COMM - PREF
    RETURN
end subroutine CORPBE
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
    ! slimmed down version of GCOR used in PW91 routines, to interpolate
    ! LSD correlation energy, as given by (10) of
    ! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
    ! K. Burke, May 11, 1996.
    
    implicit none
    
    ! Input variables
    real(dp), intent(in) :: a, a1, b1, b2, b3, b4, rtrs
    real(dp), intent(out) :: gg, ggrs
    
    ! Local variables
    real(dp) :: q0, q1, q2, q3

    ! implicit real*8 (a-h,o-z)
    q0 = -2.d0*a*(1.d0+a1*rtrs*rtrs)
    q1 = 2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
    q2 = log(1.d0+1.d0/q1)
    gg = q0*q2
    q3 = a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
    ggrs = -2.d0*a*a1*q2-q0*q3/(q1*(1.d0+q1))

    return
end subroutine gcor2

end module pbe_m
