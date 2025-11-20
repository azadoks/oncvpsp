module teter_m
   implicit none
   private
   public :: gg1cc
   public :: gp1cc
   public :: gpp1cc
contains
!{\src2tex{textfont=tt}}
!!****f* ABINIT/gg1cc
!! NAME
!! gg1cc
!!
!! FUNCTION
!! gg1cc_xx=$(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4xx^2)(1-xx^2)})^2$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx= abscisse to which gg1cc_xx is calculated
!!
!! OUTPUT
!!  gg1cc_xx= gg1cc_x(xx)
!!
!! PARENTS
!!      psp1cc
!!
!! CHILDREN
!!
!! SOURCE
subroutine gg1cc(gg1cc_xx, xx)

   !This section has been created automatically by the script Abilint (TD).
   !Do not modify the following lines by hand.

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

   !Arguments ------------------------------------
   !scalars
   real(dp), intent(in) :: xx
   real(dp), intent(out) :: gg1cc_xx

   !Local variables -------------------------------------------
   !The c s are coefficients for Taylor expansion of the analytic
   !form near xx=0, 1/2, and 1.
   !scalars
   real(dp) :: c21 = 4.d0/9.d0, c22 = -40.d0/27.d0, c23 = 20.d0/3.d0 - 16.d0*pi**2/27.d0
   real(dp) :: c24 = -4160.d0/243.d0 + 160.d0*pi**2/81.d0, c31 = 1.d0/36.d0
   real(dp) :: c32 = -25.d0/108.d0, c33 = 485.d0/432.d0 - pi**2/27.d0
   real(dp) :: c34 = -4055.d0/972.d0 + 25.d0*pi**2/81.d0
   real(dp) :: sox, yy

   ! *************************************************************************

   !Cut off beyond 3/gcut=xcccrc
   if (xx > 3.0d0) then
      gg1cc_xx = 0.0d0
      !  Take care of difficult limits near x=0, 1/2, and 1
      !else if (abs(xx)<=1.d-09) then
      !  gg1cc_xx=1.d0
   else if (abs(xx) <= 1.d-03) then
      yy = 2.0d0*pi*xx
      sox = 1.0d0 - yy**2/6.0d0 + yy**4/120.0d0 + yy**6/5040.0d0
      gg1cc_xx = (sox/((1.d0 - 4.0d0*xx**2)*(1.d0 - xx**2)))**2
   else if (abs(xx - 0.5d0) <= 1.d-04) then
      !  (this limit and next are more troublesome for numerical cancellation)
      gg1cc_xx = c21 + (xx - 0.5d0)*(c22 + (xx - 0.5d0)*(c23 + (xx - 0.5d0)*c24))
   else if (abs(xx - 1.d0) <= 1.d-04) then
      gg1cc_xx = c31 + (xx - 1.0d0)*(c32 + (xx - 1.0d0)*(c33 + (xx - 1.0d0)*c34))
   else
      !  The following is the square of the Fourier transform of a
      !  function built out of two spherical bessel functions in G
      !  space and cut off absolutely beyond gcut
      gg1cc_xx = (sin(2.0d0*pi*xx)/((2.0d0*pi*xx)* &
      &   (1.d0 - 4.0d0*xx**2)*(1.d0 - xx**2)))**2
   end if

end subroutine gg1cc
!!***
!{\src2tex{textfont=tt}}
!!****f* ABINIT/gp1cc
!! NAME
!! gp1cc
!!
!! FUNCTION
!! Derivative of gg(xx) wrt xx.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx=abscisse to which gp1cc_xx is calculated
!!
!! OUTPUT
!!  gp1cc_xx=derivative of gg(xx) wrt xx.
!!
!! NOTES
!! $ phi(x) = \frac{\sin(2\pi x)}{(2\pi x)(1-4x^2)(1-x^2)}$
!! $ gg(x)= phi(x)^2$
!! $ gp(x)= 2 * phi(x) * phi''(x)$
!! $ phi''(x)=\frac{\cos(2\pi x)-(1-15x^2+20x^4) phi(x)}{x(1-4x^2)(1-x^2)}$
!!
!!
!! PARENTS
!!      psp1cc
!!
!! CHILDREN
!!
!! SOURCE
subroutine gp1cc(gp1cc_xx, xx)

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: two_pi = 2.0d0*pi

   !Arguments ------------------------------------
   !scalars
   real(dp), intent(in) :: xx
   real(dp), intent(out) :: gp1cc_xx

   !Local variables -------------------------------------------
   !scalars
   real(dp), parameter :: c11 = 20.d0 - 8.d0*pi**2/3.d0
   real(dp), parameter :: c12 = 268.d0 - 160.d0/3.d0*pi**2 + 128.d0/45.d0*pi**4
   real(dp), parameter :: c21 = -40.d0/27.d0, c22 = 40.d0/3.d0 - 32.d0*pi**2/27.d0
   real(dp), parameter :: c23 = -4160.d0/81.d0 + 160.d0*pi**2/27.d0
   real(dp), parameter :: c24 = 157712.d0/729.d0 - 320.d0*pi**2/9.d0 + 512.d0*pi**4/405.d0
   real(dp), parameter :: c25 = -452200.d0/729.d0 + 83200.d0*pi**2/729.d0 - 1280.d0*pi**4/243.d0
   real(dp), parameter :: c31 = -25.d0/108.d0, c32 = 485.d0/216.d0 - 2.d0*pi**2/27.d0
   real(dp), parameter :: c33 = -4055.d0/324.d0 + 25.d0*pi**2/27.d0
   real(dp), parameter :: c34 = 616697.d0/11664.d0 - 485.d0*pi**2/81.d0 + 32.d0*pi**4/405.d0
   real(dp), parameter :: c35 = -2933875.d0/15552.d0 + 20275.d0*pi**2/729.d0 - 200.d0*pi**4/243.d0
   real(dp), parameter :: two_pim1 = 1.0d0/two_pi
   real(dp) :: denom, phi, phip

   ! *************************************************************************

   !Cut off beyond r=3*xcccrc is already done at the calling level
   if (xx > 1.001d0) then
      !  The part that follows will be repeated later, but written in this way,
      !  only one "if" condition is tested in most of the cases (1.001 < x < 3.0)
      denom = 1.d0/(xx*(1.d0 - 4.d0*xx**2)*(1.d0 - xx**2))
      phi = denom*sin(two_pi*xx)*two_pim1
      phip = denom*(cos(two_pi*xx) - (1.d0 - xx**2*(15.d0 - xx**2*20))*phi)
      gp1cc_xx = 2.d0*phi*phip
      !  Handle limits where denominator vanishes
   else if (abs(xx) < 1.d-03) then
      gp1cc_xx = xx*(c11 + xx**2*c12)
   else if (abs(xx - 0.5d0) <= 1.d-03) then
      gp1cc_xx = c21 + (xx - 0.5d0)*(c22 + (xx - 0.5d0)*(c23 + (xx - 0.5d0)*(c24 + (xx - 0.5d0)*c25)))
   else if (abs(xx - 1.d0) <= 1.d-03) then
      gp1cc_xx = c31 + (xx - 1.0d0)*(c32 + (xx - 1.0d0)*(c33 + (xx - 1.0d0)*(c34 + (xx - 1.0d0)*c35)))
   else
      !  Here is the repeated part ...
      denom = 1.d0/(xx*(1.d0 - 4.d0*xx**2)*(1.d0 - xx**2))
      phi = denom*sin(two_pi*xx)*two_pim1
      phip = denom*(cos(two_pi*xx) - (1.d0 - xx**2*(15.d0 - xx**2*20))*phi)
      gp1cc_xx = 2.d0*phi*phip
   end if

end subroutine gp1cc
!!***
!{\src2tex{textfont=tt}}
!!****f* ABINIT/gpp1cc
!! NAME
!! gpp1cc
!!
!! FUNCTION
!! Second derivative of gg wrt xx.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx= abscisse to which gpp1cc_xx is calculated
!!
!! OUTPUT
!!  gpp1cc_xx=second derivative of gg wrt xx.
!!
!!
!! PARENTS
!!      psp1cc
!!
!! CHILDREN
!!
!! SOURCE
subroutine gpp1cc(gpp1cc_xx, xx)

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: two_pi = 2.0d0*pi

   !Arguments ------------------------------------
   !scalars
   real(dp), intent(in) :: xx
   real(dp), intent(out) :: gpp1cc_xx

   !Local variables -------------------------------------------
   !scalars
   real(dp), parameter :: c2 = 40.d0/3.d0 - 32.d0*pi**2/27.d0
   real(dp), parameter :: c3 = -8320.d0/81.d0 + 320.d0*pi**2/27.d0
   real(dp), parameter :: c4 = 157712.d0/243.d0 - 320.d0*pi**2/3.d0 + 512.d0*pi**4/135.d0
   real(dp), parameter :: c5 = -18088.d2/729.d0 + 3328.d2*pi**2/729.d0 - 5120.d0*pi**4/243.d0
   real(dp), parameter :: c6 = 485.d0/216.d0 - 2.d0*pi**2/27.d0
   real(dp), parameter :: c7 = -4055.d0/162.d0 + 50.d0*pi**2/27.d0
   real(dp), parameter :: c8 = 616697.d0/3888.d0 - 485.d0*pi**2/27.d0 + 32.d0*pi**4/135.d0
   real(dp), parameter :: c9 = -2933875.d0/3888.d0 + 81100.d0*pi**2/729.d0 - 800.d0*pi**4/243.d0

   !Series expansion coefficients around xx=0 from Mathematica
   real(dp), parameter :: cpp0 = 20.d0 - 8.d0*pi**2/3.d00

   real(dp), parameter :: cpp2 = 804.d0 - 160.d0*pi**2 + (128.d0*pi**4)/15.d0

   real(dp), parameter :: cpp4 = 11400.d0 - 2680.d0*pi**2 &
   &                          + (640.d0*pi**4)/3.d0 - (128.d0*pi**6)/21.d0

   real(dp), parameter :: cpp6 = 8.d0*(27967275.d0 - 7182000.d0*pi**2 &
   &                          + 675360.d0*pi**4 &
   &                          - 28800.d0*pi**6 + 512.d0*pi**8)/2025.d0

   real(dp) :: t1, t10, t100, t11, t12, t120, t121, t122, t127, t138, t14, t140, t15, t152
   real(dp) :: t157, t16, t160, t17, t174, t175, t18, t19, t2, t20, t21, t23, t24, t3, t31, t33
   real(dp) :: t34, t4, t41, t42, t44, t45, t46, t5, t54, t55, t56, t57, t6, t62, t64, t65, t7
   real(dp) :: t72, t78, t79, t8, t85, t9, t93

   ! *************************************************************************

   if (xx > 3.0d0) then
      !  Cut off beyond 3/gcut=3*xcccrc
      gpp1cc_xx = 0.0d0
      !  Take care of difficult limits near xx=0, 1/2, and 1
   else if (abs(xx) <= 1.d-03) then
      gpp1cc_xx = cpp0 + cpp2*xx**2 + cpp4*xx**4 + cpp6*xx**6
   else if (abs(xx - 0.5d0) <= 1.d-04) then
      !  (this limit and next are more troublesome for numerical cancellation)
      gpp1cc_xx = c2 + (xx - 0.5d0)*(c3 + (xx - 0.5d0)*(c4 + (xx - 0.5d0)*c5))
   else if (abs(xx - 1.d0) <= 1.d-04) then
      gpp1cc_xx = c6 + (xx - 1.0d0)*(c7 + (xx - 1.0d0)*(c8 + (xx - 1.0d0)*c9))
   else

      !  Should fix up this Maple fortran later
      t1 = xx**2
      t2 = 1/t1
      t3 = 1/Pi
      t4 = 2*xx
      t5 = t4 - 1
      t6 = t5**2
      t7 = 1/t6
      t8 = t4 + 1
      t9 = t8**2
      t10 = 1/t9
      t11 = xx - 1
      t12 = t11**2
      t14 = 1/t12/t11
      t15 = xx + 1
      t16 = t15**2
      t17 = 1/t16
      t18 = Pi*xx
      t19 = sin(t18)
      t20 = cos(t18)
      t21 = t20**2
      t23 = t19*t21*t20
      t24 = t17*t23
      t31 = t19**2
      t33 = t31*t19*t20
      t34 = t17*t33
      t41 = Pi**2
      t42 = 1/t41
      t44 = 1/t16/t15
      t45 = t31*t21
      t46 = t44*t45
      t54 = 1/t1/xx
      t55 = 1/t12
      t56 = t55*t46
      t57 = t10*t56
      t62 = t9**2
      t64 = t17*t45
      t65 = t55*t64
      t72 = 1/t9/t8
      t78 = t14*t64
      t79 = t10*t78
      t85 = t12**2
      t93 = t21**2
      t100 = t31**2
      t120 = 1/t6/t5
      t121 = t55*t34
      t122 = t10*t121
      t127 = t16**2
      t138 = t6**2
      t140 = t10*t65
      t152 = t72*t65
      t157 = t7*t140
      t160 = t1**2
      t174 = t55*t24
      t175 = t10*t174
      gpp1cc_xx = 8*t2*t3*t7*t10*t14*t34 + 8*t2*t42*t7*t10*t14*t46&
      &   - 8*t2*t3*t7*t10*t14*t24 + 8*t2*t3*t7*t10*t55*t44*t33 +&
      &   6*t2*t42*t7*t10*t55/t127*t45 + 24*t2*t42/t138*t140 +&
      &   16*t54*t42*t120*t140 + 16*t2*t3*t120*t122 + 16*t2&
      &   *t42*t7*t72*t78 - 8*t2*t3*t7*t10*t55*t44*t23 - 8*t54*t3*t7*t175&
      &   + 2*t2*t7*t10*t55*t17*t100 + 2*t2*t7*t10*t55*t17*t93 +&
      &   8*t54*t42*t7*t79 + 16*t2*t42*t7*t72*t56 + 6*t2*t42*t7*t10/t85&
      &   *t64 + 24*t2*t42*t7/t62*t65 + 8*t54*t42*t7*t57 -&
      &   16*t2*t3*t7*t72*t174 + 8*t54*t3*t7*t122 - 16*t2*t3*t120*t175&
      &   + 16*t2*t42*t120*t79 + 16*t2*t42*t120*t57 + 16*t54*t42*t7*t152 +&
      &   32*t2*t42*t120*t152 + 16*t2*t3*t7*t72*t121 - 12*t2*t157 +&
      &   6/t160*t42*t157
   end if

end subroutine gpp1cc
!!***
end module teter_m
