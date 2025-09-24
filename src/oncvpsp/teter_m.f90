module teter_m
   use precision_m, only: dp
   use constants_m, only:  pi, twopi
   implicit none
   private
   public :: teter, teter_deriv_1, teter_deriv_2
contains

!> NAME
!> gg1cc
!>
!> FUNCTION
!> yy=$(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4xx^2)(1-xx^2)})^2$
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx= abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=teter(xx)
function teter(xx) result(yy)

   !Arguments ------------------------------------
   real(dp), intent(in) :: xx

   !Return value -------------------------------------------
   real(dp) :: yy

   !Local variables -------------------------------------------
   ! The cs are coefficients for Taylor expansion of the analytic form near xx=0, 1/2, and 1.
   real(dp), parameter :: c21 = 4.0_dp / 9.0_dp
   real(dp), parameter :: c22 = -40.0_dp / 27.0_dp
   real(dp), parameter :: c23 = 20.0_dp / 3.0_dp - 16.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c24 = -4160.0_dp / 243.0_dp + 160.0_dp * pi**2 / 81.0_dp
   real(dp), parameter :: c31 = 1.0_dp / 36.0_dp
   real(dp), parameter :: c32 = -25.0_dp / 108.0_dp
   real(dp), parameter :: c33 = 485.0_dp / 432.0_dp - pi**2 / 27.0_dp
   real(dp), parameter :: c34 = -4055.0_dp / 972.0_dp + 25.0_dp * pi**2 / 81.0_dp
   real(dp) :: sox
   real(dp) :: zz

   ! Take care of difficult limits near x=0, 1/2, and 1
   if (xx > 3.0_dp) then
      ! Cut off beyond 3/gcut=xcccrc
      yy = 0.0_dp
   else if (abs(xx) <= 1.0e-3_dp) then
      zz = 2.0_dp * pi * xx
      sox = 1.0_dp - zz**2 / 6.0_dp + zz**4 / 120.0_dp + zz**6 / 5040.0_dp
      yy = (sox / ((1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2)))**2
   else if (abs(xx - 0.5_dp) <= 1.0e-4_dp) then
      !  (this limit and next are more troublesome for numerical cancellation)
      yy = c21 + (xx - 0.5_dp) * (c22 + (xx - 0.5_dp) * (c23 + (xx - 0.5_dp) * c24))
   else if (abs(xx - 1.0_dp) <= 1.0e-4_dp) then
      yy = c31 + (xx - 1.0_dp) * (c32 + (xx - 1.0_dp) * (c33 + (xx - 1.0_dp) * c34))
   else
      ! The following is the square of the Fourier transform of a
      ! function built out of two spherical bessel functions in G
      ! space and cut off absolutely beyond gcut
      yy = (sin(2.0_dp * pi * xx) / ((2.0_dp * pi * xx) * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2)))**2
   end if
   return
end function teter

!> NAME
!> gp1cc
!>
!> FUNCTION
!> Derivative of teter(xx) wrt xx.
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx=abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=derivative of teter(xx) wrt xx.
!>
!> NOTES
!> $ phi(x) = \frac{\sin(2\pi x)}{(2\pi x)(1-4x^2)(1-x^2)}$
!> $ teter(x)= phi(x)^2$
!> $ gp(x)= 2 * phi(x) * phi''(x)$
!> $ phi''(x)=\frac{\cos(2\pi x)-(1-15x^2+20x^4) phi(x)}{x(1-4x^2)(1-x^2)}$
function teter_deriv_1(xx) result(yy)

   ! Arguments ------------------------------------
   real(dp), intent(in) :: xx

   ! Return value -------------------------------------------
   real(dp) :: yy

   ! Local variables -------------------------------------------
   real(dp), parameter :: c11 = 20.0_dp - 8.0_dp * pi**2 / 3.0_dp
   real(dp), parameter :: c12 = 268.0_dp - 160.0_dp / 3.0_dp * pi**2 + 128.0_dp / 45.0_dp * pi**4
   real(dp), parameter :: c21 = -40.0_dp / 27.0_dp
   real(dp), parameter :: c22 = 40.0_dp / 3.0_dp - 32.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c23 = -4160.0_dp / 81.0_dp + 160.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c24 = 157712.0_dp / 729.0_dp - 320.0_dp * pi**2 / 9.0_dp + 512.0_dp * pi**4 / 405.0_dp
   real(dp), parameter :: c25 = -452200.0_dp / 729.0_dp + 83200.0_dp * pi**2 / 729.0_dp - 1280.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: c31 = -25.0_dp / 108.0_dp
   real(dp), parameter :: c32 = 485.0_dp / 216.0_dp - 2.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c33 = -4055.0_dp / 324.0_dp + 25.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c34 = 616697.0_dp / 11664.0_dp - 485.0_dp * pi**2 / 81.0_dp + 32.0_dp * pi**4 / 405.0_dp
   real(dp), parameter :: c35 = -2933875.0_dp / 15552.0_dp + 20275.0_dp * pi**2 / 729.0_dp - 200.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: invtwopi = 1.0_dp / twopi
   real(dp) :: denom, phi, phip

   if (xx > 3.0_dp) then
      !  Cut off beyond 3/gcut=3*xcccrc
      yy = 0.0_dp
   else if (xx > 1.0010_dp) then
      ! The part that follows will be repeated later, but written in this way,
      ! only one "if" condition is tested in most of the cases (1.001 < x < 3.0)
      denom = 1.0_dp / (xx * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2))
      phi = denom * sin(twopi * xx) * invtwopi
      phip = denom * (cos(twopi * xx) - (1.0_dp - xx**2 * (15.0_dp - xx**2 * 20.0_dp)) * phi)
      yy = 2.0_dp * phi * phip
      ! Handle limits where denominator vanishes
   else if (abs(xx) < 1.0e-3_dp) then
      yy = xx * (c11 + xx**2 * c12)
   else if (abs(xx - 0.5_dp) <= 1.0e-3_dp) then
      yy = c21 + (xx - 0.5_dp) * (c22 + (xx - 0.5_dp) * (c23 + (xx - 0.5_dp) * (c24 + (xx - 0.5_dp) * c25)))
   else if (abs(xx - 1.0_dp) <= 1.0e-3_dp) then
      yy = c31 + (xx - 1.0_dp) * (c32 + (xx - 1.0_dp) * (c33 + (xx - 1.0_dp) * (c34 + (xx - 1.0_dp) * c35)))
   else
      ! Here is the repeated part ...
      denom = 1.0_dp / (xx * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2))
      phi = denom * sin(twopi * xx) * invtwopi
      phip = denom * (cos(twopi * xx) - (1.0_dp - xx**2 * (15.0_dp - xx**2 * 20.0_dp)) * phi)
      yy = 2.0_dp * phi * phip
   end if
   return
end function teter_deriv_1

!> NAME
!> gpp1cc
!>
!> FUNCTION
!> Second derivative of teter wrt xx.
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM, DRH)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx= abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=second derivative of teter wrt xx.
function teter_deriv_2(xx) result(yy)

   ! Arguments ------------------------------------
   real(dp), intent(in) :: xx

   ! Return value -------------------------------------------
   real(dp) :: yy

   ! Local variables -------------------------------------------
   real(dp), parameter :: c2 = 40.0_dp / 3.0_dp - 32.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c3 = -8320.0_dp / 81.0_dp + 320.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c4 = 157712.0_dp / 243.0_dp - 320.0_dp * pi**2 / 3.0_dp + 512.0_dp * pi**4 / 135.0_dp
   real(dp), parameter :: c5 = -18088.d2 / 729.0_dp + 3328.d2 * pi**2 / 729.0_dp - 5120.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: c6 = 485.0_dp / 216.0_dp - 2.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c7 = -4055.0_dp / 162.0_dp + 50.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c8 = 616697.0_dp / 3888.0_dp - 485.0_dp * pi**2 / 27.0_dp + 32.0_dp * pi**4 / 135.0_dp
   real(dp), parameter :: c9 = -2933875.0_dp / 3888.0_dp + 81100.0_dp * pi**2 / 729.0_dp - 800.0_dp * pi**4 / 243.0_dp
   !Series expansion coefficients around xx=0 from Mathematica
   real(dp), parameter :: cpp0 = 20.0_dp - 8.0_dp * pi**2 / 3.0_dp
   real(dp), parameter :: cpp2 = 804.0_dp - 160.0_dp * pi**2 + (128.0_dp * pi**4) / 15.0_dp
   real(dp), parameter :: cpp4 = 11400.0_dp - 2680.0_dp * pi**2 + (640.0_dp * pi**4) / 3.0_dp &
      - (128.0_dp * pi**6) / 21.0_dp
   real(dp), parameter :: cpp6 = 8.0_dp * (27967275.0_dp - 7182000.0_dp * pi**2 &
                                           + 675360.0_dp * pi**4 &
                                           - 28800.0_dp * pi**6 + 512.0_dp * pi**8) / 2025.0_dp

   real(dp) :: t1, t10, t100, t11, t12, t120, t121, t122, t127, t138, t14, t140, t15, t152
   real(dp) :: t157, t16, t160, t17, t174, t175, t18, t19, t2, t20, t21, t23, t24, t3, t31, t33
   real(dp) :: t34, t4, t41, t42, t44, t45, t46, t5, t54, t55, t56, t57, t6, t62, t64, t65, t7
   real(dp) :: t72, t78, t79, t8, t85, t9, t93

   if (xx > 3.0_dp) then
      !  Cut off beyond 3/gcut=3*xcccrc
      yy = 0.0_dp
      !  Take care of difficult limits near xx=0, 1/2, and 1
   else if (abs(xx) <= 1.0e-3_dp) then
      yy = cpp0 + cpp2 * xx**2 + cpp4 * xx**4 + cpp6 * xx**6
   else if (abs(xx - 0.5_dp) <= 1.0e-4_dp) then
      !  (this limit and next are more troublesome for numerical cancellation)
      yy = c2 + (xx - 0.5_dp) * (c3 + (xx - 0.5_dp) * (c4 + (xx - 0.5_dp) * c5))
   else if (abs(xx - 1.0_dp) <= 1.0e-4_dp) then
      yy = c6 + (xx - 1.0_dp) * (c7 + (xx - 1.0_dp) * (c8 + (xx - 1.0_dp) * c9))
   else
      !  Should fix up this Maple fortran later
      t1 = xx**2
      t2 = 1 / t1
      t3 = 1 / pi
      t4 = 2 * xx
      t5 = t4 - 1
      t6 = t5**2
      t7 = 1 / t6
      t8 = t4 + 1
      t9 = t8**2
      t10 = 1 / t9
      t11 = xx - 1
      t12 = t11**2
      t14 = 1 / t12 / t11
      t15 = xx + 1
      t16 = t15**2
      t17 = 1 / t16
      t18 = pi * xx
      t19 = sin(t18)
      t20 = cos(t18)
      t21 = t20**2
      t23 = t19 * t21 * t20
      t24 = t17 * t23
      t31 = t19**2
      t33 = t31 * t19 * t20
      t34 = t17 * t33
      t41 = pi**2
      t42 = 1 / t41
      t44 = 1 / t16 / t15
      t45 = t31 * t21
      t46 = t44 * t45
      t54 = 1 / t1 / xx
      t55 = 1 / t12
      t56 = t55 * t46
      t57 = t10 * t56
      t62 = t9**2
      t64 = t17 * t45
      t65 = t55 * t64
      t72 = 1 / t9 / t8
      t78 = t14 * t64
      t79 = t10 * t78
      t85 = t12**2
      t93 = t21**2
      t100 = t31**2
      t120 = 1 / t6 / t5
      t121 = t55 * t34
      t122 = t10 * t121
      t127 = t16**2
      t138 = t6**2
      t140 = t10 * t65
      t152 = t72 * t65
      t157 = t7 * t140
      t160 = t1**2
      t174 = t55 * t24
      t175 = t10 * t174
      yy = 8 * t2 * t3 * t7 * t10 * t14 * t34 + 8 * t2 * t42 * t7 * t10 * t14 * t46&
      &   - 8 * t2 * t3 * t7 * t10 * t14 * t24 + 8 * t2 * t3 * t7 * t10 * t55 * t44 * t33 +&
      &   6 * t2 * t42 * t7 * t10 * t55 / t127 * t45 + 24 * t2 * t42 / t138 * t140 +&
      &   16 * t54 * t42 * t120 * t140 + 16 * t2 * t3 * t120 * t122 + 16 * t2&
      &   * t42 * t7 * t72 * t78 - 8 * t2 * t3 * t7 * t10 * t55 * t44 * t23 - 8 * t54 * t3 * t7 * t175&
      &   + 2 * t2 * t7 * t10 * t55 * t17 * t100 + 2 * t2 * t7 * t10 * t55 * t17 * t93 +&
      &   8 * t54 * t42 * t7 * t79 + 16 * t2 * t42 * t7 * t72 * t56 + 6 * t2 * t42 * t7 * t10 / t85&
      &   * t64 + 24 * t2 * t42 * t7 / t62 * t65 + 8 * t54 * t42 * t7 * t57 -&
      &   16 * t2 * t3 * t7 * t72 * t174 + 8 * t54 * t3 * t7 * t122 - 16 * t2 * t3 * t120 * t175&
      &   + 16 * t2 * t42 * t120 * t79 + 16 * t2 * t42 * t120 * t57 + 16 * t54 * t42 * t7 * t152 +&
      &   32 * t2 * t42 * t120 * t152 + 16 * t2 * t3 * t7 * t72 * t121 - 12 * t2 * t157 +&
      &   6 / t160 * t42 * t157
   end if
   return
end function teter_deriv_2

end module teter_m
