# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from numpy.testing import assert_allclose
import pytest
from math import *
from chemicals.iapws import *
from fluids.numerics import assert_close, linspace, logspace
from chemicals.iapws import REGION_3A, REGION_3B, REGION_3C, REGION_3D, REGION_3E, REGION_3F, REGION_3G, REGION_3H, REGION_3I, REGION_3J, REGION_3K, REGION_3L, REGION_3M, REGION_3N, REGION_3O, REGION_3P, REGION_3Q, REGION_3R, REGION_3S, REGION_3T, REGION_3U, REGION_3V, REGION_3W, REGION_3X, REGION_3Y, REGION_3Z
from chemicals.vapor_pressure import Psat_IAPWS


### IAPWS Region 1 tests
nis1 = [0.14632971213167, -0.84548187169114, -0.37563603672040E1,
       0.33855169168385E1, -0.95791963387872, 0.15772038513228,
       -0.16616417199501E-1, 0.81214629983568E-3, 0.28319080123804E-3,
       -0.60706301565874E-3, -0.18990068218419E-1, -0.32529748770505E-1,
       -0.21841717175414E-1, -0.52838357969930E-4, -0.47184321073267E-3,
       -0.30001780793026E-3, 0.47661393906987E-4, -0.44141845330846E-5,
       -0.72694996297594E-15, -0.31679644845054E-4, -0.28270797985312E-5,
       -0.85205128120103E-9, -0.22425281908000E-5, -0.65171222895601E-6,
       -0.14341729937924E-12, -0.40516996860117E-6, -0.12734301741641E-8,
       -0.17424871230634E-9, -0.68762131295531E-18, 0.14478307828521E-19,
       0.26335781662795E-22, -0.11947622640071E-22, 0.18228094581404E-23,
       -0.93537087292458E-25]
lis1 = [0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 2., 2., 2., 2., 2., 3., 3., 3., 4., 4., 4., 5., 8., 8., 21., 23., 29., 30., 31., 32.]
lis1 = [int(i) for i in lis1]
Jis1 = [-2., -1., 0., 1., 2., 3., 4., 5., -9., -7., -1., 0., 1., 3., -3., 0., 1., 3., 17., -4., 0., 6., -5., -2., 10., -8., -11., -6., -29., -31., -38., -39., -40., -41.]
Jis1 = [int(i) for i in Jis1]


def iapws97_G_region1_naive(tau, pi):
    return sum([nis1[i]*(7.1-pi)**lis1[i]*(tau-1.222)**Jis1[i] for i in range(34)])

def iapws97_dG_dpi_region1_naive(tau, pi):
    return sum([-nis1[i]*lis1[i]*(7.1-pi)**(lis1[i]-1)*(tau-1.222)**Jis1[i] for i in range(34)])

def iapws97_d2G_d2pi_region1_naive(tau, pi):
    return sum([nis1[i]*lis1[i]*(lis1[i]-1)*(7.1-pi)**(lis1[i]-2)*(tau-1.222)**Jis1[i] for i in range(34)])

def iapws97_dG_dtau_region1_naive(tau, pi):
    return sum([nis1[i]*Jis1[i]*(7.1-pi)**lis1[i]*(tau-1.222)**(Jis1[i]-1) for i in range(34)])

def iapws97_d2G_d2tau_region1_naive(tau, pi):
    return sum([nis1[i]*Jis1[i]*(Jis1[i]-1)*(7.1-pi)**lis1[i]*(tau-1.222)**(Jis1[i]-2) for i in range(34)])

def iapws97_d2G_dpidtau_region1_naive(tau, pi):
    return sum([-nis1[i]*Jis1[i]*lis1[i]*(7.1-pi)**(lis1[i]-1)*(tau-1.222)**(Jis1[i]-1) for i in range(34)])


# Section 2 - ideal gas part
J0is2 = [0., 1., -5., -4., -3., -2., -1., 2., 3.]
J0is2 = [int(i) for i in J0is2]

def iapws97_G0_region2_naive(tau, pi):
    return log(pi) + sum( [n0is2[i]*tau**J0is2[i] for i in range(9)] )

#        self.G0 = log(pi) + sum( [n0is2[i]*tau**J0is2[i] for i in range(9)] )
#        self.dG0_dpi = 1./pi
#        self.ddG0_ddpi = -1/pi**2
#        self.dG0_dtau = sum( [n0is2[i]*J0is2[i]*tau**(J0is2[i]-1) for i in range(9)] )
#        self.ddG0_ddtau = sum( [n0is2[i]*J0is2[i]*(J0is2[i]-1)*tau**(J0is2[i]-2) for i in range(9)] )



n0is2 = [-0.96927686500217E1, 0.10086655968018E2, -0.56087911283020E-2,
        0.71452738081455E-1, -0.40710498223928, 0.14240819171444E1,
        -0.43839511319450E1, -0.28408632460772, 0.21268463753307E-1]
# Section 2 - residual part

lis2 = [1., 1., 1., 1., 1., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 4., 4., 4.,
        5., 6., 6., 6., 7., 7., 7., 8., 8., 9., 10., 10., 10., 16., 16., 18.,
        20., 20., 20., 21., 22., 23., 24., 24., 24.]
lis2 = [int(i) for i in lis2]
Jis2 = [0., 1., 2., 3., 6., 1., 2., 4., 7., 36., 0., 1., 3., 6., 35., 1., 2.,
        3., 7., 3., 16., 35., 0., 11., 25., 8., 36., 13., 4., 10., 14., 29.,
        50., 57., 20., 35., 48., 21., 53., 39., 26., 40., 58.]
Jis2 = [int(i) for i in Jis2]
nis2 = [-0.17731742473213E-2, -0.17834862292358E-1, -0.45996013696365E-1,
        -0.57581259083432E-1, -0.50325278727930E-1, -0.33032641670203E-4,
        -0.18948987516315E-3, -0.39392777243355E-2, -0.43797295650573E-1,
        -0.26674547914087E-4, 0.20481737692309E-7, 0.43870667284435E-6,
        -0.32277677238570E-4, -0.15033924542148E-2, -0.40668253562649E-1,
        -0.78847309559367E-9, 0.12790717852285E-7, 0.48225372718507E-6,
        0.22922076337661E-5, -0.16714766451061E-10, -0.21171472321355E-2,
        -0.23895741934104E2, -0.59059564324270E-17, -0.12621808899101E-5,
        -0.38946842435739E-1, 0.11256211360459E-10, -0.82311340897998E1,
        0.19809712802088E-7, 0.10406965210174E-18, -0.10234747095929E-12,
        -0.10018179379511E-8, -0.80882908646985E-10, 0.10693031879409,
        -0.33662250574171, 0.89185845355421E-24, 0.30629316876232E-12,
        -0.42002467698208E-5, -0.59056029685639E-25, 0.37826947613457E-5,
        -0.12768608934681E-14, 0.73087610595061E-28, 0.55414715350778E-16,
        -0.94369707241210E-6]

def iapws97_Gr_region2_naive(tau, pi):
    return sum([nis2[i]*pi**lis2[i]*(tau-0.5)**Jis2[i] for i in range(43)])

def iapws97_dGr_dpi_region2_naive(tau, pi):
    return sum([nis2[i]*lis2[i]*pi**(lis2[i]-1)*(tau-0.5)**Jis2[i] for i in range(43)])

def iapws97_ddGr_ddpi_region2_naive(tau, pi):
    return sum([nis2[i]*lis2[i]*(lis2[i]-1)*pi**(lis2[i]-2)*(tau-0.5)**Jis2[i] for i in range(43)])

def iapws97_dGr_dtau_region2_naive(tau, pi):
    return sum([nis2[i]*pi**lis2[i]*Jis2[i]*(tau-0.5)**(Jis2[i]-1) for i in range(43)])

def iapws97_ddGr_ddtau_region2_naive(tau, pi):
    return sum([nis2[i]*pi**lis2[i]*Jis2[i]*(Jis2[i]-1)*(tau-0.5)**(Jis2[i]-2) for i in range(43)])

def iapws97_ddGr_dpi_dtau_region2_naive(tau, pi):
    return sum([nis2[i]*lis2[i]*pi**(lis2[i]-1)*Jis2[i]*(tau-0.5)**(Jis2[i]-1) for i in range(43)])



def test_iapws97_dG_dpi_region1():
    assert_close(iapws97_dG_dpi_region1_naive(1386/277.15, 101325/16.53E6),
                 iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6), rtol=1e-14)
    
    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 0.12923271825448354, rtol=1e-14)
    
    # Point that had bad error with horner's method
    assert_close(iapws97_dG_dpi_region1(1386 / 600.15, 10001325 / 16.53E6),
                 0.09345587583404263, rtol=1e-14)
    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 iapws97_dG_dpi_region1_naive(1386/277.15, 101325/16.53E6), rtol=1e-14)

@pytest.mark.slow
def test_iapws97_dG_dpi_region1_fuzz():
    funcs_naive = [iapws97_dG_dpi_region1_naive, iapws97_G_region1_naive, iapws97_d2G_d2pi_region1_naive,
                   iapws97_dG_dtau_region1_naive, iapws97_d2G_d2tau_region1_naive, iapws97_d2G_dpidtau_region1_naive]
    funcs_fast = [iapws97_dG_dpi_region1, iapws97_G_region1, iapws97_d2G_dpi2_region1,
                  iapws97_dG_dtau_region1, iapws97_d2G_d2tau_region1, iapws97_d2G_dpidtau_region1]
    atols = [0, 1e-14, 0, 3e-15, 0.0, 1e-16]
    rtols = [2e-13, 1e-12, 3e-12, 1e-13, 1e-12, 2e-12]

    # for testing faster
#    funcs_naive = [iapws97_d2G_dpidtau_region1_naive]
#    funcs_fast = [iapws97_d2G_dpidtau_region1]
#    atols = [1e-16]
#    rtols = [2e-12]
    
    N = 500
    Ts = linspace(273.15, 623.15, N)
    def test_Ps(T, N):
        Psat = Psat_IAPWS(T)
        return logspace(log10(Psat), log10(100e6), N)
    
    for naive, fast, rtol, atol in zip(funcs_naive, funcs_fast, rtols, atols):
        for T in Ts:
            tau = 1386.0/T
            for P in test_Ps(T, N):
                pi = P/16.53E6
                assert_close(naive(tau, pi),
                             fast(tau, pi), rtol=rtol, atol=atol)


#test_iapws97_dG_dpi_region1_fuzz()

def test_iapws97_dG_dpi_region2():
    assert_close(iapws97_dGr_dpi_region2(.656, 16), -0.006292631931275252, rtol=1e-14)
    
    
    # Point that had bad error with horner's method
    assert_close(iapws97_dGr_dpi_region2(0.788009171330091, 26.87134177929712), -0.018525334158583723, rtol=1e-14)

def test_iapws97_dG_dpi_region5():
    assert_close(iapws97_dGr_dpi_region5(.5, 30.0), 0.0004009761854002751, rtol=1e-14)

def test_iapws_boundary_equations():
    assert_close(iapws97_boundary_2_3(0.623150000E3), 16529164.252621626, rtol=1e-13)
    
    assert_close(iapws97_boundary_3uv(22.3E6), 647.7996121480069, rtol=1e-14)
    
    assert_close(iapws97_boundary_3ef(40E6), 713.959399239744, rtol=1e-14)
    
    assert_close(iapws97_boundary_3cd(25E6), 649.3659208321279, rtol=1e-14)
    
    assert_close(iapws97_boundary_3gh(25E6), 656.6980572261236, rtol=2e-14)
    
    assert_close(iapws97_boundary_3ij(25E6), 660.7865756716819, rtol=1e-14)
    
    assert_close(iapws97_boundary_3jk(25E6), 668.1915358826951, rtol=1e-14)
    
    assert_close(iapws97_boundary_3mn(22.8E6), 649.6054132953997, rtol=1e-14)
    
    assert_close(iapws97_boundary_3qu(22E6), 645.6355027340121, rtol=1e-14)
    
    assert_close(iapws97_boundary_3rx(22E6), 648.2622753670172, rtol=1e-14)
    
    assert_close(iapws97_boundary_3wx(log(22.3), 1 / log(22.3)), 648.204947950734, rtol=1e-14)
    
    assert_close(iapws97_boundary_3ab(log(40), 1 / log(40)), 693.0341408296053, rtol=1e-14)
    
    assert_close(iapws97_boundary_3op(log(22.8), 1 / log(22.8)), 650.010694314133, rtol=1e-14)
    
def test_iapws97_region_3_misc():
    assert iapws97_region_3(630, 50e6) == REGION_3A
    assert iapws97_region_3(709.5013, 50e6) == REGION_3B
    assert iapws97_region_3(709.5012, 50e6) == REGION_3A
    
    assert iapws97_region_3(700.0, 30e6) == REGION_3F
    
    # CoolProp differs but http://twt.mpei.ac.ru/MCS/Worksheets/iapws/IAPWS-IF97-Region3-VPT.xmcd confirms it is C here.
    # A test with IAPWS95 shows CoolProp matches the correct answer better
    # We are right next to a transition point / huge discontinuity here.
    assert iapws97_region_3(623.1500000001, 16529164.269161053) == REGION_3C



def test_iapws97_region_full_table():
    assert iapws97_region_3(630, 50e6) == REGION_3A
    assert_close(1/iapws97_region3_rho(T=630, P=50e6), 0.001470853100110911, rtol=1e-13)
    assert iapws97_region_3(670, 80e6) == REGION_3A
    assert_close(1/iapws97_region3_rho(T=670, P=80e6), 0.0015038313585404727, rtol=1e-13)
    
    assert iapws97_region_3(710.0, 50e6) == REGION_3B
    assert_close(1/iapws97_region3_rho(T=710, P=50e6), 0.0022047285870574838, rtol=1e-13)
    assert iapws97_region_3(750.0, 80e6) == REGION_3B
    assert_close(1/iapws97_region3_rho(T=750, P=80e6), 0.0019736929401211155, rtol=1e-13)
    
    assert iapws97_region_3(630.0, 20e6) == REGION_3C
    assert_close(1/iapws97_region3_rho(T=630, P=20e6), 0.0017616964055295276, rtol=1e-13)
    assert iapws97_region_3(650.0, 30e6) == REGION_3C
    assert_close(1/iapws97_region3_rho(T=650, P=30e6), 0.0018195606165288332, rtol=1e-13)
    
    assert iapws97_region_3(656.0, 26e6) == REGION_3D
    assert_close(1/iapws97_region3_rho(T=656, P=26e6), 0.002245587720029806, rtol=1e-13)
    assert iapws97_region_3(670.0, 30e6) == REGION_3D
    assert_close(1/iapws97_region3_rho(T=670, P=30e6), 0.002506897701629579, rtol=1e-13)
    
    assert iapws97_region_3(661.0, 26e6) == REGION_3E
    assert_close(1/iapws97_region3_rho(T=661, P=26e6), 0.0029702259620031472, rtol=1e-13)
    assert iapws97_region_3(675.0, 30e6) == REGION_3E
    assert_close(1/iapws97_region3_rho(T=675, P=30e6), 0.0030046270863580073, rtol=1e-13)
    
    assert iapws97_region_3(671.0, 26e6) == REGION_3F
    assert_close(1/iapws97_region3_rho(T=671, P=26e6), 0.00501902940104471, rtol=1e-13)
    assert iapws97_region_3(690.0, 30e6) == REGION_3F
    assert_close(1/iapws97_region3_rho(T=690, P=30e6), 0.004656470141685632, rtol=1e-13)
    
    assert iapws97_region_3(649.0, 23.6e6) == REGION_3G
    assert_close(1/iapws97_region3_rho(T=649, P=23.6e6), 0.0021631983783137196, rtol=1e-13)
    assert iapws97_region_3(650.0, 24e6) == REGION_3G
    assert_close(1/iapws97_region3_rho(T=650, P=24e6), 0.0021660441609564836, rtol=1e-13)
    
    assert iapws97_region_3(652.0, 23.6e6) == REGION_3H
    assert_close(1/iapws97_region3_rho(T=652, P=23.6e6), 0.002651081406573861, rtol=1e-13)
    assert iapws97_region_3(654.0, 24e6) == REGION_3H
    assert_close(1/iapws97_region3_rho(T=654, P=24e6), 0.0029678023349397832, rtol=1e-13)
    
    assert iapws97_region_3(653.0, 23.6e6) == REGION_3I
    assert_close(1/iapws97_region3_rho(T=653, P=23.6e6), 0.003273916815935874, rtol=1e-13)
    assert iapws97_region_3(655.0, 24e6) == REGION_3I
    assert_close(1/iapws97_region3_rho(T=655, P=24e6), 0.0035503298636594843, rtol=1e-13)
    
    assert iapws97_region_3(655.0, 23.5e6) == REGION_3J
    assert_close(1/iapws97_region3_rho(T=655, P=23.5e6), 0.004545001141649382, rtol=1e-13)
    assert iapws97_region_3(660.0, 24e6) == REGION_3J
    assert_close(1/iapws97_region3_rho(T=660, P=24e6), 0.005100267703573203, rtol=1e-13)
    
    assert iapws97_region_3(660.0, 23e6) == REGION_3K
    assert_close(1/iapws97_region3_rho(T=660, P=23e6), 0.006109525996886692, rtol=1e-13)
    assert iapws97_region_3(670.0, 24e6) == REGION_3K
    assert_close(1/iapws97_region3_rho(T=670, P=24e6), 0.0064273256447015745, rtol=1e-13)
    
    assert iapws97_region_3(646.0, 22.6e6) == REGION_3L
    assert_close(1/iapws97_region3_rho(T=646, P=22.6e6), 0.0021178608506781027, rtol=1e-13)
    assert iapws97_region_3(646.0, 23e6) == REGION_3L
    assert_close(1/iapws97_region3_rho(T=646, P=23e6), 0.002062374674146725, rtol=1e-13)
    
    assert iapws97_region_3(648.6, 22.6e6) == REGION_3M
    assert_close(1/iapws97_region3_rho(T=648.6, P=22.6e6), 0.002533063780421483, rtol=1e-13)
    assert iapws97_region_3(649.3, 22.8e6) == REGION_3M
    assert_close(1/iapws97_region3_rho(T=649.3, P=22.8e6), 0.0025729717809150347, rtol=1e-13)
    
    assert iapws97_region_3(649, 22.6e6) == REGION_3N
    assert_close(1/iapws97_region3_rho(T=649, P=22.6e6), 0.0029234327109982578, rtol=1e-13)
    assert iapws97_region_3(649.7, 22.8e6) == REGION_3N
    assert_close(1/iapws97_region3_rho(T=649.7, P=22.8e6), 0.0029133114940412745, rtol=1e-13)
    
    assert iapws97_region_3(649.1, 22.6e6) == REGION_3O
    assert_close(1/iapws97_region3_rho(T=649.1, P=22.6e6), 0.003131208996006528, rtol=1e-13)
    assert iapws97_region_3(649.9, 22.8e6) == REGION_3O
    assert_close(1/iapws97_region3_rho(T=649.9, P=22.8e6), 0.003221160277936286, rtol=1e-13)
    
    assert iapws97_region_3(649.4, 22.6e6) == REGION_3P
    assert_close(1/iapws97_region3_rho(T=649.4, P=22.6e6), 0.0037155961864873133, rtol=1e-13)
    assert iapws97_region_3(650.2, 22.8e6) == REGION_3P
    assert_close(1/iapws97_region3_rho(T=650.2, P=22.8e6), 0.0036647547896187177, rtol=1e-13)
    
    assert iapws97_region_3(640, 21.1E6) == REGION_3Q
    assert_close(1/iapws97_region3_rho(T=640, P=21.1E6), 0.001970999271891958, rtol=1e-13)
    assert iapws97_region_3(643, 21.8E6) == REGION_3Q
    assert_close(1/iapws97_region3_rho(T=643, P=21.8E6), 0.002043919160913867, rtol=1e-13)
    
    assert iapws97_region_3(644, 21.1E6) == REGION_3R
    assert_close(1/iapws97_region3_rho(T=644, P=21.1E6), 0.00525100992110033, rtol=1e-13)
    assert iapws97_region_3(648, 21.8E6) == REGION_3R
    assert_close(1/iapws97_region3_rho(T=648, P=21.8E6), 0.00525684474078012, rtol=1e-13)
    
    assert iapws97_region_3(635.0, 19.1e6) == REGION_3S
    assert_close(1/iapws97_region3_rho(T=635, P=19.1e6), 0.0019328290790263667, rtol=1e-13)
    assert iapws97_region_3(638.0, 20e6) == REGION_3S
    assert_close(1/iapws97_region3_rho(T=638, P=20e6), 0.0019853872274726695, rtol=1e-13)
    
    assert iapws97_region_3(626, 17e6) == REGION_3T
    assert_close(1/iapws97_region3_rho(T=626, P=17e6), 0.008483262001139871, rtol=1e-13)
    assert iapws97_region_3(640, 20e6) == REGION_3T
    assert_close(1/iapws97_region3_rho(T=640, P=20e6), 0.006227528101006945, rtol=1e-13)
    
    assert iapws97_region_3(644.6, 21.5e6) == REGION_3U
    assert_close(1/iapws97_region3_rho(T=644.6, P=21.5e6), 0.002268366646629464, rtol=1e-13)
    assert iapws97_region_3(646.1, 22e6) == REGION_3U
    assert_close(1/iapws97_region3_rho(T=646.1, P=22e6), 0.0022963505532551556, rtol=1e-13)
    
    assert iapws97_region_3(648.6, 22.5e6) == REGION_3V
    assert_close(1/iapws97_region3_rho(T=648.6, P=22.5e6), 0.002832373260251989, rtol=1e-13)
    assert iapws97_region_3(647.9, 22.3e6) == REGION_3V
    assert_close(1/iapws97_region3_rho(T=647.9, P=22.3e6), 0.0028114244045568644, rtol=1e-13)
    
    assert iapws97_region_3(647.5, 22.15e6) == REGION_3W
    assert_close(1/iapws97_region3_rho(T=647.5, P=22.15e6), 0.003694032280598682, rtol=1e-13)
    assert iapws97_region_3(648.1, 22.3e6) == REGION_3W
    assert_close(1/iapws97_region3_rho(T=648.1, P=22.3e6), 0.0036222263053987108, rtol=1e-13)
    
    assert iapws97_region_3(648, 22.11e6) == REGION_3X
    assert_close(1/iapws97_region3_rho(T=648, P=22.11e6), 0.004528072648832488, rtol=1e-13)
    assert iapws97_region_3(649, 22.3e6) == REGION_3X
    assert_close(1/iapws97_region3_rho(T=649, P=22.3e6), 0.004556905798876878, rtol=1e-13)
    
    assert iapws97_region_3(646.84, 22e6) == REGION_3Y
    assert_close(1/iapws97_region3_rho(T=646.84, P=22e6), 0.0026983547189247956, rtol=3e-13)
    assert iapws97_region_3(647.05, 22.064e6) == REGION_3Y
    assert_close(1/iapws97_region3_rho(T=647.05, P=22.064e6), 0.0027176556481596707, rtol=3e-12)
    
    assert iapws97_region_3(646.89, 22e6) == REGION_3Z
    assert_close(1/iapws97_region3_rho(T=646.89, P=22e6), 0.003798732962152225, rtol=2e-12)
    assert iapws97_region_3(647.15, 22.064e6) == REGION_3Z
    assert_close(1/iapws97_region3_rho(T=647.15, P=22.064e6), 0.003701940009172692, rtol=2e-12)



def test_iapws97_rho(): 
    assert_close(iapws97_rho(T=330, P=8e5), 985.1049808079207)
    assert_close(iapws97_rho(T=823, P=14e6), 40.39293607288123)
    assert_close(iapws97_rho(T=2000, P=3e7), 32.11456228328856)
    assert_close(iapws97_rho(648.6, 22.5e6), 353.06081088726)

@pytest.mark.CoolProp
@pytest.mark.slow
def test_iapws97_region_3_rho_coolprop():
    from CoolProp.CoolProp import PropsSI
    Ts = linspace(623.15+1e-10, 1073.15, 1000)
    # Do some points near where more region transitions are, and then up to the limit.
    for P_lim in (25.5e6, 100e6):
        def test_Ps(T, N):
             # Do not check too low to the boundary
             # Sometimes CoolProp says a different region
            lower_P = iapws97_boundary_2_3(T)*(1+4e-6)
            if lower_P >= P_lim:
                # No valid points in region 3
                return []
            upper_P = iapws97_boundary_2_3(T)*10.0
            upper_P = min(upper_P, P_lim)
            return logspace(log10(lower_P), log10(upper_P), N)
    
        for T in Ts:
            for P in test_Ps(T, 1000):
                assert iapws97_identify_region_TP(T, P) == 3
                rho_implemented = iapws97_rho(T=T, P=P)
                rho_CoolProp = PropsSI('DMASS','T',T,'P',P,'IF97::Water')
    #            try:
                assert_close(rho_CoolProp, rho_implemented, rtol=1e-10)
    #            except:
    #                print([T, P, 1-rho_CoolProp/rho_implemented])
#test_iapws97_region_3_rho_coolprop()
    
    
@pytest.mark.CoolProp
@pytest.mark.slow
def test_iapws97_region_5_rho_coolprop():
    # Working great!
    from CoolProp.CoolProp import PropsSI
    Ts = linspace(1073.15+1e-10, 2273.15, 100)
    def test_Ps(T, N):
        return logspace(log10(1e-6), log10(50e6), N)

    for T in Ts:
        for P in test_Ps(T, 1000):
            assert iapws97_identify_region_TP(T, P) == 5
            rho_implemented = iapws97_rho(T=T, P=P)
            rho_CoolProp = PropsSI('DMASS','T',T,'P',P,'IF97::Water')
            assert_close(rho_CoolProp, rho_implemented, rtol=1e-10)


def iapws97_dGr_dpi_region2_fastest(tau, pi):
    '''Fastest implementation, maybe near possible. Horner's method in places
    has caused issues however and this has some error in some regions.
    '''
    taut = tau - 0.5
    pi2 = pi*pi
    taut2 = taut*taut
    taut3 = taut*taut2
    taut4 = taut2*taut2
    taut6 = taut4*taut2
    taut8 = taut4*taut4
    taut13 = taut6*taut4*taut3
    taut21 = taut13*taut8
    taut29 = taut21*taut8
    # 53 from 13*13*!3*!3*3
    # 57 from 
    return (pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi*(pi2*pi2*pi2*(pi2*(pi2
        *(pi*(pi*(pi*(pi*taut13*taut13*(taut13*taut*(1.32995316841867198e-15 - 0.0000226487297378903988*taut13*taut4*taut)
        + 1.75410265428146418e-27) - 2.93678005497663003e-14*taut29*taut8*taut2) + 0.0000832192847496054092*taut*taut21*taut29*taut3)
    - 1.24017662339841913e-24*taut21) + taut8*taut8*taut4*(taut13*taut2*(6.1258633752463995e-12
    - 0.0000840049353964159951*taut13) + 1.78371690710842005e-23)) - 6.05920510335077989*taut21*taut21*taut13*taut2) 
    + taut29*(1.71088510070543998*taut21 - 1.29412653835175996e-9)) + taut4*(taut6
    *(-1.00181793795109993e-8*taut4 - 1.02347470959289996e-12) + 1.04069652101739995e-18))
    + 1.78287415218792009e-7*taut13) + taut4*taut4*(9.00496908836719986e-11 - 65.8490727183984035*taut13*taut13*taut2))
    + taut6*taut4*taut*(-0.27262789705017304*taut6*taut6*taut2 - 8.83526622937069987e-6) - 4.13416950269890026e-17)
    + taut3*(taut13*(-143.374451604623999*taut13*taut6 - 0.012702883392812999) - 1.00288598706366e-10))
    + 0.0000114610381688305001*taut6*taut) + taut*(taut*(1.92901490874028006e-6*taut + 5.11628714091400033e-8)
    - 3.15389238237468004e-9)) + taut*(taut2*(taut3*(-0.122004760687946995*taut21*taut8 
    - 0.00451017736264439952) - 0.0000968330317157100001) + 1.31612001853305008e-6) + 6.14452130769269999e-8)
    + taut*(taut*(taut2*(taut3*(-0.0000533490958281740028*taut21*taut8 - 0.0875945913011459965) 
    - 0.00787855544867100029) - 0.000378979750326299998) - 0.000066065283340406)) 
    + taut*(taut*(taut*(-0.0503252787279300021*taut3 - 0.0575812590834320001) - 0.0459960136963650026) 
    - 0.0178348622923579989) - 0.00177317424732129992)


@pytest.mark.CoolProp
@pytest.mark.slow
def test_iapws97_region_2_rho_coolprop():
    from CoolProp.CoolProp import PropsSI
    P_lim = 1e-6
    Ts = linspace(273.15+1e-10,  1073.15-1e-10, 100)
    def test_Ps(T, N):
        upper_P = iapws97_boundary_2_3(T)*(1.0-1e-10)
        if T <= 623.15:
            upper_P = min(Psat_IAPWS(T)*(1.0-1e-10), upper_P)
        
        
        if upper_P < P_lim or upper_P > 100e6:
            # No valid points in region 2
            return []

        return logspace(log10(P_lim), log10(upper_P), N)

    for T in Ts:
        for P in test_Ps(T, 1000):
            assert iapws97_identify_region_TP(T, P) == 2
            rho_implemented = iapws97_rho(T=T, P=P)
            rho_CoolProp = PropsSI('DMASS','T',T,'P',P,'IF97::Water')
#            try:
            assert_close(rho_CoolProp, rho_implemented, rtol=2e-15)
#            except:
#                print([T, P, 1-rho_implemented/rho_CoolProp])


@pytest.mark.CoolProp
@pytest.mark.slow
def test_iapws97_region_1_rho_coolprop():
    from CoolProp.CoolProp import PropsSI
    Ts = linspace(273.15+1e-10,  623.15-1e-10, 500)
    def test_Ps(T, N):
        Psat = Psat_IAPWS(T)*(1+1e-10)
        return logspace(log10(Psat), log10(100e6), N)

    for T in Ts:
        for P in test_Ps(T, 500):
            assert iapws97_identify_region_TP(T, P) == 1
            rho_implemented = iapws97_rho(T=T, P=P)
            rho_CoolProp = PropsSI('DMASS','T',T,'P',P,'IF97::Water')
            try:
                assert_close(rho_CoolProp, rho_implemented, rtol=2e-13)
            except:
                print([T, P, 1-rho_implemented/rho_CoolProp])
