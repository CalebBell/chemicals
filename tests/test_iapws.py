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

import numpy as np
from numpy.testing import assert_allclose
import pytest
import chemicals
from math import *
from chemicals.iapws import *
from chemicals import iapws
from fluids.numerics import assert_close, assert_close1d, assert_close2d, linspace, logspace, derivative
from chemicals.vapor_pressure import Psat_IAPWS
from chemicals.iapws import (REGION_3A, REGION_3B, REGION_3C, REGION_3D, REGION_3E, REGION_3F, REGION_3G,
                             REGION_3H, REGION_3I, REGION_3J, REGION_3K, REGION_3L, REGION_3M, REGION_3N,
                             REGION_3O, REGION_3P, REGION_3Q, REGION_3R, REGION_3S, REGION_3T, REGION_3U,
                             REGION_3V, REGION_3W, REGION_3X, REGION_3Y, REGION_3Z)
from chemicals.iapws import (iapws_97_Trho_err_region1, iapws_97_Trho_err_region2, iapws_97_Trho_err_region5,
                             iapws_97_Prho_err_region1, iapws_97_Prho_err_region2, iapws_97_Prho_err_region5,
                             iapws_97_Prho_err_region3, iapws95_rho_err, iapws95_T_err, iapws95_d2Ar_ddelta2_delta_1)

def make_me_precise():
    import mpmath as mp
    globals()['exp'] = mp.exp
    globals()['log'] = mp.log
    globals()['sqrt'] = mp.sqrt
    mp.mp.dps = 50

    import fluids.numerics
    fluids.numerics.exp = mp.exp
    fluids.numerics.log = mp.log
    return mp

def make_me_float():
    import math
    globals()['exp'] = math.exp
    globals()['log'] = math.log
    globals()['sqrt'] = math.sqrt
    import fluids.numerics
    fluids.numerics.exp = math.exp
    fluids.numerics.log = math.log

### IAPWS Naive Functions
### Regoin 1
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

def iapws97_d2G_dpi2_region1_naive(tau, pi):
    return sum([nis1[i]*lis1[i]*(lis1[i]-1)*(7.1-pi)**(lis1[i]-2)*(tau-1.222)**Jis1[i] for i in range(34)])

def iapws97_dG_dtau_region1_naive(tau, pi):
    return sum([nis1[i]*Jis1[i]*(7.1-pi)**lis1[i]*(tau-1.222)**(Jis1[i]-1) for i in range(34)])

def iapws97_d2G_dtau2_region1_naive(tau, pi):
    return sum([nis1[i]*Jis1[i]*(Jis1[i]-1)*(7.1-pi)**lis1[i]*(tau-1.222)**(Jis1[i]-2) for i in range(34)])

def iapws97_d2G_dpidtau_region1_naive(tau, pi):
    return sum([-nis1[i]*Jis1[i]*lis1[i]*(7.1-pi)**(lis1[i]-1)*(tau-1.222)**(Jis1[i]-1) for i in range(34)])

### Region 2
# Section 2 - ideal gas part
J0is2 = [0., 1., -5., -4., -3., -2., -1., 2., 3.]
J0is2 = [int(i) for i in J0is2]

def iapws97_G0_region2_naive(tau, pi):
    return log(pi) + sum( [n0is2[i]*tau**J0is2[i] for i in range(9)] )

def iapws97_dG0_dtau_region2_naive(tau, pi):
    return sum([n0is2[i]*J0is2[i]*tau**(J0is2[i]-1) for i in range(9)])

def iapws97_d2G0_dtau2_region2_naive(tau, pi):
    return sum([n0is2[i]*J0is2[i]*(J0is2[i]-1)*tau**(J0is2[i]-2) for i in range(9)])

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

def iapws97_d2Gr_dpi2_region2_naive(tau, pi):
    return sum([nis2[i]*lis2[i]*(lis2[i]-1)*pi**(lis2[i]-2)*(tau-0.5)**Jis2[i] for i in range(43)])

def iapws97_dGr_dtau_region2_naive(tau, pi):
    return sum([nis2[i]*pi**lis2[i]*Jis2[i]*(tau-0.5)**(Jis2[i]-1) for i in range(43)])

def iapws97_d2Gr_dtau2_region2_naive(tau, pi):
    return sum([nis2[i]*pi**lis2[i]*Jis2[i]*(Jis2[i]-1)*(tau-0.5)**(Jis2[i]-2) for i in range(43)])

def iapws97_d2Gr_dpidtau_region2_naive(tau, pi):
    return sum([nis2[i]*lis2[i]*pi**(lis2[i]-1)*Jis2[i]*(tau-0.5)**(Jis2[i]-1) for i in range(43)])

### Region 3
# Section 3 - In terms of density
lis3 = [None, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
        3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
Jis3 = [None, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2,
        4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26]
nis3 = [0.10658070028513E1, -0.15732845290239E2, 0.20944396974307E2,
        -0.76867707878716E1, 0.26185947787954E1, -0.28080781148620E1,
        0.12053369696517E1, -0.84566812812502E-2, -0.12654315477714E1,
        -0.11524407806681E1, 0.88521043984318, -0.64207765181607,
        0.38493460186671, -0.85214708824206, 0.48972281541877E1,
        -0.30502617256965E1, 0.39420536879154E-1, 0.12558408424308,
        -0.27999329698710, 0.13899799569460E1, -0.20189915023570E1,
        -0.82147637173963E-2, -0.47596035734923, 0.43984074473500E-1,
        -0.44476435428739, 0.90572070719733, 0.70522450087967,
        0.10770512626332, -0.32913623258954, -0.50871062041158,
        -0.22175400873096E-1, 0.94260751665092E-1, 0.16436278447961,
        -0.13503372241348E-1, -0.14834345352472E-1, 0.57922953628084E-3,
        0.32308904703711E-2, 0.80964802996215E-4, -0.16557679795037E-3,
        -0.44923899061815E-4]

def iapws97_A_region3_naive(tau, delta):
    return nis3[0]*log(delta) + sum([nis3[i]*delta**lis3[i]*tau**Jis3[i] for i in range(1, 40)])

def iapws97_dA_ddelta_region3_naive(tau, delta):
    return (nis3[0]/delta + sum([nis3[i]*lis3[i]*delta**(lis3[i]-1)*tau**Jis3[i] for i in range(1, 40)]))

def iapws97_d2A_ddelta2_region3_naive(tau, delta):
    return (-nis3[0]/delta**2 + sum([nis3[i]*lis3[i]*(lis3[i]-1)*delta**(lis3[i]-2)*tau**Jis3[i] for i in range(1, 40)]))

def iapws97_dA_dtau_region3_naive(tau, delta):
    return (sum([nis3[i]*Jis3[i]*delta**lis3[i]*tau**(Jis3[i]-1) for i in range(1, 40)]))

def iapws97_d2A_dtau2_region3_naive(tau, delta):
    return (sum([nis3[i]*Jis3[i]*(Jis3[i]-1)*delta**lis3[i]*tau**(Jis3[i]-2) for i in range(1, 40)]))

def iapws97_d2A_ddeltadtau_region3_naive(tau, delta):
    return (sum([nis3[i]*lis3[i]*Jis3[i]*delta**(lis3[i]-1)*tau**(Jis3[i]-1) for i in range(1, 40)]))

### Region 5
# Section 5 - ideal gas part
J0is5 = [0., 1., -3., -2., -1., 2.]
n0is5 = [-0.13179983674201E2, 0.68540841634434E1, -0.24805148933466E-1,
        0.36901534980333, -0.31161318213925E1, -0.32961626538917]
# Section 5 - residual part
lis5 = [1, 1, 1, 2, 2, 3]
Jis5 = [1, 2, 3, 3, 9, 7]
nis5 = [0.15736404855259E-2, 0.90153761673944E-3, -0.50270077677648E-2,
        0.22440037409485E-5, -0.41163275453471E-5, 0.37919454822955E-7]

def iapws97_G0_region5_naive(tau, pi):
    return  log(pi) + sum( [n0is5[i]*tau**J0is5[i] for i in range(6)] )

def iapws97_dG0_dtau_region5_naive(tau, pi):
    return sum( [n0is5[i]*J0is5[i]*tau**(J0is5[i]-1) for i in range(6)] )

def iapws97_d2G0_dtau2_region5_naive(tau, pi):
    return sum( [n0is5[i]*J0is5[i]*(J0is5[i]-1)*tau**(J0is5[i]-2) for i in range(6)] )


def iapws97_Gr_region5_naive(tau, pi):
    return sum( [nis5[i]*pi**lis5[i]*tau**Jis5[i] for i in range(6)] )

def iapws97_dGr_dpi_region5_naive(tau, pi):
    return sum( [nis5[i]*lis5[i]*pi**(lis5[i]-1)*tau**Jis5[i] for i in range(6)] )

def iapws97_d2Gr_dpi2_region5_naive(tau, pi):
    return sum( [nis5[i]*lis5[i]*(lis5[i]-1)*pi**(lis5[i]-2)*tau**Jis5[i] for i in range(6)] )

def iapws97_dGr_dtau_region5_naive(tau, pi):
    return sum( [nis5[i]*pi**lis5[i]*Jis5[i]*tau**(Jis5[i]-1) for i in range(6)] )

def iapws97_d2Gr_dtau2_region5_naive(tau, pi):
    return sum( [nis5[i]*pi**lis5[i]*Jis5[i]*(Jis5[i]-1)*tau**(Jis5[i]-2) for i in range(6)] )

def iapws97_d2Gr_dpidtau_region5_naive(tau, pi):
    return sum( [nis5[i]*lis5[i]*pi**(lis5[i]-1)*Jis5[i]*tau**(Jis5[i]-1) for i in range(6)] )



def test_iapws97_G_region1():
    assert_close(iapws97_G_region1(2.0, 10.0), -2.2706245211745117, rtol=1e-12)
    assert_close(iapws97_G_region1(4.0, .01), -0.23648176259353207, rtol=1e-12)


def test_iapws97_d2G_dpi2_region1():
    assert_close(iapws97_d2G_dpi2_region1(4.8, .8), -0.0009193903700879624, rtol=1e-12)

def test_iapws97_dG_dtau_region1():
    assert_close(iapws97_dG_dtau_region1(4.8, .8), 0.12200788755868151, rtol=1e-12)

def test_iapws97_d2G_dtau2_region1():
    assert_close(iapws97_d2G_dtau2_region1(4.8, .8), -0.38995449934370313, rtol=1e-12)

def test_iapws97_d2G_dpidtau_region1():
    assert_close(iapws97_d2G_dpidtau_region1(4.8, .8), 0.02436120286448202, rtol=1e-12)

### Fast equation fuzz tests
# check that floating points are behaving nicely
# The only two constants in this world are death and floating point error.

@pytest.mark.fuzz
@pytest.mark.slow
def test_iapws97_region1_fuzz():
    funcs_naive = [iapws97_dG_dpi_region1_naive, iapws97_G_region1_naive, iapws97_d2G_dpi2_region1_naive,
                   iapws97_dG_dtau_region1_naive, iapws97_d2G_dtau2_region1_naive, iapws97_d2G_dpidtau_region1_naive]
    funcs_fast = [iapws97_dG_dpi_region1, iapws97_G_region1, iapws97_d2G_dpi2_region1,
                  iapws97_dG_dtau_region1, iapws97_d2G_dtau2_region1, iapws97_d2G_dpidtau_region1]
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

def test_iapws97_G0_region2():
    assert_close(iapws97_G0_region2(.7, .2), -8.652623310109474, rtol=1e-12)

def test_iapws97_dG0_dtau_region2():
    assert_close(iapws97_dG0_dtau_region2(.7, .2), 13.987869557187791, rtol=1e-12)

def test_iapws97_d2G0_dtau2_region2():
    assert_close(iapws97_d2G0_dtau2_region2(.7, .2), -9.417242299320637, rtol=1e-12)

def test_iapws97_Gr_region2():
    assert_close(iapws97_Gr_region2(.7, .2), -0.0015296155747292416, rtol=1e-12)

def test_iapws97_dGr_dtau_region2():
    assert_close( iapws97_dGr_dtau_region2(.7, .2), -0.00865815941508857, rtol=1e-12)

def test_iapws97_d2Gr_dtau2_region2():
    assert_close(iapws97_d2Gr_dtau2_region2(.7, .2), -0.0328162552745155, rtol=1e-12)

def test_iapws97_d2Gr_dpidtau_region2():
    assert_close(iapws97_d2Gr_dpidtau_region2(.7, .2), -0.043342202599839647, rtol=1e-12)

def test_iapws97_dGr_dpi_region2():
    assert_close(iapws97_dGr_dpi_region2(.7, .2), -0.0076523073842253275, rtol=1e-12)

def test_iapws97_d2Gr_dpi2_region2():
    assert_close(iapws97_d2Gr_dpi2_region2(.7, .2), -4.239257829323464e-05, rtol=1e-12)


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_region2_fuzz():
    funcs_naive = [iapws97_d2G0_dtau2_region2_naive, iapws97_dG0_dtau_region2_naive, iapws97_G0_region2_naive, iapws97_d2Gr_dpidtau_region2_naive, iapws97_d2Gr_dtau2_region2_naive, iapws97_dGr_dtau_region2_naive, iapws97_d2Gr_dpi2_region2_naive, iapws97_Gr_region2_naive, iapws97_dGr_dpi_region2_naive]
    funcs_fast = [iapws97_d2G0_dtau2_region2, iapws97_dG0_dtau_region2, iapws97_G0_region2,

                  iapws97_d2Gr_dpidtau_region2, iapws97_d2Gr_dtau2_region2, iapws97_dGr_dtau_region2,
                  iapws97_d2Gr_dpi2_region2, iapws97_Gr_region2, iapws97_dGr_dpi_region2]
    atols = [0, 0, 1e-14, 0, 0, 0.0, 3e-18, 0, 0, ]
    rtols = [1e-14, 1e-14, 5e-15, 2e-14, 2e-14, 2e-15, 1e-14, 2e-15, 2e-15]

    N = 200 # tested up to 2000
    P_lim = 1e-6
    Ts = linspace(273.15, 1073.15, N)
    def test_Ps(T, N):
        upper_P = iapws97_boundary_2_3(T)
        if T <= 623.15:
            upper_P = min(Psat_IAPWS(T), upper_P)

        if upper_P < P_lim or upper_P > 100e6:
            # No valid points in region 2
            return []

        return logspace(log10(P_lim), log10(upper_P), N)

    for naive, fast, rtol, atol in zip(funcs_naive, funcs_fast, rtols, atols):
#        print(fast)
        for T in Ts:
            tau = 540.0/T
            for P in test_Ps(T, N):
                pi = P/1E6
                assert_close(naive(tau, pi),
                             fast(tau, pi), rtol=rtol, atol=atol)

#test_iapws97_region2_fuzz()

def test_iapws97_A_region3():
    assert_close(iapws97_A_region3(1.1, .5), -2.2904854361532445, rtol=1e-13)

def test_iapws97_dA_ddelta_region3():
    assert_close(iapws97_dA_ddelta_region3(1.1, .5), 0.4629318841671135, rtol=1e-13)

def test_iapws97_d2A_ddelta2_region3():
    assert_close(iapws97_d2A_ddelta2_region3(1.1, .5), -2.805028150138424, rtol=1e-13)

def test_iapws97_d2A_ddeltadtau_region3():
    assert_close(iapws97_d2A_ddeltadtau_region3(1.1, .5), -6.355921133001381, rtol=1e-13)

def test_iapws97_dA_dtau_region3():
    assert_close(iapws97_dA_dtau_region3(1.1, .5), 6.572433280048887, rtol=1e-13)

def test_iapws97_d2A_dtau2_region3():
    assert_close(iapws97_d2A_dtau2_region3(1.1, .5), -23.97237328970195, rtol=1e-13)


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_region3_fuzz():
    funcs_naive = [iapws97_d2A_ddeltadtau_region3_naive, iapws97_d2A_dtau2_region3_naive, iapws97_dA_dtau_region3_naive, iapws97_d2A_ddelta2_region3_naive, iapws97_dA_ddelta_region3_naive, iapws97_A_region3_naive]
    funcs_fast = [iapws97_d2A_ddeltadtau_region3, iapws97_d2A_dtau2_region3, iapws97_dA_dtau_region3,
                  iapws97_d2A_ddelta2_region3, iapws97_dA_ddelta_region3, iapws97_A_region3]
    atols = [0, 0, 0, 1e-13, 0.0, 0, ]
    rtols = [3e-12, 1e-11, 5e-13, 1e-12, 2e-12, 5e-14]
    N = 500
    Ts = linspace(623.15, 1073.15, N)

    for naive, fast, rtol, atol in zip(funcs_naive, funcs_fast, rtols, atols):
#        print(fast)
        # Do some points near where more region transitions are, and then up to the limit.
        for P_lim in (25.5e6, 100e6):
            def test_Ps(T, N):
                 # Do not check too low to the boundary
                 # Sometimes CoolProp says a different region
                lower_P = iapws97_boundary_2_3(T)
                if lower_P >= P_lim:
                    # No valid points in region 3
                    return []
                upper_P = iapws97_boundary_2_3(T)*10.0
                upper_P = min(upper_P, P_lim)
                return logspace(log10(lower_P), log10(upper_P), N)

            for T in Ts:
                tau = 647.096/T
                for P in test_Ps(T, N):
                    rho = iapws97_rho(T, P)
                    delta = rho/322.0
#                    print(tau, delta)
                    assert_close(naive(tau, delta),
                                 fast(tau, delta), rtol=rtol, atol=atol)
#test_iapws97_region3_fuzz()


def test_iapws97_G0_region5_naive():
    assert_close(iapws97_G0_region5(.9, 1.01), -10.309167223118518, rtol=1e-13)

def test_iapws97_dG0_dtau_region5_naive():
    assert_close(iapws97_dG0_dtau_region5(.9, 1.01), 9.208884308822196, rtol=1e-13)

def test_iapws97_d2G0_dtau2_region5():
    assert_close(iapws97_d2G0_dtau2_region5(.9, 1.01), -6.337757906177518, rtol=1e-13)

def test_iapws97_Gr_region5():
    assert_close(iapws97_Gr_region5(.9, 1.01), -0.0015332877816188876, rtol=1e-13)

def test_iapws97_dGr_dpi_region5():
    assert_close(iapws97_dGr_dpi_region5(.9, 1.01), -0.0015180281714734884, rtol=1e-13)

def test_iapws97_d2Gr_dpi2_region5():
    assert_close(iapws97_d2Gr_dpi2_region5(.9, 1.01), 1.9216694490837294e-07, rtol=1e-13)

def test_iapws97_dGr_dtau_region5():
    assert_close(iapws97_dGr_dtau_region5(.9, 1.01), -0.009119973056785994, rtol=1e-13)

def test_iapws97_d2Gr_dtau2_region5():
    assert_close(iapws97_d2Gr_dtau2_region5(.9, 1.01), -0.025727469083651373, rtol=1e-13)

def test_iapws97_d2Gr_dpidtau_region5():
    assert_close(iapws97_d2Gr_dpidtau_region5(.9, 1.01), -0.009039988008632744, rtol=1e-13)

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_region5_fuzz():
    funcs_naive = [iapws97_d2G0_dtau2_region5_naive, iapws97_dG0_dtau_region5_naive,
                   iapws97_G0_region5_naive, iapws97_d2Gr_dpidtau_region5_naive,
                   iapws97_d2Gr_dtau2_region5_naive, iapws97_dGr_dtau_region5_naive,
                   iapws97_d2Gr_dpi2_region5_naive, iapws97_Gr_region5_naive,
                   iapws97_dGr_dpi_region5_naive]
    funcs_fast = [iapws97_d2G0_dtau2_region5, iapws97_dG0_dtau_region5,
                  iapws97_G0_region5,

                  iapws97_d2Gr_dpidtau_region5,
                  iapws97_d2Gr_dtau2_region5, iapws97_dGr_dtau_region5,
                  iapws97_d2Gr_dpi2_region5, iapws97_Gr_region5,
                  iapws97_dGr_dpi_region5]
    atols = [0, 0, 0, 0, 0, 0, 5e-21, 4e-17, 1e-18]
    rtols = [2e-15, 1e-15, 1e-15, 5e-15, 2e-15, 1e-14, 1e-14, 2e-14, 2e-14]

#    funcs_naive = [iapws97_Gr_region5_naive]
#    funcs_fast = [iapws97_Gr_region5]
#    atols = [4e-17]
#    rtols = [2e-14]

    N = 2000
    Ts = linspace(1073.15, 2273.15, N)
    def test_Ps(T, N):
        return logspace(log10(1e-6), log10(50e6), N)

    for naive, fast, rtol, atol in zip(funcs_naive, funcs_fast, rtols, atols):
        errs = []
        erri = 0.0
#        print(naive)
        for T in Ts:
            tau = 1000.0/T
            for P in test_Ps(T, N):
                pi = P/1E6
#                print(tau, pi)
                v0 = naive(tau, pi)
                v1 = fast(tau, pi)
                assert_close(v0, v1, rtol=rtol, atol=atol)
                error = abs(1.0 - v1/v0)
                erri += error
                errs.append(error)
#        print(naive, erri/N**2, np.std(errs), np.max(errs))
#test_iapws97_region5_fuzz()
### Fast tests

def test_iapws97_dG_dpi_region1():
    assert_close( iapws97_dG_dpi_region1(4.8, .8), 0.12341682293659642, rtol=1e-12)

    assert_close(iapws97_dG_dpi_region1_naive(1386/277.15, 101325/16.53E6),
                 iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6), rtol=1e-14)

    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 0.12923271825448354, rtol=1e-14)

    # Point that had bad error with horner's method
    assert_close(iapws97_dG_dpi_region1(1386 / 600.15, 10001325 / 16.53E6),
                 0.09345587583404263, rtol=1e-14)
    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 iapws97_dG_dpi_region1_naive(1386/277.15, 101325/16.53E6), rtol=1e-14)



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

    # Vapor pressure boundary
    assert_close(iapws97_rho(432.0135947190398, 600559.0434678708, True), 3.171654556869339)
    assert_close(iapws97_rho(432.0135947190398, 600559.0434678708, False), 908.5584542274903)

@pytest.mark.CoolProp
@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_region_3_rho_coolprop():
    from CoolProp.CoolProp import PropsSI
    Ts = linspace(623.15+1e-10, 1073.15, 100)
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
            for P in test_Ps(T, 100):
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
@pytest.mark.fuzz
def test_iapws97_region_5_rho_coolprop():
    # Working great!
    from CoolProp.CoolProp import PropsSI
    Ts = linspace(1073.15+1e-10, 2273.15, 100)
    def test_Ps(T, N):
        return logspace(log10(1e-6), log10(50e6), N)

    for T in Ts:
        for P in test_Ps(T, 100):
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
@pytest.mark.fuzz
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
        for P in test_Ps(T, 100):
            assert iapws97_identify_region_TP(T, P) == 2
            rho_implemented = iapws97_rho(T=T, P=P)
            rho_CoolProp = PropsSI('DMASS','T',T,'P',P,'IF97::Water')
#            try:
            assert_close(rho_CoolProp, rho_implemented, rtol=2e-15)
#            except:
#                print([T, P, 1-rho_implemented/rho_CoolProp])


@pytest.mark.CoolProp
@pytest.mark.slow
@pytest.mark.fuzz
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


def test_iapws97_P():
    rho = iapws97_rho(273.15, 999)
    assert_close(iapws97_P(273.15, rho), 999, rtol=1e-10)

    rho = iapws97_rho(289.47653061224486, 145.6348477501249)
    assert_close(iapws97_P(289.47653061224486, rho), 145.6348477501249, rtol=1e-10)

    iapws97_identify_region_TP(1032.9489949748788, 1702.7691722259083)
    rho = iapws97_rho(1032.9489949748788, 1702.7691722259083)
    P_calc = iapws97_P(1032.9489949748788, rho)
    assert_close(P_calc, 1702.7691722259083, rtol=1e-10)

    rho = iapws97_rho(273.9508008008008, 696.3744730627147)
    assert_close(iapws97_P(273.9508008008008, rho), 696.3744730627147, rtol=5e-9)


    rho = iapws97_rho(275.5524024024024, 749.6781874965719)
    assert_close(iapws97_P(275.5524024024024, rho), 749.6781874965719, rtol=5e-9)

    rho = iapws97_rho(1500, 20e6)
    assert_close(iapws97_P(1500, rho), 20e6, rtol=5e-9)


    T_spec = 300
    rho_spec = .02
    for i in range(15):
        P_calc = iapws97_P(T_spec, rho_spec)
        assert_close(iapws97_rho(T_spec, P_calc), rho_spec, rtol=1e-10)
        rho_spec *= .25


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_P_fuzz():
    N = 40
    Ts = linspace(273.15, 623.15, N)
    # Ts = linspace(273.15, 1073.15, N)
    Ps = logspace(log10(1e-5), log10(100e6), N)
    for T in Ts:
        for P in Ps:
            rho = iapws97_rho(T, P)
            P_calc = iapws97_P(T, rho)
            assert_close(P, P_calc, rtol=1e-9)


    # Region 1 and 2 general - Good, working great!
    N = 100
    Ts = linspace(273.15, 1073.15, N)
    Ps = logspace(log10(1e-8), log10(100e6), N)
    for T in Ts:
        for P in Ps:
            if iapws97_identify_region_TP(T, P) != 3:
                rho = iapws97_rho(T, P)
                P_calc = iapws97_P(T, rho)
                # 5e-9 is best solvers can do
                assert_close(P, P_calc, rtol=5e-9)

    # Region 5 - works great
    N = 100
    Ts = linspace(1073.15, 2273.15, N)
    Ps = logspace(log10(1e-8), log10(50e6), N)
    for T in Ts:
        for P in Ps:
            rho = iapws97_rho(T, P)
            P_calc = iapws97_P(T, rho)
            assert_close(P, P_calc, rtol=1e-9)


def test_iapws_97_Trho_err_region():
    from chemicals.iapws import iapws_97_Trho_err_region1, iapws_97_Trho_err_region2, iapws_97_Trho_err_region5
    drho_dP_num = derivative(lambda P, *args: iapws_97_Trho_err_region1(P, *args)[0], 1e5, args=(400, 940), dx=1e-1)
    rho_err, drho_dP_analytical = iapws_97_Trho_err_region1(1e5, T=400.0, rho=940)
    assert_close(drho_dP_num, drho_dP_analytical)
    assert_close(drho_dP_analytical, 5.139076806770276e-07)
    assert_close(rho_err, -2.590879183496895)

    drho_dP_num = derivative(lambda P, *args: iapws_97_Trho_err_region2(P, *args)[0], 1e5, args=(400, .5), dx=1e-1)
    rho_err, drho_dP_analytical = iapws_97_Trho_err_region2(1e5, T=400.0, rho=.5)
    assert_close(drho_dP_num, drho_dP_analytical)
    assert_close(drho_dP_analytical, 5.5373377906291985e-06,)
    assert_close(rho_err, 0.04758348314889638)

    drho_dP_num = derivative(lambda P, *args: iapws_97_Trho_err_region5(P, *args)[0], 1e6, args=(2200, 50), dx=1)
    rho_err, drho_dP_analytical = iapws_97_Trho_err_region5(1e6, T=2200, rho=50)
    assert_close(drho_dP_num, drho_dP_analytical)
    assert_close(drho_dP_analytical, 9.84028484195585e-07)
    assert_close(rho_err, -49.01554810610934)


def test_iapws_97_Prho_err_region():
    from chemicals.iapws import iapws_97_Prho_err_region3, iapws_97_Prho_err_region2, iapws_97_Prho_err_region5, iapws_97_Prho_err_region1
    drho_dP_num = derivative(lambda T, *args: iapws_97_Prho_err_region2(T, *args)[0], 400, args=(1e5, 3), dx=1e-5)
    rho_err, drho_dP_analytical = iapws_97_Prho_err_region2(400, P=1e5, rho=3)
    assert_close(drho_dP_analytical, -0.0014400334536077983)
    assert_close(rho_err, -2.4524165168511036)
    assert_close(drho_dP_num, drho_dP_analytical)

    drho_dP_num = derivative(lambda T, *args: iapws_97_Prho_err_region5(T, *args)[0], 2000, args=(1e6, 300), dx=1e-2)
    rho_err, drho_dP_analytical = iapws_97_Prho_err_region5(2000, P=1e6, rho=3)
    assert_close(drho_dP_num, drho_dP_analytical)
    assert_close(rho_err, -1.9170536728150382, rtol=1e-10)
    assert_close(drho_dP_analytical, -0.0005418228178000102, rtol=1e-10)

    drho_dP_num = derivative(lambda T, *args: iapws_97_Prho_err_region1(T, *args)[0], 300, args=(1e6, 990), dx=1e-4)
    rho_err, drho_dP_analytical = iapws_97_Prho_err_region1(300, P=1e6, rho=990)
    assert_close(drho_dP_num, drho_dP_analytical)
    assert_close(rho_err, 6.960320342238447, rtol=1e-10)
    assert_close(drho_dP_analytical, -0.2744655492373509, rtol=1e-10)

    dP_dT_numerical = derivative(lambda T, *args: iapws_97_Prho_err_region3(T, *args)[0], 620, args=(40e6, 400), dx=.01, order=15)
    P_err, dP_dT_analytical = iapws_97_Prho_err_region3(620, P=40e6, rho=400)
    assert_close(dP_dT_numerical, dP_dT_analytical)
    assert_close(dP_dT_analytical, 215319.4089751701)
    assert_close(P_err, -25505787.520883154)


def test_iapws97_T():
    # region 5
    assert_close(iapws97_T(1e7, iapws97_rho(T=1600, P=1e7)), 1600)

    # region 2 top
    assert_close(iapws97_T(60e6, iapws97_rho(T=1000, P=60e6)), 1000)

    # region 3 top
    rho = iapws97_rho(T=640, P=60e6)
    P = iapws97_P(640, rho)
    assert_close(iapws97_T(P, rho), 640)

    # region 1
    rho = iapws97_rho(T=400, P=40e6)
    assert_close(iapws97_T(40e6, rho), 400)

    # region 2 bottom
    rho = iapws97_rho(T=800, P=1e5)
    iapws97_T(1e5, rho)
    assert_close(iapws97_T(1e5, rho), 800)

    # region 1 bottom
    rho = iapws97_rho(T=300, P=1e5)
    iapws97_T(1e5, rho)
    assert_close(iapws97_T(1e5, rho), 300)

    rho = iapws97_rho(T=273.15, P=1e-08)
    iapws97_T(1e-08, rho)
    assert_close(iapws97_T(1e-08, rho), 273.15)


    rho = iapws97_rho(273.15, 0.0065793322465757635)
    assert_close(iapws97_T(0.0065793322465757635, rho), 273.15)

    rho = iapws97_rho(273.15, 673.4150657750918)
    assert_close(iapws97_T(673.4150657750918, rho), 273.15)

    iapws97_identify_region_TP(273.15, 22570197.196339723)
    rho = iapws97_rho(273.15, 22570197.196339723)
    assert_close(iapws97_T(22570197.196339723, rho), 273.15)

    iapws97_identify_region_TP(1073.15, 68926121.04349865)
    rho = iapws97_rho(1073.15, 68926121.04349865)
    assert_close(iapws97_T(68926121.04349865, rho), 1073.15)

    rho = iapws97_rho(273.9508008008008, 17030650.2925232)
    assert_close(iapws97_T(17030650.2925232, rho), 273.9508008008008)

    # region 5 border requiring calc
    rho = iapws97_rho(1073.150000000001, 34705199.859195136)
    assert_close(iapws97_T(34705199.859195136, rho), 1073.150000000001)

    # region 5 border requiring equation 2 calc
    rho = iapws97_rho(1073.15, 52396013.53002634)
    assert_close(iapws97_T(52396013.53002634, rho), 1073.15)

def test_iapws97_identify_region_TP():
    assert 1 == iapws97_identify_region_TP(432.0135947190398, 600559.0434678708)
    assert 2 == iapws97_identify_region_TP(432.0135947190398, 600559.0434678708, use_95_boundary=True)



@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws97_T_fuzz():
    # region 2 and 1
    N = 100
    Ts = linspace(273.15, 1073.15, N)
    Ps = logspace(log10(1e-8), log10(100e6), N)
    for T in Ts:
        for P in Ps:
            if iapws97_identify_region_TP(T, P) != 3:
                rho = iapws97_rho(T, P)
                T_calc = iapws97_T(P, rho)
                try:
                    # 5e-9 is best solvers can do
                    assert_close(T, T_calc, rtol=5e-9)
                except:
                    # multiple solutions
                    rho_recalc = iapws97_rho(T_calc, P)
                    assert_close(rho, rho_recalc, rtol=5e-9)
    # region 5
    N = 100
    Ts = linspace(1073.15+1e-12, 2273.15, N)
    Ps = logspace(log10(1e-8), log10(50e6), N)
    for T in Ts:
        for P in Ps:
            assert iapws97_identify_region_TP(T, P) == 5
            rho = iapws97_rho(T, P)
            T_calc = iapws97_T(P, rho)
            assert_close(T, T_calc, rtol=5e-9)
### IAPWS95 checks
Tc = 647.096
rhoc = 322.
R = 461.51805 # kJ/kg water/K
MW = 18.015268


### Ideal part functions

ni0s = [-8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315,
        1.2795, 0.96956, 0.24873]
gammais = [None, None, None, 1.28728967, 3.53734222, 7.74073708,
           9.24437796, 27.5075105]


def dAddelta_idg(tau, delta):
    '''
    >>> dAddelta_idg(1.000148377125193, 1.1118012422360248)
    0.8994413407821229
    '''
    tot = 1./delta
    return tot


def ddAdddelta_idg(tau, delta):
    '''
    >>> ddAdddelta_idg(1.000148377125193, 1.1118012422360248)
    -0.8089947255079429
    '''
    tot = -1./(delta*delta)
    return tot


def dAdtau_idg(tau, delta):
    '''
    >>> dAdtau_idg(1.000148377125193, 1.1118012422360248)
    9.803439179390017
    '''
    tot =ni0s[1] + ni0s[2]/tau
    for i in range(3, 8):
        tot += ni0s[i]*gammais[i]*(1.0/(1 - exp(-gammais[i]*tau)) - 1.0)
    return tot


def ddAddtau_idg(tau, delta):
    '''
    >>> ddAddtau_idg(1.000148377125193, 1.1118012422360248)
    -3.433163341430695
    '''
    tot = -ni0s[2]/tau**2
    for i in range(3,8):
        tot -= ni0s[i]*gammais[i]**2*exp(-gammais[i]*tau)*(1-exp(-gammais[i]*tau))**-2
    return tot


def A_idg(tau, delta):
    '''
    >>> A_idg(1.000148377125193, 1.1118012422360248)
    -1.5631960505251727
    '''
    tot = log(delta) + ni0s[0] + ni0s[1]*tau + ni0s[2]*log(tau)
    for i in range(3,8):
        tot += ni0s[i]*log(1 - exp(-gammais[i]*tau))
    return tot


### Residual part functions
cis = [None, None, None, None, None, None, None, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
       2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6, None, None, None]
dis = [1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13,
       15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12,
       3, 4, 4, 5, 14, 3, 6, 6, 6, 3, 3, 3]
tis = [-0.5, 0.875, 1., 0.5, 0.75, 0.375, 1., 4., 6., 12., 1., 5., 4., 2.,
       13., 9., 3., 4., 11., 4., 13., 1., 7., 1., 9., 10., 10., 3., 7.,
       10., 10., 6., 10., 10., 1., 2., 3., 4., 8., 6., 9., 8., 16., 22.,
       23., 23., 10., 50., 44., 46., 50., 0., 1., 4.]
nis = [0.12533547935523E-1, 0.78957634722828E1, -0.87803203303561E1,
       0.31802509345418, -0.26145533859358, -0.78199751687981E-2,
       0.88089493102134E-2, -0.66856572307965, 0.20433810950965,
       -0.66212605039687E-4, -0.19232721156002, -0.25709043003438,
       0.16074868486251, -0.40092828925807E-1, 0.39343422603254E-6,
       -0.75941377088144E-5, 0.56250979351888E-3, -0.15608652257135E-4,
       0.11537996422951E-8, 0.36582165144204E-6, -0.13251180074668E-11,
       -0.62639586912454E-9, -0.10793600908932, 0.17611491008752E-1,
       0.22132295167546, -0.40247669763528, 0.58083399985759,
       0.49969146990806E-2, -0.31358700712549E-1, -0.74315929710341,
       0.47807329915480, 0.20527940895948E-1, -0.13636435110343,
       0.14180634400617E-1, 0.83326504880713E-2, -0.29052336009585E-1,
       0.38615085574206E-1, -0.20393486513704E-1, -0.16554050063734E-2,
       0.19955571979541E-2, 0.15870308324157E-3, -0.16388568342530E-4,
       0.43613615723811E-1, 0.34994005463765E-1, -0.76788197844621E-1,
       0.22446277332006E-1, -0.62689710414685E-4, -0.55711118565645E-9,
       -0.19905718354408, 0.31777497330738, -0.11841182425981,
       -0.31306260323435E2, 0.31546140237781E2, -0.25213154341695E4,
       -0.14874640856724, 0.31806110878444]
alphas = [20., 20., 20.]
betas = [150., 150., 250., 0.3, 0.3]
gammas = [1.21, 1.21, 1.25]
epsilons = [1., 1., 1.]
ais = [3.5, 3.5]
bis = [0.85, 0.95]
Bis = [0.2, 0.2]
Cis = [28., 32.]
Dis = [700., 800.]
Ais = [0.32, 0.32]

for arr in (cis, dis, tis, nis):
    for i in range(len(arr)):
        try:
            arr[i] = float(arr[i])
        except:
            pass


def calcA_res(tau, delta):
    '''
    >>> calcA_res(647.096/647., 358./322.)
    -1.212026565041463
    '''
    phir = 0
    for i in range(7):
        phir += nis[i]*delta**dis[i]*tau**tis[i]
    for i in range(7,51):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(-delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(
        -alphas[i-51]*(delta-epsilons[i-51])**2 - betas[i-51]*(tau-gammas[i-51])**2)
    for i in range(2):
        theta = (1-tau) + Ais[i]*((delta-1.0)**2)**(1/(2*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
        Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]
        phir += nis[i+54]*Delta**bis[i]*delta*psi
    return phir




def dAddelta_res(tau, delta):
    '''
    >>> dAddelta_res(647.096/647., 358./322.)
    -0.714012024371285
    '''
    phir = 0.0
    for i in range(7):
        phir += nis[i]*dis[i]*delta**(dis[i]-1.0)*tau**tis[i]
    for i in range(7,51):
        phir += nis[i]*exp(-delta**cis[i])*(delta**(dis[i]-1.0)*tau**tis[i])*(dis[i]-cis[i]*delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(
        -alphas[i-51]*(delta-epsilons[i-51])**2.0 - betas[i-51]*(tau-gammas[i-51])**2.0)*(
        dis[i]/delta - 2.0*alphas[i-51]*(delta-epsilons[i-51]))
    for i in range(2):
        theta = (1.0-tau) + Ais[i]*((delta-1.0)**2.0)**(1.0/(2.0*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1.0)**2.0 - Dis[i]*(tau-1.0)**2.0)
        Delta = theta**2.0 + Bis[i]*((delta-1.0)**2.0)**ais[i]
        _d_psi_d_delta = d_psi_d_delta(i, tau, delta)
        _d_Delta_bd_delta = d_Delta_bd_delta(i, tau, delta)

        phir += nis[i+54]*(Delta**bis[i]*(psi + delta*_d_psi_d_delta)
        + _d_Delta_bd_delta*psi*delta)
    return phir


def ddAdddelta_res(tau, delta):
    '''Works. Don't touch it.
    >>> ddAdddelta_res(647.096/647., 358./322.)
    0.47573069564568893
    '''
    # need this for rho solver newton
    phir = 0
    for i in range(7):
        phir += nis[i]*dis[i]*(dis[i]-1)*delta**(dis[i]-2)*tau**tis[i]
    for i in range(7,51):
        phir += nis[i]*exp(-delta**cis[i])*(delta**(dis[i]-2))*tau**tis[i]*(
        (dis[i]-cis[i]*delta**cis[i])*(dis[i]-1-cis[i]*delta**cis[i]) - cis[i]**2*delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*tau**tis[i]*exp(-alphas[i-51]*(delta-epsilons[i-51])**2
        - betas[i-51]*(tau-gammas[i-51])**2)*(
        - 2*alphas[i-51]*delta**dis[i]
        + 4*alphas[i-51]**2*delta**dis[i]*(delta-epsilons[i-51])**2
        - 4*dis[i]*alphas[i-51]*delta**(dis[i]-1)*(delta-epsilons[i-51])
        + dis[i]*(dis[i]-1)*delta**(dis[i]-2))
    for i in range(2):
        theta = (1-tau) + Ais[i]*((delta-1)**2)**(1/(2*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
        Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]
        _d_psi_d_delta = d_psi_d_delta(i, tau, delta)
        _d_Delta_bd_delta = d_Delta_bd_delta(i, tau, delta)
        _d2_psi_d2_delta = d2_psi_d2_delta(i, tau, delta)
        _d2_Delta_bd2_delta = d2_Delta_bd2_delta(i, tau, delta)

        phir += nis[i+54]*(Delta**bis[i]*(2*_d_psi_d_delta + delta*_d2_psi_d2_delta)
        + 2*_d_Delta_bd_delta*(psi + delta*_d_psi_d_delta) + _d2_Delta_bd2_delta*delta*psi)
    return phir

def dddA_ddddelta_res(tau, delta):
#    return (0.00401214146653984007*delta**23*tau**10*exp(-delta**4) + 25.576954040118963*delta**21*tau**50*exp(-delta**6) - 68.6393942343940893*delta**21*tau**46*exp(-delta**6) + 42.9963516455212797*delta**21*tau**44*exp(-delta**6) - 0.0511548036983829613*delta**19*tau**10*exp(-delta**4) + 1.20336016101793195e-7*delta**18*tau**50*exp(-delta**6) - 140.673247220654275*delta**15*tau**50*exp(-delta**6) + 377.516668289167455*delta**15*tau**46*exp(-delta**6) - 236.479934050367035*delta**15*tau**44*exp(-delta**6) + 0.17001449464462573*delta**15*tau**10*exp(-delta**4) + 0.000131108546740240002*delta**15*tau**8*exp(-delta**2) + 6.26395869124539993e-10*delta**15*tau*exp(-delta) - 2.81878141106042993e-8*delta**14*tau*exp(-delta) + 1.32511800746680002e-12*delta**13*tau**13*exp(-delta) - 0.00126962466593255998*delta**13*tau**9*exp(-delta**2) - 0.00255661666143467987*delta**13*tau**8*exp(-delta**2) - 0.0159644575836327997*delta**13*tau**6*exp(-delta**2) + 3.94629397548460203e-7*delta**13*tau*exp(-delta) - 4.81344064407172886e-7*delta**12*tau**50*exp(-delta**6) - 5.1679602291205205e-11*delta**12*tau**13*exp(-delta) + 0.0132432400509872004*delta**12*tau**8*exp(-delta**2) + 0.163147892109631987*delta**12*tau**4*exp(-delta**2) - 0.308920684593648021*delta**12*tau**3*exp(-delta**2) + 0.232418688076680008*delta**12*tau**2*exp(-delta**2) - 0.0666612039045704069*delta**12*tau*exp(-delta**2) - 1.71006072270999421e-6*delta**12*tau*exp(-delta) - 0.606049487964162026*delta**11*tau**23*exp(-delta**3) + 6.20155227494462486e-10*delta**11*tau**13*exp(-delta) - 0.136914327545672065*delta**11*tau**10*exp(-delta**4) + 0.0209488069878872411*delta**11*tau**9*exp(-delta**2) + 0.0141597230479459206*delta**11*tau**8*exp(-delta**2) + 0.263413550129941165*delta**11*tau**6*exp(-delta**2) - 3.65821651442040008e-7*delta**11*tau**4*exp(-delta) + 2.07328134180476731*delta**10*tau**23*exp(-delta**3) - 0.94483814752165507*delta**10*tau**22*exp(-delta**3) - 2.27390250081302894e-9*delta**10*tau**13*exp(-delta) - 1.15379964229510002e-9*delta**10*tau**11*exp(-delta) - 0.113445075204936005*delta**10*tau**10*exp(-delta**2) - 0.198648600764808003*delta**10*tau**8*exp(-delta**2) - 2.44721838164447991*delta**10*tau**4*exp(-delta**2) + 0.0000120721144975873203*delta**10*tau**4*exp(-delta) + 4.63381026890472025*delta**10*tau**3*exp(-delta**2) - 3.48628032115020003*delta**10*tau**2*exp(-delta**2) + 0.99991805856855609*delta**10*tau*exp(-delta**2) + 142.09418911177201*delta**9*tau**50*exp(-delta**6) - 381.329967968855954*delta**9*tau**46*exp(-delta**6) + 238.868620252896022*delta**9*tau**44*exp(-delta**6) - 1.17756762454289698*delta**9*tau**16*exp(-delta**3) + 3.46139892688530034e-8*delta**9*tau**11*exp(-delta) + 1.09091480882744007*delta**9*tau**10*exp(-delta**2) - 0.0952218499449419969*delta**9*tau**9*exp(-delta**2) - 0.0216329102121396027*delta**9*tau**8*exp(-delta**2) - 1.36155784594004392*delta**9*tau**6*exp(-delta**2) - 0.000105112492718738184*delta**9*tau**4*exp(-delta) + 4.24234641574913418*delta**8*tau**23*exp(-delta**3) - 3.11525903419677031e-7*delta**8*tau**11*exp(-delta) - 2.4632454907791681*delta**8*tau**10*exp(-delta**2) + 0.804526833097472416*delta**8*tau**8*exp(-delta**2) + 9.911234445660142*delta**8*tau**4*exp(-delta**2) - 0.0000592701760150253563*delta**8*tau**4*exp(-delta) - 18.766931589064118*delta**8*tau**3*exp(-delta**2) + 14.11943530065831*delta**8*tau**2*exp(-delta**2) - 4.0496681372026524*delta**8*tau*exp(-delta**2) - 12.4396880508286038*delta**7*tau**23*exp(-delta**3) + 5.6690288851299302*delta**7*tau**22*exp(-delta**3) + 8.30735742452472082e-7*delta**7*tau**11*exp(-delta) - 5.50933111586084046*delta**7*tau**10*exp(-delta**2) + 0.114266219933930394*delta**7*tau**9*exp(-delta**2) + 0.250869605700392018*delta**7*tau**7*exp(-delta**2) + 3.16114821778658417*delta**7*tau**6*exp(-delta**2) + 0.00337146888754115935*delta**7*tau**4*exp(-delta) - 0.0399753175926447976*delta**7*tau**3*exp(-delta**2) - 0.000562509793518880044*delta**7*tau**3*exp(-delta) + 2.77441370456912092e-7*delta**6*tau**50*exp(-delta**6) + 5.88783812271448515*delta**6*tau**16*exp(-delta**3) + 25.6054990265034803*delta**6*tau**10*exp(-delta**2) - 0.834324123212193625*delta**6*tau**8*exp(-delta**2) - 10.2783172029068144*delta**6*tau**4*exp(-delta**2) - 0.00786676073759603849*delta**6*tau**4*exp(-delta) + 19.4620031293998252*delta**6*tau**3*exp(-delta**2) + 0.0118127056638964809*delta**6*tau**3*exp(-delta) - 14.6423773488308395*delta**6*tau**2*exp(-delta**2) + 4.19965584598793562*delta**6*tau*exp(-delta**2) - 6.19517254363365666*delta**5*tau**23*exp(-delta**3) - 11.9150444067814796*delta**5*tau**10*exp(-delta**2) - 1.77058361340368009*delta**5*tau**9*exp(-delta**2) + 7.59413770881439999e-6*delta**5*tau**9*exp(-delta) - 1.88152204275294022*delta**5*tau**7*exp(-delta**2) - 4.43403523352476814*delta**5*tau**6*exp(-delta**2) + 0.299814881944835965*delta**5*tau**3*exp(-delta**2) - 0.070876233983378889*delta**5*tau**3*exp(-delta) - 0.140891928070015993*delta**5*tau*exp(-delta**2) + 14.2826047990995058*delta**4*tau**23*exp(-delta**3) - 6.50888501626029026*delta**4*tau**22*exp(-delta**3) - 3.93434226032540015e-7*delta**4*tau**13*exp(-delta) - 40.8530296559261075*delta**4*tau**10*exp(-delta**2) - 0.000113912065632216001*delta**4*tau**9*exp(-delta) + 0.863488072714559962*delta**4*tau**7*exp(-delta**2) + 0.118127056638964806*delta**4*tau**3*exp(-delta) + 0.0400928289258069975*delta**4*tau**2*exp(-delta) - 14.2094189111772007*delta**3*tau**50*exp(-delta**6) + 38.1329967968856067*delta**3*tau**46*exp(-delta**6) - 23.8868620252896022*delta**3*tau**44*exp(-delta**6) - 4.97195219251445408*delta**3*tau**16*exp(-delta**3) + 4.72121071239048018e-6*delta**3*tau**13*exp(-delta) + 40.4904092746456783*delta**3*tau**10*exp(-delta**2) + 7.9676262603165604*delta**3*tau**9*exp(-delta**2) + 0.000455648262528864003*delta**3*tau**9*exp(-delta) + 3.01043526840470443*delta**3*tau**7*exp(-delta**2) + 2.46335290751376013*delta**3*tau**6*exp(-delta**2) + 161364187.786847979*delta**3*tau**4*(delta - 1.0)**3*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) - 12102314.0840135999*delta**3*tau**4*(delta - 1.0)*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) - 0.160748684862510011*delta**3*tau**4*exp(-delta) - 0.479703811111737544*delta**3*tau**3*exp(-delta**2) - 0.48111394710968397*delta**3*tau**2*exp(-delta) - 2018952.97521798406*delta**3*tau*(delta - 1.0)**3*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 151421.473141348804*delta**3*tau*(delta - 1.0)*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 0.634013676315071995*delta**3*tau*exp(-delta**2) + 2003600.66069984017*delta**3*(delta - 1.0)**3*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) - 150270.049552488024*delta**3*(delta - 1.0)*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 1.34677663992036001*delta**2*tau**23*exp(-delta**3) - 0.0000141636321371714414*delta**2*tau**13*exp(-delta) - 2.68063804302186082*delta**2*tau**10*exp(-delta**2) - 0.000455648262528864003*delta**2*tau**9*exp(-delta) - 2.59046421814367989*delta**2*tau**7*exp(-delta**2) + 0.25709043003437998*delta**2*tau**5*exp(-delta) - 36306942.2520408034*delta**2*tau**4*(delta - 1.0)**2*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) + 907673.556301020086*delta**2*tau**4*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) + 1.44673816376259001*delta**2*tau**4*exp(-delta) + 1.44334184132905197*delta**2*tau**2*exp(-delta) + 454264.419424046355*delta**2*tau*(delta - 1.0)**2*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) - 11356.6104856011589*delta**2*tau*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 0.192327211560020001*delta**2*tau*exp(-delta) - 450810.148657464015*delta**2*(delta - 1.0)**2*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 11270.2537164366004*delta**2*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) - 1.84291674827090413*delta*tau**23*exp(-delta**3) + 0.839856131130360062*delta*tau**22*exp(-delta**3) + 9.44242142478096036e-6*delta*tau**13*exp(-delta) + 0.0000662126050396869941*delta*tau**12*exp(-delta) - 8.1763823872351189*delta*tau**10*exp(-delta**2) - 5.31175084021104027*delta*tau**9*exp(-delta**2) - 0.752608817101175998*delta*tau**7*exp(-delta**2) - 0.204338109509650007*delta*tau**6*exp(-delta) - 1.54254258020627999*delta*tau**5*exp(-delta) + 1815347.11260203994*delta*tau**4*(delta - 1.0)*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) - 2.22491060444553002*delta*tau**4*exp(-delta) + 0.119925952777934386*delta*tau**3*exp(-delta**2) - 0.962227894219367941*delta*tau**2*exp(-delta) - 22713.2209712023177*delta*tau*(delta - 1.0)*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 0.211414783445121601*delta*tau - 0.422675784210048033*delta*tau*exp(-delta**2) - 1.15396326936012006*delta*tau*exp(-delta) + 26122.2492869444213*delta*(delta - 1.0)**3*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.849999999999999978*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) - 83377.8113011882378*delta*(delta - 1.0)**3*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.949999999999999956*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) - 1399.40621180059406*delta*(delta - 1.0)*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) + 3908.33490474319842*delta*(delta - 1.0)*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) - 1399.40621180059384*delta*(delta - 1.0)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.849999999999999978*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) + 3908.33490474319842*delta*(delta - 1.0)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.949999999999999956*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) + 22540.5074328732007*delta*(delta - 1.0)*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 61.0677328866124753*delta*(0.106666666666666771*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.0700000000000000622*((delta - 1.0)**2)**3.5)*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.05000000000000004)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0) - 24.9893966392963236*delta*(0.320000000000000062*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.210000000000000048*((delta - 1.0)**2)**3.5)*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.14999999999999991)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0) + 24.9893966392963165*delta*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0) - 61.0677328866124753*delta*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0) + 24.9893966392963236*delta*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*(4.2311111111111126*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.93422222222222251*((delta - 1.0)**2)**3.33333333333333348 + 7.13999999999999879*((delta - 1.0)**2)**3.5)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0) - 61.0677328866124753*delta*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*(4.72888888888889003*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 2.16177777777777802*((delta - 1.0)**2)**3.33333333333333348 + 7.98000000000000043*((delta - 1.0)**2)**3.5)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0) + 0.318061108784439994*delta*(0.106666666666666771*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.0700000000000000622*((delta - 1.0)**2)**3.5)*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*(2.24000000000000021*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.4700000000000002*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-2.04999999999999982)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**3 - 0.636122217568879988*delta*(0.106666666666666771*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.0700000000000000622*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.05000000000000004)*(4.72888888888889003*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 2.16177777777777802*((delta - 1.0)**2)**3.33333333333333348 + 7.98000000000000043*((delta - 1.0)**2)**3.5)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**3 - 0.148746408567240002*delta*(0.320000000000000062*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.210000000000000048*((delta - 1.0)**2)**3.5)*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*(2.45333333333333359*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.60999999999999988*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-2.14999999999999991)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**3 + 0.297492817134480003*delta*(0.320000000000000062*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.210000000000000048*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.14999999999999991)*(4.2311111111111126*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.93422222222222251*((delta - 1.0)**2)**3.33333333333333348 + 7.13999999999999879*((delta - 1.0)**2)**3.5)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**3 + 0.148746408567240002*delta*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.14999999999999991)*(0.746666666666666812*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.341333333333333433*((delta - 1.0)**2)**3.33333333333333348 + 1.26000000000000001*((delta - 1.0)**2)**3.5)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**3 - 0.318061108784439994*delta*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.05000000000000004)*(0.248888888888889104*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.113777777777777894*((delta - 1.0)**2)**3.33333333333333348 + 0.420000000000000373*((delta - 1.0)**2)**3.5)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**3 - 0.148746408567240002*delta*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*(5.6414814814814811*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 13.5395555555555589*((delta - 1.0)**2)**3.33333333333333348 + 35.6999999999999886*((delta - 1.0)**2)**3.5)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**3 + 0.318061108784439994*delta*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*(6.30518518518518789*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 15.132444444444447*((delta - 1.0)**2)**3.33333333333333348 + 39.9000000000000057*((delta - 1.0)**2)**3.5)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**3 - 0.046919851012788602*tau**0.375 - 3.34266711393870005e-9*tau**50*exp(-delta**6) + 0.261681694342865978*tau**16*exp(-delta**3) - 0.000198637815119060996*tau**12*exp(-delta) + 3.48500399914553993*tau**10*exp(-delta**2) + 0.647616054535919972*tau**7*exp(-delta**2) + 0.613014328528949992*tau**6*exp(-delta) + 1.54254258020627999*tau**5*exp(-delta) - 15127.8926050169994*tau**4*exp(-20*(delta - 1)**2 - 390.625*(0.8*tau - 1)**2) - 1.04120506006389002*tau**4*exp(-delta) + 189.276841426685991*tau*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) + 1.15396326936012006*tau*exp(-delta) - 1399.40621180059406*(delta - 1.0)**2*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.849999999999999978*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) + 3908.33490474319842*(delta - 1.0)**2*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.949999999999999956*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) + 49.9787932785926472*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) - 122.135465773224951*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) + 24.9893966392963165*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.849999999999999978*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2) - 61.0677328866124753*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**0.949999999999999956*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2) - 187.837561940610016*exp(-20*(delta - 1)**2 - 219.615*(0.826446280991736*tau - 1)**2) - 0.954183326353319927*(0.106666666666666771*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.0700000000000000622*((delta - 1.0)**2)**3.5)*(2.02666666666666684*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.33000000000000007*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.05000000000000004)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**2 + 0.446239225701720033*(0.320000000000000062*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 0.210000000000000048*((delta - 1.0)**2)**3.5)*(1.81333333333333346*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.19000000000000017*((delta - 1.0)**2)**3.5)*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-1.14999999999999991)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**2 - 0.446239225701720033*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.150000000000000022)*(4.2311111111111126*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 1.93422222222222251*((delta - 1.0)**2)**3.33333333333333348 + 7.13999999999999879*((delta - 1.0)**2)**3.5)*exp(-28*(delta - 1)**2 - 700*(tau - 1)**2)/(delta - 1.0)**2 + 0.954183326353319927*((-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)**2 + 0.200000000000000011*((delta - 1.0)**2)**3.5)**(-0.0500000000000000444)*(4.72888888888889003*(-tau + 0.320000000000000007*((delta - 1.0)**2)**1.66666666666666674 + 1.0)*((delta - 1.0)**2)**1.66666666666666674 + 2.16177777777777802*((delta - 1.0)**2)**3.33333333333333348 + 7.98000000000000043*((delta - 1.0)**2)**3.5)*exp(-32*(delta - 1)**2 - 800*(tau - 1)**2)/(delta - 1.0)**2)
    # Generated with sympy as un-expanded reference here; no source to compare against
    phir = 0
    for i in range(7):
        phir += delta**(dis[i] - 3)*dis[i]*nis[i]*tau**tis[i]*(dis[i]**2 - 3*dis[i] + 2)
#        phir += delta**(dis[i] - 2)*dis[i]*nis[i]*tau**tis[i]*(dis[i] - 2)*(dis[i] - 1)/delta
    for i in range(7,51):
        phir += delta**dis[i]*nis[i]*tau**tis[i]*(-3*cis[i]*delta**cis[i]*dis[i]*(dis[i] - 1) + 3*cis[i]*delta**cis[i]*dis[i]*(cis[i]*delta**cis[i] - cis[i] + 1) - cis[i]*delta**cis[i]*(cis[i]**2*delta**(2*cis[i]) - 3*cis[i]**2*delta**cis[i] + cis[i]**2 + 3*cis[i]*delta**cis[i] - 3*cis[i] + 2) + dis[i]*(dis[i]**2 - 3*dis[i] + 2))*exp(-delta**cis[i])/delta**3
#        phir += (-cis[i]*delta**cis[i]*delta**(dis[i] - 2)*nis[i]*tau**tis[i]*(-cis[i]**2*delta**cis[i] + (-cis[i]*delta**cis[i] + dis[i])*(-cis[i]*delta**cis[i] + dis[i] - 1))*exp(-delta**cis[i])/delta + delta**(dis[i] - 2)*nis[i]*tau**tis[i]*(-cis[i]**3*delta**cis[i]/delta - cis[i]**2*delta**cis[i]*(-cis[i]*delta**cis[i] + dis[i])/delta - cis[i]**2*delta**cis[i]*(-cis[i]*delta**cis[i] + dis[i] - 1)/delta)*exp(-delta**cis[i]) + delta**(dis[i] - 2)*nis[i]*tau**tis[i]*(dis[i] - 2)*(-cis[i]**2*delta**cis[i] + (-cis[i]*delta**cis[i] + dis[i])*(-cis[i]*delta**cis[i] + dis[i] - 1))*exp(-delta**cis[i])/delta)
    for i in range(51, 54):
        phir += (delta**dis[i]*nis[i]*tau**tis[i]*(-4*alphas[i-51]**2*(delta - epsilons[i-51])*(2*alphas[i-51]*(delta - epsilons[i-51])**2 - 3) + 6*alphas[i-51]*dis[i]*(2*alphas[i-51]*(delta - epsilons[i-51])**2 - 1)/delta - 6*alphas[i-51]*dis[i]*(delta - epsilons[i-51])*(dis[i] - 1)/delta**2 + dis[i]*(dis[i]**2 - 3*dis[i] + 2)/delta**3)*exp(-alphas[i-51]*(delta - epsilons[i-51])**2 - betas[i-51]*(-gammas[i-51] + tau)**2))
#        phir += (-alphas[i-51]*nis[i]*tau**tis[i]*(2*delta - 2*epsilons[i-51])*(4*alphas[i-51]**2*delta**dis[i]*(delta - epsilons[i-51])**2 - 2*alphas[i-51]*delta**dis[i] - 4*alphas[i-51]*delta**(dis[i] - 1)*dis[i]*(delta - epsilons[i-51]) + delta**(dis[i] - 2)*dis[i]*(dis[i] - 1))*exp(-alphas[i-51]*(delta - epsilons[i-51])**2 - betas[i-51]*(-gammas[i-51] + tau)**2) + nis[i]*tau**tis[i]*(4*alphas[i-51]**2*delta**dis[i]*(2*delta - 2*epsilons[i-51]) + 4*alphas[i-51]**2*delta**dis[i]*dis[i]*(delta - epsilons[i-51])**2/delta - 4*alphas[i-51]*delta**(dis[i] - 1)*dis[i] - 2*alphas[i-51]*delta**dis[i]*dis[i]/delta - 4*alphas[i-51]*delta**(dis[i] - 1)*dis[i]*(delta - epsilons[i-51])*(dis[i] - 1)/delta + delta**(dis[i] - 2)*dis[i]*(dis[i] - 2)*(dis[i] - 1)/delta)*exp(-alphas[i-51]*(delta - epsilons[i-51])**2 - betas[i-51]*(-gammas[i-51] + tau)**2))
    for i in range(2):
        phir += (2*nis[i+54]*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**bis[i]*(-2*Cis[i]**2*delta*(delta - 1)*(2*Cis[i]*(delta - 1)**2 - 3) + 6*Cis[i]*bis[i]*delta*(2*Cis[i]*(delta - 1)**2 - 1)*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])/((delta - 1)*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) - 6*Cis[i]*bis[i]*delta*(Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**2 - Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**2 + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i] - Bis[i]*ais[i]*((delta - 1)**2)**ais[i] + 2*bis[i]*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 2*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/((delta - 1)*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) - 12*Cis[i]*bis[i]*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 3*Cis[i]*(2*Cis[i]*(delta - 1)**2 - 1) + bis[i]*delta*(-3*Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**2 + 3*Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**3 + 2*Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] - 3*Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**2 + Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**3 + 4*Bis[i]*ais[i]**3*((delta - 1)**2)**ais[i] - 6*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i] + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i] + 4*bis[i]**2*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + 6*bis[i]*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])*(Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**2 - Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**2 + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i] - Bis[i]*ais[i]*((delta - 1)**2)**ais[i])/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 12*bis[i]*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 6*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])*(Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**2 - Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**2 + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i] - Bis[i]*ais[i]*((delta - 1)**2)**ais[i])/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 8*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2)/((delta - 1)**3*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) + 3*bis[i]*(Ais[i]**2*((delta - 1)**2)**(1/betas[i+3])/betas[i+3]**2 - Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3]**2 + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i] - Bis[i]*ais[i]*((delta - 1)**2)**ais[i] + 2*bis[i]*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 2*(Ais[i]*(Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1)**2)**(1/(2*betas[i+3]))/betas[i+3] + Bis[i]*ais[i]*((delta - 1)**2)**ais[i])**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/((delta - 1)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1)**2)**(1/(2*betas[i+3])) - tau + 1)**2)))*exp(-Cis[i]*(delta - 1)**2 - Dis[i]*(tau - 1)**2))
#        phir += (nis[i+54]*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**bis[i]*(-4*Cis[i]**2*delta*(delta - 1)*(2*Cis[i]*(delta - 1)**2 - 3) - 6*Cis[i]*bis[i]*delta*(delta - 1)*(2*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + 2*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - 2*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 2*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 4*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2 + bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 6*Cis[i]*bis[i]*delta*(2*Cis[i]*(delta - 1)**2 - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 12*Cis[i]*bis[i]*(delta - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 6*Cis[i]*(2*Cis[i]*(delta - 1)**2 - 1) + bis[i]*delta*(4*Ais[i]**2*(delta - 1)*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + Ais[i]**2*(2*delta - 2.0)*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) - 12*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**5) + 3*Ais[i]**2*(delta - 1)**2*(2*delta - 2.0)*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**3*(delta - 1.0)**6) - 8*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 6*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**4) + 4*Ais[i]*(delta - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) - 12*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**5) + Ais[i]*(delta - 1)**2*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**3*(delta - 1.0)**6) + 8*Bis[i]*ais[i]**3*((delta - 1)**2)**ais[i]/(delta - 1)**3 - 12*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**3 + 4*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**3 + bis[i]**2*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + 6*bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 3*bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 6*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 2*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 2*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**3/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 3*bis[i]*(2*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + 2*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - 2*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 2*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 4*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2 + bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))*exp(-Cis[i]*(delta - 1)**2 - Dis[i]*(tau - 1)**2))
    return phir


def dAdtau_res(tau, delta):
    '''
    >>> dAdtau_res(647.096/647., 358./322.)
    -3.2172250077516558
    '''
    phir = 0
    for i in range(7):
        phir += nis[i]*tis[i]*delta**dis[i]*tau**(tis[i]-1)
    for i in range(7,51):
        phir += nis[i]*tis[i]*delta**dis[i]*tau**(tis[i]-1)*exp(-delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(
        -alphas[i-51]*(delta-epsilons[i-51])**2 - betas[i-51]*(tau-gammas[i-51])**2)*(
        tis[i]/tau -2*betas[i-51]*(tau-gammas[i-51]))
    for i in range(2):
        theta = (1-tau) + Ais[i]*((delta-1)**2)**(1/(2*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
        Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]
        _d_Delta_bd_tau = d_Delta_bd_tau(i, tau, delta)
        _d_psi_d_tau = d_psi_d_tau(i, tau, delta)

        phir += nis[i+54]*delta*(_d_Delta_bd_tau*psi + Delta**bis[i]*_d_psi_d_tau)
    return phir
def ddAddtau_res(tau, delta):
    '''
    >>> ddAddtau_res(647.096/647., 358./322.)
    -9.960295065592888
    '''
    phir = 0
    for i in range(7):
        phir += nis[i]*tis[i]*(tis[i]-1)*delta**dis[i]*tau**(tis[i]-2)
    for i in range(7,51):
        phir += nis[i]*tis[i]*(tis[i]-1)*delta**dis[i]*tau**(tis[i]-2)*exp(-delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(
        -alphas[i-51]*(delta-epsilons[i-51])**2 - betas[i-51]*(tau-gammas[i-51])**2)*(
        (tis[i]/tau -2*betas[i-51]*(tau-gammas[i-51]))**2 - tis[i]/tau**2 -2*betas[i-51])
    for i in range(2):
        theta = (1-tau) + Ais[i]*((delta-1)**2)**(1/(2*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
        Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

        _d2_Delta_bd2_tau = d2_Delta_bd2_tau(i, tau, delta)
        _d_Delta_bd_tau = d_Delta_bd_tau(i, tau, delta)
        _d_psi_d_tau = d_psi_d_tau(i, tau, delta)
        _d2_psi_d2_tau = d2_psi_d2_tau(i, tau, delta)

        phir += nis[i+54]*delta*(_d2_Delta_bd2_tau*psi
        + 2*_d_Delta_bd_tau*_d_psi_d_tau + Delta**bis[i]*_d2_psi_d2_tau)
    return phir


def dAddeltatau_res(tau, delta):
    '''
    >>> dAddeltatau_res(647.096/647., 358./322.)
    -1.332147204361434
    '''
    phir = 0
    for i in range(7):
        phir += nis[i]*dis[i]*tis[i]*delta**(dis[i]-1)*tau**(tis[i]-1)
    for i in range(7,51):
        phir += nis[i]*tis[i]*delta**(dis[i]-1)*tau**(tis[i]-1)*exp(-delta**cis[i])*(dis[i] - cis[i]*delta**cis[i])
    for i in range(51, 54):
        phir += nis[i]*delta**dis[i]*tau**tis[i]*exp(
        -alphas[i-51]*(delta-epsilons[i-51])**2 - betas[i-51]*(tau-gammas[i-51])**2)*(
        dis[i]/delta - 2*alphas[i-51]*(delta-epsilons[i-51]))*(
        (tis[i]/tau - 2*betas[i-51]*(tau-gammas[i-51])))
    for i in range(2):
        theta = (1-tau) + Ais[i]*((delta-1)**2)**(1/(2*betas[i+3]))
        psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
        Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

        _d_psi_d_tau = d_psi_d_tau(i, tau, delta)
        _d2_psi_d_delta_d_tau = d2_psi_d_delta_d_tau(i, tau, delta)
        _d_Delta_bd_delta = d_Delta_bd_delta(i, tau, delta)
        _d_Delta_bd_tau = d_Delta_bd_tau(i, tau, delta)
        _d_psi_d_delta = d_psi_d_delta(i, tau, delta)
        _d2_Delta_bd_delta_d_tau = d2_Delta_bd_delta_d_tau(i, tau, delta)

        phir += nis[i+54]*(Delta**bis[i]*(_d_psi_d_tau + delta*_d2_psi_d_delta_d_tau)
        + delta*_d_Delta_bd_delta*_d_psi_d_tau + _d_Delta_bd_tau*(psi + delta*_d_psi_d_delta)
        + _d2_Delta_bd_delta_d_tau*psi*delta)
    return phir

def iapws95_d3Ar_ddeltadtau2_naive(tau, delta):
    # Not in publication
    phir = 0.0
    for i in range(7):
        phir += delta**dis[i]*dis[i]*nis[i]*tau**tis[i]*tis[i]*(tis[i] - 1)/(delta*tau**2)
    for i in range(7,51):
        phir += delta**dis[i]*nis[i]*tau**tis[i]*tis[i]*(-cis[i]*delta**cis[i]*tis[i] + cis[i]*delta**cis[i] + dis[i]*tis[i] - dis[i])*exp(-delta**cis[i])/(delta*tau**2)
    for i in range(51, 54):
        phir += delta**dis[i]*nis[i]*tau**tis[i]*(-8*alphas[i-51]*betas[i-51]**2*(delta - epsilons[i-51])*(gammas[i-51] - tau)**2 + 4*alphas[i-51]*betas[i-51]*(delta - epsilons[i-51]) - 8*alphas[i-51]*betas[i-51]*tis[i]*(delta - epsilons[i-51])*(gammas[i-51] - tau)/tau - 2*alphas[i-51]*tis[i]**2*(delta - epsilons[i-51])/tau**2 + 2*alphas[i-51]*tis[i]*(delta - epsilons[i-51])/tau**2 + 4*betas[i-51]**2*dis[i]*(gammas[i-51] - tau)**2/delta - 2*betas[i-51]*dis[i]/delta + 4*betas[i-51]*dis[i]*tis[i]*(gammas[i-51] - tau)/(delta*tau) + dis[i]*tis[i]**2/(delta*tau**2) - dis[i]*tis[i]/(delta*tau**2))*exp(-alphas[i-51]*(delta - epsilons[i-51])**2 - betas[i-51]*(-gammas[i-51] + tau)**2)
    for i in range(2):
        phir += 2*nis[i+54]*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**bis[i]*(2*Ais[i]* Dis[i]*bis[i]*delta*(2*delta - 2.0)*(tau - 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) + 2*Ais[i]*bis[i]**2*delta*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2) - 2*Ais[i]*bis[i]*delta*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2) - 4*Cis[i]* Dis[i]**2*delta*(delta - 1)*(tau - 1)**2 - 8*Cis[i]* Dis[i]*bis[i]*delta*(delta - 1)*(tau - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 2*Cis[i]* Dis[i]*delta*(delta - 1) - 4*Cis[i]*bis[i]**2*delta*(delta - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 2*Cis[i]*bis[i]*delta*(delta - 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 4*Cis[i]*bis[i]*delta*(delta - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + 2* Dis[i]**2*bis[i]*delta*(tau - 1)**2*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 2* Dis[i]**2*(tau - 1)**2 + 4* Dis[i]*bis[i]**2*delta*(tau - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 4* Dis[i]*bis[i]*delta*(tau - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 -  Dis[i]*bis[i]*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 4* Dis[i]*bis[i]*(tau - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) -  Dis[i] + 2*bis[i]**3*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**3 + bis[i]**2*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 6*bis[i]**2*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**3 + 2*bis[i]**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - bis[i]*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + 4*bis[i]*delta*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**3 + bis[i]/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 2*bis[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2)*exp(-Cis[i]*(delta - 1)**2 -  Dis[i]*(tau - 1)**2)
    return phir

def iapws95_d3Ar_ddelta2dtau_naive(tau, delta):
    # Not in publication
    phir = 0
    for i in range(7):
        phir += delta**dis[i]*dis[i]*nis[i]*tau**tis[i]*tis[i]*(dis[i] - 1)/(delta**2*tau)
    for i in range(7,51):
        phir += delta**dis[i]*nis[i]*tau**tis[i]*tis[i]*(-2*cis[i]*delta**cis[i]*dis[i] + cis[i]*delta**cis[i]*(cis[i]*delta**cis[i] - cis[i] + 1) + dis[i]*(dis[i] - 1))*exp(-delta**cis[i])/(delta**2*tau)
    for i in range(51, 54):
        phir += delta**dis[i]*nis[i]*tau**tis[i]*(4*alphas[i-51]*betas[i-51]*(gammas[i-51] - tau)*(2*alphas[i-51]*(delta - epsilons[i-51])**2 - 1) - 8*alphas[i-51]*betas[i-51]*dis[i]*(delta - epsilons[i-51])*(gammas[i-51] - tau)/delta + 2*alphas[i-51]*tis[i]*(2*alphas[i-51]*(delta - epsilons[i-51])**2 - 1)/tau - 4*alphas[i-51]*dis[i]*tis[i]*(delta - epsilons[i-51])/(delta*tau) + 2*betas[i-51]*dis[i]*(dis[i] - 1)*(gammas[i-51] - tau)/delta**2 + dis[i]*tis[i]*(dis[i] - 1)/(delta**2*tau))*exp(-alphas[i-51]*(delta - epsilons[i-51])**2 - betas[i-51]*(-gammas[i-51] + tau)**2)
    for i in range(2):
        phir += 2*nis[i+54]*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**bis[i]*(2*Ais[i]*Cis[i]*bis[i]*delta*(delta - 1)*(2*delta - 2.0)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) - Ais[i]*bis[i]*(2*delta - 2.0)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) + 4*Cis[i]*Dis[i]*bis[i]*delta*(delta - 1)*(tau - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 2*Cis[i]*Dis[i]*delta*(tau - 1)*(2*Cis[i]*(delta - 1)**2 - 1) + 4*Cis[i]*Dis[i]*(delta - 1)*(tau - 1) + 4*Cis[i]*bis[i]**2*delta*(delta - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 4*Cis[i]*bis[i]*delta*(delta - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 2*Cis[i]*bis[i]*delta*(2*Cis[i]*(delta - 1)**2 - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + 4*Cis[i]*bis[i]*(delta - 1)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - Dis[i]*bis[i]*delta*(tau - 1)*(2*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + 2*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - 2*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 2*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 4*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2 + bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - 2*Dis[i]*bis[i]*(tau - 1)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - bis[i]**2*delta*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*(2*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + 2*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - 2*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 2*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 4*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2 + bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - 2*bis[i]**2*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 - bis[i]*delta*(Ais[i]*bis[i]*(2*delta - 2.0)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) - Ais[i]*(2*delta - 2.0)*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2*(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)) + Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - Ais[i]*(2*delta - 2.0)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + Ais[i]*(delta - 1)**2*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) - bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) + bis[i]*delta*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*(2*Ais[i]**2*(delta - 1)**2*((delta - 1.0)**2)**(1/betas[i+3])/(betas[i+3]**2*(delta - 1.0)**4) + 2*Ais[i]*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) - 2*Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**3) + 2*Ais[i]*(delta - 1)**2*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]**2*(delta - 1.0)**4) + 4*Bis[i]*ais[i]**2*((delta - 1)**2)**ais[i]/(delta - 1)**2 - 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1)**2 + bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2) - (Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))**2/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2))/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2 + 2*bis[i]*(Ais[i]*(2*delta - 2.0)*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)*((delta - 1.0)**2)**(1/(2*betas[i+3]))/(betas[i+3]*(delta - 1.0)**2) + 2*Bis[i]*ais[i]*((delta - 1)**2)**ais[i]/(delta - 1))*(Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)/(Bis[i]*((delta - 1)**2)**ais[i] + (Ais[i]*((delta - 1.0)**2)**(1/(2*betas[i+3])) - tau + 1)**2)**2)*exp(-Cis[i]*(delta - 1)**2 - Dis[i]*(tau - 1)**2)

    return phir

### Derivatives of Distance Function

def d_psi_d_delta(i, tau, delta):
    '''i is either 0 or 1 for 55 or 56.
    >>> d_psi_d_delta(0, 647.096/647., 358./322.)
    -4.411951793785948
    '''
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    ans = -2*Cis[i]*(delta-1)*psi
    return ans

def d_Delta_d_delta(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d_Delta_d_delta(1, 647.096/647., 358./322.)
    3.595414062719538e-06
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))
    ans = (delta - 1)*(
    Ais[i]*theta*2/betas[i+3]*((delta-1.)**2 )**(1./(2*betas[i+3])-1.)
    + 2*Bis[i]*ais[i]*((delta-1)**2)**(ais[i]-1))
    return ans


def d_Delta_bd_delta(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d_Delta_bd_delta(1, 647.096/647., 358./322.)
    7.931159558108671e-06
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))
    Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

    _d_Delta_d_delta = (delta - 1)*(
    Ais[i]*theta*2/betas[i+3]*((delta-1.)**2 )**(1./(2*betas[i+3])-1.)
    + 2*Bis[i]*ais[i]*((delta-1)**2)**(ais[i]-1))

    ans = bis[i]*Delta**(bis[i]-1)*_d_Delta_d_delta
    return ans


def d2_psi_d2_delta(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_psi_d2_delta(1, 647.096/647., 358./322.)
    -8.581401910121393
    '''
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    ans = (2*Cis[i]*(delta-1)**2 - 1)*2*Cis[i]*psi
    return ans


def d2_Delta_d2_delta(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_Delta_d2_delta(1, 647.096/647., 358./322.)
    0.0002472143243416378
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))

    _d_Delta_d_delta = (delta - 1)*(
    Ais[i]*theta*2/betas[i+3]*((delta-1.)**2 )**(1./(2*betas[i+3])-1.)
    + 2*Bis[i]*ais[i]*((delta-1)**2)**(ais[i]-1))

    first = 1./(delta-1.)*_d_Delta_d_delta
    second = 4*Bis[i]*ais[i]*(ais[i]-1)*((delta-1)**2)**(ais[i]-2)
    third = 2*Ais[i]**2*(1./betas[i+3])**2*(((delta-1)**2)**(1./(2*betas[i+3])-1))**2
    fourth = Ais[i]*theta*4/betas[i+3]*(1./(2*betas[i+3])-1)*((delta-1)**2)**(1./(2*betas[i+3])-2)
    ans = first + (delta-1.)**2*(second + third + fourth)
    return ans


def d2_Delta_bd2_delta(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_Delta_bd2_delta(1, 647.096/647., 358./322.)
    0.0005157293089972383
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))
    Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

    _d_Delta_d_delta = (delta - 1)*(
    Ais[i]*theta*2/betas[i+3]*((delta-1.)**2 )**(1./(2*betas[i+3])-1.)
    + 2*Bis[i]*ais[i]*((delta-1)**2)**(ais[i]-1))

    first = 1./(delta-1.)*_d_Delta_d_delta
    second = 4*Bis[i]*ais[i]*(ais[i]-1)*((delta-1)**2)**(ais[i]-2)
    third = 2*Ais[i]**2*(1./betas[i+3])**2*(((delta-1)**2)**(1./(2*betas[i+3])-1))**2
    fourth = Ais[i]*theta*4/betas[i+3]*(1./(2*betas[i+3])-1)*((delta-1)**2)**(1./(2*betas[i+3])-2)
    _d2_Delta_d2_delta = first + (delta-1.)**2*(second + third + fourth)

    ans = bis[i]*(Delta**(bis[i]-1)*_d2_Delta_d2_delta +
    (bis[i]-1)*Delta**(bis[i]-2)*_d_Delta_d_delta**2)
    return ans


def d_Delta_bd_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d_Delta_bd_tau(1, 647.096/647., 358./322.)
    -0.0002958235123606516
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))
    Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]
    ans = -2*theta*bis[i]*Delta**(bis[i]-1)
    return ans


def d_psi_d_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d_psi_d_tau(1, 647.096/647., 358./322.)
    -0.159135911130251
    '''
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    ans = -2*Dis[i]*(tau-1)*psi
    return ans


def d2_psi_d2_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_psi_d2_tau(1, 647.096/647., 358./322.)
    -1072.4719549824797
    '''
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    ans = (2*Dis[i]*(tau-1)**2 -1)*2*Dis[i]*psi
    return ans


def d2_Delta_bd2_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_Delta_bd2_tau(1, 647.096/647., 358./322.)
    4.370635612303636
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1/(2*betas[i+3]))
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

    ans = 2*bis[i]*Delta**(bis[i]-1) + 4*theta**2*bis[i]*(bis[i]-1)*Delta**(bis[i]-2)
    return ans


def d2_psi_d_delta_d_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_psi_d_delta_d_tau(1, 647.096/647., 358./322.)
    1.1386619231183175
    '''
    psi = exp(-Cis[i]*(delta-1)**2 - Dis[i]*(tau-1)**2)
    ans = 4*Cis[i]*Dis[i]*(delta-1)*(tau-1)*psi
    return ans


def d2_Delta_bd_delta_d_tau(i, tau, delta):
    ''''i is either 0 or 1 for 55 or 56.
    >>> d2_Delta_bd_delta_d_tau(1, 647.096/647., 358./322.)
    -0.027232925382835605
    '''
    theta = (1-tau) + Ais[i]*((delta-1)**2)**(1./(2*betas[i+3]))
    Delta = theta**2 + Bis[i]*((delta-1)**2)**ais[i]

    _d_Delta_d_delta = (delta - 1)*(
    Ais[i]*theta*2/betas[i+3]*((delta-1.)**2 )**(1./(2*betas[i+3])-1.)
    + 2*Bis[i]*ais[i]*((delta-1)**2)**(ais[i]-1))

    first = -Ais[i]*bis[i]*2/betas[i+3]*Delta**(bis[i]-1)*(delta-1)*((delta-1)**2)**(1/(2*betas[i+3])-1)
    second = -2*theta*bis[i]*(bis[i]-1)*Delta**(bis[i]-2)*_d_Delta_d_delta
    ans = first + second
    return ans

def test_iapws95_d2A_d2deltar():
    assert_close(iapws95_d2Ar_ddelta2(3.23548, 2.652088725981779), -681188.0609390885, rtol=1e-11)

    # Test equation for delta = 1 exactly, derived with sympy limit(thing, delta, 1)
    for tau in linspace(.01, .8, 100):
        assert_close(iapws95_d2Ar_ddelta2(tau,.99999999999), iapws95_d2Ar_ddelta2(tau, 1), rtol=2e-8)
    for tau in linspace(10, 1.5, 100):
        assert_close(iapws95_d2Ar_ddelta2(tau,.99999999999), iapws95_d2Ar_ddelta2(tau, 1), rtol=2e-8)


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d2A_d2deltar_vs_naive(precise=False, allow_fail=True):
    '''Overall performs very well. 2e-10 was needed in 2000^2 points for like 1 point.
    Smaller number of points work to 1e-12. Having an absolute tolerance of 1e-15
    would also work find.
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 500
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        print(i)
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d2Ar_ddelta2(tau, delta)
            val_naive = float(ddAdddelta_res(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=2e-10)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()

#test_iapws95_d2A_d2deltar_vs_naive()
#test_iapws95_d2A_d2deltar_vs_naive(precise=True, allow_fail=False)

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d3A_d3deltar_vs_naive(precise=False, allow_fail=True):
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 100
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
#        print(i)
#        rhos = logspace(log10(1e-10), log10(1e-4), N)
        rhos = logspace(log10(1e-4), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d3Ar_ddelta3(tau, delta)
            val_naive = float(dddA_ddddelta_res(mpf(tau), mpf(delta)))
            if allow_fail:
                try:
                    assert_close(val, val_naive, rtol=1e-10)
                except:
#                    print([T, rho, abs(1.0 - val/val_naive)])
                    print([tau, delta, val, T, rho])
            rerri = abs(1.0 - val/val_naive)
#            if rerri > 1e-4:
#                print([T, rho, rerri])
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()

#test_iapws95_d3A_d3deltar_vs_naive(precise=True, allow_fail=True)

#test_iapws95_d3A_d3deltar_vs_naive(precise=False, allow_fail=True)



def test_iapws95_d3Ar_ddelta3():
    # Extremely flat function. Requires huge steps to get results.

    tau, delta = .7, 2.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 1e-3, 2.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 1e-8, 2.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 15.5, 2.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 15.5, 11.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 2.5, 1.2
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-11)

    tau, delta = 1e-4, 1e-3
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*1e-2, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-10)

    tau, delta = 0.32916247703300405, 1.7314054528190082
    assert_close(derivative(lambda delta: iapws95_d2Ar_ddelta2(tau, delta), delta, dx=delta*5e-4, order=7),
                 iapws95_d3Ar_ddelta3(tau, delta), rtol=1e-8)


    assert_close(iapws95_d3Ar_ddelta3(.5, .6), 0.06421246015370234, rtol=1e-9)

    tau, delta, T, rho = 0.5028042728122091, 3.3706405889480084e-09, 1286.9739478957888, 1.0853462696412588e-06
    num = derivative(lambda d: iapws95_d2Ar_ddelta2(tau, d), delta, dx=delta*.3, order=7)
    assert_close(iapws95_d3Ar_ddelta3(tau, delta), num, rtol=1e-7)

@pytest.mark.slow
@pytest.mark.mpmath
def test_iapws95_d3Ar_ddelta3_mpmath():
    import mpmath as mp
    mp.mp.dps = 100
    tau, delta, T, rho = 0.5028042728122091, 3.3706405889480084e-09, 1286.9739478957888, 1.0853462696412588e-06
    assert_close(iapws95_d3Ar_ddelta3(tau, delta),
                 float(dddA_ddddelta_res(mp.mpf(tau), mp.mpf(delta))))

    tau, delta, T, rho = 0.6879013719642095, 4.2599878142989604e-13, 940.6813627254511, 1.3717160762042654e-10
    assert_close(iapws95_d3Ar_ddelta3(tau, delta),
                 float(dddA_ddddelta_res(mp.mpf(tau), mp.mpf(delta))), rtol=1e-3)
#
#
#0.5809660021590505 3.308239737943376e-13 1113.8276553106205 1.065253195617767e-10 0.668894347379583 0.6689624104241152 0.00010174419888420161
#0.5809660021590505 3.5241129527138654e-13 1113.8276553106205 1.1347643707738647e-10 0.6688943474959489 0.6689624104240671 0.00010174402486240464
#0.5809660021590505 3.7540725845964094e-13 1113.8276553106205 1.208811372240044e-10 0.668894347495896 0.6689624104240159 0.00010174402486495815

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d3Ar_ddeltadtau2_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 100
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
#        print(i)
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d3Ar_ddeltadtau2(tau, delta)
            val_naive = float(iapws95_d3Ar_ddeltadtau2_naive(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=2e-10)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()
#test_iapws95_d3Ar_ddeltadtau2_vs_naive(precise=True, allow_fail=False)
#test_iapws95_d3Ar_ddeltadtau2_vs_naive(precise=False, allow_fail=False)



@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d3Ar_ddelta2dtau_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 100
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        print(i)
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d3Ar_ddelta2dtau(tau, delta)
            val_naive = float(iapws95_d3Ar_ddelta2dtau_naive(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=2e-10)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()
#test_iapws95_d3Ar_ddelta2dtau_vs_naive(precise=True, allow_fail=False)
#test_iapws95_d3Ar_ddelta2dtau_vs_naive(precise=False, allow_fail=True)



def test_iapws95_dA_ddeltar():
    rhoc_inv = (1.0/322.0)
    # Has the highest error found so far. No idea where it comes from.
    T, rho = 221.12845138055224, 1032.1609281563633
    tau = 647.096/T
    delta = rho*rhoc_inv
    assert_close(iapws95_dAr_ddelta(tau, delta),
                 dAddelta_res(tau, delta), rtol=5e-9)

    assert_close(iapws95_dAr_ddelta(.999999999999999, .999999999999999),
                 iapws95_dAr_ddelta(1,1), rtol=1e-13)

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_dA_ddeltar_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 500
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_dAr_ddelta(tau, delta)
            val_naive = float(dAddelta_res(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=5e-9)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()

#test_iapws95_dA_ddeltar_vs_naive()
#test_iapws95_dA_ddeltar_vs_naive(precise=True, allow_fail=False)


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_Ar_vs_naive(precise=False, allow_fail=True):
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 300
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_Ar(tau, delta)
            val_naive = float(calcA_res(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=2e-9)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()
#test_iapws95_Ar_vs_naive()
#test_iapws95_Ar_vs_naive(precise=True, allow_fail=False)


def test_iapws95_Ar():
    assert_close(iapws95_Ar(647.096/300.0, 999.0/322), -9.57577716026768, rtol=1e-11)

    # Point was being issue in third delta derivative
    assert_close(iapws95_Ar(0.827680930232558, 4.443026216725733e-07), -3.726849328574788e-7, rtol=1e-13)




def test_iapws95_dAr_dtau():
    assert_close(iapws95_dAr_dtau(647.096/300.0, 999.0/322),
                 -7.704333630957023, rtol=1e-11)

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_dAr_dtau_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 300
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_dAr_dtau(tau, delta)
            val_naive = float(dAdtau_res(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=1e-8)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()
#test_iapws95_dAr_dtau_vs_naive()




@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d2Ar_dtau2_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x

    errs = []
    rerr = 0
    N = 500
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d2Ar_dtau2(tau, delta)
            val_naive = float(ddAddtau_res(mpf(tau), mpf(delta)))
            if allow_fail:
                assert_close(val, val_naive, rtol=5e-10)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()
#test_iapws95_d2Ar_dtau2_vs_naive(precise=True, allow_fail=False)


def test_iapws95_d2Ar_dtau2():
    assert_close(iapws95_d2Ar_dtau2(647.096/300.0, 999.0/322),
                 -1.2616419775539731, rtol=1e-11)


@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_d2Ar_ddeltadtau_vs_naive(precise=False, allow_fail=True):
    '''
    '''
    if precise:
        mp = make_me_precise()
        mpf = mp.mpf
    else:
        mpf = lambda x: x
    errs = []
    rerr = 0
    N = 300
    Ts = linspace(200.0, 5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
#        print(i)
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d2Ar_ddeltadtau(tau, delta)
            val_naive = float(dAddeltatau_res(mpf(tau), mpf(delta)))
            assert_close(val.real, val_naive, rtol=1e-8)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
#    print(rerr/N**2, np.std(errs), np.max(errs))
    make_me_float()

#test_iapws95_d2Ar_ddeltadtau_vs_naive()
#test_iapws95_d2Ar_ddeltadtau_vs_naive(precise=True, allow_fail=False)

def test_iapws95_d2Ar_ddeltadtau():
    val = iapws95_d2Ar_ddeltadtau(647.096/300.0, 999.0/322)
    assert_close(val, -0.1984035623854279, rtol=1e-12)

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_A0_vs_naive():
    '''Can't think of a way to optimize this.
    '''
    errs = []
    rerr = 0
    N = 400
    Ts = linspace(200.0,  5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_A0(tau, delta)
            val_naive = iapws95_A0(tau, delta)
            assert_close(val, val_naive, rtol=1e-15)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
#    print(rerr/N**2, np.std(errs), np.max(errs))
#test_iapws95_A0_vs_naive()

@pytest.mark.slow
@pytest.mark.fuzz
def test_iapws95_iapws95_dA0_dtau_vs_naive():
    '''Can't think of a way to optimize this.
    '''
    errs = []
    rerr = 0
    N = 300
    Ts = linspace(200.0,  5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_dA0_dtau(tau, delta)
            val_naive = dAdtau_idg(tau, delta)
            assert_close(val, val_naive, rtol=1e-15)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
#    print(rerr/N**2, np.std(errs), np.max(errs))

@pytest.mark.slow
@pytest.mark.fuzz
def test_ddAddtau_idg_vs_naive():
    '''Can't think of a way to optimize this.
    '''
    errs = []
    rerr = 0
    N = 100
    Ts = linspace(200.0,  5000.0, N)
    rhoc_inv = (1.0/322.0)
    for i, T in enumerate(Ts):
        rhos = logspace(log10(1e-10), log10(5000), N)
        for rho in rhos:
            tau = 647.096/T
            delta = rho*rhoc_inv
            val = iapws95_d2A0_dtau2(tau, delta)
            val_naive = ddAddtau_idg(tau, delta)
            assert_close(val, val_naive, rtol=5e-15)
            rerri = abs(1.0 - val/val_naive)
            rerr += rerri
            errs.append(rerri)
#    print(rerr/N**2, np.std(errs), np.max(errs))
#test_ddAddtau_idg_vs_naive()


def test_iapws95_A0():
    ans = iapws95_A0(.5345, .575745)
    assert_close(ans, -7.3790791583143, rtol=1e-14)

def test_iapws95_dA0_dtau():
    ans = iapws95_dA0_dtau(.5345, .575745)
    assert_close(ans, 13.16120203177092, rtol=1e-14)

def test_iapws95_d2A0_dtau2():
    ans = iapws95_d2A0_dtau2(.5345, .575745)
    assert_close(ans, -14.97966801918871, rtol=1e-14)
    assert_close(ans, derivative(iapws95_dA0_dtau, .5345, args=(.575745,), dx=1e-8))

def test_iapws95_d3A0_dtau3():
    ans = iapws95_d3A0_dtau3(.5345, .575745)
    assert_close(ans, 67.50267480495313, rtol=1e-14)
    assert_close(ans, derivative(iapws95_d2A0_dtau2, .5345, args=(.575745,), dx=1e-8))

@pytest.mark.slow
@pytest.mark.CoolProp
@pytest.mark.fuzz
def test_rho_iapws95_CoolProp():
    from CoolProp.CoolProp import PropsSI
    N = 40
    Ts = linspace(273.16+1e-10,  1073.15-1e-10, N)
    Ps = logspace(log10(1e-3), log10(100e6), N)

    for T in Ts:
        for P in Ps:
            rho_implemented = iapws95_rho(T, P)
            rho_CoolProp = PropsSI('DMASS', 'T', T, 'P', P, 'water')
            assert_close(rho_implemented, rho_CoolProp, rtol=1e-10)
            # some points found to fail near Psat curve as expected.


def test_iapws97_rho_extrapolated():
    region5_highT = iapws97_rho_extrapolated(2300, 20e6)
    region5_highT_num = (iapws97_region5_rho(2273.15, 20e6)
                         + (2300-2273.15)*derivative(iapws97_region5_rho, 2273.15, args=(20e6,), dx=.1, order=3))
    assert_close(region5_highT, region5_highT_num, rtol=1e-10)

    region2_highT = iapws97_rho_extrapolated(1100, 80e6)
    region2_highT_num = (iapws97_region2_rho(1073.15, 80e6)
                         + (1100-1073.15)*derivative(iapws97_region2_rho, 1073.15, args=(80e6,), dx=.1, order=5))

    assert_close(region2_highT, region2_highT_num, rtol=1e-9)


    region1_lowT = iapws97_rho_extrapolated(200, 80e6)
    assert region1_lowT == iapws97_rho(273.15, 80e6)

def test_iapws95_P():
    assert_close(iapws95_rho(300.0, iapws95_P(300, 1000)), 1000)

def test_iapws95_T_err():
    err, derr = iapws95_T_err(300, 1000, 1e5)
    assert_close(err, 7733001.355973767, rtol=1e-11)
    assert_close(derr, 639359.0465881994, rtol=1e-11)

    assert_close(derivative(lambda T: iapws95_T_err(T, 1000, 1e5)[0], 300, dx=1e-4),
                 iapws95_T_err(300, 1000, 1e5)[1])

def test_iapws95_rho():
    '''TODO points:

    iapws95_rho(200.0, 1e9) - not solving
    '''
    assert_close(iapws95_rho(273.1600000001, 0.001), 7.932210036861784e-09, rtol=1e-8)

    assert_close(iapws95_rho(350.0, 1e6), 974.1288271329855, rtol=1e-8)
    assert_close(iapws95_rho(981.3822764554016, 171493178.34983346), 444.5570512999293)

    # Point where was starting from negative density initially.
    assert_close(iapws95_rho(2357., 97719212), 85.77393882818544, rtol=1e-9)


    # Three points CoolProp is finding the vapor root when the liquid one is NOT correct
    # initially was making the wrong call there
    assert_close(iapws95_rho(432.0135947190398, 600559.0434678708), 3.1715230968689263)

    assert_close(iapws95_rho(443.36028005610694, 796123.0461361709), 4.141564959829041)

    assert_close(iapws95_rho(485.9103500701087, 2014934.1250668736), 10.114793546282295)

    assert_close(iapws95_rho(472.89458917842654, 1546542.3293244045), 7.819904266670308)

    # Found comparing against coolprop
    assert_close(iapws95_rho(640, 20265239.648236595), 481.5275168680331, rtol=1e-10)

    # Slightly different density than CoolProp here
    assert_close(iapws95_rho(647.08, 22059526.03804436), 295.66686689744756, rtol=1e-10)


@pytest.mark.slow
@pytest.mark.CoolProp
def test_iapws95_rho_vs_Coolprop():
    from CoolProp.CoolProp import PropsSI
    assert_close(iapws95_rho(2357., 97719212),
                 PropsSI('DMASS', 'T', 2357, 'P', 97719212.0, 'water'),
                 rtol=1e-9)

def test_iapws95_T():
    assert_close(iapws95_T(P=20265239.648236595, rho=481.5275168680331), 640.0)

    # Point where converging to wrong solution
    assert_close(iapws95_T(P=1000000000, rho=669.0726669889119), 1749.5356805149595, rtol=1e-7)

    # Point where vapor pressure was failing to calculate
    rho = iapws95_rho(T=1000.0, P=10)
    assert_close(iapws95_T(P=10.0, rho=rho), 1000, rtol=1e-11)


def test_rhol_sat_IAPWS():
    assert_close(iapws92_rhol_sat(273.16), 999.7891346511478, rtol=1e-13)
    assert_close(iapws92_rhol_sat(373.1243), 958.3652337840979, rtol=1e-13)
    assert_close(iapws92_rhol_sat(647.096), 322.0, rtol=1e-13)

def test_rhog_sat_IAPWS():
    assert_close(iapws92_rhog_sat(373.1243), 0.5975863185878799, rtol=1e-13)
    assert_close(iapws92_rhog_sat(273.16), 0.004854262596261411, rtol=1e-13)
    assert_close(iapws92_rhog_sat(647.096), 322)


def test_iapws95_saturation():
    from chemicals.iapws import iapws95_sat_err_and_jac

    err, jac = iapws95_sat_err_and_jac([1.807830655839175e-05, 0.7040053821406961], 300.0)
    assert_close1d(err, [2.660272002685815e-10, 1.681859203017666e-05], rtol=1e-6)
    assert_close2d(jac, [[-2219129832.1579313, 3530.261802535087],
                         [-122750979191014.77, 5014.538087479508]], rtol=1e-7)

    vals = iapws95_saturation(300.0, xtol=1e-5)
    assert_close1d(vals, (3536.806752287503, 996.5130274681349, 0.025589673682920273))

def test_iapws95_dPsat_dT():
    assert_close(derivative(iapws95_Psat, 330.0, dx=330*1e-6),
                 iapws95_dPsat_dT(330)[0], rtol=1e-7)

    dPsat_dT, Psat = iapws95_dPsat_dT(500.0)
    assert_close(dPsat_dT, 49008.17866580053, rtol=1e-10)
    assert_close(Psat, 2639195.8717618496, rtol=1e-13)

    # Check for functional equivalence
    for T in linspace(235.0, 647.096):
        assert_close(iapws95_Psat(T), iapws95_dPsat_dT(T)[1], rtol=1e-15)

def test_iapws95_Tsat():
    # Tested with a LOT of points
    for T in linspace(235.0, 647.096, 100):
        assert_close(iapws95_Tsat(iapws95_Psat(T)), T, rtol=2e-14)

def test_iapws95_Psat():
    assert_close(iapws95_Psat(300.0), 3536.806752274638, rtol=1e-12)
    assert_close(iapws95_Psat(260.0), 222.5574677094123, rtol=1e-12)
    assert_close(iapws95_Psat(500.0), 2639195.8717618496, rtol=1e-12)
    assert_close(iapws95_Psat(620.0), 15900579.384968637, rtol=1e-12)
    assert_close(iapws95_Psat(640.0), 20265209.268028494, rtol=1e-12)
    assert_close(iapws95_Psat(645.0), 21515208.664142672, rtol=1e-12)
    assert_close(iapws95_Psat(647.0), 22038405.72692307, rtol=1e-10)
    assert iapws95_Psat(647.096) == 22064000.0 # Should be dead on

    with pytest.raises(ValueError):
        iapws95_Psat(150.0)

def test_iapws95_rhol_sat():
    assert_close(iapws95_rhol_sat(250.0), 991.1730911906317, rtol=1e-13)
    assert_close(iapws95_rhol_sat(300.0), 996.5130274681278, rtol=1e-13)
    assert_close(iapws95_rhol_sat(500.0), 831.3134495854619, rtol=1e-13)
    assert_close(iapws95_rhol_sat(620.0), 586.8776188274974, rtol=1e-13)
    assert_close(iapws95_rhol_sat(640.0), 481.52614604429033, rtol=1e-13)
    assert_close(iapws95_rhol_sat(645.0), 425.04824669763036, rtol=1e-13)
    assert_close(iapws95_rhol_sat(647.0), 357.34089197171636, rtol=1e-13)
    assert_close(iapws95_rhol_sat(647.08), 340.387972615003, rtol=1e-12)
    assert_close(iapws95_rhol_sat(647.094), 329.1865391205182, rtol=1e-12)
    assert_close(iapws95_rhol_sat(647.0955), 325.7094848352297, rtol=1e-11)
    assert_close(iapws95_rhol_sat(647.09597), 322.9332279314927, rtol=1e-11)
    assert_close(iapws95_rhol_sat(647.095998), 322.2429814792015, rtol=1e-11)
    assert_close(iapws95_rhol_sat(647.09599997), 322.02998574491903, rtol=1e-11)

    # Linear interp - may get replaced in future
    assert_close(iapws95_rhol_sat(647.0959999999), 322.00017602415505, rtol=1e-7)
    assert 322.0 == iapws95_rhol_sat(647.096)

    with pytest.raises(ValueError):
        iapws95_rhol_sat(200.0)

def test_iapws95_rhol_sat_dT():
    for T in linspace(235.0, 647.096, 100):
        assert_close(iapws95_drhol_sat_dT(T)[1], iapws95_rhol_sat(T), rtol=1e-13)

    assert_close1d(iapws95_drhol_sat_dT(277), (0.002398816135429972, 999.9249513005213), rtol=1e-13)
    assert_close1d(iapws95_drhol_sat_dT(647.09599999999), (-1759460.0473706475, 322.00001760241554), rtol=1e-13)

    # Numerical dereivatives do fail near T = 277 K
    for T in linspace(235.0+10e-4, 647.095-10e-4, 10):
        assert_close(iapws95_drhol_sat_dT(T)[0],
                     derivative(iapws95_rhol_sat, T, dx=T*4e-7, order=7), rtol=1e-4)


def test_rhog_sat_IAPWS95():
    assert_close(iapws95_rhog_sat(260.0), 0.0018552889771409127, rtol=1e-13)
    assert_close(iapws95_rhog_sat(400.0), 1.3694075410068125, rtol=1e-13)
    assert_close(iapws95_rhog_sat(600.0), 72.84231718283309, rtol=1e-13)
    assert_close(iapws95_rhog_sat(630.0), 132.8395560369342, rtol=1e-13)
    assert_close(iapws95_rhog_sat(645.0), 224.4505402883077, rtol=1e-13)
    assert_close(iapws95_rhog_sat(647.0), 286.5083958147434, rtol=1e-13)
    assert_close(iapws95_rhog_sat(647.08), 303.4596095608764, rtol=1e-13)
    assert_close(iapws95_rhog_sat(647.094), 314.7568132641434, rtol=1e-13)
    assert_close(iapws95_rhog_sat(647.0958), 319.61981063229433, rtol=2e-12)
    assert_close(iapws95_rhog_sat(647.09598), 321.2362413951111, rtol=1e-10)
    assert_close(iapws95_rhog_sat(647.095998), 321.7569904453351, rtol=1e-10)
    assert_close(iapws95_rhog_sat(647.0959998), 321.92296716781414, rtol=1e-10)
    assert_close(iapws95_rhog_sat(647.09599998), 321.9758031526337, rtol=1e-10)

    assert 322 == iapws95_rhog_sat(647.0959999999999999)
    assert 322 == iapws95_rhog_sat(647.0959999999999)
    assert_close(iapws95_rhog_sat(647.09599999999), 321.9999832063918, rtol=1e-13)

    with pytest.raises(ValueError):
        iapws95_rhog_sat(200.0)

@pytest.mark.slow
@pytest.mark.mpmath
@pytest.mark.fuzz
def test_iapws95_saturation_fits():
    import mpmath as mp
    mp.mp.dps = 50
    iapws.use_mpmath_backend()
    N = 100 # should be able to set arbitrarily high, tested to 1000
    Ts = linspace(273.15, 647.09, N)

    for T in Ts:
        rhol = iapws95_rhol_sat(T)
        P_corr = float(iapws95_Psat(T))
        Psat_mp, rhol_mp, rhog_mp = iapws95_saturation(mp.mpf(T), xtol=1e-20)
            # Everything except > 646.7 is under 2E-13
        assert_close(P_corr, float(Psat_mp), rtol=7.5e-13)
        # Almost everything is under 2e-14
        assert_close(rhol, float(rhol_mp), rtol=2e-13)


    # Test the low temperature regime with a lower tolerance - fitting issues for
    # density
    Ts = linspace(235, 273.15, N)
    for T in Ts:
        rhol = iapws95_rhol_sat(T)
        P_corr = float(iapws95_Psat(T))
        Psat_mp, rhol_mp, rhog_mp = iapws95_saturation(mp.mpf(T), xtol=1e-20)
        assert_close(P_corr, float(Psat_mp), rtol=7.5e-13)
        assert_close(rhol, float(rhol_mp), rtol=3e-11)

    iapws.reset_backend()


def test_rhog_sat_IAPWS95_vs_saturation():
    # Specific points
    Ts = [260.0, 400.0, 600.0, 630.0, 645]
    for T in Ts:
        assert_close(iapws95_saturation(T)[2],
                     iapws95_rhog_sat(T), rtol=1e-12)

    # 647 requires mpmath
@pytest.mark.slow
@pytest.mark.fuzz
def test_rhog_sat_IAPWS95_vs_saturation2():
    Ts = [260.0, 400.0, 600.0, 630.0, 645]
    import mpmath as mp
    mp.mp.dps = 50
    iapws.use_mpmath_backend()
    Ts = [647, 647.08, 647.094, 647.0959]
    for T in Ts:
        Psat_mp, rhol_mp, rhog_mp = iapws95_saturation(mp.mpf(T), xtol=1e-20)
        assert_close(float(rhog_mp), float(iapws95_rhog_sat(T)), rtol=1e-13)

    Ts = [647.0958, 647.09598, 647.095998, 647.0959998, 647.09599998]
    # Fit lower accuracy points
    for T in Ts:
        Psat_mp, rhol_mp, rhog_mp = iapws95_saturation(mp.mpf(T), xtol=1e-20)
        assert_close(float(rhog_mp), float(iapws95_rhog_sat(T)), rtol=1e-10)
    iapws.reset_backend()

@pytest.mark.slow
@pytest.mark.CoolProp
def test_rhog_sat_IAPWS95_CoolProp():
    from CoolProp.CoolProp import PropsSI
    Ts = [400.0, 600.0]
    for T in Ts:
        assert_close(iapws95_rhog_sat(T),
                     PropsSI('DMASS', 'T', T, 'Q', 1, 'water'), rtol=1e-13)

    Ts = [630.0, 645, 647]
    for T in Ts:
        assert_close(iapws95_rhog_sat(T),
                     PropsSI('DMASS', 'T', T, 'Q', 1, 'water'), rtol=1e-9)

    Ts = [647.08, 647.094,
#          647.0958
#    647.09599
          ]
    for T in Ts:
        assert_close(iapws95_rhog_sat(T),
                     PropsSI('DMASS', 'T', T, 'Q', 1, 'water'), rtol=1e-4)

#
#test_rhog_sat_IAPWS95()
#test_rhog_sat_IAPWS95_CoolProp()
#test_rhog_sat_IAPWS95_vs_saturation()

def test_iapws95_d4Ar_ddelta2dtau2():
    assert_close(iapws95_d4Ar_ddelta2dtau2(647.096/300.0, 999.0/322), -2.6564229154801002)
    assert_close(iapws95_d4Ar_ddelta2dtau2(.7, 1.2), 1.28154351717541)


def test_iapws95_A0_tau_derivatives():
    tau, delta = 647.096/300.0, 999.0/322
    together = iapws95_A0_tau_derivatives(tau, delta)
    A0 = iapws95_A0(tau, delta)
    dA0_dtau = iapws95_dA0_dtau(tau, delta)
    d2A0_dtau2 = iapws95_d2A0_dtau2(tau, delta)
    d3A0_dtau3 = iapws95_d3A0_dtau3(tau, delta)
    assert_close1d(together, (A0, dA0_dtau, d2A0_dtau2, d3A0_dtau3), rtol=1e-15)

def test_consistency_iapws95_rho_iapws95_P():
    P = 39.0693
    rho = iapws95_rho(T=235.0, P=P)
    P_check = iapws95_P(T=235.0, rho=rho)
    assert abs(1-P/P_check) < 5e-6

def test_iapws95_properties():
    expect = [996.5563403888951, 112553.33413264707, 393.06243381456477, 112653.67968858521, 4130.178615033825, 4180.639522022912, 1501.520415056628, -2.2023653545981183e-07, 0.0009207295643366906, 1.978788044482276e-08, 4.4896388297803826e-07]
    assert_close1d(iapws95_properties(T=300.0, P=1e5), expect, rtol=1e-13)


def test_iapws92_Psat():
    assert_close(iapws92_Psat(400.0), 245765.263541822, rtol=1e-13)

def test_iapws92_dPsat_dT():
    assert_close(iapws92_dPsat_dT(400.0)[0], 7483.4709410560408287, rtol=1e-12)
    assert_close(iapws92_dPsat_dT(400.0)[1], 245765.263541822, rtol=1e-12)
