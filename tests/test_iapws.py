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


def test_iapws97_dG_dpi_region1():
    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 0.12923271825448354, rtol=1e-14)


def test_iapws97_dG_dpi_region2():
    assert_close(iapws97_dGr_dpi_region2(.656, 16), -0.006292631931275252, rtol=1e-14)

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


@pytest.mark.CoolProp
@pytest.mark.slow
@pytest.mark.xfail
def test_iapws97_region_2_rho_coolprop():
    # Need to fix - either find where error is and fix there, or stop using horner's method
    # gas
    from CoolProp.CoolProp import PropsSI
    P_lim = 1e-6
    Ts = linspace(273.15+1e-10, 1073.15-1e-10, 100)
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
            try:
                assert_close(rho_CoolProp, rho_implemented, rtol=2e-15)
            except:
                print([T, P, 1-rho_implemented/rho_CoolProp])
test_iapws97_region_2_rho_coolprop()