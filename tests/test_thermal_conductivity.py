# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
SOFTWARE.'''

from numpy.testing import assert_allclose
import pytest
from fluids.numerics import assert_close, assert_close1d
from chemicals.thermal_conductivity import *
from thermo.identifiers import checkCAS
from chemicals.thermal_conductivity import k_data_Perrys_8E_2_314, k_data_Perrys_8E_2_315, k_data_VDI_PPDS_10, k_data_VDI_PPDS_9

def test_Perrys2_314_data():
    # In perry's, only 102 is used. No chemicals are missing.
    # Tmaxs all match to 5E-4. Tmins match to 1E-3.
    assert all([checkCAS(i) for i in k_data_Perrys_8E_2_314.index])
    tots_calc = [k_data_Perrys_8E_2_314[i].abs().sum() for i in [u'C1', u'C2', u'C3', u'C4', u'Tmin', u'Tmax']]
    tots = [48935634.823768869, 297.41545078799999, 421906466448.71423, 232863514627157.62, 125020.26000000001, 347743.42000000004]
    assert_close(tots_calc, tots)
    
    assert k_data_Perrys_8E_2_314.index.is_unique
    assert k_data_Perrys_8E_2_314.shape == (345, 7)


def test_Perrys2_315_data():
    # From perry's - Deuterium  ,  Nitrogen trifluoride  ,  Nitrous oxide   Silicon tetrafluoride  ,  Terephthalic acid  all have no data
    # All perry's use #100.
    # Tmins all match at 5E-4.
    # Tmaxs all match at 2E-3.
    assert all([checkCAS(i) for i in k_data_Perrys_8E_2_315.index])
    tots_calc = [k_data_Perrys_8E_2_315[i].abs().sum() for i in [u'C1', u'C2', u'C3', u'C4', u'C5', u'Tmin', u'Tmax']]
    tots = [82.001667499999996, 0.19894598900000002, 0.0065330144999999999, 0.00046928630199999995, 1.0268010799999999e-07, 70996.369999999995, 138833.41]
    assert_close1d(tots_calc, tots)
    
    assert k_data_Perrys_8E_2_315.index.is_unique
    assert k_data_Perrys_8E_2_315.shape == (340, 8)


def test_VDI_PPDS_10_data():
    '''Average deviation of 2.4% from tabulated values. Many chemicals have
    much higher deviations. 10% or more deviations:
    ['75-34-3', '107-06-2', '106-93-4', '420-46-2', '71-55-6', '79-34-5', 
    '67-72-1', '76-12-0', '76-13-1', '76-14-2', '540-54-5', '75-01-4', 
    '75-35-4', '79-01-6', '127-18-4', '462-06-6', '108-90-7', '108-86-1', 
    '108-41-8', '100-44-7', '108-93-0', '100-61-8', '121-69-7', '91-66-7']
    
    These have been checked - it appears the tabulated data is just incorrect.
    '''

    assert all([checkCAS(i) for i in k_data_VDI_PPDS_10.index])
    tots_calc = [k_data_VDI_PPDS_10[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E']]
    tots = [2.2974640014599998, 0.015556001460000001, 1.9897655000000001e-05, 6.7747269999999993e-09, 2.3260109999999999e-12]
    assert_close1d(tots_calc, tots)
    
    assert k_data_VDI_PPDS_10.index.is_unique
    assert k_data_VDI_PPDS_10.shape == (275, 6)


def test_VDI_PPDS_9_data():
    '''Average deviation of 0.71% from tabulated values. The following have 
    larger deviations
        
    ['124-18-5', '629-59-4', '629-78-7', '526-73-8', '95-63-6']
    
    These have been checked - it appears the tabulated data is just incorrect.
    '''

    assert all([checkCAS(i) for i in k_data_VDI_PPDS_9.index])
    tots_calc = [k_data_VDI_PPDS_9[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E']]
    tots = [63.458699999999993, 0.14461469999999998, 0.00042270770000000005, 1.7062660000000002e-06, 3.2715370000000003e-09]
    assert_close1d(tots_calc, tots)
    
    assert k_data_VDI_PPDS_9.index.is_unique
    assert k_data_VDI_PPDS_9.shape == (271, 6)


def test_CSP_liq():
    kl = Sheffy_Johnson(300, 47, 280)
    assert_close(kl, 0.17740150413112196)

    kl = Sato_Riedel(300, 47, 390, 520)
    assert_close(kl, 0.2103769246133769)

    kl = Lakshmi_Prasad(273.15, 100)
    assert_close(kl, 0.013664450000000009)

    kl = Gharagheizi_liquid(300, 40, 350, 1E6, 0.27)
    assert_close(kl, 0.2171113029534838)

    kl = Nicola_original(300, 142.3, 611.7, 0.49, 201853)
    assert_close(kl, 0.2305018632230984)

    kl = Nicola(300, 142.3, 611.7, 2110000.0, 0.49)
    assert_close(kl, 0.10863821554584034)
    # Not at all sure about this one

    kl = Bahadori_liquid(273.15, 170)
    assert_close(kl, 0.14274278108272603)
    
    kl = Mersmann_Kind_thermal_conductivity_liquid(400, 170.33484, 658.0, 0.000754, {'C': 12, 'H': 26})
    assert_close(kl, 0.0895271829899285)


def test_CSP_liq_dense():
    # From [2]_, for butyl acetate.
    kl_dense = DIPPR9G(515.05, 3.92E7, 579.15, 3.212E6, 7.085E-2)
    assert_close(kl_dense, 0.0864419738671184)

    kld1 = Missenard(304., 6330E5, 591.8, 41E5, 0.129)
    assert_close(kld1, 0.21983757770696569)
    # # butyl acetate
    kld2 = Missenard(515.05, 3.92E7, 579.15, 3.212E6, 7.085E-2)
    assert_close(kld2, 0.086362465280714396)


def test_CSP_gas():
    # 2-methylbutane at low pressure, 373.15 K. Mathes calculation in [1]_.
    kg =  Eucken(72.151, 135.9, 8.77E-6)
    assert_close(kg, 0.018792644287722975)

    # 2-methylbutane at low pressure, 373.15 K. Mathes calculation in [1]_.
    kg = Eucken_modified(72.151, 135.9, 8.77E-6)
    assert_close(kg, 0.023593536999201956)

    # CO, brute force tests on three  options for chemtype
    kg1 = DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear')
    assert kg1 == DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92) # No argument
    kg2 = DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='monoatomic')
    kg3 = DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='nonlinear')
    assert_allclose([kg1, kg2, kg3], [0.01813208676438415, 0.023736881470903245, 0.018625352738307743])

    with pytest.raises(ValueError):
        DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='FAIL')


    kg = Chung(T=373.15, MW=72.151, Tc=460.4, omega=0.227, Cvm=135.9, mu=8.77E-6)
    assert_close(kg, 0.023015653729496946)

    kg = Eli_Hanley(T=373.15, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267, omega=0.227, Cvm=135.9)
    assert_close(kg, 0.022479517891353377)

    kg = Eli_Hanley(T=1000, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267, omega=0.227, Cvm=135.9)
    assert_close(kg, 0.06369581356766069)

    kg = Bahadori_gas(40+273.15, 20)
    assert_close(kg, 0.031968165337873326)

    kg = Gharagheizi_gas(580., 16.04246, 111.66, 4599000.0, 0.0115478000)
    assert_close(kg, 0.09594861261873211)


def test_CSP_gas_dense():
    kgs = [Stiel_Thodos_dense(T=378.15, MW=44.013, Tc=309.6, Pc=72.4E5, Vc=97.4E-6, Zc=0.274, Vm=i, kg=2.34E-2) for i in [144E-6, 24E-6, 240E-6]]
    kgs_exp = [0.041245574404863684, 0.9158718777539487, 0.03258313269922979]
    assert_allclose(kgs, kgs_exp)


    kgs = [Eli_Hanley_dense(T=T, MW=42.081, Tc=364.9, Vc=1.81E-4, Zc=0.274, omega=0.144, Cvm=82.70, Vm=1.721E-4) for T in [473., 900]]
    kgs_exp = [0.06038475936515042, 0.08987438807653142]
    assert_allclose(kgs, kgs_exp)

    kg = Eli_Hanley_dense(700, MW=42.081, Tc=364.9, Vc=1.81E-4, Zc=0.274, omega=0.144, Cvm=82.70, Vm=3.721E-4)
    assert_allclose(kg, 0.06953791121177173)

    kg = Chung_dense(T=473., MW=42.081, Tc=364.9, Vc=184.6E-6, omega=0.142, Cvm=82.67, Vm=172.1E-6, mu=134E-7, dipole=0.4)
    assert_allclose(kg, 0.06160569232570781)


def test_DIPPR9H():
    k = DIPPR9H([0.258, 0.742], [0.1692, 0.1528])
    assert_allclose(k, 0.15657104706719646)

#    with pytest.raises(Exception):
#        DIPPR9H([0.258, 0.742], [0.1692])

def test_Filippov():
    kl = Filippov([0.258, 0.742], [0.1692, 0.1528])
    assert_allclose(kl, 0.15929167628799998)

    with pytest.raises(ValueError):
        Filippov([0.258], [0.1692, 0.1528])


def test_Lindsay_Bromley():
    kg = Lindsay_Bromley(323.15, [0.23, 0.77], [1.939E-2, 1.231E-2], [1.002E-5, 1.015E-5], [248.31, 248.93], [46.07, 50.49])
    assert_allclose(kg, 0.01390264417969313)

#    with pytest.raises(Exception):
#        Lindsay_Bromley(323.15, [0.23], [1.939E-2, 1.231E-2], [1.002E-5, 1.015E-5], [248.31, 248.93], [46.07, 50.49])


