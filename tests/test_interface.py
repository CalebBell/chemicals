# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.
"""

import pytest
import numpy as np
import pandas as pd
from fluids.numerics import assert_close, assert_close1d
from chemicals.utils import rho_to_Vm
from chemicals.interface import *
from chemicals.interface import (sigma_data_Mulero_Cachadina, sigma_data_Jasper_Lange, sigma_data_Somayajulu,
                                 sigma_data_VDI_PPDS_11, sigma_data_Somayajulu2)
from chemicals.identifiers import check_CAS

def test_sigma_IAPWS():
    assert_close(sigma_IAPWS(300.), 0.07168596252716256, rtol=1e-13)
    assert_close(sigma_IAPWS(450.), 0.04289149915650591, rtol=1e-13)
    assert_close(sigma_IAPWS(600.), 0.00837561087288565, rtol=1e-13)

def test_CSP():
    # p-dichloribenzene at 412.15 K, from DIPPR; value differs due to a slight
    # difference in method.
    sigma1 = Brock_Bird(412.15, 447.3, 685, 3.952E6)
    assert_close(sigma1, 0.02208448325192495)

    # Chlorobenzene from Poling, as compared with a % error value at 293 K.
    sigma1 = Brock_Bird(293.15, 404.75, 633.0, 4530000.0)
    assert_close(sigma1, 0.032985686413713036)

    # TODO: Find parameters where Brock Bird is negative

    # Chlorobenzene from Poling, as compared with a % error value at 293 K.
    sigma1 = Pitzer_sigma(293., 633.0, 4530000.0, 0.249)
    assert_close(sigma1, 0.03458453513446387)

    sigma1 = Sastri_Rao(293.15, 404.75, 633.0, 4530000.0)
    sigma2 = Sastri_Rao(293.15, 404.75, 633.0, 4530000.0, chemicaltype='alcohol')
    sigma3 = Sastri_Rao(293.15, 404.75, 633.0, 4530000.0, chemicaltype='acid')
    sigmas = [0.03234567739694441, 0.023255104102733407, 0.02558993464948134]
    assert_close1d([sigma1, sigma2, sigma3], sigmas)

    sigma = Zuo_Stenby(293., 633.0, 4530000.0, 0.249)
    assert_close(sigma, 0.03345569011871088)

    # 1-butanol, as compared to value in CRC Handbook of 0.02493.
    sigma = Hakim_Steinberg_Stiel(298.15, 563.0, 4414000.0, 0.59, StielPolar=-0.07872)
    assert_close(sigma, 0.021907902575190447)

    sigma_calc = Miqueu(300., 340.1, 0.000199, 0.1687)
    assert_close(sigma_calc, 0.003474100774091376)

    sigma_calc = Mersmann_Kind_sigma(298.15, 164.15, 328.25, 497.1, 3430000.0)
    assert_close(0.016744311449290426, sigma_calc)


def test_Aleem():
    st = Aleem(T=90, MW=16.04246, Tb=111.6, rhol=458.7, Hvap_Tb=510870., Cpl=2465.)
    assert_close(st, 0.01669970221165325)

    # Test the "critical" point (where the correlation predicts sigma=0)
    st_at_Tc = Aleem(T=318.8494929006085, MW=16.04246, Tb=111.6, rhol=458.7,
                     Hvap_Tb=510870.,  Cpl=2465.)
    assert_close(st_at_Tc, 0, atol=1E-12)


def test_REFPROP():
    sigma = REFPROP_sigma(298.15, 647.096, -0.1306, 2.471, 0.2151, 1.233)
    assert_close(sigma, 0.07205503890847453)


def test_Somayajulu():
    sigma = Somayajulu(300, 647.126, 232.713514, -140.18645, -4.890098)
    assert_close(sigma, 0.07166386387996757)


def test_Jasper():
    sigma = Jasper(298.15, 24, 0.0773)
    assert_close(sigma, 0.0220675)


def test_data():
    tot = sum([sigma_data_Mulero_Cachadina[i].sum() for i in sigma_data_Mulero_Cachadina.columns[1:]])
    assert_close(tot, 114350.07371931802)

    assert np.all(sigma_data_Mulero_Cachadina.columns == [u'Fluid', u'sigma0', u'n0',
                                                    u'sigma1', u'n1', u'sigma2',
                                                    u'n2', u'Tc', u'Tmin',
                                                    u'Tmax'])
    assert sigma_data_Mulero_Cachadina.shape == (115, 10)
    assert sigma_data_Mulero_Cachadina.index.is_unique

    tot = sum([sigma_data_Jasper_Lange[i].sum() for i in sigma_data_Jasper_Lange.columns[1:]])
    assert_close(tot, 343153.38953333395)

    assert np.all(sigma_data_Jasper_Lange.columns == [u'Name', u'a', u'b', u'Tmin', u'Tmax'])
    assert sigma_data_Jasper_Lange.shape == (522, 5)
    assert sigma_data_Jasper_Lange.index.is_unique

    tot = sum([sigma_data_Somayajulu[i].sum() for i in sigma_data_Somayajulu.columns[1:]])
    assert_close(tot, 38941.199955999997)

    assert np.all(sigma_data_Somayajulu.columns == [u'Chemical', u'Tt', u'Tc', u'A', u'B', u'C'])
    assert sigma_data_Somayajulu.shape == (64, 6)
    assert sigma_data_Somayajulu.index.is_unique

    tot = sum([sigma_data_Somayajulu2[i].sum() for i in sigma_data_Somayajulu2.columns[1:]])
    assert_close(tot, 39471.356771000006)
    assert np.all(sigma_data_Somayajulu2.columns == [u'Chemical', u'Tt', u'Tc', u'A', u'B', u'C'])
    assert sigma_data_Somayajulu2.shape == (64, 6)
    assert sigma_data_Somayajulu2.index.is_unique

def test_VDI_PPDS_11_data():
    """I believe there are no errors here."""
    for i in sigma_data_VDI_PPDS_11.index:
        assert check_CAS(i)

    assert sigma_data_VDI_PPDS_11.index.is_unique
    assert sigma_data_VDI_PPDS_11.shape == (272, 8)

    # Doing the sums on the arrays is faster but much uglier. Worth it?
    tots_calc = [sigma_data_VDI_PPDS_11[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E', u'Tc', u'Tm']]
    tots = [18.495069999999998, 336.69950000000006, 6.5941200000000002, 7.7347200000000003, 6.4262199999999998, 150142.28, 56917.699999999997]
    for calc, fixed in zip(tots_calc, tots):
        assert_close(calc, fixed)

def test_Weinaug_Katz():
    # sample test case
    sigma = Weinaug_Katz([5.1e-5, 7.2e-5], Vml=0.000125, Vmg=0.02011, xs=[.4, .6], ys=[.6, .4])
    assert_close(sigma, 0.06547479150776776)

    # pure component check it checks out
    Vml = rho_to_Vm(800.8088185536124, 100.15888)
    Vmg = rho_to_Vm(4.97865317223119, 100.15888)
    sigma = Weinaug_Katz([5.088443542210164e-05], Vml, Vmg, [1.0], [1.0])
    assert_close(sigma, 0.026721669606560042, rtol=1e-13)


def test_Winterfeld_Scriven_Davis():
    # The example is from [2]_; all results agree.
    # The original source has not been reviewed.
    # 16.06 mol% n-pentane, 83.94 mol% dichloromethane at 298.15 K.
    # sigmas are 15.47 and 28.77 respectively, rhos 8.61 kmol/m^3 and 15.53 kmol/m^3

    sigma = Winterfeld_Scriven_Davis([0.1606, 0.8394], [0.01547, 0.02877], [8610., 15530.])
    assert_close(sigma, 0.024967388450439817)

#    with pytest.raises(Exception):
#        Winterfeld_Scriven_Davis([0.1606, 0.8394, 0.118], [0.01547, 0.02877], [8610., 15530.])


def test_Diguilio_Teja():
    # Winterfeld_Scriven_Davis example, with Diguilio Teja. Calculated sigmas at Tbs are
    # 0.01424 and 0.02530.
    # Checked and edited 2017-01-30

    sigma = Diguilio_Teja(T=298.15, xs=[0.1606, 0.8394], sigmas_Tb=[0.01424, 0.02530], Tbs=[309.21, 312.95], Tcs=[469.7, 508.0])
    assert_close(sigma, 0.025716823875045505)

#    with pytest.raises(Exception):
#        Diguilio_Teja(T=298.15, xs=[0.1606, 0.8394, 0.118], sigmas_Tb=[0.01424, 0.02530], Tbs=[309.21, 312.95], Tcs=[469.7, 508.0])

    with pytest.raises(Exception):
         Diguilio_Teja(T=501.85, xs=[0.1606, 0.8394], sigmas_Tb=[0.01424, 0.02530], Tbs=[309.21, 312.95], Tcs=[469.7, 508.0])



def test_Meybodi_Daryasafar_Karimi():
    sigma = Meybodi_Daryasafar_Karimi(980, 760, 580, 914)
    assert_close(sigma, 0.02893598143089256)


def test_API10A32():
    from fluids.core import F2K, R2K
    assert_close(API10A32(T=F2K(60), Tc=R2K(1334), K_W=12.4), 29.577333312096968, rtol=1e-13)


def test_PPDS14():
    sigma = PPDS14(T=280, Tc=562.05, a0=0.0786269, a1=1.28646, a2=-0.112304)
    assert_close(sigma, 0.030559764256249854, rtol=1e-14)
    
def test_Watson_sigma():
    sigma = Watson_sigma(T=350.0, Tc=543.836, a1=-3.02417, a2=1.21792, a3=-5.26877e-9, a4=5.62659e-9, a5=-2.27553e-9)
    assert_close(sigma, 0.013834092660564925, rtol=1e-14)
    
def test_ISTExpansion():
    sigma = ISTExpansion(T=400.0, Tc=776.0, a1=0.037545, a2=0.0363288)
    assert_close(sigma, 0.02672100905515996, rtol=1e-13)
    
    sigma = ISTExpansion(T=400.0, Tc=776.0, a1=0.037545, a2=0.0363288, a3=1e-4, a4=1e-3, a5=1e-4)
    assert_close(sigma, 0.02679017489704166, rtol=1e-15)