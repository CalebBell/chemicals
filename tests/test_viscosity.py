# -*- coding: utf-8 -*-
r"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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

from random import uniform
import pytest
from math import log, log10
import numpy as np
import pandas as pd
from fluids.numerics import assert_close, assert_close1d, derivative
from fluids.constants import psi, atm, foot, lb
from fluids.core import R2K, F2K
from chemicals.utils import normalize, mixing_simple
from chemicals.viscosity import *
from chemicals.viscosity import (mu_data_Dutt_Prasad, mu_data_VN3, mu_data_VN2,
                                 mu_data_VN2E, mu_data_Perrys_8E_2_313, mu_data_Perrys_8E_2_312,
                                 mu_data_VDI_PPDS_7, mu_data_VDI_PPDS_8)
from chemicals.identifiers import check_CAS

### Check data integrity

def test_Dutt_Prasad_data():
    assert all([check_CAS(i) for i in mu_data_Dutt_Prasad.index])
    tots_calc = [mu_data_Dutt_Prasad[i].abs().sum() for i in ['A', 'B', 'C', 'Tmin', 'Tmax']]
    tots = [195.89260000000002, 65395.299999999996, 9849.1899999999987, 25952, 35016]
    assert_close1d(tots_calc, tots)

    assert mu_data_Dutt_Prasad.index.is_unique
    assert mu_data_Dutt_Prasad.shape == (100, 6)


def test_VN2E_data():
    assert all([check_CAS(i) for i in mu_data_VN2E.index])
    tots_calc = [mu_data_VN2E[i].abs().sum() for i in ['C', 'D', 'Tmin', 'Tmax']]
    tots = [567743298666.74878, 48.8643, 3690, 4860]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN2E.index.is_unique
    assert mu_data_VN2E.shape == (14, 6)


def test_VN3_data():
    assert all([check_CAS(i) for i in mu_data_VN3.index])
    tots_calc = [mu_data_VN3[i].abs().sum() for i in ['A', 'B', 'C', 'Tmin', 'Tmax']]
    tots = [645.18849999999998, 169572.65159999998, 50050.151870000002, 126495, 175660]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN3.index.is_unique
    assert mu_data_VN3.shape == (432, 7)

def test_VN2_data():
    assert all([check_CAS(i) for i in mu_data_VN3.index])
    tots_calc = [mu_data_VN2[i].abs().sum() for i in ['A', 'B', 'Tmin', 'Tmax']]
    tots = [674.10069999999996, 83331.98599999999, 39580, 47897]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN2.index.is_unique
    assert mu_data_VN2.shape == (135, 6)


def test_Perrys2_313_data():
    # All values calculated at Tmin and Tmax check out to at least 5E-3 precision
    # The rounding has some effect, but it is not worrying.
    assert all([check_CAS(i) for i in mu_data_Perrys_8E_2_313.index])
    tots_calc = [mu_data_Perrys_8E_2_313[i].abs().sum() for i in [u'C1', u'C2', u'C3', u'C4', u'C5', u'Tmin', u'Tmax']]
    tots = [9166.6971369999992, 615425.94497999991, 1125.5317557875198, 9.054869390623603e+34, 402.21244000000002, 72467.140000000014, 136954.85999999999]
    assert_close1d(tots_calc, tots)

    assert mu_data_Perrys_8E_2_313.index.is_unique
    assert mu_data_Perrys_8E_2_313.shape == (337, 8)


def test_Perrys2_312_data():
    # Argon, Difluoromethane, 1-Hexyne, Methylsilane, Nitrogen trifluoride,
    # Vinyl chloride all do not match on Tmax at 1E-3 - their models predict
    # ~1E-5 Pa*S, but listed values are ~1E-10 to 1E-12. Unsure of the cause.
    # All coumpounds match at 1E-3 for Tmin.

    assert all([check_CAS(i) for i in mu_data_Perrys_8E_2_312.index])
    tots_calc = [mu_data_Perrys_8E_2_312[i].abs().sum() for i in [u'C1', u'C2', u'C3', u'C4', u'Tmin', u'Tmax']]
    tots = [0.00019683902626010103, 250.10520100000002, 65862.829200000007, 191286, 74802.639999999999, 355064.37]
    assert_close1d(tots_calc, tots)

    assert mu_data_Perrys_8E_2_312.index.is_unique
    assert mu_data_Perrys_8E_2_312.shape == (345, 7)


def test_VDI_PPDS_7_data():
    assert all([check_CAS(i) for i in mu_data_VDI_PPDS_7.index])
    tots_calc = [mu_data_VDI_PPDS_7[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E']]
    tots = [507.14607000000001, 1680.7624099999998, 165461.14259999999, 46770.887000000002, 0.057384780000000003]
    assert_close1d(tots_calc, tots)

    assert mu_data_VDI_PPDS_7.index.is_unique
    assert mu_data_VDI_PPDS_7.shape == (271, 7)

def test_VDI_PPDS_8_data():
    # Coefficients for water are incorrect - obtained an average deviation of 150%!
    assert all([check_CAS(i) for i in mu_data_VDI_PPDS_8.index])
    assert mu_data_VDI_PPDS_8.index.is_unique
    assert mu_data_VDI_PPDS_8.shape == (274, 6)

    tots_calc = [mu_data_VDI_PPDS_8[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E']]
    tots = [0.00032879559999999999, 9.5561339999999995e-06, 2.8377710000000001e-09, 2.8713399999999998e-12, 2.8409200000000004e-15]
    assert_close1d(tots_calc, tots)


def test_ViswanathNatarajan():
    mu = Viswanath_Natarajan_2(348.15, -5.9719-log(100), 1007.0)
    assert_close(mu, 0.00045983686956829517)

    mu = Viswanath_Natarajan_2_exponential(298.15, 4900800,  -3.8075)
    assert_close(mu, 0.0018571903840928496)

    mu = Viswanath_Natarajan_3(298.15, -2.7173-log10(1000), -1071.18, -129.51)
    assert_close(mu, 0.0006129806445142112)

def test_PPDS9():
    mu = PPDS9(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    assert_close(mu, 0.00035091137378230684, rtol=1e-13)

    coeffs = (1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    # normal region
    mu_expect = PPDS9(400.0, *coeffs)
    dmu, mu = dPPDS9_dT(400.0, *coeffs)
    assert_close(mu, mu_expect, rtol=1e-13)
    assert_close(dmu, -3.186540635882627e-06, rtol=1e-10)
    assert_close(derivative(PPDS9, 400.0, args=coeffs, dx=1e-4), dmu, rtol=1e-9)

    mu_expect = PPDS9(5.0, *coeffs)
    dmu, mu = dPPDS9_dT(5.0, *coeffs)
    assert_close(mu, mu_expect, rtol=1e-13)
    assert_close(dmu, 1126796480623.1184, rtol=1e-10)
    assert_close(derivative(PPDS9, 5.0, args=coeffs, dx=1e-5), dmu, rtol=1e-9)

    # Check can go super low T without overflow
    coeffs = [1.20479, 0.6058, 216.325, 2.278, 8.756e-05]
    assert PPDS9(3, *coeffs) > 1e10
    assert PPDS9(2, *coeffs) > 1e10
    assert PPDS9(1, *coeffs) > 1e10

    dPPDS9_dT(2, *coeffs)
    dPPDS9_dT(1, *coeffs)
    dPPDS9_dT(3, *coeffs)

def test_Letsou_Stiel():
    # Checked 2017-03-05
    mu = Letsou_Stiel(400., 46.07, 516.25, 6.383E6, 0.6371)
    assert_close(mu, 0.0002036150875308151)

    # Butane, custom case; vs 0.000166 exp
    mu = Letsou_Stiel(298.15, 58.1222, 425.12, 3796000.0, 0.193)
    assert_close(mu, 0.00015385075092073157)


def test_Przedziecki_Sridhar():
    # Reid (1983) toluene, 383K
    mu = Przedziecki_Sridhar(383., 178., 591.8, 41E5, 316E-6, 95E-6, .263, 92.14)
    assert_close(mu, 0.0002198147995603383)


def test_Lucas():
    # methylcyclohexane
    mu = Lucas(300., 500E5, 572.2, 34.7E5, 0.236, 0, 0.00068)
    assert_close(mu, 0.0010683738499316518)

    # Psat > P:
    mu = Lucas(300., 500E5, 572.2, 34.7E5, 0.236, 501E5, 0.00068)
    assert_close(mu, 0.00068)


def test_Stiel_Thodos():
    # CCl4
    mu = Stiel_Thodos(300., 556.35, 4.5596E6, 153.8)
    assert_close(mu, 1.0408926223608723e-05)

    # Tr > 1.5
    mu = Stiel_Thodos(900., 556.35, 4.5596E6, 153.8)
    assert_close(mu, 2.899111242556782e-05)


def test_Yoon_Thodos():
    # CCl4
    mu = Yoon_Thodos(300., 556.35, 4.5596E6, 153.8)
    assert_close(mu, 1.0194885727776819e-05)


def test_viscosity_gas_Gharagheizi():
    mu = viscosity_gas_Gharagheizi(120., 190.564, 45.99E5, 16.04246)
    assert_close(mu, 5.215761625399613e-06)


def test_Lucas_gas():
    mu = Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)
    assert_close(mu, 1.7822676912698928e-05)

    mu = Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=None)
    assert_close(mu, 1.371116974367763e-05)

    mu = Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=8)
    assert_close(mu, 1.7811559961984407e-05)

    # Helium, testing Q
    mu = Lucas_gas(T=6, Tc=5.1889, Pc=226968.0, Zc=0.3014, MW=4.002602, CASRN='7440-59-7')
    assert_close(mu, 1.3042945737346396e-06)

    mu = Lucas_gas(T=150, Tc=5.1889, Pc=226968.0, Zc=0.3014, MW=4.002602, CASRN='7440-59-7')
    assert_close(mu, 1.2558477184738118e-05)


def test_Wilke():
    mu = Wilke([0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07])
    assert_close(mu, 9.701614885866193e-06, rtol=1e-10)

#    with pytest.raises(Exception):
#        Wilke([0.05], [1.34E-5, 9.5029E-6], [64.06, 46.07])

    mu = Wilke_large([0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07])
    assert_close(mu, 9.701614885866193e-06, rtol=1e-10)

    mu = Wilke_prefactored([0.05, 0.95], [1.34E-5, 9.5029E-6], *Wilke_prefactors([64.06, 46.07]))
    assert_close(mu, 9.701614885866193e-06, rtol=1e-10)

    # Large composition test
    zs = [0.10456352460469782, 0.10472506156674823, 0.10347781516834291, 1.4089716440797791e-05, 0.10488254481011455, 0.10078888107401028, 0.09902003237540975, 0.09045109410107957, 0.08642540418108867, 0.1043609364553231, 0.10129061594674436]
    mus = [1.1601665408586192e-05, 9.408370570946896e-06, 8.19709294177777e-06, 1.4314548719058091e-05, 1.5057441002481923e-05, 7.5434795308593725e-06, 7.447082353139856e-06, 7.0365592301967965e-06, 6.720364621681796e-06, 1.2157004301638695e-05, 1.3006463728382868e-05]
    MWs = [16.04246, 30.06904, 44.09562, 2.01588, 44.0095, 58.1222, 58.1222, 72.14878, 72.14878, 34.08088, 64.0638]

    # Make a large set of data, but don't actually use it anywhere; for easy of bencharmking
    zs_new = []
    mus_new = []
    MWs_new = []
    for i in range(20):
        zs_new.extend(zs)
        mus_new.extend(mus)
        MWs_new.extend(MWs)
    zs_large = normalize(zs_new)
    mus_large = mus_new
    MWs_large = MWs_new

    mu_expect = 9.282311227289109e-06
    assert_close(Wilke(zs, mus, MWs), mu_expect, rtol=1e-10)
    assert_close(Wilke_large(zs, mus, MWs), mu_expect, rtol=1e-10)
    prefactors = Wilke_prefactors(MWs)
    assert_close(Wilke_prefactored(zs, mus, *prefactors), mu_expect, rtol=1e-10)

    # Test that the prefactors really work - use a different composition
    zs_diff = normalize([.1, .2,.3, .4, .5, .6, .7, .8, .9, .10, .2])
    mu_expect = 8.238656569251283e-06
    assert_close(Wilke(zs_diff, mus, MWs), mu_expect, rtol=1e-10)
    assert_close(Wilke_large(zs_diff, mus, MWs), mu_expect, rtol=1e-10)
    assert_close(Wilke_prefactored(zs_diff, mus, *prefactors), mu_expect, rtol=1e-10)



def test_Herning_Zipperer():
    mu = Herning_Zipperer([0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07])
    assert_close(mu, 9.730630997268096e-06)

#    with pytest.raises(Exception):
#        Herning_Zipperer([0.05], [1.34E-5, 9.5029E-6], [64.06, 46.07])

    zs = [0.5, 0.25, 0.25]
    mus =  [1.78e-05, 1.12e-05, 9.35e-06]
    MWs = [28.0134, 16.043, 30.07]
    MWs_roots = [i**0.5 for i in MWs]
    mu_root = Herning_Zipperer(zs, mus, MWs, MWs_roots)
    assert_close(mu_root, 1.4174908599465168e-05, rtol=1e-12)
    mu_root = Herning_Zipperer(zs, mus, None, MWs_roots)
    assert_close(mu_root, 1.4174908599465168e-05, rtol=1e-12)

def test_Brockaw():
    mu = Brokaw(308.2, [0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07], [0.42, 0.19], [347, 432])
    assert_close(mu, 9.699085099801568e-06)

#    with pytest.raises(Exception):
#        Brokaw(308.2, [0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07], [0.42, 0.19], [347, 432])

    # Test < 0.1 MD
    mu = Brokaw(308.2, [0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07], [0.42, 0.05], [347, 432])
    assert_close(mu, 9.70504431025103e-06)


def test_round_whole_even():
    from chemicals.viscosity import _round_whole_even
    assert _round_whole_even(116.4) == 116
    assert _round_whole_even(116.6) == 117
    assert _round_whole_even(116.5) == 116
    assert _round_whole_even(115.5) == 116


def test_viscosity_index():
    # 4/5 examples from the official document, one custom example
    assert_close(92.4296472453428, viscosity_index(73.3E-6, 8.86E-6))
    assert_close(156.42348293257797, viscosity_index(22.83E-6, 5.05E-6))
    assert_close(111.30701701381422, viscosity_index(53.47E-6, 7.80E-6))
    assert_close(92.03329369797858, viscosity_index(73.5E-6, 8.86E-6))
    assert_close(192.9975428057893, viscosity_index(1000E-6, 100E-6)) # Custom
    assert 193 == viscosity_index(1000E-6, 100E-6, rounding=True) # custom, rounded
    assert 92 == viscosity_index(73.3E-6, 8.86E-6, rounding=True)
    assert None == viscosity_index(3E-6, 1.5E-6)


def test_Lorentz_Bray_Clarke():
    # Made up example
    T = 300.0
    P = 1e6
    zs = [.4, .3, .3]
    MWs = [16.04246, 30.06904, 44.09562]
    Tcs = [190.564, 305.32, 369.83]
    Pcs = [4599000.0, 4872000.0, 4248000.0]
    Vcs = [9.86e-05, 0.0001455, 0.0002]
    Vm = 0.002302491921416089

    mu = Lorentz_Bray_Clarke(T, P, Vm, zs, MWs, Tcs, Pcs, Vcs)
    assert_close(mu, 9.925488946486405e-06, rtol=1e-6)

    #  2,000 psig and 160°F.
    zs = [0.875, 0.083, 0.021, 0.006, 0.008, 0.003, 0.002, 0.001, 0.001]
    MWs = [16.04, 30.07, 44.09, 58.12, 58.12, 72.15, 72.15, 86.17, 114.00]
    Pcs = [667.8*psi, 707.8*psi, 616.3*psi, 529.1*psi, 550.7*psi, 490.4*psi, 488.6*psi, 436.9*psi, 360.6*psi]
    Tcs = [R2K(343.0), R2K(549.8), R2K(665.7), R2K(734.7), R2K(765.3), R2K(828.8), R2K(845.4), R2K(913.4), R2K(1023.9)]
    Vcs = [1.590*foot**3/lb, 2.370*foot**3/lb, 3.250*foot**3/lb, 4.208*foot**3/lb, 4.080*foot**3/lb, 4.899*foot**3/lb, 4.870*foot**3/lb, 5.929*foot**3/lb, 7.882*foot**3/lb]
    P = atm + 2000*psi
    T = F2K(160.0)

    MW = mixing_simple(zs, MWs)
    rho_mass = 6.74*lb/foot**3
    rhom = rho_mass/MW
    Vm = 1.0/rhom

    mu = Lorentz_Bray_Clarke(T, P, Vm, zs, MWs, Tcs, Pcs, Vcs)
    assert_close(mu, 1.636032602394696e-05)

def test_viscosity_converter():
    # Barbey - todo viscosity_converter(95, 'barbey', 'parlin cup #7')

    visc = viscosity_converter(8.79, 'engler', 'parlin cup #7')
    assert type(visc) is float
    assert_close(visc, 52.5)

    # seconds/degrees string inputs and capitals
    visc = viscosity_converter(8.79, 'degrees engler', 'seconds parlin cup #7')
    assert type(visc) is float
    assert_close(visc, 52.5)

    visc = viscosity_converter(8.79, '    degrees engler', 'seconds parlin cup #7    ')
    assert type(visc) is float
    assert_close(visc, 52.5)

    visc = viscosity_converter(8.79, 'Engler', 'PARLIN cup #7')
    assert type(visc) is float
    assert_close(visc, 52.5)


    visc = viscosity_converter(8.78, 'engler', 'parlin cup #7')
    assert type(visc) is float
    assert_close(visc, 52.45389001785669)

    visc = viscosity_converter(5.91, 'engler', 'parlin cup #7', True)
    assert type(visc) is float
    assert_close(visc, 39.96017612902695)

    with pytest.raises(Exception):
        # limit is 5.92, but even that fails due to float conversion
        viscosity_converter(5.91, 'engler', 'parlin cup #7', extrapolate=False)

    viscosity_converter(5.919999, 'engler', 'parlin cup #7')
    with pytest.raises(Exception):
        # too little
        viscosity_converter(5.91999, 'engler', 'parlin cup #7')

    with pytest.raises(Exception):
        viscosity_converter(8.79, 'NOTAREALSCALE', 'kinematic viscosity')
    with pytest.raises(Exception):
        viscosity_converter(8.79, 'kinematic viscosity', 'NOTAREALSCALE')

    nu = viscosity_converter(8.79, 'engler', 'kinematic viscosity')
    assert type(nu) is float
    assert_close(nu, 6.5E-5)

    with pytest.raises(Exception):
        viscosity_converter(8.79, 'pratt lambert e', 'kinematic viscosity')

    t = viscosity_converter(0.0002925, 'kinematic viscosity', 'pratt lambert g')
    assert type(nu) is float
    assert_close(t, 7.697368421052632)
    nu = viscosity_converter(7.697368421052632, 'pratt lambert g', 'kinematic viscosity', )
    assert type(nu) is float
    assert_close(nu, .0002925)

    with pytest.raises(Exception):
        viscosity_converter(0.00002925, 'kinematic viscosity', 'pratt lambert g')
    viscosity_converter(0.00002925, 'kinematic viscosity', 'pratt lambert g', True)

    # Too low on the input side
    with pytest.raises(Exception):
        viscosity_converter(6, 'pratt lambert g', 'kinematic viscosity')
    viscosity_converter(6, 'pratt lambert g', 'kinematic viscosity', True)

    nu = viscosity_converter(700, 'Saybolt Universal Seconds', 'kinematic viscosity')
    assert type(nu) is float
    assert_close(nu, 0.00015108914751515542)

    t = viscosity_converter(0.00015108914751515542, 'kinematic viscosity', 'Saybolt Universal Seconds')
    assert type(t) is float
    assert_close(t, 700)

    # Barbey custom tests
    nu = viscosity_converter(483, 'barbey', 'kinematic viscosity')*1E6
    assert type(nu) is float
    assert_close(nu, 13.1)
    nu = viscosity_converter(1.4, 'barbey', 'kinematic viscosity')*1E6
    assert type(nu) is float
    assert_close(nu, 4400)
    nu = viscosity_converter(6200, 'barbey', 'kinematic viscosity')*1E6
    assert type(nu) is float
    assert_close(nu, 1)
    nu = viscosity_converter(6300, 'barbey', 'kinematic viscosity', extrapolate=True)*1E6
    assert type(nu) is float
    assert_close(nu, 0.984, rtol=1E-3)
    # The extrapolation when barbey is known and under 1.4 is not working and goes the wrong direction
    # viscosity_converter(1.39999, 'barbey', 'kinematic viscosity', extrapolate=True)*1E6

@pytest.mark.slow
@pytest.mark.fuzz
def test_viscosity_converter_fuzz_SUS():
    for i in range(20):
        # fuzz the numerical solver for SUS a bit -increase to try harder, but
        # all efforts show the function is monotonic. It turns negative at 25.5
        # and stops working on the high side at 8000000000000. Plenty of room
        # for newton's method to converge!
        SUS = uniform(31.0, 20000.0)
        nu = viscosity_converter(SUS, 'Saybolt Universal Seconds', 'kinematic viscosity')
        SUS2 = viscosity_converter(nu, 'kinematic viscosity', 'Saybolt Universal Seconds')
        assert_close(SUS, SUS2)



def test_mu_IAPWS():
    # # Main PDF table 4
    Ts_mu20 = [298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15]
    rhos_mu20 = [998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400]
    mus20_listed = [889.735100E-6, 1437.649467E-6, 307.883622E-6, 14.538324E-6, 217.685358E-6, 32.619287E-6, 35.802262E-6, 77.430195E-6, 44.217245E-6, 47.640433E-6, 64.154608E-6]
    for i in range(len(Ts_mu20)):
        mu = mu_IAPWS(T=Ts_mu20[i], rho=rhos_mu20[i])
        assert_close(mus20_listed[i], mu, rtol=4e-8)

    assert_close(mu_IAPWS(T=647.35, rho=122.0, drho_dP=17.109308489109e-6), 2.55206768504972e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=122.0, drho_dP=17.109308489109e-6, drho_dP_Tr=2.936891667997e-6)
    assert_close(mu, 2.552067683647617e-05, rtol=1e-13)

    mu = mu_IAPWS(T=647.35, rho=222, drho_dP=175.456980972231e-6)
    assert_close(mu, 3.133759043936869e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=222, drho_dP=175.456980972231e-6, drho_dP_Tr=3.119177410324e-6)
    assert_close(mu, 3.1337589197275484e-05, rtol=1e-13)

    mu = mu_IAPWS(T=647.35, rho=272, drho_dP=1508.2800389184448e-6)
    assert_close(mu, 3.6228143762668687e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=272, drho_dP=1508.2800389184448e-6, drho_dP_Tr=2.999611040849e-6)
    assert_close(mu, 3.622814313612717e-05, rtol=1e-13)


    mu = mu_IAPWS(T=647.35, rho=322, drho_dP=1.213641949033e-2, drho_dP_Tr=2.751438963343e-6)
    assert_close(mu, 4.296157881023829e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=322, drho_dP=1.213641949033e-2)
    assert_close(mu, 4.2961578738287014e-05, rtol=1e-13)


    mu = mu_IAPWS(T=647.35, rho=372, drho_dP=1245.917204367233E-6, drho_dP_Tr=2.415440238773e-6)
    assert_close(mu, 4.568820447470762e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=372, drho_dP=1245.917204367233E-6)
    assert_close(mu, 4.5688204316365544e-05, rtol=1e-13)

    mu = mu_IAPWS(T=647.35, rho=422, drho_dP=130.393537965256e-6, drho_dP_Tr=2.046542440571e-6)
    assert_close(mu, 4.943625601494998e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=422, drho_dP=130.393537965256e-6)
    assert_close(mu, 4.943625789265867e-05, rtol=1e-13)

    mu = mu_IAPWS(T=647.35, rho=20, drho_dP=4.074182233978e-6)
    assert_close(mu, 2.3253942194448227e-05, rtol=1e-13) # rhor <= 0.310559006
    mu = mu_IAPWS(T=647.35, rho=20, drho_dP=4.074182233978e-6, drho_dP_Tr=2.352881641109e-6)
    assert_close(mu, 2.3253942194448267e-05, rtol=1e-13)


    mu = mu_IAPWS(T=647.35, rho=600, drho_dP=4.28639862246e-6, drho_dP_Tr=9.903267886625e-7)
    assert_close(mu, 7.001987553170616e-05, rtol=1e-13)
    mu = mu_IAPWS(T=647.35, rho=600, drho_dP=4.28639862246e-6)
    assert_close(mu, 7.001987566162558e-05, rtol=1e-13)

    # Test case with zero
    assert_close(mu_IAPWS(347.0, 975.5266664069043, 4.439512743107522e-07), 0.00038316609714314585, rtol=1e-13)

def test_Twu_1985():
    from chemicals.viscosity import Twu_1985_internal
    mu = Twu_1985_internal(T=609.67, Tb=1210.17, SG=0.8964)
    assert_close(mu, 9.195790397643691, rtol=1e-13)

    mu = Twu_1985(T=R2K(609.67), Tb=R2K(1210.17), rho=0.8964*999.0170824078306)
    assert_close(mu, 0.008235004218042592, rtol=1e-13)


def test_mu_air_lemmon():
    assert_close(mu_air_lemmon(300.0, 40.10292351061862), 1.853715185567247e-05, rtol=1e-13)

    # Values in a check table 5 of [1]_.
    assert round(mu_air_lemmon(100.0, 0.0), 11) == 7.09559e-6
    assert round(mu_air_lemmon(300.0, 0.0), 10) == 18.5230e-6
    assert round(mu_air_lemmon(100.0, 28e3), 9) == 107.923e-6
    assert round(mu_air_lemmon(200.0, 10e3), 10) == 21.1392e-6
    assert round(mu_air_lemmon(300.0, 5e3), 10) == 21.3241e-6
    assert round(mu_air_lemmon(132.64, 10.4e3), 10) == 17.7623e-6