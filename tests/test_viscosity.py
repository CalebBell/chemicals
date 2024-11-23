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

from math import log, log10, exp
from random import uniform
from random import Random

import pytest
from fluids.constants import atm, foot, lb, psi
from fluids.core import F2K, R2K
from fluids.numerics import assert_close, assert_close1d, derivative

from chemicals.utils import mixing_simple, normalize
from chemicals.viscosity import (
    PPDS5,
    PPDS9,
    Brokaw,
    Herning_Zipperer,
    Letsou_Stiel,
    Lorentz_Bray_Clarke,
    Lucas,
    Lucas_gas,
    Przedziecki_Sridhar,
    Stiel_Thodos,
    Twu_1985,
    Viswanath_Natarajan_2,
    Viswanath_Natarajan_2_exponential,
    Viswanath_Natarajan_3,
    Wilke,
    Wilke_large,
    Wilke_prefactored,
    Wilke_prefactors,
    Yoon_Thodos,
    dmu_Yaws_dT,
    dPPDS9_dT,
    mu_air_lemmon,
    mu_data_Dutt_Prasad,
    mu_data_Perrys_8E_2_312,
    mu_data_Perrys_8E_2_313,
    mu_data_VDI_PPDS_7,
    mu_data_VDI_PPDS_8,
    mu_data_VN2,
    mu_data_VN2E,
    mu_data_VN3,
    mu_IAPWS,
    mu_TDE,
    mu_Yaws,
    mu_Yaws_fitting_jacobian,
    viscosity_converter,
    viscosity_gas_Gharagheizi,
    viscosity_index,
    viscosity_scales, 
    viscosity_scales_linear, 
    viscosity_converter_limits,
)

### Check data integrity

def test_Dutt_Prasad_data():
    tots_calc = [mu_data_Dutt_Prasad[i].abs().sum() for i in ['A', 'B', 'C', 'Tmin', 'Tmax']]
    tots = [195.89260000000002, 65395.299999999996, 9849.1899999999987, 25952, 35016]
    assert_close1d(tots_calc, tots)

    assert mu_data_Dutt_Prasad.shape == (100, 6)


def test_VN2E_data():
    tots_calc = [mu_data_VN2E[i].abs().sum() for i in ['C', 'D', 'Tmin', 'Tmax']]
    tots = [567743298666.74878, 48.8643, 3690, 4860]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN2E.shape == (14, 6)


def test_VN3_data():
    tots_calc = [mu_data_VN3[i].abs().sum() for i in ['A', 'B', 'C', 'Tmin', 'Tmax']]
    tots = [645.18849999999998, 169572.65159999998, 50050.151870000002, 126495, 175660]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN3.shape == (432, 7)

def test_VN2_data():
    tots_calc = [mu_data_VN2[i].abs().sum() for i in ['A', 'B', 'Tmin', 'Tmax']]
    tots = [674.10069999999996, 83331.98599999999, 39580, 47897]
    assert_close1d(tots_calc, tots)

    assert mu_data_VN2.shape == (135, 6)


def test_Perrys2_313_data():
    # All values calculated at Tmin and Tmax check out to at least 5E-3 precision
    # The rounding has some effect, but it is not worrying.
    tots_calc = [mu_data_Perrys_8E_2_313[i].abs().sum() for i in ['C1', 'C2', 'C3', 'C4', 'C5', 'Tmin', 'Tmax']]
    tots = [9166.6971369999992, 615425.94497999991, 1125.5317557875198, 9.054869390623603e+34, 402.21244000000002, 72467.140000000014, 136954.85999999999]
    assert_close1d(tots_calc, tots)

    assert mu_data_Perrys_8E_2_313.shape == (337, 8)


def test_Perrys2_312_data():
    # Argon, Difluoromethane, 1-Hexyne, Methylsilane, Nitrogen trifluoride,
    # Vinyl chloride all do not match on Tmax at 1E-3 - their models predict
    # ~1E-5 Pa*S, but listed values are ~1E-10 to 1E-12. Unsure of the cause.
    # All coumpounds match at 1E-3 for Tmin.

    tots_calc = [mu_data_Perrys_8E_2_312[i].abs().sum() for i in ['C1', 'C2', 'C3', 'C4', 'Tmin', 'Tmax']]
    tots = [0.00019683902626010103, 250.10520100000002, 65862.829200000007, 191286, 74802.639999999999, 355064.37]
    assert_close1d(tots_calc, tots)

    assert mu_data_Perrys_8E_2_312.shape == (345, 7)


def test_VDI_PPDS_7_data():
    tots_calc = [mu_data_VDI_PPDS_7[i].abs().sum() for i in ['A', 'B', 'C', 'D', 'E']]
    tots = [507.14607000000001, 1680.7624099999998, 165461.14259999999, 46770.887000000002, 0.057384780000000003]
    assert_close1d(tots_calc, tots)

    assert mu_data_VDI_PPDS_7.shape == (271, 7)

def test_VDI_PPDS_8_data():
    # Coefficients for water are incorrect - obtained an average deviation of 150%!
    assert mu_data_VDI_PPDS_8.shape == (274, 6)

    tots_calc = [mu_data_VDI_PPDS_8[i].abs().sum() for i in ['A', 'B', 'C', 'D', 'E']]
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
    assert None is viscosity_index(3e-06, 1.5e-06)


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
    assert_close(visc, 52.45389001785669, rtol=2e-5)

    visc = viscosity_converter(5.91, 'engler', 'parlin cup #7', True)
    assert type(visc) is float
    assert_close(visc, 39.96017612902695, rtol=2e-5)

    with pytest.raises(Exception):
        # limit is 5.92, but even that fails due to float conversion
        viscosity_converter(5.91, 'engler', 'parlin cup #7', extrapolate=False)

    viscosity_converter(5.919999999, 'engler', 'parlin cup #7')
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

def test_viscosity_converter_spline_fits():
    # Saybolt Furol conversions
    assert_close(viscosity_converter(12.95, 'saybolt furol', 'kinematic viscosity'), 1.3099999999999998e-05)
    assert_close(viscosity_converter(100.7, 'saybolt furol', 'kinematic viscosity'), 0.00021999999999999998)
    assert_close(viscosity_converter(2000.0, 'saybolt furol', 'kinematic viscosity'), 0.004399999999999999)

    # Redwood Standard conversions
    assert_close(viscosity_converter(29.0, 'redwood standard', 'kinematic viscosity'), 1e-06)
    assert_close(viscosity_converter(592.0, 'redwood standard', 'kinematic viscosity'), 0.000154)
    assert_close(viscosity_converter(18400.0, 'redwood standard', 'kinematic viscosity'), 0.004399999999999999)

    # Redwood Admiralty conversions
    assert_close(viscosity_converter(5.1, 'redwood admiralty', 'kinematic viscosity'), 4.2999999999999995e-06)
    assert_close(viscosity_converter(64.6, 'redwood admiralty', 'kinematic viscosity'), 0.000154)
    assert_close(viscosity_converter(921.0, 'redwood admiralty', 'kinematic viscosity'), 0.0021999999999999997)

    # Engler conversions
    assert_close(viscosity_converter(1.0, 'engler', 'kinematic viscosity'), 1e-06)
    assert_close(viscosity_converter(20.45, 'engler', 'kinematic viscosity'), 0.000154)
    assert_close(viscosity_converter(584.0, 'engler', 'kinematic viscosity'), 0.004399999999999999)

    # Barbey conversions
    assert_close(viscosity_converter(6200.0, 'barbey', 'kinematic viscosity'), 1e-06)
    assert_close(viscosity_converter(40.3, 'barbey', 'kinematic viscosity'), 0.000154)
    assert_close(viscosity_converter(1.4, 'barbey', 'kinematic viscosity'), 0.004399999999999999)

    # Parlin Cup #7 conversions
    assert_close(viscosity_converter(40.0, 'parlin cup #7', 'kinematic viscosity'), 4.32e-05)
    assert_close(viscosity_converter(92.0, 'parlin cup #7', 'kinematic viscosity'), 0.00013199999999999998)
    assert_close(viscosity_converter(149.0, 'parlin cup #7', 'kinematic viscosity'), 0.00021999999999999998)

    # Parlin Cup #10 conversions
    assert_close(viscosity_converter(15.0, 'parlin cup #10', 'kinematic viscosity'), 6.5e-05)
    assert_close(viscosity_converter(108.0, 'parlin cup #10', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(860.0, 'parlin cup #10', 'kinematic viscosity'), 0.004399999999999999)

    # Parlin Cup #15 conversions
    assert_close(viscosity_converter(6.0, 'parlin cup #15', 'kinematic viscosity'), 6.5e-05)
    assert_close(viscosity_converter(24.0, 'parlin cup #15', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(203.0, 'parlin cup #15', 'kinematic viscosity'), 0.004399999999999999)

    # Parlin Cup #20 conversions
    assert_close(viscosity_converter(3.0, 'parlin cup #20', 'kinematic viscosity'), 6.5e-05)
    assert_close(viscosity_converter(9.0, 'parlin cup #20', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(70.0, 'parlin cup #20', 'kinematic viscosity'), 0.004399999999999999)

    # Ford Cup #3 conversions
    assert_close(viscosity_converter(30.0, 'ford cup #3', 'kinematic viscosity'), 6.5e-05)
    assert_close(viscosity_converter(218.0, 'ford cup #3', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(1715.0, 'ford cup #3', 'kinematic viscosity'), 0.004399999999999999)

    # Ford Cup #4 conversions
    assert_close(viscosity_converter(20.0, 'ford cup #4', 'kinematic viscosity'), 6.5e-05)
    assert_close(viscosity_converter(147.0, 'ford cup #4', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(1150.0, 'ford cup #4', 'kinematic viscosity'), 0.004399999999999999)

    # Mac Michael conversions
    assert_close(viscosity_converter(125.0, 'mac michael', 'kinematic viscosity'), 2.06e-05)
    assert_close(viscosity_converter(805.0, 'mac michael', 'kinematic viscosity'), 0.00033)
    assert_close(viscosity_converter(10500.0, 'mac michael', 'kinematic viscosity'), 0.004399999999999999)

    # Zahn Cup #1 conversions
    assert_close(viscosity_converter(38.0, 'zahn cup #1', 'kinematic viscosity'), 2.06e-05)
    assert_close(viscosity_converter(62.0, 'zahn cup #1', 'kinematic viscosity'), 5.4e-05)
    assert_close(viscosity_converter(90.0, 'zahn cup #1', 'kinematic viscosity'), 8.759999999999999e-05)

    # Zahn Cup #2 conversions
    assert_close(viscosity_converter(18.0, 'zahn cup #2', 'kinematic viscosity'), 2.06e-05)
    assert_close(viscosity_converter(46.0, 'zahn cup #2', 'kinematic viscosity'), 0.00010999999999999999)
    assert_close(viscosity_converter(88.0, 'zahn cup #2', 'kinematic viscosity'), 0.00021999999999999998)

    # Zahn Cup #3 conversions
    assert_close(viscosity_converter(22.5, 'zahn cup #3', 'kinematic viscosity'), 0.000154)
    assert_close(viscosity_converter(40.0, 'zahn cup #3', 'kinematic viscosity'), 0.00033)
    assert_close(viscosity_converter(75.0, 'zahn cup #3', 'kinematic viscosity'), 0.00066)

    # Zahn Cup #4 conversions
    assert_close(viscosity_converter(18.0, 'zahn cup #4', 'kinematic viscosity'), 0.000198)
    assert_close(viscosity_converter(41.0, 'zahn cup #4', 'kinematic viscosity'), 0.0005499999999999999)
    assert_close(viscosity_converter(77.0, 'zahn cup #4', 'kinematic viscosity'), 0.0010999999999999998)

    # Zahn Cup #5 conversions
    assert_close(viscosity_converter(13.0, 'zahn cup #5', 'kinematic viscosity'), 0.00021999999999999998)
    assert_close(viscosity_converter(43.0, 'zahn cup #5', 'kinematic viscosity'), 0.0008799999999999999)
    assert_close(viscosity_converter(96.0, 'zahn cup #5', 'kinematic viscosity'), 0.00198)

    # Demmier #1 conversions
    assert_close(viscosity_converter(1.3, 'demmier #1', 'kinematic viscosity'), 4.2999999999999995e-06)
    assert_close(viscosity_converter(55.0, 'demmier #1', 'kinematic viscosity'), 0.000176)
    assert_close(viscosity_converter(1370.0, 'demmier #1', 'kinematic viscosity'), 0.004399999999999999)

    # Demmier #10 conversions
    assert_close(viscosity_converter(1.0, 'demmier #10', 'kinematic viscosity'), 3.21e-05)
    assert_close(viscosity_converter(13.7, 'demmier #10', 'kinematic viscosity'), 0.00043999999999999996)
    assert_close(viscosity_converter(137.0, 'demmier #10', 'kinematic viscosity'), 0.004399999999999999)

    # Stormer 100G Load conversions
    assert_close(viscosity_converter(2.6, 'stormer 100g load', 'kinematic viscosity'), 7.4e-06)
    assert_close(viscosity_converter(70.0, 'stormer 100g load', 'kinematic viscosity'), 0.000198)
    assert_close(viscosity_converter(1540.0, 'stormer 100g load', 'kinematic viscosity'), 0.004399999999999999)

    # Pratt Lambert F conversions
    assert_close(viscosity_converter(7.0, 'pratt lambert f', 'kinematic viscosity'), 8.759999999999999e-05)
    assert_close(viscosity_converter(33.7, 'pratt lambert f', 'kinematic viscosity'), 0.00066)
    assert_close(viscosity_converter(234.0, 'pratt lambert f', 'kinematic viscosity'), 0.004399999999999999)

def test_viscosity_converter_linear_ones():
    assert_close(viscosity_converter(35, 'american can', 'kinematic viscosity'), 0.0001225)
    assert_close(viscosity_converter(70, 'american can', 'kinematic viscosity'), 0.000245)
    assert_close(viscosity_converter(60, 'astm 0.07', 'kinematic viscosity'), 8.4e-05)
    assert_close(viscosity_converter(120, 'astm 0.07', 'kinematic viscosity'), 0.000168)
    assert_close(viscosity_converter(25, 'astm 0.10', 'kinematic viscosity'), 0.00011999999999999999)
    assert_close(viscosity_converter(50, 'astm 0.10', 'kinematic viscosity'), 0.00023999999999999998)
    assert_close(viscosity_converter(9, 'astm 0.15', 'kinematic viscosity'), 0.00018899999999999999)
    assert_close(viscosity_converter(18, 'astm 0.15', 'kinematic viscosity'), 0.00037799999999999997)
    assert_close(viscosity_converter(5, 'astm 0.20', 'kinematic viscosity'), 0.000305)
    assert_close(viscosity_converter(10, 'astm 0.20', 'kinematic viscosity'), 0.00061)
    assert_close(viscosity_converter(4, 'astm 0.25', 'kinematic viscosity'), 0.00056)
    assert_close(viscosity_converter(8, 'astm 0.25', 'kinematic viscosity'), 0.00112)
    assert_close(viscosity_converter(10, 'a&w b', 'kinematic viscosity'), 0.000185)
    assert_close(viscosity_converter(20, 'a&w b', 'kinematic viscosity'), 0.00037)
    assert_close(viscosity_converter(12, 'a&w crucible', 'kinematic viscosity'), 0.00014039999999999997)
    assert_close(viscosity_converter(24, 'a&w crucible', 'kinematic viscosity'), 0.00028079999999999994)
    assert_close(viscosity_converter(39, 'caspers tin plate', 'kinematic viscosity'), 0.0001404)
    assert_close(viscosity_converter(78, 'caspers tin plate', 'kinematic viscosity'), 0.0002808)
    assert_close(viscosity_converter(12, 'continental can', 'kinematic viscosity'), 3.9599999999999994e-05)
    assert_close(viscosity_converter(24, 'continental can', 'kinematic viscosity'), 7.919999999999999e-05)
    assert_close(viscosity_converter(12, 'crown cork and seal', 'kinematic viscosity'), 3.9599999999999994e-05)
    assert_close(viscosity_converter(24, 'crown cork and seal', 'kinematic viscosity'), 7.919999999999999e-05)
    assert_close(viscosity_converter(24, 'murphy varnish', 'kinematic viscosity'), 7.44e-05)
    assert_close(viscosity_converter(48, 'murphy varnish', 'kinematic viscosity'), 0.0001488)
    assert_close(viscosity_converter(15, 'parlin cup #25', 'kinematic viscosity'), 0.0021)
    assert_close(viscosity_converter(30, 'parlin cup #25', 'kinematic viscosity'), 0.0042)
    assert_close(viscosity_converter(10, 'parlin cup #30', 'kinematic viscosity'), 0.0026)
    assert_close(viscosity_converter(20, 'parlin cup #30', 'kinematic viscosity'), 0.0052)
    assert_close(viscosity_converter(70, 'pratt lambert a', 'kinematic viscosity'), 4.2699999999999994e-05)
    assert_close(viscosity_converter(140, 'pratt lambert a', 'kinematic viscosity'), 8.539999999999999e-05)
    assert_close(viscosity_converter(60, 'pratt lambert b', 'kinematic viscosity'), 7.32e-05)
    assert_close(viscosity_converter(120, 'pratt lambert b', 'kinematic viscosity'), 0.0001464)
    assert_close(viscosity_converter(40, 'pratt lambert c', 'kinematic viscosity'), 9.72e-05)
    assert_close(viscosity_converter(80, 'pratt lambert c', 'kinematic viscosity'), 0.0001944)
    assert_close(viscosity_converter(25, 'pratt lambert d', 'kinematic viscosity'), 0.00012174999999999999)
    assert_close(viscosity_converter(50, 'pratt lambert d', 'kinematic viscosity'), 0.00024349999999999998)
    assert_close(viscosity_converter(15, 'pratt lambert e', 'kinematic viscosity'), 0.00014625)
    assert_close(viscosity_converter(30, 'pratt lambert e', 'kinematic viscosity'), 0.0002925)
    assert_close(viscosity_converter(7, 'pratt lambert g', 'kinematic viscosity'), 0.000266)
    assert_close(viscosity_converter(14, 'pratt lambert g', 'kinematic viscosity'), 0.000532)
    assert_close(viscosity_converter(5, 'pratt lambert h', 'kinematic viscosity'), 0.00037999999999999997)
    assert_close(viscosity_converter(10, 'pratt lambert h', 'kinematic viscosity'), 0.0007599999999999999)
    assert_close(viscosity_converter(4, 'pratt lambert i', 'kinematic viscosity'), 0.0006079999999999999)
    assert_close(viscosity_converter(8, 'pratt lambert i', 'kinematic viscosity'), 0.0012159999999999999)
    assert_close(viscosity_converter(20, 'scott', 'kinematic viscosity'), 3.2e-05)
    assert_close(viscosity_converter(40, 'scott', 'kinematic viscosity'), 6.4e-05)
    assert_close(viscosity_converter(30, 'westinghouse', 'kinematic viscosity'), 0.000102)
    assert_close(viscosity_converter(60, 'westinghouse', 'kinematic viscosity'), 0.000204)

def test_roundtrip_viscosity_conversions():
    '''Test roundtrip conversions between viscosity scales using random points
    in valid overlapping ranges, accounting for scale differences.'''
    viscosity_converter(35, 'american can', 'kinematic viscosity')
    # Get all scales that can be tested
    tabulated_scales = set(viscosity_scales.keys()) - {'kinematic viscosity'}
    linear_scales = set(viscosity_scales_linear.keys())
    all_scales = tabulated_scales | linear_scales
    
    ranges = {}
    # For linear scales, use the minimum time and a reasonable maximum; test only a few scales as they're all the same
    linear_subset = {'caspers tin plate', 'american can', 'parlin cup #30', 'a&w crucible', 'astm 0.20'}
    for scale, (coef, min_time) in ((s, viscosity_scales_linear[s]) for s in linear_subset):
        # Use 10x minimum time as a reasonable maximum
        ranges[scale] = (min_time, min_time * 10)
    
    # For tabulated scales, use the limits from viscosity_converter_limits instead of the linear scales
    for scale in tabulated_scales:
        scale_min, scale_max, nu_min, nu_max = viscosity_converter_limits[scale]
        ranges[scale] = (scale_min, scale_max)
    
    rng = Random(0)
    pts = 3
    scales = sorted(ranges.keys())
    for i, scale1 in enumerate(scales):
        for scale2 in scales[i+1:]:
            # Convert range limits to kinematic viscosity for comparison
            scale1_nu_min = viscosity_converter(ranges[scale1][0], scale1, 'kinematic viscosity', extrapolate=True)
            scale1_nu_max = viscosity_converter(ranges[scale1][1], scale1, 'kinematic viscosity', extrapolate=True)
            scale2_nu_min = viscosity_converter(ranges[scale2][0], scale2, 'kinematic viscosity', extrapolate=True)
            scale2_nu_max = viscosity_converter(ranges[scale2][1], scale2, 'kinematic viscosity', extrapolate=True)
            
            # Find overlapping range in kinematic viscosity, leaving room for the polishing
            min_nu = max(scale1_nu_min, scale2_nu_min)*1.1
            max_nu = min(scale1_nu_max, scale2_nu_max)*0.9
            
            if min_nu >= max_nu:
                continue
            
            # Generate test points in kinematic viscosity space
            # Use log-uniform distribution since viscosity ranges are often logarithmic
            log_min = log(min_nu)
            log_max = log(max_nu)
            test_nus = [exp(rng.uniform(log_min, log_max)) for _ in range(pts)]
            
            for nu in test_nus:
                in_scale_1 = viscosity_converter(nu,'kinematic viscosity', scale1)
                in_scale_2 = viscosity_converter(in_scale_1, scale1, scale2)
                roundtrip = viscosity_converter(in_scale_2, scale2, scale1)
                assert_close(in_scale_1, roundtrip, rtol=1e-9)

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


def test_mu_Yaws():
    assert_close(mu_Yaws(300.0, -6.4406-log10(1000), 1117.6, 0.0137, -0.000015465), 0.00100666120816515, rtol=1e-14)

    d = dmu_Yaws_dT(300.0, -9.4406, 1117.6, 0.0137, -0.000015465)
    d_num = derivative(lambda T: mu_Yaws(T,  -9.4406, 1117.6, 0.0137, -0.000015465), 300.0, dx=300.0*1e-7)
    assert_close(d, d_num)


    T, A, B, C, D = 300.0, -9.4406, 1117.6, 0.0137, -0.000015465
    der_num = [derivative(lambda A: mu_Yaws(T, A, B, C, D), A, dx=A*1e-5),
         derivative(lambda B: mu_Yaws(T, A, B, C, D), B, dx=B*1e-5),
         derivative(lambda C: mu_Yaws(T, A, B, C, D), C, dx=C*1e-5),
         derivative(lambda D: mu_Yaws(T, A, B, C, D), D, dx=D*1e-5)]

    der_expect = [[0.0023179230916164487, 7.726410305388163e-06, 0.6953769274849346, 208.61307824548038]]
    der_analytical = mu_Yaws_fitting_jacobian([T], A, B, C, D)
    assert_close1d(der_expect, der_analytical, rtol=1e-13)
    assert_close1d(der_analytical, [der_num], rtol=1e-7)

    # Point where overflow would occur
    kwargs = {'T': 489.2, 'A': 1.0, 'B': 1.0, 'C': 1.0, 'D': 1.0}
    mu_Yaws(**kwargs)


def test_mu_TDE():
    assert_close(mu_TDE(400.0, -14.0878, 3500.26, -678132.0, 6.17706e7), 0.00018221752814389364, rtol=1e-12)

def test_PPDS5():
    assert_close(PPDS5(T=350.0, Tc=470.008, a0=1.08003e-5, a1=0.19583, a2=0.811897), 8.096643275836458e-06, rtol=1e-13)
