"""Chemical Engineering Design Library (ChEDL).

Utilities for process modeling. Copyright (C) 2016, Caleb Bell
<Caleb.Andrew.Bell@gmail.com> Permission is hereby granted, free of charge, to
any person obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including without
limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions: The
above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS
IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from math import log, log10

import numpy as np
import pytest
from fluids.constants import R, h, k
from fluids.numerics import assert_close, assert_close1d, linspace, logspace

from chemicals.heat_capacity import (
    PPDS2,
    PPDS15,
    Cp_data_Poling,
    Cpg_statistical_mechanics,
    Cpg_statistical_mechanics_integral,
    Cpg_statistical_mechanics_integral_over_T,
    CRC_standard_data,
    Dadgostar_Shaw,
    Lastovka_Shaw,
    Lastovka_Shaw_integral,
    Lastovka_Shaw_integral_over_T,
    Lastovka_Shaw_T_for_Hm,
    Lastovka_Shaw_T_for_Sm,
    Lastovka_solid,
    Lastovka_solid_integral,
    Lastovka_solid_integral_over_T,
    PiecewiseHeatCapacity,
    Rowlinson_Bondi,
    Rowlinson_Poling,
    Shomate,
    Shomate_integral,
    Shomate_integral_over_T,
    TDE_CSExpansion,
    TRC_gas_data,
    TRCCp,
    TRCCp_integral,
    TRCCp_integral_over_T,
    Zabransky_cubic,
    Zabransky_cubic_integral,
    Zabransky_cubic_integral_over_T,
    Zabransky_quasi_polynomial,
    Zabransky_quasi_polynomial_integral,
    Zabransky_quasi_polynomial_integral_over_T,
    ZabranskyQuasipolynomial,
    ZabranskySpline,
    vibration_frequency_cm_to_characteristic_temperature,
    zabransky_dict_const_s,
    zabransky_dict_iso_s,
    zabransky_dict_sat_s,
)


def test_heat_capacity_CSP():
    # Example is for cis-2-butene at 350K from Poling. It is not consistent with
    # the example presented. The error is in the main expressesion.
    Cp1 = Rowlinson_Poling(350.0, 435.5, 0.203, 91.21)
    assert_close(Cp1, 143.80196224081436, rtol=1e-12)
    Cp2 = Rowlinson_Poling(373.28, 535.55, 0.323, 119.342)
    assert_close(Cp2, 177.6260232047048)

    # Example in VDI Heat Atlas, Example 11 in VDI D1 p 139.
    # Also formerly in ChemSep. Also in COCO.
    Cp = Rowlinson_Bondi(373.28, 535.55, 0.323, 119.342)
    assert_close(Cp, 175.3976263003074)


def test_solid_models():
    # Point at 300K from fig 5 for polyethylene, as read from the pixels of the
    # graph. 1684 was the best the pixels could give.
    Cp = Lastovka_solid(300, 0.2139)
    assert_close(Cp, 1682.063629222013)

    dH = Lastovka_solid_integral(300, 0.2139)
    assert_close(dH, 283246.1242170376)

    dS = Lastovka_solid_integral_over_T(300, 0.2139)
    assert_close(dS, 1947.5537561495557)

    Cp = Dadgostar_Shaw(355.6, 0.139)
    assert_close(Cp, 1802.5291501191516)

    # Data in article:
    alphas = [0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.154, 0.154, 0.154, 0.154, 0.154, 0.135, 0.135, 0.135, 0.135, 0.135, 0.144, 0.144, 0.144, 0.144]
    Ts = [196.42, 208.55, 220.61, 229.57, 238.45, 250.18, 261.78, 270.39, 276.2, 281.74, 287.37, 295.68, 308.45, 317.91, 330.38, 336.54, 342.66, 223.2, 227.5, 244.5, 275, 278.2, 283.3, 289.4, 295, 318.15, 333.15, 348.15, 363.15, 373.15, 246.73, 249.91, 257.54, 276.24, 298.54, 495, 497, 498.15, 500, 502, 401.61, 404.6, 407.59, 410.57]
    expecteds = [1775, 1821, 1867, 1901, 1934, 1978, 2022, 2054, 2076, 2097, 2118, 2149, 2196, 2231, 2277, 2299, 2322, 1849, 1866, 1932, 2050, 2062, 2082, 2105, 2126, 2213, 2269, 2325, 2380, 2416, 1485, 1500, 1535, 1618, 1711, 2114, 2117, 2119, 2122, 2126, 1995, 2003, 2012, 2020]

    # I guess 5 J isn't bad! Could try recalculating alphas as well.
    Calculated = [Dadgostar_Shaw(T, alpha) for T, alpha in zip(Ts, alphas)]
    assert_close1d(expecteds, Calculated, atol=5)


def test_Zabransky_quasi_polynomial():
    Cp = Zabransky_quasi_polynomial(330, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(Cp, 165.472878778683)

    H2 = Zabransky_quasi_polynomial_integral(300, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    H1 = Zabransky_quasi_polynomial_integral(200, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(H2 - H1, 14662.031376528757)

    S2 = Zabransky_quasi_polynomial_integral_over_T(300, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    S1 = Zabransky_quasi_polynomial_integral_over_T(200, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(S2-S1, 59.16999297436473) # result from quadrature, not the integral

def test_Zabransky_cubic():
    Cp = Zabransky_cubic(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    assert_close(Cp, 75.31465144297991)

    H0 = Zabransky_cubic_integral(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    assert_close(H0, 31051.690370364562)

    S0 = Zabransky_cubic_integral_over_T(298.15, 20.9634, -10.1344, 2.8253,  -0.256738)
    assert_close(S0, 24.732465342840854)

def test_Shomate_single():
    from scipy.integrate import quad
    water_low_gas_coeffs = [30.09200, 6.832514/1e3, 6.793435/1e6, -2.534480/1e9, 0.082139*1e6]
    assert_close(Shomate(500, *water_low_gas_coeffs), 35.21836175, rtol=1e-12)

    assert_close(Shomate_integral(500, *water_low_gas_coeffs), 15979.244791666666, rtol=1e-12)



    assert_close(Shomate_integral_over_T(500, *water_low_gas_coeffs), 191.00554193938726, rtol=1e-12)

    over_T_test_expect = Shomate_integral_over_T(500, *water_low_gas_coeffs) - Shomate_integral_over_T(300, *water_low_gas_coeffs)
    over_T_calc = quad(lambda T: Shomate(T, *water_low_gas_coeffs)/T, 300, 500)[0]
    assert_close(over_T_test_expect, over_T_calc)

    enthalpy_calc = quad(lambda T: Shomate(T, *water_low_gas_coeffs), 400, 500)[0]
    enthalpy_expect = Shomate_integral(500, *water_low_gas_coeffs) - Shomate_integral(400, *water_low_gas_coeffs)

    assert_close(enthalpy_calc, enthalpy_expect, rtol=1e-7)

def test_Lastovka_Shaw():
    # C64H52S2 (M = 885.2 alpha = 0.1333 mol/g
    # From figure 22, part b
    # Some examples didn't match so well; were the coefficients rounded?
    assert_close(Lastovka_Shaw(1000.0, 0.1333), 2467.113309084757)

    # Same, but try the correlation for cyclic aliphatic compounds
    assert_close(Lastovka_Shaw(1000.0, 0.1333, cyclic_aliphatic=True), 2187.36187944884)


def test_Lastovka_Shaw_integral():
    # C64H52S2 (M = 885.2 alpha = 0.1333 mol/g
    # From figure 22, part b
    # Some examples didn't match so well; were the coefficients rounded?
    assert_close(Lastovka_Shaw_integral(1000.0, 0.1333), 6615282.290516732)

    # Same, but try the correlation for cyclic aliphatic compounds
    assert_close(Lastovka_Shaw_integral(1000.0, 0.1333, cyclic_aliphatic=True), 6335530.860880815)

def test_Lastovka_Shaw_integral_over_T():
    dS = Lastovka_Shaw_integral_over_T(300.0, 0.1333)
    assert_close(dS, 3609.791928945323)

    dS = Lastovka_Shaw_integral_over_T(1000.0, 0.1333, cyclic_aliphatic=True)
    assert_close(dS, 3790.4489380423597)

def test_Lastovka_Shaw_T_for_Hm():
    sv, MW, Hm = 0.05, 12.0, 3727593.7203149837
    T = Lastovka_Shaw_T_for_Hm(Hm=Hm, MW=MW, similarity_variable=sv)
    assert_close(T, 153731.155415042)

    T = Lastovka_Shaw_T_for_Hm(Hm=55000, MW=80.0, similarity_variable=0.23)
    assert_close(T, 600.0943429567604)

    T = Lastovka_Shaw_T_for_Hm(Hm=55000, MW=80.0, similarity_variable=0.23, factor=2)
    assert_close(T, 469.688546084446)

def test_Lastovka_Shaw_T_for_Sm():
    T = Lastovka_Shaw_T_for_Sm(Sm=112.80, MW=72.151, similarity_variable=0.2356)
    assert_close(T, 603.4298291570275)

    T = Lastovka_Shaw_T_for_Sm(Sm=112.80, MW=72.151, similarity_variable=0.2356, factor=2)
    assert_close(T, 446.62661557878624)

@pytest.mark.slow
@pytest.mark.fuzz
def test_Lastovka_Shaw_T_for_Hm_fuzz():
    T_ref = 298.15
    factor = 1.0

    # Originally tested with 50 points in everything
    similarity_variables = linspace(.05, .5, 8)
    MWs = linspace(12, 1200, 8)
    Hms = [-i for i in logspace(log10(1e4), log10(1), 15)] + logspace(log10(1), log10(1e7), 15)

    for sv in similarity_variables:
        for MW in MWs:
            for Hm in Hms:
                try:
                    Lastovka_Shaw_T_for_Hm(Hm=Hm, MW=MW, similarity_variable=sv, T_ref=T_ref)
                except Exception as e:
                    if 'negative temperature' in str(e):
                        continue
    #                 print(sv, MW, Hm, e)


# This test used to check that the only error raised was NotBoundedError
# However, to ad dnumba compatibility the exception could not be caught as e.
#@pytest.mark.slow
#@pytest.mark.fuzz
#def test_Lastovka_Shaw_T_for_Sm_fuzz():
#    T_ref = 298.15
#    factor = 1.0
#
#    similarity_variables = linspace(.05, .5, 8)
#    MWs = linspace(12, 1200, 8)
#    Sms = logspace(log10(3000), log10(300), 15)
#
#    for sv in similarity_variables:
#        for MW in MWs:
#            for Sm in Sms:
#                try:
#                    Lastovka_Shaw_T_for_Sm(Sm=Sm, MW=MW, similarity_variable=sv, T_ref=T_ref)
#                except Exception as e:
#                    if 'negative temperature' in str(e):
#                        continue
#                    elif isinstance(e, NotBoundedError):
#                        continue
#                    else:
#                        raise ValueError("Could not converge with unexpected error")

def test_CRC_standard_data():
    tots_calc = [CRC_standard_data[i].abs().sum() for i in ['Hfs', 'Gfs', 'S0s', 'Cps', 'Hfl', 'Gfl', 'S0l', 'Cpl', 'Hfg', 'Gfg', 'S0g', 'Cpg']]
    tots = [628580900.0, 306298700.0, 68541.800000000003, 56554.400000000001, 265782700.0, 23685900.0, 61274.0, 88464.399999999994, 392946600.0, 121270700.0, 141558.29999999999, 33903.300000000003]
    assert_close1d(tots_calc, tots)

    assert CRC_standard_data.shape == (2470, 13)


def test_Cp_data_Poling():
    tots_calc = [Cp_data_Poling[i].abs().sum() for i in ['Tmin', 'Tmax', 'a0', 'a1', 'a2', 'a3', 'a4', 'Cpg', 'Cpl']]
    tots = [40966.0, 301000.0, 1394.792, 10.3125808, 0.024578948000000003, 3.1149673000000004e-05, 1.25391256e-08, 43400.95, 49813.16]
    assert_close1d(tots_calc, tots)

    assert Cp_data_Poling.shape == (367, 10)


def test_TRC_gas_data():
    tots_calc = [TRC_gas_data[i].abs().sum() for i in ['Tmin', 'Tmax', 'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'I', 'J', 'Hfg']]
    tots = [365114, 3142000, 7794.1999999999998, 24465781000, 1056180, 133537.068, 67639.309000000008, 156121050000, 387884, 212320, 967467.89999999991, 30371.91, 495689880.0]
    assert_close1d(tots_calc, tots)

    assert TRC_gas_data.shape == (1961, 14)

def test_TRC_gas():
    Cps = [TRCCp(T, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201.) for T in [150, 300]]
    assert_close1d(Cps, [35.58433189527471, 42.06527108097466])

    # Test a zero division error
    Cp_zero = TRCCp(298, 4.0, 32753000, 892, 76.209, -0.254, -58870000, 98, 298)
    Cp_almost_zero = TRCCp(298*(1+1e-13), 4.0, 32753000, 892, 76.209, -0.254, -58870000, 98, 298)
    assert_close(Cp_zero, Cp_almost_zero)


def test_TRCCp_integral():
    dH = TRCCp_integral(298.15, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201., 1.2)
    assert_close(dH, 10802.536262068483)
    dH = TRCCp_integral(150, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201., 1.2)
    assert_close(dH, 5071.357470491908)


def test_TRCCp_integral_over_T():
    coeffs = [4.0, 124000, 245, 50.538999999999994, -49.468999999999994, 220440000, 560, 78]
    dS = TRCCp_integral_over_T(300, *coeffs) - TRCCp_integral_over_T(200, *coeffs)
    assert_close(dS, 23.4427894111345)

def test_Zabransky_dicts():
    from chemicals.heat_capacity import (
        zabransky_dict_const_p,
        zabransky_dict_const_s,
        zabransky_dict_iso_p,
        zabransky_dict_iso_s,
        zabransky_dict_sat_p,
        zabransky_dict_sat_s,
    )
    quasi_dicts = [zabransky_dict_sat_p, zabransky_dict_const_p, zabransky_dict_iso_p]
    spline_dicts = [zabransky_dict_sat_s, zabransky_dict_const_s, zabransky_dict_iso_s]

    sums = [[4811.400000000001, 7889.5, 11323.099999999999],
            [6724.9, 11083.200000000004, 17140.75],
            [37003.999999999985, 72467.1, 91646.64999999997]]
    coeff_sums = [[1509.3503910000002, 170.63668272100003, 553.71843, 3602.9634764, 2731.1423000000004, 2505.7230399999994],
                  [394736.92543421406, 52342.25656440046, 54451.52735, 366067.89141800004, 61161.632348850006, 207335.59372000452],
                  [85568.9366422, 9692.534972905993, 13110.905983999992, 97564.75442449997, 30855.65738500001, 73289.607074896]]
    attrs = ['Tmin', 'Tmax', 'Tc']
    getattr_ = getattr
    sum_ = sum
    abs_ = abs
    for i in range(len(quasi_dicts)):
        tot_calc = [sum_([abs_(getattr_(k, j)) for k in quasi_dicts[i].values()]) for j in attrs]
        assert_close1d(tot_calc, sums[i])
        coeff_calc = sum_(np.abs([getattr_(k, 'coeffs') for k in quasi_dicts[i].values()]))
        assert_close1d(coeff_calc, coeff_sums[i])

    sums = [[108141.20000000004, 153049.59999999992, 16980025.79560638],
            [15739.999999999998, 22191.1, 60581.66851630001],
            [223632.70000000065, 291814.8999999999, 146374405.8581997]]
    get_Tmin_sum = lambda z_splines: sum_([i.Tmin for i in z_splines])
    get_Tmax_sum = lambda z_splines: sum_([i.Tmax for i in z_splines])
    get_coeffs_sum = lambda z_splines: sum_([sum_([abs_(j) for j in i.coeffs]) for i in z_splines])
    for values, spline_dict in zip(sums, spline_dicts):
        spline_tuple = tuple(spline_dict.values())
        Tmin_sums = sum_([get_Tmin_sum(i) for i in spline_tuple])
        Tmax_sums = sum_([get_Tmax_sum(i) for i in spline_tuple])
        coeffs_sums = sum_([get_coeffs_sum(i) for i in spline_tuple])
        assert_close1d([Tmin_sums, Tmax_sums, coeffs_sums], values)

def test_HeatCapacityClass_methods():
    classes = [ZabranskySpline, ZabranskyQuasipolynomial]
    for c in classes:
        assert hasattr(c, 'calculate')
        assert hasattr(c, 'calculate_integral')
        assert hasattr(c, 'calculate_integral_over_T')

def test_ZABRANSKY_SPLINE():
    from chemicals.heat_capacity import zabransky_dict_iso_s
    d = zabransky_dict_iso_s['7732-18-5']
    assert_close(d.force_calculate(645), 4521.767801978968)
    assert_close(d.calculate(400), 76.53795418514787)
    assert_close(d.force_calculate(250), 77.1082743553618)
#    assert_close1d(d.Ts, [273.6, 380.0, 590.0, 635.0, 644.6])

    # Test enthalpy integrals

    assert_close(d.force_calculate_integral(200, 270), 5505.7344456789615)
    assert_close(d.force_calculate_integral(200, 280), 6264.127552384685)
    assert_close(d.force_calculate_integral(200, 380), 13817.514452840238)
    assert_close(d.force_calculate_integral(200, 380+1E-9), 13817.514452916399)
    assert_close(d.force_calculate_integral(300, 590-1E-9), 24309.137141163086)
    assert_close(d.force_calculate_integral(300, 590+1E-9), 24309.13714137674)
    assert_close(d.force_calculate_integral(200, 635-1E-9), 39698.85970996421)
    assert_close(d.force_calculate_integral(200, 635+1E-9), 39698.85970642518)
    assert_close(d.force_calculate_integral(200, 644.6), 76304.80397369813)
    assert_close(d.force_calculate_integral(200, 645), 78093.69002487611)

    # Same test cases, flipped around
    assert_close(d.force_calculate_integral(270, 200), -5505.7344456789615)
    assert_close(d.force_calculate_integral(280, 200), -6264.127552384685)
    assert_close(d.force_calculate_integral(380, 200), -13817.514452840238)
    assert_close(d.force_calculate_integral(380+1E-9, 200), -13817.514452916399)
    assert_close(d.force_calculate_integral(590-1E-9, 300), -24309.137141163086)
    assert_close(d.force_calculate_integral(590+1E-9, 300), -24309.13714137674)
    assert_close(d.force_calculate_integral(635-1E-9, 200), -39698.85970996421)
    assert_close(d.force_calculate_integral(635+1E-9, 200), -39698.85970642518)
    assert_close(d.force_calculate_integral(644.6, 200), -76304.80397369813)
    assert_close(d.force_calculate_integral(645, 200), -78093.69002487611)

    # Test entropy integrals
    assert_close(d.force_calculate_integral_over_T(200, 270), 23.65404552882864)
    assert_close(d.force_calculate_integral_over_T(200, 280), 26.412181131648946)
    assert_close(d.force_calculate_integral_over_T(200, 380), 49.473640416141926)
    assert_close(d.force_calculate_integral_over_T(200, 380+1E-9), 49.473640416342356)
    assert_close(d.force_calculate_integral_over_T(300, 590-1E-9), 55.558387640582296)
    assert_close(d.force_calculate_integral_over_T(300, 590+1E-9), 55.558387640900335)
    assert_close(d.force_calculate_integral_over_T(200, 635-1E-9), 99.5436500019074)
    assert_close(d.force_calculate_integral_over_T(200, 635+1E-9), 99.5436500025113)
    assert_close(d.force_calculate_integral_over_T(200, 644.6), 156.74432313471482)
    assert_close(d.force_calculate_integral_over_T(200, 645), 159.51864708803043)

    # Same test cases, flipped around
    assert_close(d.force_calculate_integral_over_T(270, 200), -23.65404552882864)
    assert_close(d.force_calculate_integral_over_T(280, 200), -26.412181131648946)
    assert_close(d.force_calculate_integral_over_T(380, 200), -49.473640416141926)
    assert_close(d.force_calculate_integral_over_T(380+1E-9, 200), -49.473640416342356)
    assert_close(d.force_calculate_integral_over_T(590-1E-9, 300), -55.558387640582296)
    assert_close(d.force_calculate_integral_over_T(590+1E-9, 300), -55.558387640900335)
    assert_close(d.force_calculate_integral_over_T(635-1E-9, 200), -99.5436500019074)
    assert_close(d.force_calculate_integral_over_T(635+1E-9, 200), -99.5436500025113)
    assert_close(d.force_calculate_integral_over_T(644.6, 200), -156.74432313471482)
    assert_close(d.force_calculate_integral_over_T(645, 200), -159.51864708803043)


    # Test a chemical with only one set of coefficients
    d = zabransky_dict_iso_s['2016-57-1']
    assert_close(d.calculate(310), 375.543177681642)
    assert_close(d.force_calculate_integral(290, 340), 18857.29436766774)
    assert_close(d.force_calculate_integral_over_T(290, 340), 59.96511735461314)


def test_zabransky_dict_types():
    for d in (zabransky_dict_sat_s, zabransky_dict_iso_s, zabransky_dict_const_s):
        for v in d.values():
            assert type(v) is PiecewiseHeatCapacity


def test_PPDS2():
    Cp = PPDS2(T=350.0, Ts=462.493, C_low=4.54115, C_inf=9.96847, a1=-103.419, a2=695.484, a3=-2006.1, a4=2476.84, a5=-1186.47)
    assert_close(Cp, 136.46338956689826, rtol=1e-13)

def test_PPDS15():
    Cp = PPDS15(T=400.0, Tc=562.05, a0=0.198892, a1=24.1389, a2=-20.2301, a3=5.72481, a4=4.43613e-7, a5=-3.10751e-7)
    assert_close(Cp, 161.89831435090065, rtol=1e-14)

def test_TDE_CSExpansion():
    Cp = TDE_CSExpansion(550.0, 778.0, 0.626549, 120.705, 0.255987, 0.000381027, -3.03077e-7)
    assert_close(Cp, 328.4720426864035, rtol=1e-14)

def test_piece_wise_heat_capacity():
    class HeatCapacity:
        __slots__ = ('value', 'Tmin', 'Tmax')
        def __init__(self, value, Tmin, Tmax):
            self.value = value
            self.Tmin = Tmin
            self.Tmax = Tmax

        def calculate_integral(self, Ta, Tb):
            return self.value * (Tb - Ta)

        def calculate_integral_over_T(self, Ta, Tb):
            return self.value * log(Tb/Ta)

    models = [HeatCapacity(1, 200, 250), HeatCapacity(1, 250, 300), HeatCapacity(1, 300, 350)]

    Cp = PiecewiseHeatCapacity(models)
    # Trivial
    assert_close(0., Cp.force_calculate_integral(298.15, 298.15))
    assert_close(0., Cp.force_calculate_integral_over_T(298.15, 298.15))

    # Within bounds
    H_TminTmax = sum([i.calculate_integral(i.Tmin, i.Tmax) for i in models])
    S_TminTmax = sum([i.calculate_integral_over_T(i.Tmin, i.Tmax) for i in models])
    assert_close(H_TminTmax, Cp.force_calculate_integral(Cp.Tmin, Cp.Tmax))
    assert_close(S_TminTmax, Cp.force_calculate_integral_over_T(Cp.Tmin, Cp.Tmax))

    # Left to center
    assert_close(models[0].calculate_integral(100, 250), Cp.force_calculate_integral(100, 250))
    assert_close(models[0].calculate_integral_over_T(100, 250), Cp.force_calculate_integral_over_T(100, 250))

    # Center to right
    assert_close(models[-1].calculate_integral(300, 400), Cp.force_calculate_integral(300, 400))
    assert_close(models[-1].calculate_integral_over_T(300, 400), Cp.force_calculate_integral_over_T(300, 400))

    # Across both sides
    H = (models[0].calculate_integral(100, 200)
         + H_TminTmax
         +models[-1].calculate_integral(350, 400))
    S = (models[0].calculate_integral_over_T(100, 200)
         + S_TminTmax
         +models[-1].calculate_integral_over_T(350, 400))
    assert_close(H, Cp.force_calculate_integral(100, 400))
    assert_close(S, Cp.force_calculate_integral_over_T(100, 400))


def test_Cpg_statistical_mechanics():
    thetas = [1360, 2330, 2330, 4800, 4880, 4880]
    Cp_from_spectra = Cpg_statistical_mechanics(300.0,thetas)
    assert_close(Cp_from_spectra, 35.55983440173097, rtol=1e-12)

    test = Cpg_statistical_mechanics(300.0,thetas, linear=True)
    assert_close(test, 35.55983440173097-R/2, rtol=1e-12)

    # Test the high limit
    assert_close(Cpg_statistical_mechanics(1e100, thetas)/R, 10, rtol=1e-14)

    # Test there is a check to ensure at high temperatures numerical error does not make the value too high

    assert_close(Cpg_statistical_mechanics(1e10, thetas)/R, 9.99999999999993, rtol=1e-14)

    #Test the use of expm1 at high temperatures
    assert_close(Cpg_statistical_mechanics(1e11, thetas)/R, 10, rtol=1e-15)

    # Low temperatures

    assert 4 == Cpg_statistical_mechanics(1, thetas)/R
    assert 4 == Cpg_statistical_mechanics(10, thetas)/R
    assert 4 == Cpg_statistical_mechanics(0, thetas)/R
    assert 4 == Cpg_statistical_mechanics(1e-10, thetas)/R

    # Other Perry's sample from software as Hz
    v_scaled = [3.24, 4.97, 4.97, 9.90, 10.26, 10.26] # divided by 10^13,
    v_scaled_Hz = [vj*1e13 for vj in v_scaled]
    thetas_comp = [h*vj/k for vj in v_scaled_Hz]

    Cp_perry_stat = Cpg_statistical_mechanics(300.0,thetas_comp)
    assert_close(Cp_perry_stat, 34.89647856513431, rtol=1e-12)


    thetas_caleb_psi4_mp2_631G = [1615.6879, 2486.5201, 2486.6163, 5128.4685, 5353.9398, 5354.6923]
    thetas_caleb_psi4_mp2_6311ppG3df3pd = [1472.1887, 2385.2593, 2385.3405, 5084.3735, 5296.0091, 5296.3665]

    # Yep, it's matching
    Cp = Cpg_statistical_mechanics(298.15,thetas_caleb_psi4_mp2_631G)
    assert_close(Cp, 34.626, atol=1e-3)


def test_Cpg_statistical_mechanics_integral():
    thetas = [1360, 2330, 2330, 4800, 4880, 4880]

    # numerical = quad(Cpg_statistical_mechanics, 399, 400, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(400, thetas)-Cpg_statistical_mechanics_integral(399, thetas)
    assert_close(38.371249200803504, analytical, rtol=1e-13)


    # numerical = quad(Cpg_statistical_mechanics, 1, 400, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(400, thetas)-Cpg_statistical_mechanics_integral(1, thetas)
    assert_close(13775.685078010132, analytical, rtol=1e-13)

    # numerical = quad(Cpg_statistical_mechanics, 0, 400, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(400, thetas)-Cpg_statistical_mechanics_integral(0, thetas)
    assert_close(13808.942928482746, analytical, rtol=1e-13)

    # numerical = quad(Cpg_statistical_mechanics, 0, 1e-10, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(1e-10, thetas)-Cpg_statistical_mechanics_integral(0, thetas)
    assert_close(3.325785047261296e-09, analytical, rtol=1e-13)

    # numerical = quad(Cpg_statistical_mechanics, 1000, 2000, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(2000, thetas)-Cpg_statistical_mechanics_integral(1000, thetas)
    assert_close(65130.40200286856, analytical, rtol=1e-13)

    # numerical = quad(Cpg_statistical_mechanics, 1e9, 2e9, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(2e9, thetas)-Cpg_statistical_mechanics_integral(1e9, thetas)
    assert_close(83144626181.50354, analytical, rtol=1e-13)

    val = Cpg_statistical_mechanics_integral(1e100, thetas)
    assert val > 0

    # numerical = quad(Cpg_statistical_mechanics, 1e100, 2e100, args=(thetas,))[0]
    analytical = Cpg_statistical_mechanics_integral(2e100, thetas)-Cpg_statistical_mechanics_integral(1e100, thetas)
    assert_close(8.31446261815324e+101, analytical, rtol=1e-13)

def test_Cpg_statistical_mechanics_integral_over_T():
    thetas = [1360, 2330, 2330, 4800, 4880, 4880]

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 399, 400)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(400, thetas)-Cpg_statistical_mechanics_integral_over_T(399, thetas)
    assert_close(0.09604821735794644, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1e-2, 1)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(1, thetas)-Cpg_statistical_mechanics_integral_over_T(1e-2, thetas)
    assert_close(153.15806144652714, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1e-2, 10)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(10, thetas)-Cpg_statistical_mechanics_integral_over_T(1e-2, thetas)
    assert_close(229.73709216979074, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1e-20, 1e-10)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(1e-10, thetas)-Cpg_statistical_mechanics_integral_over_T(1e-20, thetas)
    assert_close(765.7903072326358, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1000, 2000)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(2000, thetas)-Cpg_statistical_mechanics_integral_over_T(1000, thetas)
    assert_close(44.506610727671244, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1e9, 2e9)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(2e9, thetas)-Cpg_statistical_mechanics_integral_over_T(1e9, thetas)
    assert_close(57.6314632164183, analytical, rtol=1e-13)

    # numerical = quad(lambda T: Cpg_statistical_mechanics(T, thetas)/T, 1, 400)[0]
    analytical = Cpg_statistical_mechanics_integral_over_T(400, thetas)-Cpg_statistical_mechanics_integral_over_T(1, thetas)
    assert_close(200.85926489512875, analytical, rtol=1e-13)



def test_vibration_frequency_cm_to_characteristic_temperature():
    T = vibration_frequency_cm_to_characteristic_temperature(667)
    assert_close(T, 959.6641613636505, rtol=1e-13)

    NIST_RECOMMENDED_SCALE_FREQ = 0.9365

    # acetylene
    cminvs = [458.0292,726.0247,755.0628, 2002.2779,3501.5995,3591.7215]
    things = [vibration_frequency_cm_to_characteristic_temperature(f, NIST_RECOMMENDED_SCALE_FREQ) for f in cminvs]

    Cp = Cpg_statistical_mechanics(298.15, things)
    assert_close(Cp, 46.27123343052294, rtol=1e-12)

