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

import numpy as np
import pytest
from fluids.numerics import assert_close, assert_close1d
from chemicals.heat_capacity import *
from chemicals.heat_capacity import TRC_gas_data, CRC_standard_data, Cp_data_Poling

def test_heat_capacity_CSP():
    # Example is for cis-2-butene at 350K from Poling. It is not consistent with
    # the example presented. The error is in the main expressesion.
    Cp1 = Rowlinson_Poling(350.0, 435.5, 0.203, 91.21)
    Cp2 = Rowlinson_Poling(373.28, 535.55, 0.323, 119.342)
    assert_close1d([Cp1, Cp2], [143.80194441498296, 177.62600344957252])

    # Example in VDI Heat Atlas, Example 11 in VDI D1 p 139.
    # Also formerly in ChemSep. Also in COCO.
    Cp = Rowlinson_Bondi(373.28, 535.55, 0.323, 119.342)
    assert_close(Cp, 175.39760730048116)


def test_solid_models():
    # Point at 300K from fig 5 for polyethylene, as read from the pixels of the
    # graph. 1684 was the best the pixels could give.
    Cp = Lastovka_solid(300, 0.2139)
    assert_close(Cp, 1682.063629222013)
    
    dH = Lastovka_solid_integral(300, 0.2139)
    assert_close(dH, 283246.1242170376)
    
    dS = Lastovka_solid_integral_over_T(300, 0.2139)
    assert_close(dS, 1947.553552666818)

    Cp = Dadgostar_Shaw(355.6, 0.139)
    assert_close(Cp, 1802.5291501191516)

    # Data in article:
    alphas = [0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.154, 0.154, 0.154, 0.154, 0.154, 0.135, 0.135, 0.135, 0.135, 0.135, 0.144, 0.144, 0.144, 0.144]
    Ts = [196.42, 208.55, 220.61, 229.57, 238.45, 250.18, 261.78, 270.39, 276.2, 281.74, 287.37, 295.68, 308.45, 317.91, 330.38, 336.54, 342.66, 223.2, 227.5, 244.5, 275, 278.2, 283.3, 289.4, 295, 318.15, 333.15, 348.15, 363.15, 373.15, 246.73, 249.91, 257.54, 276.24, 298.54, 495, 497, 498.15, 500, 502, 401.61, 404.6, 407.59, 410.57]
    expecteds = [1775, 1821, 1867, 1901, 1934, 1978, 2022, 2054, 2076, 2097, 2118, 2149, 2196, 2231, 2277, 2299, 2322, 1849, 1866, 1932, 2050, 2062, 2082, 2105, 2126, 2213, 2269, 2325, 2380, 2416, 1485, 1500, 1535, 1618, 1711, 2114, 2117, 2119, 2122, 2126, 1995, 2003, 2012, 2020]

    # I guess 5 J isn't bad! Could try recalculating alphas as well.
    Calculated = [Dadgostar_Shaw(T, alpha) for T, alpha in zip(Ts, alphas)]
    assert_close(expecteds, Calculated, atol=5)


def test_Zabransky_quasi_polynomial():
    Cp = Zabransky_quasi_polynomial(330, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(Cp, 165.4728226923247)

    H2 = Zabransky_quasi_polynomial_integral(300, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    H1 = Zabransky_quasi_polynomial_integral(200, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(H2 - H1, 14662.026406892926)
    
    S2 = Zabransky_quasi_polynomial_integral_over_T(300, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    S1 = Zabransky_quasi_polynomial_integral_over_T(200, 591.79, -3.12743, 0.0857315, 13.7282, 1.28971, 6.42297, 4.10989)
    assert_close(S2-S1, 59.169972919442074) # result from quadrature, not the integral


def test_Zabransky_cubic():
    Cp = Zabransky_cubic(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    assert_close(Cp, 75.31462591538555)
    
    H0 = Zabransky_cubic_integral(298.15, 20.9634, -10.1344, 2.8253, -0.256738)
    assert_close(H0, 31051.679845520584)
    
    S0 = Zabransky_cubic_integral_over_T(298.15, 20.9634, -10.1344, 2.8253,  -0.256738)
    assert_close(S0, 24.73245695987246)
    

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


def test_CRC_standard_data():
    tots_calc = [CRC_standard_data[i].abs().sum() for i in [u'Hfs', u'Gfs', u'S0s', u'Cps', u'Hfl', u'Gfl', u'S0l', 'Cpl', u'Hfg', u'Gfg', u'S0g', u'Cpg']]
    tots = [628580900.0, 306298700.0, 68541.800000000003, 56554.400000000001, 265782700.0, 23685900.0, 61274.0, 88464.399999999994, 392946600.0, 121270700.0, 141558.29999999999, 33903.300000000003]
    assert_close1d(tots_calc, tots)

    assert CRC_standard_data.index.is_unique
    assert CRC_standard_data.shape == (2470, 13)


def test_Cp_data_Poling():
    tots_calc = [Cp_data_Poling[i].abs().sum() for i in ['Tmin', 'Tmax', 'a0', 'a1', 'a2', 'a3', 'a4', 'Cpg', 'Cpl']]
    tots = [40966.0, 301000.0, 1394.7919999999999, 10.312580799999999, 0.024578948000000003, 3.1149672999999997e-05, 1.2539125599999999e-08, 43530.690000000002, 50002.459999999999]
    assert_close1d(tots_calc, tots)


    assert Cp_data_Poling.index.is_unique
    assert Cp_data_Poling.shape == (368, 10)


def test_TRC_gas_data():
    tots_calc = [TRC_gas_data[i].abs().sum() for i in ['Tmin', 'Tmax', 'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'I', 'J', 'Hfg']]
    tots = [365114, 3142000, 7794.1999999999998, 24465781000, 1056180, 133537.068, 67639.309000000008, 156121050000, 387884, 212320, 967467.89999999991, 30371.91, 495689880.0]
    assert_close1d(tots_calc, tots)

    assert TRC_gas_data.index.is_unique
    assert TRC_gas_data.shape == (1961, 14)

def test_TRC_gas():
    Cps = [TRCCp(T, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201.) for T in [150, 300]]
    assert_close1d(Cps, [35.584319834110346, 42.06525682312236])


def test_TRCCp_integral():
    dH = TRCCp_integral(298.15, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201., 1.2)
    assert_close(dH, 10802.532600592816)
    dH = TRCCp_integral(150, 4.0, 7.65E5, 720., 3.565, -0.052, -1.55E6, 52., 201., 1.2)
    assert_close(dH, 5071.355751575949)


def test_TRCCp_integral_over_T():
    coeffs = [4.0, 124000, 245, 50.538999999999994, -49.468999999999994, 220440000, 560, 78]
    dS = TRCCp_integral_over_T(300, *coeffs) - TRCCp_integral_over_T(200, *coeffs)
    assert_close(dS, 23.44278146529652)
    
def test_Zabransky_dicts():
    from chemicals.heat_capacity import (zabransky_dict_sat_s, 
                                         zabransky_dict_sat_p, 
                                         zabransky_dict_const_s, 
                                         zabransky_dict_const_p, 
                                         zabransky_dict_iso_s, 
                                         zabransky_dict_iso_p)
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
            [224963.9000000004, 293247.29999999993, 146377687.7182099]] 
    get_Tmin_sum = lambda z_splines: sum_([i.Tmin for i in z_splines])
    get_Tmax_sum = lambda z_splines: sum_([i.Tmax for i in z_splines])
    get_coeffs_sum = lambda z_splines: sum_([sum_([abs_(j) for j in i.coeffs]) for i in z_splines])
    for values, spline_dict in zip(sums, spline_dicts):
        spline_tuple = tuple(spline_dict.values())
        Tmin_sums = sum_([get_Tmin_sum(i) for i in spline_tuple])
        Tmax_sums = sum_([get_Tmax_sum(i) for i in spline_tuple])
        coeffs_sums = sum_([get_coeffs_sum(i) for i in spline_tuple])
        assert_close1d([Tmin_sums, Tmax_sums, coeffs_sums], values)
