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

from numpy.testing import assert_allclose
import pytest
import numpy as np
import pandas as pd
from fluids.numerics import assert_close, derivative, assert_close1d
from chemicals.vapor_pressure import *
from chemicals.vapor_pressure import Psat_data_WagnerMcGarry, Psat_data_AntoinePoling, Psat_data_WagnerPoling, Psat_data_AntoineExtended, Psat_data_Perrys2_8, Psat_data_VDI_PPDS_3
from chemicals.identifiers import check_CAS
from math import *

### Regression equations

def test_Wagner_original():
    # Methane, coefficients from [2]_, at 100 K
    Psat = Wagner_original(100.0, 190.53, 4596420., -6.00435, 1.1885, -0.834082, -1.22833)
    assert_allclose(Psat, 34520.44601450496)


def test_Wagner():
    # Methane, coefficients from [2]_, at 100 K.
    Psat = Wagner(100., 190.551, 4599200, -6.02242, 1.26652, -0.5707, -1.366)
    assert_allclose(Psat, 34415.00476263708)


def test_TRC_Antoine_extended():
    # Tetrafluoromethane, coefficients from [1]_, at 180 K:
    Psat = TRC_Antoine_extended(180.0, 227.51, -120., 8.95894, 510.595, -15.95, 2.41377, -93.74, 7425.9)
    assert_allclose(Psat, 706317.0898414153)

    # Test x is restricted to 0
    Psat = TRC_Antoine_extended(120.0, 227.51, -120., 8.95894, 510.595, -15.95, 2.41377, -93.74, 7425.9)
    assert_allclose(Psat, 11265.018958511126)


def test_Antoine():
    # Methane, coefficients from [1]_, at 100 K:
    Psat = Antoine(100.0, 8.7687, 395.744, -6.469)
    assert_allclose(Psat, 34478.367349639906)
    
    # Tetrafluoromethane, coefficients from [1]_, at 180 K
    Psat = Antoine(180, A=8.95894, B=510.595, C=-15.95)
    assert_allclose(Psat, 702271.0518579542)
    
    # Oxygen at 94.91 K, with coefficients from [3]_ in units of °C, mmHg, log10,
    # showing the conversion of coefficients A (mmHg to Pa) and C (°C to K)
    Psat = Antoine(94.91, 6.83706+2.1249, 339.2095, 268.70-273.15)
    assert_allclose(Psat, 162978.88655572367)

def test_Antoine_fit_extrapolate():
    T = 178.01
    A, B, C =  (24.0989474955895, 4346.793091137991, -18.96968471040141)
    Psat, dPsat_dT, d2Psat_dT2 = (0.03946094565666715, 0.006781441203850251, 0.0010801244983894853)
    
    assert_close(Psat, Antoine(T, A, B, C, base=e), rtol=1e-10)
    dPsat_dT_num = derivative(lambda T: Antoine(T, A, B, C, base=e), T, dx=T*1e-9)
    assert_close(dPsat_dT_num, dPsat_dT, rtol=1e-7)
    d2Psat_dT2_num = derivative(lambda T: Antoine(T, A, B, C, base=e), T, dx=T*1e-7, n=2, order=13)
    d2Psat_dT2_analytical = B*(B/(C + T) - 2)*exp(A - B/(C + T))/(C + T)**3
    assert_close(d2Psat_dT2_num, d2Psat_dT2, rtol=1e-4)
    assert_close(d2Psat_dT2_analytical, d2Psat_dT2, rtol=1e-10)

def test_Antoine_coeffs_from_point():
    T = 178.01
    A, B, C = (24.0989474955895, 4346.793091137991, -18.96968471040141)
    Psat = Antoine(T, A, B, C, base=exp(1))
    dPsat_dT, d2Psat_dT2 = (0.006781441203850251, 0.0010801244983894853)
    new_coeffs = Antoine_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2, base=exp(1))
    assert_close1d(new_coeffs, [A, B, C], rtol=1e-9)

def test_test_Antoine_AB():
    T = 178.01
    A, B = (27.358925161569008, 5445.569591293226)
    Psat = Antoine(T, A, B, C=0, base=e)
    dPsat_dT = B*e**(A - B/T)*log(e)/T**2
    AB = Antoine_AB_coeffs_from_point(T, Psat, dPsat_dT, base=exp(1))
    assert_close(A, AB[0])
    assert_close(B, AB[1])



def test_Antoine_AB_fit_extrapolate():
    T = 178.01
    Psat, dPsat_dT, d2Psat_dT2 = (0.03946094565666715, 0.006781441203850251, 0.0010801244983894853)
    
    A, B = (27.358925161569008, 5445.569591293226)
    C = 0.0
    assert_close(Psat, Antoine(T, A, B, C, base=e), rtol=1e-10)
    dPsat_dT_num = derivative(lambda T: Antoine(T, A, B, C, base=e), T, dx=T*1e-9)
    assert_close(dPsat_dT_num, dPsat_dT, rtol=1e-7)
    # d2Psat_dT2_num = derivative(lambda T: Antoine(T, A, B, C, base=e), T, dx=T*1e-7, n=2, order=13)
    
    d2Psat_dT2_analytical = B*(B/T - 2)*exp(A - B/T)/T**3
    
    # Second derivative does not match, but is similar - small jump
    assert_close(d2Psat_dT2, d2Psat_dT2_analytical, rtol=1e-2)

def test_DIPPR101_ABC_coeffs_from_point():
    T = 178.01
    Psat, dPsat_dT, d2Psat_dT2 = (0.03946094565666715, 0.006781441203850251, 0.0010801244983894853)
    ABC = DIPPR101_ABC_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2)
    assert_close1d(ABC, (72.47169926642722, -6744.620564969687, -7.2976291987890844))



### Data integrity
def test_WagnerMcGarry():
    sums_calc = [Psat_data_WagnerMcGarry[i].abs().sum() for i in ['A', 'B', 'C', 'D', 'Pc', 'Tc', 'Tmin']]
    sums = [1889.3027499999998, 509.57053652899992, 1098.2766456999998, 1258.0866876, 1005210819, 129293.19100000001, 68482]
    assert_allclose(sums_calc, sums)

    assert Psat_data_WagnerMcGarry.index.is_unique
    assert Psat_data_WagnerMcGarry.shape == (245, 8)
    for i in Psat_data_WagnerMcGarry.index:
        assert check_CAS(i)


def test_AntoinePoling():
    sums_calc =  [Psat_data_AntoinePoling[i].abs().sum() for i in ['A', 'B', 'C', 'Tmin', 'Tmax']]
    sums = [2959.75131, 398207.29786, 18732.24601, 86349.09, 120340.66]
    assert_allclose(sums_calc, sums)

    assert Psat_data_AntoinePoling.index.is_unique
    assert Psat_data_AntoinePoling.shape == (325, 6)
    assert all([check_CAS(i) for i in Psat_data_AntoinePoling.index])


def test_WagnerPoling():
    sums_calc =  [Psat_data_WagnerPoling[i].abs().sum() for i in ['A', 'B', 'C', 'D', 'Tmin', 'Tmax', 'Tc', 'Pc']]
    sums = [894.39071999999999, 271.76480999999995, 525.8134399999999, 538.25393000000008, 24348.006000000001, 59970.149999999994, 63016.021000000001, 357635500]
    assert_allclose(sums_calc, sums)

    assert Psat_data_WagnerPoling.index.is_unique
    assert Psat_data_WagnerPoling.shape == (104, 9)
    assert all([check_CAS(i) for i in Psat_data_WagnerPoling.index])


def test_AntoineExtended():
    sums_calc = [Psat_data_AntoineExtended[i].abs().sum() for i in ['A', 'B', 'C', 'Tc', 'to', 'n', 'E', 'F', 'Tmin', 'Tmax']]
    sums = [873.55827000000011, 107160.285, 4699.9650000000001, 47592.470000000001, 7647, 241.56537999999998, 22816.815000000002, 1646509.79, 33570.550000000003, 46510.849999999999]
    assert_allclose(sums_calc, sums)

    assert Psat_data_AntoineExtended.index.is_unique
    assert Psat_data_AntoineExtended.shape == (97, 11)
    assert all([check_CAS(i) for i in Psat_data_AntoineExtended.index])
    
def test_VDI_PPDS_3_data():
    """I believe there are no errors here.

    Average temperature deviation
    0.144% vs tabulated values.
    """
    assert all([check_CAS(i) for i in Psat_data_VDI_PPDS_3.index])
    tots_calc = [Psat_data_VDI_PPDS_3[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'Tc', u'Pc', u'Tm']]
    tots = [2171.4607300000002, 694.38631999999996, 931.3604499999999, 919.88944000000004, 150225.16000000003, 1265565000, 56957.849999999991]
    assert_allclose(tots_calc, tots)
    
    assert Psat_data_VDI_PPDS_3.index.is_unique
    assert Psat_data_VDI_PPDS_3.shape == (275, 8)

def test_Perrys2_8():
    assert all([check_CAS(i) for i in Psat_data_Perrys2_8.index])
    tots_calc = [Psat_data_Perrys2_8[i].abs().sum() for i in ['C1', 'C2', 'C3', 'C4', 'C5', 'Tmin', 'Tmax']]
    tots = [30288.457300000002, 2574584.506, 3394.9677, 1.1357374248600194, 1223.7, 70399.92, 187558.921]
    assert_allclose(tots_calc, tots)

    assert Psat_data_Perrys2_8.index.is_unique
    assert Psat_data_Perrys2_8.shape == (340, 8)

### CSP relationships
def test_boiling_critical_relation():
    P = boiling_critical_relation(347.2, 409.3, 617.1, 36E5)
    assert_allclose(P, 15209.467273093938)


def test_Lee_Kesler():
    # Example from [2]_; ethylbenzene at 347.2 K.
    # Their result is 0.132 bar.
    P = Lee_Kesler(347.2, 617.1, 36E5, 0.299)
    assert_allclose(P, 13078.694162949312)


def test_Ambrose_Walton():
    # Example from [2]_; ethylbenzene at 347.25 K.
    # Their result is 0.1329 bar.
    Psat = Ambrose_Walton(347.25, 617.15, 36.09E5, 0.304)
    assert_allclose(Psat, 13278.878504306222)

    
def test_Edalat():
    # No check data, but gives the same results as the other CSP relationships
    Psat = Edalat(347.2, 617.1, 36E5, 0.299)
    assert_allclose(Psat, 13461.273080743307)

def test_Sanjari():
    P = Sanjari(347.2, 617.1, 36E5, 0.299)
    assert_allclose(P, 13651.916109552498)

    Ts_dat = [125.45, 126.54775, 127.6455, 128.74325, 129.841, 130.93875, 132.0365, 133.13425, 134.232, 135.32975, 136.4275, 137.52525, 138.623, 139.72075, 140.8185, 141.91625, 143.014, 144.11175, 145.2095, 146.30725, 147.405, 148.50275, 149.6005, 150.69825, 151.796, 152.89375, 153.9915, 155.08925, 156.187, 157.28475, 158.3825, 159.48025, 160.578, 161.67575, 162.7735, 163.87125, 164.969, 166.06675, 167.1645, 168.26225, 169.36, 170.45775, 171.5555, 172.65325, 173.751, 174.84875, 175.9465, 177.04425, 178.142, 179.23975, 180.3375, 181.43525, 182.533, 183.63075, 184.7285, 185.82625, 186.924, 188.02175, 189.1195, 190.21725, 191.315, 192.41275, 193.5105, 194.60825, 195.706, 196.80375, 197.9015, 198.99925, 200.097, 201.19475, 202.2925, 203.39025, 204.488, 205.58575, 206.6835, 207.78125, 208.879, 209.97675, 211.0745, 212.17225, 213.27, 214.36775, 215.4655, 216.56325, 217.661, 218.75875, 219.8565, 220.95425, 222.052, 223.14975, 224.2475, 225.34525, 226.443, 227.54075, 228.6385, 229.73625, 230.834, 231.93175, 233.0295, 234.12725, 235.225, 236.32275, 237.4205, 238.51825, 239.616, 240.71375, 241.8115, 242.90925, 244.007, 245.10475, 246.2025, 247.30025, 248.398, 249.49575, 250.5935, 251.69125, 252.789, 253.88675, 254.9845, 256.08225, 257.18, 258.27775, 259.3755, 260.47325, 261.571, 262.66875, 263.7665, 264.86425, 265.962, 267.05975, 268.1575, 269.25525, 270.353, 271.45075, 272.5485, 273.64625, 274.744, 275.84175, 276.9395, 278.03725, 279.135, 280.23275, 281.3305, 282.42825, 283.526, 284.62375, 285.7215, 286.81925, 287.917, 289.01475, 290.1125, 291.21025, 292.308, 293.40575, 294.5035, 295.60125, 296.699, 297.79675, 298.8945, 299.99225, 301.09, 302.18775, 303.2855, 304.38325, 305.481, 306.57875, 307.6765, 308.77425, 309.872, 310.96975, 312.0675, 313.16525, 314.263, 315.36075, 316.4585, 317.55625, 318.654, 319.75175, 320.8495, 321.94725, 323.045, 324.14275, 325.2405, 326.33825, 327.436, 328.53375, 329.6315, 330.72925, 331.827, 332.92475, 334.0225, 335.12025, 336.218, 337.31575, 338.4135, 339.51125, 340.609, 341.70675, 342.8045, 343.90225]
    Ps_dat = [2.01857353521E-006, 0.000002517, 3.12468960653E-006, 3.86254620966E-006, 4.75480477553E-006, 5.82952636953E-006, 7.11906108993E-006, 8.66057817586E-006, 1.04966307916E-005, 1.2675775303E-005, 1.52532334366E-005, 1.82916090268E-005, 2.18616534529E-005, 2.60430842891E-005, 3.0925458133E-005, 3.66090963473E-005, 4.32060629106E-005, 5.08412011933E-005, 5.96532246339E-005, 6.97958540095E-005, 0.000081439, 9.47700963623E-005, 0.0001099952, 0.0001273406, 0.0001470541, 0.0001694061, 0.000194692, 0.0002232326, 0.0002553766, 0.0002915015, 0.0003320157, 0.00037736, 0.0004280092, 0.0004844737, 0.0005473018, 0.0006170807, 0.0006944384, 0.0007800461, 0.0008746188, 0.0009789181, 0.0010937533, 0.0012199834, 0.0013585185, 0.0015103219, 0.0016764114, 0.0018578609, 0.0020558023, 0.0022714268, 0.0025059864, 0.0027607955, 0.0030372323, 0.0033367399, 0.0036608279, 0.0040110736, 0.0043891231, 0.0047966926, 0.0052355691, 0.0057076118, 0.0062147529, 0.0067589985, 0.0073424292, 0.0079672011, 0.0086355463, 0.0093497735, 0.0101122685, 0.0109254949, 0.011791994, 0.0127143855, 0.0136953675, 0.0147377169, 0.0158442892, 0.0170180188, 0.0182619187, 0.0195790806, 0.0209726748, 0.0224459497, 0.0240022316, 0.0256449245, 0.0273775098, 0.0292035452, 0.0311266652, 0.0331505797, 0.0352790737, 0.0375160069, 0.0398653126, 0.0423309973, 0.0449171398, 0.0476278905, 0.0504674705, 0.0534401708, 0.0565503513, 0.0598024404, 0.0632009333, 0.0667503919, 0.0704554433, 0.074320779, 0.078351154, 0.0825513859, 0.0869263538, 0.0914809975, 0.0962203163, 0.1011493683, 0.1062732693, 0.111597192, 0.1171263649, 0.1228660716, 0.1288216497, 0.1349984903, 0.1414020366, 0.1480377837, 0.1549112773, 0.1620281134, 0.1693939372, 0.1770144428, 0.1848953722, 0.1930425148, 0.2014617071, 0.2101588319, 0.219139818, 0.2284106396, 0.2379773164, 0.2478459127, 0.2580225378, 0.2685133453, 0.2793245335, 0.2904623449, 0.3019330665, 0.3137430302, 0.3258986125, 0.3384062351, 0.3512723653, 0.3645035164, 0.3781062484, 0.3920871688, 0.4064529329, 0.4212102456, 0.4363658616, 0.4519265872, 0.4678992814, 0.4842908574, 0.5011082844, 0.5183585892, 0.5360488584, 0.5541862404, 0.5727779482, 0.5918312616, 0.6113535304, 0.6313521773, 0.6518347014, 0.6728086819, 0.6942817823, 0.7162617546, 0.7387564438, 0.7617737937, 0.785321852, 0.8094087764, 0.8340428416, 0.8592324462, 0.8849861205, 0.9113125351, 0.9382205104, 0.9657190263, 0.9938172336, 1.0225244661, 1.0518502533, 1.0818043358, 1.1123966802, 1.1436374973, 1.1755372609, 1.2081067294, 1.2413569686, 1.2752993784, 1.3099457208, 1.3453081525, 1.3813992602, 1.4182321006, 1.4558202448, 1.4941778289, 1.5333196103, 1.573261032, 1.614018295, 1.6556084416, 1.6980494505, 1.7413603456, 1.7855613233, 1.8306738995, 1.8767210817, 1.9237275734, 1.9717200158, 2.0207272795, 2.0707808186, 2.1219151068, 2.1741681851, 2.2275823605, 2.2822051216, 2.338090372, 2.3953001459, 2.4539070606, 2.5139977584, 2.575676075]

    AARD_calc = sum([abs(Sanjari(T, 345.0, 26.40E5, 0.3170)-P*1E6)/(P*1E6) for T, P in zip(Ts_dat, Ps_dat)])/len(Ts_dat)
    assert_allclose(AARD_calc, 0.006445800342803334)

    # Supposed to be 1.387 %, small difference
    # Functions are identical- data simply must be different.
    # Or different method sof calculating AARD. No worries.
    AARD_calc = sum([abs(Lee_Kesler(T, 345.0, 26.40E5, 0.3170)-P*1E6)/(P*1E6) for T, P in zip(Ts_dat, Ps_dat)])/len(Ts_dat)
    assert_allclose(AARD_calc, 0.01370923047231833)

    # Supposed to be 0.785 %, small difference; plus formula matches
    AARD_calc = sum([abs(Ambrose_Walton(T, 345.0, 26.40E5, 0.3170)-P*1E6)/(P*1E6) for T, P in zip(Ts_dat, Ps_dat)])/len(Ts_dat)
    assert_allclose(AARD_calc, 0.00841629399152493)

def test_Psat_IAPWS():
    Psat = Psat_IAPWS(300.)
    assert_close(Psat, 3536.5894130130105, rtol=1e-12)
    Psat = Psat_IAPWS(500.)
    assert_close(Psat, 2638897.7562732217, rtol=1e-12)
    Psat = Psat_IAPWS(600.)
    assert_close(Psat, 12344314.578376647, rtol=1e-12)
    
    # Points obtained with CoolProp
    assert_close(dPsat_IAPWS_dT(300.0), 207.88388134164325327, rtol=5e-15)
    assert_close(dPsat_IAPWS_dT(500.0), 49008.859762957346618, rtol=5e-15)
    assert_close(dPsat_IAPWS_dT(600.0), 160315.50500711359916, rtol=5e-15)
    assert_close(dPsat_IAPWS_dT(100.0), -0.094715356431800557785, rtol=5e-15)

def test_Tsat_IAPWS():
    # 3 checks from book
    Tsat = Tsat_IAPWS(1E5)
    assert_close(Tsat, 372.7559186113376, rtol=1e-13)
    Tsat = Tsat_IAPWS(1E6)
    assert_close(Tsat, 453.0356323914666, rtol=1e-13)
    Tsat = Tsat_IAPWS(1E7)
    assert_close(Tsat, 584.1494879985282, rtol=1e-13)
    
    assert_close(Tsat_IAPWS(Psat_IAPWS(553.123521)), 553.123521, rtol=1e-13)
    
    T_min = 159.77353993926621
    assert_close(Tsat_IAPWS(Psat_IAPWS(T_min)), T_min, rtol=1e-13)

def test_Psub_Clapeyron():
    Psub = Psub_Clapeyron(250.0, Tt=273.15, Pt=611.0, Hsub_t=51100.0)
    assert_allclose(Psub, 76.06457150831804)
