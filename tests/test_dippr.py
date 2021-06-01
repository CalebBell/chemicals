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

from fluids.numerics import assert_close, assert_close1d
from scipy.misc import derivative
from scipy.integrate import quad
import pytest

from chemicals.dippr import *

def test_Eqs():
    a = EQ100(300, 276370., -2090.1, 8.125, -0.014116, 0.0000093701)
    assert_close(a, 75355.81)

    a = EQ101(300, 73.649, -7258.2, -7.3037, 4.1653E-6, 2)
    assert_close(a, 3537.44834545549)

    a = EQ102(300, 1.7096E-8, 1.1146, 0, 0)
    assert_close(a, 9.860384711890639e-06)

    a = EQ104(300, 0.02222, -26.38, -16750000, -3.894E19, 3.133E21)
    assert_close(a, -1.1204179007265151)

    a = EQ106(300, 647.096, 0.17766, 2.567, -3.3377, 1.9699)
    assert_close(a, 0.07231499373541)

    a = EQ107(300., 33363., 26790., 2610.5, 8896., 1169.)
    assert_close(a, 33585.90452768923)

    a = EQ114(20, 33.19, 66.653, 6765.9, -123.63, 478.27)
    assert_close(a, 19423.948911676463)

    a = EQ116(300., 647.096, 17.863, 58.606, -95.396, 213.89, -141.26)
    assert_close(a, 55.17615446406527)

    a = EQ127(20., 3.3258E4, 3.6199E4, 1.2057E3, 1.5373E7, 3.2122E3, -1.5318E7, 3.2122E3)
    assert_close(a, 33258.0)

    # Random coefficients
    a = EQ115(300, 0.01, 0.002, 0.0003, 0.00004)
    assert_close(a, 37.02960772416336)


def test_EQ105_more():
    a = EQ105(300., 0.70824, 0.26411, 507.6, 0.27537)
    assert_close(a, 7.593170096339236)

    d1 = EQ105(300, 0.70824, 0.26411, 507.6, 0.27537, order=1)
    assert_close(d1, -0.01048317827611096, rtol=1e-13)
    assert_close(derivative(lambda T: EQ105(T, 0.70824, 0.26411, 507.6, 0.27537), 300, dx=1e-3), d1)

    d2 = EQ105(300, 0.70824, 0.26411, 507.6, 0.27537, order=2)
    assert_close(d2, -2.2118503164569428e-05, rtol=1e-13)
    assert_close(derivative(lambda T: EQ105(T, 0.70824, 0.26411, 507.6, 0.27537, order=1), 300, dx=1e-3), d2)

    d3 = EQ105(300, 0.70824, 0.26411, 507.6, 0.27537, order=3)
    assert_close(derivative(lambda T: EQ105(T, 0.70824, 0.26411, 507.6, 0.27537, order=2), 300, dx=1e-3), d3)
    
    
    # avoid complex numbers
    kwargs = {'T': 195.0, 'A': 7247.0, 'B': 0.418, 'C': 5.2, 'D': 0.24, 'order': 0}
    assert_close(EQ105(**kwargs), 17337.32057416268)


def test_EQ106_more():
    d1_analytical = EQ106(300, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=1)
    d1_numerical = derivative(lambda T: EQ106(T, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=0), 300.0, dx=1e-4)
    assert_close(d1_analytical, -0.00019544758472630726, rtol=1e-13)
    assert_close(d1_analytical, d1_numerical, rtol=1e-8)

    d2_analytical = EQ106(300, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=2)
    d2_numerical = derivative(lambda T: EQ106(T, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=1), 300.0, dx=1e-4)
    assert_close(d2_analytical, 2.113551968259993e-07, rtol=1e-13)
    assert_close(d2_analytical, d2_numerical, rtol=1e-8)

    d3_analytical = EQ106(300, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=3)
    d3_numerical = derivative(lambda T: EQ106(T, 647.096, 0.17766, 2.567, -3.3377, 1.9699, order=2), 300.0, dx=1e-4)
    assert_close(d3_analytical, -5.524739890945817e-09, rtol=1e-13)
    assert_close(d3_analytical, d3_numerical, rtol=1e-8)
    
    # check that this regression set of points does not produce an error
    overflow_kwargs = {'T': 304.747, 'Tc': 405.4, 'A': 56.743647419038744, 'B': -75.36242555958763,
                   'C': -141.1028969227863, 'D': -254.76349199392695, 'E': -442.5916844036474, 'order': 0}
    EQ106(**overflow_kwargs)
    
    # Point exactly on the critical point that was an error, needed an if statement.
    assert 0.0 == EQ106(T=473.2, Tc=473.2, **{'A': 4761730.0, 'B': -11.5565, 'C': 30.6629, 'D': -31.89366, 'E': 12.67797})


def test_EQ101_more():
    assert_close(EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=1), 208.00259945348506, rtol=1e-13)
    assert_close(derivative(lambda T: EQ101(T, 73.649, -7258.2, -7.3037, 4.1653E-6, 2), 300, dx=1e-4),
                 EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=1), rtol=1e-9)

    assert_close(EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=2), 10.64524169930617, rtol=1e-13)
    assert_close(derivative(lambda T: EQ101(T, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=1), 300, dx=1e-4),
                 EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=2), rtol=1e-9)

    assert_close(EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=3), 0.4566096458105825, rtol=1e-13)
    assert_close(derivative(lambda T: EQ101(T, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=2), 300, dx=1e-4),
                 EQ101(300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=3), rtol=1e-9)

    with pytest.raises(ValueError):
        EQ101(300, 73.649, -7258.2, -7.3037, 4.1653E-6, 2, order=1000)

def test_EQ127_more():
    # T derivative
    coeffs = (3.3258E4, 3.6199E4, 1.2057E3, 1.5373E7, 3.2122E3, -1.5318E7, 3.2122E3)
    diff_1T = derivative(EQ127, 50,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ127(50., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T, 0.000313581049006, rtol=1E-4)

    # Integral
    int_50 = EQ127(50., *coeffs, order=-1)
    int_20 = EQ127(20., *coeffs, order=-1)
    numerical_1T = quad(EQ127, 20, 50, args=coeffs)[0]
    assert_close(int_50 - int_20, numerical_1T)
    assert_close(numerical_1T, 997740.00147014)

    # Integral over T
    T_int_50 = EQ127(50., *coeffs, order=-1j)
    T_int_20 = EQ127(20., *coeffs, order=-1j)

    to_int = lambda T :EQ127(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 20, 50)[0]
    assert_close(T_int_50 - T_int_20, numerical_1_over_T)
    assert_close(T_int_50 - T_int_20, 30473.9971912935)

    with pytest.raises(Exception):
        EQ127(20., *coeffs, order=1E100)

def test_EQ115_more():
    args = (0.01, 0.002, 0.0003, 0.00004, 10.0)
    d1_analytical = EQ115(300.0, *args, order=1)
    d1_numerical = derivative(lambda T: EQ115(T, *args, order=0), 300.0, dx=1e-2)
    assert_close(d1_analytical, 0.8888181148504066, rtol=1e-13)
    assert_close(d1_analytical, d1_numerical)

    d2_analytical = EQ115(300.0, *args, order=2)
    d2_numerical = derivative(lambda T: EQ115(T, *args, order=1), 300.0, dx=1e-2)
    assert_close(d2_analytical, 0.024294699592116293, rtol=1e-13)
    assert_close(d2_analytical, d2_numerical)

    d3_analytical = EQ115(300.0, *args, order=3)
    d3_numerical = derivative(lambda T: EQ115(T, *args, order=2), 300.0, dx=1e-2)
    assert_close(d3_analytical, 0.0007252940633608988, rtol=1e-13)
    assert_close(d3_analytical, d3_numerical)
    
    # Check case avoid overflow
    kwargs = {'T': 400.05, 'A': 1.0, 'B': 1.0, 'C': 1.0, 'D': 1.0, 'E': 1.0, 'order': 0}
    EQ115(**kwargs)



def test_EQ116_more():
    # T derivative
    coeffs = (647.096, 17.863, 58.606, -95.396, 213.89, -141.26)
    diff_1T = derivative(EQ116, 50,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ116(50., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T_analytical, 0.020379262711650914)

    # Integral
    int_50 = EQ116(50., *coeffs, order=-1)
    int_20 = EQ116(20., *coeffs, order=-1)
    numerical_1T = quad(EQ116, 20, 50, args=coeffs)[0]
    assert_close(int_50 - int_20, numerical_1T)
    assert_close(int_50 - int_20, 1636.962423782701)

    # Integral over T
    T_int_50 = EQ116(50., *coeffs, order=-1j)
    T_int_20 = EQ116(20., *coeffs, order=-1j)

    to_int = lambda T :EQ116(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 20, 50)[0]
    assert_close(T_int_50 - T_int_20, numerical_1_over_T)
    assert_close(T_int_50 - T_int_20, 49.95109104018752)

    with pytest.raises(Exception):
        EQ116(20., *coeffs, order=1E100)

def test_EQ107_more():
    # T derivative
    coeffs = (33363., 26790., 2610.5, 8896., 1169.)
    diff_1T = derivative(EQ107, 250,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ107(250., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T_analytical, 1.985822265543943)

    # Integral
    int_50 = EQ107(50., *coeffs, order=-1)
    int_20 = EQ107(20., *coeffs, order=-1)
    numerical_1T = quad(EQ107, 20, 50, args=coeffs)[0]
    assert_close(int_50 - int_20, numerical_1T)
    assert_close(numerical_1T, 1000890.0)

    # Integral over T
    T_int_50 = EQ107(50., *coeffs, order=-1j)
    T_int_20 = EQ107(20., *coeffs, order=-1j)

    to_int = lambda T :EQ107(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 20, 50)[0]
    assert_close(T_int_50 - T_int_20, numerical_1_over_T)
    assert_close(T_int_50 - T_int_20, 30570.20768751744)


    with pytest.raises(Exception):
        EQ107(20., *coeffs, order=1E100)
    
    # Case that requires overflow handling
    EQ107(**{'T': 377.77777777777777, 'A': 1539249.2020718465, 'B': -46807441.804555826, 'C': -409401.9169728528, 'D': -2164118.45731599, 'E': 339.5030595758336, 'order': 0})


def test_EQ114_more():
    # T derivative
    coeffs = (33.19, 66.653, 6765.9, -123.63, 478.27)
    diff_1T = derivative(EQ114, 20,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ114(20., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T, 1135.38618941)

    # Integral
    int_50 = EQ114(30., *coeffs, order=-1)
    int_20 = EQ114(20., *coeffs, order=-1)
    numerical_1T = quad(EQ114, 20, 30, args=coeffs)[0]
    assert_close(int_50 - int_20, numerical_1T)
    assert_close(int_50 - int_20, 295697.48978888744)

#     Integral over T
    T_int_50 = EQ114(30., *coeffs, order=-1j)
    T_int_20 = EQ114(20., *coeffs, order=-1j)

    to_int = lambda T :EQ114(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 20, 30)[0]
    assert_close(T_int_50 - T_int_20, numerical_1_over_T)
    assert_close(T_int_50 - T_int_20, 11612.331762721366)

    with pytest.raises(Exception):
        EQ114(20., *coeffs, order=1E100)


def test_EQ102_more():
    # T derivative
    coeffs = (1.7096E-8, 1.1146, 1.1, 2.1)
    diff_1T = derivative(EQ102, 250,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ102(250., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T, 3.5861274167602139e-08)

    # Integral
    int_250 = EQ102(250., *coeffs, order=-1)
    int_220 = EQ102(220., *coeffs, order=-1)
    numerical_1T = quad(EQ102, 220, 250, args=coeffs)[0]
    assert_close(int_250 - int_220, numerical_1T)
    assert_close(int_250 - int_220, 0.00022428562125110119)

#     Integral over T
    T_int_250 = EQ102(250., *coeffs, order=-1j)
    T_int_220 = EQ102(220., *coeffs, order=-1j)

    to_int = lambda T :EQ102(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 220, 250)[0]
    assert_close(T_int_250 - T_int_220, numerical_1_over_T)
    assert_close(T_int_250 - T_int_220, 9.5425212178091671e-07)
#
    with pytest.raises(Exception):
        EQ102(20., *coeffs, order=1E100)
        
    # No overflow
    EQ102(T = 194.6, A = 0.0, B = 4000.0, C = 0.75, D = 0.0, order = 0)


def test_EQ100_more():
    # T derivative
    coeffs = (276370., -2090.1, 8.125, -0.014116, 0.0000093701)
    diff_1T = derivative(EQ100, 250,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ100(250., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T, -88.7187500531)

    # Integral
    int_250 = EQ100(250., *coeffs, order=-1)
    int_220 = EQ100(220., *coeffs, order=-1)
    numerical_1T = quad(EQ100, 220, 250, args=coeffs)[0]
    assert_close(int_250 - int_220, numerical_1T)
    assert_close(int_250 - int_220, 2381304.7021859996)

#     Integral over T
    T_int_250 = EQ100(250., *coeffs, order=-1j)
    T_int_220 = EQ100(220., *coeffs, order=-1j)

    to_int = lambda T :EQ100(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 220, 250)[0]
    assert_close(T_int_250 - T_int_220, numerical_1_over_T)
    assert_close(T_int_250 - T_int_220, 10152.09780143667)

    with pytest.raises(Exception):
        EQ100(20., *coeffs, order=1E100)


def test_EQ104_more():
    # T derivative
    coeffs = (0.02222, -26.38, -16750000, -3.894E19, 3.133E21)
    diff_1T = derivative(EQ104, 250,  dx=1E-3, order=21, args=coeffs)
    diff_1T_analytical = EQ104(250., *coeffs, order=1)
    assert_close(diff_1T, diff_1T_analytical, rtol=1E-3)
    assert_close(diff_1T, 0.0653824814073)

    # Integral
    int_250 = EQ104(250., *coeffs, order=-1)
    int_220 = EQ104(220., *coeffs, order=-1)
    numerical_1T = quad(EQ104, 220, 250, args=coeffs)[0]
    assert_close(int_250 - int_220, numerical_1T)
    assert_close(int_250 - int_220, -127.91851427119406)

#     Integral over T
    T_int_250 = EQ104(250., *coeffs, order=-1j)
    T_int_220 = EQ104(220., *coeffs, order=-1j)

    to_int = lambda T :EQ104(T, *coeffs)/T
    numerical_1_over_T = quad(to_int, 220, 250)[0]
    assert_close(T_int_250 - T_int_220, numerical_1_over_T)
    assert_close(T_int_250 - T_int_220, -0.5494851210308727)

    with pytest.raises(Exception):
        EQ104(20., *coeffs, order=1E100)

def test_EQ101_fitting():
    T, A, B, C, D, E = 300.0, 73.649, -7258.2, -7.3037, 4.1653E-6, 2
    der_num = [derivative(lambda A: EQ101(T, A, B, C, D, E), A, dx=A*1e-5),
                 derivative(lambda B: EQ101(T, A, B, C, D, E), B, dx=B*1e-5),
                 derivative(lambda C: EQ101(T, A, B, C, D, E), C, dx=C*1e-5),
                 derivative(lambda D: EQ101(T, A, B, C, D, E), D, dx=D*1e-5),
                 derivative(lambda E: EQ101(T, A, B, C, D, E), E, dx=E*1e-5),
              ]
    
    der_expect = [[3537.44834545549, 11.791494484851635, 20176.835877810598, 318370351.0909941, 7563.831703366002]]
    der_analytical = EQ101_fitting_jacobian([T], A, B, C, D, E)
    assert_close1d(der_analytical, [der_num])
    assert_close1d(der_analytical, der_expect, rtol=1e-13)


def test_EQ102_fitting():
    T, A, B, C, D = 300.0, 2e-6, 0.42, 900.0, -4e4
    der_num = [derivative(lambda A: EQ102(T, A, B, C, D), A, dx=A*1e-5),
                 derivative(lambda B: EQ102(T, A, B, C, D), B, dx=B*1e-5),
                 derivative(lambda C: EQ102(T, A, B, C, D), C, dx=C*1e-5),
                 derivative(lambda D: EQ102(T, A, B, C, D), D, dx=D*1e-5)]
    der_expect = [[3.08662207669995, 3.521084181393621e-05, -5.787416393812407e-09, -1.9291387979374693e-11]]
    der_analytical = EQ102_fitting_jacobian([T], A, B, C, D)
    assert_close1d(der_analytical, [der_num])
    assert_close1d(der_analytical, der_expect, rtol=1e-13)

def test_EQ105_fitting():
    T, A, B, C, D = 300., 0.70824, 0.26411, 507.6, 0.27537
    der_num = [derivative(lambda A: EQ105(T, A, B, C, D), A, dx=A*1e-5),
                 derivative(lambda B: EQ105(T, A, B, C, D), B, dx=B*1e-5),
                 derivative(lambda C: EQ105(T, A, B, C, D), C, dx=C*1e-5),
                 derivative(lambda D: EQ105(T, A, B, C, D), D, dx=D*1e-5)]
    der_expect = [[10.721182221195127, -51.225752844842916, 0.006195731841673147, -7.0661094492142285]]
    der_analytical = EQ105_fitting_jacobian([T], A, B, C, D)
    assert_close1d(der_analytical, [der_num])
    assert_close1d(der_analytical, der_expect, rtol=1e-13)


    T, A, B, C, D = 304, 3733.087734888731, 0.30552803535622014, 301.6993863907116, 0.3415888512743092
    der_num = [derivative(lambda A: EQ105(T, A, B, C, D), A, dx=A*1e-5),
                 derivative(lambda B: EQ105(T, A, B, C, D), B, dx=B*1e-5),
                 derivative(lambda C: EQ105(T, A, B, C, D), C, dx=C*1e-5),
                 derivative(lambda D: EQ105(T, A, B, C, D), D, dx=D*1e-5)]
    der_analytical = EQ105_fitting_jacobian([T], A, B, C, D)
    assert_close1d(der_analytical, [der_num])


def test_EQ106_fitting():
    T, Tc, A, B, C, D, E = 300, 647.096, 0.17766, 2.567, -3.3377, 1.9699, 0.25
    der_num = [derivative(lambda A: EQ106(T, Tc, A, B, C, D, E), A, dx=A*1e-5),
                 derivative(lambda B: EQ106(T, Tc, A, B, C, D, E), B, dx=B*1e-5),
                 derivative(lambda C: EQ106(T, Tc, A, B, C, D, E), C, dx=C*1e-5),
                 derivative(lambda D: EQ106(T, Tc, A, B, C, D, E), D, dx=D*1e-5),
                 derivative(lambda E: EQ106(T, Tc, A, B, C, D, E), E, dx=E*1e-5),
              ]
        
    der_expect = [[0.4007741423076755, -0.04435095583995359, -0.020561534535812425, -0.009532527415937865, -0.004419372434354963]]
    der_analytical = EQ106_fitting_jacobian([T], Tc, A, B, C, D, E)
    assert_close1d(der_analytical, [der_num])
    assert_close1d(der_analytical, der_expect, rtol=1e-13)
    
    # Test case with errors
    EQ106_fitting_jacobian(Ts=[466.],Tc = 466.0, A = 47700.0, B = 0.37, C = 0.0, D = 0.0, E = 0.0)

def test_EQ107_fitting():
    T, A, B, C, D, E = 250.0, 33363., 26790., 2610.5, 8896., 1169
    der_num = [derivative(lambda A: EQ107(T, A, B, C, D, E), A, dx=A*1e-5),
                 derivative(lambda B: EQ107(T, A, B, C, D, E), B, dx=B*3e-3, order=3),
                 derivative(lambda C: EQ107(T, A, B, C, D, E), C, dx=C*1e-3, order=9),
                 derivative(lambda D: EQ107(T, A, B, C, D, E), D, dx=D*1e-5, order=3),
                 derivative(lambda E: EQ107(T, A, B, C, D, E), E, dx=E*1e-5),
              ]
    
    der_expect = [[1.0, 3.7138247806865474e-07, -7.197214038362036e-05, 0.00758947296962729, -0.42452325330497365]]
    der_analytical = EQ107_fitting_jacobian([T], A, B, C, D, E)
    assert_close1d(der_analytical, [der_num])
    assert_close1d(der_analytical, der_expect, rtol=1e-13)

