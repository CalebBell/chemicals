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

import numpy as np
import pytest
from fluids.numerics import assert_close, assert_close1d, linspace

from chemicals.temperature import ITS90_68_difference, T_converter, Ts_27, Ts_48, Ts_68, Ts_76, diffs_27, diffs_48, diffs_68, diffs_76


def test_data():
    Ts_sums_calc = [np.array(i).sum() for i in [Ts_68, Ts_48, Ts_76, Ts_27]]
    Ts_sums = [186818.69999999998, 175181.39999999997, 368, 133893.09999999998]
    assert_close1d(Ts_sums_calc, Ts_sums)
    diffs_sums_calc = [abs(np.array(i)).sum() for i in [diffs_68, diffs_48, diffs_76, diffs_27]]
    diffs_sums = [46.304000000000016, 151.31800000000001, 0.038800000000000001, 411.17999999999995]
    assert_close1d(diffs_sums_calc, diffs_sums)


def test_conversion():
    # TODO actually test data points instead of covering everything in a slow test

    T2 = T_converter(500, 'ITS-68', 'ITS-48')
    assert_close(T2, 499.9470092992346)

    high_scales = ('ITS-90', 'ITS-68', 'ITS-27', 'ITS-48')

    for scale1 in high_scales:
        for scale2 in high_scales:
            T = T_converter(1000.0, scale1, scale2)
            assert_close(T_converter(T, scale2, scale1), 1000.0)

    mid_scales = ('ITS-90', 'ITS-68', 'ITS-48')

    for Ti in linspace(100, 1000, 10):
        for scale1 in mid_scales:
            for scale2 in mid_scales:
                T = T_converter(Ti, scale1, scale2)
                assert_close(T_converter(T, scale2, scale1), Ti, rtol=1e-6)

    low_scales = ('ITS-90', 'ITS-68', 'ITS-76')

    for Ti in (15, 17, 19, 21, 23, 25):
        for scale1 in low_scales:
            for scale2 in low_scales:
                T = T_converter(Ti, scale1, scale2)
                assert_close(T_converter(T, scale2, scale1), Ti)

    with pytest.raises(ValueError):
        T_converter(10, 'ITS-27', 'ITS-48')

    with pytest.raises(ValueError):
        T_converter(10, 'FAIL', 'ITS-48')

    with pytest.raises(ValueError):
        T_converter(10, 'ITS-76', 'FAIL')


def test_diff_68():
    dTs_calc = [ITS90_68_difference(i) for i in [13.7, 70, 80.5, 298.15, 1000, 1500]]


    dTs = [0, 0.006818871618271216, 0, -0.006253950277664615,
           0.01231818956580355, -0.31455]
    assert_close1d(dTs, dTs_calc)


def test_scale_boundaries():
    # ITS-76 (5K to 27K)
    T1 = T_converter(5.0, 'ITS-76', 'ITS-90')
    assert_close(T1, 4.9999)
    T2 = T_converter(27.0, 'ITS-76', 'ITS-90')
    assert_close(T2, 26.9959)

    # ITS-68 (14K to 4300K)
    T3 = T_converter(14.0, 'ITS-68', 'ITS-90')
    assert_close(T3, 13.994)
    T4 = T_converter(4300.0, 'ITS-68', 'ITS-90')
    assert_close(T4, 4297.42)

def test_fixed_points():
    # Triple point of water (273.16 K)
    T2 = T_converter(273.16, 'ITS-68', 'ITS-90')
    assert_close(T2, 273.16)
    T3 = T_converter(273.16, 'ITS-48', 'ITS-90')
    assert_close(T3, 273.16)

    # Freezing point of tin (505.078 K)
    T4 = T_converter(505.078, 'ITS-68', 'ITS-90')
    assert_close(T4, 505.037995)

    # Freezing point of zinc (692.677 K)
    T5 = T_converter(692.677, 'ITS-68', 'ITS-90')
    assert_close(T5, 692.62431)

def test_scale_overlaps():
    # Test point where ITS-48 and ITS-68 overlap (904K)
    T1_48 = T_converter(904.0, 'ITS-48', 'ITS-90')
    T1_68 = T_converter(904.0, 'ITS-68', 'ITS-90')
    assert_close(T1_48, T1_68, rtol=1e-3)

    # Test point where ITS-68 and ITS-76 overlap (20K)
    T2_76 = T_converter(20.0, 'ITS-76', 'ITS-90')
    T2_68 = T_converter(20.0, 'ITS-68', 'ITS-90')
    assert_close(T2_76, T2_68, rtol=1e-3)

def test_critical_points():
    # Critical point of water (647.096 K)
    T1 = 647.096
    T1_90_68 = T_converter(T1, 'ITS-90', 'ITS-68')
    T1_68_90 = T_converter(T1_90_68, 'ITS-68', 'ITS-90')
    assert_close(T1, T1_68_90)

    # Critical point of oxygen (154.581 K)
    T2 = 154.581
    T2_90_68 = T_converter(T2, 'ITS-90', 'ITS-68')
    T2_68_90 = T_converter(T2_90_68, 'ITS-68', 'ITS-90')
    assert_close(T2, T2_68_90)

def test_edge_cases():
    with pytest.raises(ValueError):
        T_converter(4.9, 'ITS-76', 'ITS-90')  # Below ITS-76 minimum

    with pytest.raises(ValueError):
        T_converter(27.1, 'ITS-76', 'ITS-90')  # Above ITS-76 maximum

    with pytest.raises(ValueError):
        T_converter(13.9, 'ITS-68', 'ITS-90')  # Below ITS-68 minimum

    with pytest.raises(ValueError):
        T_converter(4300.1, 'ITS-68', 'ITS-90')  # Above ITS-68 maximum

def test_linear_conversion_points():
    # ITS-76 <-> ITS-90 points
    assert_close(T_converter(5, 'ITS-76', 'ITS-90'), 4.9999)
    assert_close(T_converter(7.2, 'ITS-76', 'ITS-90'), 7.19968002423)
    assert_close(T_converter(9.4, 'ITS-76', 'ITS-90'), 9.39946065917)
    assert_close(T_converter(11.6, 'ITS-76', 'ITS-90'), 11.5992515766)
    assert_close(T_converter(13.8, 'ITS-76', 'ITS-90'), 13.7989249491)
    assert_close(T_converter(16, 'ITS-76', 'ITS-90'), 15.9986)
    assert_close(T_converter(18.2, 'ITS-76', 'ITS-90'), 18.1981594572)
    assert_close(T_converter(20.4, 'ITS-76', 'ITS-90'), 20.3976794018)
    assert_close(T_converter(22.6, 'ITS-76', 'ITS-90'), 22.5971182748)
    assert_close(T_converter(24.8, 'ITS-76', 'ITS-90'), 24.7965647013)
    assert_close(T_converter(27, 'ITS-76', 'ITS-90'), 26.9959)

    # ITS-68 <-> ITS-90 points
    assert_close(T_converter(14, 'ITS-68', 'ITS-90'), 13.994)
    assert_close(T_converter(442.6, 'ITS-68', 'ITS-90'), 442.561701802)
    assert_close(T_converter(871.2, 'ITS-68', 'ITS-90'), 871.085683799)
    assert_close(T_converter(1299.8, 'ITS-68', 'ITS-90'), 1299.5701936)
    assert_close(T_converter(1728.4, 'ITS-68', 'ITS-90'), 1727.98677691)
    assert_close(T_converter(2157, 'ITS-68', 'ITS-90'), 2156.3455346)
    assert_close(T_converter(2585.6, 'ITS-68', 'ITS-90'), 2584.66153245)
    assert_close(T_converter(3014.2, 'ITS-68', 'ITS-90'), 3012.92823042)
    assert_close(T_converter(3442.8, 'ITS-68', 'ITS-90'), 3441.14084888)
    assert_close(T_converter(3871.4, 'ITS-68', 'ITS-90'), 3869.30209144)
    assert_close(T_converter(4300, 'ITS-68', 'ITS-90'), 4297.42)

    # ITS-48 <-> ITS-90 points
    assert_close(T_converter(93.15, 'ITS-48', 'ITS-90'), 93.17)
    assert_close(T_converter(511.15, 'ITS-48', 'ITS-90'), 511.167218737)
    assert_close(T_converter(929.15, 'ITS-48', 'ITS-90'), 929.418086942)
    assert_close(T_converter(1347.15, 'ITS-48', 'ITS-90'), 1348.3450273)
    assert_close(T_converter(1765.15, 'ITS-48', 'ITS-90'), 1766.93759523)
    assert_close(T_converter(2183.15, 'ITS-48', 'ITS-90'), 2185.46620148)
    assert_close(T_converter(2601.15, 'ITS-48', 'ITS-90'), 2604.28326291)
    assert_close(T_converter(3019.15, 'ITS-48', 'ITS-90'), 3023.02890902)
    assert_close(T_converter(3437.15, 'ITS-48', 'ITS-90'), 3441.8674108)
    assert_close(T_converter(3855.15, 'ITS-48', 'ITS-90'), 3860.90135997)
    assert_close(T_converter(4273.15, 'ITS-48', 'ITS-90'), 4279.95)

    # ITS-27 <-> ITS-90 points
    assert_close(T_converter(903.15, 'ITS-27', 'ITS-90'), 903.23)
    assert_close(T_converter(1240.15, 'ITS-27', 'ITS-90'), 1241.42131767)
    assert_close(T_converter(1577.15, 'ITS-27', 'ITS-90'), 1577.47759384)
    assert_close(T_converter(1914.15, 'ITS-27', 'ITS-90'), 1912.85010146)
    assert_close(T_converter(2251.15, 'ITS-27', 'ITS-90'), 2247.4203438)
    assert_close(T_converter(2588.15, 'ITS-27', 'ITS-90'), 2581.45972044)
    assert_close(T_converter(2925.15, 'ITS-27', 'ITS-90'), 2914.74270106)
    assert_close(T_converter(3262.15, 'ITS-27', 'ITS-90'), 3247.21287016)
    assert_close(T_converter(3599.15, 'ITS-27', 'ITS-90'), 3578.81025784)
    assert_close(T_converter(3936.15, 'ITS-27', 'ITS-90'), 3909.07648984)
    assert_close(T_converter(4273.15, 'ITS-27', 'ITS-90'), 4238.05)