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
from random import Random
from chemicals.temperature import ITS90_68_difference, T_converter, Ts_27, Ts_48, Ts_68, Ts_76, diffs_27, diffs_48, diffs_68, diffs_76


def test_data():
    Ts_sums_calc = [np.array(i).sum() for i in [Ts_68, Ts_48, Ts_76, Ts_27]]
    Ts_sums = [186818.69999999998, 175181.39999999997, 368, 133893.09999999998]
    assert_close1d(Ts_sums_calc, Ts_sums)
    diffs_sums_calc = [abs(np.array(i)).sum() for i in [diffs_68, diffs_48, diffs_76, diffs_27]]
    diffs_sums = [46.304000000000016, 151.31800000000001, 0.038800000000000001, 411.17999999999995]
    assert_close1d(diffs_sums_calc, diffs_sums)


def test_conversion():
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


    # Random T90 -> ITS-76 points and round trips
    assert_close(T_converter(13.435713414167, 'ITS-90', 'ITS-76'), 13.436758011396)
    assert_close(T_converter(T_converter(13.435713414167, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 13.435713414167)
    assert_close(T_converter(14.320405585456, 'ITS-90', 'ITS-76'), 14.321565696878)
    assert_close(T_converter(T_converter(14.320405585456, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 14.320405585456)
    assert_close(T_converter(14.626918647779, 'ITS-90', 'ITS-76'), 14.628146618509)
    assert_close(T_converter(T_converter(14.626918647779, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 14.626918647779)
    assert_close(T_converter(16.987430025932, 'ITS-90', 'ITS-76'), 16.989027641309)
    assert_close(T_converter(T_converter(16.987430025932, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 16.987430025932)
    assert_close(T_converter(17.073897086401, 'ITS-90', 'ITS-76'), 17.075513283896)
    assert_close(T_converter(T_converter(17.073897086401, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 17.073897086401)
    assert_close(T_converter(18.260794273576, 'ITS-90', 'ITS-76'), 18.262647665938)
    assert_close(T_converter(T_converter(18.260794273576, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 18.260794273576)
    assert_close(T_converter(19.209670487466, 'ITS-90', 'ITS-76'), 19.211708200102)
    assert_close(T_converter(T_converter(19.209670487466, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 19.209670487466)
    assert_close(T_converter(20.734166060193, 'ITS-90', 'ITS-76'), 20.736593177144)
    assert_close(T_converter(T_converter(20.734166060193, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 20.734166060193)
    assert_close(T_converter(24.619006017206, 'ITS-90', 'ITS-76'), 24.622383782420)
    assert_close(T_converter(T_converter(24.619006017206, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 24.619006017206)
    assert_close(T_converter(26.200580731023, 'ITS-90', 'ITS-76'), 26.204440144401)
    assert_close(T_converter(T_converter(26.200580731023, 'ITS-90', 'ITS-76'), 'ITS-76','ITS-90'), 26.200580731023)

    # Random T90 -> ITS-68 points and round trips
    assert_close(T_converter(100.656051429236, 'ITS-90', 'ITS-68'), 100.647040261164)
    assert_close(T_converter(T_converter(100.656051429236, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 100.656051429236)
    assert_close(T_converter(318.460545436143, 'ITS-90', 'ITS-68'), 318.472077466932)
    assert_close(T_converter(T_converter(318.460545436143, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 318.460545436143)
    assert_close(T_converter(387.436178520803, 'ITS-90', 'ITS-68'), 387.464674690384)
    assert_close(T_converter(T_converter(387.436178520803, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 387.436178520803)
    assert_close(T_converter(2280.843626060949, 'ITS-90', 'ITS-68'), 2281.571645701955)
    assert_close(T_converter(T_converter(2280.843626060949, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 2280.843626060949)
    assert_close(T_converter(2448.638988848594, 'ITS-90', 'ITS-68'), 2449.477910998919)
    assert_close(T_converter(T_converter(2448.638988848594, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 2448.638988848594)
    assert_close(T_converter(3349.179834571059, 'ITS-90', 'ITS-68'), 3350.750661627371)
    assert_close(T_converter(T_converter(3349.179834571059, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 3349.179834571059)
    assert_close(T_converter(3407.333513222301, 'ITS-90', 'ITS-68'), 3408.961984805847)
    assert_close(T_converter(T_converter(3407.333513222301, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 3407.333513222301)
    assert_close(T_converter(3582.608658018462, 'ITS-90', 'ITS-68'), 3584.401939925675)
    assert_close(T_converter(T_converter(3582.608658018462, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 3582.608658018462)
    assert_close(T_converter(3742.872067385867, 'ITS-90', 'ITS-68'), 3744.837138119218)
    assert_close(T_converter(T_converter(3742.872067385867, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 3742.872067385867)
    assert_close(T_converter(3981.107191722345, 'ITS-90', 'ITS-68'), 3983.328923919461)
    assert_close(T_converter(T_converter(3981.107191722345, 'ITS-90', 'ITS-68'), 'ITS-68','ITS-90'), 3981.107191722345)

    # Random T90 -> ITS-48 points and round trips
    assert_close(T_converter(587.537100132141, 'ITS-90', 'ITS-48'), 587.501584795817)
    assert_close(T_converter(T_converter(587.537100132141, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 587.537100132141)
    assert_close(T_converter(692.366741369814, 'ITS-90', 'ITS-48'), 692.344588159499)
    assert_close(T_converter(T_converter(692.366741369814, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 692.366741369814)
    assert_close(T_converter(1826.436909160389, 'ITS-90', 'ITS-48'), 1824.585952661900)
    assert_close(T_converter(T_converter(1826.436909160389, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 1826.436909160389)
    assert_close(T_converter(2022.133734217255, 'ITS-90', 'ITS-48'), 2019.973788475658)
    assert_close(T_converter(T_converter(2022.133734217255, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 2022.133734217255)
    assert_close(T_converter(2274.475984915300, 'ITS-90', 'ITS-48'), 2271.978485162272)
    assert_close(T_converter(T_converter(2274.475984915300, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 2274.475984915300)
    assert_close(T_converter(2768.019869149050, 'ITS-90', 'ITS-48'), 2764.644373694928)
    assert_close(T_converter(T_converter(2768.019869149050, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 2768.019869149050)
    assert_close(T_converter(3355.761956877384, 'ITS-90', 'ITS-48'), 3351.203345457763)
    assert_close(T_converter(T_converter(3355.761956877384, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 3355.761956877384)
    assert_close(T_converter(3433.632798425905, 'ITS-90', 'ITS-48'), 3428.931655843199)
    assert_close(T_converter(T_converter(3433.632798425905, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 3433.632798425905)
    assert_close(T_converter(4041.866073267261, 'ITS-90', 'ITS-48'), 4035.688458555037)
    assert_close(T_converter(T_converter(4041.866073267261, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 4041.866073267261)
    assert_close(T_converter(4183.774670532954, 'ITS-90', 'ITS-48'), 4177.164252155940)
    assert_close(T_converter(T_converter(4183.774670532954, 'ITS-90', 'ITS-48'), 'ITS-48','ITS-90'), 4183.774670532954)

    # Random T90 -> ITS-27 points and round trips
    assert_close(T_converter(966.471627470517, 'ITS-90', 'ITS-27'), 965.716314921189)
    assert_close(T_converter(T_converter(966.471627470517, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 966.471627470517)
    assert_close(T_converter(1794.702412792593, 'ITS-90', 'ITS-27'), 1795.361363571567)
    assert_close(T_converter(T_converter(1794.702412792593, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 1794.702412792593)
    assert_close(T_converter(2440.376619569769, 'ITS-90', 'ITS-27'), 2445.677432276432)
    assert_close(T_converter(T_converter(2440.376619569769, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 2440.376619569769)
    assert_close(T_converter(2818.772407687345, 'ITS-90', 'ITS-27'), 2827.963374114178)
    assert_close(T_converter(T_converter(2818.772407687345, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 2818.772407687345)
    assert_close(T_converter(2965.912585574560, 'ITS-90', 'ITS-27'), 2976.958286576851)
    assert_close(T_converter(T_converter(2965.912585574560, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 2965.912585574560)
    assert_close(T_converter(2982.217569467931, 'ITS-90', 'ITS-27'), 2993.459134966685)
    assert_close(T_converter(T_converter(2982.217569467931, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 2982.217569467931)
    assert_close(T_converter(2984.581625145705, 'ITS-90', 'ITS-27'), 2995.851378838990)
    assert_close(T_converter(T_converter(2984.581625145705, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 2984.581625145705)
    assert_close(T_converter(3200.884407978739, 'ITS-90', 'ITS-27'), 3215.114723062496)
    assert_close(T_converter(T_converter(3200.884407978739, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 3200.884407978739)
    assert_close(T_converter(3512.317533393310, 'ITS-90', 'ITS-27'), 3531.550641197200)
    assert_close(T_converter(T_converter(3512.317533393310, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 3512.317533393310)
    assert_close(T_converter(4083.581024594283, 'ITS-90', 'ITS-27'), 4114.606650574897)
    assert_close(T_converter(T_converter(4083.581024594283, 'ITS-90', 'ITS-27'), 'ITS-27','ITS-90'), 4083.581024594283)


def test_roundtrip_conversions():
    """Test roundtrip conversions between all temperature scales using random points
    in valid overlapping ranges, accounting for differences between scales."""
    
    # Calculate the maximum differences for each scale to ITS-90
    multiplier =(1+0.2) # safety factor
    max_diff_27 = max(abs(min(diffs_27)), abs(max(diffs_27)))*multiplier
    max_diff_48 = max(abs(min(diffs_48)), abs(max(diffs_48)))*multiplier
    max_diff_68 = max(abs(min(diffs_68)), abs(max(diffs_68)))*multiplier
    max_diff_76 = max(abs(min(diffs_76)), abs(max(diffs_76)))*multiplier
    
    # Adjust ranges to account for these differences
    scales = ['ITS-68', 'ITS-48', 'ITS-76', 'ITS-27']
    ranges = {
        'ITS-68': (14.0 + max_diff_68, 4300.0 - max_diff_68),
        'ITS-48': (93.15 + max_diff_48, 4273.15 - max_diff_48),
        'ITS-76': (5.0 + max_diff_76, 27.0 - max_diff_76),
        'ITS-27': (903.15 + max_diff_27, 4273.15 - max_diff_27)
    }
    
    rng = Random(0)
    pts = 50

    for scale1 in scales:
        for scale2 in scales:
            if scale1 >= scale2:  # Skip redundant combinations
                continue
                
            min_temp = max(ranges[scale1][0], ranges[scale2][0])
            max_temp = min(ranges[scale1][1], ranges[scale2][1])
            
            if min_temp >= max_temp:
                continue
                
            test_points = [rng.uniform(min_temp, max_temp) for _ in range(pts)]
            for val in test_points:
                # Test both directions
                result1 = T_converter(T_converter(val, scale1, scale2), scale2, scale1)
                result2 = T_converter(T_converter(val, scale2, scale1), scale1, scale2)
                assert_close(result1, val, rtol=1e-13)
                assert_close(result2, val, rtol=1e-13)


def test_T_converter_range_helpful_strings():
    with pytest.raises(ValueError) as exc_info:
        T_converter(13.998, 'ITS-68', 'ITS-90')
    assert "Temperature 13.998 K is below minimum 13.999 K for ITS-68 scale" in str(exc_info.value)
    
    with pytest.raises(ValueError) as exc_info:
        T_converter(4300.0002, 'ITS-68', 'ITS-90')
    assert "Temperature 4300.0002 K is above maximum 4300.0001 K for ITS-68 scale" in str(exc_info.value)
    
    # Test ITS-76 ranges
    with pytest.raises(ValueError) as exc_info:
        T_converter(4.9998, 'ITS-76', 'ITS-90')
    assert "Temperature 4.9998 K is below minimum 4.9999 K for ITS-76 scale" in str(exc_info.value)
    
    with pytest.raises(ValueError) as exc_info:
        T_converter(27.0002, 'ITS-76', 'ITS-90')
    assert "Temperature 27.0002 K is above maximum 27.0001 K for ITS-76 scale" in str(exc_info.value)