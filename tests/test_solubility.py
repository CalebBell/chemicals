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
from fluids.numerics import assert_close, assert_close1d

from chemicals.solubility import (
    Henry_constants,
    Henry_converter,
    Henry_pressure,
    Henry_pressure_mixture,
    Tm_depression_eutectic,
    d2Henry_constants_dT2,
    dHenry_constants_dT,
    solubility_eutectic,
    solubility_parameter,
    hansen_delta_d,
    hansen_delta_d_methods,
    hansen_delta_p,
    hansen_delta_p_methods,
    hansen_delta_h,
    hansen_delta_h_methods,
    ALSHERI_HANSEN,
    HSPIPY,
    WDR_SCHRIER,
    MANUEL_RUBEN_2022,
    alsheri_hansen_data, 
    hspipy_data, 
    wdr_schrier_data,
    manuel_ruben_2022_data,
)


def test_hansen_delta():
    # Test values for ethanol (64-17-5)
    d_calc = hansen_delta_d('64-17-5')
    p_calc = hansen_delta_p('64-17-5')
    h_calc = hansen_delta_h('64-17-5')
    
    assert_close1d([d_calc, p_calc, h_calc], [15800.0, 8800.0, 19400.0])
    
    hits_d = hansen_delta_d_methods('64-17-5')
    hits_p = hansen_delta_p_methods('64-17-5')
    hits_h = hansen_delta_h_methods('64-17-5')
    
    assert hits_d == hits_p
    assert hits_p == hits_h
    
    with pytest.raises(Exception):
        hansen_delta_d('64-17-5', method='BADMETHOD')
    with pytest.raises(Exception):
        hansen_delta_p('64-17-5', method='BADMETHOD')
    with pytest.raises(Exception):
        hansen_delta_h('64-17-5', method='BADMETHOD')

    assert None is hansen_delta_d('99234-43-50-0')
    assert [] == hansen_delta_d_methods('99234-43-50-0')
    assert None is hansen_delta_p('99234-43-50-0')
    assert [] == hansen_delta_p_methods('99234-43-50-0')
    assert None is hansen_delta_h('99234-43-50-0')
    assert [] == hansen_delta_h_methods('99234-43-50-0')

    w_methods = hansen_delta_d_methods('7732-18-5')
    assert w_methods == ['MANUEL_RUBEN_2022', 'HSPIPY']
    w_methods = hansen_delta_p_methods('7732-18-5')
    assert w_methods == ['MANUEL_RUBEN_2022', 'HSPIPY']
    w_methods = hansen_delta_h_methods('7732-18-5')
    assert w_methods == ['MANUEL_RUBEN_2022', 'HSPIPY']

    w_values = [hansen_delta_d('7732-18-5', method=method) for method in w_methods]
    assert_close1d(w_values, [15500.0, 15500.0])
    w_values = [hansen_delta_p('7732-18-5', method=method) for method in w_methods]
    assert_close1d(w_values, [16000.0, 16000.0])
    w_values = [hansen_delta_h('7732-18-5', method=method) for method in w_methods]
    assert_close1d(w_values, [42300.0, 42300.0])


    test_cases = ['64-17-5', '7732-18-5', '71-43-2', '67-56-1']
    for cas in test_cases:
        d_methods = set(hansen_delta_d_methods(cas))
        p_methods = set(hansen_delta_p_methods(cas))
        h_methods = set(hansen_delta_h_methods(cas))
        assert d_methods == p_methods == h_methods, f"Inconsistent methods for {cas}"

    assert isinstance(hansen_delta_d('64-17-5'), float)
    assert isinstance(hansen_delta_p('64-17-5'), float)
    assert isinstance(hansen_delta_h('64-17-5'), float)

    assert_close(hansen_delta_d('138495-42-8'), 11600.0)  # Non-zero δD
    assert_close(hansen_delta_p('138495-42-8'), 0.0)      # Zero δP
    assert_close(hansen_delta_h('138495-42-8'), 0.0)      # Zero δH
    
def test_hansen_solubility_all_data():
    s1 = alsheri_hansen_data.index
    s2 = hspipy_data.index
    s3 = wdr_schrier_data.index
    s4 = manuel_ruben_2022_data.index
    tots = {}
    tots_exp = {
        ALSHERI_HANSEN: {'delta_d': 15302100.0, 'delta_p': 6738300.0, 'delta_h': 6748800.0},
        HSPIPY: {'delta_d': 4126200.0, 'delta_p': 1692500.0, 'delta_h': 1982300.0},
        WDR_SCHRIER: {'delta_d': 1884800.0, 'delta_p': 841900.0, 'delta_h': 982400.0},
        MANUEL_RUBEN_2022: {'delta_d': 20575260.0, 'delta_p': 8828240.0, 'delta_h': 8904890.0}
    }
    # These should match the sums of the respective series
    for s, method in zip([s1, s2, s3, s4], [ALSHERI_HANSEN, HSPIPY, WDR_SCHRIER, MANUEL_RUBEN_2022]):
        tots[method] = {
            'delta_d': sum([hansen_delta_d(i, method=method) for i in s]),
            'delta_p': sum([hansen_delta_p(i, method=method) for i in s]),
            'delta_h': sum([hansen_delta_h(i, method=method) for i in s])
        }

    for method in [ALSHERI_HANSEN, HSPIPY, WDR_SCHRIER, MANUEL_RUBEN_2022]:
        for param in ['delta_d', 'delta_p', 'delta_h']:
            assert_close(tots[method][param], tots_exp[method][param], rtol=1e-11)

    s = set()
    s.update(s1)
    s.update(s2)
    s.update(s3)
    s.update(s4)
    assert len(s) == 1294  # Correct total number of unique compounds


def test_solubility():
    # From [1]_, matching examples 1 and 2.
    x1 = solubility_eutectic(293.15, 369.4, 18640., 0, 0, 1)
    x2 = solubility_eutectic(T=260., Tm=278.68, Hm=9952., Cpl=0, Cps=0, gamma=3.0176)
    x3 = solubility_eutectic(T=260., Tm=278.68, Hm=9952., Cpl=195, Cps=60, gamma=3.0176)
    assert_close1d([x1, x2, x3], [0.20626915125512824, 0.2434007130748926, 0.2533343734537043])

    dTm1 = Tm_depression_eutectic(353.35, 19110, 0.02)
    dTm2 = Tm_depression_eutectic(353.35, 19110, M=0.4, MW=40.)
    assert_close1d([dTm1, dTm2], [1.0864598583150953, 0.8691678866520763])
    with pytest.raises(Exception):
        Tm_depression_eutectic(353.35, 19110)


def test_solubility_parameter():
    delta = solubility_parameter(T=298.2, Hvapm=26403.3, Vml=0.000116055)
    assert delta == solubility_parameter(298.2, 26403.3, 0.000116055)
    assert_close(delta, 14357.681538173534)

    assert None is solubility_parameter(T=3500.2, Hvapm=26403.3, Vml=0.000116055)


def test_Henry_converter():
    from chemicals.solubility import (
        HENRY_SCALES_BUNSEN,
        HENRY_SCALES_HBP,
        HENRY_SCALES_HBP_SI,
        HENRY_SCALES_HCC,
        HENRY_SCALES_HCP,
        HENRY_SCALES_HCP_MOLALITY,
        HENRY_SCALES_HXP,
        HENRY_SCALES_KHCC,
        HENRY_SCALES_KHPC,
        HENRY_SCALES_KHPC_SI,
        HENRY_SCALES_KHPX,
        HENRY_SCALES_SI,
    )
    test_values = [1.2E-05, 0.0012159, 0.0297475, 1.20361E-08, 0.00121956,
                   2.19707E-05, 0.0272532, 45515.2, 83333.3, 0.822436, 33.6163,
                   4611823929.1419935]
    test_scales = [HENRY_SCALES_HCP, HENRY_SCALES_HCP_MOLALITY, HENRY_SCALES_HCC,
                   HENRY_SCALES_HBP_SI, HENRY_SCALES_HBP, HENRY_SCALES_HXP,
                   HENRY_SCALES_BUNSEN, HENRY_SCALES_KHPX, HENRY_SCALES_KHPC_SI,
                   HENRY_SCALES_KHPC, HENRY_SCALES_KHCC, HENRY_SCALES_SI]
    for v, scales in zip(test_values, test_scales):
        for scale in scales:
            calc = Henry_converter(v, old_scale=scale, new_scale='Hxp', rhom=55341.9, MW=18.01528)
            # Best we can match given the digits provided
            assert_close(calc, 2.19707E-05, rtol=2e-6)
            recalc = Henry_converter(v, old_scale=scale, new_scale=scale, rhom=55341.9, MW=18.01528)
            assert_close(v, recalc, rtol=1e-14)


def test_Henry_pressure():
    H = Henry_pressure(300.0, A=15.0, B=300.0, C=.04, D=1e-3, E=1e2, F=1e-5)
    assert_close(H, 37105004.47898146)

def test_Henry_pressure_mixture():
    H = Henry_pressure_mixture([1072330.36341, 744479.751106, None], zs=[.48, .48, .04])
    assert_close(H, 893492.1611602883)

def test_Henry_constants():
    lnHenry_matrix = [[0.0, 0.0, 0.0], [22.13581843104147, 0.0, 0.0], [22.239038459475733, 0.0, 0.0]]
    Hs = Henry_constants(lnHenry_matrix, [0.8, 0.15, 0.05], [False, True, True], True)
    assert_close1d([0.0, 4106424071.093, 4552937470.331], Hs)

    # Test no error, henry goes too high
    kwargs = {'lnHenry_matrix': [[0.0, 0.0, 0.0, 10.219309168602226], [25.321910952013628, 0.0, 20.22048862340195, 18.156499431206274], [24.269833800974848, 15.365770016174015, 0.0, 18.57423773432834], [15.37228504674701, 0.0, 0.0, 0.0]], 'zs': [0.000242872042323124, 0.5796842606361365, 0.41952722776289886, 0.000545639558641481], 'henry_components': [False, True, True, False], 'skip_zero': True, 'Hs': [0.0, 0.0, 0.0, 0.0]}
    Henry_constants(**kwargs)


def test_dHenry_constants_dT():
    lnHenry_matrix = [[0.0, 0.0, 0.0], [22.13581843104147, 0.0, 0.0], [22.239038459475733, 0.0, 0.0]]
    dlnHenry_matrix_dT = [[0.0, 0.0, 0.0], [0.017113988888888904, 0.0, 0.0], [0.015461911111111101, 0.0, 0.0]]
    calc = dHenry_constants_dT(lnHenry_matrix, dlnHenry_matrix_dT, [0.8, 0.15, 0.05], [False, True, True], True)
    assert_close1d(calc, [0.0, 70277295.92576516, 70397114.46071726])

def test_d2Henry_constants_dT2():
    lnHenry_matrix = [[0.0, 0.0, 0.0], [22.13581843104147, 0.0, 0.0], [22.239038459475733, 0.0, 0.0]]
    dlnHenry_matrix_dT = [[0.0, 0.0, 0.0], [0.017113988888888904, 0.0, 0.0], [0.015461911111111101, 0.0, 0.0]]
    d2lnHenry_matrix_dT2 = [[0.0, 0.0, 0.0], [-0.0004070325925925928, 0.0, 0.0], [-0.00034016518518518524, 0.0, 0.0]]

    calc = d2Henry_constants_dT2(lnHenry_matrix, dlnHenry_matrix_dT, d2lnHenry_matrix_dT2, [0.8, 0.15, 0.05], [False, True, True], True)
    assert_close1d(calc, [0.0, -468723.574327235, -460276.89146166])
