"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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

import numpy as np
import pytest
from fluids.numerics import assert_close, assert_close1d, assert_close2d

from chemicals.utils import (
    SG,
    API_to_SG,
    Cp_minus_Cv,
    Joule_Thomson,
    Parachor,
    SG_to_API,
    Vfs_to_zs,
    Vm_to_rho,
    Watson_K,
    Z,
    d2ns_to_dn2_partials,
    d2xs_to_d2xsn1,
    dns_to_dn_partials,
    dxs_to_dn_partials,
    dxs_to_dns,
    dxs_to_dxsn1,
    hash_any_primitive,
    isentropic_exponent,
    isentropic_exponent_PT,
    isentropic_exponent_PV,
    isentropic_exponent_TV,
    isobaric_expansion,
    isothermal_compressibility,
    mix_component_flows,
    mix_multiple_component_flows,
    mixing_logarithmic,
    mixing_power,
    mixing_simple,
    molar_velocity_to_velocity,
    none_and_length_check,
    normalize,
    phase_identification_parameter,
    phase_identification_parameter_phase,
    property_mass_to_molar,
    property_molar_to_mass,
    radius_of_gyration,
    recursive_copy,
    remove_zeros,
    rho_to_Vm,
    solve_flow_composition_mix,
    speed_of_sound,
    to_num,
    v_molar_to_v,
    v_to_v_molar,
    vapor_mass_quality,
    velocity_to_molar_velocity,
    ws_to_zs,
    zs_to_Vfs,
    zs_to_ws,
    ms_to_ns,
    ns_to_ms,
    ns_to_Qls,
    Qls_to_ns,
    Qls_to_ms,
    ms_to_Qls,
    rho_to_Baume_light, Baume_light_to_rho, SG_to_Baume_light, Baume_light_to_SG, rho_to_Baume_heavy, Baume_heavy_to_rho, SG_to_Baume_heavy, Baume_heavy_to_SG,
)


def test_recursive_copy():
    import array
    from copy import deepcopy
    from decimal import Decimal
    from fractions import Fraction
    from collections import namedtuple

    test_cases = [None,
                  True, False,
                  float,
                  int,
                  str,
                 -1, 0, 1, 2, 2000, 2**65,
                  -1.1234, 0, 2e200, float("nan"),
                  1j, 1.1234j, 1234132+123.234j, -12341-1234j,
                  'a', 'b', 'asdfsdfasdf', 'asdf asdfadf',
                  Decimal(2),Fraction(2.34),
                 ]

    numpy_test_cases = [ np.complex128(1), np.complex64(1), np.float16(1),
                        np.float32(1), np.float64(1), np.int16(1), np.int32(1), np.int64(1),
                        np.int8(1), np.longlong(1), np.uint16(1), np.uint32(1), np.uint64(1),
                        np.uint8(1), np.ulonglong(1),]
    try:
        numpy_test_cases.append(np.float128(1))
    except:
        # not supported on windows https://github.com/cupy/cupy/issues/2916
        pass
    try:
        numpy_test_cases.append(np.complex256(1))
    except:
        # not supported on windows https://github.com/cupy/cupy/issues/2916
        pass


    InnerTuple = namedtuple('InnerTuple', ['a', 'b'])
    OuterTuple = namedtuple('OuterTuple', ['x', 'y'])
    named_tuple_1 = InnerTuple(1, (2, 3, 4))
    named_tuple_2 = OuterTuple(InnerTuple(1, (2, 3, 4)), 'abc')
    named_tuple_3 = OuterTuple(InnerTuple('asdfas', (2, '125dsaf', 4)), frozenset([0]))

    tuple_cases = [(12,21,34,3, None, -1, 'asd', 1j),
                   (12., (123.352, None, 4, (123.352, None, 4), [123]), ),
                   tuple(),
                   tuple(range(129)),
                   named_tuple_1,
                   named_tuple_2,
                   named_tuple_3
                  ]

    list_cases = [[2,21,34,3, None, -1, 'asd', 1j],
                   [12., (123.352, None, 4, (123.352, None, 4), [123]), ],
                   list(),
                   list(range(129)),

    ]

    set_cases = [{1,2,'a', 1j, np.float32(1)},
    ]

    frozenset_cases = [frozenset([1,2,'a', 1j, np.float32(1)]),
    ]

    numpy_cases = [np.arange(13),
                   np.arange(13, dtype=int),
                   np.arange(13, dtype=np.uint32),
                   np.arange(13, dtype=np.float32),
                   np.arange(13, dtype=np.float16),
                  ]

    other_cases = [
    ]


    byte_cases = [b'asdfasd', b'1232',
                  bytes.fromhex('2Ef0 F1f2  '),
                 ]

    test_cases += numpy_test_cases
    test_cases += tuple_cases
    test_cases += list_cases
    test_cases += numpy_cases
    test_cases += other_cases
    test_cases += set_cases
    test_cases += frozenset_cases
    test_cases += byte_cases




    for case in test_cases:
        implemented_copy = recursive_copy(case)
        original_copy = deepcopy(case)
        assert type(implemented_copy) == type(original_copy)
        assert hash_any_primitive(implemented_copy) == hash_any_primitive(original_copy)


    unhashable_cases = [range(10),
                        range(1,2,10),
                        range(1,29,10),
                        bytearray([0,1,3,2]),
                        array.array('l', [-11111111, 22222222, -33333333, 44444444]),
                       ]


    for case in unhashable_cases:
        implemented_copy = recursive_copy(case)
        original_copy = deepcopy(case)
        assert type(implemented_copy) == type(original_copy)
        assert id(implemented_copy) != case


    read_only_cases = [{'a': 3, 'b': 5}.items(),]

    for case in read_only_cases:
        implemented_copy = recursive_copy(case)
        assert type(implemented_copy) == type(case)
        assert implemented_copy is case


    # numpy case with object arrays
    original = np.array([{2**66}, {2**67}])
    copy = recursive_copy(original)
    original[1].add(1)
    assert len(copy[1]) == 1



def test_to_num():
    assert to_num(['1', '1.1', '1E5', '0xB4', '']) == [1.0, 1.1, 100000.0, '0xB4', None]


def test_remove_zeros():
    a = remove_zeros([0, 1e-9, 1], 1e-12)
    assert type(a) == list
    assert_close1d(a, [9.99999998999e-13, 9.99999998999e-10, 0.999999998999])
    b = remove_zeros(np.array([0, 1e-9, 1]), 1e-12)
    assert_close1d(b, [9.99999998999e-13, 9.99999998999e-10, 0.999999998999])
    assert type(b) == np.ndarray

    assert remove_zeros([.3, .2, .5], 1e-12) == [0.3, 0.2, 0.5]

    c = remove_zeros(np.array([.3, .2, .5]), 1e-12)
    assert_close1d(c, [.3, .2, .5])
    assert type(c) == np.ndarray


def test_none_and_length_check():
    assert True is none_and_length_check([[1, 2, 3]])

    assert True is none_and_length_check(([1, 1], [1, 1], [1, 30], [10, 0]), length=2)
    assert True is none_and_length_check(([1, 1], [1, 1], [1, 30], [10, 0]))

    assert False is none_and_length_check(([1, 1], [None, 1], [1, 30], [10, 0]), length=2)
    assert False is none_and_length_check(([None, 1], [1, 1], [1, 30], [10, 0]))
    assert False is none_and_length_check(([1, 1], [None, 1], [1, None], [10, 0]), length=2)
    assert False is none_and_length_check(([None, 1], [1, 1], [1, 30], [None, 0]))

    assert False is none_and_length_check(([1, 1, 1], [1, 1], [1, 30], [10, 0]), length=2)
    assert False is none_and_length_check(([1, 1], [1, 1, 1], [1, 30], [10, 0]))
    assert False is none_and_length_check(([1, 1, 1], [1, 1], [1, 30, 1], [10, 0]), length=2)
    assert False is none_and_length_check(([1, 1], [1, 1, 1], [1, 30], [10, 0, 1]))
    assert False is none_and_length_check(([1, 1, 1], [1, 1, 1], [1, 30, 1], [10, 0]), length=3)
    assert False is none_and_length_check(([1, 1], [1, 1, 1], [1, 30, 1], [10, 0, 1]))

    assert True is none_and_length_check(([1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=8)
    assert True is none_and_length_check(([1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]))

    assert False is none_and_length_check(([1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=9)
    assert False is none_and_length_check(([1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=7)
    assert False is none_and_length_check(([1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, None]), length=8)
    assert False is none_and_length_check(([1, 1, None, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=8)

    # Test list input instead of tuples
    assert True is none_and_length_check([[1, 1], [1, 1], [1, 30], [10, 0]], length=2)
    assert True is none_and_length_check([[1, 1], [1, 1], [1, 30], [10, 0]])

    assert True is none_and_length_check([[1, 1], [1, 1], [1, 30], [10, 0]], length=2)
    assert True is none_and_length_check([[1, 1], [1, 1], [1, 30], [10, 0]])

    # Test with numpy arrays
    assert True is none_and_length_check((np.array([1, 1, 1, 1, 1, 1, 1, 1]), np.array([1, 1, 1, 1, 1, 1, 1, 1]), [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=8)
    assert False is none_and_length_check((np.array([1, 1, 1, 1, 1, 1, 1]), np.array([1, 1, 1, 1, 1, 1, 1, 1]), [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=8)
    assert False is none_and_length_check((np.array([1, 1, 1, 1, 1, 1, 1, 7]), np.array([1, 1, 1, 1, 1, 1, 1, 1]), [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]), length=7)
    assert True is none_and_length_check(np.array([[1, 1], [1, 1], [1, 30], [10, 0]]))

    assert True is none_and_length_check(np.array([[1, 1], [1, 1], [1, 30], [10, 0]]))
    assert True is none_and_length_check(np.array([[1, 1], [1, 1], [1, 30], [10, 0]]), length=2)
    assert False is none_and_length_check(np.array([[1, 1], [1, 1], [1, 30], [10, 0]], dtype=object), length=3)
    assert False is none_and_length_check(np.array([[1, 1], [1, 1, 10], [1, 30], [10, 0]], dtype=object), length=3)


def test_property_molar_to_mass():
    assert_close(property_molar_to_mass(500, 18.015), 27754.648903691366, rtol=1e-14)

def test_property_mass_to_molar():
    assert_close(property_mass_to_molar(20.0, 18.015), 0.3603, rtol=1e-14)

def test_zs_to_ws():
    ws_calc = zs_to_ws([0.5, 0.5], [10, 20])
    ws = [0.3333333333333333, 0.6666666666666666]
    assert_close1d(ws_calc, ws)


def test_ws_to_zs():
    zs_calc = ws_to_zs([0.3333333333333333, 0.6666666666666666], [10, 20])
    zs = [0.5, 0.5]
    assert_close1d(zs_calc, zs)


def test_zs_to_Vfs():
    Vfs_calc = zs_to_Vfs([0.637, 0.363], [8.0234e-05, 9.543e-05])
    Vfs = [0.5960229712956298, 0.4039770287043703]
    assert_close1d(Vfs_calc, Vfs)


def test_Vfs_to_zs():
    zs_calc = Vfs_to_zs([0.596, 0.404], [8.0234e-05, 9.543e-05])
    zs = [0.6369779395901142, 0.3630220604098858]
    assert_close1d(zs_calc, zs)

def test_molar_velocity_conv():
    v_molar_calc = v_to_v_molar(500, 18.015)
    assert_close(v_molar_calc, 67.10998435404377, rtol=1e-12)

    v_calc = v_molar_to_v(67.10998435404377, 18.015)
    assert_close(v_calc, 500, rtol=1e-12)


def test_Z():
    Z_calc = Z(600, P=1E6, V=0.00463)
    assert_close(Z_calc, 0.9281016730797027)

def test_vapor_mass_quality():
    x = vapor_mass_quality(0.5, 60, 30)
    assert_close(x, 1/3.)

def test_Vm_to_rho():
    rho = Vm_to_rho(0.000132, 86.18)
    assert_close(rho, 652.8787878787879)


def test_rho_to_Vm():
    Vm = rho_to_Vm(652.9, 86.18)
    assert_close(Vm, 0.00013199571144126206)
    assert_close(652.9, Vm_to_rho(rho_to_Vm(652.9, 86.18), 86.18))


def test_isentropic_exponent():
    k = isentropic_exponent(33.6, 25.27)
    assert_close(k, 1.329639889196676)

    k = isentropic_exponent_TV(Cv=23.98081290153672, Vm=4.730885141495376e-05, dP_dT_V=509689.2959155567)
    assert_close(k, 2.005504495083913, rtol=1e-14)

    k = isentropic_exponent_PT(Cp=38.36583283578205, P=100000000.0, dV_dT_P=9.407705210161724e-08)
    assert_close(k, 1.32487270350443, rtol=1e-14)

    k = isentropic_exponent_PV(Cp=38.36583283578205, Cv=23.98081290153672, Vm=4.730885141495376e-05, P=100000000.0, dP_dV_T=-5417785576072.434)
    assert_close(k, 4.100576762582646, rtol=1e-14)

def test_Parachor():
    P = Parachor(100.15888, 800.8088185536124, 4.97865317223119, 0.02672166960656005)
    assert_close(P, 5.088443542210164e-05)



def test_mixing_simple():
    prop = mixing_simple([0.1, 0.9], [0.01, 0.02])
    assert_close(prop, 0.019)

def test_mixing_logarithmic():
    prop = mixing_logarithmic([0.1, 0.9], [0.01, 0.02])
    assert_close(prop, 0.01866065983073615)

def test_normalize():
    fractions_calc = normalize([3, 2, 1])
    fractions = [0.5, 0.3333333333333333, 0.16666666666666666]
    assert_close1d(fractions, fractions_calc)

def test_isobaric_expansion():
    beta = isobaric_expansion(0.000130229900873546, 1.58875261849113e-7)
    assert_close(beta, 0.0012199599384121608)

def test_isothermal_compressibility():
    kappa = isothermal_compressibility(0.000130229900873546, -2.72902118209903e-13)
    assert_close(2.09554116511916e-9, kappa)


def test_phase_identification_parameter():
    PIP = phase_identification_parameter(0.000130229900874, 582169.397484, -3.66431747236e+12, 4.48067893805e+17, -20518995218.2)
    assert_close(PIP, 11.33428990564796)

    assert 'l' == phase_identification_parameter_phase(-20518995218.2, 0.000130229900874, 582169.397484, -3.66431747236e+12, 4.48067893805e+17)
    # Artificially give a value to make it be solid
    assert 's' == phase_identification_parameter_phase(+20518995218.2)

def test_Cp_minus_Cv():
    d = Cp_minus_Cv(299, 582232.475794113, -3665180614672.253)
    assert_close(d, 27.654681381642394)


def test_speed_of_sound():
    # Matches exactly
    w = speed_of_sound(V=0.00229754, dP_dV=-3.5459e+08, Cp=153.235, Cv=132.435, MW=67.152)
    assert_close(w, 179.5868138460819)

    # No MW version
    w = speed_of_sound(V=0.00229754, dP_dV=-3.5459e+08, Cp=153.235, Cv=132.435)
    assert_close(w, 46.537593457316525)


def test_mu_JT():
    # Matches exactly
    mu_JT = Joule_Thomson(T=390, V=0.00229754, Cp=153.235, dV_dT=1.226396e-05)
    assert_close(mu_JT, 1.621956080529905e-05)

    mu_JT = Joule_Thomson(T=390, V=0.00229754, Cp=153.235, beta=0.005337865717245402)
    assert_close(mu_JT, 1.621956080529905e-05)

    with pytest.raises(Exception):
        Joule_Thomson(T=390, V=0.00229754, Cp=153.235)

def test_API_SG():
    SG = API_to_SG(60.62)
    assert_close(SG, 0.7365188423901728)

    API = SG_to_API(0.7365)
    assert_close(API, 60.62491513917175)

def test_SG():
    sg = SG(860)
    assert_close(sg, 0.8608461408159591)

def test_Watson_K():
    assert_close(Watson_K(400, .8), 11.20351186639291, rtol=1e-13)


def test_mix_component_flows():
    names, flows = mix_component_flows(['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], 1, 1, [0.5, 0.5], [0.5, 0.5])
    assert names == ['64-17-5', '67-56-1', '7732-18-5']
    assert_close1d(flows, [ 0.5,  0.5,  1.])

    names, flows = mix_component_flows(['7732-18-5', '67-56-1'], ['7732-18-5', '67-56-1'], 1, 1, [0.2, 0.8], [0.3, 0.7])
    assert names == ['7732-18-5', '67-56-1']
    assert_close1d(flows,  [0.5, 1.5])

def test_mix_multiple_component_flows():
    names, flows = mix_multiple_component_flows([['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1']], [1, 1], [[0.5, 0.5], [0.5, 0.5]])
    assert names == ['64-17-5', '67-56-1', '7732-18-5']
    assert_close1d(flows, [ 0.5,  0.5,  1.])

    args = ([['7732-18-5', '64-17-5'], ['7732-18-5', '64-17-5'], ['7732-18-5', '67-56-1'], ['7732-18-5', '67-56-1']], [1, 1, 1, 1], [[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
    names, flows = mix_multiple_component_flows(*args)
    assert names == ['64-17-5', '67-56-1', '7732-18-5']
    assert_close1d(flows, [1.0, 1.0, 2.0])

    names, flows = mix_multiple_component_flows([['7732-18-5', '64-17-5']], [1], [[0.5 , 0.5]])
    assert names == ['7732-18-5', '64-17-5']
    assert_close1d(flows, [ 0.5, 0.5,])


def test_solve_flow_composition_mix():
    # Needs fast assert_close1d
    Fs = [3600, None, None, None, None]
    zs = [None, .1, .2, None, None]
    ws = [None, None, None, .01, .02]
    MWs = [18.01528, 46.06844, 32.04186, 72.151, 142.286]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    Fs_expect = [3600, 519.3039148597746, 1038.6078297195493, 17.44015034881175, 17.687253669610733]
    zs_expect = [0.6932356751002142, 0.1, 0.2, 0.003358370666918819, 0.0034059542328670383]
    ws_expect = [0.5154077420893427, 0.19012206531421305, 0.26447019259644433, 0.01, 0.02]

    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)


    Fs = [3600, 519.3039148597746, None, None, None]
    zs = [None, None, .2, None, None]
    ws = [None, None, None, .01, .02]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [3600, 519.3039148597746, 1038.6078297195493, None, None]
    zs = [None, None, None, None, None]
    ws = [None, None, None, .01, .02]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [3600, 519.3039148597746, 1038.6078297195493, None, None]
    zs = [None, None, None, 0.003358370666918819, None]
    ws = [None, None, None, None, .02]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [None, None, None, None, 17.687253669610733]
    zs = [0.6932356751002142, .1, .2, None, None]
    ws = [None, None, None, .01, None]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [None, None, None, None, 17.687253669610733]
    zs = [0.6932356751002142, .1, .2, 0.003358370666918816, None]
    ws = [None, None, None, None, None]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [3600., 519.3039148597746, None, None, 17.687253669610733]
    zs = [None, None, .2, 0.003358370666918816, None]
    ws = [None, None, None, None, None]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

    Fs = [3600.0, 519.3039148597746, 1038.6078297195495, 17.440150348811738, 17.687253669610733]
    zs = [None, None, None, None, None]
    ws = [None, None, None, None, None]
    Fs, zs, ws = solve_flow_composition_mix(Fs, zs, ws, MWs)
    assert_close1d(Fs, Fs_expect)
    assert_close1d(ws, ws_expect)
    assert_close1d(zs, zs_expect)

@pytest.mark.fuzz
@pytest.mark.slow
def test_solve_flow_composition_mix_fuzz():
    # Solve a large-scale problem
    from random import randint, uniform
    N = 1000
    Fs = [3600]
    zs = [None]
    ws = [None]
    MWs = [18.015]
    for i in range(1, N):
        Fs.append(None)
        frac = uniform(1e-3/N, 1e-2/N)
        if randint(0, 1):
            zs.append(frac)
            ws.append(None)
        else:
            zs.append(None)
            ws.append(frac)
        MWs.append(uniform(100, 200))

    Fs_calc, zs_calc, ws_calc = solve_flow_composition_mix(Fs, zs, ws, MWs)


    """
    For 4 mass specs, 16 calcs; for 15 mass specs, 225 mw multipliy calcs; numerical solver is
    needed.

    Very fast. Can cache more on the iterations.

    Does not work without at least one mole flow spec. Mass flow specs should be converted to mole flows first.


    Derived from the following sympy code:
    from sympy import *
    F1, F2, F3, F4, F5, F6, F7, F8, F9, MW1, MW2, MW3, MW4, MW5, MW6, MW7, MW8, MW9, Ft = symbols(
        'F1, F2, F3, F4, F5, F6, F7, F8, F9, MW1, MW2, MW3, MW4, MW5, MW6, MW7, MW8, MW9, Ft')

    x2, x3, x4, x5 = symbols('x2, x3, x4, x5')
    w6, w7, w8, w9 = symbols('w6, w7, w8, w9')

    F2 = Ft*x2
    F3 = Ft*x3
    F4 = Ft*x4
    F5 = Ft*x5

    Eq1 = Eq(Ft, F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9)

    den = F1/Ft*MW1 + F2/Ft*MW2 + F3/Ft*MW3 + F4/Ft*MW4 + F5/Ft*MW5 + F6/Ft*MW6 + F7/Ft*MW7 + F8/Ft*MW8 + F9/Ft*MW9

    Eq6 = Eq(w6, F6/Ft*MW6/(den))
    Eq7 = Eq(w7, F7/Ft*MW7/(den))
    Eq8 = Eq(w8, F8/Ft*MW8/(den))
    Eq9 = Eq(w9, F9/Ft*MW9/(den))

    m_ans = solve([Eq6, Eq7, Eq8, Eq9], [F6, F7, F8, F9])
    {F7: -w7*(F1*MW1 + Ft*MW2*x2 + Ft*MW3*x3 + Ft*MW4*x4 + Ft*MW5*x5)/(MW7*(w6 + w7 + w8 + w9 - 1)),
     F6: -w6*(F1*MW1 + Ft*MW2*x2 + Ft*MW3*x3 + Ft*MW4*x4 + Ft*MW5*x5)/(MW6*(w6 + w7 + w8 + w9 - 1)),
     F9: -w9*(F1*MW1 + Ft*MW2*x2 + Ft*MW3*x3 + Ft*MW4*x4 + Ft*MW5*x5)/(MW9*(w6 + w7 + w8 + w9 - 1)),
     F8: -w8*(F1*MW1 + Ft*MW2*x2 + Ft*MW3*x3 + Ft*MW4*x4 + Ft*MW5*x5)/(MW8*(w6 + w7 + w8 + w9 - 1))}
    """


def test_dxs_to_dns():
    ans = dxs_to_dns([-0.0028, -0.00719, -0.00859], [0.7, 0.2, 0.1])
    assert_close1d(ans, [0.001457, -0.0029330000000000003, -0.004333])

    out = [0.0]*3
    ans = dxs_to_dns([-0.0028, -0.00719, -0.00859], [0.7, 0.2, 0.1], out)
    assert_close1d(ans, [0.001457, -0.0029330000000000003, -0.004333])
    assert out is ans

def test_dns_to_dn_partials():
    ans = dns_to_dn_partials([0.001459, -0.002939, -0.004334], -0.0016567)
    assert_close1d(ans, [-0.0001977000000000001, -0.0045957, -0.0059907])

    out_arr = [0.0]*3
    ans = dns_to_dn_partials([0.001459, -0.002939, -0.004334], -0.0016567, out_arr)
    assert_close1d(ans, [-0.0001977000000000001, -0.0045957, -0.0059907])
    assert out_arr is ans

def test_dxs_to_dn_partials():
    ans = dxs_to_dn_partials([-0.0026404, -0.00719, -0.00859], [0.7, 0.2, 0.1], -0.0016567)
    assert_close1d(ans, [-0.00015182, -0.00470142, -0.00610142])

    out = [0.0]*3
    ans = dxs_to_dn_partials([-0.0026404, -0.00719, -0.00859], [0.7, 0.2, 0.1], -0.0016567, out)
    assert_close1d(ans, [-0.00015182, -0.00470142, -0.00610142])
    assert ans is out

def test_dxs_to_dxsn1():
    dxsm1_calc = dxs_to_dxsn1([-2651.3181821109024, -2085.574403592012, -2295.0860830203587])
    dxsm1_expect = [-356.23209909054367, 209.51167942834672]
    assert_close1d(dxsm1_calc, dxsm1_expect, rtol=1e-12)


def test_d2xs_to_d2xsn1():
    test_data = [[-4622.597669746761, 1098.2168754455354, -1014.6525577464832, 210.01951660524486, 974.1035993187835],
                 [1098.2168754455347, -2340.178897128616, -791.9786180652402, 384.16062898478407, 2502.1075882633986],
                 [-1014.6525577464832, -791.9786180652403, -6071.991240755533, 3359.1726518884193, 1823.0297151891225],
                 [210.01951660524486, 384.16062898478407, 3359.1726518884193, -2723.623727200983, 319.21808436722444],
                 [974.1035993187835, 2502.1075882633972, 1823.0297151891225, 319.21808436722444, -5223.4014263476565]]
    d2xsn1_expect = [[-11794.206294731985, -7601.395738484303, -9035.187298602046, -6306.7035934284195],
                     [-7601.395738484303, -12567.79550000307, -10340.517347865418, -7660.566469993495],
                     [-9035.187298602046, -10340.517347865418, -14941.452097481435, -4006.4765740155844],
                     [-6306.7035934284195, -7660.566469993495, -4006.4765740155844, -8585.461322283089]]
    assert_close2d(d2xs_to_d2xsn1(test_data), d2xsn1_expect, rtol=1e-12)

def test_d2ns_to_dn2_partials():
    d2ns = [[0.152, 0.08, 0.547], [0.08, 0.674, 0.729], [0.547, 0.729, 0.131]]
    res = d2ns_to_dn2_partials(d2ns, [20.0, .124, 900.52])
    res_expect = [[40.152, 20.203999999999997, 921.067], [20.204, 0.922, 901.3729999999999], [921.067, 901.373, 1801.1709999999998]]
    assert_close2d(res, res_expect)


def test_mixing_power():
    zs = [.1, .4, .5]
    props = [.01, .015, .025]
    assert_close(mixing_power(zs, props, r=-4), 0.015110570481191988, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-3), 0.01586305104547743, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-2), 0.016718346377260586, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-1), 0.01764705882352941, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=1), 0.0195, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=2), 0.020310096011589902, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=3), 0.021001133725640605, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=4), 0.021572268349078574, rtol=1e-13)

    # Creating these was a lot of work
    assert_close(mixing_power(zs, props, r=1.0/6.0), 0.01875088400005579, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=5.0/6.0), 0.019354727234672513, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=5.0/3.0), 0.020052643726129516, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-1/6), 0.018438750245785035, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-5/6), 0.017805472008199413, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-5/3), 0.01702198311393396, rtol=1e-13)

    assert_close(mixing_power(zs, props, r=-.25), 0.018360027162617934, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-.5), 0.018122801833129563, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-.75), 0.01788479898151489, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-1.5), 0.017176405936050576, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.25), 0.01882804734220324, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.5), 0.019056901419612297, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.75), 0.01928114152959784, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=1.5), 0.005541354936265382, rtol=1e-13)
    
    assert_close(mixing_power(zs, props, r=-4/3.0), 0.017332239251910703, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-2/3), 0.017964155302880425, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-1/3), 0.01828110264599271, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=4/3), 0.0197823134832162, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=2/3), 0.01920695698168168, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=1/3), 0.01890479323567482, rtol=1e-13)

    assert_close(mixing_power(zs, props, r=0.125), 0.018712156392543106, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.375), 0.01894299926200097, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.625), 0.019169648589127595, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=0.875), 0.019391287384440178, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-0.125), 0.018478024469926933, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-0.375), 0.01824157672809132, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-0.625), 0.018003832442176182, rtol=1e-13)
    assert_close(mixing_power(zs, props, r=-0.875), 0.017765831497664065, rtol=1e-13)

    # generic case
    assert_close(mixing_power(zs, props, r=3.3423), 0.02120966068692269, rtol=1e-13)

    

def test_hash_any_primitive():
    from chemicals.utils import hash_any_primitive
    a = {'a': 1, 'b': 2}
    b = {'b': 2, 'a': 1}
    assert hash_any_primitive(a) == hash_any_primitive(b)


    assert hash_any_primitive([]) == hash_any_primitive([])

    assert hash_any_primitive([1,2]) == hash_any_primitive([1,2])

    assert not hash_any_primitive([1,2]) == hash_any_primitive([1,2,3])

    assert not hash_any_primitive([1,2,4]) == hash_any_primitive([1,2,3])


def test_radius_of_gyration():
    assert_close(radius_of_gyration(MW=32.00452, planar=False, A=6.345205205562681e-47, B=3.2663291891213418e-46, C=3.4321304373822523e-46),
                 1.507895671067054e-10, rtol=1e-13)

    calc = radius_of_gyration(A=9.972482951577335e-47, B=1.941784104125803e-45, planar=True, MW=111.01, C=0)
    assert_close(calc, 4.8859099776532495e-11, rtol=1e-13)


def test_molar_velocity():
    assert_close(molar_velocity_to_velocity(46.537593457316525, 67.152), 179.5868138460819, rtol=1e-12)

    assert_close(velocity_to_molar_velocity(179.5868138460819, 67.152), 46.537593457316525, rtol=1e-12)

def test_ms_to_ns():
    ans = ms_to_ns([0.00033915944387075575, 0.005163034654211768, 0.002485711001895458],[18.01528, 16.04246, 142.28168])
    assert_close1d(ans, [0.018826209965693327, 0.3218355946788565, 0.0174703517831351])


def test_ns_to_ms():
    ans = ns_to_ms([0.018826209965693327, 0.3218355946788565, 0.0174703517831351], [18.01528, 16.04246, 142.28168])
    assert_close1d(ans, [0.00033915944387075575, 0.005163034654211768, 0.002485711001895458])


def test_ns_to_Qls():
    ans = ns_to_Qls([0.018826209965693327, 0.3218355946788565, 0.0174703517831351], [1.8087205105724903e-05, 5.858784737690099e-05, 0.00019580845677748954])
    assert_close1d(ans, [3.405135210129374e-07, 1.885565470149901e-05, 3.4208426220155466e-06])

def test_Qls_to_ns():
    ans = Qls_to_ns([3.405135210129374e-07, 1.885565470149901e-05, 3.4208426220155466e-06] , [1.8087205105724903e-05, 5.858784737690099e-05, 0.00019580845677748954])
    assert_close1d(ans, [0.018826209965693327, 0.3218355946788565, 0.0174703517831351])


def test_ms_to_Qls():
    ans = ms_to_Qls([0.00033915944387075575, 0.005163034654211768, 0.002485711001895458], [18.01528, 16.04246, 142.28168], [1.8087205105724903e-05, 5.858784737690099e-05, 0.00019580845677748954])
    assert_close1d(ans, [3.405135210129374e-07, 1.885565470149901e-05, 3.4208426220155466e-06])
    # Qls_to_ms, ms_to_Qls

def test_Qls_to_ms():
    ans = Qls_to_ms([3.405135210129374e-07, 1.885565470149901e-05, 3.4208426220155466e-06], [18.01528, 16.04246, 142.28168], [1.8087205105724903e-05, 5.858784737690099e-05, 0.00019580845677748954])
    assert_close1d(ans, [0.00033915944387075575, 0.005163034654211768, 0.002485711001895458])


def test_Baume_light_roundtrip():
    # Test rho conversion round trip
    rho = 820.0
    Be = rho_to_Baume_light(rho)
    rho2 = Baume_light_to_rho(Be)
    assert_close(rho, rho2)
    
    # Test SG conversion round trip
    SG = 0.8209
    Be = SG_to_Baume_light(SG) 
    SG2 = Baume_light_to_SG(Be)
    assert_close(SG, SG2)

def test_Baume_heavy_roundtrip():
    # Test rho conversion round trip
    rho = 1200.0
    Be = rho_to_Baume_heavy(rho)
    rho2 = Baume_heavy_to_rho(Be)
    assert_close(rho, rho2)
    
    # Test SG conversion round trip
    SG = 1.2012
    Be = SG_to_Baume_heavy(SG)
    SG2 = Baume_heavy_to_SG(Be)
    assert_close(SG, SG2)