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
from fluids.constants import R as _R
from fluids.numerics import assert_close, assert_close1d, assert_close2d, assert_close3d, derivative, hessian, jacobian, linspace

from chemicals import rho_to_Vm
from chemicals.virial import (
    B_from_Z,
    B_to_Z,
    BVirial_Abbott,
    BVirial_Abbott_fast,
    BVirial_Abbott_mat,
    BVirial_Abbott_vec,
    BVirial_Meng,
    BVirial_Meng_mat,
    BVirial_Meng_vec,
    BVirial_mixture,
    BVirial_Oconnell_Prausnitz,
    BVirial_Oconnell_Prausnitz_mat,
    BVirial_Oconnell_Prausnitz_vec,
    BVirial_Pitzer_Curl,
    BVirial_Pitzer_Curl_fast,
    BVirial_Pitzer_Curl_mat,
    BVirial_Pitzer_Curl_vec,
    BVirial_Tsonopoulos,
    BVirial_Tsonopoulos_extended,
    BVirial_Tsonopoulos_extended_fast,
    BVirial_Tsonopoulos_extended_mat,
    BVirial_Tsonopoulos_extended_vec,
    BVirial_Tsonopoulos_fast,
    BVirial_Tsonopoulos_mat,
    BVirial_Tsonopoulos_vec,
    BVirial_Xiang,
    BVirial_Xiang_mat,
    BVirial_Xiang_vec,
    CVirial_Liu_Xiang,
    CVirial_Liu_Xiang_mat,
    CVirial_Liu_Xiang_vec,
    CVirial_mixture_Orentlicher_Prausnitz,
    CVirial_Orbey_Vera,
    CVirial_Orbey_Vera_mat,
    CVirial_Orbey_Vera_vec,
    Lee_Kesler_virial_CSP_Vcijs,
    Meng_Duan_2005_virial_CSP_kijs,
    Meng_virial_a,
    Tarakad_Danner_virial_CSP_kijs,
    Tarakad_Danner_virial_CSP_omegaijs,
    Tarakad_Danner_virial_CSP_Pcijs,
    Tarakad_Danner_virial_CSP_Tcijs,
    Z_from_virial_density_form,
    Z_from_virial_pressure_form,
    d2BVirial_mixture_dzizjs,
    d2CVirial_mixture_dT2_Orentlicher_Prausnitz,
    d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs,
    d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs,
    d2V_dzizjs_virial,
    d3BVirial_mixture_dzizjzks,
    d3CVirial_mixture_dT3_Orentlicher_Prausnitz,
    d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks,
    dBVirial_mixture_dzs,
    dCVirial_mixture_dT_Orentlicher_Prausnitz,
    dCVirial_mixture_Orentlicher_Prausnitz_dzs,
    dV_dzs_virial,
)


def test_B_To_Z():
    Z_calc = B_to_Z(-0.0015, 300.0, 1E5)
    assert_close(Z_calc, 0.9398638020957176)


def test_B_from_Z():
    B_calc = B_from_Z(0.94, 300.0, 1E5)
    assert_close(B_calc, -0.0014966032712675846)


def test_Z_from_virial_density_form():
    Z_calc = Z_from_virial_density_form(300.0, 122057.233762653, (1E-4, 1E-5, 1E-6, 1E-7))
    assert_close(Z_calc, 1.2843494052609186)

    Z_calc = Z_from_virial_density_form(300, 102031.881198762, (1e-4, 1e-5, 1e-6))
    assert_close(Z_calc, 1.0736323841544937)

    Z_calc = Z_from_virial_density_form(300, 96775.8831504971, (1e-4, 1e-5))
    assert_close(Z_calc, 1.018326089216066)

    Z_calc = Z_from_virial_density_form(300, 95396.3561037084, (1e-4,))
    assert_close(Z_calc,  1.003809998713499)

    assert_close(1, Z_from_virial_density_form(300, 95396.3561037084))

    """B-only solution, derived as follows:

    >>> B, C, D, E = symbols('B, C, D, E')
    >>> P, V, Z, R, T = symbols('P, V, Z, R, T', positive=True, real=True, nonzero=True)
    >>> rho = 1/V
    >>> to_slv = Eq(P*V/R/T, 1 + B*rho)
    >>> slns = solve(to_slv, V)
    >>> simplify(slns[1]*P/R/T)
    1/2 + sqrt(4*B*P + R*T)/(2*sqrt(R)*sqrt(T))

    To check this, simply disable the if statement and allow the numerical
    algorithm to run.
    """


def test_Z_from_virial_pressure_form():
    Z_calc = Z_from_virial_pressure_form(102919.99946855308, (4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19))
    assert_close(Z_calc, 1.00283753944)

    Z_calc = Z_from_virial_pressure_form(102847.17619188508, (4.032286555169439e-09, 1.6197059494442215e-13))
    assert_close(Z_calc, 1.00212796)

    Z_calc = Z_from_virial_pressure_form(102671.27455742132, (4.032286555169439e-09,))
    assert_close(Z_calc, 1.000414)

    Z_calc = Z_calc = Z_from_virial_pressure_form(102671.27455742132)
    assert_close(Z_calc, 1)


def test_BVirial_Pitzer_Curl():
    # doctest
    B = BVirial_Pitzer_Curl(510., 425.2, 38E5, 0.193)
    assert_close(B, -0.00020845362479301725)

    with pytest.raises(Exception):
        BVirial_Pitzer_Curl(510., 425.2, 38E5, 0.193, order=-3)

    values_expect = [-0.00020845362479301725, 1.0653775169998656e-06, -5.795710171294467e-09, 4.513533043400151e-11, -0.437891506790894, 8.720086532349054]
    values = [BVirial_Pitzer_Curl(510., 425.2, 38E5, 0.193, order) for order in [0, 1, 2, 3, -1, -2]]
    assert_close1d(values, values_expect, rtol=1e-12)



@pytest.mark.slow
@pytest.mark.sympy
def test_BVirial_Pitzer_Curl_calculus():
    from sympy import Rational, diff, integrate, lambdify, symbols

    # Derivatives check
    # Uses SymPy
    T, Tc, Pc, omega, R = symbols('T, Tc, Pc, omega, R')
    Tr = T/Tc
    B0 = Rational(1445,10000) - Rational(33,100)/Tr - Rational(1385,10000)/Tr**2 - Rational(121,10000)/Tr**3
    B1 = Rational(73,1000) + Rational(46,100)/Tr - Rational(1,2)/Tr**2 - Rational(97,1000)/Tr**3 - Rational(73,10000)/Tr**8

    # Note: scipy.misc.derivative was attempted, but found to given too
    # incorrect of derivatives for higher orders, so there is no reasons to
    # implement it. Plus, no uses have yet been found for the higher
    # derivatives/integrals. For integrals, there is no way to get the
    # indefinite integral.

    # Test points in vector form for order 1, 2, and 3
    # Use lambdify for fast evaluation
    pts = 3 # points total = pts^4
    _Ts = linspace(5,500,pts)
    _Tcs = linspace(501,900,pts)
    _Pcs = linspace(2E5, 1E7,pts)
    _omegas = linspace(-1, 10,pts)

    _Ts, _Tcs, _Pcs, _omegas = np.meshgrid(_Ts, _Tcs, _Pcs, _omegas)
    _Ts, _Tcs, _Pcs, _omegas = _Ts.ravel(), _Tcs.ravel(), _Pcs.ravel(), _omegas.ravel()
    for order in range(1,4):
        B0c = diff(B0, T, order)
        B1c = diff(B1, T, order)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)
        Bcalc2 = BVirial_Pitzer_Curl(_Ts, _Tcs, _Pcs, _omegas, order)
        assert_close1d(Bcalcs, Bcalc2)


    # Check integrals using SymPy:
    for order in range(-2, 0):
        B0c = B0
        B1c = B1
        for i in range(abs(order)):
            B0c = integrate(B0c, T)
            B1c = integrate(B1c, T)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)
        Bcalc2 = [BVirial_Pitzer_Curl(T2, Tc2, Pc2, omega2, order) for T2, Tc2, Pc2, omega2 in zip(_Ts, _Tcs, _Pcs, _omegas)]
        assert_close1d(Bcalcs, Bcalc2)

def test_BVirial_Abbott_calculus():
    B = BVirial_Abbott(510., 425.2, 38E5, 0.193)
    assert_close(B, -0.00020570185009564064)

    with pytest.raises(Exception):
        BVirial_Abbott(510., 425.2, 38E5, 0.193, order=-3)

    values_expect = [-0.00020570185009564064, 1.0392492946983827e-06, -5.9022336392448295e-09, 4.782227646523899e-11, 0.30386992442862953, 330.826226911517]
    values = [BVirial_Abbott(510., 425.2, 38E5, 0.193, order) for order in [0, 1, 2, 3, -1, -2]]
    assert_close1d(values, values_expect, rtol=1e-12)

@pytest.mark.slow
@pytest.mark.sympy
def test_BVirial_Abbott_calculus_sympy():
    from sympy import diff, integrate, lambdify, symbols

    B = BVirial_Abbott(510., 425.2, 38E5, 0.193)
    assert_close(B, -0.00020570185009564064)

    with pytest.raises(Exception):
        BVirial_Abbott(510., 425.2, 38E5, 0.193, order=-3)

    T, Tc, Pc, omega, R = symbols('T, Tc, Pc, omega, R')
    Tr = T/Tc
    B0 = 0.083 - 0.422/Tr**1.6
    B1 = 0.139 - 0.172/Tr**4.2

    # Test points in vector form for order 1, 2, and 3
    # Use lambdify for fast evaluation
    pts = 3
    _Ts = np.linspace(5,500,pts)
    _Tcs = np.linspace(501,900,pts)
    _Pcs = np.linspace(2E5, 1E7,pts)
    _omegas = np.linspace(-1, 10,pts)
    _Ts, _Tcs, _Pcs, _omegas = np.meshgrid(_Ts, _Tcs, _Pcs, _omegas)
    _Ts, _Tcs, _Pcs, _omegas = _Ts.ravel(), _Tcs.ravel(), _Pcs.ravel(), _omegas.ravel()


    for order in range(1,4):
        B0c = diff(B0, T, order)
        B1c = diff(B1, T, order)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)
        Bcalc2 = BVirial_Abbott(_Ts, _Tcs, _Pcs, _omegas, order)
        assert_close1d(Bcalcs, Bcalc2)


    # Check integrals using SymPy:
    for order in range(-2, 0):
        B0c = B0
        B1c = B1
        for i in range(abs(order)):
            B0c = integrate(B0c, T)
            B1c = integrate(B1c, T)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)

        Bcalc2 = [BVirial_Abbott(T2, Tc2, Pc2, omega2, order) for T2, Tc2, Pc2, omega2 in zip(_Ts, _Tcs, _Pcs, _omegas)]
        assert_close1d(Bcalcs, Bcalc2)

def test_BVirial_Tsonopoulos():
    B = BVirial_Tsonopoulos(510., 425.2, 38E5, 0.193)
    assert_close(B, -0.00020935295404416802)

    with pytest.raises(Exception):
        BVirial_Tsonopoulos(510., 425.2, 38E5, 0.193, order=-3)


    values_expect = [-0.00020935295404416802, 9.95742355603791e-07, -5.542344657946387e-09, 4.570351609785339e-11, -0.7019279964346002, -257.84756571017147]
    values = [BVirial_Tsonopoulos(510., 425.2, 38E5, 0.193, order) for order in [0, 1, 2, 3, -1, -2]]
    assert_close1d(values_expect, values)

@pytest.mark.slow
@pytest.mark.sympy
def test_BVirial_Tsonopoulos_calculus():
    from sympy import Rational, diff, integrate, lambdify, symbols

    T, Tc, Pc, omega, R = symbols('T, Tc, Pc, omega, R')
    Tr = T/Tc
    B0 = Rational(1445, 10000) - Rational(33,100)/Tr - Rational(1385,10000)/Tr**2 - Rational(121,10000)/Tr**3 - Rational(607,1000000)/Tr**8
    B1 = Rational(637,10000) + Rational(331,1000)/Tr**2 - Rational(423,1000)/Tr**3 - Rational(8,1000)/Tr**8

    # Test points in vector form for order 1, 2, and 3
    # Use lambdify for fast evaluation
    pts = 3
    _Ts = np.linspace(5,500,pts)
    _Tcs = np.linspace(501,900,pts)
    _Pcs = np.linspace(2E5, 1E7,pts)
    _omegas = np.linspace(-1, 10,pts)
    _Ts, _Tcs, _Pcs, _omegas = np.meshgrid(_Ts, _Tcs, _Pcs, _omegas)
    _Ts, _Tcs, _Pcs, _omegas = _Ts.ravel(), _Tcs.ravel(), _Pcs.ravel(), _omegas.ravel()


    for order in range(1,4):
        B0c = diff(B0, T, order)
        B1c = diff(B1, T, order)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)
        Bcalc2 = BVirial_Tsonopoulos(_Ts, _Tcs, _Pcs, _omegas, order)
        assert_close1d(Bcalcs, Bcalc2)

    # Check integrals using SymPy:
    for order in range(-2, 0):
        B0c = B0
        B1c = B1
        for i in range(abs(order)):
            B0c = integrate(B0c, T)
            B1c = integrate(B1c, T)
        Br = B0c + omega*B1c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)

        Bcalc2 = [BVirial_Tsonopoulos(T2, Tc2, Pc2, omega2, order) for T2, Tc2, Pc2, omega2 in zip(_Ts, _Tcs, _Pcs, _omegas)]
        assert_close1d(Bcalcs, Bcalc2)

def test_BVirial_Tsonopoulos_extended():
    B = BVirial_Tsonopoulos_extended(510., 425.2, 38E5, 0.193, species_type='normal', dipole=0)
    assert_close(B, -0.00020935295404416802)

    B = BVirial_Tsonopoulos_extended_fast(510., 425.2, 38E5, 0.193, a=0, b=0)[0]
    assert_close(B, -0.00020935295404416802)
    for i in range(1, 4):
        assert_close(BVirial_Tsonopoulos_extended(510., 425.2, 38E5, 0.193, a=0, b=0, order=i),
                     BVirial_Tsonopoulos_extended_fast(510., 425.2, 38E5, 0.193, a=0, b=0)[i], rtol=1e-13)


    B = BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469)
    assert_close(B, -9.679718337596426e-05)

    with pytest.raises(Exception):
        BVirial_Tsonopoulos_extended(510., 425.2, 38E5, 0.193, order=-3)

    # Test all of the different types
    types = ['simple', 'normal', 'methyl alcohol', 'water', 'ketone',
    'aldehyde', 'alkyl nitrile', 'ether', 'carboxylic acid', 'ester', 'carboxylic acid',
    'ester', 'alkyl halide', 'mercaptan', 'sulfide', 'disulfide', 'alkanol']

    Bs_calc = [BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608,
                                            a=0, b=0, species_type=i, dipole=0.1) for i in types]
    Bs = [-9.00253249139901e-05, -9.00253249139901e-05, -8.136808332317606e-05, -9.232253763245037e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.00558374295638e-05, -9.003498498098181e-05, -9.003498498098181e-05, -9.003498498098181e-05, -9.003498498098181e-05, -7.331249596682434e-05]
    assert_close1d(Bs_calc, Bs)

    # Test the fast call - with a group specified
    calc0 = BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608,
                                            a=0, b=0, species_type='ketone', dipole=0.1)

    calc1 = BVirial_Tsonopoulos_extended_fast(430., 405.65, 11.28E6, 0.252608,
                                            a=-0.00014477824238067583, b=0)[0]
    assert_close(calc0, calc1, rtol=1e-10)

    # Test compare a and b both being specified

    calc0 = BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608,
                                            a=-1e-5, b=-1e-7)

    calc1 = BVirial_Tsonopoulos_extended_fast(430., 405.65, 11.28E6, 0.252608,
                                            a=-1e-5, b=-1e-7)[0]
    assert_close(calc0, calc1, rtol=1e-10)



@pytest.mark.slow
@pytest.mark.sympy
def test_BVirial_Tsonopoulos_extended_calculus():
    from sympy import Rational, diff, integrate, lambdify, symbols

    # Use lambdify for fast evaluation
    pts = 3
    _Ts = np.linspace(5,500,pts)
    _Tcs = np.linspace(501,900,pts)
    _Pcs = np.linspace(2E5, 1E7,pts)
    _omegas = np.linspace(-1, 10,pts)
    _Ts, _Tcs, _Pcs, _omegas = np.meshgrid(_Ts, _Tcs, _Pcs, _omegas)
    _Ts, _Tcs, _Pcs, _omegas = _Ts.ravel(), _Tcs.ravel(), _Pcs.ravel(), _omegas.ravel()

    T, Tc, Pc, omega, R = symbols('T, Tc, Pc, omega, R')
    Tr = T/Tc
    # Attempted to optimize things - the integrate expressions for B)c and B1c are the slowest.
    B0 = Rational(1445, 10000) - Rational(33,100)/Tr - Rational(1385,10000)/Tr**2 - Rational(121,10000)/Tr**3 - Rational(607,1000000)/Tr**8
    B1 = Rational(637,10000) + Rational(331,1000)/Tr**2 - Rational(423,1000)/Tr**3 - Rational(8,1000)/Tr**8
    B2 = 1/Tr**6
    B3 = -1/Tr**8

    a = 0.1
    b = 0.2

    for order in range(1,4):
        B0c = diff(B0, T, order)
        B1c = diff(B1, T, order)
        B2c = diff(B2, T, order)
        B3c = diff(B3, T, order)

        Br = B0c + omega*B1c + a*B2c + b*B3c
        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)
        Bcalc2 = BVirial_Tsonopoulos_extended(_Ts, _Tcs, _Pcs, _omegas, order=order, a=a, b=b)
        assert_close1d(Bcalcs, Bcalc2)


    # Check integrals using SymPy:
    for order in range(-2, 0):
        B0c = B0
        B1c = B1
        B2c = B2
        B3c = B3
        for i in range(abs(order)):
            B0c = integrate(B0c, T)
            B1c = integrate(B1c, T)
            B2c = integrate(B2c, T)
            B3c = integrate(B3c, T)
        Br = B0c + omega*B1c + a*B2c + b*B3c

        BVirial = (Br*R*Tc/Pc).subs(R, _R)
        f = lambdify((T, Tc, Pc, omega), BVirial, "numpy")

        Bcalcs = f(_Ts, _Tcs, _Pcs, _omegas)


        Bcalc2 = [BVirial_Tsonopoulos_extended(T2, Tc2, Pc2, omega2, a=a, b=b, order=order) for T2, Tc2, Pc2, omega2 in zip(_Ts, _Tcs, _Pcs, _omegas)]
        assert_close1d(Bcalcs, Bcalc2)


def test_CVirial_Orbey_Vera():
    Tc, Pc, omega = 568.7, 2490000.0, 0.394
    T = 0.866*Tc
    expect = (1.9314651915020253e-07, 3.3074098632105e-11, -2.2771746296770828e-11, 7.437013303758471e-13)
    calc = CVirial_Orbey_Vera(T, Tc, Pc, omega)
    assert_close1d(expect, calc, rtol=1e-14)
    assert_close(derivative(lambda T: CVirial_Orbey_Vera(T, Tc, Pc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: CVirial_Orbey_Vera(T, Tc, Pc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: CVirial_Orbey_Vera(T, Tc, Pc, omega)[2], T, dx=T*1e-6), expect[3])

    # Vector call with out memory savings
    vec_call = CVirial_Orbey_Vera_vec(T, [Tc], [Pc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Cs_out = [0]
    dCs_out = [0]
    d2Cs_out = [0]
    d3Cs_out = [0]

    # vector call with memory savings
    CVirial_Orbey_Vera_vec(T, [Tc], [Pc],  [omega], Cs_out, dCs_out, d2Cs_out, d3Cs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Cs_out, dCs_out, d2Cs_out, d3Cs_out], rtol=1e-13)

    # matrix call
    mat_call = CVirial_Orbey_Vera_mat(T, [[Tc]], [[Pc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Cs_out = [[0]]
    dCs_out = [[0]]
    d2Cs_out = [[0]]
    d3Cs_out = [[0]]
    CVirial_Orbey_Vera_mat(T, [[Tc]], [[Pc]],  [[omega]], Cs_out, dCs_out, d2Cs_out, d3Cs_out)
    assert_close3d(expect_mat, [Cs_out, dCs_out, d2Cs_out, d3Cs_out], rtol=1e-13)




def test_CVirial_Liu_Xiang():
    T = 300
    Tc = 568.7
    Pc = 2490000.0
    Vc = 0.000492
    omega = 0.394
    expect = (-1.3720744771519268e-05, 5.131512522084089e-07, -2.031692047127711e-08, 8.713867012492766e-10)
    C = CVirial_Liu_Xiang(T=T, Tc=Tc, Pc=Pc, Vc=Vc, omega=omega)[0]
    assert_close(C, expect[0])

    assert_close(derivative(lambda T: CVirial_Liu_Xiang(T, Tc, Pc, Vc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: CVirial_Liu_Xiang(T, Tc, Pc, Vc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: CVirial_Liu_Xiang(T, Tc, Pc, Vc, omega)[2], T, dx=T*1e-6), expect[3])

    # Vector call with out memory savings
    vec_call = CVirial_Liu_Xiang_vec(T, [Tc], [Pc], [Vc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Cs_out = [0]
    dCs_out = [0]
    d2Cs_out = [0]
    d3Cs_out = [0]

    # vector call with memory savings
    CVirial_Liu_Xiang_vec(T, [Tc], [Pc], [Vc], [omega], Cs_out, dCs_out, d2Cs_out, d3Cs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Cs_out, dCs_out, d2Cs_out, d3Cs_out], rtol=1e-13)

    # matrix call
    mat_call = CVirial_Liu_Xiang_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Cs_out = [[0]]
    dCs_out = [[0]]
    d2Cs_out = [[0]]
    d3Cs_out = [[0]]
    CVirial_Liu_Xiang_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]], Cs_out, dCs_out, d2Cs_out, d3Cs_out)
    assert_close3d(expect_mat, [Cs_out, dCs_out, d2Cs_out, d3Cs_out], rtol=1e-13)


    # POint with a graph
    Tc = 647.10
    Pc = 22050e3
    omega = 0.344
    MW = 18.015
    rhoc_mass = 325.0
    Vc = rho_to_Vm(MW=MW, rho=rhoc_mass)

    graph_point = CVirial_Liu_Xiang(T=0.6*Tc, Tc=Tc, Pc=Pc, Vc=Vc, omega=omega)[0]/Vc**2
    assert_close(graph_point, -48.10297670037914)


def test_BVirial_Pitzer_Curl_fast():
    Tc = 508.1
    Pc = 4696e3
    T = 0.7*Tc
    omega = 0.308
    expect = [-0.00077399438330036, 5.418010746475681e-06, -5.971092529726755e-08, 1.0270254054853065e-09]
    assert_close(BVirial_Pitzer_Curl_fast(T, Tc, Pc, omega)[0], -0.00077399438330036, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Pitzer_Curl_fast(T, Tc, Pc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Pitzer_Curl_fast(T, Tc, Pc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Pitzer_Curl_fast(T, Tc, Pc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Pitzer_Curl_vec(T, [Tc], [Pc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Pitzer_Curl_vec(T, [Tc], [Pc], [omega], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Pitzer_Curl_mat(T, [[Tc]], [[Pc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Pitzer_Curl_mat(T, [[Tc]], [[Pc]], [[omega]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_BVirial_Abbott_fast():
    Tc = 508.1
    Pc = 4696e3
    T = 0.7*Tc
    omega = 0.308
    expect = [-0.0007717412887871863, 5.539165302002657e-06, -5.889351243581266e-08, 8.651388154761998e-10]
    assert_close(BVirial_Abbott_fast(T, Tc, Pc, omega)[0], -0.0007717412887871863, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Abbott_fast(T, Tc, Pc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Abbott_fast(T, Tc, Pc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Abbott_fast(T, Tc, Pc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Abbott_vec(T, [Tc], [Pc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Abbott_vec(T, [Tc], [Pc], [omega], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Abbott_mat(T, [[Tc]], [[Pc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Abbott_mat(T, [[Tc]], [[Pc]], [[omega]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_BVirial_Tsonopoulos_fast():
    Tc = 508.1
    Pc = 4696e3
    T = 0.7*Tc
    omega = 0.308
    expect = [-0.0007649313249027217, 5.797597679149657e-06, -7.258950480792238e-08, 1.3572606703229816e-09]
    assert_close(BVirial_Tsonopoulos_fast(T, Tc, Pc, omega)[0], -0.0007649313249027217, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_fast(T, Tc, Pc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_fast(T, Tc, Pc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_fast(T, Tc, Pc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Tsonopoulos_vec(T, [Tc], [Pc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Tsonopoulos_vec(T, [Tc], [Pc], [omega], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Tsonopoulos_mat(T, [[Tc]], [[Pc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Tsonopoulos_mat(T, [[Tc]], [[Pc]], [[omega]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_BVirial_Xiang():
    # Acetone, figure 5 acetone
    MW = 58.080
    Tc = 508.1
    Pc = 4696e3
    rhoc_mass = 268
    omega = 0.308
    T = 0.7*Tc
    Vc = rho_to_Vm(MW=MW, rho=rhoc_mass)
    # ~-4.5 expected if ans not multiplied by Vc
    expect = [-0.0008878022797789029, 7.467113021775213e-06, -9.800480205221981e-08, 1.7837228334430914e-09]
    assert_close(BVirial_Xiang(T, Tc, Pc, Vc, omega)[0], -0.0008878022797789028, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Xiang(T, Tc, Pc, Vc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Xiang(T, Tc, Pc, Vc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Xiang(T, Tc, Pc, Vc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Xiang_vec(T, [Tc], [Pc], [Vc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Xiang_vec(T, [Tc], [Pc], [Vc], [omega], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Xiang_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Xiang_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_Tarakad_Danner_virial_CSP_kijs():

    expect = [[0.0, 0.016463320918394864, 0.048781060667435705, 0.2243198001812905],
     [0.016463320918394864, 0.0, 0.11579684127861978, 0.1338961298780077],
     [0.048781060667435705, 0.11579684127861978, 0.0, 0.4001799396649175],
     [0.2243198001812905, 0.1338961298780077, 0.4001799396649175, 0.0]]


    Vcs = [0.000168, 0.000316, 5.6e-05, 0.002055]
    ans = Tarakad_Danner_virial_CSP_kijs(Vcs)
    assert_close2d(ans, expect, rtol=1e-13)

def test_Tarakad_Danner_virial_CSP_Tcijs():
    Vcs = [0.000168, 0.000316, 5.6e-05, 0.002055]
    Tcs = [514.0, 591.75, 647.14, 843.0]
    kijs = Tarakad_Danner_virial_CSP_kijs(Vcs)
    Tcijs = Tarakad_Danner_virial_CSP_Tcijs(Tcs=Tcs, kijs=kijs)
    Tcijs_expect = [[514.0, 542.4269432446305, 548.606779975124, 510.5967574676473],
     [542.4269432446305, 591.7500000000001, 547.1675300600267, 611.7203098423653],
     [548.606779975124, 547.1675300600267, 647.14, 443.0307753796196],
     [510.5967574676473, 611.7203098423653, 443.0307753796196, 842.9999999999999]]
    assert_close2d(Tcijs, Tcijs_expect, rtol=1e-12)

def test_Tarakad_Danner_virial_CSP_Pcijs():
    kijs = Tarakad_Danner_virial_CSP_kijs(Vcs=[0.000168, 0.000316])
    Tcijs = Tarakad_Danner_virial_CSP_Tcijs(Tcs=[514.0, 591.75], kijs=kijs)
    Pcijs = Tarakad_Danner_virial_CSP_Pcijs(Tcs=[514.0, 591.75], Pcs=[6137000.0, 4108000.0], Vcs=[0.000168, 0.000316], Tcijs=Tcijs)
    expect = [[6136999.999999997, 4861936.434873204], [4861936.434873204, 4107999.9999999995]]
    assert_close2d(Pcijs, expect, rtol=1e-13)

def test_Tarakad_Danner_virial_CSP_omegaijs():
    omegaijs = Tarakad_Danner_virial_CSP_omegaijs([0.635, 0.257, 0.344, 1.26])
    omegaijs_expect = [[0.635, 0.446, 0.4895, 0.9475], [0.446, 0.257, 0.3005, 0.7585], [0.4895, 0.3005, 0.344, 0.802], [0.9475, 0.7585, 0.802, 1.26]]

    assert_close2d(omegaijs, omegaijs_expect, rtol=1e-12)


def test_Meng_Duan_2005_virial_CSP_kijs():
    CASs = ['74-82-8', '74-84-0', '124-38-9', '7727-37-9', '7439-89-6']
    atomss = [{'C': 1, 'H': 4}, {'C': 2, 'H': 6}, {'C': 1, 'O': 2}, {'N': 2}, {'Fe': 1}]
    kijs = Meng_Duan_2005_virial_CSP_kijs(CASs=CASs, atomss=atomss)
    kijs_expect = [[0.0, 0.0014070591327669385, 0.04313694538361394, 0.024878043016556484, 0.0],
                   [0.0014070591327669385, 0.0, 0.08607516737052616, 0.0496414777972359, 0.0],
                   [0.04313694538361394, 0.08607516737052616, 0.0, 0.0, 0.0],
                   [0.024878043016556484, 0.0496414777972359, 0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0, 0.0, 0.0]]

    assert_close2d(kijs, kijs_expect, rtol=1e-13)

def test_Lee_Kesler_virial_CSP_Vcijs():
    Vcs = [0.000168, 0.000316, 5.6e-05, 0.002055]
    ans = Lee_Kesler_virial_CSP_Vcijs(Vcs)
    expect = [[0.00016800000000000004, 0.00023426511495004188, 0.00010196900126054562, 0.0007574916472147805], [0.00023426511495004188, 0.00031600000000000015, 0.00015044767921762743, 0.0009304209382011844], [0.00010196900126054562, 0.00015044767921762743, 5.6000000000000047e-05, 0.0005655603315853534], [0.0007574916472147805, 0.0009304209382011844, 0.0005655603315853534, 0.002055000000000001]]
    assert_close2d(expect, ans, rtol=1e-13)

def test_Meng_virial_a():
    calc = Meng_virial_a(317.4, 5870000.0, 1.85, haloalkane=True)
    assert_close(calc, -0.04493829786760545, rtol=1e-13)

    calc = Meng_virial_a(514.0, 6137000.0, 1.44, haloalkane=False)
    assert_close(calc, -0.006378416625935997, rtol=1e-13)

def test_BVirial_Meng():
    T = 388.26
    Tc = 647.1
    Pc = 22050000.0
    Vc = 5.543076e-05
    omega = 0.344
    a = 0.0
    calc = BVirial_Meng(T, Tc, Pc, Vc, omega, a)
    expect = (-0.00032436028497558636, 2.470038900338557e-06, -3.132003987118147e-08, 5.776332655071256e-10)
    assert_close1d(calc, expect, rtol=1e-13)

    assert_close(BVirial_Meng(T, Tc, Pc, Vc, omega)[0], -0.00032436028497558636, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Meng(T, Tc, Pc, Vc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Meng(T, Tc, Pc, Vc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Meng(T, Tc, Pc, Vc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Meng_vec(T, [Tc], [Pc], [Vc], [omega], [a])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Meng_vec(T, [Tc], [Pc], [Vc], [omega], [a], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Meng_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]], [[a]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Meng_mat(T, [[Tc]], [[Pc]], [[Vc]], [[omega]], [[a]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_BVirial_Tsonopoulos_extended_vec():
    T = 388.26
    Tc = 647.1
    Pc = 22050000.0
    omega = 0.344
    a = 1e-7
    b = 1e-12
    calc = BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a, b)
    expect = (-0.0003371378639247552, 2.8128391774419057e-06, -3.992508039858756e-08, 8.034145734426842e-10)
    assert_close1d(calc, expect, rtol=1e-13)

    assert_close(BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a, b)[0], -0.0003371378639247552, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a, b)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a, b)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Tsonopoulos_extended_fast(T, Tc, Pc, omega, a, b)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Tsonopoulos_extended_vec(T, [Tc], [Pc], [omega], [a], [b])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Tsonopoulos_extended_vec(T, [Tc], [Pc], [omega], [a], [b], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Tsonopoulos_extended_mat(T, [[Tc]], [[Pc]], [[omega]], [[a]], [[b]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Tsonopoulos_extended_mat(T, [[Tc]], [[Pc]], [[omega]], [[a]], [[b]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

def test_BVirial_Oconnell_Prausnitz():
    Tc = 508.1
    Pc = 4696e3
    T = 0.7*Tc
    omega = 0.308
    expect = [-0.0011699094955283378, 1.5320766743417143e-05, -3.131112087345748e-07, 8.09517379066269e-09]
    assert_close(BVirial_Oconnell_Prausnitz(T, Tc, Pc, omega)[0], -0.0011699094955283378, rtol=1e-13)
    assert_close(derivative(lambda T: BVirial_Oconnell_Prausnitz(T, Tc, Pc, omega)[0], T, dx=T*1e-6), expect[1])
    assert_close(derivative(lambda T: BVirial_Oconnell_Prausnitz(T, Tc, Pc, omega)[1], T, dx=T*1e-6), expect[2])
    assert_close(derivative(lambda T: BVirial_Oconnell_Prausnitz(T, Tc, Pc, omega)[2], T, dx=T*1e-6), expect[3])


    # Vector call with out memory savings
    vec_call = BVirial_Oconnell_Prausnitz_vec(T, [Tc], [Pc], [omega])
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, vec_call, rtol=1e-13)

    Bs_out = [0]
    dBs_out = [0]
    d2Bs_out = [0]
    d3Bs_out = [0]

    # vector call with memory savings
    BVirial_Oconnell_Prausnitz_vec(T, [Tc], [Pc], [omega], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    expect_vec = [[expect[0]], [expect[1]], [expect[2]], [expect[3]]]
    assert_close2d(expect_vec, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)

    # matrix call
    mat_call = BVirial_Oconnell_Prausnitz_mat(T, [[Tc]], [[Pc]], [[omega]])
    expect_mat = [[[expect[0]]], [[expect[1]]], [[expect[2]]], [[expect[3]]]]
    assert_close2d(expect_mat, mat_call, rtol=1e-13)

    # matrix call for memory savings
    Bs_out = [[0]]
    dBs_out = [[0]]
    d2Bs_out = [[0]]
    d3Bs_out = [[0]]
    BVirial_Oconnell_Prausnitz_mat(T, [[Tc]], [[Pc]], [[omega]], Bs_out, dBs_out, d2Bs_out, d3Bs_out)
    assert_close3d(expect_mat, [Bs_out, dBs_out, d2Bs_out, d3Bs_out], rtol=1e-13)


def test_BVirial_mixture():
    Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    zs = [.5, .3, .2]
    ans = BVirial_mixture(zs=zs, Bijs=Bijs)
    assert_close(ans, -3.19884e-05, rtol=1e-13)

def test_CVirial_mixture_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    zs = [.5, .3, .2]
    C = CVirial_mixture_Orentlicher_Prausnitz(zs, Cijs)
    assert_close(C, 2.079044009541466e-09, rtol=1e-13)

def test_dCVirial_mixture_dT_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    zs = [.5, .3, .2]
    dC = dCVirial_mixture_dT_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs)
    assert_close(dC, -7.275151799622592e-12, rtol=1e-14)

def test_d2CVirial_mixture_dT2_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    zs = [.5, .3, .2]
    ans = d2CVirial_mixture_dT2_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s)
    assert_close(ans, 6.723710778756013e-14, rtol=1e-14)

def test_d3CVirial_mixture_dT3_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    d3Cij_dT3s = [[-4.2300e-16, -7.9727e-16, -1.6962e-15], [-7.9727e-16, -1.3826e-15, -1.4525e-15], [-1.6962e-15, -1.4525e-15,  1.9786e-14]]
    zs = [.5, .3, .2]
    ans = d3CVirial_mixture_dT3_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s, d3Cij_dT3s)
    assert_close(ans, -3.735836855582578e-16, rtol=1e-14)

def test_dBVirial_mixture_dzs():
    Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    zs = [.5, .3, .2]
    calc = dBVirial_mixture_dzs(zs=zs, Bijs=Bijs)
    expect = [-3.4089e-05, -7.2301e-05, -0.00012621]
    assert_close1d(calc, expect, rtol=1e-13)

    numerical = jacobian(BVirial_mixture, zs, args=(Bijs,), perturbation=1e-7)
    assert_close1d(numerical, expect, rtol=1e-7)

    out = [0.0]*3
    calc = dBVirial_mixture_dzs(zs=zs, Bijs=Bijs, dB_dzs=out)
    assert calc is out

def test_d2BVirial_mixture_dzizjs():
    Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    zs = [.5, .3, .2]
    calc = d2BVirial_mixture_dzizjs(zs=zs, Bijs=Bijs)
    expect = [[-1.248e-05, -4.023e-05, -7.89e-05], [-4.023e-05, -8.782e-05, -0.0001292], [-7.89e-05, -0.0001292, -0.00024]]
    assert_close2d(calc, expect, rtol=1e-13)

    out = [[0.0]*3 for _ in range(3)]
    calc = d2BVirial_mixture_dzizjs(zs=zs, Bijs=Bijs, d2B_dzizjs=out)
    assert calc is out

    numerical = hessian(BVirial_mixture, zs, args=(Bijs,), perturbation=25e-5)
    assert_close1d(numerical, expect, rtol=1e-7)


def test_d3BVirial_mixture_dzizjzks():
    Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    zs = [.5, .3, .2]
    out = d3BVirial_mixture_dzizjzks(zs=zs, Bijs=Bijs)
    expect = [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
    assert_close3d(out, expect, atol=0)

    data = [[[0.0]*3 for _ in range(3)] for _ in range(3)]
    out = d3BVirial_mixture_dzizjzks(zs=zs, Bijs=Bijs, d3B_dzizjzks=data)
    assert out is data


def test_dCVirial_mixture_Orentlicher_Prausnitz_dzs():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    zs = [.5, .3, .2]
    dCs = dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs)
    expect = [5.444504708906765e-09, 6.549687763198116e-09, 7.749866726057895e-09]
    assert_close1d(dCs, expect, rtol=1e-13)

    def to_jac(zs): return CVirial_mixture_Orentlicher_Prausnitz(zs, Cijs)

    numerical = jacobian(to_jac, zs, perturbation=1e-7)
    assert_close1d(dCs, numerical, rtol=1e-7)


    # Large case to test the indexes better
    Cijs = [[1.465189176575471e-09, 1.860435325911775e-09, 2.1207286035797035e-09, 3.4042362688416274e-09, 5.466385407462434e-09, 1.3120738106399996e-09], [1.860435325911775e-09, 2.4658163172184164e-09, 3.2101245346622573e-09, 4.92501739124829e-09, 8.157718655133752e-09, 1.6899226468398887e-09], [2.1207286035797035e-09, 3.2101245346622573e-09, 4.927821768068707e-09, 7.308756958853639e-09, 1.1651843838258561e-08, 2.0186446325603855e-09], [3.4042362688416274e-09, 4.92501739124829e-09, 7.308756958853639e-09, 1.047676572661861e-08, 1.606047782005346e-08, 3.228332375603084e-09], [5.466385407462434e-09, 8.157718655133752e-09, 1.1651843838258561e-08, 1.606047782005346e-08, 1.940333988459971e-08, 5.3077543954058676e-09], [1.3120738106399996e-09, 1.6899226468398887e-09, 2.0186446325603855e-09, 3.228332375603084e-09, 5.3077543954058676e-09, 1.1706982435200877e-09]]
    zs = [.105, .095, .1, .11, .09, .5]
    dCs = dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs)
    expect = [6.723078338327948e-09, 8.31259485308936e-09, 1.0024423007706207e-08, 1.3182720227087849e-08, 1.7513764717771297e-08, 6.37515424248066e-09]
    assert_close1d(expect, dCs, rtol=1e-13)

    numerical = jacobian(to_jac, zs, perturbation=1e-7)
    assert_close1d(dCs, numerical, rtol=1e-7)

def test_d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    zs = [.5, .3, .2]
    d2Cs = d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs(zs, Cijs)
    expect =[[9.682788665539636e-09, 1.144914672502962e-08, 1.306435533767415e-08], [1.144914672502962e-08, 1.3855767429479784e-08, 1.6090359675187448e-08], [1.306435533767415e-08, 1.6090359675187448e-08, 2.070223940361239e-08]]
    assert_close2d(d2Cs, expect, rtol=1e-13)

    def to_jac(zs): return dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs)

    numerical = jacobian(to_jac, zs, perturbation=1e-7, scalar=False)
    assert_close2d(d2Cs, numerical, rtol=1e-7)


    # Large case to test the indexes better
    Cijs = [[1.465189176575471e-09, 1.860435325911775e-09, 2.1207286035797035e-09, 3.4042362688416274e-09, 5.466385407462434e-09, 1.3120738106399996e-09], [1.860435325911775e-09, 2.4658163172184164e-09, 3.2101245346622573e-09, 4.92501739124829e-09, 8.157718655133752e-09, 1.6899226468398887e-09], [2.1207286035797035e-09, 3.2101245346622573e-09, 4.927821768068707e-09, 7.308756958853639e-09, 1.1651843838258561e-08, 2.0186446325603855e-09], [3.4042362688416274e-09, 4.92501739124829e-09, 7.308756958853639e-09, 1.047676572661861e-08, 1.606047782005346e-08, 3.228332375603084e-09], [5.466385407462434e-09, 8.157718655133752e-09, 1.1651843838258561e-08, 1.606047782005346e-08, 1.940333988459971e-08, 5.3077543954058676e-09], [1.3120738106399996e-09, 1.6899226468398887e-09, 2.0186446325603855e-09, 3.228332375603084e-09, 5.3077543954058676e-09, 1.1706982435200877e-09]]
    zs = [.105, .095, .1, .11, .09, .5]
    d2Cs = d2CVirial_mixture_Orentlicher_Prausnitz_dzizjs(zs, Cijs)
    expect = [[1.071059420883899e-08, 1.2889939578666922e-08, 1.4774121572395523e-08, 1.9844462997775223e-08, 2.679252138715349e-08, 1.0050740025831605e-08], [1.2889939578666922e-08, 1.5744639953466747e-08, 1.887821938847448e-08, 2.496687544775576e-08, 3.403908283614765e-08, 1.2156619122970956e-08], [1.4774121572395523e-08, 1.887821938847448e-08, 2.394400241774843e-08, 3.129230073624375e-08, 4.2095296823327694e-08, 1.4158004743089302e-08], [1.9844462997775223e-08, 2.496687544775576e-08, 3.129230073624375e-08, 4.044962155652963e-08, 5.374127429084261e-08, 1.898903108170806e-08], [2.679252138715349e-08, 3.403908283614765e-08, 4.2095296823327694e-08, 5.374127429084261e-08, 6.604354929432281e-08, 2.5831225059285898e-08], [1.0050740025831605e-08, 1.2156619122970956e-08, 1.4158004743089302e-08, 1.898903108170806e-08, 2.5831225059285898e-08, 9.421395633868426e-09]]
    assert_close2d(expect, d2Cs, rtol=1e-13)

    numerical = jacobian(to_jac, zs, perturbation=1e-7, scalar=False)
    assert_close2d(d2Cs, numerical, rtol=1e-7)

def test_d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    zs = [.5, .3, .2]
    d3Cs = d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks(zs, Cijs)
    expect =[[[8.76000000000001e-09, 1.0187346981802932e-08, 1.1232922854993754e-08], [1.0187346981802932e-08, 1.2122397359312199e-08, 1.3593770131672463e-08], [1.1232922854993754e-08, 1.3593770131672463e-08, 1.6848814353377677e-08]], [[1.0187346981802932e-08, 1.2122397359312199e-08, 1.3593770131672463e-08], [1.2122397359312199e-08, 1.4760000000000021e-08, 1.6832843749118396e-08], [1.3593770131672463e-08, 1.6832843749118396e-08, 2.1218107423078488e-08]], [[1.1232922854993754e-08, 1.3593770131672463e-08, 1.6848814353377677e-08], [1.3593770131672463e-08, 1.6832843749118396e-08, 2.1218107423078488e-08], [1.6848814353377677e-08, 2.1218107423078488e-08, 2.956200000000003e-08]]]
    assert_close3d(d3Cs, expect, rtol=1e-13)

    def to_jac(zs): return dCVirial_mixture_Orentlicher_Prausnitz_dzs(zs, Cijs)

    numerical = hessian(to_jac, zs, perturbation=3e-5, scalar=False)
    assert_close3d(d3Cs, numerical, rtol=3e-5)


    # Large case to test the indexes better
    Cijs = [[1.465189176575471e-09, 1.860435325911775e-09, 2.1207286035797035e-09, 3.4042362688416274e-09, 5.466385407462434e-09, 1.3120738106399996e-09], [1.860435325911775e-09, 2.4658163172184164e-09, 3.2101245346622573e-09, 4.92501739124829e-09, 8.157718655133752e-09, 1.6899226468398887e-09], [2.1207286035797035e-09, 3.2101245346622573e-09, 4.927821768068707e-09, 7.308756958853639e-09, 1.1651843838258561e-08, 2.0186446325603855e-09], [3.4042362688416274e-09, 4.92501739124829e-09, 7.308756958853639e-09, 1.047676572661861e-08, 1.606047782005346e-08, 3.228332375603084e-09], [5.466385407462434e-09, 8.157718655133752e-09, 1.1651843838258561e-08, 1.606047782005346e-08, 1.940333988459971e-08, 5.3077543954058676e-09], [1.3120738106399996e-09, 1.6899226468398887e-09, 2.0186446325603855e-09, 3.228332375603084e-09, 5.3077543954058676e-09, 1.1706982435200877e-09]]
    zs = [.105, .095, .1, .11, .09, .5]
    d3Cs = d3CVirial_mixture_Orentlicher_Prausnitz_dzizjzks(zs, Cijs)
    expect = [[[8.791135059452836e-09, 1.0308422318043785e-08, 1.1248791844050504e-08, 1.5421582703357614e-08, 2.1147026843830464e-08, 8.167478419326304e-09], [1.0308422318043785e-08, 1.2261620654600401e-08, 1.3985943446213539e-08, 1.8887099304259298e-08, 2.6168516315435204e-08, 9.622719073212483e-09], [1.1248791844050504e-08, 1.3985943446213539e-08, 1.685361143642469e-08, 2.2504395584853537e-08, 3.0785384558550047e-08, 1.0665609066268141e-08], [1.5421582703357614e-08, 1.8887099304259298e-08, 2.2504395584853537e-08, 2.9710283153602206e-08, 4.011536227739754e-08, 1.4603938139341327e-08], [2.1147026843830464e-08, 2.6168516315435204e-08, 3.0785384558550047e-08, 4.011536227739754e-08, 5.003150644938515e-08, 2.0184021263543112e-08], [8.167478419326304e-09, 9.622719073212483e-09, 1.0665609066268141e-08, 1.4603938139341327e-08, 2.0184021263543112e-08, 7.578880928347833e-09]], [[1.0308422318043785e-08, 1.2261620654600401e-08, 1.3985943446213539e-08, 1.8887099304259298e-08, 2.6168516315435204e-08, 9.622719073212483e-09], [1.2261620654600401e-08, 1.4794897903310516e-08, 1.7639492188523473e-08, 2.346442970935061e-08, 3.2848633417587723e-08, 1.150048197891079e-08], [1.3985943446213539e-08, 1.7639492188523473e-08, 2.2218641495188888e-08, 2.9224178168642914e-08, 4.0393641206428843e-08, 1.3323984224128247e-08], [1.8887099304259298e-08, 2.346442970935061e-08, 2.9224178168642914e-08, 3.800425063355443e-08, 5.18477086792256e-08, 1.7970860061469278e-08], [2.6168516315435204e-08, 3.2848633417587723e-08, 4.0393641206428843e-08, 5.18477086792256e-08, 6.533650453500665e-08, 2.509574192969567e-08], [9.622719073212483e-09, 1.150048197891079e-08, 1.3323984224128247e-08, 1.7970860061469278e-08, 2.509574192969567e-08, 8.971756058880126e-09]], [[1.1248791844050504e-08, 1.3985943446213539e-08, 1.685361143642469e-08, 2.2504395584853537e-08, 3.0785384558550047e-08, 1.0665609066268141e-08], [1.3985943446213539e-08, 1.7639492188523473e-08, 2.2218641495188888e-08, 2.9224178168642914e-08, 4.0393641206428843e-08, 1.3323984224128247e-08], [1.685361143642469e-08, 2.2218641495188888e-08, 2.9566930608412272e-08, 3.845315051033282e-08, 5.247663974134904e-08, 1.6308330162363284e-08], [2.2504395584853537e-08, 2.9224178168642914e-08, 3.845315051033282e-08, 4.944503342771887e-08, 6.660141603770983e-08, 2.1749292204673606e-08], [3.0785384558550047e-08, 4.0393641206428843e-08, 5.247663974134904e-08, 6.660141603770983e-08, 8.286539109452554e-08, 2.998746118655785e-08], [1.0665609066268141e-08, 1.3323984224128247e-08, 1.6308330162363284e-08, 2.1749292204673606e-08, 2.998746118655785e-08, 1.0100421248596669e-08]], [[1.5421582703357614e-08, 1.8887099304259298e-08, 2.2504395584853537e-08, 2.9710283153602206e-08, 4.011536227739754e-08, 1.4603938139341327e-08], [1.8887099304259298e-08, 2.346442970935061e-08, 2.9224178168642914e-08, 3.800425063355443e-08, 5.18477086792256e-08, 1.7970860061469278e-08], [2.2504395584853537e-08, 2.9224178168642914e-08, 3.845315051033282e-08, 4.944503342771887e-08, 6.660141603770983e-08, 2.1749292204673606e-08], [2.9710283153602206e-08, 3.800425063355443e-08, 4.944503342771887e-08, 6.286059435971174e-08, 8.357299557997483e-08, 2.867779938135162e-08], [4.011536227739754e-08, 5.18477086792256e-08, 6.660141603770983e-08, 8.357299557997483e-08, 1.0263190467585411e-07, 3.902717277758873e-08], [1.4603938139341327e-08, 1.7970860061469278e-08, 2.1749292204673606e-08, 2.867779938135162e-08, 3.902717277758873e-08, 1.3812906337677218e-08]], [[2.1147026843830464e-08, 2.6168516315435204e-08, 3.0785384558550047e-08, 4.011536227739754e-08, 5.003150644938515e-08, 2.0184021263543112e-08], [2.6168516315435204e-08, 3.2848633417587723e-08, 4.0393641206428843e-08, 5.18477086792256e-08, 6.533650453500665e-08, 2.509574192969567e-08], [3.0785384558550047e-08, 4.0393641206428843e-08, 5.247663974134904e-08, 6.660141603770983e-08, 8.286539109452554e-08, 2.998746118655785e-08], [4.011536227739754e-08, 5.18477086792256e-08, 6.660141603770983e-08, 8.357299557997483e-08, 1.0263190467585411e-07, 3.902717277758873e-08], [5.003150644938515e-08, 6.533650453500665e-08, 8.286539109452554e-08, 1.0263190467585411e-07, 1.1642003930759836e-07, 4.905884204966277e-08], [2.0184021263543112e-08, 2.509574192969567e-08, 2.998746118655785e-08, 3.902717277758873e-08, 4.905884204966277e-08, 1.924155286926517e-08]], [[8.167478419326304e-09, 9.622719073212483e-09, 1.0665609066268141e-08, 1.4603938139341327e-08, 2.0184021263543112e-08, 7.578880928347833e-09], [9.622719073212483e-09, 1.150048197891079e-08, 1.3323984224128247e-08, 1.7970860061469278e-08, 2.509574192969567e-08, 8.971756058880126e-09], [1.0665609066268141e-08, 1.3323984224128247e-08, 1.6308330162363284e-08, 2.1749292204673606e-08, 2.998746118655785e-08, 1.0100421248596669e-08], [1.4603938139341327e-08, 1.7970860061469278e-08, 2.1749292204673606e-08, 2.867779938135162e-08, 3.902717277758873e-08, 1.3812906337677218e-08], [2.0184021263543112e-08, 2.509574192969567e-08, 2.998746118655785e-08, 3.902717277758873e-08, 4.905884204966277e-08, 1.924155286926517e-08], [7.578880928347833e-09, 8.971756058880126e-09, 1.0100421248596669e-08, 1.3812906337677218e-08, 1.924155286926517e-08, 7.024189461120531e-09]]]
    assert_close3d(expect, d3Cs, rtol=1e-13)

    numerical = hessian(to_jac, zs, perturbation=8e-5, scalar=False)
    assert_close3d(d3Cs, numerical, rtol=3e-5)

def test_dV_dzs_virial():
    calc = dV_dzs_virial(B=-5.130920247359858e-05, C=2.6627784284381213e-09, V=0.024892080086430797, dB_dzs=[-4.457911131778849e-05, -9.174964457681726e-05, -0.0001594258679841028], dC_dzs=[6.270599057032657e-09, 7.766612052069565e-09, 9.503031492910165e-09])
    expect = [-4.4510120473455416e-05, -9.181495962913208e-05, -0.00015970040988493522]
    assert_close1d(calc, expect, rtol=1e-14)


def test_d2V_dzizjs_virial():
    d2C_dzizjs = [[1.0287075724127612e-08, 1.2388277824773021e-08, 1.4298813522844275e-08], [1.2388277824773021e-08, 1.514162073913238e-08, 1.8282527232061114e-08], [1.4298813522844275e-08, 1.8282527232061114e-08, 2.3350122217403063e-08]]
    d2B_dzizjs = [[-1.0639357784985337e-05, -3.966321845899801e-05, -7.53987684376414e-05], [-3.966321845899801e-05, -8.286257232134107e-05, -0.00014128571574782375], [-7.53987684376414e-05, -0.00014128571574782375, -0.00024567752140887547]]
    dB_dzs = [-4.457911131778849e-05, -9.174964457681726e-05, -0.0001594258679841028]
    dC_dzs = [6.270599057032657e-09, 7.766612052069565e-09, 9.503031492910165e-09]
    dV_dzs = [-4.4510120473455416e-05, -9.181495962913208e-05, -0.00015970040988493522]
    calc = d2V_dzizjs_virial(B=-5.130920247359858e-05, C=2.6627784284381213e-09, V=0.024892080086430797, dB_dzs=dB_dzs, dC_dzs=dC_dzs, dV_dzs=dV_dzs, d2B_dzizjs=d2B_dzizjs, d2C_dzizjs=d2C_dzizjs)
    expect = [[-1.0426891738921789e-05, -3.965469458896068e-05, -7.570310078917325e-05],
     [-3.965469458896068e-05, -8.327011676733816e-05, -0.00014230835847955974],
     [-7.570310078917325e-05, -0.00014230835847955974, -0.000247797888252081]]

    assert_close2d(calc, expect, rtol=1e-14)


def test_d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs():
    Cijs = [[1.46e-09, 1.831e-09, 2.12e-09], [1.831e-09, 2.46e-09, 2.996e-09], [2.12e-09, 2.996e-09, 4.927e-09]]

    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    zs = [.5, .3, .2]
    calc = d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs(zs, Cijs, dCij_dTs)
    expect = [-1.574099410300064e-11, -2.2726730950177544e-11, -3.5684695311570996e-11]
    assert_close(calc, expect, rtol=1e-13)

    # larger test case

    Cijs = [[1.465189176575471e-09, 1.860435325911775e-09, 2.1207286035797035e-09, 3.4042362688416274e-09, 5.466385407462434e-09, 1.3120738106399996e-09], [1.860435325911775e-09, 2.4658163172184164e-09, 3.2101245346622573e-09, 4.92501739124829e-09, 8.157718655133752e-09, 1.6899226468398887e-09], [2.1207286035797035e-09, 3.2101245346622573e-09, 4.927821768068707e-09, 7.308756958853639e-09, 1.1651843838258561e-08, 2.0186446325603855e-09], [3.4042362688416274e-09, 4.92501739124829e-09, 7.308756958853639e-09, 1.047676572661861e-08, 1.606047782005346e-08, 3.228332375603084e-09], [5.466385407462434e-09, 8.157718655133752e-09, 1.1651843838258561e-08, 1.606047782005346e-08, 1.940333988459971e-08, 5.3077543954058676e-09], [1.3120738106399996e-09, 1.6899226468398887e-09, 2.0186446325603855e-09, 3.228332375603084e-09, 5.3077543954058676e-09, 1.1706982435200877e-09]]
    # approximate the values
    dCij_dTs = [[-i*1e-3 for i in r] for r in Cijs]
    zs = [.105, .095, .1, .11, .09, .5]
    calc = d2CVirial_mixture_Orentlicher_Prausnitz_dTdzs(zs, Cijs, dCij_dTs)
    expect = [-6.723078338327947e-12, -8.312594853089363e-12, -1.0024423007706204e-11, -1.3182720227087848e-11, -1.75137647177713e-11, -6.375154242480662e-12]
    assert_close(calc, expect, rtol=1e-13)
