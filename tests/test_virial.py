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

import pytest
import numpy as np
from fluids.numerics import linspace, assert_close, assert_close1d, assert_close2d, assert_close3d, derivative
from chemicals.virial import *
from chemicals import rho_to_Vm
from fluids.constants import R as _R
from scipy.integrate import quad

def test_B_To_Z():
    Z_calc = B_to_Z(-0.0015, 300.0, 1E5)
    assert_close(Z_calc, 0.9398638020957176)


def test_B_from_Z():
    B_calc = B_from_Z(0.94, 300.0, 1E5)
    assert_close(B_calc, -0.0014966032712675846)


def test_Z_from_virial_density_form():
    Z_calc = Z_from_virial_density_form(300.0, 122057.233762653, 1E-4, 1E-5, 1E-6, 1E-7)
    assert_close(Z_calc, 1.2843494052609186)

    Z_calc = Z_from_virial_density_form(300, 102031.881198762, 1e-4, 1e-5, 1e-6)
    assert_close(Z_calc, 1.0736323841544937)

    Z_calc = Z_from_virial_density_form(300, 96775.8831504971, 1e-4, 1e-5)
    assert_close(Z_calc, 1.018326089216066)

    Z_calc = Z_from_virial_density_form(300, 95396.3561037084, 1e-4)
    assert_close(Z_calc,  1.003809998713499)

    assert_close(1, Z_from_virial_density_form(300, 95396.3561037084))

    '''B-only solution, derived as follows:

    >>> B, C, D, E = symbols('B, C, D, E')
    >>> P, V, Z, R, T = symbols('P, V, Z, R, T', positive=True, real=True, nonzero=True)
    >>> rho = 1/V
    >>> to_slv = Eq(P*V/R/T, 1 + B*rho)
    >>> slns = solve(to_slv, V)
    >>> simplify(slns[1]*P/R/T)
    1/2 + sqrt(4*B*P + R*T)/(2*sqrt(R)*sqrt(T))

    To check this, simply disable the if statement and allow the numerical
    algorithm to run.
    '''


def test_Z_from_virial_pressure_form():
    Z_calc = Z_from_virial_pressure_form(102919.99946855308, 4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19)
    assert_close(Z_calc, 1.00283753944)

    Z_calc = Z_from_virial_pressure_form(102847.17619188508, 4.032286555169439e-09, 1.6197059494442215e-13)
    assert_close(Z_calc, 1.00212796)

    Z_calc = Z_from_virial_pressure_form(102671.27455742132, 4.032286555169439e-09)
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
    from sympy import symbols, Rational, diff, lambdify, integrate
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
def test_BVirial_Abbott_calculus():
    from sympy import symbols, Rational, diff, lambdify, integrate

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
    from sympy import symbols, Rational, diff, lambdify, integrate

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


@pytest.mark.slow
@pytest.mark.sympy
def test_BVirial_Tsonopoulos_extended_calculus():
    from sympy import symbols, Rational, diff, lambdify, integrate


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
    calc = BVirial_Meng(388.26, 647.1, 22050000.0, 5.543076e-05, 0.344)
    expect = (-0.00032436028497558636, 2.470038900338557e-06, -3.132003987118147e-08, 5.776332655071256e-10)
    assert_close1d(calc, expect, rtol=1e-13)
    
def test_BVirial_mixture():
    Bijs = [[-6.24e-06, -2.013e-05, -3.9e-05], [-2.01e-05, -4.391e-05, -6.46e-05], [-3.99e-05, -6.46e-05, -0.00012]]
    zs = [.5, .3, .2]
    ans = BVirial_mixture(zs=zs, Bijs=Bijs)
    assert_close(ans, -0.000116432, rtol=1e-13)
    
def test_CVirial_mixture_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.1207e-09], [1.83e-09, 2.46e-09, 2.996e-09], [2.120e-09, 2.996e-09, 4.927e-09]]
    zs = [.5, .3, .2]
    C = CVirial_mixture_Orentlicher_Prausnitz(zs, Cijs)
    assert_close(C, 2.0787313269445096e-09, rtol=1e-13)
    
def test_dCVirial_mixture_dT_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.1207e-09], [1.83e-09, 2.46e-09, 2.996e-09], [2.120e-09, 2.996e-09, 4.927e-09]]
    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    zs = [.5, .3, .2]
    dC = dCVirial_mixture_dT_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs)
    assert_close(dC, -7.276497863811498e-12, rtol=1e-14)
    
def test_d2CVirial_mixture_dT2_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.1207e-09], [1.83e-09, 2.46e-09, 2.996e-09], [2.120e-09, 2.996e-09, 4.927e-09]]
    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    zs = [.5, .3, .2]
    ans = d2CVirial_mixture_dT2_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s)
    assert_close(ans, 6.723627522976708e-14, rtol=1e-14)
    
def test_d3CVirial_mixture_dT3_Orentlicher_Prausnitz():
    Cijs = [[1.46e-09, 1.831e-09, 2.1207e-09], [1.83e-09, 2.46e-09, 2.996e-09], [2.120e-09, 2.996e-09, 4.927e-09]]
    dCij_dTs = [[-2.212e-12, -4.137e-12, -1.079e-11], [-4.137e-12, -7.669e-12, -1.809e-11], [-1.079e-11, -1.809e-11, -2.010e-11]]
    d2Cij_dT2s = [[ 2.6469e-14,  5.0512e-14,  1.1509e-13], [ 5.0512e-14,  9.3272e-14,  1.7836e-13], [ 1.1509e-13,  1.7836e-13, -1.4906e-13]]
    d3Cij_dT3s = [[-4.2300e-16, -7.9727e-16, -1.6962e-15], [-7.9727e-16, -1.3826e-15, -1.4525e-15], [-1.6962e-15, -1.4525e-15,  1.9786e-14]]
    zs = [.5, .3, .2]
    ans = d3CVirial_mixture_dT3_Orentlicher_Prausnitz(zs, Cijs, dCij_dTs, d2Cij_dT2s, d3Cij_dT3s)
    assert_close(ans, -3.7356470379612563e-16, rtol=1e-14)
