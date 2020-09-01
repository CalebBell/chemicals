# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from __future__ import division
from chemicals import *
import chemicals.vectorized
from math import *
from random import random
from fluids.constants import *
from fluids.numerics import assert_close, assert_close1d, assert_close2d
from numpy.testing import assert_allclose
import pytest
try:
    import numba
    import chemicals.numba
    import chemicals.numba_vectorized
except:
    numba = None
import numpy as np


def swap_funcs_and_test(names, substitutions, test):
    '''
    names : list[str]
        object names to switch out
    substitutions : list[obj]
        Objects to put in
    test : function
        Unit test to run in the file
    '''
    originals = {}
    glob = test.__globals__
    for name, sub in zip(names, substitutions):
        originals[name] = glob[name]
        glob[name] = sub
    try:
        test()
    except Exception as e:
        glob.update(originals)
        raise e
    glob.update(originals)

def mark_as_numba(func):
    func = pytest.mark.numba(func)
    func = pytest.mark.slow(func)
    func = pytest.mark.skipif(numba is None, reason="Numba is missing")(func)
    return func
    

@mark_as_numba
def test_return_1d_array():
    
    # Functions which initialize an array, and then need to return the correct value 
    N = 30
    zs = zs_orig = normalize([random() for i in range(N)])
    MWs = [random()*200 for i in range(N)]
    zs2 = np.array(zs)
    MWs2 = np.array(MWs)
    
    # Took the slightest performance hit to CPython only, 186 us original, 190 us revised
    # at 1000 elements; no performance difference < 50 compounds
    ws = zs_to_ws(zs, MWs)
    ws_np = chemicals.numba.zs_to_ws(zs2, MWs2)
    assert type(ws_np) is np.ndarray
    assert_close1d(ws, ws_np)
    
    zs = ws_to_zs(ws, MWs)
    zs_np = chemicals.numba.ws_to_zs(ws_np, MWs2)
    assert type(zs_np) is np.ndarray
    assert_close1d(zs, zs_np)
    
    # Treat MWs as Vfs; doesn't matter to math
    Vfs = zs_to_Vfs(zs, MWs)
    Vfs_np = chemicals.numba.zs_to_Vfs(zs2, MWs2)
    assert type(Vfs_np) is np.ndarray
    assert_close1d(Vfs, Vfs_np)
    
    zs = Vfs_to_zs(Vfs, MWs)
    zs_np = chemicals.numba.Vfs_to_zs(Vfs_np, MWs2)
    assert type(Vfs_np) is np.ndarray
    assert_close1d(zs, zs_np)
    
    # Functions which have a return list comprehension 
    vals = [-2651.3181821109024, -2085.574403592012, -2295.0860830203587]
    dxsn1 = chemicals.dxs_to_dxsn1(vals)
    dxsn1_np = chemicals.numba.dxs_to_dxsn1(np.array(vals))
    assert_close1d(dxsn1, dxsn1_np)
    assert type(dxsn1_np) is np.ndarray

    dxs, xs = [-0.0028, -0.00719, -0.00859], [0.7, 0.2, 0.1]
    dns = dxs_to_dns(dxs, xs)
    dns_np = chemicals.numba.dxs_to_dns(np.array(dxs), np.array(xs))
    assert type(dns_np) is np.ndarray
    assert_close1d(dns, dns_np)
    
    dns = [0.001459, -0.002939, -0.004334]
    dn_partials = dns_to_dn_partials(dns, -0.0016567)
    dn_partials_np = chemicals.numba.dns_to_dn_partials(np.array(dns), -0.0016567)
    assert type(dn_partials_np) is np.ndarray
    assert_close1d(dn_partials_np, dn_partials)

    dxs = [-0.0026404, -0.00719, -0.00859]
    xs = [0.7, 0.2, 0.1]
    F = -0.0016567
    dn_partials = dxs_to_dn_partials(dxs, xs, F)
    dn_partials_np = chemicals.numba.dxs_to_dn_partials(np.array(dxs), np.array(xs), F)
    assert_close1d(dn_partials, dn_partials_np)
    assert type(dn_partials_np) is np.ndarray


@mark_as_numba
def test_return_2d_array():
    d2xs = [[0.152, 0.08, 0.547], [0.08, 0.674, 0.729], [0.547, 0.729, 0.131]]
    xs = [0.7, 0.2, 0.1]
    dxdn_partials = d2xs_to_dxdn_partials(d2xs, xs)
    a, b = np.array(d2xs), np.array(xs)
    dxdn_partials_np = chemicals.numba.d2xs_to_dxdn_partials(a, b)
    
    assert type(dxdn_partials_np) is np.ndarray
    assert_close1d(dxdn_partials, dxdn_partials_np)



@mark_as_numba
def test_mixing_simple():
    a = np.array([1,2])
    b = np.array([.1, .2])
    tot = chemicals.numba.mixing_simple(a, b)
    assert_close(tot, 0.5, rtol=1e-14)
    
    
    a = np.array([.1, .9])
    b = np.array([.01, .02])
    val = chemicals.numba.mixing_logarithmic(a, b)
    assert_close(val, 0.01866065983073615, rtol=1e-13)
    
    
@mark_as_numba
def test_dippr_correlations():
    orders = (0, 1, -1, -1j)
    args = (20, 33.19, 66.653, 6765.9, -123.63, 478.27)
    for i in orders:
        assert_close(chemicals.numba.EQ114(*args, order=i), chemicals.numba.EQ114(*args, order=i), rtol=1e-13)

    args = (300, 276370., -2090.1, 8.125, -0.014116, 0.0000093701)
    for i in orders:
        assert_close(chemicals.numba.EQ100(*args, order=i), chemicals.numba.EQ100(*args, order=i), rtol=1e-13)
        
    # EQ102 -  numba-scipy does not support complex numbers so this does not work in numba

    args = (300., 647.096, 17.863, 58.606, -95.396, 213.89, -141.26)
    for i in orders:
        assert_close(chemicals.numba.EQ116(*args, order=i), chemicals.numba.EQ116(*args, order=i), rtol=1e-13)
    
    args = (20., 3.3258E4, 3.6199E4, 1.2057E3, 1.5373E7, 3.2122E3, -1.5318E7, 3.2122E3)
    for i in orders:
        assert_close(chemicals.numba.EQ127(*args, order=i), chemicals.numba.EQ127(*args, order=i), rtol=1e-13)

    args = (300., 33363., 26790., 2610.5, 8896., 1169)
    for i in orders:
        assert_close(chemicals.numba.EQ107(*args, order=i), chemicals.numba.EQ107(*args, order=i), rtol=1e-13)

        args = (300.0, 0.02222, -26.38, -16750000, -3.894E19, 3.133E21)
    for i in orders:
        assert_close(chemicals.numba.EQ104(*args, order=i), chemicals.numba.EQ104(*args, order=i), rtol=1e-13)

@mark_as_numba
def test_thermal_conductivity_misc():
    assert_close(chemicals.numba.Bahadori_liquid(273.15, 170),
                 Bahadori_liquid(273.15, 170))

    assert_close(chemicals.numba.Missenard(304., 6330E5, 591.8, 41E5, 0.129),
                 chemicals.Missenard(304., 6330E5, 591.8, 41E5, 0.129))

    assert_close(chemicals.numba.DIPPR9H(np.array([0.258, 0.742]), np.array([0.1692, 0.1528])),
                 DIPPR9H([0.258, 0.742], [0.1692, 0.1528]))
    
    assert_close(chemicals.numba.Filippov(np.array([0.258, 0.742]), np.array([0.1692, 0.1528])),
                 Filippov([0.258, 0.742], [0.1692, 0.1528]))
    
    assert_close(chemicals.numba.DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear'),
                 chemicals.DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear'))

    assert_close(chemicals.numba.Eli_Hanley(T=373.15, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267, omega=0.227, Cvm=135.9),
                 chemicals.Eli_Hanley(T=373.15, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267, omega=0.227, Cvm=135.9))

    assert_close(chemicals.numba.Eli_Hanley_dense(T=473., MW=42.081, Tc=364.9, Vc=1.81E-4, Zc=0.274, omega=0.144, Cvm=82.70, Vm=1.721E-4),
                 chemicals.Eli_Hanley_dense(T=473., MW=42.081, Tc=364.9, Vc=1.81E-4, Zc=0.274, omega=0.144, Cvm=82.70, Vm=1.721E-4))

    assert_close(chemicals.numba.Chung_dense(T=473., MW=42.081, Tc=364.9, Vc=184.6E-6, omega=0.142, Cvm=82.67, Vm=172.1E-6, mu=134E-7, dipole=0.4),
                 chemicals.Chung_dense(T=473., MW=42.081, Tc=364.9, Vc=184.6E-6, omega=0.142, Cvm=82.67, Vm=172.1E-6, mu=134E-7, dipole=0.4))

    # Does not work - atom input
#    chemicals.numba.Mersmann_Kind_thermal_conductivity_liquid(400, 170.33484, 658.0, 0.000754, {'C': 12, 'H': 26})

@mark_as_numba
def test_viscosity_misc():
    # Has a min, if statement
    args = (300., 500E5, 572.2, 34.7E5, 0.236, 0, 0.00068)
    ans = chemicals.numba.Lucas(*args)
    ans_base = chemicals.viscosity.Lucas(*args)
    assert_close(ans, ans_base, rtol=1e-14)
    
    # There is a dict lokup but it is not always needed
    new = Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)
    fast = chemicals.numba.Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)
    assert_close(new, fast, rtol=1e-12)
    
    # Test the dict lookup has been turned into a couple if statements - not suitable for large
    # tables but for three elements it is just as fast as a dict lookup
    kwargs = dict(T=6, Tc=5.1889, Pc=226968.0, Zc=0.3014, MW=4.002602, CASRN='7440-59-7')
    assert_close(chemicals.numba.Lucas_gas(**kwargs), Lucas_gas(**kwargs), rtol=1e-14)
    
    # A couple of points with Herning-Sipperer; works fine
    zs = np.array([0.5, 0.25, 0.25]*10)
    mus = np.array([1.78e-05, 1.12e-05, 9.35e-06]*10)
    MWs = np.array([28.0134, 16.043, 30.07]*10)
    
    fast = chemicals.numba.Herning_Zipperer(zs, mus, MWs)
    base = chemicals.Herning_Zipperer(zs.tolist(), mus.tolist(), MWs.tolist())
    assert_close(fast, base, rtol=1e-14)
    
    
    # Function calling other functions
    n = 1
    zs = np.array([.4, .3, .3]*n)
    MWs = np.array([16.04246, 30.06904, 44.09562]*n)
    Tcs = np.array([190.564, 305.32, 369.83]*n)
    Pcs = np.array([4599000.0, 4872000.0, 4248000.0]*n)
    Vcs = np.array([9.86e-05, 0.0001455, 0.0002]*n)
    mu = chemicals.numba.Lorentz_Bray_Clarke(T=300.0, P=1e6, Vm=0.0023025, zs=zs, MWs=MWs, Tcs=Tcs, Pcs=Pcs, Vcs=Vcs)
    assert_close(mu, 9.925488160761484e-06, rtol=1e-14)
    
    # Viscosity index - works beautifully
    assert_close(chemicals.numba.viscosity_index(73.3E-6, 8.86E-6, rounding=False),
                 chemicals.viscosity_index(73.3E-6, 8.86E-6, rounding=False), rtol=1e-14)
    
    assert_close(chemicals.numba.viscosity_index(73.3E-6, 8.86E-6, rounding=True),
                 chemicals.viscosity_index(73.3E-6, 8.86E-6, rounding=True), rtol=1e-14)
    
@mark_as_numba
def test_interface_misc():
    
    # Tested quite a bit with numba/PyPy
    # At first numba had 3x the speed, but then I made the optimizations by hand
    # I knew were possible. Their speed is about equal after, with a slight edge up
    # by numba with large arrays
    n = 1
    xs = np.array([0.1606, 0.8394]*n)
    xs /= sum(xs)
    sigmas = np.array([0.01547, 0.02877]*n)
    rhoms = np.array([8610., 15530.]*n)
    xs2, sigmas2, rhoms2 = xs.tolist(), sigmas.tolist(), rhoms.tolist()
    assert_close(chemicals.numba.Winterfeld_Scriven_Davis(xs, sigmas, rhoms), 
                 Winterfeld_Scriven_Davis(xs2, sigmas2, rhoms2))
    
    
    n = 1
    xs = np.array([0.1606, 0.8394]*n)
    sigmas_Tb = np.array([0.01424, 0.02530]*n)
    Tbs = np.array([309.21, 312.95]*n)
    Tcs = np.array([469.7, 508.0]*n)
    assert_close(chemicals.Diguilio_Teja(T=298.15, xs=xs,sigmas_Tb=sigmas_Tb, Tbs=Tbs, Tcs=Tcs),
    chemicals.numba.Diguilio_Teja(T=298.15, xs=xs,sigmas_Tb=sigmas_Tb, Tbs=Tbs, Tcs=Tcs), rtol=1e-12)
    
    # Exception is correctly raised with numba
    with pytest.raises(ValueError):
        chemicals.numba.Diguilio_Teja(T=1000, xs=xs,sigmas_Tb=sigmas_Tb, Tbs=Tbs, Tcs=Tcs)
        
        
@mark_as_numba
def test_virial():
    Z = chemicals.numba.Z_from_virial_pressure_form(102919.99946855308, 4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19)
    assert_close(Z, 1.00283753944, rtol=1e-13)
    
#    # Takes 8 seconds to compile. Fun!
#    assert_close(chemicals.numba.BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469),
#                 chemicals.BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469),
#                rtol=1e-13)

@mark_as_numba
def test_phase_change():
    # Function had some duplicated powers; numba was optimizing them on me anyway
    # Had list-in-list constants being indexed. I thought that would take a lot of time
    # but instead removing it only saved 25%, and ~8% in CPython, and zilch in PyPy.
    # PyPy takes 19% of the time numba does here, numba has a high overhead.
    assert_close(chemicals.numba.MK(553.15, 751.35, 0.302),
             chemicals.MK(553.15, 751.35, 0.302), rtol=1e-12)


@mark_as_numba
def test_vapor_pressure():
    # PyPy 75 ns, CPython 2470 ns, numba 214 ns
    assert_close(chemicals.numba.dPsat_IAPWS_dT(300.), 
                 chemicals.dPsat_IAPWS_dT(300.), rtol=1e-14)
    
    
    Psats_vec_expect = [34478.367349639906, 33596697.716487624, 109799836.81382856, 179376011.49286702, 234627689.09298804]
    Ts = np.linspace(100, 1000, 5)
    Psats_calc = chemicals.numba_vectorized.Antoine(Ts, 8.7687, 395.744, -6.469, 10)
    assert_close(Psats_calc, Psats_vec_expect, rtol=1e-11)
    

@mark_as_numba
def test_temperature():
    # Note also the last four decimals are different!
    # 494 us numba, 388 us PyPy, 1740 us CPython
    assert_close(chemicals.numba.ITS90_68_difference(1000.),
                 chemicals.ITS90_68_difference(1000.0), rtol=1e-12)
    
    # Probably never going to work
#    chemicals.numba.T_converter(500, 'ITS-68', 'ITS-48')

@mark_as_numba
def test_critical():
    assert_close(chemicals.numba.Li(np.array([0.6449, 0.2359, 0.1192]), np.array([425.12, 469.7, 507.6]),np.array([0.000255, 0.000313, 0.000371])),
             Li([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [0.000255, 0.000313, 0.000371]), rtol=1e-13)

    assert_close(chemicals.numba.Chueh_Prausnitz_Tc(np.array([0.6449, 0.2359, 0.1192]), np.array([425.12, 469.7, 507.6]),
        np.array([0.000255, 0.000313, 0.000371]), np.array([[0, 1.92681, 6.80358],
        [1.92681, 0, 1.89312], [ 6.80358, 1.89312, 0]])),
    Chueh_Prausnitz_Tc([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6],
        [0.000255, 0.000313, 0.000371], [[0, 1.92681, 6.80358],
         [1.92681, 0, 1.89312], [ 6.80358, 1.89312, 0]]), rtol=1e-13)

    zs = np.array([0.6449, 0.2359, 0.1192])
    Tcs = np.array([425.12, 469.7, 507.6])
    Aijs = np.array([[0, 1.2503, 1.516], [0.799807, 0, 1.23843], [0.659633, 0.807474, 0]])
    
    
    assert_close(chemicals.numba.Grieves_Thodos(zs, Tcs, Aijs),
                 Grieves_Thodos(zs, Tcs, Aijs), rtol=1e-12)
    
    Aijs = np.array([[0, 1.174450, 1.274390], [0.835914, 0, 1.21038], [0.746878, 0.80677, 0]])
    assert_close(chemicals.numba.modified_Wilson_Tc(zs, Tcs, Aijs),
                 modified_Wilson_Tc(zs, Tcs, Aijs), rtol=1e-12)
    
    
    assert_close(chemicals.numba.Chueh_Prausnitz_Vc(np.array([0.4271, 0.5729]), np.array([0.000273, 0.000256]), np.array([[0, 5.61847], [5.61847, 0]])),
                 Chueh_Prausnitz_Vc([0.4271, 0.5729], [0.000273, 0.000256], [[0, 5.61847], [5.61847, 0]]), rtol=1e-13)

    assert_close(chemicals.numba.modified_Wilson_Vc(np.array([0.4271, 0.5729]), np.array([0.000273, 0.000256]), np.array([[0, 0.6671250], [1.3939900, 0]])),
                 modified_Wilson_Vc([0.4271, 0.5729], [0.000273, 0.000256], [[0, 0.6671250], [1.3939900, 0]]), rtol=1e-13)

    # Not working yet: Ihmels, Meissner, Grigoras, critical_surface_methods
    # Maybe a future numba update will make this work. 

@mark_as_numba
def test_volume():
    assert_close(chemicals.numba.Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245),
                chemicals.Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245))
    assert_close(chemicals.numba.COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532),
                 chemicals.COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532))

    assert_close(chemicals.numba.Bhirud_normal(280.0, 469.7, 33.7E5, 0.252),
                 Bhirud_normal(280.0, 469.7, 33.7E5, 0.252))

    assert_close(chemicals.numba.SNM0(121, 150.8, 7.49e-05, -0.004),
                 SNM0(121, 150.8, 7.49e-05, -0.004))
    
    assert_close(chemicals.numba.SNM0(121, 150.8, 7.49e-05, -0.004, -0.03259620),
                 SNM0(121, 150.8, 7.49e-05, -0.004, -0.03259620))

    
    kwargs = dict(T=405.45, Tb=239.82, Tc=405.45, Pc=111.7*101325, M=17.03, dipole=None)
    assert_close(chemicals.numba.Campbell_Thodos(**kwargs), 
                 Campbell_Thodos(**kwargs))

    # Test a slow one
    # 81.2 us orig, then 67.6 after optimizations in CPython
    # numba: 2.25 µs, PYPY: 1.31; numba with numpy: 4 us
        
    N = 100
    xs = [0.4576, 0.5424]*N
    MWs = [32.04, 18.01]*N
    Tcs = [512.58, 647.29]*N
    Pcs = [8.096E6, 2.209E7]*N
    Zrs = [0.2332, 0.2374]*N
    
    xs2 = np.array(xs)
    MWs2 = np.array(MWs)
    Tcs2 = np.array(Tcs)
    Pcs2 = np.array(Pcs)
    Zrs2 = np.array(Zrs)
    
    orig = Rackett_mixture(T=298., xs=xs, MWs=MWs, Tcs=Tcs, Pcs=Pcs, Zrs=Zrs)
    
    new = chemicals.numba.Rackett_mixture(T=298., xs=xs2, MWs=MWs2, Tcs=Tcs2, Pcs=Pcs2, Zrs=Zrs2)
    assert_close(orig, new)
    
    
    
    # Test COSTALD_mixture - even slower
    # timing after optimization at 200 elements - 1.49 m CPython, 27.1 µs numba, 63.5 µs PyPy3, 71.4 us PyPy2
    T = 300.0
    N = 15
    xs = normalize([0.4576, 0.5424]*N)
    Tcs = [512.58, 647.29]*N
    Vcs =  [0.000117, 5.6e-05]*N
    omegas = [0.559,0.344]*N
    
    xs2 = np.array(xs)
    Tcs2 = np.array(Tcs)
    Vcs2 = np.array(Vcs)
    omegas2 = np.array(omegas)
    assert_close(COSTALD_mixture(xs, T, Tcs, Vcs, omegas),
                 chemicals.numba.COSTALD_mixture(xs2, T, Tcs2, Vcs2, omegas2))

@mark_as_numba
def test_solbility():
    assert_close(Henry_converter(1.2e-5, old_scale='Hcp', new_scale='SI', rhom=55344.59,  MW=18.01528),
                 chemicals.numba.Henry_converter(1.2e-5, old_scale='Hcp', new_scale='SI', rhom=55344.59,  MW=18.01528))

@mark_as_numba
def test_refractivity():
    assert_close(brix_to_RI(5.8), chemicals.numba.brix_to_RI(5.8))

@mark_as_numba
def test_rachford_rice():
    n = 10
    zs = np.array([0.5, 0.3, 0.2]*n)
    Ks = np.array([1.685, 0.742, 0.532]*n)
    
    assert_close(chemicals.numba.Rachford_Rice_flash_error(0.5, zs=zs, Ks=Ks),
                 Rachford_Rice_flash_error(0.5, zs=zs, Ks=Ks))

    zs = np.array([0.5, 0.3, 0.2])
    Ks = np.array([1.685, 0.742, 0.532])
    VF_new, xs_new, ys_new = chemicals.numba.flash_inner_loop(zs=zs, Ks=Ks)
    VF, xs, ys = flash_inner_loop(zs=zs.tolist(), Ks=Ks.tolist())
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)

@mark_as_numba
def test_Rachford_Rice_solutionN():
    ns = [0.204322076984, 0.070970999150, 0.267194323384, 0.296291964579, 0.067046080882, 0.062489248292, 0.031685306730]
    Ks_y = [1.23466988745, 0.89727701141, 2.29525708098, 1.58954899888, 0.23349348597, 0.02038108640, 1.40715641002]
    Ks_z = [1.52713341421, 0.02456487977, 1.46348240453, 1.16090546194, 0.24166289908, 0.14815282572, 14.3128010831]
    ns2, Ks2, betas2 = np.array(ns), np.array([Ks_y, Ks_z]), np.array([.1, .6])
    betas_new, zs_new = chemicals.numba.Rachford_Rice_solutionN(ns2, Ks2, betas2)
    betas, zs = Rachford_Rice_solutionN(ns, [Ks_y, Ks_z], [.1, .6])
    assert_close1d(betas, betas_new, rtol=1e-14)
    assert_close2d(zs, zs_new, rtol=1e-14)

@mark_as_numba
def test_Rachford_Rice_solution2():
    ns = [0.204322076984, 0.070970999150, 0.267194323384, 0.296291964579, 0.067046080882, 0.062489248292, 0.031685306730]
    Ks_y = [1.23466988745, 0.89727701141, 2.29525708098, 1.58954899888, 0.23349348597, 0.02038108640, 1.40715641002]
    Ks_z = [1.52713341421, 0.02456487977, 1.46348240453, 1.16090546194, 0.24166289908, 0.14815282572, 14.3128010831]
    ns2, Ksy2, Ksz2 = np.array(ns), np.array(Ks_y), np.array(Ks_z)
    
    beta0_new, beta1_new, z0_new, z1_new, z2_new = chemicals.numba.Rachford_Rice_solution2(ns2, Ksy2, Ksz2, beta_y=.1, beta_z=.6)
    beta0, beta1, z0, z1, z2 = Rachford_Rice_solution2(ns, Ks_y, Ks_z, beta_y=.1, beta_z=.6)
    assert_close(beta0_new, beta0)
    assert_close(beta1_new, beta1)
    assert_close1d(z0, z0_new)
    assert_close1d(z1, z1_new)
    assert_close1d(z2, z2_new)
    

@mark_as_numba
def test_rachford_rice_polynomial():
    zs, Ks = [.4, .6], [2, .5]
    VF_new, xs_new, ys_new = chemicals.numba.Rachford_Rice_solution_polynomial(np.array(zs), np.array(Ks))
    VF, xs, ys = Rachford_Rice_solution_polynomial(zs, Ks)
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)

    zs = [0.5, 0.3, 0.2]
    Ks = [1.685, 0.742, 0.532]
    VF_new, xs_new, ys_new = chemicals.numba.Rachford_Rice_solution_polynomial(np.array(zs), np.array(Ks))
    VF, xs, ys = Rachford_Rice_solution_polynomial(zs, Ks)
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)
    
    zs = [0.2, 0.3, 0.4, 0.1]
    Ks = [2.5250, 0.7708, 1.0660, 0.2401]
    VF_new, xs_new, ys_new = chemicals.numba.Rachford_Rice_solution_polynomial(np.array(zs), np.array(Ks))
    VF, xs, ys = Rachford_Rice_solution_polynomial(zs, Ks)
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)
    
    zs = [0.2, 0.3, 0.4, 0.05, 0.05]
    Ks = [2.5250, 0.7708, 1.0660, 0.2401, 0.3140]
    VF_new, xs_new, ys_new = chemicals.numba.Rachford_Rice_solution_polynomial(np.array(zs), np.array(Ks))
    VF, xs, ys = Rachford_Rice_solution_polynomial(zs, Ks)
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)
    
    # 6 and higher use generic routine
    zs = [0.05, 0.10, 0.15, 0.30, 0.30, 0.10]
    Ks = [6.0934, 2.3714, 1.3924, 1.1418, 0.6457, 0.5563]
    VF_new, xs_new, ys_new = chemicals.numba.Rachford_Rice_solution_polynomial(np.array(zs), np.array(Ks))
    VF, xs, ys = Rachford_Rice_solution_polynomial(zs, Ks)
    assert_close(VF, VF_new)
    assert_close1d(xs, xs_new)
    assert_close1d(ys, ys_new)


@mark_as_numba
def test_lazy_loading():
    # Numba interfers with to_num
    # The data_reader functions are not part of the public API so are not converted
    chemicals.numba.heat_capacity.zabransky_dicts
    chemicals.numba.heat_capacity.CRC_standard_data
    
    assert 'jitclass' in str(chemicals.numba.heat_capacity.ZabranskySpline)
    assert 'jitclass' in str(chemicals.numba.heat_capacity.ZabranskyQuasipolynomial)
    assert 'jitclass' in str(chemicals.numba.heat_capacity.zabransky_dict_iso_s['2016-57-1'].models[0])

    
@mark_as_numba
def test_safety_functions():
    import test_safety
    swap_funcs_and_test(['NFPA_30_classification'], 
                        [chemicals.numba.NFPA_30_classification], 
                        test_safety.test_NFPA_30_classification)
    
        
