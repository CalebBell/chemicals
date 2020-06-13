# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.'''

from __future__ import division
from chemicals import *
import chemicals.vectorized
from math import *
from random import random
from fluids.constants import *
from fluids.numerics import assert_close, assert_close1d
import pytest
try:
    import numba
    import chemicals.numba
    import chemicals.numba_vectorized
except:
    numba = None
import numpy as np



@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
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


@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_return_2d_array():
    d2xs = [[0.152, 0.08, 0.547], [0.08, 0.674, 0.729], [0.547, 0.729, 0.131]]
    xs = [0.7, 0.2, 0.1]
    dxdn_partials = d2xs_to_dxdn_partials(d2xs, xs)
    a, b = np.array(d2xs), np.array(xs)
    dxdn_partials_np = chemicals.numba.d2xs_to_dxdn_partials(a, b)
    
    assert type(dxdn_partials_np) is np.ndarray
    assert_close1d(dxdn_partials, dxdn_partials_np)



@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_mixing_simple():
    a = np.array([1,2])
    b = np.array([.1, .2])
    tot = chemicals.numba.mixing_simple(a, b)
    assert_close(tot, 0.5, rtol=1e-14)
    
    
    a = np.array([.1, .9])
    b = np.array([.01, .02])
    val = chemicals.numba.mixing_logarithmic(a, b)
    assert_close(val, 0.01866065983073615, rtol=1e-13)

@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")    
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

@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")    
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
    
@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")    
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
        
        
#@pytest.mark.numba
#@pytest.mark.skipif(numba is None, reason="Numba is missing")
#def test_virial():
#    # Takes 8 seconds to compile. Fun!
#    assert_close(chemicals.numba.BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469),
#                 chemicals.BVirial_Tsonopoulos_extended(430., 405.65, 11.28E6, 0.252608, a=0, b=0, species_type='ketone', dipole=1.469),
#                rtol=1e-13)

@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_phase_change():
    # Function had some duplicated powers; numba was optimizing them on me anyway
    # Had list-in-list constants being indexed. I thought that would take a lot of time
    # but instead removing it only saved 25%, and ~8% in CPython, and zilch in PyPy.
    # PyPy takes 19% of the time numba does here, numba has a high overhead.
    assert_close(chemicals.numba.MK(553.15, 751.35, 0.302),
             chemicals.MK(553.15, 751.35, 0.302), rtol=1e-12)


@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_vapor_pressure():
    # PyPy 75 ns, CPython 2470 ns, numba 214 ns
    assert_close(chemicals.numba.dPsat_IAPWS_dT(300.), 
                 chemicals.dPsat_IAPWS_dT(300.), rtol=1e-14)
    
    
    Psats_vec_expect = [34478.367349639906, 33596697.716487624, 109799836.81382856, 179376011.49286702, 234627689.09298804]
    Ts = np.linspace(100, 1000, 5)
    Psats_calc = chemicals.numba_vectorized.Antoine(Ts, 8.7687, 395.744, -6.469, 10)
    assert_close(Psats_calc, Psats_vec_expect, rtol=1e-11)
    

@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_temperature():
    # Note also the last four decimals are different!
    # 494 us numba, 388 us PyPy, 1740 us CPython
    assert_close(chemicals.numba.ITS90_68_difference(1000.),
                 chemicals.ITS90_68_difference(1000.0), rtol=1e-12)
    
    # Probably never going to work
#    chemicals.numba.T_converter(500, 'ITS-68', 'ITS-48')
    
@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_volume():
    assert_close(chemicals.numba.Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245),
                chemicals.Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245))
    assert_close(chemicals.numba.COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532),
                 chemicals.COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532))

    assert_close(chemicals.numba.Bhirud_normal(280.0, 469.7, 33.7E5, 0.252),
                 Bhirud_normal(280.0, 469.7, 33.7E5, 0.252))
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
    # timing after optimization at 200 elements - 1.49 m CPython, 27.1 µs numba, 63.5 µs PyPy
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
    
@pytest.mark.numba
@pytest.mark.skipif(numba is None, reason="Numba is missing")
def test_rachford_rice():
    n = 10
    zs = np.array([0.5, 0.3, 0.2]*n)
    Ks = np.array([1.685, 0.742, 0.532]*n)
    
    assert_close(chemicals.numba.Rachford_Rice_flash_error(0.5, zs=zs, Ks=Ks),
                 Rachford_Rice_flash_error(0.5, zs=zs, Ks=Ks))
