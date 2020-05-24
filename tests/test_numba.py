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
from fluids.constants import *
from fluids.numerics import assert_close, assert_close1d
import pytest
try:
    import numba
    import chemicals.numba
#    import chemicals.numba_vectorized
except:
    numba = None
import numpy as np

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
    
@pytest.mark.skipif(numba is None, reason="Numba is missing")    
def test_visc_misc():
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