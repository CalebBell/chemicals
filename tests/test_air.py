# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell
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
import chemicals
from chemicals.air import *
from math import *
from fluids.numerics import assert_close, assert_close1d, assert_close2d, linspace, logspace, derivative

from chemicals.air import TAU_MAX_EXP_87


def test_lemmon2000_A0():
    assert_close(lemmon2000_air_A0(0.36842, 0.5), -17.026123512818458, rtol=1e-15)
    
    # if statement test
    assert_close(lemmon2000_air_A0(TAU_MAX_EXP_87*(1-1e-15), .5), 
                 lemmon2000_air_A0(TAU_MAX_EXP_87*(1+1e-15), .5), rtol=1e-15)
    
    # Value computed with SymPy
    assert_close(lemmon2000_air_A0(300, .5), -3.54214522586700104, rtol=5e-15)
    
    
    # Check points that may fail
    lemmon2000_air_A0(1e20, 1e20)
    lemmon2000_air_A0(1e-10, 1e-10)
    
    
    
def test_lemmon2000_air_dA0_dtau():
    assert_close(derivative(lambda tau: lemmon2000_air_A0(tau, 13000/11183.9), 126.192/200.0, dx=1e-7),
                 lemmon2000_air_dA0_dtau(126.192/200.0, 13000/11183.9), rtol=1e-9)
    
    assert_close(lemmon2000_air_dA0_dtau(0.36842, 0.5), 6.764336610288353, rtol=1e-14)
    
    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_A0(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_dA0_dtau(rat, .5))
    
    
def test_lemmon2000_air_d2A0_dtau2():
    assert_close(lemmon2000_air_d2A0_dtau2(126.192/200.0, 13000/11183.9), -6.260482844274295, rtol=1e-13)

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_dA0_dtau(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d2A0_dtau2(rat, .5))
    
    
def test_lemmon2000_air_d3A0_dtau3():
    assert_close(lemmon2000_air_d3A0_dtau3(0.36842, .5), 102.9144884392338, rtol=2e-12)
    assert_close(derivative(lambda tau: lemmon2000_air_d2A0_dtau2(tau, .5), 0.36842, dx=1e-7),
                 lemmon2000_air_d3A0_dtau3(0.36842, .5))

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0, 2.5, 2.0, 1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d2A0_dtau2(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d3A0_dtau3(rat, .5))
    
def test_lemmon2000_air_d4A0_dtau4():
    assert_close(lemmon2000_air_d4A0_dtau4(126.192/200.0, 13000/11183.9), -94.8155327278803, rtol=1e-13)
    
    assert_close(derivative(lambda tau: lemmon2000_air_d3A0_dtau3(tau, .5), 0.36842, dx=4e-7),
                 lemmon2000_air_d4A0_dtau4(0.36842, .5))

    for rat in (1000.0, 100.0, 10.0, 5.0, 3.0,2.5, 2.0,  1.5, 1.0, .5, .2, .1, .01, .001):
        assert_close(derivative(lambda tau: lemmon2000_air_d3A0_dtau3(tau, .5), rat, dx=rat*1e-5),
                     lemmon2000_air_d4A0_dtau4(rat, .5))
        
    
def test_lemmon2000_air_Ar():
    assert_close(lemmon2000_air_Ar(0.36842, 0.15880050154579475), 0.004798812280624336, rtol=1e-13)