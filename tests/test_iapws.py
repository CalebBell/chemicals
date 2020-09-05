# -*- coding: utf-8 -*-
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

from numpy.testing import assert_allclose
import pytest
from math import *
from chemicals.iapws import *
from fluids.numerics import assert_close

def test_iapws97_dG_dpi_region1():
    assert_close(iapws97_dG_dpi_region1(1386/277.15, 101325/16.53E6),
                 0.12923271825448354, rtol=1e-14)


def test_iapws97_dG_dpi_region2():
    assert_close(iapws97_dGr_dpi_region2(.656, 16), -0.006292631931275252, rtol=1e-14)

def test_iapws97_dG_dpi_region5():
    assert_close(iapws97_dGr_dpi_region5(.5, 30.0), 0.0004009761854002751, rtol=1e-14)

def test_iapws_boundary_equations():
    assert_close(iapws97_boundary_2_3(0.623150000E3), 16529164.252621626, rtol=1e-13)
    
    assert_close(iapws97_boundary_3uv(22.3E6), 647.7996121480069, rtol=1e-14)
    
    assert_close(iapws97_boundary_3ef(40E6), 713.959399239744, rtol=1e-14)
    
    assert_close(iapws97_boundary_3cd(25E6), 649.3659208321279, rtol=1e-14)
    
    assert_close(iapws97_boundary_3gh(25E6), 656.6980572261236, rtol=2e-14)
    
    assert_close(iapws97_boundary_3ij(25E6), 660.7865756716819, rtol=1e-14)
    
    assert_close(iapws97_boundary_3jk(25E6), 668.1915358826951, rtol=1e-14)
    
    assert_close(iapws97_boundary_3mn(22.8E6), 649.6054132953997, rtol=1e-14)
    
    assert_close(iapws97_boundary_3qu(22E6), 645.6355027340121, rtol=1e-14)
    
    assert_close(iapws97_boundary_3rx(22E6), 648.2622753670172, rtol=1e-14)
    
    assert_close(iapws97_boundary_3wx(log(22.3), 1 / log(22.3)), 648.204947950734, rtol=1e-14)
    
    assert_close(iapws97_boundary_3ab(log(40), 1 / log(40)), 693.0341408296053, rtol=1e-14)
    
    assert_close(iapws97_boundary_3op(log(22.8), 1 / log(22.8)), 650.010694314133, rtol=1e-14)
    
