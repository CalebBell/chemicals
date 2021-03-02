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
import types
from fluids.numerics import assert_close, assert_close1d
import pytest
import chemicals
from chemicals.units import *

def assert_pint_allclose(value, magnitude, units, rtol=1e-7, atol=0):
    assert_close(value.to_base_units().magnitude, magnitude, rtol=rtol, atol=atol)
    if type(units) != dict:
        units = dict(units.dimensionality)
    assert dict(value.dimensionality) == units


def test_misc_units():
    Hf = LHV_from_HHV(-726024.0*u.J/u.mol, 2)
    assert_pint_allclose(Hf, -638001.008, {'[length]': 2, '[mass]': 1, '[substance]': -1, '[time]': -2})

def test_speed_of_sound():
    v = speed_of_sound(V=0.00229754*u.m**3/u.mol, dP_dV=-3.5459e+08*u("Pa*mol/m^3"), Cp=153.235*u("J/mol/K"), Cv=132.435*u("J/mol/K"), MW=67.152*u("g/mol"))

    assert_pint_allclose(v, 179.5868138460819, {'[length]': 1, '[time]': -1})

    v = speed_of_sound(V=0.00229754*u.m**3/u.mol, dP_dV=-3.5459e+08*u("Pa*mol/m^3"), Cp=153.235*u("J/mol/K"), Cv=132.435*u("J/mol/K"))

    assert_pint_allclose(v, 46.537593457316525, {'[length]': 1, '[mass]': 0.5, '[substance]': -0.5, '[time]': -1})


def test_heat_capacity_variable_units():
    S = Lastovka_solid_integral_over_T(300*u.K, 0.2139*u.mol/u.g)
    assert_pint_allclose(S, 1947.5537561495564, {'[length]': 2, '[temperature]': -1, '[time]': -2})

    S = Lastovka_solid_integral_over_T(300*u.K, 0.2139*u.mol/u.g, MW=48*u.g/u.mol)
    assert_pint_allclose(S, 93.4825802951787, {'[length]': 2, '[mass]': 1, '[substance]': -1, '[temperature]': -1, '[time]': -2})
