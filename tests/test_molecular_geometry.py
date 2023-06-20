"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2022 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
from fluids.numerics import assert_close

from chemicals import molecular_geometry
from chemicals.molecular_geometry import RG, RG_methods, linear, linear_methods


def test_molecular_geometry_methods():

    assert len(molecular_geometry.radius_of_gyration_data_psi4_2022a['RG']) > 200

def test_RG():
    val = RG(CASRN='64-17-5')
    assert_close(val, 2.225708e-10, rtol=1e-1)


    methods = RG_methods(CASRN='64-17-5')
    methods_fixed = ['PSI4_2022A', 'CHEMSEP']
    assert methods == methods_fixed

    with pytest.raises(Exception):
        RG(CASRN='78-78-4', method='FAIL')


def test_linear():
    assert linear(CASRN='64-17-5') is False


    assert linear_methods(CASRN='64-17-5')  == ['PSI4_2022A']
