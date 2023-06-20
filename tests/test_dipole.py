"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from chemicals import dipole
from chemicals.dipole import dipole_moment, dipole_moment_methods


def test_dipole_moment_methods():

    tot = dipole.dipole_data_Poling['dipole_moment'].sum()
    assert_close(tot, 248.59999999999999)

    tot = dipole.dipole_data_CCDB['dipole_moment'].sum()
    assert_close(tot, 629.69)

    tot = dipole.dipole_data_Muller['dipole_moment'].sum()
    assert_close(tot, 420.05190108045235)

def test_dipole():
    d = dipole_moment(CASRN='64-17-5')
    assert_close(d, 1.44)

    d = dipole_moment(CASRN='75-52-5', method='POLING')
    assert_close(d, 3.1)

    d = dipole_moment(CASRN='56-81-5', method='MULLER')
    assert_close(d, 4.21)

    methods = dipole_moment_methods(CASRN='78-78-4')
    methods_fixed = ['CCCBDB', 'MULLER', 'POLING']
    assert methods == methods_fixed

    with pytest.raises(Exception):
        dipole_moment(CASRN='78-78-4', method='FAIL')


    # Just test one thing with psi4
    # git can track everything
    val = dipole_moment(CASRN='64-17-5', method='PSI4_2022A')
    assert val > 1
