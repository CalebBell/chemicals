# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from fluids.numerics import assert_close
import pytest
import pandas as pd
from chemicals.acentric import *
from chemicals.critical import critical_data_PSRKR4, critical_data_PassutDanner, critical_data_Yaws


@pytest.mark.fuzz
@pytest.mark.slow
def test_acentric_fuzz():
    sources = [critical_data_PSRKR4, critical_data_PassutDanner, critical_data_Yaws]

    CASs = set()
    [CASs.update(list(k.index.values)) for k in sources]

    # Use the default method for each chemical in this file
    omegas = [omega(i) for i in CASs] # This is quite slow
    assert len([i for i in omegas if i is not None]) == 6306

    omegas_default_sum = sum(abs(i) for i in omegas if i is not None)
    assert_close(omegas_default_sum, 3244.876299999988)

def test_acentric_main():
    omega_calc = omega('629-92-5', method='PSRK')
    assert_close(omega_calc, 0.8486, rtol=1e-13)
    assert_close(omega('629-92-5'), 0.8486, rtol=1e-13)

    omega_calc = omega('629-92-5', method='PD')
    assert_close(omega_calc, 0.8271, rtol=1e-13)

    omega_calc = omega('629-92-5', method='YAWS')
    assert_close(omega_calc, 0.852, rtol=1e-13)

    methods = omega_methods('74-98-6')
    assert methods == ['PSRK', 'PD', 'YAWS']

    # Error handling
    assert None == omega(CASRN='BADCAS')

@pytest.mark.skip_types
def test_acentric_main_bad_methods():
    with pytest.raises(ValueError):
        omega(CASRN='98-01-1', method='BADMETHOD')
    with pytest.raises(ValueError):
        omega(CASRN='98-01-1', method=1)
    with pytest.raises(ValueError):
        omega(CASRN='98-01-1', method=omega)




def test_acentric_correlation():
    omega = LK_omega(425.6, 631.1, 32.1E5)
    assert_close(omega, 0.32544249926397856)


def test_omega_definition():
    omega = omega_definition(999542.0, 22048320.0)
    assert_close(omega, 0.3435744558761711, rtol=1e-12)

def test_Stiel_polar_factor():
    factor = Stiel_polar_factor(Psat=169745.0, Pc=22048321.0, omega=0.344)
    assert_close(0.02322146744772713, factor)

