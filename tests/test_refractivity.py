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

from math import isnan
from fluids.numerics import assert_close, assert_close1d
import pytest
import pandas as pd
import numpy as np
from chemicals.refractivity import *
from chemicals.refractivity import RI_data_CRC_organic

def test_refractivity_CRC():
    assert_close(RI_data_CRC_organic['RI'].sum(), 6602.78821)
    assert_close(RI_data_CRC_organic['RIT'].sum(), 1304152.35)

@pytest.mark.slow
@pytest.mark.fuzz
def test_refractivity_all_answers():
    vals = [RI(i) for i in  RI_data_CRC_organic.index.values]
    RI_sum = sum([v[0] for v in vals])
    T_sum = sum([v[1] for v in vals if not isnan(v[1])])
    assert len(vals) == 4490
    assert type(vals[0][0]) is float
    assert type(vals[0][1]) is float
    
    assert_close(RI_sum, 6602.78821, rtol=1e-10)
    assert_close(T_sum, 1304152.35, rtol=1e-10)
    
    
def test_refractivity_general():
    vals = RI(CASRN='64-17-5')
    assert type(vals) is tuple
    assert_close1d(vals, (1.3611, 293.15))

    vals = RI_methods(CASRN='64-17-5')
    assert vals == ['CRC']
    assert RI_data_CRC_organic.index.is_unique
    assert RI_data_CRC_organic.shape == (4490, 2)    
    assert RI_methods(CASRN='6400000-17-5') == []

    with pytest.raises(Exception):
        RI(CASRN='64-17-5', method='FAIL')



def test_polarizability_from_RI():
    # Ethanol, with known datum RI and Vm
    alpha = polarizability_from_RI(1.3611, 5.8676E-5)
    assert_close(alpha, 5.147658123614415e-30)
    # Experimental value is 5.112 Angstrom^3 from cccbdb, http://cccbdb.nist.gov/polcalccomp2.asp?method=55&basis=9
    # Very good comparison.


def test_molar_refractivity_from_RI():
    # Ethanol, with known datum RI and Vm
    Rm = molar_refractivity_from_RI(1.3611, 5.8676E-5)
    assert_close(Rm, 1.2985217089649597e-05)
    # Confirmed with a value of 12.5355 cm^3/mol in http://rasayanjournal.co.in/vol-4/issue-4/38.pdf


def test_RI_from_molar_refractivity():
    RI = RI_from_molar_refractivity(1.2985e-5, 5.8676E-5)
    assert_close(RI, 1.3610932757685672)
    # Same value backwards

    assert_close(RI_from_molar_refractivity(molar_refractivity_from_RI(1.3611, 5.8676E-5), 5.8676E-5), 1.3611)


def test_RI_IAPWS():
    assert_close(RI_IAPWS(298.15, 997.047435, 0.5893), 1.3328581926471605, rtol=1e-12)


def test_RI_to_brix():
    assert_close(RI_to_brix(1.33299), 0.0)
    assert_close(RI_to_brix(1.532), 95.)
    assert_close(RI_to_brix(1.341452), 5.8)
    assert_close(RI_to_brix(1.3), -23.069930069929285)
    
    
def test_brix_to_RI():
    assert_close(brix_to_RI(5.8), 1.341452)
    assert_close(brix_to_RI(0), 1.33299)
    assert_close(brix_to_RI(95.0), 1.532)
    assert_close(brix_to_RI(RI_to_brix(1.3)), 1.3)
