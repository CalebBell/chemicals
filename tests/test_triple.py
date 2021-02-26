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

from fluids.numerics import assert_close, assert_close1d
import pytest
import pandas as pd
from chemicals.triple import *
from chemicals.triple import triple_data_Staveley


def test_data():
    Tt_sum = triple_data_Staveley['Tt68'].sum()
    assert_close(Tt_sum, 31251.845000000001)

    Pt_sum = triple_data_Staveley['Pt'].sum()
    assert_close(Pt_sum, 1886624.8374376972)

    Pt_uncertainty_sum = triple_data_Staveley['Pt_uncertainty'].sum()
    assert_close(Pt_uncertainty_sum, 138.65526315789461)

    assert triple_data_Staveley.index.is_unique
    assert triple_data_Staveley.shape == (189, 5)


def test_Tt():
    Tt1_calc = Tt('7664-41-7')
    Tt1 = 195.48
    Tt2_calc = Tt('74-82-8', method='MELTING')
    Tt2 = 90.75
    Tt3_calc = Tt('74-82-8')
    Tt3 = 90.69
    assert_close1d([Tt1_calc, Tt2_calc, Tt3_calc], [Tt1, Tt2, Tt3])

    m = Tt_methods('7439-90-9')
    assert m == ['STAVELEY', 'MELTING']
    assert None == Tt('72433223439-90-9')
    with pytest.raises(Exception):
        Tt('74-82-8', method='BADMETHOD')


@pytest.mark.slow
def test_Tt_fuzz():
    Tt_sum = sum([Tt(i) for i in triple_data_Staveley.index])
    assert_close(Tt_sum, 31251.845000000001)
    Tt_sum2 = pd.Series([Tt(i, method='MELTING') for i in triple_data_Staveley.index]).sum()
    assert_close(Tt_sum2, 28778.196000000004)


def test_Pt():
    Pt1_calc = Pt('7664-41-7')
    Pt1 = 6079.5
    assert_close(Pt1_calc, Pt1)

    m = Pt_methods('7664-41-7')
    assert m == ['STAVELEY']
    assert None == Pt('72433223439-90-9')
    with pytest.raises(Exception):
        Pt('74-82-8', method='BADMETHOD')

@pytest.mark.slow
def test_Pt_fuzz():
    Pt_sum = sum([Pt(i) for i in triple_data_Staveley.index if pd.notnull(triple_data_Staveley.at[i, 'Pt'])])
    assert_close(Pt_sum, 1886624.8374376972)
