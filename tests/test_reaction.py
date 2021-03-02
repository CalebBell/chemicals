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

from fluids.numerics import assert_close, assert_close1d, assert_close2d
import pytest
import pandas as pd

from chemicals.identifiers import check_CAS
from chemicals.reaction import *
from chemicals.heat_capacity import TRC_gas_data, CRC_standard_data
from chemicals.reaction import Hfg_API_TDB_data, Hfg_ATcT_data, Hfl_ATcT_data, Hfg_S0g_YAWS_data

def test_API_TDB_data():
    assert Hfg_API_TDB_data.index.is_unique
    assert Hfg_API_TDB_data['Hfg'].abs().sum() == 101711260
    assert Hfg_API_TDB_data.shape == (571, 2)
    assert all([check_CAS(i) for i in Hfg_API_TDB_data.index])


def test_ATcT_l():
    assert Hfl_ATcT_data.index.is_unique
    assert Hfl_ATcT_data.shape == (34,5)
    assert all([check_CAS(i) for i in Hfl_ATcT_data.index])
    tots_calc = [Hfl_ATcT_data[i].abs().sum() for i in ['Hfl_0K', 'Hfl', 'uncertainty']]
    tots = [2179500.0, 6819443, 19290]
    assert_close1d(tots_calc, tots)


def test_Hfg_ATcT_data():
    assert Hfg_ATcT_data.index.is_unique
    assert Hfg_ATcT_data.shape == (595, 5)
    assert all([check_CAS(i) for i in Hfg_ATcT_data.index])
    tots_calc = [Hfg_ATcT_data[i].abs().sum() for i in ['Hfg_0K', 'Hfg', 'uncertainty']]
    tots = [300788330, 300592764, 829204]
    assert_close1d(tots_calc, tots)

def test_Hfg_API_TDB_data():
    assert_close(Hfg('7732-18-5', method='API_TDB_G'), -241820.0)

    assert Hfg_methods('7732-18-5') == ['ATCT_G', 'CRC', 'API_TDB_G', 'TRC', 'YAWS']

    assert None == Hfg('98-00-1')

    with pytest.raises(Exception):
        Hf('98-00-0', method='BADMETHOD')

@pytest.mark.slow
def test_Hfg_API_TDB_data_fuzz():
    tot = sum([abs(Hfg(i, method='API_TDB_G')) for i in Hfg_API_TDB_data.index])
    assert_close(tot, 101711260.0)


def test_Hfl():
    Hfs = [Hfl('67-56-1'), Hfl('67-56-1', method='ATCT_L')]
    assert_close1d(Hfs, [-238400.0]*2)

    assert Hfl_methods('67-56-1') == ['ATCT_L', 'CRC']
    assert None == Hfl('98-00-1')

    tot = sum([abs(Hfl(i)) for i in Hfl_ATcT_data.index])
    assert_close(tot, 6819443.0)

    with pytest.raises(Exception):
        Hfl('98-00-0', method='BADMETHOD')


def test_Hfg():
    # default method ATCT_G
    assert_close(Hfg('7732-18-5'), -241822.0)

    Hfs = [Hfg('67-56-1', method=i) for i in Hfg_all_methods]
    assert_close1d(Hfs, [-200700., -190100., -201000., -200900.])

    assert Hfg_methods('67-56-1') == ['ATCT_G', 'CRC', 'API_TDB_G','TRC', 'YAWS']
    assert_close(-211800.0, Hfg('98-00-0'))

    with pytest.raises(Exception):
        Hfg('98-00-0', method='BADMETHOD')

def test_Hfs():
    assert_close(Hfs('101-81-5'), 71500)
    assert_close(Hfs('101-81-5', method='CRC'), 71500)
    assert ['CRC'] == Hfs_methods('101-81-5')



@pytest.mark.fuzz
@pytest.mark.slow
def test_Hfg_all_values():
    tot1 = sum([abs(Hfg(i, method='TRC')) for i in TRC_gas_data.index[pd.notnull(TRC_gas_data['Hfg'])]])
    assert_close(tot1, 495689880.0)

    tot2 = sum([abs(Hfg(i, method='ATCT_G')) for i in Hfg_ATcT_data.index])
    assert_close(tot2, 300592764.0)

    tot3 = sum([abs(Hfg(i, method='YAWS')) for i in Hfg_S0g_YAWS_data.index[pd.notnull(Hfg_S0g_YAWS_data['Hfg'])]])
    assert_close(tot3, 1545148533.0)

    tot4 = sum([abs(Hfg(i, method='CRC')) for i in CRC_standard_data.index[pd.notnull(CRC_standard_data['Hfg'])]])
    assert_close(tot4, 392946600.0)

def test_S0g():
    S0s = [S0g('7732-18-5', method=i) for i in S0g_all_methods]
    assert_close1d(S0s, [188.8, 188.84])

    assert S0g_methods('67-56-1') == ['CRC', 'YAWS']

    assert_close(239.9, S0g('67-56-1'))

    with pytest.raises(Exception):
        S0g('98-00-0', method='BADMETHOD')

@pytest.mark.fuzz
@pytest.mark.slow
def test_S0g_all_values():
    tot3 = sum([abs(S0g(i, method='YAWS')) for i in Hfg_S0g_YAWS_data.index[pd.notnull(Hfg_S0g_YAWS_data['S0g'])]])
    assert_close(tot3, 2691892.382999995)

    tot4 = sum([abs(S0g(i, method='CRC')) for i in CRC_standard_data.index[pd.notnull(CRC_standard_data['S0g'])]])
    assert_close(tot4, 141558.30000000008)


def test_S0s():
    assert_close(S0s('7439-93-2'), 29.1) # Lithium
    assert_close(S0s('7439-93-2', method='CRC'), 29.1)

    methods = S0s_methods('7439-93-2')
    assert methods == ['CRC']

def test_S0l():
    assert_close(S0l('7439-97-6'), 75.9) # Lithium
    assert_close(S0l('7439-97-6', method='CRC'), 75.9)

    methods = S0l_methods('7439-97-6')
    assert methods == ['CRC']

def test_Gibbs_formation():
    Gf =  Gibbs_formation(-285830.0, 69.91,  [0.0, 0.0], [130.571, 205.147], [1.0, .5])
    assert_close(Gf, -237161.633825)

    Gf = Gibbs_formation(-241818, 188.825,  [0.0, 0], [130.571, 205.147], [1.0, .5])
    assert_close(Gf, -228604.141075)

    Gf = Gibbs_formation(-648980, 297.713, [0.0, 0.0, 0.0], [5.74, 152.206, 202.789], [1, .5, 1.5])
    assert_close(Gf, -622649.329975)


def test_Hf_basis_converter():
    assert_close(Hf_basis_converter(44018.0, Hf_liq=-285830.0), -241812)

    assert_close(Hf_basis_converter(44018, Hf_gas=-241812.0), -285830)

    with pytest.raises(ValueError):
        Hf_basis_converter(44018, Hf_liq=None)
    with pytest.raises(ValueError):
        Hf_basis_converter(2000, Hf_gas=None, Hf_liq=None)
    with pytest.raises(ValueError):
        Hf_basis_converter(Hvapm=-1, Hf_liq=1)
    with pytest.raises(ValueError):
        Hf_basis_converter(Hvapm=None, Hf_liq=1)

def test_entropy_formation():
    Sf = entropy_formation(Hf=-74520.0, Gf=-50490.0)
    assert_close(Sf, -80.59701492537314)

    Sf = entropy_formation(Hf=-241818, Gf=-228572)
    assert_close(Sf, -44.427301693778304)



def test_balance_stoichiometry():
    test_cases = [
    [[{'Hg': 1, 'O': 1}, {u'Hg': 1}, {'O': 2}], [True, False, False], [2.0, 2.0, 1.0]],
    [[{'Cl': 2}, {'C': 3, 'H': 6}, {'C': 3, 'Cl': 1, 'H': 5}, {'Cl': 1, 'H': 1}],
      [True, True, False, False, False],
      [1, 1, 1, 1]],
    [[{'Al': 1}, {'H': 1, 'N': 1, 'O': 3}, {'Al': 1, 'N': 3, 'O': 9}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
      [True, True, False, False, False],
      [1.0, 4.0, 1.0, 1.0, 2.0]],
    [[{'Fe': 1}, {'O': 2}, {'Fe':2, 'O': 3}], [True, True, False], [4.0, 3.0, 2.0]],
    [[{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False], [4.0, 5.0, 4.0, 6.0]],
    [[{'O': 2}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}, {'C': 6, 'H': 14}], [True, False, False, True], [19.0, 14.0, 12.0, 2.0]],

    ]

    for atomss, statuses, products in test_cases:
        assert_close1d(balance_stoichiometry(stoichiometric_matrix(atomss, statuses)), products)


def test_stoichiometric_matrix():
    res = stoichiometric_matrix([{'Mg': 1, 'O': 1}, {'Mg': 1}, {'O': 2}], [True, False, False])
    assert_close2d([[1, -1, 0], [1, 0, -2]], res)