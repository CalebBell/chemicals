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

import pandas as pd
import pytest
from fluids.numerics import assert_close, assert_close1d, assert_close2d, numpy as np

from chemicals.heat_capacity import CRC_standard_data, TRC_gas_data
from chemicals.reaction import (
    Gibbs_formation,
    Hf_basis_converter,
    Hfg,
    Hfg_all_methods,
    Hfg_API_TDB_data,
    Hfg_ATcT_data,
    Hfg_methods,
    Hfg_S0g_YAWS_data,
    Hfl,
    Hfl_ATcT_data,
    Hfl_methods,
    Hfs,
    Hfs_methods,
    S0g,
    S0g_all_methods,
    S0g_methods,
    S0l,
    S0l_methods,
    S0s,
    S0s_methods,
    balance_stoichiometry,
    entropy_formation,
    standard_formation_reaction,
    stoichiometric_matrix,
    round_to_significant,
)


def test_API_TDB_data():
    assert Hfg_API_TDB_data['Hfg'].abs().sum() == 101711260
    assert Hfg_API_TDB_data.shape == (571, 2)


def test_ATcT_l():
    assert Hfl_ATcT_data.shape == (34,5)
    tots_calc = [Hfl_ATcT_data[i].abs().sum() for i in ['Hfl_0K', 'Hfl', 'uncertainty']]
    tots = [2179500.0, 6819443, 19290]
    assert_close1d(tots_calc, tots)


def test_Hfg_ATcT_data():
    assert Hfg_ATcT_data.shape == (595, 5)
    tots_calc = [Hfg_ATcT_data[i].abs().sum() for i in ['Hfg_0K', 'Hfg', 'uncertainty']]
    tots = [300788330, 300592764, 829204]
    assert_close1d(tots_calc, tots)

def test_Hfg_API_TDB_data():
    assert_close(Hfg('7732-18-5', method='API_TDB_G'), -241820.0)

    assert Hfg_methods('7732-18-5') == ['ATCT_G', 'CRC', 'API_TDB_G', 'WEBBOOK', 'TRC', 'JANAF', 'YAWS']

    assert None is Hfg('98-00-1')

    with pytest.raises(Exception):
        Hfg('98-00-0', method='BADMETHOD')

@pytest.mark.slow
def test_Hfg_API_TDB_data_fuzz():
    tot = sum([abs(Hfg(i, method='API_TDB_G')) for i in Hfg_API_TDB_data.index])
    assert_close(tot, 101711260.0)


def test_Hfl():
    Hfs = [Hfl('67-56-1'), Hfl('67-56-1', method='ATCT_L')]
    assert_close1d(Hfs, [-238400.0]*2)

    assert Hfl_methods('67-56-1') == ['ATCT_L', 'CRC', 'WEBBOOK']
    assert None is Hfl('98-00-1')

    tot = sum([abs(Hfl(i)) for i in Hfl_ATcT_data.index])
    assert_close(tot, 6819443.0)

    with pytest.raises(Exception):
        Hfl('98-00-0', method='BADMETHOD')


def test_Hfg():
    # default method ATCT_G
    assert_close(Hfg('7732-18-5'), -241822.0)

    Hfs = [Hfg('67-56-1', method=i) for i in Hfg_all_methods]
    assert_close1d(Hfs, [-200700.0, -190100.0, -201000.0, -205000.0, None, -200900.0, -216200.0])

    assert Hfg_methods('67-56-1') == ['ATCT_G', 'CRC', 'API_TDB_G', 'WEBBOOK', 'TRC', 'YAWS', 'JOBACK']
    assert_close(-211800.0, Hfg('98-00-0'))

    with pytest.raises(Exception):
        Hfg('98-00-0', method='BADMETHOD')

def test_Hfs():
    assert_close(Hfs('101-81-5'), 71500)
    assert_close(Hfs('101-81-5', method='CRC'), 71500)
    assert ['CRC', 'WEBBOOK'] == Hfs_methods('101-81-5')



@pytest.mark.fuzz
@pytest.mark.slow
def test_Hfg_all_values():
    tot1 = sum([abs(Hfg(i, method='TRC')) for i in TRC_gas_data.index[pd.notnull(TRC_gas_data['Hfg'])]])
    assert_close(tot1, 495689880.0)

    tot2 = sum([abs(Hfg(i, method='ATCT_G')) for i in Hfg_ATcT_data.index])
    assert_close(tot2, 300592764.0)

    tot3 = sum([abs(Hfg(i, method='YAWS')) for i in Hfg_S0g_YAWS_data.index[pd.notnull(Hfg_S0g_YAWS_data['Hfg'])]])
    assert_close(tot3, 1544220403.0)

    tot4 = sum([abs(Hfg(i, method='CRC')) for i in CRC_standard_data.index[pd.notnull(CRC_standard_data['Hfg'])]])
    assert_close(tot4, 392946600.0)

def test_S0g():
    S0s = [S0g('7732-18-5', method=i) for i in S0g_all_methods]
    assert_close1d(S0s, [188.8, 188.83842, 188.834, 188.84])

    assert S0g_methods('67-56-1') == ['CRC', 'YAWS']

    assert_close(239.9, S0g('67-56-1'))

    with pytest.raises(Exception):
        S0g('98-00-0', method='BADMETHOD')

@pytest.mark.fuzz
@pytest.mark.slow
def test_S0g_all_values():
    tot3 = sum([abs(S0g(i, method='YAWS')) for i in Hfg_S0g_YAWS_data.index[pd.notnull(Hfg_S0g_YAWS_data['S0g'])]])
    assert_close(tot3, 2690113.4130000058)

    tot4 = sum([abs(S0g(i, method='CRC')) for i in CRC_standard_data.index[pd.notnull(CRC_standard_data['S0g'])]])
    assert_close(tot4, 141558.30000000008)


def test_S0s():
    assert_close(S0s('7439-93-2'), 29.1) # Lithium
    assert_close(S0s('7439-93-2', method='CRC'), 29.1)

    methods = S0s_methods('7439-93-2')
    assert methods == ['CRC', 'WEBBOOK']

def test_S0l():
    assert_close(S0l('7439-97-6'), 75.9) # Lithium
    assert_close(S0l('7439-97-6', method='CRC'), 75.9)

    methods = S0l_methods('7439-97-6')
    assert methods == ['CRC', 'WEBBOOK']

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


def check_reaction_balance(matrix, products_calc, atol=1e-13):
    """Check that coefficients satisfy the stoichiometric matrix equation within tolerance."""
    result = np.array(matrix) @ np.array(products_calc)
    assert_close1d(result, [0.0]*len(result), atol=atol)

ill_conditioned_stoich_test_cases = [
    (
        [[{'C': 100000, 'H': 200000, 'N': 1}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}],
         [True, True, False, False, False],
         [1, 150001, 100000, 100000, 1]],
        "C100000H200000N + O2 = CO2 + H2O + NO2"
    ),
    (
        [[{'C': 100000000, 'H': 200000000, 'N': 1}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}],
         [True, True, False, False, False],
         [1, 150000001, 100000000, 100000000, 1]],
        "C100000000H200000000N + O2 = CO2 + H2O + NO2"
    ),
    (
        [[{'C': 50000, 'H': 100000, 'O': 25000}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
         [True, True, False, False],
         [1.0, 62500.0, 50000.0, 50000.0]],
        "C50000H100000O25000 + O2 = CO2 + H2O"
    ),
    (
        [[{'Fe': 1000, 'C': 2000, 'N': 2000}, {'H': 2, 'O': 2}, {'Fe': 2, 'O': 3}, {'H': 1, 'C': 1, 'N': 1}, {'H': 2, 'O': 1}],
         [True, True, False, False, False],
         [1.0, 500.0, 500.0, 2000.0, -500.0]],
        "Fe1000(CN)2000 + H2O2 = Fe2O3 + HCN + H2O"
    ),
    (
        [[{'Au': 1000, 'Cu': 1000}, {'H': 1, 'N': 1, 'O': 3}, {'Au': 1}, {'Cu': 1, 'N': 2, 'O': 6}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
         [True, True, False, False, False, False],
         [3, 8000, 3000, 3000, 2000, 4000]],
        "Au1000Cu1000 + HNO3 = Au + Cu(NO3)2 + NO + H2O"
    ),
    (
        [[{'H': 100000, 'S': 50000, 'O': 200000}, {'S': 1, 'O': 3}, {'H': 2, 'O': 1}],
         [True, False, False],
         [1.0, 50000.0, 50000.0]],
        "H100000(SO4)50000 = SO3 + H2O"
    ),
    (
        [[{'C': 50000, 'H': 100000, 'N': 10000, 'O': 20000, 'S': 1000}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}, {'S': 1, 'O': 2}],
         [True, True, False, False, False, False],
         [1.0, 76000.0, 50000.0, 50000.0, 10000.0, 1000.0]],
        "C50000H100000N10000O20000S1000 + O2 = CO2 + H2O + NO2 + SO2"
    ),
    (
        [[{'C': 12594, 'H': 25422, 'O': 5004}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
         [True, True, False, False],
         [2.0, 32895.0, 25188.0, 25422.0]],
        "PG5 combustion"
    ),
    (
        [[{'C': 597000000, 'H': 744000002, 'N': 225000000, 'O': 390000001, 'P': 60000000}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}, {'P': 4, 'O': 10}],
         [True, True, False, False, False, False],
         [1.0, 888000000.0, 597000000.0, 372000001.0, 225000000.0, 15000000.0]],
        "Y chromosome DNA combustion: C597000000H744000002N225000000O390000001P60000000 + O2 = CO2 + H2O + NO2 + P4O10"
    ),
    (
        [[{'C': 597000000, 'H': 744000002, 'N': 225000000, 'O': 390000001, 'P': 60000000}, {'H': 2, 'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 2}, {'P': 2, 'O': 5}],
         [True, True, False, False, False, False],
         [1.0, 1326000000.0, 597000000.0, 1698000000.0, 112500000.0, 30000000.0]],
        "Y chromosome with H2O2: C597000000H744000002N225000000O390000001P60000000 + H2O2 = CO2 + H2O + N2 + P2O5"
    ),
    (
        [[{'C': 597000000, 'H': 744000002, 'N': 225000000, 'O': 390000001, 'P': 60000000}, {'H': 2, 'S': 1, 'O': 4}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}, {'H': 1, 'P': 1, 'O': 4}, {'S': 1, 'O': 2}],
         [True, True, False, False, False, False, False],
         [1.0, 1836000000.0, 597000000.0, 2178000000.0, 225000000.0, 60000000.0, 1836000000.0]],
        "Y chromosome with H2SO4: C597000000H744000002N225000000O390000001P60000000 + H2SO4 = CO2 + H2O + NO2 + HPO4 + SO2"
    ),
    (
        [[{'C': 597000000, 'H': 744000002, 'N': 225000000, 'O': 390000001, 'P': 60000000}, {'O': 2}, {'C': 1, 'O': 1}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}, {'P': 4, 'O': 10}],
         [True, True, False, False, False, False],
         [1.0, 589500000.0, 597000000.0, 372000001.0, 225000000.0, 15000000.0]],
        "Y chromosome incomplete combustion: C597000000H744000002N225000000O390000001P60000000 + O2 = CO + H2O + NO2 + P4O10"
    ),
    (
        [[{'C': 597000000, 'H': 744000002, 'N': 225000000, 'O': 390000001, 'P': 60000000}, {'Cl': 2}, {'C': 1, 'Cl': 4}, {'H': 1, 'Cl': 1}, {'N': 1, 'Cl': 3}, {'P': 1, 'Cl': 5}, {'O': 2}],
         [True, True, False, False, False, False, False],
         [2.0, 4107000002.0, 1194000000.0, 1488000004.0, 450000000.0, 120000000.0, 390000001.0]],
        "Y chromosome chlorination: C597000000H744000002N225000000O390000001P60000000 + Cl2 = CCl4 + HCl + NCl3 + PCl5 + O2"
    )
]
@pytest.mark.parametrize("test_case,test_name", ill_conditioned_stoich_test_cases)
def test_balance_stoichiometry_ill_conditioned(test_case, test_name):
    """Test stoichiometry balancing for ill-conditioned reactions."""
    atomss, statuses, products = test_case
    matrix = stoichiometric_matrix(atomss, statuses)
    products_calc = balance_stoichiometry(matrix, rounding=11)
    assert_close1d(products_calc, products)
    atol = sum(abs(v) for r in matrix for v in r)*1e-10
    check_reaction_balance(matrix, products_calc, atol=atol)

stoich_test_cases = [
    [[{'Hg': 1, 'O': 1}, {'Hg': 1}, {'O': 2}], [True, False, False], [2.0, 2.0, 1.0]],
    [[{'Cl': 2}, {'C': 3, 'H': 6}, {'C': 3, 'Cl': 1, 'H': 5}, {'Cl': 1, 'H': 1}],
      [True, True, False, False, False],
      [1, 1, 1, 1]],
    [[{'Al': 1}, {'H': 1, 'N': 1, 'O': 3}, {'Al': 1, 'N': 3, 'O': 9}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
      [True, True, False, False, False],
      [1.0, 4.0, 1.0, 1.0, 2.0]],
    [[{'Fe': 1}, {'O': 2}, {'Fe':2, 'O': 3}], [True, True, False], [4.0, 3.0, 2.0]],
    [[{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False], [4.0, 5.0, 4.0, 6.0]],
    [[{'O': 2}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}, {'C': 6, 'H': 14}], [True, False, False, True], [19.0, 14.0, 12.0, 2.0]],

    # Examples from
    #Smith, William R., and Ronald W. Missen. “Using Mathematica and Maple To Obtain Chemical Equations.” Journal of Chemical Education 74, no. 11 (November 1, 1997): 1369. https://doi.org/10.1021/ed074p1369.

    # Complex redox system test case
    # [Cr(N2H4CO)6]4[Cr(CN)6]3 + KMnO4 + H2SO4 → K2Cr2O7 + MnSO4 + CO2 + KNO3 + K2SO4 + H2O
    [[
        # B = [Cr(N2H4CO)6]4[Cr(CN)6]3
        {'Cr': 7, 'N': 66, 'H': 96, 'C': 42, 'O': 24},
        # KMnO4
        {'K': 1, 'Mn': 1, 'O': 4},
        # H2SO4
        {'H': 2, 'S': 1, 'O': 4},
        # K2Cr2O7
        {'K': 2, 'Cr': 2, 'O': 7},
        # MnSO4
        {'Mn': 1, 'S': 1, 'O': 4},
        # CO2
        {'C': 1, 'O': 2},
        # KNO3
        {'K': 1, 'N': 1, 'O': 3},
        # K2SO4
        {'K': 2, 'S': 1, 'O': 4},
        # H2O
        {'H': 2, 'O': 1}
    ],
    [True, True, True, False, False, False, False, False, False],
    [10.0, 1176.0, 1399.0, 35.0, 1176.0, 420.0, 660.0, 223.0, 1879.0]],  # Expected coefficients

    # Reaction 1: H+ + OH- = H2O
    [
        [
            {'H': 1, 'p': 1},          # H+
            {'H': 1, 'O': 1, 'p': -1}, # OH-
            {'H': 2, 'O': 1}           # H2O
        ],
        [True, True, False],
        [1.0, 1.0, 1.0]
    ],

    # Reaction 2: OH- + NO+ = NO2- + H+
    [
        [
            {'H': 1, 'O': 1, 'p': -1}, # OH-
            {'N': 1, 'O': 1, 'p': 1},  # NO+
            {'N': 1, 'O': 2, 'p': -1}, # NO2-
            {'H': 1, 'p': 1}           # H+
        ],
        [True, True, False, False],
        [1.0, 1.0, 1.0, 1.0]
    ],

    # Reaction 3: 2OH- + 2NO+ = N2O3 + H+
    [
        [
            {'H': 1, 'O': 1, 'p': -1}, # OH-
            {'N': 1, 'O': 1, 'p': 1},  # NO+
            {'N': 2, 'O': 3},          # N2O3
            {'H': 1, 'p': 1}           # H+
        ],
        [True, True, False, False],
        [1.0, 2.0, 1.0, 1.0]
    ],

    # Reaction 4: OH- + NO+ = HNO2
    [
        [
            {'H': 1, 'O': 1, 'p': -1}, # OH-
            {'N': 1, 'O': 1, 'p': 1},  # NO+
            {'H': 1, 'N': 1, 'O': 2}   # HNO2
        ],
        [True, True, False],
        [1.0, 1.0, 1.0]
    ],

    # Reaction 5: OH- + NO+ + Tl+ = TlNO2 + H+
    [
        [
            {'H': 1, 'O': 1, 'p': -1}, # OH-
            {'N': 1, 'O': 1, 'p': 1},  # NO+
            {'Tl': 1, 'p': 1},         # Tl+
            {'Tl': 1, 'N': 1, 'O': 2}, # TlNO2
            {'H': 1, 'p': 1}           # H+
        ],
        [True, True, True, False, False],
        [1.0, 1.0, 1.0, 1.0, 1.0]
    ],

    # Herndon, William C. “On Balancing Chemical Equations: Past and Present.” Journal of Chemical Education 74, no. 11 (November 1997): 1359. https://doi.org/10.1021/ed074p1359.
    [
        [
            {'H': 1, 'p': 1},             # H+
            {'p': -1},                    # e-
            {'Fe': 1, 'Cl': 3},           # FeCl3
            {'H': 1, 'I': 1, 'O': 3},     # HIO3
            {'Fe': 1, 'I': 2},            # FeI2
            {'H': 1, 'Cl': 1},            # HCl
            {'H': 2, 'O': 1}              # H2O
        ],
        [True, True, True, True, False, False, False],
        [13.0, 13.0, 1.0, 2.0, 1.0, 3.0, 6.0]
    ],

    # from http://myweb.astate.edu/mdraganj/BalanceEqn.html

    # C3H8 + O2 = CO2 + H2O
    [[{'C': 3, 'H': 8}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 5.0, 3.0, 4.0]],

    # Al2(SO3)3 + NaOH = Na2SO3 + Al(OH)3
    [[{'Al': 2, 'S': 3, 'O': 9}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'S': 1, 'O': 3}, {'Al': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [1.0, 6.0, 3.0, 2.0]],

    # Al2O3 + Fe = Fe3O4 + Al
    [[{'Al': 2, 'O': 3}, {'Fe': 1}, {'Fe': 3, 'O': 4}, {'Al': 1}],
    [True, True, False, False],
    [4.0, 9.0, 3.0, 8.0]],

    # KClO3 = KCl + O2
    [[{'K': 1, 'Cl': 1, 'O': 3}, {'K': 1, 'Cl': 1}, {'O': 2}],
    [True, False, False],
    [2.0, 2.0, 3.0]],

    # NH4NO3 = N2O + H2O
    [[{'N': 2, 'H': 4, 'O': 3}, {'N': 2, 'O': 1}, {'H': 2, 'O': 1}],
    [True, False, False],
    [1.0, 1.0, 2.0]],

    # NaHCO3 = Na2CO3 + H2O + CO2
    [[{'Na': 1, 'H': 1, 'C': 1, 'O': 3}, {'Na': 2, 'C': 1, 'O': 3}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    [True, False, False, False],
    [2.0, 1.0, 1.0, 1.0]],

    # P4O10 + H2O = H3PO4
    [[{'P': 4, 'O': 10}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}],
    [True, True, False],
    [1.0, 6.0, 4.0]],

    # Be2C + H2O = Be(OH)2 + CH4
    [[{'Be': 2, 'C': 1}, {'H': 2, 'O': 1}, {'Be': 1, 'O': 2, 'H': 2}, {'C': 1, 'H': 4}],
    [True, True, False, False],
    [1.0, 4.0, 2.0, 1.0]],

    # Al + H2SO4 = Al2(SO4)3 + H2
    [[{'Al': 1}, {'H': 2, 'S': 1, 'O': 4}, {'Al': 2, 'S': 3, 'O': 12}, {'H': 2}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 3.0]],

    # S + HNO3 = H2SO4 + NO2 + H2O
    [[{'S': 1}, {'H': 1, 'N': 1, 'O': 3}, {'H': 2, 'S': 1, 'O': 4}, {'N': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 6.0, 1.0, 6.0, 2.0]],

    # NH3 + CuO = Cu + N2 + H2O
    [[{'N': 1, 'H': 3}, {'Cu': 1, 'O': 1}, {'Cu': 1}, {'N': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [2.0, 3.0, 3.0, 1.0, 3.0]],

    # Cu + HNO3 = Cu(NO3)2 + NO + H2O
    [[{'Cu': 1}, {'H': 1, 'N': 1, 'O': 3}, {'Cu': 1, 'N': 2, 'O': 6}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [3.0, 8.0, 3.0, 2.0, 4.0]],

    # from https://www.chembuddy.com/balancing-stoichiometry-questions-balancing-questions

    # K4Fe(CN)6 + KMnO4 + H2SO4 = KHSO4 + Fe2(SO4)3 + MnSO4 + HNO3 + CO2 + H2O
    [[{'K': 4, 'Fe': 1, 'C': 6, 'N': 6}, {'K': 1, 'Mn': 1, 'O': 4}, {'H': 2, 'S': 1, 'O': 4}, {'K': 1, 'H': 1, 'S': 1, 'O': 4}, {'Fe': 2, 'S': 3, 'O': 12}, {'Mn': 1, 'S': 1, 'O': 4}, {'H': 1, 'N': 1, 'O': 3}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, True, False, False, False, False, False, False],
    [10.0, 122.0, 299.0, 162.0, 5.0, 122.0, 60.0, 60.0, 188.0]],

    # K3Fe(SCN)6 + Na2Cr2O7 + H2SO4 = Fe(NO3)3 + Cr2(SO4)3 + CO2 + H2O + Na2SO4 + KNO3
    [[{'K': 3, 'Fe': 1, 'S': 6, 'C': 6, 'N': 6}, {'Na': 2, 'Cr': 2, 'O': 7}, {'H': 2, 'S': 1, 'O': 4}, {'Fe': 1, 'N': 3, 'O': 9}, {'Cr': 2, 'S': 3, 'O': 12}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'Na': 2, 'S': 1, 'O': 4}, {'K': 1, 'N': 1, 'O': 3}],
    [True, True, True, False, False, False, False, False, False],
    [1.0, 16.0, 58.0, 1.0, 16.0, 6.0, 58.0, 16.0, 3.0]],

    # K4Fe(SCN)6 + K2Cr2O7 + H2SO4 = Fe2(SO4)3 + Cr2(SO4)3 + CO2 + H2O + K2SO4 + KNO3
    [[{'K': 4, 'Fe': 1, 'S': 6, 'C': 6, 'N': 6}, {'K': 2, 'Cr': 2, 'O': 7}, {'H': 2, 'S': 1, 'O': 4}, {'Fe': 2, 'S': 3, 'O': 12}, {'Cr': 2, 'S': 3, 'O': 12}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'K': 2, 'S': 1, 'O': 4}, {'K': 1, 'N': 1, 'O': 3}],
    [True, True, True, False, False, False, False, False, False],
    [6.0, 97.0, 355.0, 3.0, 97.0, 36.0, 355.0, 91.0, 36.0]],

    # Ca10F2(PO4)6 + H2SO4 = Ca(H2PO4)2 + CaSO4 + HF
    [[{'Ca': 10, 'F': 2, 'P': 6, 'O': 24}, {'H': 2, 'S': 1, 'O': 4}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}, {'Ca': 1, 'S': 1, 'O': 4}, {'H': 1, 'F': 1}],
    [True, True, False, False, False],
    [1.0, 7.0, 3.0, 7.0, 2.0]],

    # Ca5F(PO4)3 + H2SO4 = Ca(H2PO4)2 + CaSO4 + HF
    [[{'Ca': 5, 'F': 1, 'P': 3, 'O': 12}, {'H': 2, 'S': 1, 'O': 4}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}, {'Ca': 1, 'S': 1, 'O': 4}, {'H': 1, 'F': 1}],
    [True, True, False, False, False],
    [2.0, 7.0, 3.0, 7.0, 2.0]],

    # UO2(NO3)2.6H2O = UO3 + NO2 + O2 + H2O
    [[{'U': 1, 'O': 10.8, 'N': 2.6, 'H': 2}, {'U': 1, 'O': 3}, {'N': 1, 'O': 2}, {'O': 2}, {'H': 2, 'O': 1}],
    [True, False, False, False, False],
    [5.0, 5.0, 13.0, 4.0, 5.0]],

    # S8 + O2 = SO3
    [[{'S': 8}, {'O': 2}, {'S': 1, 'O': 3}],
    [True, True, False],
    [1.0, 12.0, 8.0]],

    # NH3 + NO = N2 + H2O
    [[{'N': 1, 'H': 3}, {'N': 1, 'O': 1}, {'N': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [4.0, 6.0, 5.0, 6.0]],

    # HClO4 + P4O10 = H3PO4 + Cl2O7
    [[{'H': 1, 'Cl': 1, 'O': 4}, {'P': 4, 'O': 10}, {'H': 3, 'P': 1, 'O': 4}, {'Cl': 2, 'O': 7}],
    [True, True, False, False],
    [12.0, 1.0, 4.0, 6.0]],

    # Au + KCN + O2 + H2O = K[Au(CN)2] + KOH
    [[{'Au': 1}, {'K': 1, 'C': 1, 'N': 1}, {'O': 2}, {'H': 2, 'O': 1}, {'K': 1, 'Au': 1, 'C': 2, 'N': 2}, {'K': 1, 'O': 1, 'H': 1}],
    [True, True, True, True, False, False],
    [4.0, 8.0, 1.0, 2.0, 4.0, 4.0]],

    # CO2 + H2O = C6H12O6 + O2
    [[{'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'C': 6, 'H': 12, 'O': 6}, {'O': 2}],
    [True, True, False, False],
    [6.0, 6.0, 1.0, 6.0]],

    # V2O5 + Al = Al2O3 + V
    [[{'V': 2, 'O': 5}, {'Al': 1}, {'Al': 2, 'O': 3}, {'V': 1}],
    [True, True, False, False],
    [3.0, 10.0, 5.0, 6.0]],

    # FeS2 + O2 = Fe2O3 + SO2
    [[{'Fe': 1, 'S': 2}, {'O': 2}, {'Fe': 2, 'O': 3}, {'S': 1, 'O': 2}],
    [True, True, False, False],
    [4.0, 11.0, 2.0, 8.0]],

    # Si2H3 + O2 = SiO2 + H2O
    [[{'Si': 2, 'H': 3}, {'O': 2}, {'Si': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [4.0, 11.0, 8.0, 6.0]],

    # P4 + H2O = H3PO4 + H2
    [[{'P': 4}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}, {'H': 2}],
    [True, True, False, False],
    [1.0, 16.0, 4.0, 10.0]],

    # H2S + Cl2 = S8 + HCl
    [[{'H': 2, 'S': 1}, {'Cl': 2}, {'S': 8}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [8.0, 8.0, 1.0, 16.0]],

    # C4H10 + O2 = CO2 + H2O
    [[{'C': 4, 'H': 10}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 13.0, 8.0, 10.0]],

    # C6H6 + O2 = CO2 + H2O
    [[{'C': 6, 'H': 6}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 15.0, 12.0, 6.0]],

    # C7H16 + O2 = CO2 + H2O
    [[{'C': 7, 'H': 16}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 11.0, 7.0, 8.0]],

    # Ca3(PO4)2 + SiO2 + C = CaSiO3 + P4 + CO
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'Si': 1, 'O': 2}, {'C': 1}, {'Ca': 1, 'Si': 1, 'O': 3}, {'P': 4}, {'C': 1, 'O': 1}],
    [True, True, True, False, False, False],
    [2.0, 6.0, 10.0, 6.0, 1.0, 10.0]],

    # C10H16 + Cl2 = C + HCl
    [[{'C': 10, 'H': 16}, {'Cl': 2}, {'C': 1}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 8.0, 10.0, 16.0]],

    # C7H6O2 + O2 = CO2 + H2O
    [[{'C': 7, 'H': 6, 'O': 2}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 15.0, 14.0, 6.0]],

    # C7H10N + O2 = CO2 + H2O + NO2
    [[{'C': 7, 'H': 10, 'N': 1}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 1, 'O': 2}],
    [True, True, False, False, False],
    [2.0, 21.0, 14.0, 10.0, 2.0]],

    # KNO3 + C12H22O11 = N2 + CO2 + H2O + K2CO3
    [[{'K': 1, 'N': 1, 'O': 3}, {'C': 12, 'H': 22, 'O': 11}, {'N': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'K': 2, 'C': 1, 'O': 3}],
    [True, True, False, False, False, False],
    [48.0, 5.0, 24.0, 36.0, 55.0, 24.0]],

    # V2O5 + HCl = VOCl3 + H2O
    [[{'V': 2, 'O': 5}, {'H': 1, 'Cl': 1}, {'V': 1, 'O': 1, 'Cl': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 6.0, 2.0, 3.0]],

    # Cr(OH)3 + H2SO4 = Cr2(SO4)3 + H2O
    [[{'Cr': 1, 'O': 3, 'H': 3}, {'H': 2, 'S': 1, 'O': 4}, {'Cr': 2, 'S': 3, 'O': 12}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 6.0]],

    # Hg(OH)2 + H3PO4 = Hg3(PO4)2 + H2O
    [[{'Hg': 1, 'O': 2, 'H': 2}, {'H': 3, 'P': 1, 'O': 4}, {'Hg': 3, 'P': 2, 'O': 8}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [3.0, 2.0, 1.0, 6.0]],

    # Fe + H2O = Fe3O4 + H2
    [[{'Fe': 1}, {'H': 2, 'O': 1}, {'Fe': 3, 'O': 4}, {'H': 2}],
    [True, True, False, False],
    [3.0, 4.0, 1.0, 4.0]],

    # Ca3P2 + H2O = Ca(OH)2 + PH3
    [[{'Ca': 3, 'P': 2}, {'H': 2, 'O': 1}, {'Ca': 1, 'O': 2, 'H': 2}, {'P': 1, 'H': 3}],
    [True, True, False, False],
    [1.0, 6.0, 3.0, 2.0]],

    # H2SO4 + Al(OH)3 = Al2(SO4)3 + H2O
    [[{'H': 2, 'S': 1, 'O': 4}, {'Al': 1, 'O': 3, 'H': 3}, {'Al': 2, 'S': 3, 'O': 12}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [3.0, 2.0, 1.0, 6.0]],

    # Al(NO3)3 + Na2CO3 = Al2(CO3)3 + NaNO3
    [[{'Al': 1, 'N': 3, 'O': 9}, {'Na': 2, 'C': 1, 'O': 3}, {'Al': 2, 'C': 3, 'O': 9}, {'Na': 1, 'N': 1, 'O': 3}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 6.0]],

    # K2MnO4 + H2SO4 = KMnO4 + MnO2 + K2SO4 + H2O
    [[{'K': 2, 'Mn': 1, 'O': 4}, {'H': 2, 'S': 1, 'O': 4}, {'K': 1, 'Mn': 1, 'O': 4}, {'Mn': 1, 'O': 2}, {'K': 2, 'S': 1, 'O': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False, False, False],
    [3.0, 2.0, 2.0, 1.0, 2.0, 2.0]],

    # Na3AsO3 + H2S = As2S3 + NaOH
    [[{'Na': 3, 'As': 1, 'O': 3}, {'H': 2, 'S': 1}, {'As': 2, 'S': 3}, {'Na': 1, 'O': 1, 'H': 1}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 6.0]],

    # Mg3N2 + H2O = Mg(OH)2 + NH3
    [[{'Mg': 3, 'N': 2}, {'H': 2, 'O': 1}, {'Mg': 1, 'O': 2, 'H': 2}, {'N': 1, 'H': 3}],
    [True, True, False, False],
    [1.0, 6.0, 3.0, 2.0]],

    # Fe3O4 + H2 = Fe + H2O
    [[{'Fe': 3, 'O': 4}, {'H': 2}, {'Fe': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 4.0, 3.0, 4.0]],

    # C2H2 + O2 = CO2 + H2O
    [[{'C': 2, 'H': 2}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 5.0, 4.0, 2.0]],

    # (NH4)2Cr2O7 = NH3 + H2O + Cr2O3 + O2
    [[{'N': 2, 'H': 8, 'Cr': 2, 'O': 7}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}, {'Cr': 2, 'O': 3}, {'O': 2}],
    [True, False, False, False, False],
    [2.0, 4.0, 2.0, 2.0, 3.0]],

    # C3H8 + O2 = CO2 + H2O
    [[{'C': 3, 'H': 8}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 5.0, 3.0, 4.0]],

    # As + NaOH = Na3AsO3 + H2
    [[{'As': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 3, 'As': 1, 'O': 3}, {'H': 2}],
    [True, True, False, False],
    [2.0, 6.0, 2.0, 3.0]],

    # H3BO3 + Na2CO3 = Na2B4O7 + CO2 + H2O
    [[{'H': 3, 'B': 1, 'O': 3}, {'Na': 2, 'C': 1, 'O': 3}, {'Na': 2, 'B': 4, 'O': 7}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [4.0, 1.0, 1.0, 1.0, 6.0]],

    # Al + HCl = AlCl3 + H2
    [[{'Al': 1}, {'H': 1, 'Cl': 1}, {'Al': 1, 'Cl': 3}, {'H': 2}],
    [True, True, False, False],
    [2.0, 6.0, 2.0, 3.0]],

    # V2O5 + Ca = CaO + V
    [[{'V': 2, 'O': 5}, {'Ca': 1}, {'Ca': 1, 'O': 1}, {'V': 1}],
    [True, True, False, False],
    [1.0, 5.0, 5.0, 2.0]],

    # H3BO3 = H4B6O11 + H2O
    [[{'H': 3, 'B': 1, 'O': 3}, {'H': 4, 'B': 6, 'O': 11}, {'H': 2, 'O': 1}],
    [True, False, False],
    [6.0, 1.0, 7.0]],

    # Na2B4O7 + HCl + H2O = NaCl + H3BO3
    [[{'Na': 2, 'B': 4, 'O': 7}, {'H': 1, 'Cl': 1}, {'H': 2, 'O': 1}, {'Na': 1, 'Cl': 1}, {'H': 3, 'B': 1, 'O': 3}],
    [True, True, True, False, False],
    [1.0, 2.0, 5.0, 2.0, 4.0]],

    # Pb + Na + C2H5Cl = Pb(C2H5)4 + NaCl
    [[{'Pb': 1}, {'Na': 1}, {'C': 2, 'H': 5, 'Cl': 1}, {'Pb': 1, 'C': 8, 'H': 20}, {'Na': 1, 'Cl': 1}],
    [True, True, True, False, False],
    [1.0, 4.0, 4.0, 1.0, 4.0]],

    # C2H3Cl + O2 = CO2 + H2O + HCl
    [[{'C': 2, 'H': 3, 'Cl': 1}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}, {'H': 1, 'Cl': 1}],
    [True, True, False, False, False],
    [2.0, 5.0, 4.0, 2.0, 2.0]],

    # Ca3(PO4)2 + SiO2 = P4O10 + CaSiO3
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'Si': 1, 'O': 2}, {'P': 4, 'O': 10}, {'Ca': 1, 'Si': 1, 'O': 3}],
    [True, True, False, False],
    [2.0, 6.0, 1.0, 6.0]],

    # Se + NaOH = Na2Se + Na2SeO3 + H2O
    [[{'Se': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'Se': 1}, {'Na': 2, 'Se': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [3.0, 6.0, 2.0, 1.0, 3.0]],

    # Al + NaOH + H2O = NaAl(OH)4 + H2
    [[{'Al': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}, {'Na': 1, 'Al': 1, 'O': 4, 'H': 4}, {'H': 2}],
    [True, True, True, False, False],
    [2.0, 2.0, 6.0, 2.0, 3.0]],

    # K3AsO4 + H2S = As2S5 + KOH + H2O
    [[{'K': 3, 'As': 1, 'O': 4}, {'H': 2, 'S': 1}, {'As': 2, 'S': 5}, {'K': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [2.0, 5.0, 1.0, 6.0, 2.0]],

    # I2 + HNO3 = HIO3 + NO2 + H2
    [[{'I': 2}, {'H': 1, 'N': 1, 'O': 3}, {'H': 1, 'I': 1, 'O': 3}, {'N': 1, 'O': 2}, {'H': 2}],
    [True, True, False, False, False],
    [1.0, 6.0, 2.0, 6.0, 2.0]],

    # Al + NH4ClO4 = Al2O3 + AlCl3 + NO + H2O
    [[{'Al': 1}, {'N': 1, 'H': 4, 'Cl': 1, 'O': 4}, {'Al': 2, 'O': 3}, {'Al': 1, 'Cl': 3}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False, False],
    [3.0, 3.0, 1.0, 1.0, 3.0, 6.0]],

    # FeS + O2 = Fe2O3 + SO2
    [[{'Fe': 1, 'S': 1}, {'O': 2}, {'Fe': 2, 'O': 3}, {'S': 1, 'O': 2}],
    [True, True, False, False],
    [4.0, 7.0, 2.0, 4.0]],

    # Ca3(PO4)2 + C = Ca3P2 + CO
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'C': 1}, {'Ca': 3, 'P': 2}, {'C': 1, 'O': 1}],
    [True, True, False, False],
    [1.0, 8.0, 1.0, 8.0]],

    # H2SO4 + HI = H2S + I2 + H2O
    [[{'H': 2, 'S': 1, 'O': 4}, {'H': 1, 'I': 1}, {'H': 2, 'S': 1}, {'I': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 8.0, 1.0, 4.0, 4.0]],

    # U3O8 + HNO3 = UO2(NO3)2 + NO2 + H2O
    [[{'U': 3, 'O': 8}, {'H': 1, 'N': 1, 'O': 3}, {'U': 1, 'O': 8, 'N': 2}, {'N': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 8.0, 3.0, 2.0, 4.0]],

    # (NH4)3AsS4 + HCl = As2S5 + H2S + NH4Cl
    [[{'N': 3, 'H': 12, 'As': 1, 'S': 4}, {'H': 1, 'Cl': 1}, {'As': 2, 'S': 5}, {'H': 2, 'S': 1}, {'N': 1, 'H': 4, 'Cl': 1}],
    [True, True, False, False, False],
    [2.0, 6.0, 1.0, 3.0, 6.0]],

    # Pb3(VO4)2.PbCl2 + HCl = VO2Cl + PbCl2 + H2O
    [[{'Pb': 4, 'V': 2, 'O': 8, 'Cl': 2}, {'H': 1, 'Cl': 1}, {'V': 1, 'O': 2, 'Cl': 1}, {'Pb': 1, 'Cl': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 8.0, 2.0, 4.0, 4.0]],

    # NH3 + O2 = NO + H2O
    [[{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [4.0, 5.0, 4.0, 6.0]],

    # Hg2CrO4 = Cr2O3 + Hg + O2
    [[{'Hg': 2, 'Cr': 1, 'O': 4}, {'Cr': 2, 'O': 3}, {'Hg': 1}, {'O': 2}],
    [True, False, False, False],
    [4.0, 2.0, 8.0, 5.0]],

    # Al4C3 + H2O = CH4 + Al(OH)3
    [[{'Al': 4, 'C': 3}, {'H': 2, 'O': 1}, {'C': 1, 'H': 4}, {'Al': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [1.0, 12.0, 3.0, 4.0]],

    # CaHPO4(H2O)2 + NaOH + H2O = Na2HPO4(H2O)12 + Ca(OH)2
    [[{'Ca': 1, 'H': 5, 'P': 1, 'O': 6}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}, {'Na': 2, 'H': 25, 'P': 1, 'O': 16}, {'Ca': 1, 'O': 2, 'H': 2}],
    [True, True, True, False, False],
    [1.0, 2.0, 10.0, 1.0, 1.0]],

    # FeC2O4(H2O)2 + H2C2O4 + H2O2 + K2C2O4 = K3Fe(C2O4)3(H2O)3
    [[{'Fe': 1, 'C': 2, 'O': 6, 'H': 4}, {'H': 2, 'C': 2, 'O': 4}, {'H': 2, 'O': 2}, {'K': 2, 'C': 2, 'O': 4}, {'K': 3, 'Fe': 1, 'C': 6, 'O': 15, 'H': 6}],
    [True, True, True, True, False],
    [2.0, 1.0, 1.0, 3.0, 2.0]],

    # MgNH4AsO4(H2O)6 = Mg2As2O7 + NH3 + H2O
    [[{'Mg': 1, 'N': 1, 'H': 16, 'As': 1, 'O': 10}, {'Mg': 2, 'As': 2, 'O': 7}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [2.0, 1.0, 2.0, 13.0]],


    # POCl3 + H2O = H3PO4 + HCl
    [[{'P': 1, 'O': 1, 'Cl': 3}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 3.0]],

    # C2H5OH + O2 = CO + H2O
    [[{'C': 2, 'H': 6, 'O': 1}, {'O': 2}, {'C': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 2.0, 3.0]],

    # UO2 + HF = UF4 + H2O
    [[{'U': 1, 'O': 2}, {'H': 1, 'F': 1}, {'U': 1, 'F': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 4.0, 1.0, 2.0]],

    # Ag2S + KCN = KAg(CN)2 + K2S
    [[{'Ag': 2, 'S': 1}, {'K': 1, 'C': 1, 'N': 1}, {'K': 1, 'Ag': 1, 'C': 2, 'N': 2}, {'K': 2, 'S': 1}],
    [True, True, False, False],
    [1.0, 4.0, 2.0, 1.0]],

    # C + SiO2 + Cl2 = SiCl4 + CO
    [[{'C': 1}, {'Si': 1, 'O': 2}, {'Cl': 2}, {'Si': 1, 'Cl': 4}, {'C': 1, 'O': 1}],
    [True, True, True, False, False],
    [2.0, 1.0, 2.0, 1.0, 2.0]],

    # PCl3 + H2O = H3PO3 + HCl
    [[{'P': 1, 'Cl': 3}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 3}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 3.0]],

    # H3PO4 + (NH4)2MoO4 + HNO3 = (NH4)3PO4(MoO3)12 + NH4NO3 + H2O
    [[{'H': 3, 'P': 1, 'O': 4}, {'N': 2, 'H': 8, 'Mo': 1, 'O': 4}, {'H': 1, 'N': 1, 'O': 3}, {'N': 3, 'H': 12, 'P': 1, 'O': 40, 'Mo': 12}, {'N': 2, 'H': 4, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, True, False, False, False],
    [1.0, 12.0, 21.0, 1.0, 21.0, 12.0]],

    # MnO2 + HCl = MnCl2 + H2O + Cl2
    [[{'Mn': 1, 'O': 2}, {'H': 1, 'Cl': 1}, {'Mn': 1, 'Cl': 2}, {'H': 2, 'O': 1}, {'Cl': 2}],
    [True, True, False, False, False],
    [1.0, 4.0, 1.0, 2.0, 1.0]],

    # Fe2O3 + C = CO + Fe
    [[{'Fe': 2, 'O': 3}, {'C': 1}, {'C': 1, 'O': 1}, {'Fe': 1}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 2.0]],

    # PCl5 + P2O5 = POCl3
    [[{'P': 1, 'Cl': 5}, {'P': 2, 'O': 5}, {'P': 1, 'O': 1, 'Cl': 3}],
    [True, True, False],
    [3.0, 1.0, 5.0]],

    # FeO + H3PO4 = Fe3(PO4)2 + H2O
    [[{'Fe': 1, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}, {'Fe': 3, 'P': 2, 'O': 8}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [3.0, 2.0, 1.0, 3.0]],

    # Ca(NO3)2 = CaO + NO2 + O2
    [[{'Ca': 1, 'N': 2, 'O': 6}, {'Ca': 1, 'O': 1}, {'N': 1, 'O': 2}, {'O': 2}],
    [True, False, False, False],
    [2.0, 2.0, 4.0, 1.0]],

    # Fe + O2 = Fe2O3
    [[{'Fe': 1}, {'O': 2}, {'Fe': 2, 'O': 3}],
    [True, True, False],
    [4.0, 3.0, 2.0]],

    # Fe2O3 + H2 = Fe + H2O
    [[{'Fe': 2, 'O': 3}, {'H': 2}, {'Fe': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 3.0, 2.0, 3.0]],

    # FeSO4 + K3Fe(CN)6 = Fe3(Fe(CN)6)2 + K2SO4
    [[{'Fe': 1, 'S': 1, 'O': 4}, {'K': 3, 'Fe': 1, 'C': 6, 'N': 6}, {'Fe': 5, 'C': 12, 'N': 12}, {'K': 2, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [3.0, 2.0, 1.0, 3.0]],

    # NH3 + O2 = HNO2 + H2O
    [[{'N': 1, 'H': 3}, {'O': 2}, {'H': 1, 'N': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 3.0, 2.0, 2.0]],

    # Al + O2 = Al2O3
    [[{'Al': 1}, {'O': 2}, {'Al': 2, 'O': 3}],
    [True, True, False],
    [4.0, 3.0, 2.0]],

    # BaCl2 + Al2(SO4)3 = BaSO4 + AlCl3
    [[{'Ba': 1, 'Cl': 2}, {'Al': 2, 'S': 3, 'O': 12}, {'Ba': 1, 'S': 1, 'O': 4}, {'Al': 1, 'Cl': 3}],
    [True, True, False, False],
    [3.0, 1.0, 3.0, 2.0]],

    # Fe2(SO4)3 + Ba(NO3)2 = BaSO4 + Fe(NO3)3
    [[{'Fe': 2, 'S': 3, 'O': 12}, {'Ba': 1, 'N': 2, 'O': 6}, {'Ba': 1, 'S': 1, 'O': 4}, {'Fe': 1, 'N': 3, 'O': 9}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 2.0]],

    # Au2S3 + H2 = Au + H2S
    [[{'Au': 2, 'S': 3}, {'H': 2}, {'Au': 1}, {'H': 2, 'S': 1}],
    [True, True, False, False],
    [1.0, 3.0, 2.0, 3.0]],

    # Au + HCl + HNO3 = AuCl3 + NO + H2O
    [[{'Au': 1}, {'H': 1, 'Cl': 1}, {'H': 1, 'N': 1, 'O': 3}, {'Au': 1, 'Cl': 3}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, True, False, False, False],
    [1.0, 3.0, 1.0, 1.0, 1.0, 2.0]],

    # NiS + O2 = NiO + SO2
    [[{'Ni': 1, 'S': 1}, {'O': 2}, {'Ni': 1, 'O': 1}, {'S': 1, 'O': 2}],
    [True, True, False, False],
    [2.0, 3.0, 2.0, 2.0]],

    # Al + FeO = Al2O3 + Fe
    [[{'Al': 1}, {'Fe': 1, 'O': 1}, {'Al': 2, 'O': 3}, {'Fe': 1}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 3.0]],

    # C2H5OH + O2 = CO2 + H2O
    [[{'C': 2, 'H': 6, 'O': 1}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 3.0, 2.0, 3.0]],

    # Na2O2 + H2O = NaOH + O2
    [[{'Na': 2, 'O': 2}, {'H': 2, 'O': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'O': 2}],
    [True, True, False, False],
    [2.0, 2.0, 4.0, 1.0]],

    # Fe2O3 + CO = Fe + CO2
    [[{'Fe': 2, 'O': 3}, {'C': 1, 'O': 1}, {'Fe': 1}, {'C': 1, 'O': 2}],
    [True, True, False, False],
    [1.0, 3.0, 2.0, 3.0]],

    # Pb(NO3)2 = PbO + NO2 + O2
    [[{'Pb': 1, 'N': 2, 'O': 6}, {'Pb': 1, 'O': 1}, {'N': 1, 'O': 2}, {'O': 2}],
    [True, False, False, False],
    [2.0, 2.0, 4.0, 1.0]],

    # Al2(SO4)3 + Ca(OH)2 = CaSO4 + Al(OH)3
    [[{'Al': 2, 'S': 3, 'O': 12}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'S': 1, 'O': 4}, {'Al': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 2.0]],

    # Au2O3 = Au + O2
    [[{'Au': 2, 'O': 3}, {'Au': 1}, {'O': 2}],
    [True, False, False],
    [2.0, 4.0, 3.0]],

    # Ca3(PO4)2 + H2SO4 = CaSO4 + H3PO4
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'H': 2, 'S': 1, 'O': 4}, {'Ca': 1, 'S': 1, 'O': 4}, {'H': 3, 'P': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 2.0]],

    # SiCl4 + H2O = H4SiO4 + HCl
    [[{'Si': 1, 'Cl': 4}, {'H': 2, 'O': 1}, {'H': 4, 'Si': 1, 'O': 4}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 4.0, 1.0, 4.0]],

    # Ca + AlCl3 = CaCl2 + Al
    [[{'Ca': 1}, {'Al': 1, 'Cl': 3}, {'Ca': 1, 'Cl': 2}, {'Al': 1}],
    [True, True, False, False],
    [3.0, 2.0, 3.0, 2.0]],

    # FeCl3 + Ca(OH)2 = CaCl2 + Fe(OH)3
    [[{'Fe': 1, 'Cl': 3}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'Cl': 2}, {'Fe': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [2.0, 3.0, 3.0, 2.0]],

    # Al2O3 + C + N2 = AlN + CO
    [[{'Al': 2, 'O': 3}, {'C': 1}, {'N': 2}, {'Al': 1, 'N': 1}, {'C': 1, 'O': 1}],
    [True, True, True, False, False],
    [1.0, 3.0, 1.0, 2.0, 3.0]],

    # NO + NaOH = NaNO2 + H2O + N2O
    [[{'N': 1, 'O': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 1, 'N': 1, 'O': 2}, {'H': 2, 'O': 1}, {'N': 2, 'O': 1}],
    [True, True, False, False, False],
    [4.0, 2.0, 2.0, 1.0, 1.0]],

    # Pb3O4 + HNO3 = Pb(NO3)2 + PbO2 + H2O
    [[{'Pb': 3, 'O': 4}, {'H': 1, 'N': 1, 'O': 3}, {'Pb': 1, 'N': 2, 'O': 6}, {'Pb': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 4.0, 2.0, 1.0, 2.0]],

    # CuSO4 + KCN = CuCN + K2SO4 + C2N2
    [[{'Cu': 1, 'S': 1, 'O': 4}, {'K': 1, 'C': 1, 'N': 1}, {'Cu': 1, 'C': 1, 'N': 1}, {'K': 2, 'S': 1, 'O': 4}, {'C': 2, 'N': 2}],
    [True, True, False, False, False],
    [2.0, 4.0, 2.0, 2.0, 1.0]],

    # KO2 + CO2 = K2CO3 + O2
    [[{'K': 1, 'O': 2}, {'C': 1, 'O': 2}, {'K': 2, 'C': 1, 'O': 3}, {'O': 2}],
    [True, True, False, False],
    [4.0, 2.0, 2.0, 3.0]],

    # P4O6 = P4 + P2O4
    [[{'P': 4, 'O': 6}, {'P': 4}, {'P': 2, 'O': 4}],
    [True, False, False],
    [4.0, 1.0, 6.0]],

    # P4O10 + H2O = H3PO4
    [[{'P': 4, 'O': 10}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}],
    [True, True, False],
    [1.0, 6.0, 4.0]],

    # Al + KOH + H2O = KAlO2 + H2
    [[{'Al': 1}, {'K': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}, {'K': 1, 'Al': 1, 'O': 2}, {'H': 2}],
    [True, True, True, False, False],
    [2.0, 2.0, 2.0, 2.0, 3.0]],

    # Fe + H2O + O2 = Fe2O3(H2O)
    [[{'Fe': 1}, {'H': 2, 'O': 1}, {'O': 2}, {'Fe': 2, 'O': 4, 'H': 2}],
    [True, True, True, False],
    [4.0, 2.0, 3.0, 2.0]],

    # H3PO4 + HCl = PCl5 + H2O
    [[{'H': 3, 'P': 1, 'O': 4}, {'H': 1, 'Cl': 1}, {'P': 1, 'Cl': 5}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 5.0, 1.0, 4.0]],

    # MnO2 + KOH + O2 = K2MnO4 + H2O
    [[{'Mn': 1, 'O': 2}, {'K': 1, 'O': 1, 'H': 1}, {'O': 2}, {'K': 2, 'Mn': 1, 'O': 4}, {'H': 2, 'O': 1}],
    [True, True, True, False, False],
    [2.0, 4.0, 1.0, 2.0, 2.0]],

    # K2CO3 + C + N2 = KCN + CO
    [[{'K': 2, 'C': 1, 'O': 3}, {'C': 1}, {'N': 2}, {'K': 1, 'C': 1, 'N': 1}, {'C': 1, 'O': 1}],
    [True, True, True, False, False],
    [1.0, 4.0, 1.0, 2.0, 3.0]],

    # PCl5 + H2O = H3PO4 + HCl
    [[{'P': 1, 'Cl': 5}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 4}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 4.0, 1.0, 5.0]],

    # P4O6 + H2O = H3PO3
    [[{'P': 4, 'O': 6}, {'H': 2, 'O': 1}, {'H': 3, 'P': 1, 'O': 3}],
    [True, True, False],
    [1.0, 6.0, 4.0]],

    # Al(OH)3 + H2SO4 = Al2(SO4)3 + H2O
    [[{'Al': 1, 'O': 3, 'H': 3}, {'H': 2, 'S': 1, 'O': 4}, {'Al': 2, 'S': 3, 'O': 12}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 6.0]],

    # Fe2(SO4)3 + KOH = K2SO4 + Fe(OH)3
    [[{'Fe': 2, 'S': 3, 'O': 12}, {'K': 1, 'O': 1, 'H': 1}, {'K': 2, 'S': 1, 'O': 4}, {'Fe': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [1.0, 6.0, 3.0, 2.0]],

    # Bi(NO3)3 + H2S = Bi2S3 + HNO3
    [[{'Bi': 1, 'N': 3, 'O': 9}, {'H': 2, 'S': 1}, {'Bi': 2, 'S': 3}, {'H': 1, 'N': 1, 'O': 3}],
    [True, True, False, False],
    [2.0, 3.0, 1.0, 6.0]],

    # MgNH4PO4 = Mg2P2O7 + NH3 + H2O
    [[{'Mg': 1, 'N': 1, 'H': 4, 'P': 1, 'O': 4}, {'Mg': 2, 'P': 2, 'O': 7}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [2.0, 1.0, 2.0, 1.0]],

    # H3PO4 + Ca(OH)2 = Ca(H2PO4)2 + H2O
    [[{'H': 3, 'P': 1, 'O': 4}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [2.0, 1.0, 1.0, 2.0]],

    # Ag(NH3)2Cl + HNO3 = NH4NO3 + AgCl
    [[{'Ag': 1, 'N': 2, 'H': 6, 'Cl': 1}, {'H': 1, 'N': 1, 'O': 3}, {'N': 2, 'H': 4, 'O': 3}, {'Ag': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 2.0, 2.0, 1.0]],

    # CaS + H2O = Ca(HS)2 + Ca(OH)2
    [[{'Ca': 1, 'S': 1}, {'H': 2, 'O': 1}, {'Ca': 1, 'H': 2, 'S': 2}, {'Ca': 1, 'O': 2, 'H': 2}],
    [True, True, False, False],
    [2.0, 2.0, 1.0, 1.0]],

    # Cu + CO2 + O2 + H2O = CuCO3Cu(OH)2
    [[{'Cu': 1}, {'C': 1, 'O': 2}, {'O': 2}, {'H': 2, 'O': 1}, {'Cu': 2, 'C': 1, 'O': 5, 'H': 2}],
    [True, True, True, True, False],
    [2.0, 1.0, 1.0, 1.0, 1.0]],

    # (NH4)2BeF4 = BeF2 + NH3 + HF
    [[{'N': 2, 'H': 8, 'Be': 1, 'F': 4}, {'Be': 1, 'F': 2}, {'N': 1, 'H': 3}, {'H': 1, 'F': 1}],
    [True, False, False, False],
    [1.0, 1.0, 2.0, 2.0]],

    # Sn(OH)2 + NaOH = Na2SnO2 + H2O
    [[{'Sn': 1, 'O': 2, 'H': 2}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'Sn': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # NH4VO3 = V2O5 + NH3 + H2O
    [[{'N': 1, 'H': 4, 'V': 1, 'O': 3}, {'V': 2, 'O': 5}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [2.0, 1.0, 2.0, 1.0]],

    # H3AsO3 = As2O3 + H2O
    [[{'H': 3, 'As': 1, 'O': 3}, {'As': 2, 'O': 3}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 3.0]],

    # NaCl + H2SO4 = Na2SO4 + HCl
    [[{'Na': 1, 'Cl': 1}, {'H': 2, 'S': 1, 'O': 4}, {'Na': 2, 'S': 1, 'O': 4}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [2.0, 1.0, 1.0, 2.0]],

    # Fe(OH)3 = Fe2O3 + H2O
    [[{'Fe': 1, 'O': 3, 'H': 3}, {'Fe': 2, 'O': 3}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 3.0]],

    # As2O5 + H2O = H3AsO4
    [[{'As': 2, 'O': 5}, {'H': 2, 'O': 1}, {'H': 3, 'As': 1, 'O': 4}],
    [True, True, False],
    [1.0, 3.0, 2.0]],

    # NaOH + Cl2 = NaCl + NaClO + H2O
    [[{'Na': 1, 'O': 1, 'H': 1}, {'Cl': 2}, {'Na': 1, 'Cl': 1}, {'Na': 1, 'Cl': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [2.0, 1.0, 1.0, 1.0, 1.0]],

    # VO2Cl + NH4OH = NH4VO3 + NH4Cl + H2O
    [[{'V': 1, 'O': 2, 'Cl': 1}, {'N': 1, 'H': 5, 'O': 1}, {'N': 1, 'H': 4, 'V': 1, 'O': 3}, {'N': 1, 'H': 4, 'Cl': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 1.0, 1.0]],

    # B2O3 + H2O = H3BO3
    [[{'B': 2, 'O': 3}, {'H': 2, 'O': 1}, {'H': 3, 'B': 1, 'O': 3}],
    [True, True, False],
    [1.0, 3.0, 2.0]],

    # CH4 + O2 = CO2 + H2O
    [[{'C': 1, 'H': 4}, {'O': 2}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # SiH4 + O2 = SiO2 + H2O
    [[{'Si': 1, 'H': 4}, {'O': 2}, {'Si': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # TiCl4 + Mg = MgCl2 + Ti
    [[{'Ti': 1, 'Cl': 4}, {'Mg': 1}, {'Mg': 1, 'Cl': 2}, {'Ti': 1}],
    [True, True, False, False],
    [1.0, 2.0, 2.0, 1.0]],

    # Pb(OH)2 + NaOH = Na2PbO2 + H2O
    [[{'Pb': 1, 'O': 2, 'H': 2}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'Pb': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # Si + NaOH + H2O = Na2SiO3 + H2
    [[{'Si': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}, {'Na': 2, 'Si': 1, 'O': 3}, {'H': 2}],
    [True, True, True, False, False],
    [1.0, 2.0, 1.0, 1.0, 2.0]],

    # Si + S8 = Si2S4
    [[{'Si': 1}, {'S': 8}, {'Si': 2, 'S': 4}],
    [True, True, False],
    [4.0, 1.0, 2.0]],

    # CaS2 + O2 = CaS2O3
    [[{'Ca': 1, 'S': 2}, {'O': 2}, {'Ca': 1, 'S': 2, 'O': 3}],
    [True, True, False],
    [2.0, 3.0, 2.0]],

    # Na2SnO3 + H2S = SnS2 + NaOH + H2O
    [[{'Na': 2, 'Sn': 1, 'O': 3}, {'H': 2, 'S': 1}, {'Sn': 1, 'S': 2}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 2.0, 1.0]],

    # Na2S2 + O2 = Na2S2O3
    [[{'Na': 2, 'S': 2}, {'O': 2}, {'Na': 2, 'S': 2, 'O': 3}],
    [True, True, False],
    [2.0, 3.0, 2.0]],

    # (NH4)2Cr2O7 = Cr2O3 + N2 + H2O
    [[{'N': 2, 'H': 8, 'Cr': 2, 'O': 7}, {'Cr': 2, 'O': 3}, {'N': 2}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [1.0, 1.0, 1.0, 4.0]],

    # HCl + K2CO3 = KCl + H2O + CO2
    [[{'H': 1, 'Cl': 1}, {'K': 2, 'C': 1, 'O': 3}, {'K': 1, 'Cl': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    [True, True, False, False, False],
    [2.0, 1.0, 2.0, 1.0, 1.0]],

    # KClO3 = KCl + O2
    [[{'K': 1, 'Cl': 1, 'O': 3}, {'K': 1, 'Cl': 1}, {'O': 2}],
    [True, False, False],
    [2.0, 2.0, 3.0]],

    # Zn + NaOH + H2O = Na2Zn(OH)4 + H2
    [[{'Zn': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2, 'O': 1}, {'Na': 2, 'Zn': 1, 'O': 4, 'H': 4}, {'H': 2}],
    [True, True, True, False, False],
    [1.0, 2.0, 2.0, 1.0, 1.0]],

    # Na2CO3 + HCl = NaCl + H2O + CO2
    [[{'Na': 2, 'C': 1, 'O': 3}, {'H': 1, 'Cl': 1}, {'Na': 1, 'Cl': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    [True, True, False, False, False],
    [1.0, 2.0, 2.0, 1.0, 1.0]],

    # Ca(OH)2 + P4O10 + H2O = Ca(H2PO4)2
    [[{'Ca': 1, 'O': 2, 'H': 2}, {'P': 4, 'O': 10}, {'H': 2, 'O': 1}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}],
    [True, True, True, False],
    [2.0, 1.0, 2.0, 2.0]],

    # CaS + H2O + CO2 = Ca(HCO3)2 + H2S
    [[{'Ca': 1, 'S': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}, {'Ca': 1, 'H': 2, 'C': 2, 'O': 6}, {'H': 2, 'S': 1}],
    [True, True, True, False, False],
    [1.0, 2.0, 2.0, 1.0, 1.0]],

    # Sn(OH)4 + NaOH = Na2SnO3 + H2O
    [[{'Sn': 1, 'O': 4, 'H': 4}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'Sn': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 3.0]],

    # Na + H2O = NaOH + H2
    [[{'Na': 1}, {'H': 2, 'O': 1}, {'Na': 1, 'O': 1, 'H': 1}, {'H': 2}],
    [True, True, False, False],
    [2.0, 2.0, 2.0, 1.0]],

    # Ca3(PO4)2 + SiO2 = CaSiO3 + P2O5
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'Si': 1, 'O': 2}, {'Ca': 1, 'Si': 1, 'O': 3}, {'P': 2, 'O': 5}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 1.0]],

    # FeCl3 + NH4OH = Fe(OH)3 + NH4Cl
    [[{'Fe': 1, 'Cl': 3}, {'N': 1, 'H': 5, 'O': 1}, {'Fe': 1, 'O': 3, 'H': 3}, {'N': 1, 'H': 4, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 3.0]],

    # H3PO3 = H3PO4 + PH3
    [[{'H': 3, 'P': 1, 'O': 3}, {'H': 3, 'P': 1, 'O': 4}, {'P': 1, 'H': 3}],
    [True, False, False],
    [4.0, 3.0, 1.0]],

    # AlCl3 + AgNO3 = AgCl + Al(NO3)3
    [[{'Al': 1, 'Cl': 3}, {'Ag': 1, 'N': 1, 'O': 3}, {'Ag': 1, 'Cl': 1}, {'Al': 1, 'N': 3, 'O': 9}],
    [True, True, False, False],
    [1.0, 3.0, 3.0, 1.0]],

    # KOH + AlCl3 = KCl + Al(OH)3
    [[{'K': 1, 'O': 1, 'H': 1}, {'Al': 1, 'Cl': 3}, {'K': 1, 'Cl': 1}, {'Al': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [3.0, 1.0, 3.0, 1.0]],

    # H2SO4 + NaHCO3 = Na2SO4 + CO2 + H2O
    [[{'H': 2, 'S': 1, 'O': 4}, {'Na': 1, 'H': 1, 'C': 1, 'O': 3}, {'Na': 2, 'S': 1, 'O': 4}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 2.0, 2.0]],

    # SiO2 + HF = SiF4 + H2O
    [[{'Si': 1, 'O': 2}, {'H': 1, 'F': 1}, {'Si': 1, 'F': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 4.0, 1.0, 2.0]],

    # CaCN2 + H2O = CaCO3 + NH3
    [[{'Ca': 1, 'C': 1, 'N': 2}, {'H': 2, 'O': 1}, {'Ca': 1, 'C': 1, 'O': 3}, {'N': 1, 'H': 3}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 2.0]],

    # HCl + HNO3 = NOCl + Cl2 + H2O
    [[{'H': 1, 'Cl': 1}, {'H': 1, 'N': 1, 'O': 3}, {'N': 1, 'O': 1, 'Cl': 1}, {'Cl': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [3.0, 1.0, 1.0, 1.0, 2.0]],

    # KClO3 = KClO4 + KCl
    [[{'K': 1, 'Cl': 1, 'O': 3}, {'K': 1, 'Cl': 1, 'O': 4}, {'K': 1, 'Cl': 1}],
    [True, False, False],
    [4.0, 3.0, 1.0]],

    # P4 + O2 = P2O5
    [[{'P': 4}, {'O': 2}, {'P': 2, 'O': 5}],
    [True, True, False],
    [1.0, 5.0, 2.0]],

    # P4O10 + HCl = POCl3 + HPO3
    [[{'P': 4, 'O': 10}, {'H': 1, 'Cl': 1}, {'P': 1, 'O': 1, 'Cl': 3}, {'H': 1, 'P': 1, 'O': 3}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 3.0]],

    # Sb + O2 = Sb4O6
    [[{'Sb': 1}, {'O': 2}, {'Sb': 4, 'O': 6}],
    [True, True, False],
    [4.0, 3.0, 1.0]],

    # NH4Cl + Ca(OH)2 = CaCl2 + NH3 + H2O
    [[{'N': 1, 'H': 4, 'Cl': 1}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'Cl': 2}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [2.0, 1.0, 1.0, 2.0, 2.0]],

    # KBr + Al(ClO4)3 = AlBr3 + KClO4
    [[{'K': 1, 'Br': 1}, {'Al': 1, 'Cl': 3, 'O': 12}, {'Al': 1, 'Br': 3}, {'K': 1, 'Cl': 1, 'O': 4}],
    [True, True, False, False],
    [3.0, 1.0, 1.0, 3.0]],

    # AgNO3 + FeCl3 = Fe(NO3)3 + AgCl
    [[{'Ag': 1, 'N': 1, 'O': 3}, {'Fe': 1, 'Cl': 3}, {'Fe': 1, 'N': 3, 'O': 9}, {'Ag': 1, 'Cl': 1}],
    [True, True, False, False],
    [3.0, 1.0, 1.0, 3.0]],

    # Ca3(PO4)2 + H3PO4 = Ca(H2PO4)2
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'H': 3, 'P': 1, 'O': 4}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}],
    [True, True, False],
    [1.0, 4.0, 3.0]],

    # (CuOH)2CO3 = CuO + CO2 + H2O
    [[{'Cu': 2, 'O': 5, 'H': 2, 'C': 1}, {'Cu': 1, 'O': 1}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # SrBr2 + (NH4)2CO3 = SrCO3 + NH4Br
    [[{'Sr': 1, 'Br': 2}, {'N': 2, 'H': 8, 'C': 1, 'O': 3}, {'Sr': 1, 'C': 1, 'O': 3}, {'N': 1, 'H': 4, 'Br': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 2.0]],

    # H2O2 = H2O + O2
    [[{'H': 2, 'O': 2}, {'H': 2, 'O': 1}, {'O': 2}],
    [True, False, False],
    [2.0, 2.0, 1.0]],

    # Ca(ClO3)2 = CaCl2 + O2
    [[{'Ca': 1, 'Cl': 2, 'O': 6}, {'Ca': 1, 'Cl': 2}, {'O': 2}],
    [True, False, False],
    [1.0, 1.0, 3.0]],

    # PCl5 + H2O = POCl3 + HCl
    [[{'P': 1, 'Cl': 5}, {'H': 2, 'O': 1}, {'P': 1, 'O': 1, 'Cl': 3}, {'H': 1, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 2.0]],

    # Zn + KOH = K2ZnO2 + H2
    [[{'Zn': 1}, {'K': 1, 'O': 1, 'H': 1}, {'K': 2, 'Zn': 1, 'O': 2}, {'H': 2}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # Al2O3 + Na2CO3 = NaAlO2 + CO2
    [[{'Al': 2, 'O': 3}, {'Na': 2, 'C': 1, 'O': 3}, {'Na': 1, 'Al': 1, 'O': 2}, {'C': 1, 'O': 2}],
    [True, True, False, False],
    [1.0, 1.0, 2.0, 1.0]],

    # PCl5 + KNO2 = NOCl + POCl3 + KCl
    [[{'P': 1, 'Cl': 5}, {'K': 1, 'N': 1, 'O': 2}, {'N': 1, 'O': 1, 'Cl': 1}, {'P': 1, 'O': 1, 'Cl': 3}, {'K': 1, 'Cl': 1}],
    [True, True, False, False, False],
    [1.0, 1.0, 1.0, 1.0, 1.0]],

    # Zn + HCl = ZnCl2 + H2
    [[{'Zn': 1}, {'H': 1, 'Cl': 1}, {'Zn': 1, 'Cl': 2}, {'H': 2}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # BeO + C + Cl2 = BeCl2 + CO
    [[{'Be': 1, 'O': 1}, {'C': 1}, {'Cl': 2}, {'Be': 1, 'Cl': 2}, {'C': 1, 'O': 1}],
    [True, True, True, False, False],
    [1.0, 1.0, 1.0, 1.0, 1.0]],

    # AgBr + Na2S2O3 = Na3Ag(S2O3)2 + NaBr
    [[{'Ag': 1, 'Br': 1}, {'Na': 2, 'S': 2, 'O': 3}, {'Na': 3, 'Ag': 1, 'S': 4, 'O': 6}, {'Na': 1, 'Br': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # N2 + O2 = N2O
    [[{'N': 2}, {'O': 2}, {'N': 2, 'O': 1}],
    [True, True, False],
    [2.0, 1.0, 2.0]],

    # BeSO4 + NH4OH = Be(OH)2 + (NH4)2SO4
    [[{'Be': 1, 'S': 1, 'O': 4}, {'N': 1, 'H': 5, 'O': 1}, {'Be': 1, 'O': 2, 'H': 2}, {'N': 2, 'H': 8, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # Cu(CN)2 = CuCN + C2N2
    [[{'Cu': 1, 'C': 2, 'N': 2}, {'Cu': 1, 'C': 1, 'N': 1}, {'C': 2, 'N': 2}],
    [True, False, False],
    [2.0, 2.0, 1.0]],

    # SiC + Cl2 = SiCl4 + C
    [[{'Si': 1, 'C': 1}, {'Cl': 2}, {'Si': 1, 'Cl': 4}, {'C': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # NH3 + O2 = HNO3 + H2O
    [[{'N': 1, 'H': 3}, {'O': 2}, {'H': 1, 'N': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # Fe2(C2O4)3 = FeC2O4 + CO2
    [[{'Fe': 2, 'C': 6, 'O': 12}, {'Fe': 1, 'C': 2, 'O': 4}, {'C': 1, 'O': 2}],
    [True, False, False],
    [1.0, 2.0, 2.0]],

    # H2 + O2 = H2O
    [[{'H': 2}, {'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False],
    [2.0, 1.0, 2.0]],

    # K + Br2 = KBr
    [[{'K': 1}, {'Br': 2}, {'K': 1, 'Br': 1}],
    [True, True, False],
    [2.0, 1.0, 2.0]],

    # CO + O2 = CO2
    [[{'C': 1, 'O': 1}, {'O': 2}, {'C': 1, 'O': 2}],
    [True, True, False],
    [2.0, 1.0, 2.0]],

    # HNO2 + O2 = HNO3
    [[{'H': 1, 'N': 1, 'O': 2}, {'O': 2}, {'H': 1, 'N': 1, 'O': 3}],
    [True, True, False],
    [2.0, 1.0, 2.0]],

    # O2 = O3
    [[{'O': 2}, {'O': 3}],
    [True, False],
    [3.0, 2.0]],

    # NaHCO3 = Na2CO3 + CO2 + H2O
    [[{'Na': 1, 'H': 1, 'C': 1, 'O': 3}, {'Na': 2, 'C': 1, 'O': 3}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [2.0, 1.0, 1.0, 1.0]],

    # CO2 + NH3 = OC(NH2)2 + H2O
    [[{'C': 1, 'O': 2}, {'N': 1, 'H': 3}, {'O': 1, 'C': 1, 'N': 2, 'H': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # Xe + F2 = XeF6
    [[{'Xe': 1}, {'F': 2}, {'Xe': 1, 'F': 6}],
    [True, True, False],
    [1.0, 3.0, 1.0]],

    # MnS + HCl = H2S + MnCl2
    [[{'Mn': 1, 'S': 1}, {'H': 1, 'Cl': 1}, {'H': 2, 'S': 1}, {'Mn': 1, 'Cl': 2}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # CaC2 + H2O = C2H2 + Ca(OH)2
    [[{'Ca': 1, 'C': 2}, {'H': 2, 'O': 1}, {'C': 2, 'H': 2}, {'Ca': 1, 'O': 2, 'H': 2}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # ClO2 + H2O = HClO2 + HClO3
    [[{'Cl': 1, 'O': 2}, {'H': 2, 'O': 1}, {'H': 1, 'Cl': 1, 'O': 2}, {'H': 1, 'Cl': 1, 'O': 3}],
    [True, True, False, False],
    [2.0, 1.0, 1.0, 1.0]],

    # CuSO4 + KCN = Cu(CN)2 + K2SO4
    [[{'Cu': 1, 'S': 1, 'O': 4}, {'K': 1, 'C': 1, 'N': 1}, {'Cu': 1, 'C': 2, 'N': 2}, {'K': 2, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # NaOH + FeSO4 = Na2SO4 + Fe(OH)2
    [[{'Na': 1, 'O': 1, 'H': 1}, {'Fe': 1, 'S': 1, 'O': 4}, {'Na': 2, 'S': 1, 'O': 4}, {'Fe': 1, 'O': 2, 'H': 2}],
    [True, True, False, False],
    [2.0, 1.0, 1.0, 1.0]],

    # Ca(OH)2 + H3PO4 = CaHPO4 + H2O
    [[{'Ca': 1, 'O': 2, 'H': 2}, {'H': 3, 'P': 1, 'O': 4}, {'Ca': 1, 'H': 1, 'P': 1, 'O': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 2.0]],

    # PbCrO4 + HNO3 = Pb(NO3)2 + H2CrO4
    [[{'Pb': 1, 'Cr': 1, 'O': 4}, {'H': 1, 'N': 1, 'O': 3}, {'Pb': 1, 'N': 2, 'O': 6}, {'H': 2, 'Cr': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # HgO = Hg + O2
    [[{'Hg': 1, 'O': 1}, {'Hg': 1}, {'O': 2}],
    [True, False, False],
    [2.0, 2.0, 1.0]],

    # (CN)2 + NaOH = NaCN + NaOCN + H2O
    [[{'C': 2, 'N': 2}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 1, 'C': 1, 'N': 1}, {'Na': 1, 'O': 1, 'C': 1, 'N': 1}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 1.0, 1.0]],

    # BaCO3 + HNO3 = Ba(NO3)2 + CO2 + H2O
    [[{'Ba': 1, 'C': 1, 'O': 3}, {'H': 1, 'N': 1, 'O': 3}, {'Ba': 1, 'N': 2, 'O': 6}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 1.0, 1.0]],

    # H3AsO4 = As2O5 + H2O
    [[{'H': 3, 'As': 1, 'O': 4}, {'As': 2, 'O': 5}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 3.0]],

    # CaO + C = CaC2 + CO
    [[{'Ca': 1, 'O': 1}, {'C': 1}, {'Ca': 1, 'C': 2}, {'C': 1, 'O': 1}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 1.0]],

    # Zn(OH)2 + NaOH = Na2ZnO2 + H2O
    [[{'Zn': 1, 'O': 2, 'H': 2}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 2, 'Zn': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # HNO3 + P2O5 = N2O5 + HPO3
    [[{'H': 1, 'N': 1, 'O': 3}, {'P': 2, 'O': 5}, {'N': 2, 'O': 5}, {'H': 1, 'P': 1, 'O': 3}],
    [True, True, False, False],
    [2.0, 1.0, 1.0, 2.0]],

    # UF4 + Mg = MgF2 + U
    [[{'U': 1, 'F': 4}, {'Mg': 1}, {'Mg': 1, 'F': 2}, {'U': 1}],
    [True, True, False, False],
    [1.0, 2.0, 2.0, 1.0]],

    # Mn2O3 + Al = Al2O3 + Mn
    [[{'Mn': 2, 'O': 3}, {'Al': 1}, {'Al': 2, 'O': 3}, {'Mn': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # MnO2 + K2CO3 + KNO3 = K2MnO4 + KNO2 + CO2
    [[{'Mn': 1, 'O': 2}, {'K': 2, 'C': 1, 'O': 3}, {'K': 1, 'N': 1, 'O': 3}, {'K': 2, 'Mn': 1, 'O': 4}, {'K': 1, 'N': 1, 'O': 2}, {'C': 1, 'O': 2}],
    [True, True, True, False, False, False],
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],

    # AlN + H2O = NH3 + Al(OH)3
    [[{'Al': 1, 'N': 1}, {'H': 2, 'O': 1}, {'N': 1, 'H': 3}, {'Al': 1, 'O': 3, 'H': 3}],
    [True, True, False, False],
    [1.0, 3.0, 1.0, 1.0]],

    # Ca3(PO4)2 + H2SO4 = CaSO4 + Ca(H2PO4)2
    [[{'Ca': 3, 'P': 2, 'O': 8}, {'H': 2, 'S': 1, 'O': 4}, {'Ca': 1, 'S': 1, 'O': 4}, {'Ca': 1, 'H': 4, 'P': 2, 'O': 8}],
    [True, True, False, False],
    [1.0, 2.0, 2.0, 1.0]],

    # S + N2O = SO2 + N2
    [[{'S': 1}, {'N': 2, 'O': 1}, {'S': 1, 'O': 2}, {'N': 2}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # N2 + H2 = NH3
    [[{'N': 2}, {'H': 2}, {'N': 1, 'H': 3}],
    [True, True, False],
    [1.0, 3.0, 2.0]],

    # CaCO3 + HCl = CaCl2 + H2O + CO2
    [[{'Ca': 1, 'C': 1, 'O': 3}, {'H': 1, 'Cl': 1}, {'Ca': 1, 'Cl': 2}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    [True, True, False, False, False],
    [1.0, 2.0, 1.0, 1.0, 1.0]],

    # As2O3 + H2O = H3AsO3
    [[{'As': 2, 'O': 3}, {'H': 2, 'O': 1}, {'H': 3, 'As': 1, 'O': 3}],
    [True, True, False],
    [1.0, 3.0, 2.0]],

    # Be(OH)2 + NH4HF2 = (NH4)2BeF4 + H2O
    [[{'Be': 1, 'O': 2, 'H': 2}, {'N': 1, 'H': 5, 'F': 2}, {'N': 2, 'H': 8, 'Be': 1, 'F': 4}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 2.0]],

    # NaOH + Zn(NO3)2 = NaNO3 + Zn(OH)2
    [[{'Na': 1, 'O': 1, 'H': 1}, {'Zn': 1, 'N': 2, 'O': 6}, {'Na': 1, 'N': 1, 'O': 3}, {'Zn': 1, 'O': 2, 'H': 2}],
    [True, True, False, False],
    [2.0, 1.0, 2.0, 1.0]],

    # NaH2PO4 = NaPO3 + H2O
    [[{'Na': 1, 'H': 2, 'P': 1, 'O': 4}, {'Na': 1, 'P': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, False, False],
    [1.0, 1.0, 1.0]],

    # H2CO3 = H2O + CO2
    [[{'H': 2, 'C': 1, 'O': 3}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    [True, False, False],
    [1.0, 1.0, 1.0]],

    # BaSO4 + H2SO4 = Ba(HSO4)2
    [[{'Ba': 1, 'S': 1, 'O': 4}, {'H': 2, 'S': 1, 'O': 4}, {'Ba': 1, 'H': 2, 'S': 2, 'O': 8}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # CaCO3 = CaO + CO2
    [[{'Ca': 1, 'C': 1, 'O': 3}, {'Ca': 1, 'O': 1}, {'C': 1, 'O': 2}],
    [True, False, False],
    [1.0, 1.0, 1.0]],

    # CaO + H2O = Ca(OH)2
    [[{'Ca': 1, 'O': 1}, {'H': 2, 'O': 1}, {'Ca': 1, 'O': 2, 'H': 2}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # H2SO3 = H2O + SO2
    [[{'H': 2, 'S': 1, 'O': 3}, {'H': 2, 'O': 1}, {'S': 1, 'O': 2}],
    [True, False, False],
    [1.0, 1.0, 1.0]],

    # H3PO4 + Ca(OH)2 = CaHPO4(H2O)2
    [[{'H': 3, 'P': 1, 'O': 4}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'H': 5, 'P': 1, 'O': 6}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # NaPO3 + CuO = NaCuPO4
    [[{'Na': 1, 'P': 1, 'O': 3}, {'Cu': 1, 'O': 1}, {'Na': 1, 'Cu': 1, 'P': 1, 'O': 4}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # SO3 + H2O = H2SO4
    [[{'S': 1, 'O': 3}, {'H': 2, 'O': 1}, {'H': 2, 'S': 1, 'O': 4}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # Be(OH)2 = BeO + H2O
    [[{'Be': 1, 'O': 2, 'H': 2}, {'Be': 1, 'O': 1}, {'H': 2, 'O': 1}],
    [True, False, False],
    [1.0, 1.0, 1.0]],

    # BaO + H2O = Ba(OH)2
    [[{'Ba': 1, 'O': 1}, {'H': 2, 'O': 1}, {'Ba': 1, 'O': 2, 'H': 2}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # Na2SO3 + S = Na2S2O3
    [[{'Na': 2, 'S': 1, 'O': 3}, {'S': 1}, {'Na': 2, 'S': 2, 'O': 3}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # SO2 + H2O = H2SO3
    [[{'S': 1, 'O': 2}, {'H': 2, 'O': 1}, {'H': 2, 'S': 1, 'O': 3}],
    [True, True, False],
    [1.0, 1.0, 1.0]],

    # Li2O + H2O = LiOH
    [[{'Li': 2, 'O': 1}, {'H': 2, 'O': 1}, {'Li': 1, 'O': 1, 'H': 1}],
    [True, True, False],
    [1.0, 1.0, 2.0]],

    # Na2HPO4 = Na4P2O7 + H2O
    [[{'Na': 2, 'H': 1, 'P': 1, 'O': 4}, {'Na': 4, 'P': 2, 'O': 7}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # H4As2O7 = As2O5 + H2O
    [[{'H': 4, 'As': 2, 'O': 7}, {'As': 2, 'O': 5}, {'H': 2, 'O': 1}],
    [True, False, False],
    [1.0, 1.0, 2.0]],

    # CaC2 + N2 = CaCN2 + C
    [[{'Ca': 1, 'C': 2}, {'N': 2}, {'Ca': 1, 'C': 1, 'N': 2}, {'C': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Mg(OH)2 = (MgOH)2O + H2O
    [[{'Mg': 1, 'O': 2, 'H': 2}, {'Mg': 2, 'O': 3, 'H': 2}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # HAsO3 = As2O5 + H2O
    [[{'H': 1, 'As': 1, 'O': 3}, {'As': 2, 'O': 5}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # KHSO4 = K2S2O7 + H2O
    [[{'K': 1, 'H': 1, 'S': 1, 'O': 4}, {'K': 2, 'S': 2, 'O': 7}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # H3PO4 = H4P2O7 + H2O
    [[{'H': 3, 'P': 1, 'O': 4}, {'H': 4, 'P': 2, 'O': 7}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # NaCl + NH4HCO3 = NaHCO3 + NH4Cl
    [[{'Na': 1, 'Cl': 1}, {'N': 1, 'H': 5, 'C': 1, 'O': 3}, {'Na': 1, 'H': 1, 'C': 1, 'O': 3}, {'N': 1, 'H': 4, 'Cl': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # HAsO2 = As2O3 + H2O
    [[{'H': 1, 'As': 1, 'O': 2}, {'As': 2, 'O': 3}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # UO3 + H2 = UO2 + H2O
    [[{'U': 1, 'O': 3}, {'H': 2}, {'U': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # CdSO4 + H2S = CdS + H2SO4
    [[{'Cd': 1, 'S': 1, 'O': 4}, {'H': 2, 'S': 1}, {'Cd': 1, 'S': 1}, {'H': 2, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # HIO3 = I2O5 + H2O
    [[{'H': 1, 'I': 1, 'O': 3}, {'I': 2, 'O': 5}, {'H': 2, 'O': 1}],
    [True, False, False],
    [2.0, 1.0, 1.0]],

    # Ca(HCO3)2 = CaCO3 + CO2 + H2O
    [[{'Ca': 1, 'H': 2, 'C': 2, 'O': 6}, {'Ca': 1, 'C': 1, 'O': 3}, {'C': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, False, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # FeS + H2SO4 = H2S + FeSO4
    [[{'Fe': 1, 'S': 1}, {'H': 2, 'S': 1, 'O': 4}, {'H': 2, 'S': 1}, {'Fe': 1, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # (NH4)2SO4 + CaCO3 = (NH4)2CO3 + CaSO4
    [[{'N': 2, 'H': 8, 'S': 1, 'O': 4}, {'Ca': 1, 'C': 1, 'O': 3}, {'N': 2, 'H': 8, 'C': 1, 'O': 3}, {'Ca': 1, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Hg2CO3 = Hg + HgO + CO2
    [[{'Hg': 2, 'C': 1, 'O': 3}, {'Hg': 1}, {'Hg': 1, 'O': 1}, {'C': 1, 'O': 2}],
    [True, False, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # CaSO4 = CaS + O2
    [[{'Ca': 1, 'S': 1, 'O': 4}, {'Ca': 1, 'S': 1}, {'O': 2}],
    [True, False, False],
    [1.0, 1.0, 2.0]],

    # BeF2 + Mg = MgF2 + Be
    [[{'Be': 1, 'F': 2}, {'Mg': 1}, {'Mg': 1, 'F': 2}, {'Be': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Mg + N2 = Mg3N2
    [[{'Mg': 1}, {'N': 2}, {'Mg': 3, 'N': 2}],
    [True, True, False],
    [3.0, 1.0, 1.0]],

    # SiO2 + Ca(OH)2 = CaSiO3 + H2O
    [[{'Si': 1, 'O': 2}, {'Ca': 1, 'O': 2, 'H': 2}, {'Ca': 1, 'Si': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # K2O + H2O = KOH
    [[{'K': 2, 'O': 1}, {'H': 2, 'O': 1}, {'K': 1, 'O': 1, 'H': 1}],
    [True, True, False],
    [1.0, 1.0, 2.0]],

    # C + H2O = CO + H2
    [[{'C': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 1}, {'H': 2}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Ca(OH)2 + CO2 = Ca(HCO3)2
    [[{'Ca': 1, 'O': 2, 'H': 2}, {'C': 1, 'O': 2}, {'Ca': 1, 'H': 2, 'C': 2, 'O': 6}],
    [True, True, False],
    [1.0, 2.0, 1.0]],

    # N2O3 + H2O = HNO2
    [[{'N': 2, 'O': 3}, {'H': 2, 'O': 1}, {'H': 1, 'N': 1, 'O': 2}],
    [True, True, False],
    [1.0, 1.0, 2.0]],

    # SiO2 + Na2CO3 = Na2SiO3 + CO2
    [[{'Si': 1, 'O': 2}, {'Na': 2, 'C': 1, 'O': 3}, {'Na': 2, 'Si': 1, 'O': 3}, {'C': 1, 'O': 2}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # BaO2 + H2SO4 = BaSO4 + H2O2
    [[{'Ba': 1, 'O': 2}, {'H': 2, 'S': 1, 'O': 4}, {'Ba': 1, 'S': 1, 'O': 4}, {'H': 2, 'O': 2}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Na2Cr2O7 + S = Cr2O3 + Na2SO4
    [[{'Na': 2, 'Cr': 2, 'O': 7}, {'S': 1}, {'Cr': 2, 'O': 3}, {'Na': 2, 'S': 1, 'O': 4}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Ca(OH)2 + CO2 = CaCO3 + H2O
    [[{'Ca': 1, 'O': 2, 'H': 2}, {'C': 1, 'O': 2}, {'Ca': 1, 'C': 1, 'O': 3}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Fe2O3 + SiO2 = Fe2Si2O7
    [[{'Fe': 2, 'O': 3}, {'Si': 1, 'O': 2}, {'Fe': 2, 'Si': 2, 'O': 7}],
    [True, True, False],
    [1.0, 2.0, 1.0]],

    # CO2 + NH3 + H2O = NH4HCO3
    [[{'C': 1, 'O': 2}, {'N': 1, 'H': 3}, {'H': 2, 'O': 1}, {'N': 1, 'H': 5, 'C': 1, 'O': 3}],
    [True, True, True, False],
    [1.0, 1.0, 1.0, 1.0]],

    # Na2O + H2O = NaOH
    [[{'Na': 2, 'O': 1}, {'H': 2, 'O': 1}, {'Na': 1, 'O': 1, 'H': 1}],
    [True, True, False],
    [1.0, 1.0, 2.0]],

    # NH4NO3 = N2O + H2O
    [[{'N': 2, 'H': 4, 'O': 3}, {'N': 2, 'O': 1}, {'H': 2, 'O': 1}],
    [True, False, False],
    [1.0, 1.0, 2.0]],

    # N2O5 + H2O = HNO3
    [[{'N': 2, 'O': 5}, {'H': 2, 'O': 1}, {'H': 1, 'N': 1, 'O': 3}],
    [True, True, False],
    [1.0, 1.0, 2.0]],

    # CaS + H2O = Ca(OH)2 + H2S
    [[{'Ca': 1, 'S': 1}, {'H': 2, 'O': 1}, {'Ca': 1, 'O': 2, 'H': 2}, {'H': 2, 'S': 1}],
    [True, True, False, False],
    [1.0, 2.0, 1.0, 1.0]],

    # Al(OH)3 + NaOH = NaAlO2 + H2O
    [[{'Al': 1, 'O': 3, 'H': 3}, {'Na': 1, 'O': 1, 'H': 1}, {'Na': 1, 'Al': 1, 'O': 2}, {'H': 2, 'O': 1}],
    [True, True, False, False],
    [1.0, 1.0, 1.0, 2.0]],

]

@pytest.mark.parametrize("test_case", stoich_test_cases)
def test_balance_stoichiometry(test_case):
    atomss, statuses, products = test_case
    for settings in [{'rounding': 9, 'allow_fractional': False},
                     {'rounding': 16, 'allow_fractional': True}]:
        matrix = stoichiometric_matrix(atomss, statuses)
        products_calc = balance_stoichiometry(matrix, **settings)
        if not settings['allow_fractional']:
            # when we allow fractions we stil have valid ratios but they do not match the hardcoded answers
            assert_close1d(products_calc, products)
        check_reaction_balance(matrix, products_calc, atol=1e-11)


def test_round_to_significant():
    """Test the round_to_significant function with various cases. TODO: Move to fluids.numerics."""
    # Test cases as (input, significant_digits, expected_output)
    test_cases = [
        (1234.567, 3, 1230.0),           # Normal positive number rounding
        (0.004567, 2, 0.0046),           # Small positive number
        (-9876.543, 4, -9877.0),         # Negative number rounding
        (1.2345e5, 2, 1.2e5),            # Scientific notation input
        (1.2345e-3, 3, 0.00123),         # Small number in scientific notation
        (0.0, 3, 0.0),                   # Zero input
        (999.999, 2, 1000.0),            # Rounding up to next order of magnitude
        (0.0001234, 2, 0.00012),         # Small positive number near zero
        (-0.0009876, 3, -0.000988),      # Small negative number near zero
        (123456789.0, 5, 123460000.0),   # Large number rounding
        (1.0, 1, 1.0),                   # Simple case with one significant digit
        (-1.0, 1, -1.0),                 # Simple negative number
        (0.12345, 3, 0.123),             # Decimal number rounding
        (-0.12345, 3, -0.123),           # Negative decimal number rounding
        (999.5, 3, 1000.0),              # Rounding with carryover
        (-0.0054321, 2, -0.0054),        # Small negative number with rounding down

        # Rounding to 5 significant digits
        (123456.789123, 5, 123460.0),        # Normal positive number
        (-987654.321987, 5, -987650.0),      # Normal negative number
        (0.000123456789, 5, 0.00012346),     # Small positive number
        (-0.000987654321, 5, -0.00098765),   # Small negative number
        (1234567890.12345, 5, 1234600000.0), # Large positive number
        (-9876543210.54321, 5, -9876500000.0), # Large negative number
        (1.23456789123e5, 5, 1.2346e5),      # Scientific notation positive
        (-9.87654321987e-3, 5, -0.0098765),  # Scientific notation negative
        (99999.9999999, 5, 100000.0),        # Rounding up
        (-99999.9999999, 5, -100000.0),      # Rounding up negative

        # Rounding to 6 significant digits
        (123456.789123, 6, 123457.0),
        (-987654.321987, 6, -987654.0),
        (0.000123456789, 6, 0.000123457),
        (-0.000987654321, 6, -0.000987654),
        (1234567890.12345, 6, 1234570000.0),
        (-9876543210.54321, 6, -9876540000.0),
        (1.23456789123e5, 6, 1.23457e5),
        (-9.87654321987e-3, 6, -0.00987654),
        (123456.1234567, 6, 123456.0),
        (-123456.1234567, 6, -123456.0),

        # Rounding to 7 significant digits
        (123456.789123, 7, 123456.8),
        (-987654.321987, 7, -987654.3),
        (0.000123456789, 7, 0.0001234568),
        (-0.000987654321, 7, -0.0009876543),
        (1234567890.12345, 7, 1234568000.0),
        (-9876543210.54321, 7, -9876543000.0),
        (1.23456789123e5, 7, 1.234568e5),
        (-9.87654321987e-3, 7, -0.009876543),
        (123456.7890123, 7, 123456.8),
        (-123456.7890123, 7, -123456.8),

        # Rounding to 8 significant digits
        (123456.789123, 8, 123456.79),
        (-987654.321987, 8, -987654.32),
        (0.000123456789, 8, 0.00012345679),
        (-0.000987654321, 8, -0.00098765432),
        (1234567890.12345, 8, 1234567900.0),
        (-9876543210.54321, 8, -9876543200.0),
        (1.23456789123e5, 8, 1.2345679e5),
        (-9.87654321987e-3, 8, -0.0098765432),
        (1234567.8901234, 8, 1234567.9),
        (-1234567.8901234, 8, -1234567.9),

        # Rounding to 9 significant digits
        (123456.789123, 9, 123456.789),
        (-987654.321987, 9, -987654.322),
        (0.000123456789, 9, 0.000123456789),
        (-0.000987654321, 9, -0.000987654321),
        (1234567890.12345, 9, 1234567890.0),
        (-9876543210.54321, 9, -9876543210.0),
        (1.23456789123e5, 9, 1.23456789e5),
        (-9.87654321987e-3, 9, -0.00987654322),
        (12345678.9012345, 9, 12345678.9),
        (-12345678.9012345, 9, -12345678.9),

        # Rounding to 10 significant digits
        (123456.789123, 10, 123456.7891),
        (-987654.321987, 10, -987654.3220),
        (0.000123456789, 10, 0.0001234567890),
        (-0.000987654321, 10, -0.0009876543210),
        (1234567890.12345, 10, 1234567890.0),
        (-9876543210.54321, 10, -9876543211.0),
        (1.23456789123e5, 10, 1.234567891e5),
        (-9.87654321987e-3, 10, -0.009876543220),
        (123456789.012345, 10, 123456789.0),
        (-123456789.012345, 10, -123456789.0),    
    ]

    # Perform tests
    for i, (x, significant_digits, expected) in enumerate(test_cases):
        result = round_to_significant(x, significant_digits)
        assert result == expected


def test_stoichiometric_matrix():
    res = stoichiometric_matrix([{'Mg': 1, 'O': 1}, {'Mg': 1}, {'O': 2}], [True, False, False])
    assert_close2d([[1, -1, 0], [1, 0, -2]], res)


def test_standard_formation_reaction():
    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 3, 'H': 8})
    assert coeff_test == [3.0, 4.0]
    assert reactant_atomss_test == [{'C': 1}, {'H': 2}]

    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1})
    assert coeff_test == [6.0, 7.0, 1.0, 2.0, 2.0]
    assert reactant_atomss_test == [{'C': 1}, {'H': 2}, {'N': 2}, {'O': 2}, {'S': 1}]

    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 6, 'H': 7, 'B': 1, 'O': 2})
    assert coeff_test == [12.0, 7.0, 2.0, 2.0]
    assert reactant_atomss_test == [{'C': 1}, {'H': 2}, {'B': 1}, {'O': 2}]

    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 4, 'H': 12, 'Si': 1})
    assert coeff_test == [4.0, 6.0, 1.0]
    assert reactant_atomss_test ==  [{'C': 1}, {'H': 2}, {'Si': 1}]

    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 12, 'H': 10, 'Cl': 1, 'P': 1})
    assert coeff_test == [24.0, 10.0, 1.0, 2.0]
    assert reactant_atomss_test ==  [{'C': 1}, {'H': 2}, {'Cl': 2}, {'P': 1}]

    reactant_coeff, coeff_test, reactant_atomss_test = standard_formation_reaction({'C': 2, 'H': 4, 'Br': 1, 'F': 1})
    assert coeff_test == [4.0, 4.0, 1.0, 1.0]
    assert reactant_atomss_test == [{'C': 1}, {'H': 2}, {'Br': 2}, {'F': 2}]
