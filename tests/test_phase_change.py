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

import pytest

from fluids.numerics import assert_close, assert_close1d
from chemicals.phase_change import *
from chemicals.phase_change import (Hvap_data_CRC, Hfus_data_CRC,
                                    Hvap_data_Gharagheizi, Hsub_data_Gharagheizi,
                                    Tb_data_Yaws, Tm_ON_data,
                                    phase_change_data_Perrys2_150,
                                    phase_change_data_Alibakhshi_Cs,
                                    phase_change_data_VDI_PPDS_4)
from chemicals.miscdata import CRC_inorganic_data, CRC_organic_data
from chemicals.identifiers import check_CAS


def test_Watson():
    Hvap = Watson(T=320, Hvap_ref=43908, T_ref=300.0, Tc=647.14)
    assert_close(Hvap, 42928.990094915454, rtol=1e-12)

    # Supercritical
    Hvap = Watson(T=480, **{'Hvap_ref': 862.4034055086013, 'T_ref': 469.6530300000001, 'Tc': 469.7, 'exponent': 0.48746055545202976})
    assert Hvap == 0.0


def test_Clapeyron():
    Hvap = Clapeyron(294.0, 466.0, 5.55E6)
    assert_close(Hvap, 26512.36357131963)

    # Test at 1/2 bar, sZ=0.98
    Hvap = Clapeyron(264.0, 466.0, 5.55E6, 0.98, 5E4)
    assert_close(Hvap, 23370.947571814384)


def test_Pitzer():
    Hvap = Pitzer(452, 645.6, 0.35017)
    assert_close(Hvap, 36696.749078320056)


def test_SMK():
    Hvap = SMK(553.15, 751.35, 0.302)
    assert_close(Hvap, 39866.18999046232)


def test_MK():
    # Problem in article for SMK function.
    Hv1 = MK(553.15, 751.35, 0.302)
    # data in [1]_., should give 26.43 KJ/mol
    Hv2 = MK(298.15, 469.69, 0.2507)
    assert_close(Hv1, 38728.00667307733, rtol=1e-12)
    assert_close(Hv2, 25940.988533726406, rtol=1e-12)


def test_Velasco():
    Hv1 = Velasco(553.15, 751.35, 0.302)
    Hv2 = Velasco(333.2, 476.0, 0.5559)
    assert_close(Hv1, 39524.251054691274, rtol=1e-12)
    assert_close(Hv2, 33299.428636069264, rtol=1e-12)


def test_Riedel():
    # same problem as in Perry's examples
    Hv1 = Riedel(294.0, 466.0, 5.55E6)
    # Pyridine, 0.0% err vs. exp: 35090 J/mol; from Poling [2]_.
    Hv2 = Riedel(388.4, 620.0, 56.3E5)
    assert_close(Hv1, 26828.59040728512, rtol=1e-12)
    assert_close(Hv2, 35089.80179000598, rtol=1e-12)


def test_Chen():
    Hv1 = Chen(294.0, 466.0, 5.55E6)
    assert_close(Hv1, 26705.902558030946)


def test_Liu():
    Hv1 = Liu(294.0, 466.0, 5.55E6)
    assert_close(Hv1, 26378.575260517395)


def test_Vetere():
    Hv1 = Vetere(294.0, 466.0, 5.55E6)
    assert_close(Hv1, 26363.43895706672)


def test_Hvap_CRC_data():

    HvapTb_tot = Hvap_data_CRC['HvapTb'].sum()
    assert_close(HvapTb_tot, 30251890.0)

    Hvap298_tot = Hvap_data_CRC['Hvap298'].sum()
    assert_close(Hvap298_tot, 29343710.0)

    Tb_tot = Hvap_data_CRC['Tb'].sum()
    assert_close(Tb_tot, 407502.95600000001)

    assert Hvap_data_CRC.index.is_unique
    assert Hvap_data_CRC.shape == (926, 5)

    assert all([check_CAS(i) for i in list(Hvap_data_CRC.index)])


def test_Hfus_CRC_data():
    Hfus_total = Hfus_data_CRC['Hfus'].sum()
    assert_close(Hfus_total, 29131241)
    assert Hfus_data_CRC.index.is_unique
    assert Hfus_data_CRC.shape == (1112, 3)

    assert all([check_CAS(i) for i in list(Hfus_data_CRC.index)])


def test_Hfus():
    assert_close(Hfus('462-06-6', method='CRC'), 11310.0, rtol=1e-12)
    assert_close(Hfus('462-06-6'), 11310.0, rtol=1e-12)
    assert_close(Hfus(CASRN='75-07-0'), 2310.0)
    assert Hfus(CASRN='75000-07-0') is None

    assert Hfus_methods('7732-18-5') == ['CRC']

def test_Gharagheizi_Hvap_data():
    # 51 CAS number DO NOT validate
    Hvap298_tot = Hvap_data_Gharagheizi['Hvap298'].sum()
    assert_close(Hvap298_tot, 173584900)

    assert Hvap_data_Gharagheizi.index.is_unique
    assert Hvap_data_Gharagheizi.shape == (2730, 2)


def test_Gharagheizi_Hsub_data():
    tots = [Hsub_data_Gharagheizi[i].sum() for i in ['Hsub', 'error']]
    assert_close(tots[0], 130537650)
    assert_close(tots[1], 1522960.0)

    assert Hsub_data_Gharagheizi.index.is_unique
    assert Hsub_data_Gharagheizi.shape == (1241, 3)


def test_Yaws_Tb_data():
    tot = Tb_data_Yaws.sum()
    assert_close(tot, 6631287.51)

    assert Tb_data_Yaws.index.is_unique
    assert Tb_data_Yaws.shape == (13461, 1)

@pytest.mark.slow
def test_Yaws_Tb_CAS_valid():
    assert all([check_CAS(i) for i in Tb_data_Yaws.index])

def test_Tm_ON_data():
    tot = Tm_ON_data.sum()
    assert_close(tot, 4059989.425)

    assert Tm_ON_data.shape == (11549, 1)
    assert Tm_ON_data.index.is_unique

@pytest.mark.slow
def test_Tm_ON_data_CAS_valid():
    assert all([check_CAS(i) for i in Tm_ON_data.index])


def test_Perrys2_150_data():
    # rtol=2E-4 for Tmin; only helium-4 needs a higher tolerance
    # Everything hits 0 at Tmax except Difluoromethane, methane, and water;
    # those needed their Tmax adjusted to their real Tc.
    # C1 is divided by 1000, to give units of J/mol instead of J/kmol
    # Terephthalic acid removed, was a constant value only.

    assert all([check_CAS(i) for i in phase_change_data_Perrys2_150.index])
    tots_calc = [phase_change_data_Perrys2_150[i].abs().sum() for i in [u'Tc', u'C1', u'C2', u'C3', u'C4', u'Tmin', u'Tmax']]
    tots = [189407.42499999999, 18617223.739999998, 174.34494000000001, 112.51209900000001, 63.894040000000004, 70810.849999999991, 189407.005]
    assert_close1d(tots_calc, tots)

    assert phase_change_data_Perrys2_150.index.is_unique
    assert phase_change_data_Perrys2_150.shape == (344, 8)


def test_Alibakhshi_Cs_data():
    # Oops, a bunch of these now-lonely coefficients have an invalid CAS...
    # assert all([check_CAS(i) for i in phase_change_data_Alibakhshi_Cs.index])
    tots_calc = [phase_change_data_Alibakhshi_Cs[i].abs().sum() for i in [u'C']]
    tots = [28154.361500000003]
    assert_close1d(tots_calc, tots)

    assert phase_change_data_Alibakhshi_Cs.index.is_unique
    assert phase_change_data_Alibakhshi_Cs.shape == (1890, 2)


def test_VDI_PPDS_4_data():
    """I believe there are no errors here."""
    assert all([check_CAS(i) for i in phase_change_data_VDI_PPDS_4.index])
    tots_calc = [phase_change_data_VDI_PPDS_4[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'E', u'Tc', u'MW']]
    tots = [1974.2929800000002, 2653.9399000000003, 2022.530649, 943.25633100000005, 3124.9258610000002, 150142.28, 27786.919999999998]
    assert_close1d(tots_calc, tots)

    assert phase_change_data_VDI_PPDS_4.index.is_unique
    assert phase_change_data_VDI_PPDS_4.shape == (272, 8)

@pytest.mark.slow
@pytest.mark.fuzz
def test_Tb_all_values():
    s1 = CRC_inorganic_data.index[CRC_inorganic_data['Tb'].notnull()]
    s2 = CRC_organic_data.index[CRC_organic_data['Tb'].notnull()]
    s3 = Tb_data_Yaws.index

    tots = []
    tots_exp = [639213.2310000042, 2280667.079999829, 6631287.510000873]
    # These should match the sums of the respective series
    for s, method in zip([s1, s2, s3], ['CRC_INORG', 'CRC_ORG', 'YAWS']):
        tots.append(sum([Tb(i, method=method) for i in s]))
    assert_close1d(tots, tots_exp, rtol=1e-11)

    s = set(); s.update(s1); s.update(s2); s.update(s3)
    assert len(s) == 13868

def test_Tb():
    # CRC_inorg, CRC org, Yaws
    Tbs_calc = Tb('993-50-0'), Tb('626-94-8'), Tb('7631-99-4')
    Tbs = [399.15, 412.15, 653.15]
    assert_close1d(Tbs, Tbs_calc)

    hits = [Tb_methods(i) for i in ['993-50-0', '626-94-8', '7631-99-4']]
    assert hits == [['CRC_INORG'], ['CRC_ORG'], ['YAWS']]


    with pytest.raises(Exception):
        Tb('993-50-0', method='BADMETHOD')

    assert None == Tb('9923443-50-0')
    assert [] == Tb_methods('9923443-50-0')


    w_methods = Tb_methods('7732-18-5')
    assert w_methods == ['CRC_INORG', 'YAWS']

    Tbs = [Tb('7732-18-5', method=i) for i in w_methods]
    assert_close1d(Tbs, [373.124, 373.15])


@pytest.mark.slow
@pytest.mark.fuzz
def test_Tm_all_values():
    s1 = CRC_inorganic_data.index[CRC_inorganic_data['Tm'].notnull()]
    s2 = CRC_organic_data.index[CRC_organic_data['Tm'].notnull()]
    s3 = Tm_ON_data.index
    tots = []
    tots_exp = [1543322.6125999668, 2571284.480399755, 4059989.4249993376]
    # These should match the sums of the respective series
    for s, method in zip([s1, s2, s3], ['CRC_INORG', 'CRC_ORG', 'OPEN_NTBKM']):
        tots.append(sum([Tm(i, method=method) for i in s]))
    assert_close1d(tots, tots_exp, rtol=1e-11)

    s = set(); s.update(s1); s.update(s2); s.update(s3)
    assert len(s) == 14723

def test_Tm():
    # Open notebook, CRC organic, CRC inorg
    Tms_calc = Tm('996-50-9'), Tm('999-78-0'), Tm('993-50-0')
    Tms = [263.15, 191.15, 274.15]
    assert_close1d(Tms, Tms_calc)

    hits = [Tm_methods(i) for i in ['996-50-9', '999-78-0', '993-50-0']]
    assert hits == [['OPEN_NTBKM'], ['CRC_ORG'], ['CRC_INORG']]

    with pytest.raises(Exception):
        Tm('993-50-0', method='BADMETHOD')

    assert  Tm('9923443-50-0') is None
    assert [] == Tm_methods('9923443-50-0')


    w_methods = Tm_methods('7732-18-5')
    assert w_methods == ['OPEN_NTBKM', 'CRC_INORG']

    Tms = [Tm('7732-18-5', method=i) for i in w_methods]
    assert_close1d(Tms, [273.15, 273.15])


def test_Alibakhshi():
    Hvap = Alibakhshi(T=320.0, Tc=647.14, C=-16.7171)
    assert_close(Hvap, 41961.30490225752, rtol=1e-13)

def test_PPDS12():
    Hvap = PPDS12(300.0, 591.75, 4.60584, 13.97224, -10.592315, 2.120205, 4.277128)
    assert_close(Hvap, 37948.76862035927, rtol=1e-13)
