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

import numpy as np
import pandas as pd
import pytest
from fluids.numerics import assert_close, assert_close1d

from chemicals.environment import (
    GWP,
    IPCC_1995_100YR_GWP,
    IPCC_2007_20YR_GWP,
    IPCC_2007_100YR_GWP,
    IPCC_2007_500YR_GWP,
    IPCC_2014_20YR_GWP,
    IPCC_2014_100YR_GWP,
    ODP,
    GWP_all_methods,
    GWP_methods,
    IPCC_2007_GWPs,
    IPCC_2014_GWPs,
    ODP_all_methods,
    ODP_data,
    ODP_methods,
    logP,
    logP_data_CRC,
    logP_data_Syrres,
    logP_methods,
    IPCC_2021_20YR_GWP, IPCC_2021_100YR_GWP, IPCC_2021_500YR_GWP, IPCC_2021_GWPs,
    GTP, GTP_methods, GTP_all_methods, IPCC_2014_20YR_GTP, IPCC_2014_50YR_GTP, IPCC_2014_100YR_GTP, IPCC_2021_50YR_GTP, IPCC_2021_100YR_GTP,
)


def test_CFC_11_in_all_GWP_methods():
    assert set(GWP_methods('75-69-4')) == set(GWP_all_methods)

def test_IPCC_2007_GWPs():
    dat_calc = [IPCC_2007_GWPs[i].sum() for i in ['Lifetime, years', 'Radiative efficiency, W/m^2/ppb', 'SAR 100yr', '20yr GWP', '100yr GWP', '500yr GWP']]
    dat = [85518.965000000011, 17.063414000000002, 128282.0, 288251, 274671.70000000001, 269051.29999999999]
    assert_close1d(dat_calc, dat)

def test_IPCC_2014_data():
    dat_calc = [IPCC_2014_GWPs[i].sum() for i in ['Lifetime, years', 'Radiative efficiency, W/m^2/ppb', '20yr GWP', '100yr GWP', '20yr GTP', '50yr GTP', '100yr GTP', '20yr AGWP', '100yr AGWP', '20yr AGTP', '50yr AGTP', '100yr AGTP']]
    dat = [99873.62999999999, 51.46000000000001, 545139.8269, 402141.69450999994, 491538.24423, 371372.361614, 330862.17449599993, 1.3603687693e-08, 3.6891094493e-08, 3.3635677188e-10, 2.2900033490900003e-10, 1.8091880086000002e-10]
    assert_close1d(dat_calc, dat)


def test_only_removed_GWPs():
    old = set(IPCC_2007_GWPs.index)
    new = set(IPCC_2014_GWPs.index)
    # dimethyl ether is the only chemical in 4e but not 5e; it has a GWP of 1
    assert old.difference(new) == {'115-10-6'}

def test_GWP():
    GWP1_calc = GWP(CASRN='74-82-8')
    assert_close(GWP1_calc, 81.2) # methane 20 year

    GWP1_calc = GWP(CASRN='74-82-8', method='IPCC (2014) 20yr')
    assert_close(GWP1_calc, 84.0) # methane 20 year

    GWP2_calc = GWP(CASRN='74-82-8', method='IPCC (1995) 100yr')
    assert_close(GWP2_calc, 21.0)

    GWP_available = GWP_methods(CASRN='56-23-5')
    assert set(GWP_available) == set(GWP_all_methods)

    with pytest.raises(Exception):
        GWP(CASRN='74-82-8', method='BADMETHOD')

    assert GWP('7732-18-5', method=None) is None
    assert GWP_methods('14882353275-98-3') == []
    assert type(GWP(CASRN='74-82-8')) is float

    assert_close(GWP(CASRN='74-82-8', method='IPCC (2021) 500yr'), 7.95)
    assert_close(GWP(CASRN='74-82-8', method='IPCC (2021) 100yr'), 27.9)

@pytest.mark.slow
@pytest.mark.fuzz
def test_GWP_all_values():
    values_1995 = [GWP(i, method=IPCC_1995_100YR_GWP) for i in IPCC_2007_GWPs.index]
    sum_1995_100 = sum(filter(lambda x: x is not None, values_1995))
    assert_close(sum_1995_100, 128282.0)

    sum_2007_100 = sum([GWP(i, method=IPCC_2007_100YR_GWP) for i in IPCC_2007_GWPs.index])
    assert_close(sum_2007_100, 274671.7)

    sum_2007_20 = sum([GWP(i, method=IPCC_2007_20YR_GWP) for i in IPCC_2007_GWPs.index])
    assert_close(sum_2007_20, 288251.0)

    sum_2007_500 = sum([GWP(i, method=IPCC_2007_500YR_GWP) for i in IPCC_2007_GWPs.index])
    assert_close(sum_2007_500, 269051.3)

    sum_2014_20 = sum([GWP(i, method=IPCC_2014_20YR_GWP) for i in IPCC_2014_GWPs.index])
    assert_close(sum_2014_20, 545139.8269)

    sum_2014_100 = sum([GWP(i, method=IPCC_2014_100YR_GWP) for i in IPCC_2014_GWPs.index])
    assert_close(sum_2014_100, 402141.69451)

    sum_2021_20 = sum([GWP(i, method=IPCC_2021_20YR_GWP) for i in IPCC_2021_GWPs.index])
    assert_close(sum_2021_20, 617406.6410000003)

    sum_2021_100 = sum([GWP(i, method=IPCC_2021_100YR_GWP) for i in IPCC_2021_GWPs.index])
    assert_close(sum_2021_100, 504313.68200000003)

    sum_2021_500 = sum([GWP(i, method=IPCC_2021_500YR_GWP) for i in IPCC_2021_GWPs.index])
    assert_close(sum_2021_500, 415507.3559999999)

def test_logP_data():
    tot = np.abs(logP_data_CRC['logP']).sum()
    assert_close(tot, 1216.99)

    tot = np.abs(logP_data_Syrres['logP']).sum()
    assert_close(tot, 25658.06)


def test_logP():
    vals = logP('67-56-1'), logP('124-18-5'), logP('7732-18-5')
    assert_close1d(vals, [-0.74, 6.25, -1.38])

    with pytest.raises(Exception):
        logP(CASRN='74-82-8', method='BADMETHOD')

    logP_available = logP_methods('110-54-3')
    assert logP_available == ['CRC', 'SYRRES']

    assert logP('1124321250-54-3') is None

@pytest.mark.fuzz
@pytest.mark.slow
def test_logP_all_values():
    tot_CRC = np.sum(np.abs(np.array([logP(i, method='CRC') for i in logP_data_CRC.index])))
    assert_close(tot_CRC, 1216.99)

    tot_SYRRES = np.sum(np.abs(np.array([logP(i, method='SYRRES') for i in logP_data_Syrres.index])))
    assert_close(tot_SYRRES, 25658.060000000001)

def test_ODP_data():

    dat_calc = [ODP_data[i].sum() for i in ['ODP2 Max', 'ODP2 Min', 'ODP1 Max', 'ODP1 Min', 'ODP2 Design', 'ODP1 Design', 'Lifetime']]
    dat = [77.641999999999996, 58.521999999999998, 64.140000000000001, 42.734000000000002, 63.10509761272651, 47.809027930358717, 2268.1700000000001]
    assert_close1d(dat_calc, dat)


def test_ODP():
    V1 = ODP(CASRN='460-86-6')
    V2 = ODP(CASRN='76-14-2', method='ODP2 Max')
    V3 = ODP(CASRN='76-14-2', method='ODP1 Max')
    assert_close1d([V1, V2, V3], [7.5, 0.58, 1.0])

    assert ODP(CASRN='148875-98-3', method='ODP2 string') == '0.2-2.1'

    methods = ['ODP2 Max', 'ODP1 Max', 'ODP2 logarithmic average', 'ODP1 logarithmic average', 'ODP2 Min', 'ODP1 Min', 'ODP2 string', 'ODP1 string']

    assert methods == ODP_methods(CASRN='148875-98-3')

    with pytest.raises(Exception):
        ODP(CASRN='148875-98-3', method='BADMETHOD')

    assert ODP(CASRN='14882353275-98-3') is None

    assert ODP_methods('14882353275-98-3') == []

@pytest.mark.slow
@pytest.mark.fuzz
def test_ODP_all_values():
    dat_calc = [pd.to_numeric(pd.Series([ODP(i, method=j) for i in ODP_data.index]), errors='coerce').sum() for j in ODP_all_methods]

    dat = [77.641999999999996, 64.140000000000001, 63.10509761272651, 47.809027930358717, 58.521999999999998, 42.734000000000002, 54.342000000000006, 38.280000000000001]
    assert_close1d(dat_calc, dat, rtol=1e-12)

def test_GTP_stuff():
    assert_close(GTP('10024-97-2', method=IPCC_2021_50YR_GTP), 290)
    assert_close(GTP('10024-97-2', method=IPCC_2021_100YR_GTP), 233)

    assert_close(GTP('10024-97-2', method=IPCC_2014_20YR_GTP), 277)
    assert_close(GTP('10024-97-2', method=IPCC_2014_50YR_GTP), 282.0)
    assert_close(GTP('10024-97-2', method=IPCC_2014_100YR_GTP), 234.0)

    assert frozenset(GTP_methods('10024-97-2')) == frozenset(GTP_all_methods)

    assert GTP('50-00-0') is None