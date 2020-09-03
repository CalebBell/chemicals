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

from numpy.testing import assert_allclose
import pytest
import pandas as pd
import numpy as np
from fluids.core import F2K
from chemicals.identifiers import check_CAS
from chemicals.safety import *
from chemicals.utils import normalize
from fluids.numerics import assert_close, assert_close1d

SUZUKI = 'Suzuki (1994)'
CROWLLOUVAR = 'Crowl and Louvar (2001)'
IEC = 'IEC 60079-20-1 (2010)'
NFPA = 'NFPA 497 (2008)'
IARC = 'International Agency for Research on Cancer'
NTP = 'National Toxicology Program 13th Report on Carcinogens'
UNLISTED = 'Unlisted'
COMBINED = 'Combined'
ONTARIO = 'Ontario Limits'
NONE = 'None'


# Not passing due to differences in pandas parsing versions
# TODO clean up file
@pytest.mark.xfail
def test_OntarioExposureLimits():
    from chemicals.safety import _OntarioExposureLimits
    pts = [_OntarioExposureLimits[i]["TWA (ppm)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 41047.08621213534)

    pts = [_OntarioExposureLimits[i]["TWA (mg/m^3)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 108342.92212201601)

    pts = [_OntarioExposureLimits[i]["STEL (ppm)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 44849.91366780729)

    pts = [_OntarioExposureLimits[i]["STEL (mg/m^3)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 95303.402815878886)

    pts = [_OntarioExposureLimits[i]["Ceiling (ppm)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 1140.6482916385789)

    pts = [_OntarioExposureLimits[i]["Ceiling (mg/m^3)"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 6093.2716177993389)

    pts = [_OntarioExposureLimits[i]["Skin"] for i in _OntarioExposureLimits.keys()]
    tot = pd.DataFrame(pts)[0].sum()
    assert_allclose(tot, 236)


def test_IARC_data():
    assert IARC_data.index.is_unique
    assert all([check_CAS(i) for i in IARC_data.index])
    assert IARC_data.shape == (863, 4)

    dict_exp = {11: 76, 1: 75, 3: 438, 12: 274}
    dict_calc = IARC_data['group'].value_counts().to_dict()
    assert dict_exp == dict_calc


def test_NTP_data():
    dict_exp = {1: 36, 2: 191}
    dict_calc = NTP_data['Listing'].value_counts().to_dict()
    assert dict_exp == dict_calc


def test_Carcinogen():
    expected = {NTP: NTP_codes[2], IARC: IARC_codes[3]}
    assert Carcinogen('61-82-5') == expected

    expected = {NTP: NTP_codes[1], IARC: IARC_codes[1]}
    assert Carcinogen('71-43-2') == expected

    expected = {NTP: UNLISTED, IARC: UNLISTED}
    assert Carcinogen('7732-18-5') == expected

    assert Carcinogen_methods('71-43-2') == [IARC, NTP]
    assert Carcinogen('71-43-2', method=NTP) == NTP_codes[1]

    with pytest.raises(Exception):
        Carcinogen('71-43-2', method='BADMETHOD')
        
    # New item 2020
    assert Carcinogen('100-00-5')[IARC] == IARC_codes[12]
    
    # Trichloroethylene added to 14th ed of NTP as known carcinogen.
    expected = {NTP: NTP_codes[1], IARC: IARC_codes[1]}
    assert Carcinogen('79-01-6') == expected
    
    # Cobalt, added in 14th ed of NTP
    expected = {NTP: NTP_codes[2], IARC: IARC_codes[11]}
    assert Carcinogen('7440-48-4') == expected
    
def test_Skin():
    assert Skin('108-94-1')
    assert not Skin('1395-21-7')

def test_safety_predictions():

    # Pentane, 1.5 % LFL in literature
    LFL1 = Suzuki_LFL(-3536600)
    # Methylamine, 4.9% LFL in article
    LFL2 = Suzuki_LFL(-1085100)

    assert_allclose([LFL1, LFL2], [0.014276107095811811, 0.043977077259217984])

    # Pentane, literature 7.8% UFL
    UFL1 = Suzuki_UFL(-3536600)
    # Methylamine, literature 20.6% UFL
    UFL2 = Suzuki_UFL(-1085100)
    assert_allclose([UFL1, UFL2], [0.0831119493052, 0.17331479619669998])

def test_safety_atom_predictions():
    # Pentane, 1.5 % LFL in literature
    LFL1 = Crowl_Louvar_LFL({'H': 12, 'C': 5})
    # Hexane, example from [1]_, lit. 1.2 %
    LFL2 = Crowl_Louvar_LFL({'H': 14, 'C': 6})
    # Example with atoms not considered
    LFL3 = Crowl_Louvar_LFL({'H': 14, 'C': 6, 'Na': 100})
    # Formaldehyde
    LFL4 = Crowl_Louvar_LFL({'H': 2, 'C': 1, 'O': 1})

    assert_allclose([LFL1, LFL2, LFL3, LFL4], [0.014073694984646879,
                    0.011899610558199915, 0.011899610558199915,
                    0.09548611111111112])

    assert Crowl_Louvar_LFL({'H': 14, 'C': 0}) == None

    # Pentane, 7.8 % LFL in literature
    UFL1 = Crowl_Louvar_UFL({'H': 12, 'C': 5})
    # Hexane, example from [1]_, lit. 7.5 %
    UFL2 = Crowl_Louvar_UFL({'H': 14, 'C': 6})
    # Example with atoms not considered
    UFL3 = Crowl_Louvar_UFL({'H': 14, 'C': 6, 'Na': 100})
    UFL4 = Crowl_Louvar_UFL({'H': 2, 'C': 1, 'O': 1})

    assert_allclose([UFL1, UFL2, UFL3, UFL4], [0.08955987717502559, 0.07572479446127219, 0.07572479446127219, 0.607638888888889])

    assert Crowl_Louvar_UFL({'H': 14, 'C': 0}) == None



def test_NFPA_497_2008():
    tots_calc = [NFPA_2008_data[i].sum() for i in ['T_flash', 'T_autoignition', 'LFL', 'UFL']]
    tots = [52112.200000000019, 127523.69999999998, 4.2690000000000001, 29.948999999999998]
    assert_allclose(tots_calc, tots)

    assert NFPA_2008_data.index.is_unique
    assert NFPA_2008_data.shape == (231, 5)


def test_IEC_2010():
    tots_calc = [IEC_2010_data[i].sum() for i in ['T_flash', 'T_autoignition', 'LFL', 'UFL']]
    tots = [83054.500000000015, 199499.94999999998, 6.4436999999999998, 40.034999999999997]
    assert_allclose(tots_calc, tots)

    assert IEC_2010_data.index.is_unique
    assert IEC_2010_data.shape == (327, 5)

def test_DIPPR_SERAT():
    assert_allclose(DIPPR_SERAT_data['T_flash'].sum(), 285171.13471004181)

    assert DIPPR_SERAT_data.index.is_unique
    assert DIPPR_SERAT_data.shape == (870, 2)



def test_Tflash():
    T1 = T_flash('8006-61-9', method=NFPA)
    T2 = T_flash('71-43-2')
    T3 = T_flash('71-43-2', method=IEC)

    Ts = [227.15, 262.15, 262.15]
    assert_allclose([T1, T2, T3], Ts)

    methods = T_flash_methods('110-54-3')
    assert methods == list(T_flash_all_methods)
    
    tot1 = pd.Series([T_flash(i) for i in IEC_2010_data.index]).sum()
    tot2 = pd.Series([T_flash(i) for i in NFPA_2008_data.index]).sum()
    tot3 = pd.Series([T_flash(i) for i in DIPPR_SERAT_data.index]).sum()
    assert_allclose([tot1, tot2, tot3], [86127.510323478724, 59397.72151504083, 286056.08653090859])

    tot_default = pd.Series([T_flash(i) for i in set([*IEC_2010_data.index, *NFPA_2008_data.index, *DIPPR_SERAT_data.index])]).sum()
    assert_allclose(tot_default, 324881.68653090857)
    
    assert None == T_flash(CASRN='132451235-2151234-1234123')

    with pytest.raises(Exception):
        T_flash(CASRN='8006-61-9', method='BADMETHOD')


def test_Tautoignition():
    T1 = T_autoignition('8006-61-9', method=NFPA)
    T2 = T_autoignition('71-43-2')
    T3 = T_autoignition('71-43-2', method=IEC)

    Ts = [553.15, 771.15, 771.15]
    assert_allclose([T1, T2, T3], Ts)

    methods = T_autoignition_methods('8006-61-9')
    assert methods == list(T_autoignition_all_methods)

    tot_default = pd.Series([T_autoignition(i) for i in set([*IEC_2010_data.index, *NFPA_2008_data.index])]).sum()
    assert_allclose(tot_default, 229841.29999999993)

    assert None == T_autoignition(CASRN='132451235-2151234-1234123')

    with pytest.raises(Exception):
        T_autoignition(CASRN='8006-61-9', method='BADMETHOD')


def test_LFL():
    LFL1 = LFL(CASRN='8006-61-9')
    LFL2 = LFL(CASRN='71-43-2', method=NFPA)
    LFL3 = LFL(CASRN='71-43-2', method=IEC)
    LFL4 = LFL(CASRN='71-43-2')
    LFL5 = LFL(Hc=-764464.0)
    LFL6 = LFL(atoms={'H': 4, 'C': 1, 'O': 1})
    LFLs = [0.014, 0.012, 0.012, 0.012, 0.05870183749384112, 0.06756756756756757]
    assert_allclose([LFL1, LFL2, LFL3, LFL4, LFL5, LFL6], LFLs)

    methods = LFL_methods(CASRN='71-43-2', Hc=-764464, atoms={'H': 4, 'C': 1, 'O': 1})
    assert methods == list(LFL_all_methods)

    tot_default = pd.Series([LFL(CASRN=i) for i in set([*IEC_2010_data.index, *NFPA_2008_data.index])]).sum()
    assert_allclose(tot_default, 7.0637000000000008)

    assert None == LFL(CASRN='132451235-2151234-1234123')

    with pytest.raises(Exception):
        LFL(CASRN='8006-61-9', method='BADMETHOD')

def test_LFL_ISO_10156_2017():
    # Example 1
    zs = [.8, .2]
    LFLs = [.044, 0.024]
    CASs = ['74-82-8', '74-84-0']
    res = LFL_ISO_10156_2017(zs, LFLs, CASs)
    assert_close(res, 0.037714285714285714, rtol=1e-13)
    
    # Example 2
    zs = [.48, .5+.79*.02, .02*.21]
    LFLs = [0.04, None, None]
    CASs = ['1333-74-0', '7727-37-9', '7782-44-7']
    res = LFL_ISO_10156_2017(zs, LFLs, CASs)
    assert_close(res, 0.08333333333333333, rtol=1e-13)
    
    # Example 3
    zs = [.4, .6]
    LFLs = [.044, None]
    CASs = ['74-82-8', '124-38-9']
    res = LFL_ISO_10156_2017(zs, LFLs, CASs)
    assert_close(res, 0.11379707112970712, rtol=1e-13)
    
    # Example 4
    zs = [.15, .15, .3, .35+.05*.79, .05*.21]
    LFLs = [.04, .044, None, None, None]
    CASs = ['1333-74-0', '74-82-8', '124-38-9', '7727-37-9', '7782-44-7']
    res = LFL_ISO_10156_2017(zs, LFLs, CASs)
    assert_close(res, 0.14273722742907632, rtol=1e-13)
    
def test_UFL():
    UFL1 = UFL(CASRN='8006-61-9')
    UFL2 = UFL(CASRN='71-43-2', method=NFPA)
    UFL3 = UFL(CASRN='71-43-2', method=IEC)
    UFL4 = UFL(CASRN='71-43-2')
    UFL5 = UFL(Hc=-764464.0)
    UFL6 = UFL(atoms={'H': 4, 'C': 1, 'O': 1})
    UFLs = [0.076, 0.078, 0.086, 0.086, 0.1901523455253683, 0.4299754299754299]
    assert_allclose([UFL1, UFL2, UFL3, UFL4, UFL5, UFL6], UFLs)

    methods = UFL_methods(CASRN='71-43-2', Hc=-764464, atoms={'H': 4, 'C': 1, 'O': 1})
    assert methods == list(UFL_all_methods)

    tot_default = pd.Series([UFL(CASRN=i) for i in set([*IEC_2010_data.index, *NFPA_2008_data.index])]).sum()
    assert_allclose(tot_default, 46.364000000000004)

    assert None == UFL(CASRN='132451235-2151234-1234123')

    with pytest.raises(Exception):
        UFL(CASRN='8006-61-9', method='BADMETHOD')


def test_unit_conv_TLV():
    mgm3 = ppmv_to_mgm3(1, 40)
    assert_allclose(mgm3, 1.6349617809430446)

    ppmv = mgm3_to_ppmv(1.635, 40)
    assert_allclose(ppmv, 1.0000233761164334)

def test_fire_mixing():
    LFL = fire_mixing(ys=normalize([0.0024, 0.0061, 0.0015]), FLs=[.012, .053, .031])
    assert_close(LFL, 0.02751172136637642, rtol=1e-13)
    
    UFL = fire_mixing(ys=normalize([0.0024, 0.0061, 0.0015]), FLs=[.075, .15, .32])
    assert_close(UFL, 0.12927551844869378, rtol=1e-13)


def test_NFPA_30_classification():
    assert NFPA_30_classification(253.15, 283.55) == 'IA' # ethylene oxide
    assert NFPA_30_classification(253.15, Psat_100F=268062) == 'IA' # ethylene oxide
    
    assert NFPA_30_classification(227.15, 249.05) == 'IA' # methyl chloride
    assert NFPA_30_classification(227.15, Psat_100F=812201) == 'IA' # methyl chloride
    
    assert NFPA_30_classification(233.15, 309.21) == 'IA' # pentane
    assert NFPA_30_classification(233.15, Psat_100F=107351) == 'IA' # pentane
    assert NFPA_30_classification(233.15, Psat_100F=101325) == 'IA' # pentane fake point to trigger border
    assert NFPA_30_classification(233.15, Psat_100F=101324.99999) == 'IB' # pentane fake point to trigger border
    assert NFPA_30_classification(233.15, Tb=310.92777777777) == 'IA' # pentane fake point under border
    assert NFPA_30_classification(233.15, Tb=310.92777777777777) == 'IB' # pentane fake point above border
    
    assert NFPA_30_classification(253.15, 329.23) == 'IB' # acetone
    assert NFPA_30_classification(253.15, Psat_100F=51979) == 'IB' # acetone
    
    assert NFPA_30_classification(262.15, 353.23) == 'IB' # benzene
    assert NFPA_30_classification(262.15, Psat_100F=22215) == 'IB' # benzene
    
    assert NFPA_30_classification(308.15, 390.75) == 'IC' # butyl alcohol
    assert NFPA_30_classification(308.15, Psat_100F=2158) == 'IC' # butyl alcohol
    assert NFPA_30_classification(308.15) == 'IC' # butyl alcohol
    
    assert NFPA_30_classification(F2K(100)) == 'II' # made up
    assert NFPA_30_classification(F2K(140)*(1-1e-13)) == 'II' # made up
    assert NFPA_30_classification(F2K(120), Tb=1e100) == 'II' # made up
    assert NFPA_30_classification(F2K(120), Tb=-1e100) == 'II' # made up
    assert NFPA_30_classification(F2K(120), Psat_100F=-1e100) == 'II' # made up
    assert NFPA_30_classification(F2K(120), Psat_100F=1e100) == 'II' # made up
    
    assert NFPA_30_classification(F2K(140)) == 'IIIA' # made up
    assert NFPA_30_classification(F2K(200)*(1-1e-13)) == 'IIIA' # made up
    assert NFPA_30_classification(F2K(170), Tb=1e100) == 'IIIA' # made up
    assert NFPA_30_classification(F2K(170), Tb=-1e100) == 'IIIA' # made up
    assert NFPA_30_classification(F2K(170), Psat_100F=-1e100) == 'IIIA' # made up
    assert NFPA_30_classification(F2K(170), Psat_100F=1e100) == 'IIIA' # made up
    
    assert NFPA_30_classification(F2K(200)) == 'IIIB' # made up
    assert NFPA_30_classification(F2K(200)*(1+1e-13)) == 'IIIB' # made up
    assert NFPA_30_classification(F2K(300), Tb=1e100) == 'IIIB' # made up
    assert NFPA_30_classification(F2K(3000), Tb=-1e100) == 'IIIB' # made up
    assert NFPA_30_classification(F2K(300000), Psat_100F=-1e100) == 'IIIB' # made up
    assert NFPA_30_classification(F2K(300000000), Psat_100F=1e100) == 'IIIB' # made up

    with pytest.raises(ValueError):
        NFPA_30_classification(253.15)