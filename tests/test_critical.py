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
import pandas as pd
from fluids.numerics import assert_close, assert_close1d
from chemicals.critical import *
from chemicals.critical import (critical_data_IUPAC,
                                critical_data_Matthews,
                                critical_data_CRC, 
                                critical_data_PSRKR4, 
                                critical_data_Yaws, 
                                critical_data_PassutDanner)

def test_data_IUPAC():
    MW_sum = critical_data_IUPAC['MW'].sum()
    assert_close(MW_sum,122998.43799999992)

    Tc_sum = critical_data_IUPAC['Tc'].sum()
    assert_close(Tc_sum, 462157.51300000004)

    Pc_sum = critical_data_IUPAC['Pc'].sum()
    assert_close(Pc_sum, 2063753000.0)

    Vc_sum = critical_data_IUPAC['Vc'].sum()
    assert_close(Vc_sum, 0.19953190000000001)

    Zc_sum = critical_data_IUPAC['Zc'].sum()
    assert_close(Zc_sum, 109.892)

    assert critical_data_IUPAC.shape == (810, 7)
    assert critical_data_IUPAC.index.is_unique


def test_data_Matthews():
    MW_sum = critical_data_Matthews['MW'].sum()
    assert_close(MW_sum, 19541.760399999999)

    Tc_sum = critical_data_Matthews['Tc'].sum()
    assert_close(Tc_sum, 65343.210900000005)

    Pc_sum = critical_data_Matthews['Pc'].sum()
    assert_close(Pc_sum, 579365204.25)

    Vc_sum = critical_data_Matthews['Vc'].sum()
    assert_close(Vc_sum, 0.014921000000000002)

    Zc_sum = critical_data_Matthews['Zc'].sum()
    assert_close(Zc_sum, 12.141000000000002)

    assert critical_data_Matthews.shape == (120, 6)
    assert critical_data_Matthews.index.is_unique
    
def test_data_CRC():
    Tc_sum = critical_data_CRC['Tc'].sum()
    assert_close(Tc_sum, 514092.75)

    Pc_sum = critical_data_CRC['Pc'].sum()
    assert_close(Pc_sum, 2700259000.0)

    Vc_sum = critical_data_CRC['Vc'].sum()
    assert_close(Vc_sum, 0.38876929999999998)

    Zc_sum = critical_data_CRC['Zc'].sum()
    assert_close(Zc_sum, 207.98663028416496, 1e-6)

    assert critical_data_CRC.shape == (861, 8)
    assert critical_data_CRC.index.is_unique

    Tc_error_sum = critical_data_CRC['Tc_error'].sum()
    assert_close(Tc_error_sum, 2444.24)

    Pc_error_sum = critical_data_CRC['Pc_error'].sum()
    assert_close(Pc_error_sum, 1.2587e+08)

    Vc_error_sum = critical_data_CRC['Vc_error'].sum()
    assert_close(Vc_error_sum, 0.014151)


def test_data_PSRKR4():
    Tc_sum = critical_data_PSRKR4['Tc'].sum()
    assert_close(Tc_sum, 597984.0)

    Pc_sum = critical_data_PSRKR4['Pc'].sum()
    assert_close(Pc_sum, 3708990509)

    Vc_sum = critical_data_PSRKR4['Vc'].sum()
    assert_close(Vc_sum, 0.40726849999999998)

    Zc_sum = critical_data_PSRKR4['Zc'].sum()
    assert_close(Zc_sum, 251.29839643655527, 1e-6)

    omega_sum = critical_data_PSRKR4['omega'].sum()
    assert_close(omega_sum, 410.50560000000002)

    assert critical_data_PSRKR4.shape == (995, 6)
    assert critical_data_PSRKR4.index.is_unique


def test_data_PassutDanner():
    Tc_sum = critical_data_PassutDanner['Tc'].sum()
    assert_close(Tc_sum, 111665.28333333334)

    Pc_sum = critical_data_PassutDanner['Pc'].sum()
    assert_close(Pc_sum, 579756767.55527318)

    omega_sum = critical_data_PassutDanner['omega'].sum()
    assert_close(omega_sum, 65.567000000000007)

    assert critical_data_PassutDanner.shape == (192, 4)
    assert critical_data_PassutDanner.index.is_unique


def test_data_Yaws():
    Tc_sum = critical_data_Yaws['Tc'].sum()
    assert_close(Tc_sum, 5862006.9500000002)

    Pc_sum = critical_data_Yaws['Pc'].sum()
    assert_close(Pc_sum, 62251189000.0)

    Vc_sum = critical_data_Yaws['Vc'].sum()
    assert_close(Vc_sum, 4.65511199)

    Zc_sum = critical_data_Yaws['Zc'].sum()
    assert_close(Zc_sum, 1859.6389176846883, 1e-6)

    omega_sum = critical_data_Yaws['omega'].sum()
    assert_close(omega_sum, 3170.3041999999996)

    assert critical_data_Yaws.shape == (7549, 6)
    assert critical_data_Yaws.index.is_unique


def test_relationships():
    Vc = Ihmels(Tc=599.4, Pc=1.19E6)
    Pc = Ihmels(Tc=599.4, Vc=0.00109273333333)
    Tc = Ihmels(Vc=0.00109273333333, Pc=1.19E6)
    assert_close1d([Tc, Pc, Vc], [599.3999999981714, 1190000.0000037064, 0.0010927333333333334])
    with pytest.raises(Exception):
        Ihmels(Tc=559.4)

    Vc = Meissner(Tc=599.4, Pc=1.19E6)
    Pc = Meissner(Tc=599.4, Vc=0.0010695726588235296)
    Tc = Meissner(Vc=0.0010695726588235296, Pc=1.19E6)
    assert_close1d([Tc, Pc, Vc], [599.4000000000001, 1190000.0, 0.0010695726588235296])
    with pytest.raises(Exception):
        Meissner(Tc=559.4)

    Vc = Grigoras(Tc=599.4, Pc=1.19E6)
    Pc = Grigoras(Tc=599.4, Vc=0.00134532)
    Tc = Grigoras(Vc=0.00134532, Pc=1.19E6)
    assert_close1d([Tc, Pc, Vc], [599.4, 1190000.0, 0.00134532])
    with pytest.raises(Exception):
        Grigoras(Tc=559.4)

    Vc1 = critical_surface(Tc=599.4, Pc=1.19E6, method='IHMELS')
    Vc2 = critical_surface(Tc=599.4, Pc=1.19E6, method='MEISSNER')
    Vc3 = critical_surface(Tc=599.4, Pc=1.19E6, method='GRIGORAS')
    assert_close1d([Vc1, Vc2, Vc3],
                    [0.0010927333333333334, 0.0010695726588235296, 0.00134532])

    methods = critical_surface_methods(Tc=599.4, Pc=1.19E6)
    methods.sort()
    methods_listed = ['IHMELS', 'MEISSNER', 'GRIGORAS']
    methods_listed.sort()
    assert methods == methods_listed
    with pytest.raises(ValueError):
        critical_surface()
    with pytest.raises(Exception):
        critical_surface(Tc=599.4, Pc=1.19E6, method='FAIL')
        
    assert [] == critical_surface_methods(Tc=100)

@pytest.mark.slow
def test_Tc_all_values():
    sources = [critical_data_IUPAC, critical_data_Matthews, critical_data_CRC, critical_data_PSRKR4, critical_data_PassutDanner, critical_data_Yaws]
    CASs = set()
    [CASs.update(set(k.index.values)) for k in sources]

    # Use the default method for each chemical in this file
    Tcs = [Tc(i) for i in CASs]
    Tcs_default_sum = pd.Series(Tcs).sum()
    assert_close(Tcs_default_sum, 6053723.896122222)

def test_Tc():
    Tc_val = Tc(CASRN='64-17-5')
    assert_close(514.0, Tc_val)
    assert type(Tc_val) is float
    
    assert_close(516.2, Tc(CASRN='64-17-5', method='PSRK'))
    assert Tc(CASRN='64-17-5', method='MATTHEWS') is None

    assert_close(647.3, Tc(CASRN='7732-18-5', method='PSRK'))

    assert_close(126.2, Tc(CASRN='7727-37-9', method='MATTHEWS'))

    methods = Tc_methods(CASRN='98-01-1')
    assert methods == ['IUPAC', 'PSRK', 'YAWS']

    # Error handling
    assert Tc(CASRN='BADCAS') is None
    
    with pytest.raises(Exception):
        Tc(CASRN='98-01-1', method='BADMETHOD')

def test_Pc():

    assert_close(6137000.0, Pc(CASRN='64-17-5'))

    assert_close(22048321.0, Pc(CASRN='7732-18-5', method='PSRK'))

    assert_close(3394387.5, Pc(CASRN='7727-37-9', method='MATTHEWS'))

    methods = Pc_methods(CASRN='98-01-1')
    assert methods == ['IUPAC', 'PSRK', 'YAWS']

    # Error handling
    assert None == Pc(CASRN='BADCAS')
    with pytest.raises(Exception):
        Pc(CASRN='98-01-1', method='BADMETHOD')


@pytest.mark.slow
def test_Pc_all_values():
    sources = [critical_data_IUPAC, critical_data_Matthews, critical_data_CRC, critical_data_PSRKR4, critical_data_PassutDanner, critical_data_Yaws]
    CASs = set()
    [CASs.update(set(k.index.values)) for k in sources]

    # Use the default method for each chemical in this file
    Pcs = [Pc(i) for i in CASs]
    Pcs_default_sum = pd.Series(Pcs).sum()
    assert_close(Pcs_default_sum, 63159160396.183258)

def test_Vc():
    assert_close(0.000168, Vc(CASRN='64-17-5'))

    assert_close(5.600e-05, Vc(CASRN='7732-18-5', method='PSRK'))

    assert_close(8.950e-05, Vc(CASRN='7727-37-9', method='MATTHEWS'))

    methods = Vc_methods(CASRN='98-01-1')
    assert methods == ['PSRK', 'YAWS']

    # Error handling
    assert None == Vc(CASRN='BADCAS')
    with pytest.raises(Exception):
        Vc(CASRN='98-01-1', method='BADMETHOD')

@pytest.mark.slow
def test_Vc_all_values():
    sources = [critical_data_IUPAC, critical_data_Matthews, critical_data_CRC, critical_data_PSRKR4, critical_data_Yaws]
    CASs = set()
    [CASs.update(set(k.index.values)) for k in sources]

    # Use the default method for each chemical in this file
    Vcs = [Vc(i) for i in CASs]
    Vcs_default_sum = pd.Series(Vcs).sum()
    assert_close(Vcs_default_sum, 4.7955320200000005)

def test_Zc():
    assert_close(0.241, Zc(CASRN='64-17-5'))

    assert_close(0.22941602891834947, Zc(CASRN='7732-18-5', method='PSRK'))

    assert_close(0.29, Zc(CASRN='7727-37-9', method='MATTHEWS'))

    methods = Zc_methods(CASRN='98-01-1')
    assert methods == ['PSRK', 'YAWS']

    # Error handling
    assert None == Zc(CASRN='BADCAS')
    with pytest.raises(Exception):
        Zc(CASRN='98-01-1', method='BADMETHOD')

@pytest.mark.slow
def test_Zc_all_values():
    sources = [critical_data_IUPAC, critical_data_Matthews, critical_data_CRC, critical_data_PSRKR4, critical_data_Yaws]
    CASs = set()
    [CASs.update(set(k.index.values)) for k in sources]

    # Use the default method for each chemical in this file
    Zcs = [Zc(i) for i in CASs]
    Zcs_default_sum = pd.Series(Zcs).sum()
    assert_close(Zcs_default_sum, 1930.7388004412558, 1e-6)



def test_Mersmann_Kind_predictor():
    test_atoms = {'C': 10, 'H': 22}
    
    Vc_pred = Mersmann_Kind_predictor(test_atoms)
    assert_close(Vc_pred, 0.0005851859052024729)
    
    with pytest.raises(Exception):
        Mersmann_Kind_predictor( {'C': 10, 'H': 22, 'NOTANATOM': 100})


def test_mixing_Tc():
    # Nitrogen-Argon 50/50 mixture
    Tcm =  Li([0.5, 0.5], [126.2, 150.8], [8.95e-05, 7.49e-05])
    assert_close(Tcm, 137.40766423357667)

    # example is from [2]_, for:
    # butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    # Its result is identical to that calculated in the article.
    Tcm = Li([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [0.000255, 0.000313, 0.000371])
    assert_close(Tcm, 449.68261498555444)

    # User is not responsible for validating their own input for speed
#    with pytest.raises(Exception):
#        Li([0.2359, 0.1192], [425.12, 469.7, 507.6], [0.000255, 0.000313, 0.000371])

    # First example is for an acetone/n-butane 50/50 mixture. No point is
    # available to compare the calculated value with, but it is believed
    # correct.
    Tcm = Chueh_Prausnitz_Tc([0.5, 0.5], [508.1, 425.12], [0.000213, 0.000255], [[0, -14.2619], [-14.2619, 0]])
    assert_close(Tcm, 457.01862919555555)

    # 2rd example is from [2]_, for:
    # butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    # Its result is identical to that calculated in the article.
    Tcm = Chueh_Prausnitz_Tc([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [0.000255, 0.000313, 0.000371], [[0, 1.92681, 6.80358], [1.92681, 0, 1.89312], [ 6.80358, 1.89312, 0]])
    assert_close(Tcm, 450.1225764723492)
    # 3rd example is from [2]_, for ethylene, Benzene, ethylbenzene. This is
    # similar to but not identical to the result from the article. The
    # experimental point is 486.9 K.
    Tcm = Chueh_Prausnitz_Tc([0.5, 0.447, .053], [282.34, 562.05, 617.15], [0.0001311, 0.000256, 0.000374], [[0, 37.9570, 0], [37.9570, 0, 4.2459], [0, 4.2459, 0]])
    assert_close(Tcm, 475.3154572323848)

    # User is not responsible for validating their own input for speed
#    with pytest.raises(Exception):
#        Chueh_Prausnitz_Tc([0.447, .053], [282.34, 562.05, 617.15], [0.0001311, 0.000256, 0.000374], [[0, 37.9570, 0], [37.9570, 0, 4.2459], [0, 4.2459, 0]])



    # First example is for an acetone/n-butane 50/50 mixture. No point is
    # available to compare the calculated value with, but it is believed
    # correct. Parameters here are from [2]_.
    Tcm = Grieves_Thodos([0.5, 0.5], [508.1, 425.12], [[0, 0.7137], [1.6496, 0]])
    assert_close(Tcm, 456.9398283342622)

    # 2rd example is from [2]_, for:
    # butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    # Its result is identical to that calculated in the article.
    Tcm = Grieves_Thodos([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [[0, 1.2503, 1.516], [0.799807, 0, 1.23843], [0.659633, 0.807474, 0]])
    assert_close(Tcm, 450.1839618758971)

    with pytest.raises(Exception):
        Grieves_Thodos([0.5, 0.5], [508.1, None], [[0, 0.7137], [1.6496, 0]])


    # First example is Acetone/butane 50/50 mixture.
    Tcm = modified_Wilson_Tc([0.5, 0.5], [508.1, 425.12],  [[0, 0.8359], [0, 1.1963]])
    assert_close(Tcm, 456.59176287162256)

    # 2rd example is from [2]_, for:
    # butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    # Its result is identical to that calculated in the article.
    Tcm = modified_Wilson_Tc([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [[0, 1.174450, 1.274390], [0.835914, 0, 1.21038], [0.746878, 0.80677, 0]])
    assert_close(Tcm, 450.0305966823031)

    # User is not responsible for validating their own input for speed
#    with pytest.raises(Exception):
#        Tcm = modified_Wilson_Tc([0.5, 0.5], [508.1, None],  [[0, 0.8359], [0, 1.1963]])


def test_mixing_Chueh_Prausnitz_Vc():
    Vcm = Chueh_Prausnitz_Vc([0.4271, 0.5729], [0.000273, 0.000256], [[0, 5.61847], [5.61847, 0]])
    assert_close(Vcm, 0.00026620503424517445)

    # User is not responsible for validating their own input for speed
#    with pytest.raises(Exception):
#        Chueh_Prausnitz_Vc([0.4271], [0.000273, 0.000256], [[0, 5.61847], [5.61847, 0]])

def test_mixing_modified_Wilson_Vc():
    Vcm = modified_Wilson_Vc([0.4271, 0.5729], [0.000273, 0.000256], [[0, 0.6671250], [1.3939900, 0]])
    assert_close(Vcm, 0.00026643350327068809)

    # User is not responsible for validating their own input for speed
#    with pytest.raises(Exception):
#        modified_Wilson_Vc([0.4271], [0.000273, 0.000256], [[0, 0.6671250], [1.3939900, 0]])


def test_third_property():
    with pytest.raises(Exception):
        third_property('141-62-8')
    with pytest.raises(Exception):
        third_property('1410-62-8', V=True)
        
    assert_close(third_property('110-15-6', V=True), 0.00039809186906019007)
    assert_close(third_property('110-15-6', P=True), 6095016.233766234)
    assert_close(third_property('110-15-6', T=True), 658.410835214447)


def test_Hekayati_Raeissi():
    # Tc, Pc estimation with estimated Vc
    ans = Hekayati_Raeissi(MW=92.13842, V_sat=0.00010685961260743776, Pc=4108000.0)
    assert_close1d(ans, (599.7951343525503, 4108000.0, 0.0003149083374202421), rtol=1e-10)
    ans = Hekayati_Raeissi(MW=92.13842, V_sat=0.00010685961260743776, Tc=599.7951343525503)
    assert_close1d(ans, (599.7951343525503, 4108000.0, 0.0003149083374202421), rtol=1e-10)
    
    ans = Hekayati_Raeissi(MW=92.13842, V_sat=0.00010685961260743776, Pc=4108000.0, Vc=0.000316)
    assert_close1d(ans, (600.331367806081, 4108000.0, 0.000316), rtol=1e-10)
    ans = Hekayati_Raeissi(MW=92.13842, V_sat=0.00010685961260743776, Tc=600.331367806081, Vc=0.000316)
    assert_close1d(ans, (600.331367806081, 4108000.0, 0.000316), rtol=1e-10)
    
    ans = Hekayati_Raeissi(MW=92.13842, Pc=4108000.0, Vc=0.000316)
    assert_close1d(ans, (600.331367806081, 4108000.0, 0.000316), rtol=1e-10)
    ans = Hekayati_Raeissi(MW=92.13842, Tc=600.331367806081, Vc=0.000316)
    assert_close1d(ans, (600.331367806081, 4108000.0, 0.000316), rtol=1e-10)
    
    # All inputs provided test
    ans = Hekayati_Raeissi(MW=92.13842, Tc=600., Vc=0.000316, Pc=1e7)
    assert_close1d(ans, (600.0, 10000000.0, 0.000316), rtol=1e-7)
    
    with pytest.raises(ValueError):
        Hekayati_Raeissi(MW=92.13842, Pc=4108000.0)

    with pytest.raises(ValueError):
        Hekayati_Raeissi(MW=92.13842, Vc=0.000316)
