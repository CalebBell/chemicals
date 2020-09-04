# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>

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
import numpy as np
import pytest
import pandas as pd
from chemicals.volume import *
from chemicals.volume import (rho_data_COSTALD, rho_data_SNM0, rho_data_Perry_8E_105_l, rho_data_CRC_inorg_l,
                              rho_data_VDI_PPDS_2, rho_data_CRC_inorg_l_const, rho_data_CRC_inorg_s_const, 
                              rho_data_CRC_virial)
from chemicals.utils import Vm_to_rho
from chemicals.identifiers import check_CAS

def Yen_Woods_saturation():
    V1_calc = Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245)
    V1 = 1.7695330765295693e-05
    V2_calc = Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.27)
    V2 = 1.8750391558570308e-05
    assert_allclose([V1_calc, V2_calc], [V1, V2])

def test_Rackett():
    V1_calc = Rackett(300, 647.14, 22048320.0, 0.23)
    V1 = 1.640447929033162e-05
    V2_calc = Rackett(272.03889, 369.83, 4248000.0, 0.2763)
    V2 = 8.299225005462148e-05
    assert_allclose([V1_calc, V2_calc], [V1, V2])

def test_Yamada_Gunn():
    V1_calc = Yamada_Gunn(300, 647.14, 22048320.0, 0.245)
    assert_allclose(V1_calc, 2.188284384699659e-05)

def test_Townsend_Hales():
    V1_calc = Townsend_Hales(300, 647.14, 55.95E-6, 0.3449)
    assert_allclose(V1_calc, 1.8007361992619923e-05)

def test_Bhirud_normal():
    V1_calc = Bhirud_normal(280.0, 469.7, 33.7E5, 0.252)
    V1 = 0.0001124965784251429
    V2_calc = Bhirud_normal(469.7*.99, 469.7, 33.7E5, 0.252)
    V2 = 0.00021992874232614203
    assert_allclose([V1_calc, V2_calc], [V1, V2])

    # Test above Tc, where interpolation table fails
    with pytest.raises(Exception):
        Bhirud_normal(500, 469.7, 33.7E5, 0.252)

def test_COSTALD():
    V1_calc = COSTALD(298, 647.13, 55.95E-6, 0.3449)
    V1 = 1.8133760480018036e-05
    # Propane, from API Handbook example; I have used their exact values,
    # but they rounded each step, getting 530.1
    V2_calc = COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532)
    V2 = 8.315466172295678e-05
    assert_allclose([V1_calc, V2_calc], [V1, V2])

    rho_ex = Vm_to_rho(COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532), 44.097)
    assert_allclose(rho_ex, 530.3009967969841)

def test_Campbell_Thodos():
    # Argon, with constants from [1]_ Table II and compared with listed
    # calculated result for critical volume in Table III. Tabulated s, lambda,
    # alpha, and beta are also a match.
    V1_calc = Campbell_Thodos(150.65, 87.28, 150.65, 48.02*101325, 39.948, 0.0)
    V1 = 7.538927923761386e-05

    # Water, with constants from [1]_ Table II and compared with listed
    # calculated result for critical volume in Table V. Tabulated s, lambda,
    # alpha, and beta are also a match. Deviation of 0.1% is due to author's
    # rearrangement of the formula.
    V2_calc = Campbell_Thodos(T=647.3, Tb=373.15, Tc=647.3, Pc=218.3*101325, MW=18.015, dipole=1.85, has_hydroxyl=True)
    V2 = 5.4787019341952454e-05

    # Ammonia, with constants from [1]_ Table II and compared with listed
    # calculated result for critical volume in Table IV. Tabulated s, lambda,
    # alpha, and beta are also a match. Deviation of 0.1% is due to author's
    # rearrangement of the formula.
    V3_calc = Campbell_Thodos(T=405.45, Tb=239.82, Tc=405.45, Pc=111.7*101325, MW=17.03, dipole=1.47)
    V3 = 7.347366126245346e-05
    assert_allclose([V1_calc, V2_calc, V3_calc], [V1, V2, V3])

def test_SNMO():
    # No examples for this model have been found, but it is simple and well
    # understood.
    V1_calc = SNM0(121, 150.8, 7.49e-05, -0.004)
    V1 = 3.4402256402733416e-05
    V2_calc = SNM0(121, 150.8, 7.49e-05, -0.004, -0.03259620)
    V2 = 3.493288100008123e-05
    assert_allclose([V1_calc, V2_calc], [V1, V2])

def test_volume_CSP_dense():
    V = COSTALD_compressed(303., 9.8E7, 85857.9, 466.7, 3640000.0, 0.281, 0.000105047)
    assert_allclose(V, 9.287482879788506e-05)


def test_CRC_inorganic():
#    # Lithium Sulfate:
    rho1_calc = CRC_inorganic(1133.15, 2003.0, 0.407, 1133.15)
    rho1 = 2003.0
    rho2_calc = CRC_inorganic(1405, 2003.0, 0.407, 1133.15)
    rho2 = 1892.35705
    assert_allclose([rho1_calc, rho2_calc], [rho1, rho2])

    # Tin tetrachloride
    rho = CRC_inorganic(300, 2370.0, 2.687, 239.08)
    assert_allclose(rho, 2206.30796)

def test_COSTALD_parameters():
    assert all([check_CAS(i) for i in rho_data_COSTALD.index])
    assert rho_data_COSTALD.index.is_unique
    tots_calc = [rho_data_COSTALD[i].sum() for i in ['omega_SRK', 'Vchar', 'Z_RA']]
    tots = [72.483900000000006, 0.086051663333333334, 49.013500000000001]
    assert_allclose(tots_calc, tots)

def test_SN0_data():
    assert rho_data_SNM0.index.is_unique
    assert all([check_CAS(i) for i in rho_data_SNM0.index])
    tot = rho_data_SNM0['delta_SRK'].abs().sum()
    assert_allclose(tot, 2.0715134)


def test_Perry_l_data():
    assert rho_data_Perry_8E_105_l.index.is_unique
    assert all([check_CAS(i) for i in rho_data_Perry_8E_105_l.index])

    tots_calc = [rho_data_Perry_8E_105_l[i].sum() for i in ['C1', 'C2', 'C3', 'C4', 'Tmin', 'Tmax']]
    tots = [376364.41000000003, 89.676429999999996, 189873.32999999999, 96.68741, 71151.899999999994, 189873.32999999999]
    assert_allclose(tots_calc, tots)


def test_VDI_PPDS_2_data():
    """Plenty of interesting errors here. The chemicals 463-58-1, 75-44-5,
    75-15-0, 7446-11-9, 2551-62-4 do not match the tabulated data. They are all
    in the same section, so a mixup was probably made there. The errors versus
    the tabulated data are very large.

    Note this table needed to have Tc and MW added to it as well, from the same
    source.
    """
    assert all([check_CAS(i) for i in rho_data_VDI_PPDS_2.index])
    tots_calc = [rho_data_VDI_PPDS_2[i].abs().sum() for i in [u'A', u'B', u'C', u'D', u'Tc', u'rhoc', u'MW']]
    tots = [208878.27130000002, 117504.59450000001, 202008.99950000001, 85280.333600000013, 150142.28, 97269, 27786.919999999998]
    assert_allclose(tots_calc, tots)
    
    assert rho_data_VDI_PPDS_2.index.is_unique
    assert rho_data_VDI_PPDS_2.shape == (272, 8)

def test_CRC_inorg_l_data2():
    tots_calc = [rho_data_CRC_inorg_l[i].abs().sum() for i in ['rho', 'k', 'Tm', 'Tmax']]
    tots = [882131, 181.916, 193785.09499999997, 233338.04999999996]
    assert_allclose(tots_calc, tots)

    assert rho_data_CRC_inorg_l.index.is_unique
    assert all([check_CAS(i) for i in rho_data_CRC_inorg_l.index])


def test_CRC_const_inorg_l():
    assert rho_data_CRC_inorg_l_const.index.is_unique
    assert all([check_CAS(i) for i in rho_data_CRC_inorg_l_const.index])

    tot_calc = rho_data_CRC_inorg_l_const['Vm'].sum()
    tot = 0.01106122489849834
    assert_allclose(tot_calc, tot)

def test_CRC_const_inorg_s():
    tot = rho_data_CRC_inorg_s_const['Vm'].sum()
    assert_allclose(tot, 0.13528770767318143)
    assert rho_data_CRC_inorg_s_const.index.is_unique
    assert all([check_CAS(i) for i in rho_data_CRC_inorg_s_const.index])


def test_CRC_virial_poly():
    assert rho_data_CRC_virial.index.is_unique
    assert all([check_CAS(i) for i in rho_data_CRC_virial.index])
    tots_calc = [rho_data_CRC_virial[i].abs().sum() for i in ['a1', 'a2', 'a3', 'a4', 'a5']]
    tots = [146559.69999999998, 506997.70000000001, 619708.59999999998, 120772.89999999999, 4483]
    assert_allclose(tots_calc, tots)


def test_solids_CSP():
    V = Goodman(200, 243.225, 0.00023585)
    assert_allclose(0.0002053665090860923, V)



# More gases:
def test_ideal_gas():
    assert_allclose(ideal_gas(298.15, 101325.), 0.024465403697038125)


def test_Amgat():
    Vl = Amgat([0.5, 0.5], [4.057e-05, 5.861e-05])
    assert_allclose(Vl, 4.9590000000000005e-05)

#    with pytest.raises(Exception):
#        Amgat([0.5], [4.057e-05, 5.861e-05])


def test_Rackett_mixture():
    Vl = Rackett_mixture(T=298., xs=[0.4576, 0.5424], MWs=[32.04, 18.01], Tcs=[512.58, 647.29], Pcs=[8.096E6, 2.209E7], Zrs=[0.2332, 0.2374])
    assert_allclose(Vl, 2.6252894930056885e-05)

#    with pytest.raises(Exception):
#        Rackett_mixture(T=298., xs=[0.4576], MWs=[32.04, 18.01], Tcs=[512.58, 647.29], Pcs=[8.096E6, 2.209E7], Zrs=[0.2332, 0.2374])

def test_COSTALD_mixture():
    Vl = COSTALD_mixture([0.4576, 0.5424], 298.,  [512.58, 647.29],[0.000117, 5.6e-05], [0.559,0.344] )
    assert_allclose(Vl, 2.706588773271354e-05)

#    with pytest.raises(Exception):
#        COSTALD_mixture([0.4576, 0.5424], 298.,  [512.58],[0.000117, 5.6e-05], [0.559,0.344] )
