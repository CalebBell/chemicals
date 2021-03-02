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
from chemicals.lennard_jones import *
from chemicals.lennard_jones import LJ_data_Magalhaes
from fluids.constants import k


def test_LJ_data():
    # Two instances of 174899-66-2 were present;
    # the apparently more common one, [Bmim][CF 3SO 3], was kept.
    tot = LJ_data_Magalhaes['epsilon'].abs().sum()
    assert_close(tot, 187099.82029999999)

    tot = LJ_data_Magalhaes['sigma'].abs().sum()
    assert_close(tot, 1995.8174799999997)

    assert LJ_data_Magalhaes.index.is_unique
    assert LJ_data_Magalhaes.shape == (322, 3)


def test_molecular_diameter_CSP():
    # Example from DIPPR 1983 Manual for Predicting Chemical Process Design
    # Data, Data Prediction Manual, page 10-B2, Procedure 10B.
    # Example is shown only for Monofluorobenzene and hydrogen.
    # The Monofluorobenzene example is adapted here to all methods.

    # All functions have had their units and formulas checked and are relatively
    # similar in their outputs.

    sigma = sigma_Flynn(0.000268)
    assert_close(sigma, 5.2506948422196285)

    sigma = sigma_Bird_Stewart_Lightfoot_critical_2(560.1, 4550000)
    assert_close(sigma, 5.658657684653222)

    sigma = sigma_Bird_Stewart_Lightfoot_critical_1(0.000268)
    assert_close(sigma, 5.422184116631474)

    sigma = sigma_Bird_Stewart_Lightfoot_boiling(0.0001015)
    assert_close(sigma, 5.439018856944655)

    sigma = sigma_Bird_Stewart_Lightfoot_melting(8.8e-05)
    assert_close(sigma, 5.435407341351406)

    sigma = sigma_Stiel_Thodos(0.271E-3, 0.265)
    assert_close(sigma, 5.94300853971033)

    sigma = sigma_Tee_Gotoh_Steward_1(560.1, 4550000)
    assert_close(sigma, 5.48402779790962)

    sigma = sigma_Tee_Gotoh_Steward_2(560.1, 4550000, 0.245)
    assert_close(sigma, 5.412104867264477)

    sigma = sigma_Silva_Liu_Macedo(560.1, 4550000)
    assert_close(sigma, 5.164483998730177)
    assert sigma_Silva_Liu_Macedo(1084, 3.84E5) is None

@pytest.mark.slow
@pytest.mark.fuzz
def test_stockmayer_all_lookups():
    # Use the default method for each chemical in this file
    Stockmayers = sum([Stockmayer(CASRN=i) for i in LJ_data_Magalhaes.index])
    assert_close(Stockmayers, 187099.82029999999)

def test_stockmayer_function():

    assert_close(1291.41, Stockmayer(CASRN='64-17-5'))

    methods = Stockmayer_methods(Tm=178.075, Tb=341.87, Tc=507.6, Zc=0.2638, omega=0.2975, CASRN='110-54-3')
    assert tuple(methods) == Stockmayer_all_methods

    values_calc = [Stockmayer(Tm=178.075, Tb=341.87, Tc=507.6, Zc=0.2638, omega=0.2975, CASRN='110-54-3', method=i) for i in methods]
    values = [434.76, 427.33156230000003, 273.54201582027196, 318.10801442820025, 390.85200000000003, 392.8824, 393.15049999999997, 341.90399999999994]
    assert_close1d(values_calc, values)

    # Error handling
    assert None == Stockmayer(CASRN='BADCAS')

    with pytest.raises(Exception):
        Stockmayer(CASRN='98-01-1', method='BADMETHOD')

@pytest.mark.slow
@pytest.mark.fuzz
def test_molecular_diameter_all_values():
    # Use the default method for each chemical in this file
    MDs_sum = sum([molecular_diameter(CASRN=i) for i in LJ_data_Magalhaes.index])
    assert_close(MDs_sum, 1995.8174799999997)

def test_molecular_diameter_function():

    assert_close(4.23738, molecular_diameter(CASRN='64-17-5'))

    methods = molecular_diameter_methods(Tc=507.6, Pc=3025000.0, Vc=0.000368, Zc=0.2638, omega=0.2975, Vm=0.000113, Vb=0.000140, CASRN='110-54-3')
    assert tuple(sorted(methods)) == tuple(sorted(molecular_diameter_all_methods))

    values_calc = [molecular_diameter(Tc=507.6, Pc=3025000.0, Vc=0.000368, Zc=0.2638, omega=0.2975, Vm=0.000113, Vb=0.000140, CASRN='110-54-3', method=i) for i in methods]
    values = [5.61841, 5.989061939666203, 5.688003783388763, 6.27423491655056, 6.080607912773406, 6.617051217297049, 5.960764840627408, 6.0266865190488215, 6.054448122758386, 5.9078666913304225]
    assert_close1d(values_calc, values)

    # Error handling
    assert None == molecular_diameter(CASRN='BADCAS')

    with pytest.raises(Exception):
        molecular_diameter(CASRN='98-01-1', method='BADMETHOD')


def test_stockmayer():
    epsilon_k = epsilon_Flynn(560.1)
    assert_close(epsilon_k, 345.2984087011443)

    epsilon_k = epsilon_Bird_Stewart_Lightfoot_critical(560.1)
    assert_close(epsilon_k, 431.27700000000004)

    epsilon_k = epsilon_Bird_Stewart_Lightfoot_boiling(357.85)
    assert_close(epsilon_k, 411.5275)

    epsilon_k = epsilon_Bird_Stewart_Lightfoot_melting(231.15)
    assert_close(epsilon_k, 443.808)

    epsilon_k = epsilon_Stiel_Thodos(358.5, 0.265)
    assert_close(epsilon_k, 196.3755830305783)

    epsilon_k = epsilon_Tee_Gotoh_Steward_1(560.1)
    assert_close(epsilon_k, 433.5174)

    epsilon_k = epsilon_Tee_Gotoh_Steward_2(560.1, 0.245)
    assert_close(epsilon_k, 466.55125785)


def test_collision_integral():
    lss = [(1,1), (1,2), (1,3), (2,2), (2,3), (2,4), (2,5), (2,6), (4,4)]

    omega = collision_integral_Neufeld_Janzen_Aziz(100, 1, 1)
    assert_close(omega, 0.516717697672334)

    # Calculated points for each (l,s) at 0.3 and 100, nothing to use for comparison

    omegas_03_calc = [collision_integral_Neufeld_Janzen_Aziz(0.3, l, s) for (l,s) in lss]
    omegas_03 = [2.6501763610977873, 2.2569797090694452, 1.9653064255878911, 2.8455432577357387, 2.583328280612125, 2.3647890577440056, 2.171747217180612, 2.001780269583615, 2.5720716122665257]
    assert_close1d(omegas_03_calc, omegas_03)

    with pytest.raises(Exception):
        collision_integral_Neufeld_Janzen_Aziz(0.3, l=8, s=22)

    # Example 1.6 p 21 in Benitez, Jamie - 2009 - Principles and Modern Applications of Mass Transfer Operations
    omega = collision_integral_Neufeld_Janzen_Aziz(1.283)
    assert_close(omega, 1.282, atol=0.0001)


    # More accurate formulation
    omega = collision_integral_Kim_Monroe(400, 1, 1)
    assert_close(omega, 0.4141818082392228)

    # Series of points listed in [1]_ for values of l and s for comparison with
    # the  Hirschfelder et al correlation
    omegas_03_calc = [collision_integral_Kim_Monroe(0.3, l, s) for (l,s) in lss]
    omegas_03 = [2.65, 2.2568, 1.9665, 2.8436, 2.5806, 2.3623, 2.1704, 2.0011, 2.571]
    assert_close1d(omegas_03_calc, omegas_03, rtol=1e-4)

    omegas_400_calc = [collision_integral_Kim_Monroe(400, l, s) for (l,s) in lss]
    omegas_400 = [0.41418, 0.3919, 0.37599, 0.47103, 0.45228, 0.43778, 0.42604, 0.41622, 0.4589]
    assert_close1d(omegas_400_calc, omegas_400, rtol=1e-4)

    with pytest.raises(Exception):
        collision_integral_Kim_Monroe(0.3, l=8, s=22)

    all_keys = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 3), (3, 4), (3, 5), (4, 4)]
    for i in range(len(all_keys)):
        collision_integral_Neufeld_Janzen_Aziz(2.125, l=all_keys[i][0], s=all_keys[i][1])

    all_keys = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 3), (3, 4), (3, 5), (4, 4)]
    omegas = []
    omegas_expect = [1.0533593020147094, 0.9364963413584679, 0.8755022331873744, 0.839027248432364, 0.8122702217656526, 0.7919953217020422, 0.7756030101007225, 1.1517125076228243, 1.0555472726724568, 0.9970058623517577, 0.9584301480609076, 0.9308189354078076, 0.9996578276482491, 0.9450853163684426, 0.9076763140230476, 1.0544584236379315]
    for i in range(len(all_keys)):
        omega = collision_integral_Neufeld_Janzen_Aziz(2.125, l=all_keys[i][0], s=all_keys[i][1])
        omegas.append(omega)
    assert_close1d(omegas, omegas_expect, rtol=1e-12)


def test_T_star():
    const = 308.43
    Tst1 = T_star(T=318.2, epsilon_k=const)
    Tst2 = T_star(T=318.2, epsilon=const*k)
    assert_close1d([Tst1, Tst2], [1.0316765554582887]*2)

    with pytest.raises(Exception):
        T_star(300)
