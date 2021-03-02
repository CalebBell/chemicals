# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL).

Utilities for process modeling. Copyright (C) 2016, Caleb Bell
<Caleb.Andrew.Bell@gmail.com> Permission is hereby granted, free of charge, to
any person obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including without
limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions: The
above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS
IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import pytest
from chemicals.combustion import *
from fluids.numerics import assert_close, assert_close1d
from chemicals.combustion import fuel_air_third_spec_solver, is_combustible


def test_combustion_stoichiometry():
    # methane
    stoichiometry_calc = combustion_stoichiometry({'C': 1, 'H':4})
    assert len(stoichiometry_calc) == 3
    assert_close(stoichiometry_calc['O2'], -2.0)
    assert_close(stoichiometry_calc['CO2'], 1.0)
    assert_close(stoichiometry_calc['H2O'], 2.0)

    # buckyballs!
    C60_products = combustion_stoichiometry({'C': 60})
    assert_close(C60_products['CO2'], 60)
    assert_close(C60_products['O2'], -60)

    # test with its own oxygen included
    products = combustion_stoichiometry({'C': 4, 'O': 8})
    assert_close(products['CO2'], 4)
    # O2 coefficient is 0., so it shouldn't be included just like H2O
    with pytest.raises(KeyError): products['O2']

    products = combustion_stoichiometry({'N':2})
    assert_close(products['N2'], 1)
    with pytest.raises(KeyError): products['O2']

    # test combustion of everything except N2
    products = combustion_stoichiometry({'C':2.0, 'H': 3, 'P': 1, 'I': 2, 'Cl': 3, 'F': 2, 'Br': 1, 'S': 5, 'O': 1})
    assert_close(products['CO2'], 2)
    assert_close(products['Br2'], 0.5)
    assert_close(products['I2'], 1)
    assert_close(products['HCl'], 3)
    assert_close(products['HF'], 2)
    assert_close(products['SO2'], 5)
    assert 'N2' not in products
    assert_close(products['P4O10'], 0.25)
    assert_close(products['H2O'], -1.0)
    assert_close(products['O2'], -7.25)

    # test missing elements
    # For many elements, it would be feasible to change the behavior to include additional reactions;
    # I.e. lithium burns to LiOH

    products = combustion_stoichiometry({'C':2, 'H': 3, 'Li': 3}, missing_handling='elemental')
    assert_close(products['Li'], 3)

    # test ash handling
    products = combustion_stoichiometry({'C':2, 'H': 3, 'Li': 3}, missing_handling='ash')
    assert_close(products['CO2'], 2.0)
    assert_close(products['H2O'], 1.5)
    assert_close(products['O2'], -2.75)
    assert_close(products['Ash'], 20.823)

    products = combustion_stoichiometry({'C':2, 'H': 3, 'Li': 3}, MW=47.86822, missing_handling='ash')
    assert_close(products['CO2'], 2.0)
    assert_close(products['H2O'], 1.5)
    assert_close(products['O2'], -2.75)
    assert_close(products['Ash'], 20.823)

    with pytest.raises(ValueError):
        combustion_stoichiometry({'C':2, 'H': 3, 'Li': 3}, MW=47.86822, missing_handling='BADMETHOD')

def test_HHV_stoichiometry():
    # Methane gas
    H_calc = HHV_stoichiometry({'O2': -2.0, 'CO2': 1.0, 'H2O': 2.0}, -74520.0)
    assert_close(H_calc, -890604.0)

    H_Nicotinic_acid = -2730506.5
    H_calc = HHV_stoichiometry(combustion_stoichiometry({'H': 5, 'C': 6, 'O': 2, 'N': 1}), Hf=-344900)
    assert_close(H_calc, H_Nicotinic_acid)

    H_methanol = -726024.0
    H_calc = HHV_stoichiometry(combustion_stoichiometry({'H': 4.0, 'C': 1, 'O': 1}), Hf=-239100)
    assert_close(H_calc, H_methanol)

    # Custom example of compound, to show all lines used
    H_custom = -7090590.5
    H_calc = HHV_stoichiometry(combustion_stoichiometry({'C': 10, 'H': 5.0, 'N': 3, 'O': 5, 'S': 2, 'Br': 8,
    'I':3, 'Cl': 3, 'F':2, 'P': 3}), Hf=-300000)
    assert_close(H_calc, H_custom)

    # Custom example, with different Hf for each compound
    H_custom = -7090575.75
    H_calc = HHV_stoichiometry(combustion_stoichiometry({'C': 10, 'H': 5, 'N': 3, 'O': 5, 'S': 2, 'Br': 8,
    'I':3, 'Cl': 3, 'F':2, 'P': 3}), Hf=-300000, Hf_chemicals={'H2O': -285824, 'CO2': -393473, 'SO2': -296801,
                                                'Br2': 30881, 'I2': 62416, 'HCl': -92172, 'HF': -272710,
                                                'P4O10': -3009941, 'O2': 0.0, 'N2': 0.0})
    assert_close(H_calc, H_custom)


def test_combustion():

    ### Lower level functions ###


    # Dulong dry bituminous coal
    H_calc = HHV_modified_Dulong({'C': 0.716, 'H': 0.054, 'S': 0.016, 'N': 0.016, 'O': 0.093, 'Ash': 0.105})
    assert_close(H_calc, -304.0395)

    # LHV methanol
    H_calc = LHV_from_HHV(-726024.0, 2.0)
    assert_close(H_calc, -638001.008)

    # Methanol given formula
    cd = combustion_data('CH3OH', Hf=-239100)
    assert_close(cd.HHV, -726024.0)

    # LHV of methanol
    assert_close(cd.LHV, -638001.008)

    ### Combustion data (high level function) ###

    # Methanol given atoms (selects Stoichiometric method given Hf)
    H_Nicotinic_acid = -2730506.5
    H_calc = combustion_data({'H': 5, 'C': 6, 'O': 2, 'N': 1}, Hf=-344900).HHV
    assert_close(H_calc, H_Nicotinic_acid)

    # Custome example to select Dulong method when no Hf is given
    H_custom = -3004.589626
    H_calc = combustion_data({'C': 0.716, 'H': 0.054, 'S': 0.016, 'N': 0.016, 'O': 0.01}).HHV
    assert_close(H_calc, H_custom)

    # Custom example to test Dulong method
    H_calc = combustion_data({'C': 0.716, 'H': 0.054, 'S': 0.016, 'N': 0.016, 'O': 0.01}, method='Dulong').HHV
    assert_close(H_calc, H_custom)

    # Example with Stoichiometric method selected
    H_methanol = -726024.0
    H_calc = combustion_data({'H': 4, 'C': 1, 'O': 1}, Hf=-239100.0, method='Stoichiometry').HHV
    assert_close(H_calc, H_methanol)

    # Custom example of compound, to show all lines used
    H_custom = -7090590.5
    H_calc = combustion_data({'C': 10.0, 'H': 5.0, 'N': 3, 'O': 5, 'S': 2, 'Br': 8,
    'I':3, 'Cl': 3, 'F':2, 'P': 3}, Hf=-300000.0).HHV
    assert_close(H_calc, H_custom)

    with pytest.raises(ValueError):
        combustion_data({'H': 5, 'C': 6, 'O': 2, 'N': 1}, Hf=-344900.0, method='BADMETHOD')

    with pytest.raises(ValueError):
        combustion_data(['H']) # Bad type

    with pytest.raises(ValueError):
        combustion_data('C6H12O6', method='Dulong') # too much oxygen

    with pytest.raises(ValueError):
        combustion_data('CH3OH', {'H': 4, 'C': 1, 'O': 1})

    with pytest.raises(ValueError):
        combustion_data()


def test_air_fuel_ratio_solver():
    Vm_air = 0.024936627188566596
    Vm_fuel = 0.024880983160354486
    MW_air = 28.850334
    MW_fuel = 17.86651
    n_fuel = 5.0
    n_air = 25.0

    strings = ['mole', 'mass', 'volume']
    ratios = [n_air/n_fuel, MW_air/MW_fuel*5.0, Vm_air/Vm_fuel*5.0]

    ans_expect = [n_air, n_fuel]
    ans_expect_full = ans_expect + ratios

    for ratio, s in zip(ratios, strings):
        for air_spec, fuel_spec in zip([None, n_air], [n_fuel, None]):
            ans = air_fuel_ratio_solver(ratio=ratio, Vm_air=Vm_air, Vm_fuel=Vm_fuel,
                                        MW_air=MW_air,  MW_fuel=MW_fuel, n_air=air_spec,
                                        n_fuel=fuel_spec, basis=s)
            assert_close1d(ans, ans_expect_full)

    assert is_combustible('1333-74-0', {'H':2})

def assert_comb_dict_equal(calc, expect):
    for k, v_expect in expect.items():
        v_calc = calc[k]
        assert_close(v_expect, v_calc)

def fuel_air_spec_solver_checker(inputs, ans, func=fuel_air_spec_solver):
    # One intensive variable, one extensive variable
    for n_spec, n_name in zip(['n_fuel', 'n_air', 'n_out'], [ans['n_fuel'], ans['n_air'], ans['n_out']]):
        for i_spec, i_name in zip(['O2_excess', 'frac_out_O2', 'frac_out_O2_dry'], [ans['O2_excess'], ans['frac_out_O2'], ans['frac_out_O2_dry']]):
            d2 = inputs.copy()
            d2.update({n_spec: n_name, i_spec: i_name})

            calc = func(**d2)
            assert_comb_dict_equal(calc, ans)


    for n_spec, n_name in zip(['n_fuel', 'n_air'], [ans['n_fuel'], ans['n_air']]):
        for ratio_name, ratio in zip(['mass', 'mole', 'volume'], [ans['mass_ratio'], ans['mole_ratio'], ans['volume_ratio']]):
            kwargs = inputs.copy()
            kwargs.update({n_spec: n_name, 'ratio': ratio, 'ratio_basis': ratio_name})

            calc = func(**kwargs)
            assert_comb_dict_equal(calc, ans)

#     Two extensive variables check
    for n_fuel in (None, ans['n_fuel']):
        for n_air in (None, ans['n_air']):
            for n_out in (None, ans['n_out']):
                if sum(i is not None for i in (n_fuel, n_air, n_out)) == 2:

                    calc = func(n_fuel=n_fuel, n_air=n_air, n_out=n_out, **inputs)
                    assert_comb_dict_equal(calc, ans)


def test_fuel_air_spec_solver():
    ans_N7_messy_fuel = {'O2_excess': 0.5813397129186602,
 'frac_out_O2': 0.0712025316455696,
 'frac_out_O2_dry': 0.08118672947779891,
 'mass_ratio': 24.83317358818291,
 'mole_ratio': 16.0,
 'n_air': 80.0,
 'n_fuel': 5.0,
 'n_out': 85.32000000000001,
 'ns_out': [63.325, 6.074999999999999, 0.0, 0.0, 0.0, 10.492500000000001, 5.4275],
 'volume_ratio': 16.032156677022826,
 'zs_out': [0.742205813408345, 0.0712025316455696, 0.0, 0.0, 0.0, 0.12297819971870605, 0.06361345522737927]}

    inputs_N7_messy_fuel = {'CASs': ['7727-37-9',  '7782-44-7',  '74-82-8',  '74-84-0',  '74-98-6',  '7732-18-5',  '124-38-9'],
     'MW_air': 28.793413510000008,
     'MW_fuel': 18.551580390000005,
     'Vm_air': 0.024932453821680217,
     'Vm_fuel': 0.024882445274415996,
     'zs_air': [.79, .205, 0.0, 0.0, 0.0, .0045, .0005],
     'atomss': [{'N': 2},  {'O': 2},  {'C': 1, 'H': 4},  {'C': 2, 'H': 6},  {'C': 3, 'H': 8},  {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
     'zs_fuel': [.025, .025, .85, .07, .029, .0005, .0005]}

    all_inputs = [inputs_N7_messy_fuel]
    all_ans = [ans_N7_messy_fuel]

    for inputs, ans in zip(all_inputs, all_ans):
        fuel_air_spec_solver_checker(inputs, ans)

def test_fuel_air_third_spec_solver():

    inputs_N7_messy_fuel = {'CASs': ['7727-37-9',  '7782-44-7',  '74-82-8',  '74-84-0',  '74-98-6',  '7732-18-5',  '124-38-9'],
         'MW_air': 28.793413510000008,
         'MW_fuel': 18.551580390000005,
         'MW_third': 22.594160765550242,
         'Vm_air': 0.024932453821680217,
         'Vm_fuel': 0.024882445274415996,
         'Vm_third': 0.024820446149354414,
         'n_third': 1.0,
         'zs_air': [.79, .205, 0.0, 0.0, 0.0, .0045, .0005],
         'zs_third': [0.1, 0.005, 0.5, 0.39, 0.0, 0.005, 0.0],
         'zs_fuel': [.025, .025, .85, .07, .029, .0005, .0005],
         'atomss': [{'N': 2},  {'O': 2},  {'C': 1, 'H': 4},  {'C': 2, 'H': 6},  {'C': 3, 'H': 8},  {'H': 2, 'O': 1}, {'C': 1, 'O': 2}],
    }

    ans_N7_messy_fuel = {'O2_excess': 0.289894654701522,
     'frac_out_O2': 0.04294053054383636,
     'frac_out_O2_dry': 0.0503063746233793,
     'mass_ratio': 19.969067102685425,
     'mole_ratio': 13.333333333333334,
     'n_air': 80.0,
     'n_fuel': 5.0,
     'n_out': 86.51500000000001,
     'ns_out': [63.42500000000001,  3.7150000000000034,  0.0,  0.0,  0.0,  12.667499999999999,  6.707499999999999],
     'volume_ratio': 13.365681067247072,
     'zs_out': [0.7331098653412703,  0.04294053054383636,  0.0,  0.0,  0.0, 0.14641969600647284,  0.07752990810842048]}

    inputs_N7_messy_fuel_reactivities = inputs_N7_messy_fuel.copy()
    inputs_N7_messy_fuel_reactivities['reactivities'] = [True, True, True, False, False, True, True]

    ans_N7_messy_fuel_reactivities = {'O2_excess': 0.74,
     'frac_out_O2': 0.08174418604651168,
     'frac_out_O2_dry': 0.09233901421863204,
     'mass_ratio': 19.969067102685425,
     'mole_ratio': 13.333333333333337,
     'n_air': 80.00000000000003,
     'n_fuel': 5.0,
     'n_out': 86.00000000000003,
     'ns_out': [63.425, 7.03, 0.0, 0.74, 0.145, 9.8675, 4.7925],
     'volume_ratio': 13.365681067247069,
     'zs_out': [0.7375000000000002, 0.08174418604651168, 0.0, 0.008604651162790695, 0.001686046511627906, 0.11473837209302319, 0.055726744186046485]}

    all_inputs = [inputs_N7_messy_fuel, inputs_N7_messy_fuel_reactivities]
    all_ans = [ans_N7_messy_fuel, ans_N7_messy_fuel_reactivities]

    for func in (fuel_air_third_spec_solver, combustion_spec_solver):
        for inputs, ans in zip(all_inputs, all_ans):
            fuel_air_spec_solver_checker(inputs, ans, func=func)
