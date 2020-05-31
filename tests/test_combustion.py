# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.'''
import pytest
from chemicals.combustion import *
from fluids.numerics import assert_close, assert_close1d

def test_combustion():
    
    ### Lower level functions ###
    
    stoichiometry_calc = tuple(combustion_stoichiometry({'C': 1, 'H':4}).values())
    assert_close1d(stoichiometry_calc, [-2.0, 1.0, 2.0])
    
    # Methane gas
    H_calc = HHV_stoichiometry({'O2': -2.0, 'CO2': 1, 'H2O': 2.0}, -74520.0)
    assert_close(H_calc, -890604.0)
    
    # Dulong dry bituminous coal
    H_calc = HHV_modified_Dulong({'C': 0.716, 'H': 0.054, 'S': 0.016, 'N': 0.016, 'O': 0.093, 'Ash': 0.105})
    assert_close(H_calc, -304.0395)
    
    # LHV methanol
    H_calc = LHV_from_HHV(-726024.0, 2)
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
    H_calc = combustion_data({'H': 4, 'C': 1, 'O': 1}, Hf=-239100, method='Stoichiometry').HHV
    assert_close(H_calc, H_methanol)

    # Custom example of compound, to show all lines used
    H_custom = -7090590.5
    H_calc = combustion_data({'C': 10, 'H': 5, 'N': 3, 'O': 5, 'S': 2, 'Br': 8,
    'I':3, 'Cl': 3, 'F':2, 'P': 3}, Hf=-300000).HHV
    assert_close(H_calc, H_custom)

    with pytest.raises(ValueError):
        combustion_data({'H': 5, 'C': 6, 'O': 2, 'N': 1}, Hf=-344900, method='BADMETHOD')