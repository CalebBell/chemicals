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
from fluids.numerics import assert_close

from chemicals.elements import *
from chemicals.elements import periodic_table

def test_molecular_weight():
    MW_calc = molecular_weight({'H': 12, 'C': 20, 'O': 5})
    MW = 332.30628
    assert_allclose(MW_calc, MW)

    MW_calc = molecular_weight({'C': 32, 'Cu': 1, 'H': 12, 'Na': 4, 'S': 4, 'O': 12, 'N': 8})
    MW = 984.24916
    assert_allclose(MW_calc, MW)

    with pytest.raises(Exception):
        molecular_weight({'H': 12, 'C': 20, 'FAIL': 5})

    assert_close(molecular_weight({'H': 11, 'T': 1, 'C': 20, 'O': 5}), 334.3143892)
    assert_close(molecular_weight({'H': 11, 'D': 1, 'C': 20, 'O': 5}), 333.312442)


def test_mass_fractions():
    mfs_calc = mass_fractions({'H': 12, 'C': 20, 'O': 5})
    mfs = {'H': 0.03639798802478244, 'C': 0.7228692758981262, 'O': 0.24073273607709128}
    assert_allclose(sorted(mfs.values()), sorted(mfs_calc.values()))

    mfs_calc = mass_fractions({'C': 32, 'Cu': 1, 'H': 12, 'Na': 4, 'S': 4, 'O': 12, 'N': 8})
    mfs = {'C': 0.39049299264832493, 'H': 0.012288839545466314, 'O': 0.19506524140696244, 'N': 0.11384678245496294, 'S': 0.13031253183899086, 'Na': 0.09343069187887344, 'Cu': 0.06456292022641909}
    assert_allclose(sorted(mfs.values()), sorted(mfs_calc.values()))

    # Fail two tests, one without MW, and one with MW
    with pytest.raises(Exception):
        mass_fractions({'FAIL': 12, 'C': 20, 'O': 5})
    with pytest.raises(Exception):
        mass_fractions({'FAIL': 12, 'C': 20, 'O': 5}, MW=120)


def test_atom_fractions():
    fractions_calc = atom_fractions({'H': 12, 'C': 20, 'O': 5})
    fractions = {'H': 0.32432432432432434, 'C': 0.5405405405405406, 'O': 0.13513513513513514}
    assert_allclose(sorted(fractions_calc.values()), sorted(fractions.values()))


def test_similarity_variable():
    sim1 = similarity_variable({'H': 32, 'C': 15})
    sim2 = similarity_variable({'H': 32, 'C': 15}, 212.41458)
    assert_allclose([sim1, sim2], [0.2212654140784498]*2)


def test_elements_data():
    tots_calc = [sum([getattr(i, att) for i in periodic_table if not getattr(i, att) is None]) for att in
    ['number', 'MW', 'period', 'group', 'AReneg', 'rcov', 'rvdw', 'maxbonds', 'elneg', 'ionization', 'elaffinity', 'electrons', 'protons']]
    tots_exp = [7021, 17285.2137652, 620, 895, 109.91, 144.3100000000001, 179.4300000000001, 94, 163.27000000000007, 816.4238999999999, 67.50297235000001, 7021, 7021]
    assert_allclose(tots_calc, tots_exp)

def test_misc_elements():
    assert periodic_table['H'].InChI == 'H' # 'InChI=1S/
    assert periodic_table['H'].smiles == '[H]'
    assert periodic_table[1].smiles == '[H]'
    assert periodic_table['1'].smiles == '[H]'

    assert periodic_table['Fm'].InChI_key == 'MIORUQGGZCBUGO-UHFFFAOYSA-N'

    assert periodic_table.Na is periodic_table['Na']

    assert len(periodic_table) == 118

    with pytest.raises(AttributeError):
        periodic_table.adamantium

    with pytest.raises(KeyError):
        periodic_table['BadElement']

    assert periodic_table.H.CAS_standard != periodic_table.H.CAS

    assert not 'BadElement' in periodic_table


    # Found some disagreement between sources about these - test the borders and
    # a few random ones to write expected results in the tests as well as the code.
    assert periodic_table.H.block == 's'
    assert periodic_table.He.block == 's'
    assert periodic_table.Fr.block == 's'

    assert periodic_table.Sc.block == 'd'
    assert periodic_table.Zn.block == 'd'
    assert periodic_table.Cn.block == 'd'

    assert periodic_table.Ce.block == 'f'
    assert periodic_table.Lu.block == 'f'
    assert periodic_table.Th.block == 'f'
    assert periodic_table.Lr.block == 'f'

    # Some categorize this as `d`
    assert periodic_table.La.block == 'f'
    assert periodic_table.Ac.block == 'f'

    d_block_border_elements = ['Sc', 'Y', 'Ti', 'Zr','Hf', 'Rf', 'Zn', 'Cd', 'Hg', 'Cn']
    for ele in d_block_border_elements:
        assert periodic_table[ele].block == 'd'

    p_block_border_elements = ['B', 'Ne', 'Al', 'Ar', 'Ga', 'Kr', 'In', 'Xe', 'Tl', 'Rn', 'Nh', 'Og']
    for ele in p_block_border_elements:
        assert periodic_table[ele].block == 'p'

    for i in periodic_table:
        str(i)
        i.__repr__()

    # Check that the monatomic elements all have a standard CAS number available
    assert periodic_table.Br.CAS_standard == '7726-95-6' # Br2
    assert periodic_table.Br.CAS == '10097-32-2' # monatomic

    assert periodic_table.I.CAS == '14362-44-8' # An earlier bug had 20461-54-5 which is the ionic form, not the monatomic form
    assert periodic_table.I.CAS_standard == '7553-56-2'

    assert (periodic_table.O.CAS, periodic_table.O.CAS_standard) == ('17778-80-2', '7782-44-7')
    assert (periodic_table.H.CAS, periodic_table.H.CAS_standard) == ('12385-13-6', '1333-74-0')
    assert (periodic_table.N.CAS, periodic_table.N.CAS_standard) == ('17778-88-0', '7727-37-9')
    assert (periodic_table.F.CAS, periodic_table.F.CAS_standard) == ('14762-94-8', '7782-41-4')
    assert (periodic_table.Cl.CAS, periodic_table.Cl.CAS_standard) == ('22537-15-1', '7782-50-5')

def test_Hill_formula():
    Hill_formulas = {'ClNa': {'Na': 1, 'Cl': 1}, 'BrI': {'I': 1, 'Br': 1},
                    'CCl4': {'C': 1, 'Cl': 4}, 'CH3I': {'I': 1, 'H': 3, 'C': 1},
                    'C2H5Br': {'H': 5, 'C': 2, 'Br': 1}, 'H2O4S': {'H': 2, 'S': 1, 'O': 4}}

    for formula, atoms in Hill_formulas.items():
        assert formula == atoms_to_Hill(atoms)


def test_simple_formula_parser():
    formulas = ['CO2',
                'H20OCo2NaClH4P4',
                'C.234O2',
                'O2C.234',
                'C555.234O.0000000000000000000000000000000001']
    results = [{'C': 1, 'O': 2},
               {'P': 4, 'Co': 2, 'Cl': 1, 'H': 24, 'Na': 1, 'O': 1},
               {'O': 2, 'C': 0.234},
               {'O': 2, 'C': 0.234},
               {'C': 555.234, 'O': 1e-34}
               ]

    for f in [simple_formula_parser, nested_formula_parser]:
        for formula, result in zip(formulas, results):
            assert f(formula) == result

def test_nested_formula_parser():
    with pytest.raises(ValueError):
        nested_formula_parser('Adamantium(NH3)4.0001+2')

    # repeat elements
    res = nested_formula_parser('Pd(NH3)4.0001Na(NH3)2+2')
    assert res == {'Pd': 1, 'N': 6.0001, 'H': 18.0003, 'Na': 1}

    res = nested_formula_parser('Pd(NH3)4.0001Na(NH3H2)2+2')
    assert res == {'Pd': 1, 'N': 6.0001, 'H': 22.0003, 'Na': 1}

def test_charge_from_formula():
    assert charge_from_formula('Br3-') == -1
    assert charge_from_formula('Br3-1') == -1
    assert charge_from_formula('Br3-2') == -2
    assert charge_from_formula('Br3-3') == -3
    assert charge_from_formula('Br3+') == 1
    assert charge_from_formula('Br3+1') == 1
    assert charge_from_formula('Br3+2') == 2
    assert charge_from_formula('Br3+3') == 3
    assert charge_from_formula('Br3') == 0
    assert charge_from_formula('Br3--') == -2
    assert charge_from_formula('Br3(--)') == -2
    assert charge_from_formula('Br3++') == 2
    assert charge_from_formula('Br3(++)') == 2

    assert charge_from_formula('Br3(-)') == -1
    assert charge_from_formula('Br3(-1)') == -1
    assert charge_from_formula('Br3(-2)') == -2
    assert charge_from_formula('Br3(-3)') == -3
    assert charge_from_formula('Br3(+)') == 1
    assert charge_from_formula('Br3(+1)') == 1
    assert charge_from_formula('Br3(+2)') == 2
    assert charge_from_formula('Br3(+3)') == 3

    with pytest.raises(ValueError):
        charge_from_formula('Br3(-+)')


def test_serialize_formula():
    assert serialize_formula('Pd(NH3)4+3') == 'H12N4Pd+3'
    assert 'H12N4Pd' == serialize_formula('Pd(NH3)4+0')
    assert 'H12N4Pd+' == serialize_formula('Pd(NH3)4+1')
    assert 'H12N4Pd-' == serialize_formula('Pd(NH3)4-1')
    assert 'H12N4Pd-5' == serialize_formula('Pd(NH3)4-5')


def test_mixture_atomic_composition_ordered():
    ns, names = mixture_atomic_composition_ordered([{'O': 2}, {'N': 1, 'O': 2}, {'C': 1, 'H': 4}], [0.95, 0.025, .025])
    assert names == ['H', 'C', 'N', 'O']
    assert_allclose(ns, [0.1, 0.025, 0.025, 1.95], rtol=1e-12)


def test_atom_matrix():
    atomss = [{'C': 1, 'H': 4}, {'C': 2, 'H': 6}, {'N': 2}, {'O': 2}, {'H': 2, 'O': 1}, {'C': 1, 'O': 2}]
    OCH_expect = [[0.0, 1, 4],
     [0.0, 2, 6],
     [0.0, 0.0, 0.0],
     [2, 0.0, 0.0],
     [1, 0.0, 2],
     [2, 1, 0.0]]

    default_expect = [[4, 1, 0.0, 0.0],
     [6, 2, 0.0, 0.0],
     [0.0, 0.0, 2, 0.0],
     [0.0, 0.0, 0.0, 2],
     [2, 0.0, 0.0, 1],
     [0.0, 1, 0.0, 2]]
    default = atom_matrix(atomss)

    assert_allclose(default, default_expect, rtol=1e-12)

    OCH = atom_matrix(atomss, ['O', 'C', 'H'])
    assert_allclose(OCH, OCH_expect, rtol=1e-12)