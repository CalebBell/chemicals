# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016-2020, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
from numpy.testing import assert_allclose
import pytest
from chemicals.identifiers import *
from chemicals.elements import periodic_table, nested_formula_parser, serialize_formula, molecular_weight
import os
from chemicals.identifiers import ChemicalMetadataDB, folder, pubchem_db
from chemicals.identifiers import common_mixtures

# Force the whole db to load
try:
    CAS_from_any('asdfadsfasdfasdf')
except:
    pass

def test_dippr_list():
    dippr_set = dippr_compounds()
    # TODO CASs formulas
    assert 12916928773 == sum([CAS_to_int(i) for i in dippr_set])
    assert all([check_CAS(i) for i in dippr_set])


@pytest.mark.slow
@pytest.mark.online
def test_dippr_2016_matched_meta():
    df2 = pd.read_excel('https://www.aiche.org/sites/default/files/docs/pages/dippr_compound_list_2016.xlsx')
    names = df2['Name'].tolist()
    CASs = df2['CASN'].tolist()
    for i, CAS in enumerate(CASs):
        if CAS in ['16462-44-5', '75899-69-3']:
            # CELLOBIOSE (not the latest CAS, now is 528-50-7)
            # TRIPROPYLENE GLYCOL MONOETHYL ETHER (cannot find structure)
            continue
        assert CAS_from_any(CAS) == CAS

    # TODO names?


@pytest.mark.slow
def test_Matthews_critical_names():
    from chemicals.critical import critical_data_Matthews
    for CAS, name in zip(critical_data_Matthews.index, critical_data_Matthews['Chemical']):
        assert CAS_from_any(CAS) == CAS
    # might be worth doing
#        try:
#            assert CAS_from_any(name) == CAS
#        except:
#            print(CAS, name)

# TODO uncomment when electrochem is back
#@pytest.mark.slow
#def test_Laliberte_metadata_identifiers():
#    from thermo.electrochem import _Laliberte_Density_ParametersDict, _Laliberte_Viscosity_ParametersDict, _Laliberte_Heat_Capacity_ParametersDict
#    lalib = _Laliberte_Density_ParametersDict.copy()
#    lalib.update(_Laliberte_Viscosity_ParametersDict)
#    lalib.update(_Laliberte_Heat_Capacity_ParametersDict)
#    
#    for CAS, d in lalib.items():
#        c = None
#        formula = d['Formula']
#        name = d['Name']
#        if formula not in set(['HCHO2', 'CH3CH2OH', 'HCH3CO2']):
#            assert CAS_from_any(formula) == CAS
##        try:
##            CAS_from_any(name)
##        except:
##            print(name)
#

@pytest.mark.slow
def test_pubchem_dict():
    assert all([check_CAS(i.CASs) for i in pubchem_db.CAS_index.values()])

@pytest.mark.xfail
def test_database_formulas():
    # Failures are thing slike 3He, C2D4Br2, C14H18N3NaO10[99Tc], [1H]I
    # The fix here is adding an isotope db and making the formula parser handle isotopes as well.
    # This worked until isotopes were added to formulas
    assert all([i.formula == serialize_formula(i.formula) for i in pubchem_db.CAS_index.values()])

def test_organic_user_db():
    db = ChemicalMetadataDB(elements=False,
                            main_db=None,
                            user_dbs=[os.path.join(folder, 'chemical identifiers example user db.tsv')])
    for CAS, d in  db.CAS_index.items():
        assert CAS_from_any(d.CASs) == d.CASs
    # Check something was loaded
    assert len(db.CAS_index) > 100

    # Check smiles are unique / can lookup by smiles
    for smi, d in db.smiles_index.items():
        if not smi:
            continue
        assert CAS_from_any('smiles=' + smi) == d.CASs

    # Check formula is formatted right
    assert all([i.formula == serialize_formula(i.formula) for i in db.CAS_index.values()])

    # Check CAS validity
    assert all([check_CAS(i.CASs) for i in db.CAS_index.values()])

    # MW checker
    for i in db.CAS_index.values():
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_allclose(mw_calc, i.MW, atol=0.05)


    for CAS, d in db.CAS_index.items():
        assert CAS_from_any('InChI=1S/' + d.InChI) == int_to_CAS(CAS)
        
    for CAS, d in db.CAS_index.items():
        assert CAS_from_any('InChIKey=' + d.InChI_key) == int_to_CAS(CAS)

    # Test the pubchem ids which aren't -1
    for CAS, d in db.CAS_index.items():
        if d.pubchemid != -1:
            assert CAS_from_any('PubChem=' + str(d.pubchemid)) == int_to_CAS(CAS)

    CAS_lenth = len(db.CAS_index)
    assert CAS_lenth == len(db.smiles_index)
    assert CAS_lenth == len(db.InChI_index)
    assert CAS_lenth == len(db.InChI_key_index)


def test_inorganic_db():
    db = ChemicalMetadataDB(elements=False,
                            main_db=None,
                            user_dbs=[os.path.join(folder, 'Inorganic db.tsv')])

    # Check CAS lookup
    for CAS, d in  db.CAS_index.items():
        assert CAS_from_any(d.CASs) == d.CASs

    # Try ro check formula lookups
    for formula, d in  db.formula_index.items():
        if formula in set(['H2MgO2', 'F2N2']):
            # Formulas which are not unique by design
            continue
        assert CAS_from_any(formula) == d.CASs
    
    # Check smiles are unique / can lookup by smiles
    for smi, d in db.smiles_index.items():
        if not smi:
            continue
        assert CAS_from_any('smiles=' + smi) == d.CASs

    # Check formula is formatted right
    assert all([i.formula == serialize_formula(i.formula) for i in db.CAS_index.values()])

    # Check CAS validity
    assert all([check_CAS(i.CASs) for i in db.CAS_index.values()])

    # MW checker
    for i in db.CAS_index.values():
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_allclose(mw_calc, i.MW, atol=0.05)
    

def test_mixture_from_any():
    with pytest.raises(Exception):
        mixture_from_any(['water', 'methanol'])
    with pytest.raises(Exception):
        mixture_from_any('NOTAMIXTURE')
        
    for name in ['Air', 'air', u'Air', ['air']]:
        assert mixture_from_any(name) == common_mixtures['Air']

    names = ['R-401A ', ' R-401A ', 'R401A ', 'r401a', 'r-401A', 'refrigerant-401A', 'refrigerant 401A']
    for name in names:
        assert mixture_from_any(name) == common_mixtures['R401A']
        
    assert mixture_from_any('R512A') == common_mixtures['R512A']
    assert mixture_from_any([u'air']) == common_mixtures['Air']

def test_IDs_to_CASs():
    expect = ['811-97-2', '75-37-6']
    assert IDs_to_CASs('R512A') == expect
    assert IDs_to_CASs(['R512A']) == expect
    assert IDs_to_CASs(['norflurane', '1,1-difluoroethane']) == expect
    
    assert IDs_to_CASs(['norflurane']) == ['811-97-2']
    assert IDs_to_CASs('norflurane') == ['811-97-2']
    
def test_search_chemical():
    hit0 = search_chemical('water') 
    hit1 = search_chemical('water') 
    assert hit0 is hit1
    
    with pytest.raises(ValueError):
        # Not that smart/weird
        search_chemical('(oxidane)')
        
    assert search_chemical('water').charge == 0
    
    
def test_CAS_from_any():
    assert CAS_from_any('7732-18-5 ') == '7732-18-5'
    assert CAS_from_any('   7732  -18-5 ') == '7732-18-5'
    # direct in dictionary case
    assert CAS_from_any('water') == '7732-18-5'

    # Not in the main dict, but a synonym cas exists in the database
    assert CAS_from_any('136-16-3') == '582-36-5'

    assert CAS_from_any('inchi=1s/C2H6O/c1-2-3/h3H,2H2,1H3') == '64-17-5'
    assert CAS_from_any('inchi=1/C2H6O/c1-2-3/h3H,2H2,1H3') == '64-17-5'
    assert CAS_from_any('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3') == '64-17-5'
    assert CAS_from_any('InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3') == '64-17-5'
    assert CAS_from_any('InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N') == '64-17-5'
    assert CAS_from_any('inchikey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N') == '64-17-5'
    assert CAS_from_any(' inchikey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N') == '64-17-5'

    assert CAS_from_any('InChI=1S/C6H15N/c1-5-6(2)7(3)4/h6H,5H2,1-4H3') == '921-04-0'
    
    assert CAS_from_any('pubchem=702') == '64-17-5'
    
    assert CAS_from_any('oxidane') == '7732-18-5'
    
    assert CAS_from_any('CCCCCCCCCC') == '124-18-5'
    assert CAS_from_any('SMILES=CCCCCCCCCC') == '124-18-5'

    assert CAS_from_any('S') == '7704-34-9'
    assert CAS_from_any('O') == '17778-80-2'

    assert CAS_from_any('InChiKey=QVGXLLKOCUKJST-UHFFFAOYSA-N') == '17778-80-2'
    # Just because it's an element does not mean the CAS number refers to the 
    # monatomic form unfortunately - this is the CAS for Monooxygen
    
    assert CAS_from_any('1') == '12385-13-6'
    

    assert CAS_from_any('HC2O4-') == '920-52-5'
    
    assert CAS_from_any('water (H2O)') == '7732-18-5'
    
    # Test charge interpretation
    assert CAS_from_any('Ca+2') == '14127-61-8'
    assert CAS_from_any('Ca++') == '14127-61-8'
    assert CAS_from_any('Ca(++)') == '14127-61-8'
    assert CAS_from_any('Ca(+2)') == '14127-61-8'
    assert CAS_from_any('Ca(2+)') == '14127-61-8'

    # Unknown inchi
    with pytest.raises(Exception):
        CAS_from_any('InChI=1S/C13H14N2O2S/c18-13-15-14-12(17-13)8-16-11-6-5-9-3-1-2-4-10(9)7-11/h5-7H,1-4,8H2,(H,15,18)')
    with pytest.raises(Exception):
        CAS_from_any('InChI=1/C13H14N2O2S/c18-13-15-14-12(17-13)8-16-11-6-5-9-3-1-2-4-10(9)7-11/h5-7H,1-4,8H2,(H,15,18)')
    with pytest.raises(Exception):
        CAS_from_any('InChIKey=QHHWJJGJTYTIPM-UHFFFAOYSA-N')
    # unknown pubchem
    with pytest.raises(Exception):
        CAS_from_any('pubchem=902100')
    # unknown CAS
    with pytest.raises(Exception):
        CAS_from_any('1411769-41-9')
        
    with pytest.raises(Exception):
        # This was parsed as Cerium for a little while
        CAS_from_any('Cellulose')
        
        
def test_periodic_table_variants():
    """Do a lookup in the periodic table and compare vs CAS_from_any."""
    ids = [periodic_table._CAS_to_elements, periodic_table._name_to_elements, periodic_table._symbol_to_elements]
    failed_CASs = []
    for thing in ids:
        for i in thing.keys():
            try:
                CAS_from_any(i)
            except:
                failed_CASs.append(periodic_table[i].name)
    assert 0 == len(set(failed_CASs))
    
    # Check only the 5 known diatomics have a diff case
    failed_CASs = []
    for thing in ids:
        for i in thing.keys():
            try:
                assert CAS_from_any(i) == periodic_table[i].CAS
            except:
                failed_CASs.append(periodic_table[i].name)
    assert set(['Chlorine', 'Fluorine', 'Hydrogen', 'Nitrogen', 'Oxygen']) == set(failed_CASs)     
    
    
    for CAS, d in periodic_table._CAS_to_elements.items():
        assert CAS_from_any(d.smiles) == CAS
        
    for CAS, d in periodic_table._CAS_to_elements.items():
        assert CAS_from_any('SMILES=' + d.smiles) == CAS
        
    for CAS, d in periodic_table._CAS_to_elements.items():
        assert CAS_from_any('InChI=1S/' + d.InChI) == CAS
        
    for CAS, d in periodic_table._CAS_to_elements.items():
        assert CAS_from_any('InChIKey=' + d.InChI_key) == CAS
        
        
    fail = 0
    for CAS, d in periodic_table._CAS_to_elements.items():
        
        if d.PubChem != None:
            assert CAS_from_any('PubChem=' + str(d.PubChem)) == CAS
        else:
            fail += 1
    assert fail == 9
    # 111 - 118 aren't in pubchem
    
    
def test_fake_CAS_numbers():
    """File generated with :

    known = []
    for i in reversed(range(100000)):
        s = "20{0:0>5}000-00-0".format(i)
        if check_CAS(s):
            known.append(s+'\t\n')
    f = open('Fake CAS Registry.tsv', 'w')
    f.writelines(known)
    f.close()
    """
    # TODO
    

@pytest.mark.slow
def test_db_vs_ChemSep():
    """The CAS numbers are checked, as are most of the chemical formulas. Some
    chemical structural formulas aren't supported by the current formula parser
    and are ignored; otherwise it is a very effective test.

    DO NOT TRY TO OPTimizE THis FUNCTION - IT HAS ALREADY BEEN TRIED AND
    FAILED AT. THE TIME IS ONLY TAKEN py the PARSE function.

    EVEN THAT HAS BEEN REDUCED By 80% by using cElementTree instead of
    ElementTree.
    """
    
    import xml.etree.cElementTree as ET
    folder = os.path.join(os.path.dirname(__file__), 'Data')

    tree = ET.parse(os.path.join(folder, 'chemsep1.xml'))
    root = tree.getroot()

    data = {}
    for child in root:
        CAS, name, smiles, formula = None, None, None, None
        for i in child:
            tag = i.tag
            if CAS is None and tag == 'CAS':
                CAS = i.attrib['value']
            elif name is None and tag == 'CompoundID':
                name = i.attrib['value']
            elif smiles is None and tag == 'Smiles':
                smiles = i.attrib['value']
            elif formula is None and tag == 'StructureFormula':
                formula = i.attrib['value']
        
#        CAS = [i.attrib['value'] if  ][0]
#        name = [i.attrib['value'] for i in child if i.tag ][0]
#        smiles = [i.attrib['value'] for i in child if i.tag == ]
#        formula = [i.attrib['value'] for i in child if i.tag == 'StructureFormula'][0]
        
        try:
            if '-' in formula:
                formula = None
            else:
                formula = serialize_formula(formula)
        except:
            pass
        if smiles:
            smiles = smiles[0]
        else:
            smiles = None
        
        data[CAS] = {'name': name, 'smiles': smiles, 'formula': formula}        
    
    for CAS, d in data.items():
        hit = pubchem_db.search_CAS(CAS)
        assert hit.CASs == CAS

    for CAS, d in data.items():
        assert CAS_from_any(CAS) == CAS

    for CAS, d in data.items():
        f = d['formula']
        if f is None or f == '1,4-COOH(C6H4)COOH' or d['name'] == 'Air':
            continue
        assert pubchem_db.search_CAS(CAS).formula == f

    # In an ideal world, the names would match too but ~22 don't. Adding more synonyms
    # might help.
    # Some of them are straight disagreements however
#    for CAS, d in data.items():
#        try:
#            assert CAS_from_any(d['name']) == CAS
#        except:
#            print(CAS, d['name'])
##

    # In an ideal world we could also validate against their smiles
    # but that's proving difficult due to things like 1-hexene - 
    # is it 'CCCCC=C' or 'C=CCCCC'?
#test_db_vs_ChemSep() 
    
    


def test_CAS2int():
    assert CAS_to_int('7704-34-9') == 7704349

    with pytest.raises(Exception):
        CAS_to_int(7704349)

def test_int2CAS():
    assert int_to_CAS(7704349) == '7704-34-9'

    with pytest.raises(Exception):
        CAS_to_int(7704349.0)

def test_sorted_CAS_key():
    expect = ('64-17-5', '98-00-0', '108-88-3', '7732-18-5')
    res = sorted_CAS_key(['7732-18-5', '64-17-5', '108-88-3', '98-00-0'])
    assert res == expect
    res = sorted_CAS_key(['108-88-3', '98-00-0', '7732-18-5', '64-17-5'])
    assert res == expect
    
    invalid_CAS_expect = ('641', '98-00-0', '108-88-3', '7732-8-5')
    invalid_CAS_test = sorted_CAS_key(['7732-8-5', '641', '108-88-3', '98-00-0'])
    assert invalid_CAS_expect == invalid_CAS_test
    
