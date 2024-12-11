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

import os
from math import isnan

import pandas as pd
import pytest
import re
from fluids.numerics import assert_close
from chemicals.miscdata import heos_data
import json
from chemicals.heat_capacity import Cp_data_Poling
from chemicals.elements import molecular_weight, nested_formula_parser, periodic_table, serialize_formula
from chemicals.identifiers import (
    CAS_from_any,
    CAS_to_int,
    ChemicalMetadataDB,
    ChemicalMetadata,
    IDs_to_CASs,
    check_CAS,
    common_mixtures,
    dippr_compounds,
    folder,
    int_to_CAS,
    mixture_from_any,
    pubchem_db,
    search_chemical,
    sorted_CAS_key,
    FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS,
)

# Force the whole db to load
try:
    CAS_from_any('asdfadsfasdfasdf')
except:
    pass

def test_dippr_list():
    dippr_set = dippr_compounds()
    # TODO CASs formulas
    assert 12916928773 == sum([CAS_to_int(i) for i in dippr_set])
    for i in dippr_set:
        assert check_CAS(i)


def clean_dippr_spreadsheet_name(name):
    """
    Removes the _YEARReview pattern from chemical names.
    """
    pattern = r'_\d{4}Review$'
    return re.sub(pattern, '', name)

@pytest.mark.slow
@pytest.mark.online
def test_dippr_2016_matched_meta():
    df2 = pd.read_excel('https://www.aiche.org/sites/default/files/docs/pages/dippr_compound_list_2016.xlsx')
    names = df2['Name'].tolist()
    CASs = df2['CASN'].tolist()
    for i, CAS in enumerate(CASs):
        if CAS in ('16462-44-5', '75899-69-3', '132259-10-0', '55505-26-5'):
            # CELLOBIOSE (not the latest CAS, now is 528-50-7)
            # TRIPROPYLENE GLYCOL MONOETHYL ETHER (cannot find structure)
            # 55505-26-5 is now 25339-17-7
            continue
        try:
            if isnan(CAS):
                continue
        except:
            pass
        assert CAS_from_any(CAS) == CAS


    excluded_cas = {
        '16462-44-5',  # CELLOBIOSE
        '75899-69-3',  # TRIPROPYLENE GLYCOL MONOETHYL ETHER: can't find a structure
        '132259-10-0', # AIR: not a pure chemical
        '51177-04-9',  # COBALT CARBIDE: can't find a source
        '142-47-2',    # MONOSODIUM GLUTAMATE: https://commonchemistry.cas.org/detail?cas_rn=142-47-2&search=Monosodium%20glutamate I'm pretty sure
        '14762-55-1',  # HELIUM-3 : need isotope work
        '1314-80-3',   # PHOSPHORUS PENTASULFIDE: can't find a source
        '6303-21-5',   # HYPOPHOSPHOROUS ACID: can't find a source
        '26762-93-6',  # m-DIISOPROPYLBENZENE HYDROPEROXIDE: can't find a source
        '30025-38-8',  # DIPROPYLENE GLYCOL MONOETHYL ETHER: can't find a source
        '19295-81-9'   # 1,2-BENZENEDICARBOXYLIC ACID, HEPTYL, NONYL ESTER: can't find a source
    }
    
    for _, row in df2.iterrows():
        cas = row['CASN']
        if cas in excluded_cas or pd.isna(cas):
            continue
        formula = row['Formula']
        db_chem = search_chemical(cas)
        assert serialize_formula(formula) == serialize_formula(db_chem.formula)
    
    CASs_names_do_not_match = set(['2958-75-0', '1008-17-9', '103554-13-8', '3055-14-9', '74338-98-0', '163702-07-6', '355-37-3', '13465-77-5', '577-11-7', '14868-53-2', '70807-90-8', '10049-60-2', '371-78-8', '7331-52-4', '132259-10-0', '7723-14-0', '7784-30-7', '2141-58-4', '2141-59-5', '1634-09-9', '12075-68-2', '16747-50-5', '871-28-3', '2876-53-1', '1190-76-7', '2366-36-1', '6362-80-7', '768-00-3', '767-99-7', '66325-11-9', '6145-31-9', '2150-02-9', '10441-57-3', '5451-92-3', '6163-64-0', '13286-92-5', '628-87-5', '7327-60-8', '1333-74-0', '5756-43-4', '14290-92-7', '94023-15-1', '38433-80-6', '78448-33-6', '57-88-5', '1185-39-3', '21282-97-3', '1589-47-5', '132739-31-2', '2039-93-2', '7439-15-8', '1961-96-2', '694-92-8', '15403-89-1', '94-60-0', '26762-93-6', '4454-05-1', '616-02-4', '2432-74-8', '94-60-0', '13511-13-2', '78024-33-6', '1691-17-4', '37143-54-7', '77-68-9', '54839-24-6', '29911-27-1', '2465-32-9', '15798-64-8', '39972-78-6', '106-20-7', '75899-69-3', '42448-85-1', '19295-81-9', '65185-88-8', '19177-04-9', '106538-38-9', '112-47-0', '66032-51-7', '112-47-0', '17851-27-3', '51526-06-8', '51655-57-3', '1191-87-3', '1243297-10-0', '22663-61-2', '565-48-0', '58797-58-3', '32970-45-9', '33021-02-2', '4565-32-6', '6165-55-5', '14814-09-6', '3006-96-0', '108-59-8', '35112-74-4', '24615-84-7', '111-91-1', '95-96-5', '3228-03-3', '3228-02-2', '628-08-0', '6737-11-7', '3822-68-2', '16587-40-9', '544-02-5', '52458-04-5', '100524-60-5', '7784-34-1', '2471-08-1', '7351-61-3', '99-75-2', '1551-32-2', '17890-53-8', '64001-06-5', '15507-13-8', '60956-33-4', '99172-63-1', '837-08-1', '34885-03-5', '7403-22-7', '1196-81-2', '26158-99-6', '23305-64-8', '25961-89-1', '13349-10-5', '10410-35-2', '3603-45-0', '2530-83-8', '122-52-1', '2752-17-2', '871518-84-2', '857237-25-3', '15890-40-1', '1559-81-5', '13556-58-6', '31283-14-4', '1964-45-0', '13463-40-6', '148462-57-1', '79808-30-3', '21482-12-2', '132739-31-2',
                        '55505-26-5' # see previous comment
                        ])
    
    # Not investigated
    names_lead_to_different_CAS = {'13598-36-2', # https://commonchemistry.cas.org/detail?cas_rn=10294-56-1  https://commonchemistry.cas.org/detail?cas_rn=13598-36-2
                            '14808-60-7', # I used silica not quartz
                    '2687-91-4', '117-81-7', '1345-25-1', '7803-62-5', '118-93-4', '706-31-0', '16462-44-5', '3319-31-1', '873-66-5', '1344-28-1', '18328-90-0', '7726-95-6', '16219-75-3', '21460-36-6', '7722-76-1', '872-50-4', '4050-45-7'}

    for _, row in df2.iterrows():
        cas = row['CASN']
        if cas in CASs_names_do_not_match or pd.isna(cas):
            continue
        name = clean_dippr_spreadsheet_name(row['Name'])
        if not name:
            continue
        # check we find a chemical anyway TODO check the CAS matches
        db_chem = search_chemical(name)
        if cas not in names_lead_to_different_CAS:
            assert db_chem.CASs == cas, name

def test_silica():
    chemical = search_chemical('silica')
    # Test basic properties
    assert chemical.pubchemid == 24261
    assert chemical.CASs == '7631-86-9'
    assert chemical.formula == 'O2Si'
    assert_close(chemical.MW, 60.0843)
    assert chemical.smiles == 'O=[Si]=O'
    assert chemical.InChI == 'O2Si/c1-3-2'
    assert chemical.InChI_key == 'VYPSYNLAJGMNEJ-UHFFFAOYSA-N'
    assert chemical.common_name == 'silica'
    assert chemical.iupac_name == 'dioxosilane'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'silicon dioxide',
        'SiO2',
        'silicon(IV) oxide',
        'silicon oxide',
        'silicic anhydride',
        'silica, amorphous',
        'synthetic amorphous silica',
        'silicon dioxide (amorphous)'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7631-86-9', f"Failed to find silica using term: {term}"
    
    identifier_searches = {
        'formula': 'O2Si',
        'smiles': 'O=[Si]=O',
        'inchi': 'InChI=1S/O2Si/c1-3-2',
        'inchikey': 'InChIKey=VYPSYNLAJGMNEJ-UHFFFAOYSA-N',
        'pubchem': 'pubchem=24261',
        'cas': '7631-86-9'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7631-86-9', f"Failed to find silica using {search_type}: {identifier}"


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
#    from thermo.electrochem import rho_dict_Laliberte, mu_dict_Laliberte, Cp_dict_Laliberte
#    lalib = rho_dict_Laliberte.copy()
#    lalib.update(mu_dict_Laliberte)
#    lalib.update(Cp_dict_Laliberte)
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
    for i in pubchem_db:
        assert check_CAS(i.CASs)

@pytest.mark.skip
@pytest.mark.xfail
def test_database_formulas():
    # Failures are thing slike 3He, C2D4Br2, C14H18N3NaO10[99Tc], [1H]I
    # The fix here is adding an isotope db and making the formula parser handle isotopes as well.
    # This worked until isotopes were added to formulas
    for i in pubchem_db:
        assert i.formula == serialize_formula(i.formula)

def test_organic_user_db():
    db = ChemicalMetadataDB(elements=False,
                            main_db=None,
                            user_dbs=[os.path.join(folder, 'chemical identifiers example user db.tsv')])
    pref_file = os.path.join(folder, 'organic_preferences.json')
    if os.path.exists(pref_file):
        with open(pref_file) as f:
            preferences = json.load(f)
            preferred_CAS = set(preferences.get('preferred_cas', []))
            unpreferred_CAS = set(preferences.get('unpreferred_cas', []))


    for CAS, d in  db.CAS_index.items():
        assert CAS_from_any(d.CASs) == d.CASs
    # Check something was loaded
    assert len([i for i in db]) > 100

    # Check smiles are unique / can lookup by smiles
    # Check the smiles
    errors = []
    for CAS, d in db.CAS_index.items():
        smiles = d.smiles
        if smiles:
            result = CAS_from_any('smiles='+smiles)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by smiles {smiles} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check formula is formatted right
    for i in db:
        assert i.formula == serialize_formula(i.formula)

    # Check CAS validity
    for i in db:
        assert check_CAS(i.CASs)

    # MW checker
    for i in db:
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_close(mw_calc, i.MW, atol=0.05)


    # Check the inchi
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI = d.InChI
        if InChI:
            result = CAS_from_any('InChI=1S/'+InChI)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI {InChI} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the InChI_key
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI_key = d.InChI_key
        if InChI_key:
            result = CAS_from_any('InChIKey='+InChI_key)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI_key {InChI_key} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the pubchem
    errors = []
    for CAS, d in db.CAS_index.items():
        pubchemid = d.pubchemid
        if pubchemid != -1:
            result = CAS_from_any('PubChem='+str(pubchemid))
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by pubchemid {pubchemid} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    CAS_lenth = len([v for v in db])
    assert CAS_lenth == len(db.smiles_index)
    assert CAS_lenth == len(db.InChI_index)
    assert CAS_lenth == len(db.InChI_key_index)

    # Check the name
    # # TODO get this working
    # errors = []
    # for CAS, d in db.CAS_index.items():
    #     common_name = d.common_name
    #     result = CAS_from_any(common_name)
    #     if result != d.CASs:
    #         errors.append(f"CAS {result} retrieved by common_name {common_name} on compound {d.CASs}")
    # if errors:
    #     raise ValueError('\n'.join(errors))

def test_inorganic_db():
    db = ChemicalMetadataDB(elements=False,
                            main_db=None,
                            user_dbs=[os.path.join(folder, 'Inorganic db.tsv')])

    pref_file = os.path.join(folder, 'inorganic_preferences.json')
    if os.path.exists(pref_file):
        with open(pref_file) as f:
            preferences = json.load(f)
            preferred_CAS = set(preferences.get('preferred_cas', []))
            unpreferred_CAS = set(preferences.get('unpreferred_cas', []))
            
    # Check CAS lookup
    for CAS, d in db.CAS_index.items():
        assert CAS_from_any(d.CASs) == d.CASs

    # Check formula lookups
    errors = []
    for CAS, d in db.CAS_index.items():
        formula = d.formula
        result = CAS_from_any(formula)
        if result != d.CASs:
            if int_to_CAS(CAS) not in unpreferred_CAS:
                errors.append(f"CAS {result} retrieved by formula {formula} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the smiles
    errors = []
    for CAS, d in db.CAS_index.items():
        smiles = d.smiles
        if smiles:
            result = CAS_from_any('smiles='+smiles)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by smiles {smiles} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the inchi
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI = d.InChI
        if InChI:
            result = CAS_from_any('InChI=1S/'+InChI)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI {InChI} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the InChI_key
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI_key = d.InChI_key
        if InChI_key:
            result = CAS_from_any('InChIKey='+InChI_key)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI_key {InChI_key} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the pubchem
    errors = []
    for CAS, d in db.CAS_index.items():
        pubchemid = d.pubchemid
        if pubchemid != -1:
            result = CAS_from_any('PubChem='+str(pubchemid))
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by pubchemid {pubchemid} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the name
    errors = []
    for CAS, d in db.CAS_index.items():
        common_name = d.common_name
        result = CAS_from_any(common_name)
        if result != d.CASs:
            errors.append(f"CAS {result} retrieved by common_name {common_name} on compound {d.CASs}")
    if errors:
        raise ValueError('\n'.join(errors))

    # IUPAC names are not unique, e.g Bromine (atom) or (Br2)

    # Check formula is formatted right
    for i in db:
        assert i.formula == serialize_formula(i.formula)

    # Check CAS validity
    for i in db:
        assert check_CAS(i.CASs)

    # MW checker
    for i in db:
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_close(mw_calc, i.MW, atol=0.05)

def test_ion_dbs():
    db = ChemicalMetadataDB(elements=False,
                            main_db=None,
                            user_dbs=[os.path.join(folder, 'Anion db.tsv'), os.path.join(folder, 'Cation db.tsv')])

    # Load anion preferences
    anion_pref_file = os.path.join(folder, 'anion_preferences.json')
    if os.path.exists(anion_pref_file):
        with open(anion_pref_file) as f:
            anion_preferences = json.load(f)
            preferred_CAS = set(anion_preferences.get('preferred_cas', []))
            unpreferred_CAS = set(anion_preferences.get('unpreferred_cas', []))
    
    # Load cation preferences
    cation_pref_file = os.path.join(folder, 'cation_preferences.json')
    if os.path.exists(cation_pref_file):
        with open(cation_pref_file) as f:
            cation_preferences = json.load(f)
            # Update existing sets with cation preferences
            preferred_CAS.update(cation_preferences.get('preferred_cas', []))
            unpreferred_CAS.update(cation_preferences.get('unpreferred_cas', []))
                
    # Check CAS lookup
    for CAS, d in db.CAS_index.items():
        assert CAS_from_any(d.CASs) == d.CASs

    # Check formula lookups
    errors = []
    for CAS, d in db.CAS_index.items():
        formula = d.formula
        result = CAS_from_any(formula)
        if result != d.CASs:
            if int_to_CAS(CAS) not in unpreferred_CAS:
                errors.append(f"CAS {result} retrieved by formula {formula} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the smiles
    errors = []
    for CAS, d in db.CAS_index.items():
        smiles = d.smiles
        if smiles:
            result = CAS_from_any('smiles='+smiles)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by smiles {smiles} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the inchi
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI = d.InChI
        if InChI:
            result = CAS_from_any('InChI=1S/'+InChI)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI {InChI} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the InChI_key
    errors = []
    for CAS, d in db.CAS_index.items():
        InChI_key = d.InChI_key
        if InChI_key:
            result = CAS_from_any('InChIKey='+InChI_key)
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by InChI_key {InChI_key} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the pubchem
    errors = []
    for CAS, d in db.CAS_index.items():
        pubchemid = d.pubchemid
        if pubchemid != -1:
            result = CAS_from_any('PubChem='+str(pubchemid))
            if result != d.CASs:
                if int_to_CAS(CAS) not in unpreferred_CAS:
                    errors.append(f"CAS {result} retrieved by pubchemid {pubchemid} on compound {d.CASs} not marked as unpreferred")
    if errors:
        raise ValueError('\n'.join(errors))

    # Check the name
    errors = []
    for CAS, d in db.CAS_index.items():
        common_name = d.common_name
        result = CAS_from_any(common_name)
        if result != d.CASs:
            errors.append(f"CAS {result} retrieved by common_name {common_name} on compound {d.CASs}")
    if errors:
        raise ValueError('\n'.join(errors))

    # IUPAC names are not unique, e.g Bromine (atom) or (Br2)

    # Check formula is formatted right
    for i in db:
        assert i.formula == serialize_formula(i.formula)

    # Check CAS validity
    for i in db:
        assert check_CAS(i.CASs)

    # MW checker
    for i in db:
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_close(mw_calc, i.MW, atol=0.05)


def test_mixture_from_any():
    with pytest.raises(Exception):
        mixture_from_any(['water', 'methanol'])
    with pytest.raises(Exception):
        mixture_from_any('NOTAMIXTURE')

    for name in ['Air', 'air', 'Air', ['air']]:
        assert mixture_from_any(name) == common_mixtures['Air']

    names = ['R-401A ', ' R-401A ', 'R401A ', 'r401a', 'r-401A', 'refrigerant-401A', 'refrigerant 401A']
    for name in names:
        assert mixture_from_any(name) == common_mixtures['R401A']

    assert mixture_from_any('R512A') == common_mixtures['R512A']
    assert mixture_from_any(['air']) == common_mixtures['Air']

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
    assert hit0 == hit1

    with pytest.raises(ValueError):
        # Not that smart/weird
        search_chemical('(oxidane)')

    assert search_chemical('water').charge == 0

def test_search_chemical_2024_update_metadata_organic_names():
    assert search_chemical('112-95-8').common_name == "eicosane"
    assert search_chemical('103554-13-8').common_name == "naphthalene, 2-hexyldecahydro-"
    assert search_chemical('105-76-0').common_name == "dibutyl maleate"
    assert search_chemical('107-83-5').common_name == "2-Methylpentane"
    assert search_chemical('1120-25-8').common_name == "methyl palmitoleate"
    assert search_chemical('132739-31-2').common_name == "dipropylene glycol mono-tert-butyl ether"
    assert search_chemical('142-16-5').common_name == "bis(2-ethylhexyl) maleate"
    assert search_chemical('147-81-9').common_name == "arabinose"
    assert search_chemical('163702-07-6').common_name == "methyl nonafluorobutyl ether"
    assert search_chemical('16587-32-9').common_name == "benzo[b]thiophene, 2-propyl-"
    assert search_chemical('25498-49-1').common_name == "tripropylene glycol monomethyl ether"
    assert search_chemical('30025-38-8').common_name == "dipropylene glycol monoethyl ether"
    assert search_chemical('3458-28-4').common_name == "mannose"
    assert search_chemical('34590-94-8').common_name == "dipropylene glycol monomethyl ether"
    assert search_chemical('463-72-9').common_name == "carbamic chloride"
    assert search_chemical('498-23-7').common_name == "citraconic acid"
    assert search_chemical('50-99-7').common_name == "glucose"
    assert search_chemical('528-50-7').common_name == "cellobiose"
    assert search_chemical('57-48-7').common_name == "fructose"
    assert search_chemical('59-23-4').common_name == "galactose"
    assert search_chemical('75-28-5').common_name == "isobutane"
    assert search_chemical('75-94-5').common_name == "vinyltrichlorosilane"
    assert search_chemical('87-99-0').common_name == "xylitol"
    assert search_chemical('88917-22-0').common_name == "dipropylene glycol methyl ether acetate"
    assert search_chemical('9005-53-2').common_name == "lignin"

    assert search_chemical('C30').common_name == 'triacontane'
    assert search_chemical('n-C30').common_name == 'triacontane'
    assert search_chemical('nC30').common_name == 'triacontane'

def test_search_chemical_2024_update_metadata_organic_synonyms():
    # CAS: 74-98-6 (propane)
    assert "74-98-6" not in search_chemical('74-98-6').synonyms
    assert "PROPANE" not in search_chemical('74-98-6').synonyms
    assert "R 280" in search_chemical('74-98-6').synonyms

    # CAS: 74-82-8 (methane)
    assert "74-82-8" not in search_chemical('74-82-8').synonyms

    # CAS: 74-84-0 (ethane)
    assert "74-84-0" not in search_chemical('74-84-0').synonyms
    assert "ETHANE" not in search_chemical('74-84-0').synonyms
    assert "R 170 (HYDROCARBON)" not in search_chemical('74-84-0').synonyms
    assert "R 170 (hydrocarbon)" in search_chemical('74-84-0').synonyms

    # CAS: 7194-86-7 (nonatriacontane)
    assert "7194-86-7" not in search_chemical('7194-86-7').synonyms
    assert "SNXOSZZZNFRFNZ-UHFFFAOYSA-N" not in search_chemical('7194-86-7').synonyms

    # CAS: 7194-85-6 (octatriacontane)
    assert "7194-85-6" not in search_chemical('7194-85-6').synonyms
    assert "BVKCQBBZBGYNOP-UHFFFAOYSA-N" not in search_chemical('7194-85-6').synonyms
    assert "N-OCTATRIACONTANE" not in search_chemical('7194-85-6').synonyms
    assert "Octatriacontane" not in search_chemical('7194-85-6').synonyms
    assert "n-octatriacontane" in search_chemical('7194-85-6').synonyms

    # CAS: 638-67-5 (tricosane)
    assert "638-67-5" not in search_chemical('638-67-5').synonyms
    assert "TRICOSANE" not in search_chemical('638-67-5').synonyms

    # CAS: 638-68-6 (triacontane)
    assert "638-68-6" not in search_chemical('638-68-6').synonyms
    assert "TRIACONTANE" not in search_chemical('638-68-6').synonyms

    # CAS: 646-31-1 (tetracosane)
    assert "646-31-1" not in search_chemical('646-31-1').synonyms
    assert "TETRACOSANE" not in search_chemical('646-31-1').synonyms

    # CAS: 629-97-0 (docosane)
    assert "629-97-0" not in search_chemical('629-97-0').synonyms
    assert "DOCOSANE" not in search_chemical('629-97-0').synonyms
    assert "PARAFOL 22-95" not in search_chemical('629-97-0').synonyms
    assert "parafol 22-95" in search_chemical('629-97-0').synonyms

    # CAS: 57-48-7 (fructose)
    assert "D-Fructose" not in search_chemical('57-48-7').synonyms
    assert "d-fructose" in search_chemical('57-48-7').synonyms

    # CAS: 110-54-3 (hexane)
    assert "110-54-3" not in search_chemical('110-54-3').synonyms
    assert "HEXANE" not in search_chemical('110-54-3').synonyms

def test_search_chemical_2024_update_metadata_inorganic_names():
    # Testing for name updates on common inorganic compounds
    assert search_chemical('10025-74-8').common_name == "dysprosium chloride"
    assert search_chemical('10034-96-5').common_name == "manganese sulfate monohydrate"
    assert search_chemical('10043-01-3').common_name == "aluminum sulfate"
    assert search_chemical('10043-52-4').common_name == "calcium chloride"
    assert search_chemical('10102-44-0').common_name == "nitrogen dioxide"
    assert search_chemical('10361-37-2').common_name == "barium chloride"
    assert search_chemical('12018-01-8').common_name == "chromium oxide"
    assert search_chemical('12136-45-7').common_name == "potassium oxide"
    assert search_chemical('1313-82-2').common_name == "disodium sulfide"
    assert search_chemical('13709-49-4').common_name == "yttrium fluoride"
    assert search_chemical('21645-51-2').common_name == "aluminum hydroxide"
    assert search_chemical('10043-52-4').common_name == "calcium chloride"
    assert search_chemical('7791-12-0').common_name == "thallium monochloride"
    assert search_chemical('7785-87-7').common_name == "manganese sulfate"
    assert search_chemical('10099-59-9').common_name == "lanthanum nitrate"
    assert search_chemical('14457-87-5').common_name == "cerium bromide"


def test_search_chemical_2024_update_metadata_inorganic_synonyms():
    # Testing synonym updates for select inorganic chemicals

    # CAS: 10025-73-7 (chromium chloride (CrCl3))
    assert "CHROMIC CHLORIDE (CRCL3)" not in search_chemical('10025-73-7').synonyms
    assert "CHROMIUM TRICHLORIDE (CRCL3)" not in search_chemical('10025-73-7').synonyms
    assert "chromium chloride (CrCl3)" in search_chemical('10025-73-7').synonyms

    # CAS: 10031-25-1 (chromium bromide (CrBr3))
    assert "UZDWIWGMKWZEPE-UHFFFAOYSA-K" not in search_chemical('10031-25-1').synonyms
    assert "chromium bromide (CrBr3)" in search_chemical('10031-25-1').synonyms

    # CAS: 10043-01-3 (aluminum sulfate)
    assert "ALUMINUM SULFATE" not in search_chemical('10043-01-3').synonyms
    assert "aluminum sulfate" in search_chemical('10043-01-3').synonyms

    # CAS: 10102-44-0 (nitrogen dioxide)
    assert "nitrogen oxide (NO2)" in search_chemical('10102-44-0').synonyms

    # CAS: 10257-55-3 (calcium sulfite)
    assert "CALCIUM SULFITE (1:1)" not in search_chemical('10257-55-3').synonyms
    assert "sulfurous acid, calcium salt (1:1)" in search_chemical('10257-55-3').synonyms

    # CAS: 12045-63-5 (titanium diboride)
    assert "Titanium boride (TiB2)" not in search_chemical('12045-63-5').synonyms
    assert "titanium diboride" in search_chemical('12045-63-5').synonyms

    # CAS: 12136-45-7 (potassium oxide)
    assert "potassium oxide (K2O)" in search_chemical('12136-45-7').synonyms
    assert "CHWRSCGUEQEHOH-UHFFFAOYSA-N" not in search_chemical('12136-45-7').synonyms

    # CAS: 13765-26-9 (gadolinium fluoride)
    assert "Gadolinium fluoride (GdF3)" not in search_chemical('13765-26-9').synonyms
    assert "gadolinium fluoride (GdF3)" in search_chemical('13765-26-9').synonyms

    # CAS: 20427-58-1 (zinc hydroxide)
    assert "ZINC HYDROXIDE" not in search_chemical('20427-58-1').synonyms
    assert "zinc hydroxide (Zn(OH)2)" in search_chemical('20427-58-1').synonyms


    # CAS: 7757-79-1 (potassium nitrate)
    assert "nitric acid potassium salt (1:1)" in search_chemical('7757-79-1').synonyms
    assert "POTASSIUM NITRATE" not in search_chemical('7757-79-1').synonyms

    # CAS: 7646-85-7 (zinc chloride)
    assert "zinc chloride (ZnCl2)" in search_chemical('7646-85-7').synonyms

    # CAS: 7789-40-4 (thallium bromide (TlBr))
    assert "thallium bromide (TlBr)" in search_chemical('7789-40-4').synonyms

    # CAS: 7784-30-7 (aluminum phosphate)
    assert "ALUMINIUM PHOSPHATE" not in search_chemical('7784-30-7').synonyms


def test_search_chemical_2024_update_metadata_ion_names():
    assert search_chemical('51-92-3').common_name == "tetramethylammonium"
    assert search_chemical('11062-77-4').common_name == "superoxide"
    assert search_chemical('12184-88-2').common_name == "hydride"
    assert search_chemical('13408-62-3').common_name == "ferricyanide"
    assert search_chemical('13408-63-4').common_name == "ferrocyanide"
    assert search_chemical('13907-47-6').common_name == "dichromate"
    assert search_chemical('13981-20-9').common_name == "metavanadate"
    assert search_chemical('14000-31-8').common_name == "pyrophosphate"
    assert search_chemical('14124-67-5').common_name == "selenite"
    assert search_chemical('14124-68-6').common_name == "selenate"
    assert search_chemical('14127-68-5').common_name == "triphosphate"
    assert search_chemical('14265-44-2').common_name == "phosphate"
    assert search_chemical('14265-45-3').common_name == "sulfite"
    assert search_chemical('14280-30-9').common_name == "hydroxide"
    assert search_chemical('14333-13-2').common_name == "permanganate"
    assert search_chemical('14333-24-5').common_name == "perrhenate"
    assert search_chemical('14337-12-3').common_name == "tetrachloroaurate"
    assert search_chemical('14343-69-2').common_name == "azide"
    assert search_chemical('14380-62-2').common_name == "hypobromite"
    assert search_chemical('14797-55-8').common_name == "nitrate"
    assert search_chemical('14797-65-0').common_name == "nitrite"
    assert search_chemical('14808-79-8').common_name == "sulfate"
    assert search_chemical('14866-68-3').common_name == "chlorate"
    assert search_chemical('14900-04-0').common_name == "triiodide"
    assert search_chemical('14996-02-2').common_name == "hydrogen sulfate"
    assert search_chemical('14998-27-7').common_name == "chlorite"
    assert search_chemical('15181-46-1').common_name == "bisulfite"
    assert search_chemical('15454-31-6').common_name == "iodate"
    assert search_chemical('15541-45-4').common_name == "bromate"
    assert search_chemical('15584-04-0').common_name == "arsenate"
    assert search_chemical('16833-27-5').common_name == "oxide"
    assert search_chemical('16887-00-6').common_name == "chloride"
    assert search_chemical('16971-29-2').common_name == "borohydride"
    assert search_chemical('16984-48-8').common_name == "fluoride"
    assert search_chemical('17084-08-1').common_name == "fluorosilicate"
    assert search_chemical('17306-35-3').common_name == "arsenenite"
    assert search_chemical('18851-77-9').common_name == "nitride"
    assert search_chemical('20461-54-5').common_name == "iodide"
    assert search_chemical('20561-39-1').common_name == "astatide"
    assert search_chemical('23134-05-6').common_name == "disulfite"
    assert search_chemical('24959-67-9').common_name == "bromide"
    assert search_chemical('3812-32-6').common_name == "carbonate"

    # maybe these will be changed in the future but it's what CAS says for now

    assert search_chemical('14701-22-5').common_name == "Ni2+"
    assert search_chemical('15158-11-9').common_name == "Cu2+"
    assert search_chemical('16065-83-1').common_name == "Cr3+"
    assert search_chemical('16065-91-1').common_name == "Au3+"
    assert search_chemical('16096-89-2').common_name == "La3+"
    assert search_chemical('17341-24-1').common_name == "Li1+"
    assert search_chemical('18923-26-7').common_name == "Ce3+"
    assert search_chemical('18923-27-8').common_name == "Yb3+"
    assert search_chemical('22537-40-2').common_name == "Y3+"
    assert search_chemical('22537-48-0').common_name == "Cd2+"
    assert search_chemical('22541-12-4').common_name == "Ba2+"
    assert search_chemical('22541-18-0').common_name == "Eu3+"
    assert search_chemical('22541-19-1').common_name == "Gd3+"
    assert search_chemical('22541-53-3').common_name == "Co2+"
    assert search_chemical('22541-63-5').common_name == "Co3+"


def test_search_chemical_2024_update_metadata_ion_synonyms():
    # CAS: 14701-21-4 (silver ion)
    assert "FOIXSVOLVBLSDH-UHFFFAOYSA-N" not in search_chemical('14701-21-4').synonyms
    assert "SILVER ION" not in search_chemical('14701-21-4').synonyms
    assert "Silver(1+) ion" not in search_chemical('14701-21-4').synonyms
    assert "silver, ion (Ag1+)" in search_chemical('14701-21-4').synonyms

    # CAS: 15158-11-9 (copper(II) ion)
    assert "Copper(2+)" not in search_chemical('15158-11-9').synonyms
    assert "Copper(II) cation" not in search_chemical('15158-11-9').synonyms
    assert "Cupric ion" not in search_chemical('15158-11-9').synonyms
    assert "copper, ion (Cu2+)" in search_chemical('15158-11-9').synonyms

    # CAS: 17493-86-6 (copper(I) ion)
    assert "VMQMZMRVKUZKQL-UHFFFAOYSA-N" not in search_chemical('17493-86-6').synonyms
    assert "Cuprous ion" not in search_chemical('17493-86-6').synonyms
    assert "copper, ion (Cu1+)" in search_chemical('17493-86-6').synonyms

    # CAS: 11062-77-4 (superoxide)
    assert "HYPEROXIDE" not in search_chemical('11062-77-4').synonyms
    assert "Superoxide radical anion" not in search_chemical('11062-77-4').synonyms
    assert "superoxide (8CI,9CI)" in search_chemical('11062-77-4').synonyms

    # CAS: 16833-27-5 (oxide)
    assert "Oxygen(2-)" not in search_chemical('16833-27-5').synonyms

    # CAS: 16984-48-8 (fluoride)
    assert "FLUORINE ION(F1-)" not in search_chemical('16984-48-8').synonyms
    assert "Fluoride ion" not in search_chemical('16984-48-8').synonyms

    # CAS: 18496-25-8 (sulfide)
    assert "SULFUR, ION (S2-)" not in search_chemical('18496-25-8').synonyms
    assert "Sulfide(2-)" not in search_chemical('18496-25-8').synonyms
    assert "Sulphide" not in search_chemical('18496-25-8').synonyms
    assert "UCKMPCXJQFINFW-UHFFFAOYSA-N" not in search_chemical('18496-25-8').synonyms

    # CAS: 14913-52-1 (neodymium ion)
    assert "Neodymium(3+)" not in search_chemical('14913-52-1').synonyms
    assert "neodymium, ion (Nd3+)" in search_chemical('14913-52-1').synonyms

    # CAS: 15543-40-5 (zirconium ion)
    assert "ZR4+" not in search_chemical('15543-40-5').synonyms
    assert "Zirconium(4+)" not in search_chemical('15543-40-5').synonyms
    assert "zirconium, ion (Zr4+)" in search_chemical('15543-40-5').synonyms

    # CAS: 25215-10-5 (guanidinium)
    assert "GUANIDINIUM" not in search_chemical('25215-10-5').synonyms
    assert "ZRALSGWEFCBTJO-UHFFFAOYSA-O" not in search_chemical('25215-10-5').synonyms

def test_synonym_formula_capitalization_inorganic():
    # Possible synonyms will be added/removed - just remove the test if that happens.
    # Obviously check to make sure it was actually added/removed from the database though

    # Test that formulas in synonyms are properly capitalized
    
    # Aluminum phosphate and related compounds
    assert "aluminum phosphate (al(po4))" not in search_chemical('7784-30-7').synonyms
    assert "aluminum phosphate (Al(PO4))" in search_chemical('7784-30-7').synonyms
    
    # Sodium compounds
    assert "sodium borate (na2(bo2)2)" not in search_chemical('7775-19-1').synonyms
    assert "sodium borate (Na2(BO2)2)" in search_chemical('7775-19-1').synonyms
    assert "sodium diphosphate (na4(p2o7))" not in search_chemical('7722-88-5').synonyms
    assert "sodium diphosphate (Na4(P2O7))" in search_chemical('7722-88-5').synonyms
    
    # Deuterium compounds
    assert "chloro((2)h)" not in search_chemical('7698-05-7').synonyms
    assert "chloro((2)H)" in search_chemical('7698-05-7').synonyms
    
    # Phosphorus compounds
    assert "[p(o)oh]" not in search_chemical('6303-21-5').synonyms
    assert "[P(O)OH]" in search_chemical('6303-21-5').synonyms
    
    # Lead compounds
    assert "lead orthophosphate (pb3(po4)2)" not in search_chemical('7446-27-7').synonyms
    assert "lead orthophosphate (Pb3(PO4)2)" in search_chemical('7446-27-7').synonyms
    
    # Metal cyanides
    assert "cobalt cyanide (co(cn)2)" not in search_chemical('542-84-7').synonyms
    assert "copper cyanide (cu(cn))" not in search_chemical('544-92-3').synonyms
    assert "copper cyanide (Cu(CN))" in search_chemical('544-92-3').synonyms
    
    # Mercury compounds
    assert "mercury thiocyanate (hg(scn)2)" not in search_chemical('592-85-8').synonyms
    assert "mercury thiocyanate (Hg(SCN)2)" in search_chemical('592-85-8').synonyms
    
    # Silver compounds
    assert "silver cyanide (ag(cn))" not in search_chemical('506-64-9').synonyms
    assert "silver cyanide (Ag2(CN)2)" in search_chemical('506-64-9').synonyms
    
    # Iron compounds
    assert "ferrous hydroxide (fe(oh)2)" not in search_chemical('18624-44-7').synonyms
    assert "ferrous hydroxide (Fe(OH)2)" in search_chemical('18624-44-7').synonyms
    assert "iron phosphate (fe3(po4)2)" not in search_chemical('14940-41-1').synonyms
    assert "iron phosphate (Fe3(PO4)2)" in search_chemical('14940-41-1').synonyms
    
    # Lithium compounds
    assert "lithium azide (li(n3))" not in search_chemical('19597-69-4').synonyms
    assert "lithium azide (Li(N3))" in search_chemical('19597-69-4').synonyms
    
    # Cobalt compounds
    assert "cobalt nitrate (co(no3)3)" not in search_chemical('15520-84-0').synonyms
    assert "cobalt nitrate (Co(NO3)3)" in search_chemical('15520-84-0').synonyms
    
    # Silicon compounds
    assert "sodium fluosilicate (na2(sif6))" not in search_chemical('16893-85-9').synonyms
    assert "sodium fluosilicate (Na2(SiF6))" in search_chemical('16893-85-9').synonyms
    assert "talc (mg3h2(sio3)4)" not in search_chemical('14807-96-6').synonyms
    
    # Nickel compounds
    assert "nickel selenate (ni(seo4))" not in search_chemical('15060-62-5').synonyms
    assert "nickel selenate (Ni(SeO4))" in search_chemical('15060-62-5').synonyms
    
    # Strontium compounds
    assert "strontium thiosulfate (sr(s2o3))" not in search_chemical('15123-90-7').synonyms
    assert "strontium thiosulfate (Sr(S2O3))" in search_chemical('15123-90-7').synonyms
    
    # Manganese compounds
    assert "manganese molybdate(VI) (mnna(po4))" not in search_chemical('14013-15-1').synonyms
    assert "manganese molybdate(VI) (MnNa(PO4))" in search_chemical('14013-15-1').synonyms
    assert "manganese carbonyl (mn2(co)10)" not in search_chemical('10170-69-1').synonyms
    assert "manganese carbonyl (Mn2(CO)10)" in search_chemical('10170-69-1').synonyms
    
    # Zinc compounds
    assert "willemite (zn2(sio4))" not in search_chemical('14374-77-7').synonyms
    assert "willemite (Zn2(SiO4))" in search_chemical('14374-77-7').synonyms
    
    # Hafnium compounds
    assert "hafnium silicate (hf(sio4))" not in search_chemical('13870-13-8').synonyms
    assert "hafnium silicate (Hf(SiO4))" in search_chemical('13870-13-8').synonyms
    
    # Phosphorus compounds
    assert "(p(oh)3)" not in search_chemical('13598-36-2').synonyms
    assert "(P(OH)3)" in search_chemical('13598-36-2').synonyms
    
    # Tungsten compounds
    assert "silver tungstate(VI) (ag2(wo4))" not in search_chemical('13465-93-5').synonyms
    assert "silver tungstate(VI) (Ag2(WO4))" in search_chemical('13465-93-5').synonyms
    
    # Dysprosium compounds
    assert "dysprosium nitrate (dy(no3)3)" not in search_chemical('10143-38-1').synonyms
    assert "dysprosium nitrate (Dy(NO3)3)" in search_chemical('10143-38-1').synonyms
    
    # Cadmium compounds
    assert "cadmium nitrate (cd(no3)2) tetrahydrate" not in search_chemical('10022-68-1').synonyms
    assert "cadmium nitrate (Cd(NO3)2) tetrahydrate" in search_chemical('10022-68-1').synonyms

def test_synonym_formula_capitalization_anions():
    # Superoxide
    assert "Oxygen anion (O2-)" not in search_chemical('11062-77-4').synonyms
    assert "Oxygen anion radical (O2-)" not in search_chemical('11062-77-4').synonyms
    assert "Oxygen ion (O2-)" not in search_chemical('11062-77-4').synonyms
    assert "Superoxide ion (O2-)" not in search_chemical('11062-77-4').synonyms
    assert "oxygen anion (O2-)" in search_chemical('11062-77-4').synonyms
    assert "oxygen anion radical (O2-)" in search_chemical('11062-77-4').synonyms
    assert "oxygen ion (O2-)" in search_chemical('11062-77-4').synonyms
    assert "superoxide ion (O2-)" in search_chemical('11062-77-4').synonyms
    
    # Oxide
    assert "Oxide (O2-)" not in search_chemical('16833-27-5').synonyms
    assert "Oxygen, ion (O2-)" not in search_chemical('16833-27-5').synonyms
    assert "Oxygen, ion(O2-)" not in search_chemical('16833-27-5').synonyms
    assert "oxide (O2-)" in search_chemical('16833-27-5').synonyms
    assert "oxygen, ion (O2-)" in search_chemical('16833-27-5').synonyms
    assert "oxygen, ion(O2-)" in search_chemical('16833-27-5').synonyms
    
    # Azide
    assert "Azide (N3-)" not in search_chemical('14343-69-2').synonyms
    assert "Nitrogen ion (N3-)" not in search_chemical('14343-69-2').synonyms
    assert "Trinitrogen ion (N3-)" not in search_chemical('14343-69-2').synonyms
    assert "azide (N3-)" in search_chemical('14343-69-2').synonyms
    assert "nitrogen ion (N3-)" in search_chemical('14343-69-2').synonyms
    assert "trinitrogen ion (N3-)" in search_chemical('14343-69-2').synonyms
    
    # Tribromide
    assert "Bromide (Br3-)" not in search_chemical('14522-80-6').synonyms
    assert "bromide (Br3-)" in search_chemical('14522-80-6').synonyms
    
    # Nitrate
    assert "Nitrate (No3-)" not in search_chemical('14797-55-8').synonyms
    assert "nitrate (No3-)" in search_chemical('14797-55-8').synonyms
    
    # Nitride
    assert "Ammonia, ion (N3-)" not in search_chemical('18851-77-9').synonyms
    assert "Nitrogen, ion (N3-)" not in search_chemical('18851-77-9').synonyms
    assert "ammonia, ion (N3-)" in search_chemical('18851-77-9').synonyms
    assert "nitrogen, ion (N3-)" in search_chemical('18851-77-9').synonyms
    
    # Fluoride
    assert "Fluoride ion (F-)" not in search_chemical('16984-48-8').synonyms
    assert "Fluoride ion(F-)" not in search_chemical('16984-48-8').synonyms
    assert "Fluorine ion(F1-)" not in search_chemical('16984-48-8').synonyms
    assert "fluoride ion (F-)" in search_chemical('16984-48-8').synonyms
    assert "fluoride ion(F-)" in search_chemical('16984-48-8').synonyms
    assert "fluorine ion(F1-)" in search_chemical('16984-48-8').synonyms

def test_synonym_formula_capitalization_cations():
    # Mercury(2+)
    assert "Mercury (Hg2+)" not in search_chemical('14302-87-5').synonyms
    assert "Mercury cation (Hg2+)" not in search_chemical('14302-87-5').synonyms
    assert "Mercury ion (Hg2+)" not in search_chemical('14302-87-5').synonyms
    assert "mercury (Hg2+)" in search_chemical('14302-87-5').synonyms
    assert "mercury cation (Hg2+)" in search_chemical('14302-87-5').synonyms
    assert "mercury ion (Hg2+)" in search_chemical('14302-87-5').synonyms

    # Lead(2+)
    assert "Lead (Pb2+)" not in search_chemical('14280-50-3').synonyms
    assert "Lead ion (Pb2+)" not in search_chemical('14280-50-3').synonyms
    assert "lead (Pb2+)" in search_chemical('14280-50-3').synonyms
    assert "lead ion (Pb2+)" in search_chemical('14280-50-3').synonyms

    # Silver(1+)
    assert "Silver ion (Ag+)" not in search_chemical('14701-21-4').synonyms
    assert "silver ion (Ag+)" in search_chemical('14701-21-4').synonyms

    # Copper(2+)
    assert "Copper ion (Cu++)" not in search_chemical('15158-11-9').synonyms
    assert "Cupric ion (Cu2+)" not in search_chemical('15158-11-9').synonyms
    assert "copper ion (Cu++)" in search_chemical('15158-11-9').synonyms
    assert "cupric ion (Cu2+)" in search_chemical('15158-11-9').synonyms

    # Manganese(2+)
    assert "Manganese (Mn2+)" not in search_chemical('16397-91-4').synonyms
    assert "Manganese Ion (Mn2+)" not in search_chemical('16397-91-4').synonyms
    assert "Manganese cation (Mn2+)" not in search_chemical('16397-91-4').synonyms
    assert "manganese (Mn2+)" in search_chemical('16397-91-4').synonyms
    assert "manganese cation (Mn2+)" in search_chemical('16397-91-4').synonyms
    assert "manganese ion (Mn2+)" in search_chemical('16397-91-4').synonyms

    # Erbium(3+)
    assert "Erbium (Er3+)" not in search_chemical('18472-30-5').synonyms
    assert "erbium (Er3+)" in search_chemical('18472-30-5').synonyms

    # Cerium(3+)
    assert "Cerium (Ce3+)" not in search_chemical('18923-26-7').synonyms
    assert "cerium (Ce3+)" in search_chemical('18923-26-7').synonyms

    # Iron(3+)
    assert "Iron (Fe3+)" not in search_chemical('20074-52-6').synonyms
    assert "Iron, ion (Fe3+) (8CI,9CI)" not in search_chemical('20074-52-6').synonyms
    assert "iron (Fe3+)" in search_chemical('20074-52-6').synonyms
    assert "iron, ion (Fe3+) (8CI,9CI)" in search_chemical('20074-52-6').synonyms

    # Gold(1+)
    assert "Gold (Au1+)" not in search_chemical('20681-14-5').synonyms
    assert "gold (Au1+)" in search_chemical('20681-14-5').synonyms

    # Potassium(1+)
    assert "Potassium (K+)" not in search_chemical('24203-36-9').synonyms
    assert "Potassium ion (K+)" not in search_chemical('24203-36-9').synonyms
    assert "Potassium ion (K1+)" not in search_chemical('24203-36-9').synonyms
    assert "Potassium, ion (K1+) (8CI,9CI)" not in search_chemical('24203-36-9').synonyms
    assert "potassium (K+)" in search_chemical('24203-36-9').synonyms
    assert "potassium ion (K+)" in search_chemical('24203-36-9').synonyms
    assert "potassium ion (K1+)" in search_chemical('24203-36-9').synonyms
    assert "potassium, ion (K1+) (8CI,9CI)" in search_chemical('24203-36-9').synonyms

def test_synonym_formula_capitalization_organic():
    # N-[(Phenylmethoxy)carbonyl]glycyl-L-leucine
    assert "N-(N-((Phenylmethoxy)carbonyl)glycyl)-l-leucine" not in search_chemical('1421-69-8').synonyms
    assert "N-(N-((phenylmethoxy)carbonyl)glycyl)-l-leucine" in search_chemical('1421-69-8').synonyms

    # Diisopropyl phosphite
    assert "isopropyl phosphite ((c3h7o)2(ho)p)" not in search_chemical('1809-20-7').synonyms
    assert "isopropyl phosphite ((C3H7O)2(HO)P)" in search_chemical('1809-20-7').synonyms

    # Fulminic acid
    assert "[c(h)no]" not in search_chemical('506-85-4').synonyms
    assert "[ch(no)]" not in search_chemical('506-85-4').synonyms
    assert "[C(H)NO]" in search_chemical('506-85-4').synonyms
    assert "[CH(NO)]" in search_chemical('506-85-4').synonyms

    # N-(tert-Butyloxycarbonyl)-D-methionine
    assert "N-boc-(d)-methionine" not in search_chemical('5241-66-7').synonyms
    assert "N-boc-(D)-methionine" in search_chemical('5241-66-7').synonyms

    # D-Penicillamine disulfide
    assert "3,3′-dithiobis[d-valine]" not in search_chemical('20902-45-8').synonyms
    assert "3,3′-dithiobis[D-valine]" in search_chemical('20902-45-8').synonyms

    # Stannous stearate
    assert "tin stearate (sn(c18h35o2)2)" not in search_chemical('6994-59-8').synonyms
    assert "tin stearate (Sn(C18H35O2)2)" in search_chemical('6994-59-8').synonyms

    # Loratadine
    assert "11-[N-(ethoxycarbonyl)-4-piperidylidene]-8-chloro-6,11-dihydro-5H-benzo-[5,6]cyclohepta[1,2-b]pyridine" not in search_chemical('79794-75-5').synonyms
    assert "11-[N-(ethoxycarbonyl)-4-piperidylidene]-8-chloro-6,11-dihydro-5h-benzo-[5,6]cyclohepta[1,2-b]pyridine" in search_chemical('79794-75-5').synonyms

def test_avoiding_homopolymers_copolymers_polymer_fragments():
    # these can take the place of a regular chemical, very bad
    assert search_chemical('oxetane').CASs == '503-30-0'
    assert search_chemical('hex-1-ene').CASs == '592-41-6'
    assert search_chemical('benzene-1,4-dicarbonitrile').CASs == '623-26-7'
    assert search_chemical('1,3-dioxolane').CASs == '646-06-0'
    assert search_chemical('4-methylpent-1-ene').CASs == '691-37-2'
    assert search_chemical('1,6-diisocyanatohexane').CASs == '822-06-0'
    assert search_chemical('ethenyl-dimethyl-phenylsilane').CASs == '1125-26-4'
    assert search_chemical('1,4-dioxane-2,3-dione').CASs == '3524-70-7'
    assert search_chemical('buta-1,3-diene').CASs == '106-99-0'
    assert search_chemical('prop-1-ene').CASs == '115-07-1'
    assert search_chemical('oxacyclohexadecan-2-one').CASs == '106-02-5'
    assert search_chemical('propanal').CASs == '123-38-6'
    assert search_chemical('cyclooctene').CASs == '931-88-4'
    assert search_chemical('1-chloro-1,2,2-trifluoroethene').CASs == '79-38-9'
    assert search_chemical('oxolan-2-one').CASs == '96-48-0'
    assert search_chemical('1-isocyanato-4-[(4-isocyanatophenyl)methyl]benzene').CASs == '101-68-8'
    assert search_chemical('ethenyl(trimethyl)silane').CASs == '754-05-2'
    assert search_chemical('1,5-cyclooctadiene').CASs == '111-78-4'
    assert search_chemical('1,2,3,4,5-pentadeuterio-6-(1,2,2-trideuterioethenyl)benzene').CASs == '19361-62-7'
    assert search_chemical('Tributyltin fluoride').CASs == '1983-10-4'
    assert search_chemical('503-30-0').CASs == '503-30-0'
    assert search_chemical('592-41-6').CASs == '592-41-6'
    assert search_chemical('623-26-7').CASs == '623-26-7'
    assert search_chemical('646-06-0').CASs == '646-06-0'
    assert search_chemical('691-37-2').CASs == '691-37-2'
    assert search_chemical('822-06-0').CASs == '822-06-0'
    assert search_chemical('1125-26-4').CASs == '1125-26-4'
    assert search_chemical('3524-70-7').CASs == '3524-70-7'
    assert search_chemical('106-99-0').CASs == '106-99-0'
    assert search_chemical('115-07-1').CASs == '115-07-1'
    assert search_chemical('106-02-5').CASs == '106-02-5'
    assert search_chemical('123-38-6').CASs == '123-38-6'
    assert search_chemical('931-88-4').CASs == '931-88-4'
    assert search_chemical('79-38-9').CASs == '79-38-9'
    assert search_chemical('96-48-0').CASs == '96-48-0'
    assert search_chemical('101-68-8').CASs == '101-68-8'
    assert search_chemical('754-05-2').CASs == '754-05-2'
    assert search_chemical('111-78-4').CASs == '111-78-4'
    assert search_chemical('19361-62-7').CASs == '19361-62-7'
    assert search_chemical('1983-10-4').CASs == '1983-10-4'


def test_issue_28_thermo():
    issue_28 = search_chemical('16949-15-8')
    assert issue_28.formula == 'BH4Li'
    assert issue_28.common_name == 'lithium borohydride'
    assert_close(molecular_weight(nested_formula_parser(issue_28.formula)), issue_28.MW)
    assert issue_28.pubchemid == 4148881
    with pytest.raises(Exception):
        search_chemical('pubchem=20722760')


def test_issue_45_chemicals():
    chemical = search_chemical('trichlorosilane')
    assert chemical.pubchemid == 24811
    assert chemical.CAS == 10025782
    assert chemical.formula == 'Cl3HSi'
    assert_close(chemical.MW, 135.45244)
    assert chemical.smiles == 'Cl[SiH](Cl)Cl'
    assert chemical.InChI == 'Cl3HSi/c1-4(2)3/h4H'
    assert chemical.InChI_key == 'ZDHXKXAHOVTTAH-UHFFFAOYSA-N'
    assert chemical.iupac_name == 'trichlorosilane'
    assert chemical.common_name == 'trichlorosilane'
    chemical = search_chemical('dichlorosilane')
    assert chemical.pubchemid == 61330
    assert chemical.CAS == 4109960
    assert chemical.formula == 'Cl2H2Si'
    assert_close(chemical.MW, 101.00738)
    assert chemical.smiles == 'Cl[SiH2]Cl'
    assert chemical.InChI == 'Cl2H2Si/c1-3-2/h3H2'
    assert chemical.InChI_key == 'MROCJMGDEKINLD-UHFFFAOYSA-N'
    assert chemical.iupac_name == 'dichlorosilane'
    assert chemical.common_name == 'dichlorosilane'

    chemical = search_chemical('chlorosilane')
    assert chemical.pubchemid == 61622
    assert chemical.CAS == 13465786 
    assert chemical.formula == 'ClH3Si'
    assert_close(chemical.MW, 66.56232)
    assert chemical.smiles == '[SiH3]Cl'
    assert chemical.iupac_name == 'chlorosilane'
    assert chemical.common_name == 'chlorosilane'


    # Vanadium oxytrichloride
    chemical = search_chemical("vanadium oxytrichloride")
    assert chemical.CAS == 7727186
    assert chemical.formula == "Cl3OV"
    assert chemical.pubchemid is not None
    
    # Methyldichlorosilane 
    chemical = search_chemical("methyldichlorosilane")
    assert chemical.CAS == 75547
    assert chemical.formula == "CH4Cl2Si"
    assert chemical.pubchemid is not None
    
    # Trimethoxysilane
    chemical = search_chemical("trimethoxysilane")
    assert chemical.CAS == 2487903
    assert chemical.formula == "C3H10O3Si"
    assert chemical.pubchemid is not None
    
    # 2-hydroxybut-3-enenitrile
    chemical = search_chemical("2-hydroxybut-3-enenitrile")
    assert chemical.CAS == 5809596
    assert chemical.formula == "C4H5NO"
    assert chemical.pubchemid is not None
    
    # Chlorodimethylsilane
    chemical = search_chemical("chlorodimethylsilane")
    assert chemical.CAS == 1066359
    assert chemical.formula == "C2H7ClSi"
    assert chemical.pubchemid is not None
    
    # Dimethylsilane
    chemical = search_chemical("dimethylsilane")
    assert chemical.CAS == 1111746
    assert chemical.formula == "C2H8Si"
    assert chemical.pubchemid is not None
    
    # Methylsilane
    chemical = search_chemical("methylsilane")
    assert chemical.CAS == 992949
    assert chemical.formula == "CH6Si"
    assert chemical.pubchemid is not None
    
    # Trimethylsilane
    chemical = search_chemical("trimethylsilane")
    assert chemical.CAS == 993077
    assert chemical.formula == "C3H10Si"
    assert chemical.pubchemid is not None
    
    # Silylsilane
    chemical = search_chemical("silylsilane")
    assert chemical.CAS == 1590870
    assert chemical.formula == "H6Si2"
    assert chemical.pubchemid is not None
    
    # Disilylsilane
    chemical = search_chemical("disilylsilane")
    assert chemical.CAS == 7783268
    assert chemical.formula == "H8Si3"
    assert chemical.pubchemid is not None
    
    # Monosodium glutamate
    chemical = search_chemical("monosodium glutamate")
    assert chemical.CAS == 142472
    assert chemical.formula == "C5H9NNaO4"
    assert chemical.pubchemid is not None

def test_diborane():
   # Test basic name lookup
   chemical = search_chemical("diborane")
   assert chemical.pubchemid == 12544637
   assert chemical.CAS == 19287457
   assert chemical.formula == "B2H6"
   # smiles confirmed, pubchem has the wrong structure https://en.wikipedia.org/wiki/Diborane
   assert chemical.smiles == "[BH2]1[H][BH2][H]1"
   assert chemical.CAS == 19287457
   chemical = search_chemical("Diboron hexahydride") 
   assert chemical.CAS == 19287457

def test_sodium_borohydride():
   chemical = search_chemical('sodium borohydride')
   assert chemical.pubchemid == 4311764
   assert chemical.CAS == 16940662
   assert chemical.formula == 'BH4Na' 
   assert_close(chemical.MW, 37.83253)
   assert chemical.smiles == '[BH4-].[Na+]'
   assert chemical.common_name == 'sodium borohydride'
   
   def test_acetylene():
    chemical = search_chemical('acetylene')
    assert chemical.pubchemid == 6326
    assert chemical.CAS == 74862
    assert chemical.formula == 'C2H2'
    assert_close(chemical.MW, 26.03728)
    assert chemical.smiles == 'C#C'
    assert chemical.InChI == 'C2H2/c1-2/h1-2H'
    assert chemical.InChI_key == 'HSFWRNGVRCDJHI-UHFFFAOYSA-N'
    assert chemical.common_name == 'acetylene'
    assert chemical.iupac_name == 'ethyne'

def test_water():
    chemical = search_chemical('water')
    assert chemical.pubchemid == 962
    assert chemical.CAS == 7732185
    assert chemical.formula == 'H2O'
    assert_close(chemical.MW, 18.01528)
    assert chemical.smiles == 'O'
    assert chemical.InChI == 'H2O/h1H2'
    assert chemical.InChI_key == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    assert chemical.iupac_name == 'oxidane'
    assert chemical.common_name == 'water'
    
    # Check that the problematic string with wrong separators isn't in synonyms, https://github.com/CalebBell/chemicals/issues/29
    bad_string = 'caustic soda liquid;aquafina;distilled water;hydrogen oxide (h2o);ultrexii ultrapure;'
    assert bad_string not in chemical.synonyms

def test_atomic_oxygen():
    chemical = search_chemical('atomic oxygen')
    # Test basic properties
    assert chemical.pubchemid == 159832
    assert chemical.CASs == '17778-80-2'
    assert chemical.formula == 'O'
    assert_close(chemical.MW, 15.9994)
    assert chemical.smiles == '[O]'
    assert chemical.InChI == 'O'
    assert chemical.InChI_key == 'QVGXLLKOCUKJST-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic oxygen'
    assert chemical.iupac_name == 'atomic oxygen'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'oxygen atom',
        'monooxygen',
        'ATOMIC OXYGEN',
        'OXYGEN ATOM'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '17778-80-2', f"Failed to find atomic oxygen using term: {term}"

    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'O',
        'smiles': '[O]',
        'inchi': 'InChI=1S/O',
        'inchikey': 'InChIKey=QVGXLLKOCUKJST-UHFFFAOYSA-N',
        'pubchem': 'pubchem=159832',
        'cas': '17778-80-2'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '17778-80-2', f"Failed to find atomic oxygen using {search_type}: {identifier}"
    
def test_molecular_oxygen():
    chemical = search_chemical('oxygen')
    # Test basic properties
    assert chemical.pubchemid == 977
    assert chemical.CASs == '7782-44-7'
    assert chemical.formula == 'O2'
    assert_close(chemical.MW, 31.9988)
    assert chemical.smiles == 'O=O'
    assert chemical.InChI == 'O2/c1-2'
    assert chemical.InChI_key == 'MYMOFIZGZYHOMD-UHFFFAOYSA-N'
    assert chemical.common_name == 'oxygen'
    assert chemical.iupac_name == 'molecular oxygen'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'dioxygen',
        'molecular oxygen',
        'oxygen gas',
        'O2',
        'oxygen-16',
        'liquid oxygen'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7782-44-7', f"Failed to find molecular oxygen using term: {term}"

    identifier_searches = {
        'formula': 'O2',
        'smiles': 'O=O',
        'inchi': 'InChI=1S/O2/c1-2',
        'inchikey': 'InChIKey=MYMOFIZGZYHOMD-UHFFFAOYSA-N',
        'pubchem': 'pubchem=977',
        'cas': '7782-44-7'
    }
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7782-44-7', f"Failed to find molecular oxygen using {search_type}: {identifier}"

def test_atomic_nitrogen():
    chemical = search_chemical('atomic nitrogen')
    # Test basic properties
    assert chemical.pubchemid == 57370662
    assert chemical.CASs == '17778-88-0'
    assert chemical.formula == 'N'
    assert_close(chemical.MW, 14.0067)
    assert chemical.smiles == '[N]'
    assert chemical.InChI == 'N'
    assert chemical.InChI_key == 'QJGQUHMNIGDVPM-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic nitrogen'
    assert chemical.iupac_name == 'atomic nitrogen'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'nitrogen atom',
        'mononitrogen',
        'ATOMIC NITROGEN',
        'NITROGEN ATOM',
        'nitrogen, atomic'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '17778-88-0', f"Failed to find atomic nitrogen using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'N',
        'smiles': '[N]',
        'inchi': 'InChI=1S/N',
        'inchikey': 'InChIKey=QJGQUHMNIGDVPM-UHFFFAOYSA-N',
        'pubchem': 'pubchem=57370662',
        'cas': '17778-88-0'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '17778-88-0', f"Failed to find atomic nitrogen using {search_type}: {identifier}"
    
def test_molecular_nitrogen():
    chemical = search_chemical('nitrogen')
    # Test basic properties
    assert chemical.pubchemid == 947
    assert chemical.CASs == '7727-37-9'
    assert chemical.formula == 'N2'
    assert_close(chemical.MW, 28.0134)
    assert chemical.smiles == 'N#N'
    assert chemical.InChI == 'N2/c1-2'
    assert chemical.InChI_key == 'IJGRMHOSHXDMSA-UHFFFAOYSA-N'
    assert chemical.common_name == 'nitrogen'
    assert chemical.iupac_name == 'molecular nitrogen'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'diatomic nitrogen',
        'molecular nitrogen',
        'nitrogen gas',
        'N2',
        'dinitrogen',
        'nitrogen molecule',
        'NITROGEN',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7727-37-9', f"Failed to find molecular nitrogen using term: {term}"
    
    identifier_searches = {
        'formula': 'N2',
        'smiles': 'N#N',
        'inchi': 'InChI=1S/N2/c1-2',
        'inchikey': 'InChIKey=IJGRMHOSHXDMSA-UHFFFAOYSA-N',
        'pubchem': 'pubchem=947',
        'cas': '7727-37-9'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7727-37-9', f"Failed to find molecular nitrogen using {search_type}: {identifier}"

def test_atomic_fluorine():
    chemical = search_chemical('atomic fluorine')
    # Test basic properties
    assert chemical.pubchemid == 5360525
    assert chemical.CASs == '14762-94-8'
    assert chemical.formula == 'F'
    assert_close(chemical.MW, 18.998403)
    assert chemical.smiles == '[F]'
    assert chemical.InChI == 'F'
    assert chemical.InChI_key == 'YCKRFDGAMUMZLT-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic fluorine'
    assert chemical.iupac_name == 'atomic fluorine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'fluorine atom',
        'monofluorine',
        'ATOMIC FLUORINE',
        'FLUORINE ATOM',
        'fluorine, atomic',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '14762-94-8', f"Failed to find atomic fluorine using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'F',
        'smiles': '[F]',
        'inchi': 'InChI=1S/F',
        'inchikey': 'InChIKey=YCKRFDGAMUMZLT-UHFFFAOYSA-N',
        'pubchem': 'pubchem=5360525',
        'cas': '14762-94-8'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '14762-94-8', f"Failed to find atomic fluorine using {search_type}: {identifier}"
    
def test_molecular_fluorine():
    chemical = search_chemical('fluorine')
    # Test basic properties
    assert chemical.pubchemid == 24524
    assert chemical.CASs == '7782-41-4'
    assert chemical.formula == 'F2'
    assert_close(chemical.MW, 37.996806)
    assert chemical.smiles == 'FF'
    assert chemical.InChI == 'F2/c1-2'
    assert chemical.InChI_key == 'PXGOKWXKJXAPGV-UHFFFAOYSA-N'
    assert chemical.common_name == 'fluorine'
    assert chemical.iupac_name == 'molecular fluorine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'diatomic fluorine',
        'molecular fluorine',
        'fluorine gas',
        'F2',
        'difluorine',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7782-41-4', f"Failed to find molecular fluorine using term: {term}"
    
    identifier_searches = {
        'formula': 'F2',
        'smiles': 'FF',
        'inchi': 'InChI=1S/F2/c1-2',
        'inchikey': 'InChIKey=PXGOKWXKJXAPGV-UHFFFAOYSA-N',
        'pubchem': 'pubchem=24524',
        'cas': '7782-41-4'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7782-41-4', f"Failed to find molecular fluorine using {search_type}: {identifier}"

def test_atomic_chlorine():
    chemical = search_chemical('atomic chlorine')
    # Test basic properties
    assert chemical.pubchemid == 5360523
    assert chemical.CASs == '22537-15-1'
    assert chemical.formula == 'Cl'
    assert_close(chemical.MW, 35.453)
    assert chemical.smiles == '[Cl]'
    assert chemical.InChI == 'Cl'
    assert chemical.InChI_key == 'ZAMOUSCENKQFHK-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic chlorine'
    assert chemical.iupac_name == 'atomic chlorine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'chlorine atom',
        'monochlorine',
        'ATOMIC CHLORINE',
        'CHLORINE ATOM',
        'chlorine, atomic',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '22537-15-1', f"Failed to find atomic chlorine using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'Cl',
        'smiles': '[Cl]',
        'inchi': 'InChI=1S/Cl',
        'inchikey': 'InChIKey=ZAMOUSCENKQFHK-UHFFFAOYSA-N',
        'pubchem': 'pubchem=5360523',
        'cas': '22537-15-1'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '22537-15-1', f"Failed to find atomic chlorine using {search_type}: {identifier}"
    
def test_molecular_chlorine():
    chemical = search_chemical('chlorine')
    # Test basic properties
    assert chemical.pubchemid == 24526
    assert chemical.CASs == '7782-50-5'
    assert chemical.formula == 'Cl2'
    assert_close(chemical.MW, 70.906)
    assert chemical.smiles == 'ClCl'
    assert chemical.InChI == 'Cl2/c1-2'
    assert chemical.InChI_key == 'KZBUYRJDOAKODT-UHFFFAOYSA-N'
    assert chemical.common_name == 'chlorine'
    assert chemical.iupac_name == 'molecular chlorine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'diatomic chlorine',
        'molecular chlorine',
        'Cl2',
        'dichlorine',
        'chlorine molecule',
        'chlorine mol.',
        'liquid chlorine',
        'chlorinum',
        'bertholite'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7782-50-5', f"Failed to find molecular chlorine using term: {term}"
        
    identifier_searches = {
        'formula': 'Cl2',
        'smiles': 'ClCl',
        'inchi': 'InChI=1S/Cl2/c1-2',
        'inchikey': 'InChIKey=KZBUYRJDOAKODT-UHFFFAOYSA-N',
        'pubchem': 'pubchem=24526',
        'cas': '7782-50-5'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7782-50-5', f"Failed to find molecular chlorine using {search_type}: {identifier}"

def test_atomic_iodine():
    chemical = search_chemical('atomic iodine')
    # Test basic properties
    assert chemical.pubchemid == 5360629
    assert chemical.CASs == '14362-44-8'
    assert chemical.formula == 'I'
    assert_close(chemical.MW, 126.90447)
    assert chemical.smiles == '[I]'
    assert chemical.InChI == 'I'
    assert chemical.InChI_key == 'ZCYVEMRRCGMTRW-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic iodine'
    assert chemical.iupac_name == 'atomic iodine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'iodine atom',
        'monoiodine',
        'ATOMIC IODINE',
        'IODINE ATOM',
        'iodine, atomic',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '14362-44-8', f"Failed to find atomic iodine using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'I',
        'smiles': '[I]',
        'inchi': 'InChI=1S/I',
        'inchikey': 'InChIKey=ZCYVEMRRCGMTRW-UHFFFAOYSA-N',
        'pubchem': 'pubchem=5360629',
        'cas': '14362-44-8'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '14362-44-8', f"Failed to find atomic iodine using {search_type}: {identifier}"
    
def test_molecular_iodine():
    chemical = search_chemical('iodine')
    # Test basic properties
    assert chemical.pubchemid == 807
    assert chemical.CASs == '7553-56-2'
    assert chemical.formula == 'I2'
    assert_close(chemical.MW, 253.80894)
    assert chemical.smiles == 'II'
    assert chemical.InChI == 'I2/c1-2'
    assert chemical.InChI_key == 'PNDPGZBMCMUPRI-UHFFFAOYSA-N'
    assert chemical.common_name == 'iodine'
    assert chemical.iupac_name == 'molecular iodine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'diatomic iodine',
        'molecular iodine',
        'I2',
        'diiodine',
        'iodine molecule',
        'iodium',
        'iodine (II)',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7553-56-2', f"Failed to find molecular iodine using term: {term}"
        identifier_searches = {
        'formula': 'I2',
        'smiles': 'II',
        'inchi': 'InChI=1S/I2/c1-2',
        'inchikey': 'InChIKey=PNDPGZBMCMUPRI-UHFFFAOYSA-N',
        'pubchem': 'pubchem=807',
        'cas': '7553-56-2'
    }
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7553-56-2', f"Failed to find molecular iodine using {search_type}: {identifier}"

def test_atomic_bromine():
    chemical = search_chemical('atomic bromine')
    # Test basic properties
    assert chemical.pubchemid == 5360770
    assert chemical.CASs == '10097-32-2'
    assert chemical.formula == 'Br'
    assert_close(chemical.MW, 79.904)
    assert chemical.smiles == '[Br]'
    assert chemical.InChI == 'Br'
    assert chemical.InChI_key == 'WKBOTKDWSSQWDR-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic bromine'
    assert chemical.iupac_name == 'atomic bromine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'bromine atom',
        'monobromine',
        'ATOMIC BROMINE',
        'BROMINE ATOM',
        'bromine, atomic',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '10097-32-2', f"Failed to find atomic bromine using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'Br',
        'smiles': '[Br]',
        'inchi': 'InChI=1S/Br',
        'inchikey': 'InChIKey=WKBOTKDWSSQWDR-UHFFFAOYSA-N',
        'pubchem': 'pubchem=5360770',
        'cas': '10097-32-2'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '10097-32-2', f"Failed to find atomic bromine using {search_type}: {identifier}"
    
def test_molecular_bromine():
    chemical = search_chemical('bromine')
    # Test basic properties
    assert chemical.pubchemid == 24408
    assert chemical.CASs == '7726-95-6'
    assert chemical.formula == 'Br2'
    assert_close(chemical.MW, 159.808)
    assert chemical.smiles == 'BrBr'
    assert chemical.InChI == 'Br2/c1-2'
    assert chemical.InChI_key == 'GDTBXPJZTBHREO-UHFFFAOYSA-N'
    assert chemical.common_name == 'bromine'
    assert chemical.iupac_name == 'molecular bromine'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'diatomic bromine',
        'molecular bromine',
        'bromine molecule',
        'Br2',
        'dibromine',
        'brom',
        'brome',
        'bromo',
        'bromium'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7726-95-6', f"Failed to find molecular bromine using term: {term}"
    
    identifier_searches = {
        'formula': 'Br2',
        'smiles': 'BrBr',
        'inchi': 'InChI=1S/Br2/c1-2',
        'inchikey': 'InChIKey=GDTBXPJZTBHREO-UHFFFAOYSA-N',
        'pubchem': 'pubchem=24408',
        'cas': '7726-95-6'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7726-95-6', f"Failed to find molecular bromine using {search_type}: {identifier}"

def test_parahydrogen():
    chemical = search_chemical('parahydrogen')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099490000-00-0'
    assert chemical.formula == 'H2'
    assert_close(chemical.MW, 2.01588)  # molecular weight of H2
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'parahydrogen'
    assert chemical.iupac_name == 'parahydrogen'
    # Test some key synonyms
    assert search_chemical('p-H2') == chemical
    assert search_chemical('pH2') == chemical
    assert search_chemical('para-hydrogen') == chemical

def test_orthohydrogen():
    chemical = search_chemical('orthohydrogen')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099479000-00-0'
    assert chemical.formula == 'H2'
    assert_close(chemical.MW, 2.01588)
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'orthohydrogen'
    assert chemical.iupac_name == 'orthohydrogen'
    # Test some key synonyms
    assert search_chemical('o-H2') == chemical
    assert search_chemical('ortho-H2') == chemical
    assert search_chemical('ortho-hydrogen') == chemical

def test_normal_hydrogen():
    chemical = search_chemical('normal hydrogen')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099474000-00-0'
    assert chemical.formula == 'H2'
    assert_close(chemical.MW, 2.01588)
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'normal hydrogen'
    assert chemical.iupac_name == 'normal hydrogen'
    # Test some key synonyms
    assert search_chemical('n-H2') == chemical
    assert search_chemical('nH2') == chemical
    assert '75:25 ortho:para hydrogen mixture' in chemical.synonyms


def test_equilibrium_hydrogen():
    chemical = search_chemical('hydrogen')
    assert chemical.pubchemid == 783
    assert chemical.CASs == '1333-74-0'
    assert chemical.formula == 'H2'
    assert_close(chemical.MW, 2.01588)
    assert chemical.smiles == '[HH]'
    assert chemical.InChI == 'H2/h1H'
    assert chemical.InChI_key == 'UFHFLCQGNIYNRP-UHFFFAOYSA-N'
    assert chemical.common_name == 'hydrogen'
    assert chemical.iupac_name == 'molecular hydrogen'
    # Test some key synonyms from the provided list
    assert search_chemical('dihydrogen') == chemical
    assert search_chemical('hydrogen molecule') == chemical
    assert search_chemical('equilibrium hydrogen') == chemical

def test_atomic_hydrogen():
    chemical = search_chemical('atomic hydrogen')
    # Test basic properties
    assert chemical.pubchemid == 5362549
    assert chemical.CASs == '12385-13-6'
    assert chemical.formula == 'H'
    assert_close(chemical.MW, 1.00794)
    assert chemical.smiles == '[H]'
    assert chemical.InChI == 'H'
    assert chemical.InChI_key == 'YZCKVEUIGOORGS-UHFFFAOYSA-N'
    assert chemical.common_name == 'atomic hydrogen'
    assert chemical.iupac_name == 'atomic hydrogen'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'hydrogen atom',
        'protium',
        'monatomic hydrogen'
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '12385-13-6', f"Failed to find atomic hydrogen using term: {term}"

def test_equilibrium_deuterium():
    chemical = search_chemical('deuterium')
    assert chemical.pubchemid == 24523
    assert chemical.CASs == '7782-39-0'
    assert chemical.formula == 'D2'
    assert_close(chemical.MW, 4.028204)
    assert chemical.smiles == '[2H][2H]'
    assert chemical.InChI == 'H2/h1H/i1+1D'
    assert chemical.InChI_key == 'UFHFLCQGNIYNRP-VVKOMZTBSA-N'
    assert chemical.common_name == 'deuterium'
    assert chemical.iupac_name == 'deuterium'
    # Test some key synonyms
    assert search_chemical('dideuterium') == chemical
    assert search_chemical('deuterium molecule') == chemical
    assert search_chemical('heavy hydrogen') == chemical
    assert search_chemical('hydrogen-2') == chemical

def test_paradeuterium():
    chemical = search_chemical('paradeuterium')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099458000-00-0'
    assert chemical.formula == 'D2'
    assert_close(chemical.MW, 4.028204)  # molecular weight of D2
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'paradeuterium'
    assert chemical.iupac_name == 'paradeuterium'
    # Test some key synonyms
    assert search_chemical('p-D2') == chemical
    assert search_chemical('pD2') == chemical
    assert search_chemical('para-deuterium') == chemical

def test_orthodeuterium():
    chemical = search_chemical('orthodeuterium')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099453000-00-0'
    assert chemical.formula == 'D2'
    assert_close(chemical.MW, 4.028204)
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'orthodeuterium'
    assert chemical.iupac_name == 'orthodeuterium'
    # Test some key synonyms
    assert search_chemical('o-D2') == chemical
    assert search_chemical('ortho-D2') == chemical
    assert search_chemical('ortho-deuterium') == chemical

def test_normal_deuterium():
    chemical = search_chemical('normal deuterium')
    assert chemical.pubchemid == -1
    assert chemical.CASs == '2099437000-00-0'
    assert chemical.formula == 'D2'
    assert_close(chemical.MW, 4.028204)
    assert chemical.smiles == ''
    assert chemical.InChI == ''
    assert chemical.InChI_key == ''
    assert chemical.common_name == 'normal deuterium'
    assert chemical.iupac_name == 'normal deuterium'
    # Test some key synonyms
    assert search_chemical('n-D2') == chemical
    assert search_chemical('nD2') == chemical
    assert '2:1 ortho:para deuterium mixture' in chemical.synonyms

def test_helium():
    chemical = search_chemical('helium')
    # Test basic properties
    assert chemical.pubchemid == 23987
    assert chemical.CASs == '7440-59-7'
    assert chemical.formula == 'He'
    assert_close(chemical.MW, 4.002602)
    assert chemical.smiles == '[He]'
    assert chemical.InChI == 'He'
    assert chemical.InChI_key == 'SWQJXJOGLNCZEY-UHFFFAOYSA-N'
    assert chemical.common_name == 'helium'
    assert chemical.iupac_name == 'helium'
    
    # Test the most likely search terms a user would try
    common_searches = [
        'helium',
        'atomic helium',
        'helium-4',
        'p-helium',
        'o-helium',
        'helio',
        '[He]',
        '494798-31-1',
    ]
    for term in common_searches:
        found_chemical = search_chemical(term)
        assert found_chemical.CASs == '7440-59-7', f"Failed to find helium using term: {term}"
    
    # Test identifier-based searches
    identifier_searches = {
        'formula': 'He',
        'smiles': '[He]',
        'inchi': 'InChI=1S/He',
        'inchikey': 'InChIKey=SWQJXJOGLNCZEY-UHFFFAOYSA-N',
        'pubchem': 'pubchem=23987',
        'cas': '7440-59-7'
    }
    
    for search_type, identifier in identifier_searches.items():
        found_chemical = search_chemical(identifier)
        assert found_chemical.CASs == '7440-59-7', f"Failed to find helium using {search_type}: {identifier}"
    

        
def test_absence_of_air_as_a_compound():
    # this compound has a CAS of air, but it's nonsense https://pubchem.ncbi.nlm.nih.gov/compound/195130
    with pytest.raises(Exception):
        chemical = search_chemical('132259-10-0')

        

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


    assert CAS_from_any('99685-  96-8') == '99685-96-8'
    assert CAS_from_any('99685-  96-    8') == '99685-96-8'
    assert CAS_from_any('    99685-96-8') == '99685-96-8'
    assert CAS_from_any('    99685  -  96    -    8   ') == '99685-96-8'


def test_periodic_table_variants():


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

        if d.PubChem is not None:
            assert CAS_from_any('PubChem=' + str(d.PubChem)) == CAS, (d, d.PubChem, CAS)
        else:
            fail += 1
    assert fail == 9
    # 111 - 118 aren't in pubchem

    # Direct checks on specific elements
    assert search_chemical('bromine').formula == 'Br2'
    assert search_chemical('iodine').formula == 'I2'

    assert search_chemical('Oxygen').formula == 'O2'
    assert search_chemical('nitrogen').formula == 'N2'
    assert search_chemical('fluorine').formula == 'F2'
    assert search_chemical('hydrogen').formula == 'H2'
    assert search_chemical('chlorine').formula == 'Cl2'


    assert search_chemical('monatomic bromine').formula == 'Br'
    assert search_chemical('monatomic iodine').formula == 'I'

    assert search_chemical('monatomic Oxygen').formula == 'O'
    assert search_chemical('monatomic nitrogen').formula == 'N'
    assert search_chemical('monatomic fluorine').formula == 'F'
    assert search_chemical('monatomic hydrogen').formula == 'H'
    assert search_chemical('monatomic chlorine').formula == 'Cl'

    # Check they can be looked up by their specific CAS number also
    assert search_chemical(search_chemical('monatomic bromine').CASs).formula == 'Br'
    assert search_chemical(search_chemical('monatomic iodine').CASs).formula == 'I'

    assert search_chemical(search_chemical('monatomic Oxygen').CASs).formula == 'O'
    assert search_chemical(search_chemical('monatomic nitrogen').CASs).formula == 'N'
    assert search_chemical(search_chemical('monatomic fluorine').CASs).formula == 'F'
    assert search_chemical(search_chemical('monatomic hydrogen').CASs).formula == 'H'
    assert search_chemical(search_chemical('monatomic chlorine').CASs).formula == 'Cl'

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
    assert {'Chlorine', 'Fluorine', 'Hydrogen', 'Nitrogen', 'Oxygen', 'Bromine', 'Iodine'} == set(failed_CASs)


def test_nested_brackets_with_chemical_with_brackets():
    assert search_chemical('1-(1-methylethyl)-4- methylbenzene (p-cymene)').CASs =='99-87-6'

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
    FAILED AT. THE TIME IS ONLY TAKEN by the PARSE function.

    EVEN THAT HAS BEEN REDUCED By 80% by using cElementTree instead of
    ElementTree.
    """
    import xml.etree.ElementTree as ET
    folder = os.path.join(os.path.dirname(__file__), '..', 'chemicals', 'Misc')

    tree = ET.parse(os.path.join(folder, 'ChemSep8.32.xml'))
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

        if CAS == '132259-10-0':
            # CAS for Air
            continue
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
# test_db_vs_ChemSep()

def test_heos_data_CASs():
    for CAS in heos_data.index.tolist():
        obj = search_chemical(CAS)
        assert obj.CASs == CAS
    for CAS, name in zip(heos_data.index.tolist(), heos_data['name'].tolist()):
        assert search_chemical(name).CASs == CAS


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



def test_ChemicalMetadata_basics():
    from copy import deepcopy
    import pickle
    obj = ChemicalMetadata(pubchemid=712, CAS=50000, formula='CH2O', MW=30.02598, smiles='C=O', InChI='CH2O/c1-2/h1H2', 
                    InChI_key='WSFSSNUMVMOOMR-UHFFFAOYSA-N', iupac_name='methanal', common_name='formaldehyde', 
                    synonyms=['methanal', 'formaldehyde', 'formalin', 'methanal', 'methylene oxide'])
    assert eval(obj.__repr__()) == obj
    assert hash(eval(obj.__repr__())) == hash(obj)

    # Test equality with identical object
    assert obj == obj
    
    # Test hash consistency
    assert hash(obj) == hash(obj)
    
    # Test deepcopy
    obj_copy = deepcopy(obj)
    assert obj == obj_copy
    assert hash(obj) == hash(obj_copy)
    
    # Test pickling
    obj_pickle = pickle.loads(pickle.dumps(obj))
    assert obj_pickle == obj
    assert hash(obj_pickle) == hash(obj)

    # Test hash consistency after accessing computed property
    hash_before = hash(obj)
    _ = obj.charge
    hash_after = hash(obj)
    assert hash_before == hash_after


def test_ChemicalMetadata_inequality():
    base_obj = ChemicalMetadata(
        pubchemid=712, 
        CAS=50000, 
        formula='CH2O', 
        MW=30.02598, 
        smiles='C=O', 
        InChI='CH2O/c1-2/h1H2', 
        InChI_key='WSFSSNUMVMOOMR-UHFFFAOYSA-N', 
        iupac_name='methanal', 
        common_name='formaldehyde', 
        synonyms=['methanal', 'formaldehyde']
    )
    
    # Test inequality with modified attributes
    modifications = {
        'pubchemid': 713,
        'CAS': 50001,
        'formula': 'CH3O',
        'MW': 30.02599,
        'smiles': 'CO',
        'InChI': 'different_inchi',
        'InChI_key': 'different_key',
        'iupac_name': 'different_name',
        'common_name': 'different_common_name',
        'synonyms': ['different', 'synonyms']
    }
    
    for attr, new_value in modifications.items():
        kwargs = {k: getattr(base_obj, k) for k in base_obj.__slots__ if k != '_charge'}
        kwargs[attr] = new_value
        different_obj = ChemicalMetadata(**kwargs)
        assert base_obj != different_obj
        assert hash(base_obj) != hash(different_obj)

def test_formula_search_exceptions():
    """
    Test that FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS contains exactly the 
    formulas that need special handling due to SMILES/formula interpretation conflicts.
    """
    actual_conflicts = set()
    for chemical in pubchem_db:
        if not chemical.formula:
            continue
        formula = chemical.formula
        
        try:
            # Skip single elements as they're handled by periodic table
            if formula in periodic_table:
                continue
            
            # Find both interpretations
            formula_result = pubchem_db.search_formula(formula)
            smiles_result = pubchem_db.search_smiles(formula)
            
            # If both exist and are different, we have a conflict
            if (formula_result and smiles_result and 
                formula_result.CAS != smiles_result.CAS):
                actual_conflicts.add(formula)
                
        except Exception:
            continue
    
    # Test that our exceptions list matches the actual conflicts
    assert FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS == actual_conflicts, \
        f"""FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS needs updating:
        Missing: {actual_conflicts - FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS}
        Unnecessary: {FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS - actual_conflicts}"""
    
    # Verify each exception actually needs special handling
    for formula in FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS:
        formula_result = pubchem_db.search_formula(formula)
        smiles_result = pubchem_db.search_smiles(formula)
        
        assert formula_result is not None, \
            f"Exception '{formula}' has no valid formula interpretation"
        assert smiles_result is not None, \
            f"Exception '{formula}' has no valid SMILES interpretation"
        assert formula_result.CAS != smiles_result.CAS, \
            f"Exception '{formula}' has same interpretation for formula and SMILES"

def test_Poling_databank_lookup_names():
    for CAS, name in zip(Cp_data_Poling.index.tolist(), Cp_data_Poling['Chemical'].tolist()):
        assert search_chemical(CAS).CASs == CAS
    for CAS, name in zip(Cp_data_Poling.index.tolist(), Cp_data_Poling['Chemical'].tolist()):
        assert search_chemical(name).CASs == CAS
