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
from fluids.numerics import assert_close

from chemicals.elements import molecular_weight, nested_formula_parser, periodic_table, serialize_formula
from chemicals.identifiers import (
    CAS_from_any,
    CAS_to_int,
    ChemicalMetadataDB,
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


@pytest.mark.slow
@pytest.mark.online
def test_dippr_2016_matched_meta():
    df2 = pd.read_excel('https://www.aiche.org/sites/default/files/docs/pages/dippr_compound_list_2016.xlsx')
    names = df2['Name'].tolist()
    CASs = df2['CASN'].tolist()
    for i, CAS in enumerate(CASs):
        if CAS in ('16462-44-5', '75899-69-3'):
            # CELLOBIOSE (not the latest CAS, now is 528-50-7)
            # TRIPROPYLENE GLYCOL MONOETHYL ETHER (cannot find structure)
            continue
        try:
            if isnan(CAS):
                continue
        except:
            pass
        assert CAS_from_any(CAS) == CAS


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
    for i in pubchem_db.CAS_index.values():
        assert check_CAS(i.CASs)

@pytest.mark.skip
@pytest.mark.xfail
def test_database_formulas():
    # Failures are thing slike 3He, C2D4Br2, C14H18N3NaO10[99Tc], [1H]I
    # The fix here is adding an isotope db and making the formula parser handle isotopes as well.
    # This worked until isotopes were added to formulas
    for i in pubchem_db.CAS_index.values():
        assert i.formula == serialize_formula(i.formula)

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
    for i in db.CAS_index.values():
        assert i.formula == serialize_formula(i.formula)

    # Check CAS validity
    for i in db.CAS_index.values():
        assert check_CAS(i.CASs)

    # MW checker
    for i in db.CAS_index.values():
        formula = serialize_formula(i.formula)
        atoms = nested_formula_parser(formula, check=False)
        mw_calc = molecular_weight(atoms)
        assert_close(mw_calc, i.MW, atol=0.05)


    for CAS, d in db.CAS_index.items():
        if d.InChI:
            assert CAS_from_any('InChI=1S/' + d.InChI) == int_to_CAS(CAS)

    for CAS, d in db.CAS_index.items():
        if d.InChI_key:
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
        if formula in {'H2MgO2', 'F2N2'}:
            # Formulas which are not unique by design
            continue
        assert CAS_from_any(formula) == d.CASs

    # Check smiles are unique / can lookup by smiles
    for smi, d in db.smiles_index.items():
        if not smi:
            continue
        assert CAS_from_any('smiles=' + smi) == d.CASs

    # Check formula is formatted right
    for i in db.CAS_index.values():
        assert i.formula == serialize_formula(i.formula)

    # Check CAS validity
    for i in db.CAS_index.values():
        assert check_CAS(i.CASs)

    # MW checker
    for i in db.CAS_index.values():
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
    assert hit0 is hit1

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
    assert search_chemical('15922-78-8').common_name == "sodium pyrithione"
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
    assert "n-Octatriacontane" in search_chemical('7194-85-6').synonyms

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
    assert "Parafol 22-95" in search_chemical('629-97-0').synonyms

    # CAS: 57-48-7 (fructose)
    assert "D-Fructose" not in search_chemical('57-48-7').synonyms
    assert "d-Fructose" in search_chemical('57-48-7').synonyms

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
    assert "Chromium chloride (CrCl3)" in search_chemical('10025-73-7').synonyms

    # CAS: 10031-25-1 (chromium bromide (CrBr3))
    assert "UZDWIWGMKWZEPE-UHFFFAOYSA-K" not in search_chemical('10031-25-1').synonyms
    assert "Chromium bromide (CrBr3)" in search_chemical('10031-25-1').synonyms

    # CAS: 10043-01-3 (aluminum sulfate)
    assert "ALUMINUM SULFATE" not in search_chemical('10043-01-3').synonyms
    assert "Aluminum sulfate" in search_chemical('10043-01-3').synonyms

    # CAS: 10102-44-0 (nitrogen dioxide)
    assert "nitrogen oxide (NO2)" in search_chemical('10102-44-0').synonyms

    # CAS: 10257-55-3 (calcium sulfite)
    assert "CALCIUM SULFITE (1:1)" not in search_chemical('10257-55-3').synonyms
    assert "sulfurous acid, calcium salt (1:1)" in search_chemical('10257-55-3').synonyms

    # CAS: 12045-63-5 (titanium diboride)
    assert "Titanium boride (TiB2)" not in search_chemical('12045-63-5').synonyms
    assert "Titanium diboride" in search_chemical('12045-63-5').synonyms

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
    assert {'Chlorine', 'Fluorine', 'Hydrogen', 'Nitrogen', 'Oxygen', 'Bromine', 'Iodine'} == set(failed_CASs)



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
            assert CAS_from_any('PubChem=' + str(d.PubChem)) == CAS
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

