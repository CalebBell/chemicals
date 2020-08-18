# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

from __future__ import division

__all__ = ['checkCAS', 'CAS_from_any', 'search_chemical', 
           'mixture_from_any', 'cryogenics', 'dippr_compounds',
           'get_pubchem_db']

import os
from io import open
from chemicals.utils import PY37, source_path, os_path_join, can_load_data, to_num, CAS2int, int2CAS
from chemicals.elements import (periodic_table, homonuclear_elemental_gases,
                             charge_from_formula, serialize_formula)

folder = os_path_join(source_path, 'Identifiers')

def checkCAS(CASRN):
    """Checks if a CAS number is valid. Returns False if the parser cannot parse
    the given string.

    Parameters
    ----------
    CASRN : string
        A three-piece, dash-separated set of numbers

    Returns
    -------
    result : bool
        Boolean value if CASRN was valid. If parsing fails, return False also.

    Notes
    -----
    Check method is according to Chemical Abstract Society. However, no lookup
    to their service is performed; therefore, this function cannot detect
    false positives.

    Function also does not support additional separators, apart from '-'.
    
    CAS numbers up to the series 1 XXX XXX-XX-X are now being issued.
    
    A long can hold CAS numbers up to 2 147 483-64-7

    Examples
    --------
    >>> checkCAS('7732-18-5')
    True
    >>> checkCAS('77332-18-5')
    False
    """
    try:
        check = CASRN[-1] # Don't store the int - it is not necessary and is slower
        
        productsum = 0
        i = 1
        for num in CASRN.replace('-', '')[:-1][::-1]:
            productsum += i*int(num)
            i += 1
        return productsum % 10 == int(check)
    except:
        return False


class ChemicalMetadata(object):
    __slots__ = ['pubchemid', 'formula', 'MW', 'smiles', 'InChI', 'InChI_key',
                 'iupac_name', 'common_name', 'all_names', 'CAS', '_charge']
    def __repr__(self):
        return ('<ChemicalMetadata, name=%s, formula=%s, smiles=%s, MW=%g>'
                %(self.common_name, self.formula, self.smiles, self.MW))
        
    @property
    def charge(self):
        """Charge of the species as an integer.

        Computed as a property as most species do not have a charge and so
        storing it would be a waste of memory.
        """
        try:
            return self._charge
        except AttributeError:
            self._charge = charge_from_formula(self.formula)
            return self._charge
        
    @property
    def CASs(self):
        return int2CAS(self.CAS)
    
    def __init__(self, pubchemid, CAS, formula, MW, smiles, InChI, InChI_key,
                 iupac_name, common_name, all_names):
        self.pubchemid = pubchemid
        self.CAS = CAS
        self.formula = formula
        self.MW = MW
        self.smiles = smiles
        self.InChI = InChI
        
        self.InChI_key = InChI_key
        self.iupac_name = iupac_name
        self.common_name = common_name
        self.all_names = all_names
        

class ChemicalMetadataDB(object):
    loaded_main_db = False
    def __init__(self, 
                 elements=True,
                 main_db=os.path.join(folder, 'chemical identifiers pubchem large.tsv'),
                 user_dbs=[os.path.join(folder, 'chemical identifiers pubchem small.tsv'),
                           os.path.join(folder, 'chemical identifiers example user db.tsv'),
                           os.path.join(folder, 'Cation db.tsv'),
                           os.path.join(folder, 'Anion db.tsv'),
                           os.path.join(folder, 'Inorganic db.tsv')]):
        # TODO: delay creation of indexes, as most people won't be searching with all of them.
        
        
        self.pubchem_index = {}
        self.smiles_index = {}
        self.InChI_index = {}
        self.InChI_key_index = {}
        self.name_index = {}
        self.CAS_index = {}
        self.formula_index = {}

        self.main_db = main_db
        self.user_dbs = user_dbs
        self.elements = elements

#        self.load(self.main_db)
        for db in self.user_dbs:
            self.load(db)
        self.load_elements()
        
    def load_elements(self):
        if not self.elements:
            return None
        for ele in periodic_table:
            
            CAS = int(ele.CAS.replace('-', '')) # Store as int for easier lookup
            all_names = [ele.name.lower()]
            
            obj = ChemicalMetadata(pubchemid=ele.PubChem, CAS=CAS, 
                                   formula=ele.symbol, MW=ele.MW, smiles=ele.smiles,
                                   InChI=ele.InChI, InChI_key=ele.InChI_key,
                                   iupac_name=ele.name.lower(), 
                                   common_name=ele.name.lower(), 
                                   all_names=all_names)
            
            
            if ele.InChI_key in self.InChI_key_index:
                if ele.number not in homonuclear_elemental_gases:
                    obj_old = self.InChI_key_index[ele.InChI_key]
                    for name in obj_old.all_names:
                        self.name_index[name] = obj    
            
            self.InChI_key_index[ele.InChI_key] = obj
            self.CAS_index[CAS] = obj
            self.pubchem_index[ele.PubChem] = obj
            self.smiles_index[ele.smiles] = obj
            self.InChI_index[ele.InChI] = obj
            if ele.number in homonuclear_elemental_gases:
                for name in all_names:
                    self.name_index['monatomic ' + name] = obj    
            else:
                for name in all_names:
                    self.name_index[name] = obj    
            self.formula_index[obj.formula] = obj


    def load(self, file_name):
        f = open(file_name, encoding='utf-8')
        for line in f:
            # This is effectively the documentation for the file format of the file
            values = line.rstrip('\n').split('\t')
            (pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name) = values[0:9]
            CAS = int(CAS.replace('-', '')) # Store as int for easier lookup
            
            all_names = values[7:]
            pubchemid = int(pubchemid)

            obj = ChemicalMetadata(pubchemid, CAS, formula, float(MW), smiles,
                                    InChI, InChI_key, iupac_name, common_name, 
                                    all_names)
            
            # Lookup indexes
            self.CAS_index[CAS] = obj
            self.pubchem_index[pubchemid] = obj
            self.smiles_index[smiles] = obj
            self.InChI_index[InChI] = obj
            self.InChI_key_index[InChI_key] = obj
            for name in all_names:
                self.name_index[name] = obj
            self.formula_index[obj.formula] = obj
                    
        f.close()
        
    @property
    def can_autoload(self):
        return (not self.loaded_main_db and self.main_db is not None)
        
    def autoload_next(self):
        self.load(self.main_db)
        for db in self.user_dbs:
            self.load(db)
        self.load_elements()
        self.loaded_main_db = True
        return True
        
    def _search_autoload(self, identifier, index, autoload=True):
        if index:
            if identifier in index:
                return index[identifier]
            else:
                if autoload and self.can_autoload:
                    self.autoload_next()
                    return self._search_autoload(identifier, index, autoload)
        return False
    
    def search_pubchem(self, pubchem, autoload=True):
        if type(pubchem) != int:
            pubchem = int(pubchem)
        return self._search_autoload(pubchem, self.pubchem_index, autoload=autoload)
        
    def search_CAS(self, CAS, autoload=True):
        if type(CAS) != int:
            CAS = CAS2int(CAS)
        return self._search_autoload(CAS, self.CAS_index, autoload=autoload)

    def search_smiles(self, smiles, autoload=True):
        return self._search_autoload(smiles, self.smiles_index, autoload=autoload)

    def search_InChI(self, InChI, autoload=True):
        return self._search_autoload(InChI, self.InChI_index, autoload=autoload)

    def search_InChI_key(self, InChI_key, autoload=True):
        return self._search_autoload(InChI_key, self.InChI_key_index, autoload=autoload)

    def search_name(self, name, autoload=True):
        return self._search_autoload(name, self.name_index, autoload=autoload)
    
    def search_formula(self, formula, autoload=True):
        return self._search_autoload(formula, self.formula_index, autoload=autoload)


chemical_search_cache = {}

def CAS_from_any(ID, autoload=False, cache=True):
    """Wrapper around `search_chemical` which returns the CAS number of the
    found chemical directly.

    Parameters
    ----------
    ID : str
        One of the name formats described above

    Returns
    -------
    CASRN : string
        A three-piece, dash-separated set of numbers

    Notes
    -----
    An exception is raised if the name cannot be identified. The PubChem 
    database includes a wide variety of other synonyms, but these may not be
    present for all chemcials. See `search_chemical` for more details.

    Examples
    --------
    >>> CAS_from_any('water')
    '7732-18-5'
    >>> CAS_from_any('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
    '64-17-5'
    >>> CAS_from_any('CCCCCCCCCC')
    '124-18-5'
    >>> CAS_from_any('InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    '64-17-5'
    >>> CAS_from_any('pubchem=702')
    '64-17-5'
    >>> CAS_from_any('O') # only elements can be specified by symbol
    '17778-80-2'
    """
    return search_chemical(ID, autoload=False, cache=True).CASs

def search_chemical(ID, autoload=False, cache=True):
    """Looks up metadata about a chemical by searching and testing for the input
    string being any of the following types of chemical identifiers:

    * Name, in IUPAC form or common form or a synonym registered in PubChem
    * InChI name, prefixed by 'InChI=1S/' or 'InChI=1/'
    * InChI key, prefixed by 'InChIKey='
    * PubChem CID, prefixed by 'PubChem='
    * SMILES (prefix with 'SMILES=' to ensure smiles parsing; ex.
      'C' will return Carbon as it is an element whereas the SMILES 
      interpretation for 'C' is methane)
    * CAS number (obsolete numbers may point to the current number)    

    If the input is an ID representing an element, the following additional 
    inputs may be specified as 
    
    * Atomic symbol (ex 'Na')
    * Atomic number (as a string)

    Parameters
    ----------
    ID : str
        One of the name formats described above

    Returns
    -------
    chemical_metadata : ChemicalMetadata
        A class containing attributes which describe the chemical's metadata,
        [-]

    Notes
    -----
    An exception is raised if the name cannot be identified. The PubChem 
    database includes a wide variety of other synonyms, but these may not be
    present for all chemcials.

    Examples
    --------
    >>> search_chemical('water')
    <ChemicalMetadata, name=water, formula=H2O, smiles=O, MW=18.0153>
    >>> search_chemical('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> search_chemical('CCCCCCCCCC')
    <ChemicalMetadata, name=DECANE, formula=C10H22, smiles=CCCCCCCCCC, MW=142.286>
    >>> search_chemical('InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> search_chemical('pubchem=702')
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> search_chemical('O') # only elements can be specified by symbol
    <ChemicalMetadata, name=oxygen, formula=O, smiles=[O], MW=15.9994>
    """
    if cache and ID in chemical_search_cache:
        return chemical_search_cache[ID]
    if not _pubchem_db_loaded: get_pubchem_db()
    
    ID_arg = ID
    ID = ID.strip()
    ID_lower = ID.lower()
    if ID in periodic_table:
        if periodic_table[ID].number not in homonuclear_elemental_gases:
            return pubchem_db.search_CAS(periodic_table[ID].CAS)
        else:
            for i in [periodic_table._symbol_to_elements, 
                      periodic_table._number_to_elements,
                      periodic_table._CAS_to_elements]:
                if i == periodic_table._number_to_elements:
                    if int(ID in i):
                        obj = pubchem_db.search_CAS(periodic_table[int(ID)].CAS)
                        if cache:
                            chemical_search_cache[ID_arg] = obj
                        return obj
                else:
                    if ID in i:
                        obj = pubchem_db.search_CAS(periodic_table[ID].CAS)
                        if cache:
                            chemical_search_cache[ID_arg] = obj
                        return obj
    if checkCAS(ID):
        CAS_lookup = pubchem_db.search_CAS(ID, autoload)
        if CAS_lookup:
            if cache:
                chemical_search_cache[ID_arg] = CAS_lookup
            return CAS_lookup
        # handle the case of synonyms
        CAS_alternate_loopup = pubchem_db.search_name(ID, autoload)
        if CAS_alternate_loopup:
            if cache:
                chemical_search_cache[ID_arg] = CAS_alternate_loopup
            return CAS_alternate_loopup
            
        if not autoload:
            return search_chemical(ID, autoload=True)
        raise ValueError('A valid CAS number (%s) was recognized, but is not in the database' %(ID))
        
        
    
    ID_len = len(ID)
    if ID_len > 9:
        inchi_search = False
        # normal upper case is 'InChI=1S/'
        if ID_lower[0:9] == 'inchi=1s/':
            inchi_search = ID[9:]
        elif ID_lower[0:8] == 'inchi=1/':
            inchi_search = ID[8:]
        if inchi_search:
            inchi_lookup = pubchem_db.search_InChI(inchi_search, autoload)
            if inchi_lookup:
                if cache:
                    chemical_search_cache[ID_arg] = inchi_lookup
                return inchi_lookup
            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError('A valid InChI name (%s) was recognized, but it is not in the database' %(inchi_search))
        if ID_lower[0:9] == 'inchikey=':
            inchi_key_lookup = pubchem_db.search_InChI_key(ID[9:], autoload)
            if inchi_key_lookup:
                if cache:
                    chemical_search_cache[ID_arg] = inchi_key_lookup
                return inchi_key_lookup
            else:
                if not autoload:
                    obj = search_chemical(ID, autoload=True)
                    if cache:
                        chemical_search_cache[ID_arg] = obj
                    return obj
                raise ValueError('A valid InChI Key (%s) was recognized, but it is not in the database' %(inchi_key_lookup))
    if ID_len > 8:
        if ID_lower[0:8] == 'pubchem=':
            pubchem_lookup = pubchem_db.search_pubchem(ID[8:], autoload)
            if pubchem_lookup:
                if cache:
                    chemical_search_cache[ID_arg] = pubchem_lookup
                return pubchem_lookup
                
            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError('A PubChem integer (%s) identifier was recognized, but it is not in the database.' %(ID[8:]))
    if ID_len > 7:
        if ID_lower[0:7] == 'smiles=':
            smiles_lookup = pubchem_db.search_smiles(ID[7:], autoload)
            if smiles_lookup:
                if cache:
                    chemical_search_cache[ID_arg] = smiles_lookup
                return smiles_lookup
            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError('A SMILES identifier (%s) was recognized, but it is not in the database.' %(ID[7:]))

    # Try the smiles lookup anyway
    # Parsing SMILES is an option, but this is faster
    # Pybel API also prints messages to console on failure
    smiles_lookup = pubchem_db.search_smiles(ID, autoload)
    if smiles_lookup:
        if cache:
            chemical_search_cache[ID_arg] = smiles_lookup
        return smiles_lookup
    
    try:
        formula_query = pubchem_db.search_formula(serialize_formula(ID), autoload)
        if formula_query and type(formula_query) == ChemicalMetadata:
            if cache:
                chemical_search_cache[ID_arg] = formula_query
            return formula_query
    except:
        pass
    
    # Try a direct lookup with the name - the fastest
    name_lookup = pubchem_db.search_name(ID, autoload)
    if name_lookup:
        if cache:
            chemical_search_cache[ID_arg] = name_lookup
        return name_lookup
    
#     Permutate through various name options
    ID_no_space = ID.replace(' ', '')
    ID_no_space_dash = ID_no_space.replace('-', '')
    
    for name in [ID, ID_no_space, ID_no_space_dash]:
        for name2 in [name, name.lower()]:
            name_lookup = pubchem_db.search_name(name2, autoload)
            if name_lookup:
                if cache:
                    chemical_search_cache[ID_arg] = name_lookup
                return name_lookup
    
    if ID[-1] == ')' and '(' in ID:#
        # Try to matck in the form 'water (H2O)'
        first_identifier, second_identifier = ID[0:-1].split('(', 1)
        try:
            CAS1 = search_chemical(first_identifier)
            CAS2 = search_chemical(second_identifier)
            assert CAS1 == CAS2
            CAS = CAS1
            if cache:
                chemical_search_cache[ID_arg] = CAS
            return CAS
        except:
            pass
        
    if not autoload:
        return search_chemical(ID, autoload=True)
            
    raise ValueError('Chemical name (%s) not recognized' %(ID))





### DIPPR Database, chemical list only
# Obtained via the command:
# list(pd.read_excel('http://www.aiche.org/sites/default/files/docs/pages/sponsor_compound_list-2014.xlsx')['Unnamed: 2'])[2:]
# This is consistently faster than creating a list and then making the set.
def dippr_compounds():
    """Loads and returns a set of compounds known in the DIPPR database. This
    can be useful for knowing if a chemical is of industrial relevance.

    Returns
    -------
    dippr_compounds : set([str])
        A set of CAS numbers from the 2014 edition of the DIPPR database.
    """
    dippr_compounds = set()
    with open(os.path.join(folder, 'dippr_2014.csv')) as f:
        dippr_compounds.update(f.read().split('\n'))
    return dippr_compounds

class CommonMixtureMetadata(object):
    __slots__ = ['name', 'CASs', 'N', 'source', 'names', 'ws', 'zs',
                 'synonyms']
    def __repr__(self):
        return ('<MixtureMetadata, name=%s, N=%s, CASs=%s, ws=%s, zs=%s>'
                %(self.name, self.N, self.CASs, self.ws, self.zs))

    def __init__(self, name, CASs, N, source, names, ws, zs,
                 synonyms):
        self.name = name
        self.CASs = CASs
        self.N = N
        self.source = source
        self.names = names
        self.ws = ws
        self.zs = zs
        self.synonyms = synonyms


def mixture_from_any(ID):
    """Search by string for a mixture in the included common mixture database.
    The database primarily contains refrigerant blends. The variable
    `common_mixtures` contains all loaded entries.

    Parameters
    ----------
    ID : str
        A string or 1-element list containing the name which may represent a
        mixture.

    Returns
    -------
    mixture : CommonMixtureMetadata
        Object containing basic mixture information

    Notes
    -----
    White space, '-', and upper case letters are removed in the search.

    Examples
    --------
    >>> mixture_from_any('R512A')
    <MixtureMetadata, name=R512A, N=2, CASs=['811-97-2', '75-37-6'], ws=[0.05, 0.95], zs=[0.032949, 0.96705]>
    >>> mixture_from_any(['air'])
    <MixtureMetadata, name=Air, N=3, CASs=['7727-37-9', '7440-37-1', '7782-44-7'], ws=[0.7557, 0.0127, 0.2316], zs=[0.7812, 0.0092, 0.2096]>
    """
    if not mixture_composition_loaded:
        load_mixture_composition()
    if type(ID) == list:
        if len(ID) == 1:
            ID = ID[0]
        else:
            raise ValueError('If the input is a list, the list must contain only one item.')
    ID = ID.lower().strip()
    for i in (ID, ID.replace(' ', ''), ID.replace('-', '')):
        try:
            return common_mixtures_by_synonym[i]
        except KeyError:
            pass
    raise ValueError('Mixture name not recognized')


def IDs_to_CASs(IDs):
    CASs = None
    if hasattr(IDs, 'strip') or (isinstance(IDs, list) and len(IDs) == 1):
        try:
            # Assume the name was a pre-defined mixture
            mixname = mixture_from_any(IDs)
            _d = _MixtureDict[mixname]
            CASs = _d["CASs"]
        except:
            if hasattr(IDs, 'strip'):
                CASs = [IDs]
    if CASs is None:
        CASs = [CAS_from_any(ID) for ID in IDs]
    return CASs

cryogenics = {'132259-10-0': 'Air', '7440-37-1': 'Argon', '630-08-0':
'carbon monoxide', '7782-39-0': 'deuterium', '7782-41-4': 'fluorine',
'7440-59-7': 'helium', '1333-74-0': 'hydrogen', '7439-90-9': 'krypton',
'74-82-8': 'methane', '7440-01-9': 'neon', '7727-37-9': 'nitrogen',
'7782-44-7': 'oxygen', '7440-63-3': 'xenon'}

_pubchem_db_loaded = False
def get_pubchem_db():
    """Helper function to delay the creation of the pubchem_db object.

    This avoids loading the database when it is not needed.
    """
    global _pubchem_db_loaded, pubchem_db
    if _pubchem_db_loaded:
        return pubchem_db
    else:
        pubchem_db = ChemicalMetadataDB()
    _pubchem_db_loaded = True
    return pubchem_db

mixture_composition_loaded = False
global common_mixtures_by_synonym, common_mixtures

def load_mixture_composition():
    global mixture_composition_loaded, common_mixtures_by_synonym, common_mixtures
    common_mixtures = {}
    common_mixtures_by_synonym = {}
    with open(os.path.join(folder, 'Mixtures Compositions.tsv')) as f:
        """Read in a dict of 90 or so mixutres, their components, and synonyms.

        Small errors in mole fractions not adding to 1 are known. Errors in
        adding mass fraction are less common, present at the 5th decimal. Mass
        basis is assumed for all mixtures.
        """
        next(f)
        for line in f:
            values = to_num(line.strip('\n').strip('\t').split('\t'))
            name, source, N = values[0:3]
            N = int(N)
            CASs, names, ws, zs = values[3:3+N], values[3+N:3+2*N], values[3+2*N:3+3*N], values[3+3*N:3+4*N]
            synonyms = values[3+4*N:]
            if synonyms:
                synonyms = [i.lower() for i in synonyms]
            synonyms.append(name.lower())
            obj = CommonMixtureMetadata(name=name, CASs=CASs, N=N, source=source,
                                        names=names, ws=ws, zs=zs, synonyms=synonyms)
            common_mixtures[name] = obj
            
            for syn in synonyms:
                common_mixtures_by_synonym[syn] = obj
    mixture_composition_loaded = True


if PY37:
    def __getattr__(name):
        if name == 'pubchem_db':
            return get_pubchem_db()
        elif name == 'common_mixtures' or name == 'common_mixtures_by_synonym':
            load_mixture_composition()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        get_pubchem_db()
        load_mixture_composition()