"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

This module contains a database of metadata on ~70000 chemicals from the PubChem
datase. It contains comprehensive feature for searching the metadata.
It also includes a small database of common mixture compositions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Search Functions
----------------
.. autofunction:: chemicals.identifiers.CAS_from_any
.. autofunction:: chemicals.identifiers.MW
.. autofunction:: chemicals.identifiers.search_chemical
.. autofunction:: chemicals.identifiers.IDs_to_CASs

CAS Number Utilities
--------------------
.. autofunction:: chemicals.identifiers.check_CAS
.. autofunction:: chemicals.identifiers.CAS_to_int
.. autofunction:: chemicals.identifiers.int_to_CAS
.. autofunction:: chemicals.identifiers.sorted_CAS_key

Database Objects
----------------
There is an object used to represent a chemical's metadata, an object used to
represent a common mixture's composition, and an object used to hold the
mixture metadata.

.. autoclass:: chemicals.identifiers.ChemicalMetadata
.. autoclass:: chemicals.identifiers.CommonMixtureMetadata
.. autoclass:: chemicals.identifiers.ChemicalMetadataDB
.. autofunction:: chemicals.identifiers.get_pubchem_db

Chemical Groups
---------------
It is convenient to tag some chemicals with labels like "refrigerant", or in
a certain database or not. The following chemical groups are available.

.. autodata:: chemicals.identifiers.cryogenics
.. autodata:: chemicals.identifiers.common_commercial_gases
.. autofunction:: chemicals.identifiers.dippr_compounds
"""


__all__ = ['check_CAS', 'CAS_from_any', 'MW', 'search_chemical',
           'mixture_from_any', 'cryogenics', 'common_commercial_gases', 'dippr_compounds', 'IDs_to_CASs',
           'get_pubchem_db', 'CAS_to_int', 'sorted_CAS_key', 'int_to_CAS']

import os

from chemicals.elements import charge_from_formula, homonuclear_elements_CASs_set, periodic_table, serialize_formula, nested_formula_parser
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path, to_num

folder = os_path_join(source_path, 'Identifiers')

@mark_numba_incompatible
def check_CAS(CASRN):
    """Checks if a CAS number is valid. Returns False if the parser cannot parse
    the given string.

    Parameters
    ----------
    CASRN : str
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
    >>> check_CAS('7732-18-5')
    True
    >>> check_CAS('77332-18-5')
    False
    """
    try:
        if CASRN.count('-') != 2:
            return False
        if CASRN[-2] != '-' or CASRN[-5] != '-':
            return False
        check = CASRN[-1] # Don't store the int - it is not necessary and is slower

        productsum = 0
        chars = CASRN.replace('-', '')[:-1]
        i = len(chars)
        for num in chars:
            productsum += i*int(num)
            i -= 1
        return productsum % 10 == int(check)
    except:
        return False

@mark_numba_incompatible
def CAS_to_int(CASRN):
    r'''Converts CAS number of a compounds from a string to an int. This is
    helpful when storing large amounts of CAS numbers, as their strings take up
    more memory than their numerical representational. All CAS numbers fit into
    64 bit ints.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    CASRN : int
        CASRN [-]

    Notes
    -----
    Accomplishes conversion by removing dashes only, and then converting to an
    int. An incorrect CAS number will change without exception.

    Examples
    --------
    >>> CAS_to_int('7704-34-9')
    7704349
    '''
    return int(CASRN.replace('-', ''))

@mark_numba_incompatible
def int_to_CAS(CASRN):
    r'''Converts CAS number of a compounds from an int to an string. This is
    helpful when dealing with int CAS numbers.

    Parameters
    ----------
    CASRN : int
        CASRN [-]

    Returns
    -------
    CASRN : str
        CASRN [-]

    Notes
    -----
    Handles CAS numbers with an unspecified number of digits. Does not work on
    floats.

    Examples
    --------
    >>> int_to_CAS(7704349)
    '7704-34-9'
    '''
    CASRN = str(CASRN)
    return CASRN[:-3]+'-'+CASRN[-3:-1]+'-'+CASRN[-1]

@mark_numba_incompatible
def sorted_CAS_key(CASs):
    r'''Takes a list of CAS numbers as strings, and returns a tuple of the same
    CAS numbers, sorted from smallest to largest. This is very convenient for
    obtaining a unique hash of a set of compounds, so as to see if two
    groups of compounds are the same.

    Parameters
    ----------
    CASs : list[str]
        CAS numbers as strings [-]

    Returns
    -------
    CASs_sorted : tuple[str]
        Sorted CAS numbers from lowest (first) to highest (last) [-]

    Notes
    -----
    Does not check CAS numbers for validity.

    Examples
    --------
    >>> sorted_CAS_key(['7732-18-5', '64-17-5', '108-88-3', '98-00-0'])
    ('64-17-5', '98-00-0', '108-88-3', '7732-18-5')
    '''
    int_CASs = [CAS_to_int(i) for i in CASs]
    return tuple(CAS for _, CAS in sorted(zip(int_CASs, CASs)))

class ChemicalMetadata:
    """Class for storing metadata on chemicals.

    Attributes
    ----------
    pubchemid : int
        Identification number on pubchem database; access their information
        online at https://pubchem.ncbi.nlm.nih.gov/compound/<pubchemid>
        [-]
    formula : str
        Formula of the compound; in the same format as
        :obj:`chemicals.elements.serialize_formula` generates, [-]
    MW : float
        Molecular weight of the compound as calculated with the standard
        atomic abundances; consistent with the element weights in
        :obj:`chemicals.elements.periodic_table`, [g/mol]
    smiles : str
        SMILES identification string, [-]
    InChI : str
        InChI identification string as given in pubchem (there can be multiple
        valid InChI strings for a compound), [-]
    InChI_key : str
        InChI key identification string (meant to be unique to a compound), [-]
    iupac_name : str
        IUPAC name as given in pubchem, [-]
    common_name : str
        Common name as given in pubchem, [-]
    synonyms : list[str]
        List of synonyms of the compound, [-]
    CAS : int
        CAS number of the compound; stored as an int for memory efficiency, [-]
    """

    __slots__ = ('pubchemid', 'formula', 'MW', 'smiles', 'InChI', 'InChI_key',
                 'iupac_name', 'common_name', 'synonyms', 'CAS', '_charge')

    def __repr__(self):
        return (f"ChemicalMetadata(pubchemid={self.pubchemid!r}, "
                f"CAS={self.CAS!r}, "
                f"formula={self.formula!r}, "
                f"MW={self.MW!r}, "
                f"smiles={self.smiles!r}, "
                f"InChI={self.InChI!r}, "
                f"InChI_key={self.InChI_key!r}, "
                f"iupac_name={self.iupac_name!r}, "
                f"common_name={self.common_name!r}, "
                f"synonyms={self.synonyms!r})")

    def __eq__(self, other):
        if self is other:
            return True
        return (self.pubchemid == other.pubchemid and
                self.CAS == other.CAS and
                self.formula == other.formula and
                self.MW == other.MW and
                self.smiles == other.smiles and
                self.InChI == other.InChI and
                self.InChI_key == other.InChI_key and
                self.iupac_name == other.iupac_name and
                self.common_name == other.common_name and
                self.synonyms == other.synonyms)

    def __hash__(self):
        return hash((self.pubchemid,
                    self.CAS,
                    self.formula,
                    self.MW,
                    self.smiles, 
                    self.InChI,
                    self.InChI_key,
                    self.iupac_name,
                    self.common_name,
                    tuple(self.synonyms)))

    def __str__(self):
        return (f'<ChemicalMetadata, name={self.common_name}, formula={self.formula}, smiles={self.smiles}, MW={self.MW:g}>')

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
        """CAs number of the compound as a string.
        """
        return int_to_CAS(self.CAS)

    def __init__(self, pubchemid, CAS, formula, MW, smiles, InChI, InChI_key,
                 iupac_name, common_name, synonyms):
        self.pubchemid = pubchemid
        self.CAS = CAS
        self.formula = formula
        self.MW = MW
        self.smiles = smiles
        self.InChI = InChI

        self.InChI_key = InChI_key
        self.iupac_name = iupac_name
        self.common_name = common_name
        self.synonyms = synonyms

PUBCHEM_LARGE_DB_NAME = 'chemical identifiers pubchem large.tsv'
PUBCHEM_SMALL_DB_NAME = 'chemical identifiers pubchem small.tsv'
PUBCHEM_EXAMPLE_DB_NAME = 'chemical identifiers example user db.tsv'
PUBCHEM_CATION_DB_NAME = 'Cation db.tsv'
PUBCHEM_ANION_DB_NAME = 'Anion db.tsv'
PUBCHEM_IONORGANIC_DB_NAME = 'Inorganic db.tsv'

# Preference file constants
ANION_PREFERENCES_FILE = 'anion_preferences.json'
CATION_PREFERENCES_FILE = 'cation_preferences.json'
INORGANIC_PREFERENCES_FILE = 'inorganic_preferences.json'
ORGANIC_PREFERENCES_FILE = 'organic_preferences.json'

def load_chemical_preferences():
    """Loads preferred and unpreferred CAS numbers from preference files in the 
    Identifiers folder.
        
    Returns
    -------
    preferred_CAS : set
        Set of preferred CAS numbers from all loaded files
    unpreferred_CAS : set 
        Set of unpreferred CAS numbers from all loaded files
        
    Notes
    -----

    Each file should contain 'preferred_cas' and 'unpreferred_cas' lists.
    Missing files are skipped silently.
    
    Examples
    --------
    >>> preferred, unpreferred = load_chemical_preferences()
    """
    import json
    preferred_CAS = set()
    unpreferred_CAS = set()
    
    # Files to load
    files = [
        ANION_PREFERENCES_FILE,
        CATION_PREFERENCES_FILE, 
        INORGANIC_PREFERENCES_FILE,
        ORGANIC_PREFERENCES_FILE
    ]
    
    # Load each file and update the sets
    for filename in files:
        pref_file = os.path.join(folder, filename)
        if os.path.exists(pref_file):
            with open(pref_file) as f:
                preferences = json.load(f)
                preferred_CAS.update(preferences.get('preferred_cas', []))
                unpreferred_CAS.update(preferences.get('unpreferred_cas', []))
    
    return preferred_CAS, unpreferred_CAS


class ChemicalMetadataDB:
    '''Object which holds the main database of chemical metadata.

    .. warning:: To allow the `chemicals` to grow and improve, the details of
       this class may change in the future without notice!

    '''

    loaded_main_db = False
    def __init__(self,
                 elements=True,
                 main_db=os_path_join(folder, PUBCHEM_LARGE_DB_NAME),
                 user_dbs=[os_path_join(folder, PUBCHEM_SMALL_DB_NAME),
                           os_path_join(folder, PUBCHEM_EXAMPLE_DB_NAME),
                           os_path_join(folder, PUBCHEM_CATION_DB_NAME),
                           os_path_join(folder, PUBCHEM_ANION_DB_NAME),
                           os_path_join(folder, PUBCHEM_IONORGANIC_DB_NAME)]):
                    
        '''Construct the database from its parameters, loading all of the files in
        `user_dbs`, the periodic table, and defering loading of `main_db`
        as it is very large until a search doesn't find a chemical in the smaller
        database.
        '''
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

        for db in self.user_dbs:
            self.load(db)
        self.load_elements()

    def load_elements(self):
        '''Load elements into the indexes.
        '''
        if not self.elements:
            return None

        InChI_key_index, CAS_index, pubchem_index = self.InChI_key_index, self.CAS_index, self.pubchem_index
        smiles_index, InChI_index, formula_index = self.smiles_index, self.InChI_index, self.formula_index
        name_index = self.name_index

        for ele in periodic_table:
            if ele.CAS in homonuclear_elements_CASs_set:
                continue
            CAS = int(ele.CAS.replace('-', '')) # Store as int for easier lookup
            ele_lower_name = ele.name.lower()
            obj_old = CAS_index.get(CAS, None)
            synonyms = [ele_lower_name]

            obj = ChemicalMetadata(pubchemid=ele.PubChem, CAS=CAS,
                                   formula=ele.symbol, MW=ele.MW, smiles=ele.smiles,
                                   InChI=ele.InChI, InChI_key=ele.InChI_key,
                                   iupac_name=ele_lower_name,
                                   common_name=ele_lower_name,
                                   synonyms=synonyms)
            if obj_old is not None:
                synonyms = obj_old.synonyms
                obj.synonyms = synonyms

            InChI_key_index[obj.InChI_key] = obj
            CAS_index[obj.CAS] = obj
            pubchem_index[obj.pubchemid] = obj
            smiles_index[obj.smiles] = obj
            InChI_index[obj.InChI] = obj
            for name in synonyms:
                name_index[name] = obj
            formula_index[obj.formula] = obj


    def load(self, file_name):
        '''Load a particular file into the indexes.
        '''
        f = open(file_name, encoding='utf-8')
        for line in f:
            # This is effectively the documentation for the file format of the file
            values = line.rstrip('\n').split('\t')
            (pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name) = values[0:9]
            CAS = int(CAS.replace('-', '')) # Store as int for easier lookup

            synonyms = values[7:] # include iupac name and common name in this without having to duplicate the names in the file
            pubchemid = int(pubchemid)

            obj = ChemicalMetadata(pubchemid, CAS, formula, float(MW), smiles,
                                    InChI, InChI_key, iupac_name, common_name,
                                    synonyms)

            # Lookup indexes
            self.CAS_index[CAS] = obj
            self.pubchem_index[pubchemid] = obj
            self.smiles_index[smiles] = obj
            self.InChI_index[InChI] = obj
            self.InChI_key_index[InChI_key] = obj
            for name in synonyms:
                self.name_index[name] = obj
                self.name_index[name.lower()] = obj
            self.formula_index[obj.formula] = obj

        f.close()

    def __iter__(self):
        if not self.finished_loading:
            self.autoload_main_db()
        return iter(i for i in self.CAS_index.values())

    @property
    def finished_loading(self):
        '''Whether or not the database has loaded the main database.
        '''
        return not (not self.loaded_main_db and self.main_db is not None)

    def finish_loading(self):
        '''Complete loading the main database, if it has not been fully loaded.
        '''
        if not self.finished_loading:
            self.autoload_main_db()

    def autoload_main_db(self):
        '''Load the main database when needed.
        '''
        if os.path.exists(self.main_db):
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
                if autoload and not self.finished_loading:
                    self.autoload_main_db()
                    return self._search_autoload(identifier, index, autoload)
        return False

    def search_pubchem(self, pubchem, autoload=True):
        '''Search for a chemical by its pubchem number. Accepts strings or ints.
        '''
        return self._search_autoload(int(pubchem), self.pubchem_index, autoload=autoload)

    def search_CAS(self, CAS, autoload=True):
        '''Search for a chemical by its CAS number. Accepts strings or ints.
        '''
        if type(CAS) != int:
            CAS = CAS_to_int(CAS)
        return self._search_autoload(CAS, self.CAS_index, autoload=autoload)

    def search_smiles(self, smiles, autoload=True):
        '''Search for a chemical by its smiles string.
        '''
        return self._search_autoload(smiles, self.smiles_index, autoload=autoload)

    def search_InChI(self, InChI, autoload=True):
        '''Search for a chemical by its InChI string.
        '''
        return self._search_autoload(InChI, self.InChI_index, autoload=autoload)

    def search_InChI_key(self, InChI_key, autoload=True):
        '''Search for a chemical by its InChI key.
        '''
        return self._search_autoload(InChI_key, self.InChI_key_index, autoload=autoload)

    def search_name(self, name, autoload=True):
        '''Search for a chemical by its name.
        '''
        return self._search_autoload(name, self.name_index, autoload=autoload)

    def search_formula(self, formula, autoload=True):
        '''Search for a chemical by its serialized formula.
        '''
        return self._search_autoload(formula, self.formula_index, autoload=autoload)
import sqlite3
class ChemicalMetadataDiskDB:
    """SQLite-backed version of ChemicalMetadataDB with preferred ordering support"""
    
    def __init__(self, db_path=os_path_join(folder, 'metadata.db')):
        """Initialize connection to the SQLite database
        
        Parameters
        ----------
        db_path : str or Path
            Path to the SQLite database file
        """
        self.db_path = db_path
        self._conn = sqlite3.connect(db_path)
        self._conn.row_factory = sqlite3.Row  # Allow column name access
        
    def _row_to_metadata(self, row):
        """Convert a database row to a ChemicalMetadata object"""
        if row is None:
            return None
        
        synonyms = [row['iupac_name'], row['common_name']]
        if row['raw_synonyms']:
            synonyms.extend(row['raw_synonyms'].split('\t'))
        
        return ChemicalMetadata(
            pubchemid=row['pubchemid'],
            CAS=row['cas'],
            formula=row['formula'],
            MW=row['mw'],
            smiles=row['smiles'],
            InChI=row['inchi'],
            InChI_key=row['inchi_key'],
            iupac_name=row['iupac_name'],
            common_name=row['common_name'],
            synonyms=synonyms
        )
    
    def search_pubchem(self, pubchem, autoload=True):
        """Search for a chemical by its pubchem number"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE pubchemid = ? ORDER BY preferred DESC LIMIT 1",
            (int(pubchem),)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_CAS(self, CAS, autoload=True):
        """Search for a chemical by its CAS number"""
        if isinstance(CAS, str):
            CAS = int(CAS.replace('-', ''))
        
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE cas = ? ORDER BY preferred DESC LIMIT 1",
            (CAS,)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_smiles(self, smiles, autoload=True):
        """Search for a chemical by its SMILES string"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE smiles = ? ORDER BY preferred DESC LIMIT 1",
            (smiles,)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_InChI(self, InChI, autoload=True):
        """Search for a chemical by its InChI string"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE inchi = ? ORDER BY preferred DESC LIMIT 1",
            (InChI,)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_InChI_key(self, InChI_key, autoload=True):
        """Search for a chemical by its InChI key"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE inchi_key = ? ORDER BY preferred DESC LIMIT 1",
            (InChI_key,)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_name(self, name, autoload=True):
        """Search for a chemical by its name"""
        cur = self._conn.cursor()
        cur.execute("""
            SELECT c.* FROM chemicals c
            JOIN chemical_synonyms cs ON c.cas = cs.cas
            WHERE cs.synonym = ?
            ORDER BY c.preferred DESC LIMIT 1
        """, (name,))
        return self._row_to_metadata(cur.fetchone())
    
    def search_formula(self, formula, autoload=True):
        """Search for a chemical by its formula"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE formula = ? ORDER BY preferred DESC LIMIT 1",
            (formula,)
        )
        return self._row_to_metadata(cur.fetchone())
    
    def search_pubchem(self, pubchem, autoload=True):
        """Search for a chemical by its pubchem number"""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT * FROM chemicals WHERE pubchemid = ? ORDER BY preferred DESC LIMIT 1",
            (int(pubchem),)
        )
        return self._row_to_metadata(cur.fetchone())
            

    def __iter__(self):
        """Iterate over all chemicals in the database"""
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM chemicals ORDER BY preferred DESC")
        while True:
            batch = cur.fetchmany(1000)  # Process in batches for memory efficiency
            if not batch:
                break
            for row in batch:
                yield self._row_to_metadata(row)
    
    def __len__(self):
        """Return the total number of chemicals in the database"""
        cur = self._conn.cursor()
        cur.execute("SELECT COUNT(*) FROM chemicals")
        return cur.fetchone()[0]

    @property
    def CAS_index(self):
        """Build and return a dictionary mapping CAS numbers to ChemicalMetadata objects.
        
        Returns
        -------
        dict
            Dictionary with CAS numbers (integers) as keys and ChemicalMetadata objects as values
        """
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM chemicals ORDER BY preferred DESC")
        cas_dict = {}
        for row in cur:
            metadata = self._row_to_metadata(row)
            if metadata.CAS not in cas_dict:  # Only keep the most preferred entry for each CAS
                cas_dict[metadata.CAS] = metadata
        return cas_dict


    @property
    def smiles_index(self):
        """Build and return a dictionary mapping SMILES strings to ChemicalMetadata objects.
        
        Returns
        -------
        dict
            Dictionary with SMILES strings as keys and ChemicalMetadata objects as values
        """
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM chemicals ORDER BY preferred DESC")
        smiles_dict = {}
        for row in cur:
            metadata = self._row_to_metadata(row)
            if metadata.smiles and metadata.smiles not in smiles_dict:  # Only keep the most preferred entry
                smiles_dict[metadata.smiles] = metadata
        return smiles_dict

    @property
    def InChI_index(self):
        """Build and return a dictionary mapping InChI strings to ChemicalMetadata objects.
        
        Returns
        -------
        dict
            Dictionary with InChI strings as keys and ChemicalMetadata objects as values
        """
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM chemicals ORDER BY preferred DESC")
        inchi_dict = {}
        for row in cur:
            metadata = self._row_to_metadata(row)
            if metadata.InChI and metadata.InChI not in inchi_dict:  # Only keep the most preferred entry
                inchi_dict[metadata.InChI] = metadata
        return inchi_dict

    @property
    def InChI_key_index(self):
        """Build and return a dictionary mapping InChI keys to ChemicalMetadata objects.
        
        Returns
        -------
        dict
            Dictionary with InChI keys as keys and ChemicalMetadata objects as values
        """
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM chemicals ORDER BY preferred DESC")
        inchi_key_dict = {}
        for row in cur:
            metadata = self._row_to_metadata(row)
            if metadata.InChI_key and metadata.InChI_key not in inchi_key_dict:  # Only keep the most preferred entry
                inchi_key_dict[metadata.InChI_key] = metadata
        return inchi_key_dict

    @property
    def finished_loading(self):
        """Always returns True as database is pre-loaded"""
        return True

    def close(self):
        """Explicitly close the database connection"""
        self._conn.close()

@mark_numba_incompatible
def CAS_from_any(ID, autoload=False, cache=True):
    """Wrapper around `search_chemical` which returns the CAS number of the
    found chemical directly.

    Parameters
    ----------
    ID : str
        One of the name formats described by `search_chemical`, [-]
    autoload : bool, optional
        Whether to load new chemical databanks during the search if a hit is
        not immediately found, [-]
    cache : bool, optional
        Whether or not to cache the search for faster lookup in subsequent
        queries, [-]

    Returns
    -------
    CASRN : str
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
    return search_chemical(ID, autoload=autoload, cache=cache).CASs

@mark_numba_incompatible
def MW(ID, autoload=False, cache=True):
    """Wrapper around `search_chemical` which returns the molecular weight of the
    found chemical directly.

    Parameters
    ----------
    ID : str
        One of the name formats described by `search_chemical`
    autoload : bool, optional
        Whether to load new chemical databanks during the search if a hit is
        not immediately found, [-]
    cache : bool, optional
        Whether or not to cache the search for faster lookup in subsequent
        queries, [-]

    Returns
    -------
    MW : float
        Molecular weight of chemical, [g/mol]

    Notes
    -----
    An exception is raised if the name cannot be identified. The PubChem
    database includes a wide variety of other synonyms, but these may not be
    present for all chemcials. See `search_chemical` for more details.

    Examples
    --------
    >>> MW('water')
    18.01528
    >>> MW('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
    46.06844
    >>> MW('CCCCCCCCCC')
    142.28168
    >>> MW('InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    46.06844
    >>> MW('pubchem=702')
    46.06844
    >>> MW('O') # only elements can be specified by symbol
    15.9994

    """
    return search_chemical(ID, autoload=autoload, cache=cache).MW

chemical_search_cache = {}
chemical_search_cache_max_size = 200

FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS = {
    'CO',  # Carbon monoxide vs methanol
    'CS',  # Carbon monosulfide vs methanethiol
    'NO',  # Nitric oxide vs hydroxylamine
    'CNO'  # NCO radical vs n-methylhydroxylamine
}
@mark_numba_incompatible
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
    autoload : bool, optional
        Whether to load new chemical databanks during the search if a hit is
        not immediately found, [-]
    cache : bool, optional
        Whether or not to cache the search for faster lookup in subsequent
        queries, [-]

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
    >>> print(search_chemical('water'))
    <ChemicalMetadata, name=water, formula=H2O, smiles=O, MW=18.0153>
    >>> print(search_chemical('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'))
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> print(search_chemical('CCCCCCCCCC'))
    <ChemicalMetadata, name=decane, formula=C10H22, smiles=CCCCCCCCCC, MW=142.282>
    >>> print(search_chemical('InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N'))
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> print(search_chemical('pubchem=702'))
    <ChemicalMetadata, name=ethanol, formula=C2H6O, smiles=CCO, MW=46.0684>
    >>> print(search_chemical('O')) # only elements can be specified by symbol
    <ChemicalMetadata, name=atomic oxygen, formula=O, smiles=[O], MW=15.9994>
    """
    if cache and ID in chemical_search_cache:
        return chemical_search_cache[ID]
    if not _pubchem_db_loaded: get_pubchem_db()  # pragma: no cover
    hit = _search_chemical(ID, autoload)
    if cache:
        if len(chemical_search_cache) > chemical_search_cache_max_size:
            # invalidate cache by time - first entry is removed relying on
            # dict ordering new in Python 3.7
            chemical_search_cache.pop(next(chemical_search_cache.keys().__iter__()))
        chemical_search_cache[ID] = hit
    return hit

def _search_chemical(ID, autoload):
    ID_arg = ID
    ID = ID.strip()
    ID_lower = ID.lower()
    if ID in periodic_table:
        """Special handling for homonuclear elements. Search '1'> H, 'H'> H, monotomic CAS > H
        but "Hydrogen"> H2.
        pubchem_db does not contain atomic numbers, so searching in the periodic table is necessary.
        """
        if (ID in periodic_table._symbol_to_elements or ID in periodic_table._number_to_elements
            or ID in periodic_table._CAS_to_elements):
            obj = pubchem_db.search_CAS(periodic_table[ID].CAS)
        else:
            obj = pubchem_db.search_CAS(periodic_table[ID].CAS_standard)
        return obj
    if check_CAS(ID):
        CAS_lookup = pubchem_db.search_CAS(ID, autoload)
        if CAS_lookup:
            return CAS_lookup
        # handle the case of synonyms
        CAS_alternate_loopup = pubchem_db.search_name(ID, autoload)
        if CAS_alternate_loopup:
            return CAS_alternate_loopup

        if not autoload:
            return search_chemical(ID, autoload=True)
        raise ValueError(f'A valid CAS number ({ID}) was recognized, but is not in the database')



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
                return inchi_lookup
            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError(f'A valid InChI name ({inchi_search}) was recognized, but it is not in the database')
        if ID_lower[0:9] == 'inchikey=':
            inchi_key_lookup = pubchem_db.search_InChI_key(ID[9:], autoload)
            if inchi_key_lookup:
                return inchi_key_lookup
            else:
                if not autoload:
                    obj = search_chemical(ID, autoload=True)
                    return obj
                raise ValueError(f'A valid InChI Key ({ID[9:]}) was recognized, but it is not in the database')
    if ID_len > 8:
        if ID_lower[0:8] == 'pubchem=':
            pubchem_lookup = pubchem_db.search_pubchem(ID[8:], autoload)
            if pubchem_lookup:
                return pubchem_lookup

            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError(f'A PubChem integer ({ID[8:]}) identifier was recognized, but it is not in the database.')
    if ID_len > 7:
        if ID_lower[0:7] == 'smiles=':
            smiles_lookup = pubchem_db.search_smiles(ID[7:], autoload)
            if smiles_lookup:
                return smiles_lookup
            else:
                if not autoload:
                    return search_chemical(ID, autoload=True)
                raise ValueError(f'A SMILES identifier ({ID[7:]}) was recognized, but it is not in the database.')

    # Try the smiles lookup by default for all except a few hardcoded cases that the smiles/smarts intersect
    if ID in FORMULA_SEARCH_BEFORE_SMILES_EXCEPTIONS:
        return pubchem_db.search_formula(ID, autoload)
    smiles_lookup = pubchem_db.search_smiles(ID, autoload)
    if smiles_lookup:
        return smiles_lookup

    try:
        serialized_formula = serialize_formula(ID)
        formula_query = pubchem_db.search_formula(serialized_formula, autoload)
        if formula_query and type(formula_query) == ChemicalMetadata:
            try:
                # If we found something after serializing the formula, check it is in fact a formula we were given
                # nested_formula_parser(ID, check=True)
                return formula_query
            except:
                pass
    except:
        pass

    # Try a direct lookup with the name - the fastest
    name_lookup = pubchem_db.search_name(ID, autoload)
    if name_lookup:
        return name_lookup

#     Permutate through various name options
    ID_no_space = ID.replace(' ', '')
    ID_no_space_dash = ID_no_space.replace('-', '')

    for name in [ID, ID_no_space, ID_no_space_dash]:
        for name2 in [name, name.lower()]:
            name_lookup = pubchem_db.search_name(name2, autoload)
            if name_lookup:
                return name_lookup
            # if we have a CAS number with an accidental space
            if check_CAS(name2):
                CAS_lookup = pubchem_db.search_CAS(name2, autoload)
                if CAS_lookup:
                    return CAS_lookup

    if ID[-1] == ')' and '(' in ID:
        # Try to match in the form 'water (H2O)'
        # first_identifier, second_identifier = ID[0:-1].split('(', 1)
        last_open = ID.rindex('(')
        first_identifier = ID[:last_open].strip()
        second_identifier = ID[last_open + 1:].rstrip(')')

        try:
            CAS1 = search_chemical(first_identifier, autoload)
            CAS2 = search_chemical(second_identifier, autoload)
            if CAS1 == CAS2:
                return CAS1
        except:
            pass

    if not autoload:
        return _search_chemical(ID, autoload=True)

    raise ValueError(f'Chemical name ({ID}) not recognized')





### DIPPR Database, chemical list only
# Obtained via the command:
# list(pd.read_excel('http://www.aiche.org/sites/default/files/docs/pages/sponsor_compound_list-2014.xlsx')['Unnamed: 2'])[2:]
# This is consistently faster than creating a list and then making the set.
@mark_numba_incompatible
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

class CommonMixtureMetadata:
    """Class for storing metadata on predefined chemical mixtures.

    Attributes
    ----------
    name : str
        Name of the mixture, [-]
    source : str
        Source of the mixture composition, [-]
    N : int
        Number of chemicals in the mixture, [-]
    CASs : list[str]
        CAS numbers of the mixture, [-]
    ws : list[float]
        Mass fractions of chemicals in the mixture, [-]
    zs : list[float]
        Mole fractions of chemicals in the mixture, [-]
    names : list[str]
        List of names of the chemicals in the mixture, [-]
    synonyms : list[str]
        List of synonyms of the mixture which can also be used to look it up,
        [-]
    """

    __slots__ = ['name', 'CASs', 'N', 'source', 'names', 'ws', 'zs',
                 'synonyms']
    def __repr__(self):
        return (f'<MixtureMetadata, name={self.name}, N={self.N}, CASs={self.CASs}, ws={self.ws}, zs={self.zs}>')

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


@mark_numba_incompatible
def mixture_from_any(ID):
    """Search by string for a mixture in the included common mixture database.
    The database primarily contains refrigerant blends. The variable
    `common_mixtures` contains all loaded entries.

    Parameters
    ----------
    ID : list[str] or str
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
    if not mixture_composition_loaded:  # pragma: no cover
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

@mark_numba_incompatible
def IDs_to_CASs(IDs):
    """Find the CAS numbers for multiple chemicals names at once. Also supports
    having a string input which is a common mixture name in the database.
    An error will be raised if any of the chemicals cannot be found.


    Parameters
    ----------
    IDs : list[str] or str
        A string or 1-element list containing the name which may represent a
        mixture.

    Returns
    -------
    CASs : list[str]
        CAS numbers of found chemicals, [-]

    Notes
    -----
    White space, '-', and upper case letters are removed in the search.

    Examples
    --------
    >>> IDs_to_CASs('R512A')
    ['811-97-2', '75-37-6']
    >>> IDs_to_CASs(['norflurane', '1,1-difluoroethane'])
    ['811-97-2', '75-37-6']
    """
    if hasattr(IDs, 'strip') or (isinstance(IDs, list) and len(IDs) == 1):
        try:
            # Assume the name was a pre-defined mixture
            mixname = mixture_from_any(IDs)
            return mixname.CASs
        except:
            if hasattr(IDs, 'strip'): # It it one chemical?
                return [CAS_from_any(IDs)]
    return [CAS_from_any(ID) for ID in IDs]

cryogenics = {
  '132259-10-0': 'Air',
  '7440-37-1':  'Argon',
  '630-08-0':   'carbon monoxide',
  '7782-39-0':  'deuterium',
  '7782-41-4':  'fluorine',
  '7440-59-7':  'helium',
  '1333-74-0':  'hydrogen',
  '7439-90-9':  'krypton',
  '74-82-8':    'methane',
  '7440-01-9':  'neon',
  '7727-37-9':  'nitrogen',
  '7782-44-7':  'oxygen',
  '7440-63-3':  'xenon'
}

common_commercial_gases = {"7440-37-1": "Argon", "124-38-9": "Carbon Dioxide", "7440-59-7":
      "Helium", "7440-01-9": "Neon", "7727-37-9": "Nitrogen",
      "7440-63-3": "Xenon", "10102-43-9": "Nitric Oxide", "10102-44-0":
      "Nitrogen Dioxide", "7782-44-7": "Oxygen", "132259-10-0": "Air",
      "7439-90-9": "krypton", "10043-92-2": "radon", "7732-18-5":
      "water", "7782-50-5": "chlorine", "7782-41-4": "fluorine"}



_pubchem_db_loaded = False
@mark_numba_incompatible
def get_pubchem_db():
    """Helper function to delay the creation of the pubchem_db object.

    This avoids loading the database when it is not needed.
    """
    global _pubchem_db_loaded, pubchem_db
    if _pubchem_db_loaded:  # pragma: no cover
        return pubchem_db
    else:
        pubchem_db = ChemicalMetadataDB()
        # pubchem_db = ChemicalMetadataDiskDB()
    _pubchem_db_loaded = True
    return pubchem_db

mixture_composition_loaded = False

@mark_numba_incompatible
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
        elif name in ('common_mixtures', 'common_mixtures_by_synonym'):
            load_mixture_composition()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")  # pragma: no cover
else:  # pragma: no cover
    if can_load_data:
        get_pubchem_db()
        load_mixture_composition()
