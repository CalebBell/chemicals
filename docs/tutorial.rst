Chemicals tutorial
==================

Importing
---------

Chemicals can be imported as a standalone library, or all of its functions
and classes may be imported with star imports:

>>> import numpy as np
>>> import chemicals # Good practice
>>> from chemicals import * # Bad practice but convenient

All functions are available from either the main chemicals module or the 
submodule; i.e. both chemicals.Antoine and 
chemicals.vapor_pressure.Antoine are valid ways of accessing a function.

Design philosophy
-----------------
Like all libraries, this was developed to scratch my own itches.

The bulk of this library's API is considered quite stable; enhancements to
functions and classes will still happen, and default methods when using a generic
correlation interface may change to newer and more accurate correlations as
they are published and reviewed.

All functions are designed to accept inputs in base SI units. However, any
set of consistent units given to a function will return a consistent result.
The user is directed to unit conversion libraries such as
`pint <https://github.com/hgrecco/pint>`_ to perform unit conversions if they
prefer not to work in SI units. The tutorial for using it with chemicals is
at :doc:`chemicals.units <chemcials.units>`.

There are two ways to use numpy arrays with chemicals. The easiest way to use numpy is a `vectorized` module,
which wraps all of the chemicals functions with np.vectorize. Instead of importing
from chemicals, the user can import from :doc:`chemicals.vectorized <chemicals.vectorized>`:

>>> from chemicals.vectorized import *
>>> Antoine(np.linspace(100, 200, 5), A=8.95894, B=510.595, C=-15.95)
array([7.65674361e+02, 1.89116754e+04, 1.41237759e+05, 5.60609191e+05,
       1.53010431e+06])

Inputs do not need to be numpy arrays; they can be any iterable:

>>> import chemicals.vectorized
>>> chemicals.vectorized.Tc(['108-88-3', '7732-18-5'])
array([591.75, 647.14])


It is possible to switch back and forth between the namespaces with a subsequent
import:

>>> from chemicals import *

The second way is `Numba <https://github.com/numba/numba>`_. This
optional dependency provides the speed you expect from NumPy arrays -
or better. In some cases, much better. The tutorial for using it
is at :doc:`chemicals.numba <chemicals.numba>`, but in general use it the same way but
with a different import.

>>> import chemicals.numba_vectorized # doctest: +SKIP

Note that numba can also be used to speed up scalar calculations without numpy.

>>> import chemicals.numba # doctest: +SKIP

Working with Elements
---------------------
Chemicals contains a periodic table.

>>> from chemicals import *
>>> periodic_table.Na
<Element Sodium (Na), number 11, MW=22.98977>
>>> periodic_table.U.MW
238.02891
>>> periodic_table['Th'].CAS
'7440-29-1'
>>> periodic_table.lead.protons
82
>>> periodic_table['7440-57-5'].symbol
'Au'
>>> len(periodic_table)
118
>>> 'gold' in periodic_table
True
>>> periodic_table.He.protons, periodic_table.He.neutrons, periodic_table.He.electrons # Standard number of protons, neutrons, electrons
(2, 2, 2)
>>> periodic_table.He.phase # Phase of the element in the standard state
'g'
>>> periodic_table.He.Hf # Heat of formation in standard state in J/mol - by definition 0
0.0
>>> periodic_table.He.S0 # Absolute entropy (J/(mol*K) in standard state - non-zero)
126.2
>>> periodic_table.Kr.block, periodic_table.Kr.period, periodic_table.Kr.group
('p', 4, 18)
>>> periodic_table.Rn.InChI
'Rn'
>>> periodic_table.Rn.smiles
'[Rn]'
>>> periodic_table.Pu.number
94
>>> periodic_table.Pu.PubChem
23940
>>> periodic_table.Bi.InChI_key
'JCXGWMGPZLAOME-UHFFFAOYSA-N'


The periodic table is a singleton of the periodic table class :py:class:`~.PeriodicTable`.
Each attribute accessed is a reference to an element object :py:class:`~.Element`.
The elements are the basic building blocks of every chemical.

Working with Chemical Identifiers
---------------------------------
Chemicals comes with a large library of chemical identifiers.
Chemicals has various ways of searching through its database.
There are a number of different support chemical identifiers as well.

**CAS numbers** - These are the primary identifiers in Chemicals. A CAS number uniquely identifies a chemical molecule. 7732-18-5 is the CAS number for water. Sometimes, it also identifies the phase of the chemical. `7440-44-0 <https://commonchemistry.cas.org/detail?cas_rn=7440-44-0>`_ is the CAS number for carbon in general, but `7782-42-5 <https://commonchemistry.cas.org/detail?cas_rn=7782-42-5>`_  is the CAS number for graphite and `7782-40-3 <https://commonchemistry.cas.org/detail?cas_rn=7782-40-3>`_ is the CAS number for diamond. Note that because these are assigned by people, mistakes are made and often multiple CAS numbers point to the same compound. Common Chemistry lists 57 "retired" CAS numbers which point to the element carbon. The CAS numbers in Chemicals come mostly from PubChem as there was no Common Chemistry project back then.

**PubChem IDs** - These are the identifiers for each compound in the PubChem database. Most of the metadata in Chemicals came from PubChem. `962 <https://pubchem.ncbi.nlm.nih.gov/compound/962>`_ is the Pubchem identifier for water. Each entry in PubChem comes with a structure. Sometimes structures are found to be duplicates of each other and entries are merged; these identifiers are assigned automatically by the NIH.

**Smiles** - These are actual chemicals structures, rendered into easily readable text. Multiple smiles strings can represent the same compound; they are not unique. Both "C(=O)=O" and "O=C=O" are valid SMILES strings for identifying CO2. Programs like `rdkit <https://www.rdkit.org/>`_ can create a computational representation of the molecule from a SMILES string. To solve this duplication issue, a concept of a canonical SMILES string was developed which is supposed to be unique, but in general is not reliable at all and only consistent within the same molecular modeling software. There is in general no organization which controls this format, but a there is an effort in the open source community to standardize the format called `opensmiles <http://opensmiles.org/>`_

**Chemical Formula** - These are what every student is taught in chemistry class. H2O is the formula for water. Is OH2 also a valid formula? Yes. There is a convention called the Hill convention (implemented in chemicals as :py:func:`~.atoms_to_Hill` which specified the H2O is how the formula should be written. Not all formulas, especially inorganic formulas or older formulas, follow this convention. Formulas are in general NOT unique. Even simple formulas which seem like there should only be one compound with that formula are often duplicated; carbonic acid and performic acid both have the formula "CH2O3". Searching Chemical's databases with a formula is a common mistake by users. While you can do it and you may get a match, there is no guarantee the match you wanted was found. The following snippet of code counts the number of compounds with the same formula as asprin; illustrating why searching by formula is a bad idea.

>>> from chemicals.identifiers import pubchem_db
>>> len(list(i for i in pubchem_db if i.formula == 'C9H8O4'))
20

**Chemical name** - Anyone can call a chemical by any name, so predictably names are a mess. A large number of names were retrieved from PubChem, and form the basis for searches by name in Chemicals. Only one chemical hit will be found for each name search. There is an effort by IUPAC to systematically generate names for each chemical structure, called `OPSIN <https://opsin.ch.cam.ac.uk/>`_. Most chemicals in Chemicals have a correct, associated IUPAC name retrieved from PubChem. There are in the range of a million names that can be looked by in Chemicals.

**InChI** - Short for the IUPAC International Chemical Identifier, these are programmatically derived strings which represent a compound. A non-profit was established to maintain a software package to manage this format; it is not like SMILES where lots of software implement the format. There contain all the information required to form a structure. There is a variant which is truly unique per compound; this is what is in Chemicals. They have more features than SMILES strings. "C6H14/c1-3-5-6-4-2/h3-6H2,1-2H3" is a sample string, for n-hexane. This is the best possible type of an identifier for a chemical. These can get to be quite long for complex structures.

**InChI key** - A 27-character hash of the unique InChI identifier. These are also in Chemicals and generated by the same InChI software. These were intended to be unique, and easy to search for as search engines don't search for InChI strings well. Some collisions have been detected. 'VLKZOEOYAKHREP-UHFFFAOYSA-N' is the InChI key for n-hexane as an example.

The main interface for looking up a chemical from one of these identifying markers is :py:func:`~.search_chemical`. The search can be performed with any of the following input forms:

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

Some sample queries illustrating the topic:

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


Each of those queries returns a :obj:`ChemicalMetadata <chemicals.identifiers.ChemicalMetadata>` object. The object holds the chemical metadata. It is an almost unbearable task to assemble a chemical property database. Making a database of chemical metadata is only slightly easier. The chemical metadata database doesn't have any information whatsoever about about any chemical properties; only information about the chemical structure and those identifiers mentioned above. Each of those identifiers is an attribute of the returned object.

>>> water = search_chemical('water')
>>> (water.pubchemid, water.formula, water.smiles, water.InChI, water.InChI_key, water.CASs)
(962, 'H2O', 'O', 'H2O/h1H2', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N', '7732-18-5')
>>> water.common_name, water.iupac_name, len(water.synonyms)
('water', 'oxidane', 89)


