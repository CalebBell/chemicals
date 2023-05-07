"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>
Copyright (C) 2020 Yoel Rene Cortes-Pena <yoelcortes@gmail.com>

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

This module contains lookup functions for the property dipole moment.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Lookup Functions
----------------
.. autofunction:: chemicals.dipole.dipole_moment
.. autofunction:: chemicals.dipole.dipole_moment_methods
.. autodata:: chemicals.dipole.dipole_moment_all_methods

"""
__all__ = ['dipole_moment',
           'dipole_moment_methods',
           'dipole_moment_all_methods']

from chemicals import data_reader as dr
from chemicals.data_reader import (
    data_source,
    database_constant_lookup,
    list_available_methods_from_df_dict,
    register_df_source,
    retrieve_any_from_df_dict,
    retrieve_from_df_dict,
)
from chemicals.miscdata import PSI4_2022A
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

# %% Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')

CCCBDB = 'CCCBDB'
MULLER = 'MULLER'
POLING = 'POLING'

register_df_source(folder, 'Poling Dipole.csv')
register_df_source(folder, 'cccbdb.nist.gov Dipoles.csv')
register_df_source(folder, 'Muller Supporting Info Dipoles.csv')
register_df_source(folder, 'psi4_dipoles.tsv')


_dipole_data_loaded = False
@mark_numba_incompatible
def _load_dipole_data():
    global dipole_data_CCDB, dipole_data_Muller, dipole_data_Poling, dipole_sources
    dipole_data_CCDB = data_source('cccbdb.nist.gov Dipoles.csv')
    dipole_data_Muller = data_source('Muller Supporting Info Dipoles.csv')
    dipole_data_Poling = data_source('Poling Dipole.csv')
    dipole_data_psi4_2022a = data_source('psi4_dipoles.tsv')
    dipole_sources = {
        CCCBDB: dipole_data_CCDB,
        MULLER: dipole_data_Muller,
        POLING: dipole_data_Poling,
        PSI4_2022A: dipole_data_psi4_2022a,
    }

if PY37:
    def __getattr__(name):
        if name in ('dipole_data_Poling', 'dipole_data_CCDB', 'dipole_data_Muller', 'dipole_data_psi4_2022a'):
            _load_dipole_data()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else: # pragma: no cover
    if can_load_data:
        _load_dipole_data()

# %% Dipole moment functions

dipole_moment_all_methods = (CCCBDB, MULLER, POLING, PSI4_2022A)
"""Tuple of method name keys. See the `dipole` for the actual references"""

@mark_numba_incompatible
def dipole_moment_methods(CASRN):
    """Return all methods available to obtain the dipole moment for the desired
    chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the dipole moment with the given
        inputs.

    See Also
    --------
    dipole_moment
    """
    if not _dipole_data_loaded: _load_dipole_data()
    return list_available_methods_from_df_dict(dipole_sources, CASRN, 'dipole_moment')

@mark_numba_incompatible
def dipole_moment(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's dipole moment.
    Lookup is based on CASRNs. Will automatically select a data source to use
    if no method is provided; returns None if the data is not available.

    Preferred source is 'CCCBDB'. Considerable variation in reported data has
    found.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    dipole : float
        Dipole moment, [debye]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'CCCBDB', 'MULLER', or
        'POLING', 'PSI4_2022A'. All valid values are also held in the list
        `dipole_all_methods`.

    Notes
    -----
    A total of three sources are available for this function. They are:

        * 'CCCBDB', a series of critically evaluated data for compounds in
          [1]_, intended for use in predictive modeling.
        * 'MULLER', a collection of data in a
          group-contribution scheme in [2]_.
        * 'POLING', in the appendix in [3].
        * 'PSI4_2022A', values computed using the Psi4 version 1.3.2 quantum
          chemistry software, with initialized positions from rdkit's EmbedMolecule
          method, the basis set 6-31G** and the method mp2 [4]_.

    This function returns dipole moment in units of Debye. This is actually
    a non-SI unit; to convert to SI, multiply by 3.33564095198e-30 and its
    units will be in ampere*second^2 or equivalently and more commonly given,
    coulomb*second. The constant is the result of 1E-21/c, where c is the
    speed of light.

    Examples
    --------
    >>> dipole_moment(CASRN='64-17-5')
    1.44

    See Also
    --------
    dipole_moment_methods

    References
    ----------
    .. [1] NIST Computational Chemistry Comparison and Benchmark Database
       NIST Standard Reference Database Number 101 Release 17b, September 2015,
       Editor: Russell D. Johnson III http://cccbdb.nist.gov/
    .. [2] Muller, Karsten, Liudmila Mokrushina, and Wolfgang Arlt. "Second-
       Order Group Contribution Method for the Determination of the Dipole
       Moment." Journal of Chemical & Engineering Data 57, no. 4 (April 12,
       2012): 1231-36. doi:10.1021/je2013395.
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [4] Turney, Justin M., Andrew C. Simmonett, Robert M. Parrish, Edward G.
       Hohenstein, Francesco A. Evangelista, Justin T. Fermann, Benjamin J.
       Mintz, et al. "Psi4: An Open-Source Ab Initio Electronic Structure
       Program." WIREs Computational Molecular Science 2, no. 4 (2012): 556-65.
       https://doi.org/10.1002/wcms.93.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'dipole_moment')
        if found: return val
    if not _dipole_data_loaded: _load_dipole_data()
    if method:
        return retrieve_from_df_dict(dipole_sources, CASRN, 'dipole_moment', method)
    else:
        return retrieve_any_from_df_dict(dipole_sources, CASRN, 'dipole_moment')
