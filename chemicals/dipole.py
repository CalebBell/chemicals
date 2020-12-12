# -*- coding: utf-8 -*-
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

from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

# %% Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')

CCCBDB = 'CCCBDB'
MULLER = 'MULLER'
POLING = 'POLING'

register_df_source(folder, 'Poling Dipole.csv')
register_df_source(folder, 'cccbdb.nist.gov Dipoles.csv')
register_df_source(folder, 'Muller Supporting Info Dipoles.csv')


_dipole_data_loaded = False
def _load_dipole_data():
    global dipole_data_CCDB, dipole_data_Muller, dipole_data_Poling, dipole_sources
    dipole_data_CCDB = data_source('cccbdb.nist.gov Dipoles.csv')
    dipole_data_Muller = data_source('Muller Supporting Info Dipoles.csv')
    dipole_data_Poling = data_source('Poling Dipole.csv')
    dipole_sources = {
        CCCBDB: dipole_data_CCDB,
        MULLER: dipole_data_Muller,
        POLING: dipole_data_Poling,
    }

if PY37:
    def __getattr__(name):
        if name in ('dipole_data_Poling', 'dipole_data_CCDB', 'dipole_data_Muller'):
            _load_dipole_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else: # pragma: no cover
    if can_load_data:
        _load_dipole_data()

# %% Dipole moment functions

dipole_moment_all_methods = (CCCBDB, MULLER, POLING)
'''Tuple of method name keys. See the `dipole` for the actual references'''

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
    return list_available_methods_from_df_dict(dipole_sources, CASRN, 'Dipole')

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
        'POLING'. All valid values are also held in the list `dipole_all_methods`.

    Notes
    -----
    A total of three sources are available for this function. They are:

        * 'CCCBDB', a series of critically evaluated data for compounds in
          [1]_, intended for use in predictive modeling.
        * 'MULLER', a collection of data in a
          group-contribution scheme in [2]_.
        * 'POLING', in the appendix in [3].

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
    '''
    if not _dipole_data_loaded: _load_dipole_data()
    if method:
        return retrieve_from_df_dict(dipole_sources, CASRN, 'Dipole', method)
    else:
        return retrieve_any_from_df_dict(dipole_sources, CASRN, 'Dipole')
