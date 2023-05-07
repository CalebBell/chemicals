"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2022 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRGNGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRGGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARGSING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This module contains various radius of gyration functions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Lookup Functions
----------------
.. autofunction:: chemicals.molecular_geometry.RG
.. autofunction:: chemicals.molecular_geometry.RG_methods
.. autodata:: chemicals.molecular_geometry.RG_all_methods
.. autofunction:: chemicals.molecular_geometry.linear
.. autofunction:: chemicals.molecular_geometry.linear_methods
.. autodata:: chemicals.molecular_geometry.linear_all_methods

"""

__all__ = ['RG', 'RG_methods', 'RG_all_methods',
           'linear', 'linear_methods', 'linear_all_methods']

from chemicals import data_reader as dr
from chemicals.data_reader import (
    data_source,
    database_constant_lookup,
    list_available_methods_from_df_dict,
    register_df_source,
    retrieve_any_from_df_dict,
    retrieve_from_df_dict,
)
from chemicals.miscdata import CHEMSEP, PSI4_2022A
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

# Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')
register_df_source(folder, 'psi4_radius_of_gyrations.tsv')
register_df_source(folder, 'chemsep_radius_of_gyrations.tsv')

register_df_source(folder, 'psi4_linear.tsv')

_RG_data_loaded = False
@mark_numba_incompatible
def _load_RG_data():
    global _RG_data_loaded, radius_of_gyration_data_psi4_2022a, linear_data_psi4_2022a
    global RG_sources, linear_sources, radius_of_gyration_data_chemsep
    radius_of_gyration_data_psi4_2022a = data_source('psi4_radius_of_gyrations.tsv')
    radius_of_gyration_data_chemsep = data_source('chemsep_radius_of_gyrations.tsv')
    linear_data_psi4_2022a = data_source('psi4_linear.tsv')
    RG_sources = {
        PSI4_2022A: radius_of_gyration_data_psi4_2022a,
        CHEMSEP: radius_of_gyration_data_chemsep,
    }

    linear_sources = {
        PSI4_2022A: linear_data_psi4_2022a,
    }

if PY37:
    def __getattr__(name):
        if name in ('radius_of_gyration_data_psi4_2022a',
                    'radius_of_gyration_data_chemsep',
                    'linear_data_psi4_2022a',
                    'RG_sources',
                    'linear_sources'):
            _load_RG_data()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:
    if can_load_data:
        _load_RG_data()

#  Refractive index functions

RG_all_methods = (PSI4_2022A, CHEMSEP)
"""Tuple of method name keys. See the `RG` for the actual references"""

@mark_numba_incompatible
def RG_methods(CASRN):
    """Return all methods available to obtain the radius of gyration for the
    desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the RG with the given
        inputs.

    See Also
    --------
    RG
    """
    if not _RG_data_loaded: _load_RG_data()
    return list_available_methods_from_df_dict(RG_sources, CASRN, 'RG')

@mark_numba_incompatible
def RG(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's radius of gyration.
    Lookup is based on CASRNs. Will automatically select a data source
    to use if no method is provided; returns None if the data is not available.

    Function has data for approximately 670 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    RG : float
        Radius of gyration, [m]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        RG_methods

    Notes
    -----
    The available sources are as follows:

        * 'PSI4_2022A', values computed using the Psi4 version 1.3.2 quantum
          chemistry software, with initialized positions from rdkit's EmbedMolecule
          method, the basis set 6-31G** and the method mp2 [1]_.
        * 'CHEMSEP', from the databank included and distributed with the licence
          notice ChemSep v8.1 pure component data - Copyright (c) Harry Kooijman
          and Ross Taylor (2018) - http://www.perlfoundation.org/artistic_license_2_0.
          A small portion of the data is used.

    Examples
    --------
    >>> RG(CASRN='64-17-5')
    2.225e-10

    References
    ----------
    .. [1] Turney, Justin M., Andrew C. Simmonett, Robert M. Parrish, Edward G.
       Hohenstein, Francesco A. Evangelista, Justin T. Fermann, Benjamin J.
       Mintz, et al. "Psi4: An Open-Source Ab Initio Electronic Structure
       Program." WIREs Computational Molecular Science 2, no. 4 (2012): 556-65.
       https://doi.org/10.1002/wcms.93.
    .. [2] Kooijman, Harry A., and Ross Taylor. The ChemSep Book. Books on
       Demand Norderstedt, Germany, 2000.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'RG')
        if found: return val
    if not _RG_data_loaded: _load_RG_data()
    if method:
        value = retrieve_from_df_dict(RG_sources, CASRN, 'RG', method)
    else:
        value = retrieve_any_from_df_dict(RG_sources, CASRN, 'RG')
    return value



linear_all_methods = (PSI4_2022A,)
"""Tuple of method name keys. See the `linear` for the actual references"""

@mark_numba_incompatible
def linear_methods(CASRN):
    """Return all methods available to obtain whether or not the
    desired chemical is linear.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the linear with the given
        inputs.

    See Also
    --------
    linear
    """
    if not _RG_data_loaded: _load_RG_data()
    return list_available_methods_from_df_dict(linear_sources, CASRN, 'linear')

@mark_numba_incompatible
def linear(CASRN, method=None):
    r'''This function handles the retrieval whether or not a chemical is linear.
    Lookup is based on CASRNs. Will automatically select a data source
    to use if no method is provided; returns None if the data is not available.

    Function has data for approximately 300 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    linear : bool
        Whether or not the chemical is linear, [-]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        linear_methods

    Notes
    -----
    The available sources are as follows:

        * 'PSI4_2022A', values computed using the Psi4 version 1.3.2 quantum
          chemistry software, with initialized positions from rdkit's EmbedMolecule
          method, the basis set 6-31G** and the method mp2 [1]_.

    .. warning::
        This function does not yet have a reliable data source.

    Examples
    --------
    >>> linear(CASRN='64-17-5')
    False

    References
    ----------
    .. [1] Turney, Justin M., Andrew C. Simmonett, Robert M. Parrish, Edward G.
       Hohenstein, Francesco A. Evangelista, Justin T. Fermann, Benjamin J.
       Mintz, et al. "Psi4: An Open-Source Ab Initio Electronic Structure
       Program." WIREs Computational Molecular Science 2, no. 4 (2012): 556-65.
       https://doi.org/10.1002/wcms.93.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'linear')
        if found: return bool(val)
    if not _RG_data_loaded: _load_RG_data()
    if method:
        value = retrieve_from_df_dict(linear_sources, CASRN, 'linear', method)
    else:
        value = retrieve_any_from_df_dict(linear_sources, CASRN, 'linear')
    return bool(value)
