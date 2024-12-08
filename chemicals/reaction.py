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

This module contains lookup functions enthalpies and standard entropies of
formation. Lookup functions are availa for the liquid, solid, and gas states.
A compound may be in more than one lookup function.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Solid Heat of Formation
-----------------------
.. autofunction:: chemicals.reaction.Hfs
.. autofunction:: chemicals.reaction.Hfs_methods
.. autodata:: chemicals.reaction.Hfs_all_methods

Liquid Heat of Formation
------------------------
.. autofunction:: chemicals.reaction.Hfl
.. autofunction:: chemicals.reaction.Hfl_methods
.. autodata:: chemicals.reaction.Hfl_all_methods

Gas Heat of Formation
---------------------
.. autofunction:: chemicals.reaction.Hfg
.. autofunction:: chemicals.reaction.Hfg_methods
.. autodata:: chemicals.reaction.Hfg_all_methods

Solid Absolute Entropy
----------------------
.. autofunction:: chemicals.reaction.S0s
.. autofunction:: chemicals.reaction.S0s_methods
.. autodata:: chemicals.reaction.S0s_all_methods

Liquid Absolute Entropy
-----------------------
.. autofunction:: chemicals.reaction.S0l
.. autofunction:: chemicals.reaction.S0l_methods
.. autodata:: chemicals.reaction.S0l_all_methods

Gas Absolute Entropy
--------------------
.. autofunction:: chemicals.reaction.S0g
.. autofunction:: chemicals.reaction.S0g_methods
.. autodata:: chemicals.reaction.S0g_all_methods

Utility Functions
-----------------
.. autofunction:: chemicals.reaction.Gibbs_formation
.. autofunction:: chemicals.reaction.entropy_formation
.. autofunction:: chemicals.reaction.Hf_basis_converter

Chemical Reactions
------------------
.. autofunction:: chemicals.reaction.balance_stoichiometry
.. autofunction:: chemicals.reaction.stoichiometric_matrix
.. autofunction:: chemicals.reaction.stoichiometry_mass_to_molar
.. autofunction:: chemicals.reaction.stoichiometry_molar_to_mass
.. autofunction:: chemicals.reaction.stoichiometry_MW_error
.. autofunction:: chemicals.reaction.standard_formation_reaction
"""

__all__ = ['Hfg', 'Hfl', 'Hfs', 'S0g', 'S0l', 'S0s',
           'Hfl_methods', 'Hfg_methods', 'Hfs_methods',
           'S0l_methods', 'S0g_methods', 'S0s_methods',
           'Hfl_all_methods', 'Hfg_all_methods', 'Hfs_all_methods',
           'S0l_all_methods', 'S0g_all_methods', 'S0s_all_methods',
           'Gibbs_formation', 'entropy_formation', 'Hf_basis_converter',
           'balance_stoichiometry', 'stoichiometric_matrix',
           'stoichiometry_molar_to_mass', 'stoichiometry_mass_to_molar',
           'standard_formation_reaction', 'stoichiometry_MW_error']

from math import ceil, log10, floor

from chemicals import data_reader as dr
from chemicals import heat_capacity, miscdata
from chemicals.data_reader import (
    data_source,
    database_constant_lookup,
    list_available_methods_from_df_dict,
    register_df_source,
    retrieve_any_from_df_dict,
    retrieve_from_df_dict,
)
from chemicals.elements import periodic_table, simple_formula_parser
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

# %% Register data sources and lazy load them
CRC = 'CRC'
YAWS = 'YAWS'
API_TDB_G = 'API_TDB_G'
ATCT_L = 'ATCT_L'
ATCT_G = 'ATCT_G'
TRC = 'TRC'

folder = os_path_join(source_path, 'Reactions')
register_df_source(folder, 'API TDB Albahri Hf (g).tsv')
register_df_source(folder, 'ATcT 1.112 (g).tsv')
register_df_source(folder, 'ATcT 1.112 (l).tsv')
register_df_source(folder, 'Yaws Hf S0 (g).tsv')
register_df_source(folder, 'JANAF_1998.tsv')
_reaction_data_loaded = False
def _load_reaction_data():
    global Hfg_API_TDB_data, Hfg_ATcT_data, Hfl_ATcT_data, Hfg_S0g_YAWS_data
    global Hfg_sources, Hfl_sources, Hfs_sources
    global S0g_sources, S0l_sources, S0s_sources
    global _reaction_data_loaded
    Hfg_API_TDB_data = data_source('API TDB Albahri Hf (g).tsv')
    Hfg_ATcT_data = data_source('ATcT 1.112 (g).tsv')
    Hfl_ATcT_data = data_source('ATcT 1.112 (l).tsv')
    Hfg_S0g_YAWS_data = data_source('Yaws Hf S0 (g).tsv')
    JANAF_1998_data = data_source('JANAF_1998.tsv')
    _reaction_data_loaded = True
    S0g_sources = {
        CRC: heat_capacity.CRC_standard_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
        miscdata.JANAF: JANAF_1998_data,
        YAWS: Hfg_S0g_YAWS_data,
    }
    S0l_sources = {
        CRC: heat_capacity.CRC_standard_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
        miscdata.JANAF: JANAF_1998_data,
    }
    S0s_sources = {
        CRC: heat_capacity.CRC_standard_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
    }
    Hfg_sources = {
        ATCT_G: Hfg_ATcT_data,
        CRC: heat_capacity.CRC_standard_data,
        API_TDB_G: Hfg_API_TDB_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
        TRC: heat_capacity.TRC_gas_data,
        miscdata.JANAF: JANAF_1998_data,
        YAWS: Hfg_S0g_YAWS_data,
        miscdata.JOBACK: miscdata.joback_predictions,
    }
    Hfl_sources = {
        ATCT_L: Hfl_ATcT_data,
        CRC: heat_capacity.CRC_standard_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
        miscdata.JANAF: JANAF_1998_data,
    }
    Hfs_sources = {
        CRC: heat_capacity.CRC_standard_data,
        miscdata.WEBBOOK: miscdata.webbook_data,
    }

if PY37:
    def __getattr__(name):
        if name in ('Hfg_API_TDB_data', 'Hfg_ATcT_data',
                    'Hfl_ATcT_data', 'Hfg_S0g_YAWS_data', 'JANAF_1998_data',
                    'Hfg_sources', 'Hfl_sources', 'Hfs_sources',
                    'S0g_sources', 'S0l_sources', 'S0s_sources'):
            _load_reaction_data()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:
    if can_load_data:
        _load_reaction_data()


# %% Lookup functions

# TODO: more data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3692305/
# has dippr standard heats of formation, about 55% of the database

Hfs_all_methods = (CRC, miscdata.WEBBOOK)
"""Tuple of method name keys. See the `Hfs` for the actual references"""

@mark_numba_incompatible
def Hfs_methods(CASRN):
    """Return all methods available to obtain the solid-phase heat of
    formation for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Hfs with the given
        inputs.

    See Also
    --------
    Hfs
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(Hfs_sources, CASRN, 'Hfs')

@mark_numba_incompatible
def Hfs(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's solid/crystaline
    standard phase heat of formation. The lookup is based on CASRNs. Will
    automatically select a data source to use if no method is provided; returns
    None if the data is not available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Hfs : float
        Solid standard-state heat of formation, [J/mol]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        Hfs_methods

    Notes
    -----
    Sources are:

        * 'CRC', from the CRC handbook (1360 values) [1]_
        * 'WEBBOOK' (2000 values) [2]_

    Examples
    --------
    >>> Hfs('101-81-5') # Diphenylmethane
    71500.0

    See Also
    --------
    Hfs_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [2] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'Hfs')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfs_sources, CASRN, 'Hfs', method)
    else:
        return retrieve_any_from_df_dict(Hfs_sources, CASRN, 'Hfs')

Hfl_all_methods = (ATCT_L, CRC, miscdata.WEBBOOK, miscdata.JANAF)
"""Tuple of method name keys. See the `Hfl` for the actual references"""

@mark_numba_incompatible
def Hfl_methods(CASRN):
    """Return all methods available to obtain the standard liquid-state heat
    of formation for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Hfl with the given
        inputs.

    See Also
    --------
    Hfl
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(Hfl_sources, CASRN, 'Hfl')

@mark_numba_incompatible
def Hfl(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's liquid standard
    phase heat of formation. The lookup is based on CASRNs. Will automatically
    select a data source to use if no method is provided; returns None if
    the data is not available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Hfl : float
        Liquid standard-state heat of formation, [J/mol]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `Hfl_all_methods`.

    Notes
    -----
    Sources are:

        * 'ATCT_L', the Active Thermochemical Tables version 1.112. [1]_
        * 'CRC', from the CRC handbook (1360 values) [2]_
        * 'WEBBOOK' (2000 values) [3]_

    Examples
    --------
    >>> Hfl('67-56-1')
    -238400.0

    See Also
    --------
    Hfl_methods

    References
    ----------
    .. [1] Ruscic, Branko, Reinhardt E. Pinzon, Gregor von Laszewski, Deepti
       Kodeboyina, Alexander Burcat, David Leahy, David Montoy, and Albert F.
       Wagner. "Active Thermochemical Tables: Thermochemistry for the 21st
       Century." Journal of Physics: Conference Series 16, no. 1
       (January 1, 2005): 561. doi:10.1088/1742-6596/16/1/078.
    .. [2] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [3] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'Hfl')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfl_sources, CASRN, 'Hfl', method)
    else:
        return retrieve_any_from_df_dict(Hfl_sources, CASRN, 'Hfl')

Hfg_all_methods = (ATCT_G, TRC, CRC, miscdata.WEBBOOK, miscdata.JANAF, YAWS, miscdata.JOBACK)
"""Tuple of method name keys. See the `Hfg` for the actual references"""

@mark_numba_incompatible
def Hfg_methods(CASRN):
    """Return all methods available to obtain the gas phase heat of formation
    for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Hfg with the given
        inputs.

    See Also
    --------
    Hfg
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(Hfg_sources, CASRN, 'Hfg')

@mark_numba_incompatible
def Hfg(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's gas heat of
    formation. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Hfg : float
        Ideal gas phase heat of formation, [J/mol]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        Hfg_methods

    Notes
    -----
    Function has data for approximately 8700 chemicals. Sources are:

        * 'ATCT_G', the Active Thermochemical Tables version 1.112 (600 values) [1]_
        * 'TRC', from a 1994 compilation (1750 values) [2]_
        * 'CRC', from the CRC handbook (1360 values) [3]_
        * 'WEBBOOK', a NIST resource [6]_ containing mostly experimental
          and averaged values
        * 'JANAF', the 1998 JANAF values online
        * 'JOBACK', an estimation method for organic substances in [5]_
        * 'YAWS', a large compillation of values, mostly estimated (5000 values) [4]_

    'TRC' data may have come from computational procedures, for example petane
    is off by 30%.

    Examples
    --------
    >>> Hfg('67-56-1')
    -200700.0
    >>> Hfg('67-56-1', method='YAWS')
    -200900.0
    >>> Hfg('67-56-1', method='CRC')
    -201000.0
    >>> Hfg('67-56-1', method='TRC')
    -190100.0

    See Also
    --------
    Hfg_methods

    References
    ----------
    .. [1] Ruscic, Branko, Reinhardt E. Pinzon, Gregor von Laszewski, Deepti
       Kodeboyina, Alexander Burcat, David Leahy, David Montoy, and Albert F.
       Wagner. "Active Thermochemical Tables: Thermochemistry for the 21st
       Century." Journal of Physics: Conference Series 16, no. 1
       (January 1, 2005): 561. doi:10.1088/1742-6596/16/1/078.
    .. [2] Frenkel`, M. L, Texas Engineering Experiment Station, and
       Thermodynamics Research Center. Thermodynamics of Organic Compounds in
       the Gas State. College Station, Tex.: Thermodynamics Research Center,
       1994.
    .. [3] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [4] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    .. [5] Joback, K.G., and R.C. Reid. "Estimation of Pure-Component
       Properties from Group-Contributions." Chemical Engineering
       Communications 57, no. 1-6 (July 1, 1987): 233-43.
       doi:10.1080/00986448708960487.
    .. [6] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'Hfg')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfg_sources, CASRN, 'Hfg', method)
    else:
        return retrieve_any_from_df_dict(Hfg_sources, CASRN, 'Hfg')

S0s_all_methods = (CRC, miscdata.WEBBOOK)
"""Tuple of method name keys. See the `S0s` for the actual references"""

@mark_numba_incompatible
def S0s_methods(CASRN):
    """Return all methods available to obtain the absolute entropy of the
    compound in the solid phase for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the S0s with the given
        inputs.

    See Also
    --------
    S0s
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(S0s_sources, CASRN, 'S0s')

@mark_numba_incompatible
def S0s(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's absolute
    entropy at a reference temperature of 298.15 K and pressure of 1 bar,
    in the solid state. Lookup is based on CASRNs. Will automatically select a
    data source to use if no method is provided; returns None if the data is not
    available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    S0s : float
        Ideal gas standard absolute entropy of compound, [J/mol/K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        `S0s_all_methods`.

    Notes
    -----
    Sources are:

        * 'CRC' [1]_ from the CRC handbook (1360 values)
        * 'WEBBOOK', a NIST resource [2]_ containing mostly experimental
          and averaged values

    Examples
    --------
    >>> S0s('7439-93-2') # Lithium
    29.1

    See Also
    --------
    S0s_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [2] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'S0s')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(S0s_sources, CASRN, 'S0s', method)
    else:
        return retrieve_any_from_df_dict(S0s_sources, CASRN, 'S0s')

S0l_all_methods = (CRC, miscdata.WEBBOOK, miscdata.JANAF)
"""Tuple of method name keys. See the `S0l` for the actual references"""

@mark_numba_incompatible
def S0l_methods(CASRN):
    """Return all methods available to obtain the absolute entropy for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the S0l with the given
        inputs.

    See Also
    --------
    S0l
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(S0l_sources, CASRN, 'S0l')

@mark_numba_incompatible
def S0l(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's absolute
    entropy at a reference temperature of 298.15 K and pressure of 1 bar,
    in the liquid state.

    Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    S0l : float
        Ideal gas standard absolute entropy of compound, [J/mol/K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `S0l_all_methods`.

    Notes
    -----
    Sources are:

        * 'CRC', from the CRC handbook

    Examples
    --------
    >>> S0l('7439-97-6') # Mercury
    75.9

    See Also
    --------
    S0l_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'S0l')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(S0l_sources, CASRN, 'S0l', method)
    else:
        return retrieve_any_from_df_dict(S0l_sources, CASRN, 'S0l')

S0g_all_methods = (CRC, miscdata.WEBBOOK, miscdata.JANAF, YAWS)
"""Tuple of method name keys. See the `S0g` for the actual references"""

@mark_numba_incompatible
def S0g_methods(CASRN):
    """Return all methods available to obtain the S0g for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the S0g with the given
        inputs.

    See Also
    --------
    S0g
    """
    if not _reaction_data_loaded: _load_reaction_data()
    return list_available_methods_from_df_dict(S0g_sources, CASRN, 'S0g')

@mark_numba_incompatible
def S0g(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's absolute
    entropy at a reference temperature of 298.15 K and pressure of 1 bar,
    in the ideal gas state.

    Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    S0g : float
        Ideal gas standard absolute entropy of compound, [J/mol/K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `S0g_all_methods`

    Notes
    -----
    Function has data for approximately 5400 chemicals. Sources are:

        * 'CRC', from the CRC handbook (520 values)
        * 'YAWS', a large compillation of values, mostly estimated (4890 values)
        * 'WEBBOOK', a NIST resource [3]_ containing mostly experimental
          and averaged values

    Examples
    --------
    >>> S0g('67-56-1')
    239.9
    >>> S0g('67-56-1', method='YAWS')
    239.88

    See Also
    --------
    S0g_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [2] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    .. [3] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'S0g')
        if found: return val
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(S0g_sources, CASRN, 'S0g', method)
    else:
        return retrieve_any_from_df_dict(S0g_sources, CASRN, 'S0g')


# %% Converter functions

def Hf_basis_converter(Hvapm, Hf_liq=None, Hf_gas=None):
    r'''This function converts a liquid or gas enthalpy of formation to the
    other. This is useful, as thermodynamic packages often work with ideal-
    gas as the reference state and require ideal-gas enthalpies of formation.

    Parameters
    ----------
    Hvapm : float
        Molar enthalpy of vaporization of compound at 298.15 K or (unlikely)
        the reference temperature, [J/mol]
    Hf_liq : float, optional
        Enthalpy of formation of the compound in its liquid state, [J/mol]
    Hf_gas : float, optional
        Enthalpy of formation of the compound in its ideal-gas state, [J/mol]

    Returns
    -------
    Hf_calc : float, optional
        Enthalpy of formation of the compound in the other state to the one
        provided, [J/mol]

    Notes
    -----

    Examples
    --------
    Calculate the ideal-gas enthalpy of formation for water, from its standard-
    state (liquid) value:

    >>> Hf_basis_converter(44018, Hf_liq=-285830)
    -241812

    Calculate the standard-state (liquid) enthalpy of formation for water, from
    its ideal-gas value:

    >>> Hf_basis_converter(44018, Hf_gas=-241812)
    -285830
    '''
    if Hf_liq is None and Hf_gas is None:
        raise ValueError("Provide either a liquid or a gas enthalpy of formation")
    if Hvapm is None or Hvapm < 0.0:
        raise ValueError("Enthalpy of formation unknown or zero")
    if Hf_liq is None:
        return Hf_gas - Hvapm
    else:
        return Hf_liq + Hvapm

def Gibbs_formation(dHf, S0_abs, dHfs_std, S0_abs_elements, coeffs_elements,
                    T_ref=298.15):
    r'''This function calculates the Gibbs free energy of formation of a
    compound, from its constituent elements.

    The calculated value will be for a "standard-state" value if `dHf` and
    `S0_abs` are provided in the standard state; or it will be in an
    "ideal gas" basis if they are both for an ideal gas. For compounds which
    are gases at STP, the two values are the same.

    Parameters
    ----------
    dHf : float
        Molar enthalpy of formation of the created compound, [J/mol]
    S0_abs : float
        Absolute molar entropy of the created compound at the reference
        temperature, [J/mol/K]
    dHfs_std : list[float]
        List of standard molar enthalpies of formation of all elements used in
        the formation of the created compound, [J/mol]
    S0_abs_elements : list[float]
        List of standard absolute molar entropies at the reference temperature
        of all elements used in the formation of the created compound,
        [J/mol/K]
    coeffs_elements : list[float]
        List of coefficients for each compound (i.e. 1 for C, 2 for H2 if the
        target is methane), in the same order as `dHfs_std` and
        `S0_abs_elements`, [-]
    T_ref : float, optional
        The standard state temperature, default 298.15 K; few values are
        tabulated at other temperatures, [-]

    Returns
    -------
    dGf : float
        Gibbs free energy of formation for the created compound, [J/mol]

    Notes
    -----
    Be careful for elements like Bromine - is the tabulated value for Br2 or
    Br?

    Examples
    --------
    Calculate the standard-state Gibbs free energy of formation for water,
    using water's standard state heat of formation and absolute entropy
    at 298.15 K:

    >>> Gibbs_formation(-285830, 69.91,  [0, 0], [130.571, 205.147], [1, .5])
    -237161.633825

    Calculate the ideal-gas state Gibbs free energy of formation for water,
    using water's ideal-gas state heat of formation and absolute entropy
    at 298.15 K as a gas:

    >>> Gibbs_formation(-241818, 188.825,  [0, 0], [130.571, 205.147], [1, .5])
    -228604.141075

    Calculate the Gibbs free energy of formation for CBrF3 (it is a gas at STP,
    so its standard-state and ideal-gas state values are the same) at 298.15 K:

    >>> Gibbs_formation(-648980, 297.713, [0, 0, 0], [5.74, 152.206, 202.789], [1, .5, 1.5])
    -622649.329975

    Note in the above calculation that the Bromine's `S0` and `Hf` are for Br2;
    and that the value for Bromine as a liquid, which is its standard state,
    is used.

    References
    ----------
    .. [1] "Standard Gibbs Free Energy of Formation Calculations Chemistry
       Tutorial." Accessed March, 2019. https://www.ausetute.com.au/gibbsform.html.
    '''
    N = len(coeffs_elements)
    dH = dHf
    dS = S0_abs
    for i in range(N):
        dH -= dHfs_std[i]*coeffs_elements[i]
        dS -= S0_abs_elements[i]*coeffs_elements[i]
    return dH - T_ref*dS

def entropy_formation(Hf, Gf, T_ref=298.15):
    r'''This function calculates the entropy of formation of a
    compound, from its constituent elements.

    The calculated value will be for a "standard-state" value if `Hf` and
    `Gf` are provided in the standard state; or it will be in an
    "ideal gas" basis if they are both for an ideal gas. For compounds which
    are gases at STP, the two values are the same.

    Parameters
    ----------
    Hf : float
        Molar enthalpy of formation of the compound, [J/mol]
    Gf : float
        Molar Gibbs free energy of formation of the compound, [J/mol]
    T_ref : float, optional
        The standard state temperature, default 298.15 K; few values are
        tabulated at other temperatures, [-]

    Returns
    -------
    S0 : float
        Entropy of formation of the compound, [J/mol/K]

    Notes
    -----

    Examples
    --------
    Entropy of formation of methane:

    >>> entropy_formation(Hf=-74520, Gf=-50490)
    -80.59701492537314

    Entropy of formation of water in ideal gas state:

    >>> entropy_formation(Hf=-241818, Gf=-228572)
    -44.427301693778304
    '''
    return (Hf - Gf)/T_ref


# %% Stoichiometry functions

@mark_numba_incompatible
def stoichiometric_matrix(atomss, reactants):
    r'''This function calculates a stoichiometric matrix of reactants and
    stoichiometric matrix, as required by a solver to compute the reation
    coefficients.

    Parameters
    ----------
    atomss : list[dict[(str, float)]]
        A list of dictionaties of (element, element_count) pairs for each
        chemical, [-]
    reactants : list[bool]
        List of booleans indicating whether each chemical is a reactant (True)
        or a product (False), [-]

    Returns
    -------
    matrix : list[list[float]]
        Chemical reaction matrix for further processing; rows contain element
         counts of each compound, and the columns represent each chemical, [-]

    Notes
    -----
    The rows of the matrix contain the element counts of each compound,
    and the columns represent each chemical.

    Examples
    --------
    MgO2 -> Mg + 1/2 O2
    (k=1)

    >>> stoichiometric_matrix([{'Mg': 1, 'O': 1}, {'Mg': 1}, {'O': 2}], [True, False, False])
    [[1, -1, 0], [1, 0, -2]]


    Cl2 + propylene -> allyl chloride + HCl

    >>> stoichiometric_matrix([{'Cl': 2}, {'C': 3, 'H': 6}, {'C': 3, 'Cl': 1, 'H': 5}, {'Cl': 1, 'H': 1}], [True, True, False, False, False])
    [[0, 3, -3, 0], [2, 0, -1, -1], [0, 6, -5, -1]]


    Al + 4HNO3 -> Al(NO3)3 + NO + 2H2O
    (k=1)

    >>> stoichiometric_matrix([{'Al': 1}, {'H': 1, 'N': 1, 'O': 3}, {'Al': 1, 'N': 3, 'O': 9}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False, False])
    [[1, 0, -1, 0, 0], [0, 1, 0, 0, -2], [0, 1, -3, -1, 0], [0, 3, -9, -1, -1]]

    4Fe + 3O2 -> 2(Fe2O3)
    (k=2)

    >>> stoichiometric_matrix([{'Fe': 1}, {'O': 2}, {'Fe':2, 'O': 3}], [True, True, False])
    [[1, 0, -2], [0, 2, -3]]


    4NH3 + 5O2 -> 4NO + 6(H2O)
    (k=4)

    >>> stoichiometric_matrix([{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False])
    [[3, 0, 0, -2], [1, 0, -1, 0], [0, 2, -1, -1]]


    No unique solution:
    C2H5NO2 + C3H7NO3 + 2C6H14N4O2 + 3C5H9NO2 + 2C9H11NO2 -> 8H2O + C50H73N15O11

    >>> stoichiometric_matrix([{'C': 2, 'H': 5, 'N': 1, 'O': 2}, {'C': 3, 'H': 7, 'N': 1, 'O': 3}, {'C': 6, 'H': 14, 'N': 4, 'O': 2}, {'C': 5, 'H': 9, 'N': 1, 'O': 2}, {'C': 9, 'H': 11, 'N': 1, 'O': 2}, {'H': 2, 'O': 1}, {'C': 50, 'H': 73, 'N': 15, 'O': 11}], [True, True, True, True, True, False, False])
    [[2, 3, 6, 5, 9, 0, -50], [5, 7, 14, 9, 11, -2, -73], [1, 1, 4, 1, 1, 0, -15], [2, 3, 2, 2, 2, -1, -11]]

    References
    ----------
    .. [1] Sen, S. K., Hans Agarwal, and Sagar Sen. "Chemical Equation
       Balancing: An Integer Programming Approach." Mathematical and Computer
       Modelling 44, no. 7 (October 1, 2006): 678-91.
       https://doi.org/10.1016/j.mcm.2006.02.004.
    .. [2] URAVNOTE, NOVOODKRITI PARADOKSI V. TEORIJI, and ENJA KEMIJSKIH
       REAKCIJ. "New Discovered Paradoxes in Theory of Balancing Chemical
       Reactions." Materiali in Tehnologije 45, no. 6 (2011): 503-22.
    '''
    n_compounds = len(atomss)
    elements = set()
    for atoms in atomss:
        elements.update(atoms.keys())
    elements = list(elements)
    elements.sort() # Ensure reproducibility
    n_elements = len(elements)


    matrix = [[0]*n_compounds for _ in range(n_elements)]
    element_to_row = {ele: matrix[idx] for idx, ele in enumerate(elements)}
    for i, atoms in enumerate(atomss):
        if reactants[i]:
            for k, v in atoms.items():
                element_to_row[k][i] = v
        else:
            for k, v in atoms.items():
                element_to_row[k][i] = -v
    return matrix

def round_to_significant(x, significant_digits):
    if x == 0:
        return 0.0
    magnitude = floor(log10(abs(x)))
    scale = 10 ** (significant_digits - 1 - magnitude)
    return round(x * scale) / scale

def check_reaction_balance(matrix, coeffs, atol=1e-13):
    """Check that coefficients satisfy the stoichiometric matrix equation within tolerance."""
    result = [sum(coeff * row[i] for i, coeff in enumerate(coeffs)) 
              for row in matrix]
    return all(abs(x) <= atol for x in result)

def floats_to_ints(float_list, matrix, max_denominator=1000):
    """
    Convert a list of floats to integers until we find a solution that balances.
    All floats are one or larger. The chemical equation is assumed to be reasonable.
    The SVD has already solved the problem, but there is a little numerical noise
    we need to clean up.
    
    Parameters:
    - float_list: List of floats to convert
    - matrix: Stoichiometric matrix to verify balance
    - max_denominator: Maximum scaling factor to consider
    
    Returns:
    - A list of integers scaled from the original floats
    """
    for D in range(1, max_denominator + 1):
        # Calculate rounded integers
        # in practice this works extremely well and rarely goes above 10
        # it is extremely fast compared to a Fraction/Decimal approach
        rounded = [int(round(D * x)) for x in float_list]
        # Check if these coefficients actually balance the reaction
        if check_reaction_balance(matrix, rounded):
            return rounded
    return float_list  # If we still can't find a solution, return original floats


def balance_stoichiometry(matrix, rounding=9, allow_fractional=False):
    r'''This function balances a chemical reaction.

    Parameters
    ----------
    matrix : list[list[float]]
        Chemical reaction matrix for further processing; rows contain element
         counts of each compound, and the columns represent each chemical, [-]
    rounding : int
        Roughly the number of digits of rounding to apply to the answer. As matrix
        routines are used, there is some noise; if this number is too high, the
        coefficients may become very large numberes, which are still in a correct
        ratio to each other, but are extremely ugly, [-]
    allow_fractional : bool
        Whether or not to allow the answers to be fractions, or to force them to
        integers. Setting this to True speeds up the calculation, and allows
        setting rounding arbitrarily high, [-]

    Returns
    -------
    coefficients : list[float]
        Balanced coefficients; all numbers are positive, [-]

    Notes
    -----
    Balance the reaction 4 NH3 + 5 O2 = 4 NO + 6 H2O, without knowing the
    coefficients:

    >>> matrix = stoichiometric_matrix([{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False])
    >>> matrix
    [[3, 0, 0, -2], [1, 0, -1, 0], [0, 2, -1, -1]]
    >>> balance_stoichiometry(matrix)
    [4.0, 5.0, 4.0, 6.0]
    >>> balance_stoichiometry(matrix, allow_fractional=True)
    [1.0, 1.25, 1.0, 1.5]

    The behavior of this function for inputs which do not have a unique
    solution is undefined.

    This algorithm may suffer from floating point issues. If you believe there
    is an error in the result, please report your reaction to the developers.
    This function has a comprehensive test suite and extra test cases can be added to it.

    References
    ----------
    .. [1] Sen, S. K., Hans Agarwal, and Sagar Sen. "Chemical Equation
       Balancing: An Integer Programming Approach." Mathematical and Computer
       Modelling 44, no. 7 (October 1, 2006): 678-91.
       https://doi.org/10.1016/j.mcm.2006.02.004.
    .. [2] URAVNOTE, NOVOODKRITI PARADOKSI V. TEORIJI, and ENJA KEMIJSKIH
       REAKCIJ. "New Discovered Paradoxes in Theory of Balancing Chemical
       Reactions." Materiali in Tehnologije 45, no. 6 (2011): 503-22.
    .. [3] Risteski, Ice B. "A New Approach to Balancing Chemical Equations."
       SIAM Problems & Solutions, 2007, 1-10.
    .. [4] Smith, William R., and Ronald W. Missen. "Using Mathematica and 
       Maple To Obtain Chemical Equations." Journal of Chemical Education
       74, no. 11 (November 1, 1997): 1369. https://doi.org/10.1021/ed074p1369.
    '''
    from fluids.numerics import null_space
    null_vectors = null_space(matrix, rcond=None)
    
    if not null_vectors or len(null_vectors[0]) == 0:
        raise ValueError("No solution found")
    
    # Take the first null vector (assuming unique solution)
    d = [row[0] for row in null_vectors]
    min_value_inv = 1.0/min(d, key=abs)
    d = [i*min_value_inv for i in d]
    d = [round_to_significant(v, rounding) for v in d]
    if not allow_fractional:
        d = floats_to_ints(d, matrix)
    d = [float(v) for v in d]
    return d

def stoichiometry_molar_to_mass(coefficients, MWs):
    r'''This function translates molar stoichiometric
    coefficients (most commonly used) into less commonly
    used mass-based stoichiometric coefficients.

    Parameters
    ----------
    coefficients : list[float]
        Molar balanced stoichiometric coefficients; all numbers are positive, [-]
    MWs : list[float]
        Molecular weights of all species in reaction ordered in
        the same way as the coefficients, [g/mol]

    Returns
    -------
    mass_coefficients : list[float]
        Mass-based balanced coefficients; all numbers are positive, [-]

    Notes
    -----
    Note that mass-based reactions are usually not normalized to integers.
    Mass-based coefficients are used with components that don't have well defined formulas.

    Calculate the mass based coefficients for the reaction 4 NH3 + 5 O2 = 4 NO + 6 H2O:

    >>> matrix = stoichiometric_matrix([{'N': 1, 'H': 3}, {'O': 2}, {'N': 1, 'O': 1}, {'H': 2, 'O': 1}], [True, True, False, False])
    >>> coeffs = balance_stoichiometry(matrix)
    >>> stoichiometry_molar_to_mass(coeffs, [17.03052, 31.9988, 30.0061, 18.01528])
    [68.12208, 159.994, 120.0244, 108.09168]
    '''
    return [c*MW for c, MW in zip(coefficients, MWs)]

def stoichiometry_mass_to_molar(mass_coefficients, MWs):
    r'''This function translates mass stoichiometric coefficients into the
    more commonly used mole-based stoichiometric coefficients.

    Parameters
    ----------
    mass_coefficients : list[float]
        Mass-based balanced coefficients; all numbers are positive, [-]
    MWs : list[float]
        Molecular weights of all species in reaction ordered in
        the same way as the coefficients, [g/mol]

    Returns
    -------
    coefficients : list[float]
        Molar balanced stoichiometric coefficients; all numbers are positive, [-]

    Notes
    -----
    >>> stoichiometry_mass_to_molar([68.12208, 159.994, 120.0244, 108.09168], [17.03052, 31.9988, 30.0061, 18.01528])
    [4.0, 5.0, 4.0, 6.0]
    '''
    return [c/MW for c, MW in zip(mass_coefficients, MWs)]

def stoichiometry_MW_error(coefficients, MWs, reactants):
    r'''This function calculates the molecular weight imbalance
    of a reaction given the coefficients and molecular weights of
    the involved components, and their statuses as reactants or product.

    Parameters
    ----------
    coefficients : list[float]
        Molar balanced stoichiometric coefficients; all numbers are positive, [-]
    MWs : list[float]
        Molecular weights of all species in reaction ordered in
        the same way as the coefficients, [g/mol]
    reactants : list[bool]
        List of booleans indicating whether each chemical is a reactant (True)
        or a product (False), [-]

    Returns
    -------
    MW_error : float
        The molecular weight error, [g/mol]

    Notes
    -----
    A very small value may be returned for a properly balanced
    equation because of floating-point error.

    >>> stoichiometry_MW_error([4.0, 5.0, 4.0, 6.0], [17.03052, 31.9988, 30.0061, 18.01528], [True, True, False, False])
    0.0
    '''
    reactant_MW = 0.0
    product_MW = 0.0
    for coeff, MW, stat in zip(coefficients, MWs, reactants):
        if stat:
            reactant_MW += coeff*MW
        else:
            product_MW += coeff*MW
    return reactant_MW - product_MW

def standard_formation_reaction(atoms):
    r'''This function calculates the standard reaction to reduce a chemical
    compound to its standard state elements. Any hydrogen in the compound
    is transformed to H2; oxygen to O2; carbon to graphite (single C), calcium
    to Ca, etc.

    Parameters
    ----------
    atoms : dict[(str, float)]
        A dictionary of (element, element_count) pairs for the reacting
        compound, [-]

    Returns
    -------
    reactant_coeff : float
        The coefficient of the reactant; for compounds like CO that do not
        divide evenly, this will be something other than 1 [-]
    elemental_counts : list[float]
        Balanced coefficients of each of the products, [-]
    product_atomss : list[dict[(str, float)]]
        A list of dictionaries of the elements produced, and how many atoms
        of each element are in one unit of the element in its
        standard form. Each dictionary contains a single key:value, with the key
        being the element and the value being either 1 or 2 depending on the
        standard state [-]

    Examples
    --------
    Methane

    >>> standard_formation_reaction({'C': 1, 'H': 4})
    (1.0, [1.0, 2.0], [{'C': 1}, {'H': 2}])

    Carbon monoxide

    >>> standard_formation_reaction({'C': 1, 'O': 1})
    (2.0, [2.0, 1.0], [{'C': 1}, {'O': 2}])

    Methylamine

    >>> standard_formation_reaction({'C': 1, 'H': 5, 'N': 1})
    (2.0, [2.0, 5.0, 1.0], [{'C': 1}, {'H': 2}, {'N': 2}])

    '''
    product_atomss = []
    reactants = []
    for atom in atoms:
        ele = periodic_table[atom]
        ele_atoms = simple_formula_parser(ele.formula_standard)
        product_atomss.append(ele_atoms)
        reactants.append(True)

    atoms_to_process = product_atomss + [atoms]
    reactants.append(False)
    matrix = stoichiometric_matrix(atoms_to_process, reactants)
    coeffs = balance_stoichiometry(matrix)
    reactant_coeff = coeffs[-1]
    elemental_counts = coeffs[:-1]
    return reactant_coeff, elemental_counts, product_atomss
