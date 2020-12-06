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
"""

__all__ = ['Hfg', 'Hfl', 'Hfs', 'S0g', 'S0l', 'S0s',
           'Hfl_methods', 'Hfg_methods', 'Hfs_methods',
           'S0l_methods', 'S0g_methods', 'S0s_methods',
           'Hfl_all_methods', 'Hfg_all_methods', 'Hfs_all_methods',
           'S0l_all_methods', 'S0g_all_methods', 'S0s_all_methods',
           'Gibbs_formation', 'entropy_formation', 'Hf_basis_converter',
           'balance_stoichiometry', 'stoichiometric_matrix']

from chemicals.utils import ceil, log10, PY37, source_path, os_path_join, can_load_data
from chemicals import heat_capacity
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)


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
    _reaction_data_loaded = True
    S0g_sources = {
        CRC: heat_capacity.CRC_standard_data,
        YAWS: Hfg_S0g_YAWS_data,
    }
    S0l_sources = {
        CRC: heat_capacity.CRC_standard_data,
    }
    S0s_sources = {
        CRC: heat_capacity.CRC_standard_data,
    }
    Hfg_sources = {
        ATCT_G: Hfg_ATcT_data,
        CRC: heat_capacity.CRC_standard_data,
        API_TDB_G: Hfg_API_TDB_data,
        TRC: heat_capacity.TRC_gas_data,
        YAWS: Hfg_S0g_YAWS_data,
    }
    Hfl_sources = {
        ATCT_L: Hfl_ATcT_data,
        CRC: heat_capacity.CRC_standard_data,
    }
    Hfs_sources = {
        CRC: heat_capacity.CRC_standard_data,
    }

if PY37:
    def __getattr__(name):
        if name in ('Hfg_API_TDB_data', 'Hfg_ATcT_data',
                    'Hfl_ATcT_data', 'Hfg_S0g_YAWS_data',
                    'Hfg_sources', 'Hfl_sources', 'Hfs_sources',
                    'S0g_sources', 'S0l_sources', 'S0s_sources'):
            _load_reaction_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_reaction_data()


# %% Lookup functions

# TODO: more data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3692305/
# has dippr standard heats of formation, about 55% of the database

Hfs_all_methods = (CRC,)
'''Tuple of method name keys. See the `Hfs` for the actual references'''

def Hfs_methods(CASRN):
    """Return all methods available to obtain the Hfs for the desired chemical.

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

        * 'CRC', from the CRC handbook (1360 values)

    Examples
    --------
    >>> Hfs('101-81-5') # Diphenylmethane
    71500.0

    See Also
    --------
    Hfs_methods

    References
    ----------
    .. [1] Ruscic, Branko, Reinhardt E. Pinzon, Gregor von Laszewski, Deepti
       Kodeboyina, Alexander Burcat, David Leahy, David Montoy, and Albert F.
       Wagner. "Active Thermochemical Tables: Thermochemistry for the 21st
       Century." Journal of Physics: Conference Series 16, no. 1
       (January 1, 2005): 561. doi:10.1088/1742-6596/16/1/078.
    '''
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfs_sources, CASRN, 'Hfs', method)
    else:
        return retrieve_any_from_df_dict(Hfs_sources, CASRN, 'Hfs')

Hfl_all_methods = (ATCT_L, CRC)
'''Tuple of method name keys. See the `Hfl` for the actual references'''

def Hfl_methods(CASRN):
    """Return all methods available to obtain the Hfl for the desired chemical.

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

        * 'ATCT_L', the Active Thermochemical Tables version 1.112.
        * 'CRC', from the CRC handbook (1360 values)

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
    '''
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfl_sources, CASRN, 'Hfl', method)
    else:
        return retrieve_any_from_df_dict(Hfl_sources, CASRN, 'Hfl')

Hfg_all_methods = (ATCT_G, TRC, CRC, YAWS)
'''Tuple of method name keys. See the `Hfg` for the actual references'''

def Hfg_methods(CASRN):
    """Return all methods available to obtain the Hfg for the desired chemical.

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

        * 'ATCT_G', the Active Thermochemical Tables version 1.112 (600 values)
        * 'TRC', from a 1994 compilation (1750 values)
        * 'CRC', from the CRC handbook (1360 values)
        * 'YAWS', a large compillation of values, mostly estimated (5000 values)

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
    .. [2] Frenkelʹ, M. L, Texas Engineering Experiment Station, and
       Thermodynamics Research Center. Thermodynamics of Organic Compounds in
       the Gas State. College Station, Tex.: Thermodynamics Research Center,
       1994.
    .. [3] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    .. [4] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(Hfg_sources, CASRN, 'Hfg', method)
    else:
        return retrieve_any_from_df_dict(Hfg_sources, CASRN, 'Hfg')

S0s_all_methods = (CRC,)
'''Tuple of method name keys. See the `S0s` for the actual references'''

def S0s_methods(CASRN):
    """Return all methods available to obtain the S0s for the desired chemical.

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

        * 'CRC', from the CRC handbook (1360 values)

    Examples
    --------
    >>> S0s('7439-93-2') # Lithium
    29.1

    See Also
    --------
    S0s_methods

    '''
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(S0s_sources, CASRN, 'S0s', method)
    else:
        return retrieve_any_from_df_dict(S0s_sources, CASRN, 'S0s')

S0l_all_methods = (CRC,)
'''Tuple of method name keys. See the `S0l` for the actual references'''

def S0l_methods(CASRN):
    """Return all methods available to obtain the S0l for the desired chemical.

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
    if not _reaction_data_loaded: _load_reaction_data()
    if method:
        return retrieve_from_df_dict(S0l_sources, CASRN, 'S0l', method)
    else:
        return retrieve_any_from_df_dict(S0l_sources, CASRN, 'S0l')

S0g_all_methods = (CRC, YAWS)
'''Tuple of method name keys. See the `S0g` for the actual references'''

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
    '''
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
    elements = sorted(list(elements)) # Ensure reproducibility
    n_elements = len(elements)

    matrix = [[0]*n_compounds for _ in range(n_elements)]
    for i, atoms in enumerate(atomss):
        for k, v in atoms.items():
            if not reactants[i]:
                v = -v
            matrix[elements.index(k)][i] = v
    return matrix

def balance_stoichiometry(matrix, rounding=9, allow_fractional=False):
    r'''This function balances a chemical reaction.

    Parameters
    ----------
    matrix : list[list[float]]
        Chemical reaction matrix for further processing; rows contain element
         counts of each compound, and the columns represent each chemical, [-]

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

    This algorithm relies on `scipy`.
    The behavior of this function for inputs which do not have a unique
    solution is undefined.

    This algorithm may suffer from floating point issues. If you believe there
    is an error in the result, please report your reaction to the developers.

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
    import scipy.linalg
    done = scipy.linalg.null_space(matrix)
    if len(done[0]) > 1:
        raise ValueError("No solution")
    d = done[:, 0].tolist()

    min_value_inv = 1.0/min(d)
    d = [i*min_value_inv for i in d]

    if not allow_fractional:
        from fractions import Fraction
        max_denominator = 10**rounding
        fs = [Fraction(x).limit_denominator(max_denominator=max_denominator) for x in d]
        all_denominators = set([i.denominator for i in fs])
        if 1 in all_denominators:
            all_denominators.remove(1)

        for den in sorted(list(all_denominators), reverse=True):
            fs = [num*den for num in fs]
            if all(i.denominator == 1 for i in fs):
                break

        # May have gone too far
        return [float(i) for i in fs]
#        done = False
#        for i in range(100):
#            for c in d:
#                ratio = c.as_integer_ratio()[1]
#                if ratio != 1:
#                    d = [di*ratio for di in d]
#                    break
#                done = True
#            if done:
#                break
#
#        d_as_int = [int(i) for i in d]
#        for i, j in zip(d, d_as_int):
#            if i != j:
#                raise ValueError("Could not find integer coefficients (%s, %s)" %(i, j))
#        return d_as_int
    else:
        d = [round(i, rounding + int(ceil(log10(abs(i))))) for i in d]
        return d


