# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell
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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This module contains lookup functions for three important environmental 
properties - Global Warming Potential, Ozone Depletion Potential, and 
octanol-water partition coefficient.


For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Global Warming Potential
------------------------
.. autofunction:: chemicals.environment.GWP
.. autofunction:: chemicals.environment.GWP_methods
.. autodata:: chemicals.environment.GWP_all_methods

Ozone Depletion Potential
-------------------------
.. autofunction:: chemicals.environment.ODP
.. autofunction:: chemicals.environment.ODP_methods
.. autodata:: chemicals.environment.ODP_all_methods

Octanol-Water Partition Coefficient
-----------------------------------
.. autofunction:: chemicals.environment.logP
.. autofunction:: chemicals.environment.logP_methods
.. autodata:: chemicals.environment.logP_all_methods

"""

__all__ = ['GWP', 'ODP', 'logP', 
           'GWP_all_methods', 'ODP_all_methods', 'logP_all_methods',
           'GWP_methods', 'ODP_methods', 'logP_methods']

import os
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df,
                                   retrieve_any_from_df,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

### Register data sources and lazy load them

folder = os_path_join(source_path, 'Environment')
register_df_source(folder, 'Official Global Warming Potentials.tsv')
register_df_source(folder, 'Ozone Depletion Potentials.tsv')
register_df_source(folder, 'CRC logP table.tsv')
register_df_source(folder, 'Syrres logP data.csv.gz',
                   csv_kwargs={'compression': 'gzip'})

_GWP_ODP_data_loaded = False
def _load_GWP_ODP_data():
    global _GWP_ODP_data_loaded, GWP_data, ODP_data, _GWP_keys_by_method
    global _GWP_keys_by_method, _ODP_keys_by_method
    GWP_data = data_source('Official Global Warming Potentials.tsv')
    ODP_data = data_source('Ozone Depletion Potentials.tsv')
    _GWP_ODP_data_loaded = True
    _GWP_keys_by_method = {
        'IPCC (2007) 100yr' : '100yr GWP',
        'IPCC (2007) 100yr-SAR': 'SAR 100yr',
        'IPCC (2007) 20yr': '20yr GWP',
        'IPCC (2007) 500yr': '500yr GWP',
    }
    _ODP_keys_by_method = {
        'ODP2 Max': 'ODP2 Max',
        'ODP1 Max': 'ODP1 Max',
        'ODP2 logarithmic average': 'ODP2 Design',
        'ODP1 logarithmic average': 'ODP1 Design',
        'ODP2 Min': 'ODP2 Min',
        'ODP1 Min': 'ODP1 Min',
        'ODP2 string': 'ODP2',
        'ODP1 string': 'ODP1',
    }
    
_logP_data_loaded = False
def _load_logP_data():
    global _logP_data_loaded, logP_data_CRC, logP_data_Syrres, logP_sources
    logP_data_CRC = data_source('CRC logP table.tsv')
    logP_data_Syrres = data_source('Syrres logP data.csv.gz')
    _logP_data_loaded = True
    logP_sources = {
        'CRC': logP_data_CRC,
        'SYRRES': logP_data_Syrres,
    }

if PY37:
    def __getattr__(name):
        if name in ('GWP_data', 'ODP_data'):
            _load_GWP_ODP_data()
            return globals()[name]
        elif name in ('logP_data_CRC', 'logP_data_Syrres'):
            _load_logP_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:  # pragma: no cover
    if can_load_data:
        _load_GWP_ODP_data()
        _load_logP_data()

IPCC100 = 'IPCC (2007) 100yr'
IPCC100SAR = 'IPCC (2007) 100yr-SAR'
IPCC20 = 'IPCC (2007) 20yr'
IPCC500 = 'IPCC (2007) 500yr'
GWP_all_methods = (IPCC100, IPCC100SAR, IPCC20, IPCC500)
'''Tuple of method name keys. See the `GWP` for the actual references'''


### Environmental data functions

def GWP_methods(CASRN):
    """Return all methods available to obtain GWP for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain GWP with the given inputs.

    See Also
    --------
    GWP
    """
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if CASRN in GWP_data.index:
        import pandas as pd
        return [method for method, key in _GWP_keys_by_method.items()
                if not pd.isnull(GWP_data.at[CASRN, key])]
    else:
        return []

def GWP(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's Global Warming
    Potential, relative to CO2. Lookup is based on CASRNs. Will automatically
    select a data source to use if no method is provided; returns None if the
    data is not available.

    Returns the GWP for the 100yr outlook by default.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    GWP : float
        Global warming potential, [(impact/mass chemical)/(impact/mass CO2)]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are IPCC (2007) 100yr',
        'IPCC (2007) 100yr-SAR', 'IPCC (2007) 20yr', and 'IPCC (2007) 500yr'. 
        All valid values are also held in the variable `GWP_all_methods`.

    Notes
    -----
    All data is from [1]_, the official source. Several chemicals are available
    in [1]_ are not included here as they do not have a CAS.
    Methods are 'IPCC (2007) 100yr', 'IPCC (2007) 100yr-SAR',
    'IPCC (2007) 20yr', and 'IPCC (2007) 500yr'.

    Examples
    --------
    Methane, 100-yr outlook

    >>> GWP(CASRN='74-82-8')
    25.0

    See Also
    --------
    GWP_methods

    References
    ----------
    .. [1] IPCC. "2.10.2 Direct Global Warming Potentials - AR4 WGI Chapter 2:
       Changes in Atmospheric Constituents and in Radiative Forcing." 2007.
       https://www.ipcc.ch/publications_and_data/ar4/wg1/en/ch2s2-10-2.html.
    '''
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if method:
        key = _GWP_keys_by_method[method]
        return retrieve_from_df(GWP_data, CASRN, key)
    else:
        return retrieve_any_from_df(GWP_data, CASRN, _GWP_keys_by_method.values())

### Ozone Depletion Potentials

ODP2MAX = 'ODP2 Max'
ODP2MIN = 'ODP2 Min'
ODP2STR = 'ODP2 string'
ODP2LOG = 'ODP2 logarithmic average'
ODP1MAX = 'ODP1 Max'
ODP1MIN = 'ODP1 Min'
ODP1STR = 'ODP1 string'
ODP1LOG = 'ODP1 logarithmic average'
ODP_all_methods = (ODP2MAX, ODP1MAX, ODP2LOG, ODP1LOG,
                   ODP2MIN, ODP1MIN, ODP2STR, ODP1STR)
'''Tuple of method name keys. See the `ODP` for the actual references'''

def ODP_methods(CASRN):
    """Return all methods available to obtain ODP for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain ODP with the given inputs.

    See Also
    --------
    ODP
    """
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if CASRN in ODP_data.index:
        import pandas as pd
        return [method for method, key in _ODP_keys_by_method.items()
                if not pd.isnull(ODP_data.at[CASRN, key])]
    else:
        return []

def ODP(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's Ozone Depletion
    Potential, relative to CFC-11 (trichlorofluoromethane). Lookup is based on
    CASRNs. Will automatically select a data source to use if no method is
    provided; returns None if the data is not available.

    Returns the ODP of a chemical according to [2]_ when a method is not
    specified. If a range is provided in [2]_, the highest value is returned.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    ODP : float or str
        Ozone Depletion potential, [(impact/mass chemical)/(impact/mass CFC-11)];
        if method selected has `string` in it, this will be returned as a
        string regardless of if a range is given or a number

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'ODP2 Max', 'ODP2 Min', 
        'ODP2 string', 'ODP2 logarithmic average', and methods for older values
        are 'ODP1 Max', 'ODP1 Min', 'ODP1 string', and 'ODP1 logarithmic average'.
        All valid values are also held in the list ODP_methods.

    Notes
    -----
    Values are tabulated only for a small number of halogenated hydrocarbons,
    responsible for the largest impact. The original values of ODP as defined
    in the Montreal Protocol are also available, as methods with the `ODP1`
    prefix.

    All values are somewhat emperical, as actual reaction rates of chemicals
    with ozone depend on temperature which depends on latitude, longitude,
    time of day, weather, and the concentrations of other pollutants.

    All data is from [1]_. Several mixtures listed in [1]_ are not included
    here as they are not pure species.
    Methods for values in [2]_ are 'ODP2 Max', 'ODP2 Min', 'ODP2 string',
    'ODP2 logarithmic average',  and methods for older values are 'ODP1 Max',
    'ODP1 Min', 'ODP1 string', and 'ODP1 logarithmic average'.

    Examples
    --------
    Dichlorotetrafluoroethane, according to [2]_.

    >>> ODP(CASRN='76-14-2')
    0.58

    References
    ----------
    .. [1] US EPA, OAR. "Ozone-Depleting Substances." Accessed April 26, 2016.
       https://www.epa.gov/ozone-layer-protection/ozone-depleting-substances.
    .. [2] WMO (World Meteorological Organization), 2011: Scientific Assessment
       of Ozone Depletion: 2010. Global Ozone Research and Monitoring
       Project-Report No. 52, Geneva, Switzerland, 516 p.
       https://www.wmo.int/pages/prog/arep/gaw/ozone_2010/documents/Ozone-Assessment-2010-complete.pdf
    '''
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if method:
        key = _ODP_keys_by_method[method]
        return retrieve_from_df(ODP_data, CASRN, key)
    else:
        return retrieve_any_from_df(ODP_data, CASRN, _ODP_keys_by_method.values())

### log P

SYRRES = 'SYRRES'
CRC = 'CRC'
logP_all_methods = (SYRRES, CRC)
'''Tuple of method name keys. See the `logP` for the actual references'''

def logP_methods(CASRN):
    """Return all methods available to obtain logP for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain logP with the given inputs.

    See Also
    --------
    logP
    """
    if not _logP_data_loaded: _load_logP_data()
    return list_available_methods_from_df_dict(logP_sources, CASRN, 'logP')

def logP(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's octanol-water
    partition coefficient. Lookup is based on CASRNs. Will automatically
    select a data source to use if no method is provided; returns None if the
    data is not available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    logP : float
        Octanol-water partition coefficient, [-]
        
    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'SYRRES', or 'CRC', 
        All valid values are also held in the list logP_methods.

    Notes
    -----
    .. math::
        \log P_{ oct/wat} = \log\left(\frac{\left[{solute}
        \right]_{ octanol}^{un-ionized}}{\left[{solute}
        \right]_{ water}^{ un-ionized}}\right)

    Examples
    --------
    >>> logP('67-56-1')
    -0.74

    References
    ----------
    .. [1] Syrres. 2006. KOWWIN Data, SrcKowData2.zip.
       http://esc.syrres.com/interkow/Download/SrcKowData2.zip
    .. [2] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    if not _logP_data_loaded: _load_logP_data()
    if method:
        return retrieve_from_df_dict(logP_sources, CASRN, 'logP', method) 
    else:
        return retrieve_any_from_df_dict(logP_sources, CASRN, 'logP') 
