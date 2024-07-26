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

Global Temperature Potential
----------------------------
.. autofunction:: chemicals.environment.GTP
.. autofunction:: chemicals.environment.GTP_methods
.. autodata:: chemicals.environment.GTP_all_methods


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
           'GWP_methods', 'ODP_methods', 'logP_methods',
           'GTP', 'GTP_methods', 'GTP_all_methods']
from chemicals import data_reader as dr
from chemicals import miscdata
from chemicals.data_reader import (
    data_source,
    database_constant_lookup,
    list_available_methods_from_df,
    list_available_methods_from_df_dict,
    register_df_source,
    retrieve_any_from_df,
    retrieve_any_from_df_dict,
    retrieve_from_df,
    retrieve_from_df_dict,
)
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

### Register data sources and lazy load them
folder = os_path_join(source_path, 'Environment')
register_df_source(folder, 'Official Global Warming Potentials 2007.tsv')
register_df_source(folder, 'Official Global Warming Potentials 2014.tsv')
register_df_source(folder, 'Official Global Warming Potentials 2021.tsv', index_col=1)
register_df_source(folder, 'Ozone Depletion Potentials.tsv')
register_df_source(folder, 'CRC logP table.tsv')
register_df_source(folder, 'Syrres logP data.csv.gz',
                   csv_kwargs={'compression': 'gzip'})

IPCC_2007_100YR_GWP = 'IPCC (2007) 100yr'
IPCC_1995_100YR_GWP = 'IPCC (1995) 100yr'
IPCC_2007_20YR_GWP = 'IPCC (2007) 20yr'
IPCC_2007_500YR_GWP = 'IPCC (2007) 500yr'

IPCC_2014_20YR_GWP = 'IPCC (2014) 20yr'
IPCC_2014_100YR_GWP = 'IPCC (2014) 100yr'

IPCC_2021_20YR_GWP = 'IPCC (2021) 20yr'
IPCC_2021_100YR_GWP = 'IPCC (2021) 100yr'
IPCC_2021_500YR_GWP = 'IPCC (2021) 500yr'

IPCC_2014_20YR_GTP = 'IPCC (2014) 20yr'
IPCC_2014_50YR_GTP = 'IPCC (2014) 50yr'
IPCC_2014_100YR_GTP = 'IPCC (2014) 100yr'

IPCC_2021_50YR_GTP = 'IPCC (2021) 50yr'
IPCC_2021_100YR_GTP = 'IPCC (2021) 100yr'

GWP_all_methods = (IPCC_2014_100YR_GWP, IPCC_2014_20YR_GWP,
                   IPCC_2007_100YR_GWP, IPCC_2007_20YR_GWP, IPCC_2007_500YR_GWP,
                   IPCC_1995_100YR_GWP,
                   IPCC_2021_20YR_GWP, IPCC_2021_100YR_GWP, IPCC_2021_500YR_GWP)
"""Tuple of method name keys. See the `GWP` for the actual references"""

GTP_all_methods = (IPCC_2014_20YR_GTP, IPCC_2014_50YR_GTP, IPCC_2014_100YR_GTP, IPCC_2021_50YR_GTP, IPCC_2021_100YR_GTP)
"""Tuple of method name keys. See the `GTP` for the actual references"""

_GWP_ODP_data_loaded = False
@mark_numba_incompatible
def _load_GWP_ODP_data():
    global _GWP_ODP_data_loaded, IPCC_2007_GWPs, IPCC_2014_GWPs, IPCC_2021_GWPs, ODP_data
    global _IPCC_2007_GWP_keys_by_method, _IPCC_2014_GWP_keys_by_method, _IPCC_2021_GWP_keys_by_method
    global _IPCC_2014_GTP_keys_by_method, _IPCC_2021_GTP_keys_by_method
    global _ODP_keys_by_method
    IPCC_2007_GWPs = data_source('Official Global Warming Potentials 2007.tsv')
    IPCC_2014_GWPs = data_source('Official Global Warming Potentials 2014.tsv')
    IPCC_2021_GWPs = data_source('Official Global Warming Potentials 2021.tsv')

    ODP_data = data_source('Ozone Depletion Potentials.tsv')
    _GWP_ODP_data_loaded = True
    _IPCC_2007_GWP_keys_by_method = {
        IPCC_2007_20YR_GWP: '20yr GWP',
        IPCC_2007_100YR_GWP : '100yr GWP',
        IPCC_1995_100YR_GWP: 'SAR 100yr',
        IPCC_2007_500YR_GWP: '500yr GWP',
    }
    _IPCC_2014_GWP_keys_by_method = {
        IPCC_2014_20YR_GWP : '20yr GWP',
        IPCC_2014_100YR_GWP: '100yr GWP',
    }
    _IPCC_2021_GWP_keys_by_method = {
        IPCC_2021_20YR_GWP: '20yr GWP',
        IPCC_2021_100YR_GWP: '100yr GWP',
        IPCC_2021_500YR_GWP: '500yr GWP',
    }

    _IPCC_2014_GTP_keys_by_method = {
        IPCC_2014_20YR_GTP : '20yr GTP',
        IPCC_2014_50YR_GTP : '50yr GTP',
        IPCC_2014_100YR_GTP: '100yr GTP',
    }

    _IPCC_2021_GTP_keys_by_method = {
        IPCC_2021_50YR_GTP: '50yr GTP',
        IPCC_2021_100YR_GTP: '100yr GTP',
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
@mark_numba_incompatible
def _load_logP_data():
    global _logP_data_loaded, logP_data_CRC, logP_data_Syrres, logP_sources
    logP_data_CRC = data_source('CRC logP table.tsv')
    logP_data_Syrres = data_source('Syrres logP data.csv.gz')
    _logP_data_loaded = True
    logP_sources = {
        'CRC': logP_data_CRC,
        'SYRRES': logP_data_Syrres,
        miscdata.WIKIDATA: miscdata.wikidata_data
    }

if PY37:
    def __getattr__(name):
        if name in ('IPCC_2007_GWPs', 'IPCC_2014_GWPs', 'IPCC_2021_GWPs', 'ODP_data'):
            _load_GWP_ODP_data()
            return globals()[name]
        elif name in ('logP_data_CRC', 'logP_data_Syrres'):
            _load_logP_data()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:  # pragma: no cover
    if can_load_data:
        _load_GWP_ODP_data()
        _load_logP_data()


### Environmental data functions

@mark_numba_incompatible
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
    methods_4e = list_available_methods_from_df(IPCC_2007_GWPs, CASRN, _IPCC_2007_GWP_keys_by_method)
    methods_5e = list_available_methods_from_df(IPCC_2014_GWPs, CASRN, _IPCC_2014_GWP_keys_by_method)
    methods_6e = list_available_methods_from_df(IPCC_2021_GWPs, CASRN, _IPCC_2021_GWP_keys_by_method)
    return methods_4e + methods_5e + methods_6e

@mark_numba_incompatible
def GTP_methods(CASRN):
    """Return all methods available to obtain GTP for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain GTP with the given inputs.

    See Also
    --------
    GTP
    """
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    methods_5e = list_available_methods_from_df(IPCC_2014_GWPs, CASRN, _IPCC_2014_GTP_keys_by_method)
    methods_6e = list_available_methods_from_df(IPCC_2021_GWPs, CASRN, _IPCC_2021_GTP_keys_by_method)
    return methods_5e + methods_6e

@mark_numba_incompatible
def GWP(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's Global Warming
    Potential, relative to CO2. Lookup is based on CASRNs.

    There are four sources of data:

        * IPCC Sixth Assessment Report (AR5) from 2021 [4]_
        * IPCC Fifth Assessment Report (AR5) from 2014 [2]_
        * IPCC Fourth Assessment Report (AR4) from 2007 [1]_
        * IPCC Second Assesment Report or (SAR) from 1995 [1]_

    This function returns the GWP for the 20yr outlook from the AR6 by default.

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
        The method name to use. Accepted methods are ('IPCC (2021) 100yr',
        'IPCC (2021) 20yr', 'IPCC (2021) 500yr', 'IPCC (2014) 100yr',
        'IPCC (2014) 20yr', 'IPCC (2007) 100yr', 'IPCC (2007) 20yr',
        'IPCC (2007) 500yr', 'IPCC (1995) 100yr').
        All valid values are also held in the variable `GWP_all_methods`.

    Notes
    -----
    "Fossil methane" is included in the IPCC reports to take into account
    different isotopic composition, but as that has the same CAS number it
    is not included in this function.

    Six of the entries in [2]_ are actually duplicates; the entries
    with data similar to more recent data [3]_ were prefered.

    Examples
    --------
    Methane, 20-yr outlook AR6

    >>> GWP(CASRN='74-82-8')
    81.2

    Methane, specifying the default method explicitly (this is recommended
    the default data source may be updated in the future)

    >>> GWP(CASRN='74-82-8', method='IPCC (2014) 100yr')
    28.0

    Methane, 20-year values from 1995 and 2007

    >>> (GWP(CASRN='74-82-8', method='IPCC (1995) 100yr'), GWP(CASRN='74-82-8', method='IPCC (2007) 100yr'))
    (21.0, 25.0)

    See Also
    --------
    GWP_methods

    References
    ----------
    .. [1] IPCC. "2.10.2 Direct Global Warming Potentials - AR4 WGI Chapter 2:
       Changes in Atmospheric Constituents and in Radiative Forcing." 2007.
       https://www.ipcc.ch/publications_and_data/ar4/wg1/en/ch2s2-10-2.html.
    .. [2] IPCC. "Climate Change 2013: The Physical Science Basis. - AR5 WGI Chapter 8:
       Anthropogenic and Natural Radiative Forcing." 2013.
       https://www.ipcc.ch/site/assets/uploads/2018/02/WG1AR5_Chapter08_FINAL.pdf
    .. [3] Hodnebrog, Ø., B. Aamaas, J. S. Fuglestvedt, G. Marston, G. Myhre, C.
       J. Nielsen, M. Sandstad, K. P. Shine, and T. J. Wallington. "Updated
       Global Warming Potentials and Radiative Efficiencies of Halocarbons and
       Other Weak Atmospheric Absorbers." Reviews of Geophysics 58, no. 3
       (2020): e2019RG000691. https://doi.org/10.1029/2019RG000691.
    .. [4] Masson-Delmotte, Valérie, Panmao Zhai, Anna Pirani, Sarah L. Connors,
       Clotilde Péan, Sophie Berger, Nada Caud, Yang Chen, Leah Goldfarb, and Melissa
       I. Gomis. "Climate Change 2021: The Physical Science Basis." Contribution of
       Working Group I to the Sixth Assessment Report of the Intergovernmental
       Panel on Climate Change 2 (2021): 24.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'GWP')
        if found: return val
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if method:
        if method in _IPCC_2021_GWP_keys_by_method:
            key, df = _IPCC_2021_GWP_keys_by_method[method], IPCC_2021_GWPs
        if method in _IPCC_2014_GWP_keys_by_method:
            key, df = _IPCC_2014_GWP_keys_by_method[method], IPCC_2014_GWPs
        elif method in _IPCC_2007_GWP_keys_by_method:
            key, df = _IPCC_2007_GWP_keys_by_method[method], IPCC_2007_GWPs
        return retrieve_from_df(df, CASRN, key)
    else:
        try:
            return retrieve_any_from_df(IPCC_2021_GWPs, CASRN, _IPCC_2021_GWP_keys_by_method.values())
        except:
            try:
                return retrieve_any_from_df(IPCC_2014_GWPs, CASRN, _IPCC_2014_GWP_keys_by_method.values())
            except:
                try:
                    return retrieve_any_from_df(IPCC_2007_GWPs, CASRN, _IPCC_2007_GWP_keys_by_method.values())
                except:
                    return None

@mark_numba_incompatible
def GTP(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's Global Temperature
    Potential, relative to CO2. Lookup is based on CASRNs.

    There are two sources of data:

        * IPCC Sixth Assessment Report (AR5) from 2021 [2]_
        * IPCC Fifth Assessment Report (AR5) from 2014 [1]_

    This function returns the GTP for the 50yr outlook from the AR6 by default.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    GTP : float
        Global temperature potential, [(impact/mass chemical)/(impact/mass CO2)]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are ('IPCC (2021) 100yr',
        'IPCC (2021) 50yr', 'IPCC (2014) 100yr',
        'IPCC (2014) 20yr', 'IPCC (2014) 50yr').
        All valid values are also held in the variable `GTP_all_methods`.

    Notes
    -----

    Examples
    --------
    Methane, 50-yr outlook AR6

    >>> GTP(CASRN='74-82-8')
    11.0


    See Also
    --------
    GTP_methods

    References
    ----------
    .. [1] IPCC. "Climate Change 2013: The Physical Science Basis. - AR5 WGI Chapter 8:
       Anthropogenic and Natural Radiative Forcing." 2013.
       https://www.ipcc.ch/site/assets/uploads/2018/02/WG1AR5_Chapter08_FINAL.pdf
    .. [2] Masson-Delmotte, Valérie, Panmao Zhai, Anna Pirani, Sarah L. Connors,
       Clotilde Péan, Sophie Berger, Nada Caud, Yang Chen, Leah Goldfarb, and Melissa
       I. Gomis. "Climate Change 2021: The Physical Science Basis." Contribution of
       Working Group I to the Sixth Assessment Report of the Intergovernmental
       Panel on Climate Change 2 (2021): 24.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'GTP')
        if found: return val
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if method:
        if method in _IPCC_2021_GTP_keys_by_method:
            key, df = _IPCC_2021_GTP_keys_by_method[method], IPCC_2021_GWPs
        elif method in _IPCC_2014_GTP_keys_by_method:
            key, df = _IPCC_2014_GTP_keys_by_method[method], IPCC_2014_GWPs
        return retrieve_from_df(df, CASRN, key)
    else:
        try:
            return retrieve_any_from_df(IPCC_2021_GWPs, CASRN, _IPCC_2021_GTP_keys_by_method.values())
        except:
            try:
                return retrieve_any_from_df(IPCC_2014_GWPs, CASRN, _IPCC_2014_GTP_keys_by_method.values())
            except:
                return None

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
"""Tuple of method name keys. See the `ODP` for the actual references"""

@mark_numba_incompatible
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
    return list_available_methods_from_df(ODP_data, CASRN, _ODP_keys_by_method)

@mark_numba_incompatible
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
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'ODP')
        if found: return val
    if not _GWP_ODP_data_loaded: _load_GWP_ODP_data()
    if method:
        key = _ODP_keys_by_method[method]
        return retrieve_from_df(ODP_data, CASRN, key)
    else:
        return retrieve_any_from_df(ODP_data, CASRN, _ODP_keys_by_method.values())

### log P

SYRRES = 'SYRRES'
CRC = 'CRC'
logP_all_methods = (SYRRES, CRC, miscdata.WIKIDATA)
"""Tuple of method name keys. See the `logP` for the actual references"""

@mark_numba_incompatible
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

@mark_numba_incompatible
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
        The method name to use. Accepted methods are 'SYRRES', 'CRC', and 'WIKIDATA'.
        All valid values are also held in the list logP_methods.

    Notes
    -----
    Although matimatically this could be expressed with a logarithm in any
    base, reported values are published using a  base 10 logarithm.

    .. math::
        \log_{10} P_{ oct/wat} = \log_{10}\left(\frac{\left[{solute}
        \right]_{ octanol}^{un-ionized}}{\left[{solute}
        \right]_{ water}^{ un-ionized}}\right)

    Examples
    --------
    >>> logP('67-56-1')
    -0.74
    >>> logP('100-66-3', 'WIKIDATA')
    2.11

    References
    ----------
    .. [1] Syrres. 2006. KOWWIN Data, SrcKowData2.zip.
       http://esc.syrres.com/interkow/Download/SrcKowData2.zip
    .. [2] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'logP')
        if found: return val
    if not _logP_data_loaded: _load_logP_data()
    if method:
        return retrieve_from_df_dict(logP_sources, CASRN, 'logP', method)
    else:
        return retrieve_any_from_df_dict(logP_sources, CASRN, 'logP')
