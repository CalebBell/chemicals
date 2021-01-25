# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

This module contains functions for lookup the following properties for a
chemical:

* Short-term Exposure Limit (STEL)
* Time-Weighted Average Exposure Limit (TWA)
* Celing limit for working exposure
* Whether a chemicals is absorbed thorough human skin
* Whether a chemical is a carcinogen, suspected of being a carcinogen, or has
  been identified as unlikely to be a carcinogen

* Flash point
* Auto ignition point
* Lower flammability limit
* Upper flammability limit

In addition, several estimation methods for chemicals without flammability
limits are provided and for calculating the flammability limits of mixtures.

This module also contains several utility functions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Short-term Exposure Limit
-------------------------
.. autofunction:: chemicals.safety.STEL
.. autofunction:: chemicals.safety.STEL_methods
.. autodata:: chemicals.safety.STEL_all_methods

Time-Weighted Average Exposure Limit
------------------------------------
.. autofunction:: chemicals.safety.TWA
.. autofunction:: chemicals.safety.TWA_methods
.. autodata:: chemicals.safety.TWA_all_methods

Ceiling Limit
-------------
.. autofunction:: chemicals.safety.Ceiling
.. autofunction:: chemicals.safety.Ceiling_methods
.. autodata:: chemicals.safety.Ceiling_all_methods

Skin Absorbance
---------------
.. autofunction:: chemicals.safety.Skin
.. autofunction:: chemicals.safety.Skin_methods
.. autodata:: chemicals.safety.Skin_all_methods

Carcinogenicity
---------------
.. autofunction:: chemicals.safety.Carcinogen
.. autofunction:: chemicals.safety.Carcinogen_methods
.. autodata:: chemicals.safety.Carcinogen_all_methods

Flash Point
-----------
.. autofunction:: chemicals.safety.T_flash
.. autofunction:: chemicals.safety.T_flash_methods
.. autodata:: chemicals.safety.T_flash_all_methods

Autoignition Point
------------------
.. autofunction:: chemicals.safety.T_autoignition
.. autofunction:: chemicals.safety.T_autoignition_methods
.. autodata:: chemicals.safety.T_autoignition_all_methods

Lower Flammability Limit
------------------------
.. autofunction:: chemicals.safety.LFL
.. autofunction:: chemicals.safety.LFL_methods
.. autodata:: chemicals.safety.LFL_all_methods
.. autofunction:: chemicals.safety.Suzuki_LFL
.. autofunction:: chemicals.safety.Crowl_Louvar_LFL
.. autofunction:: chemicals.safety.LFL_ISO_10156_2017

Upper Flammability Limit
------------------------
.. autofunction:: chemicals.safety.UFL
.. autofunction:: chemicals.safety.UFL_methods
.. autodata:: chemicals.safety.UFL_all_methods
.. autofunction:: chemicals.safety.Suzuki_UFL
.. autofunction:: chemicals.safety.Crowl_Louvar_UFL

Mixture Flammability Limit
--------------------------
.. autofunction:: chemicals.safety.fire_mixing

Utility Methods
---------------
.. autofunction:: chemicals.safety.ppmv_to_mgm3
.. autofunction:: chemicals.safety.mgm3_to_ppmv
.. autofunction:: chemicals.safety.NFPA_30_classification

"""

__all__ = ('ppmv_to_mgm3', 'mgm3_to_ppmv',
           'NTP_codes', 'IARC_codes',
           'Skin_all_methods',  'Ceiling_all_methods', 'STEL_all_methods',
           'TWA_all_methods',
           'TWA_methods', 'TWA', 'STEL', 'STEL_methods', 'Ceiling', 'Ceiling_methods',
           'Skin', 'Skin_methods', 'Carcinogen_methods', 'Carcinogen_all_methods',
           'Carcinogen', 'T_flash_all_methods', 'T_flash_methods',
           'T_flash', 'T_autoignition_methods', 'T_autoignition_all_methods',
           'T_autoignition', 'LFL_methods', 'LFL_all_methods',
           'LFL', 'UFL_methods', 'UFL_all_methods', 'UFL', 'fire_mixing',
           'Suzuki_LFL', 'Suzuki_UFL',
           'Crowl_Louvar_LFL', 'Crowl_Louvar_UFL', 'LFL_ISO_10156_2017',
           'NFPA_30_classification')

import os
from fluids.core import F2K
from chemicals.utils import source_path, R, none_and_length_check, normalize, PY37, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

### Utilities

def ppmv_to_mgm3(ppmv, MW, T=298.15, P=101325.0):
    r'''
    Converts a concentration in ppmv to units of mg/m^3. Used in
    industrial toxicology.

    .. math::
        \frac{mg}{m^3} = \frac{ppmv\cdot P}{RT}\cdot \frac{MW}{1000}

    Parameters
    ----------
    ppmv : float
        Concentration of a component in a gas mixure [parts per million,
        volumetric]
    MW : float
        Molecular weight of the trace gas [g/mol]
    T : float, optional
        Temperature of the gas at which the ppmv is reported, [K]
    P : float, optional
        Pressure of the gas at which the ppmv is reported, [Pa]

    Returns
    -------
    mgm3 : float
        Concentration of a substance in an ideal gas mixture [mg/m^3]

    Notes
    -----
    The term P/(RT)/1000 converts to 0.040874 at STP. Its inverse is reported
    as 24.45 in [1]_.

    Examples
    --------
    >>> ppmv_to_mgm3(1.0, 40.0)
    1.6349617809430446

    References
    ----------
    .. [1] ACGIH. Industrial Ventilation: A Manual of Recommended Practice,
       23rd Edition. American Conference of Governmental and Industrial
       Hygenists, 2004.
    '''
    parts = ppmv*1E-6
    n = parts*P/(R*T)
    mgm3 = MW*n*1000  # mol toxin /m^3 * g/mol toxis * 1000 mg/g
    return mgm3

def mgm3_to_ppmv(mgm3, MW, T=298.15, P=101325.):
    r'''
    Converts a concentration in  mg/m^3 to units of ppmv. Used in
    industrial toxicology.

    .. math::
        ppmv = \frac{1000RT}{MW\cdot P} \cdot \frac{mg}{m^3}

    Parameters
    ----------
    mgm3 : float
        Concentration of a substance in an ideal gas mixture [mg/m^3]
    MW : float
        Molecular weight of the trace gas [g/mol]
    T : float, optional
        Temperature of the gas at which the ppmv is reported, [K]
    P : float, optional
        Pressure of the gas at which the ppmv is reported, [Pa]

    Returns
    -------
    ppmv : float
        Concentration of a component in a gas mixure [parts per million,
        volumetric]

    Notes
    -----
    The term P/(RT)/1000 converts to 0.040874 at STP. Its inverse is reported
    as 24.45 in [1]_.

    Examples
    --------
    >>> mgm3_to_ppmv(1.635, 40.0)
    1.0000233761164334

    References
    ----------
    .. [1] ACGIH. Industrial Ventilation: A Manual of Recommended Practice,
       23rd Edition. American Conference of Governmental and Industrial
       Hygenists, 2004.
    '''
    n = mgm3/MW/1000.
    parts = n*R*T/P
    ppm = parts/1E-6
    return ppm

### Data


NTP_codes = {1: 'Known', 2: 'Reasonably Anticipated'}
IARC_codes = {1: 'Carcinogenic to humans (1)',
              11: 'Probably carcinogenic to humans (2A)',  # 2A
              12: 'Possibly carcinogenic to humans (2B)',  # 2B
              3: 'Not classifiable as to its carcinogenicity to humans (3)',
              4: 'Probably not carcinogenic to humans (4)'}

folder = os_path_join(source_path, 'Safety')
register_df_source(folder, 'NFPA 497 2008.tsv')
register_df_source(folder, 'IS IEC 60079-20-1 2010.tsv')
register_df_source(folder, 'DIPPR T_flash Serat.csv')
register_df_source(folder, 'National Toxicology Program Carcinogens.tsv')
register_df_source(folder, 'IARC Carcinogen Database.tsv')
_safety_data_loaded = False


IEC = 'IEC 60079-20-1 (2010)'
NFPA = 'NFPA 497 (2008)'
SERAT = 'Serat DIPPR (2017)'

SUZUKI = 'Suzuki (1994)'
CROWLLOUVAR = 'Crowl and Louvar (2001)'

def _load_safety_data():
    global Ontario_exposure_limits_dict, NFPA_2008_data, IEC_2010_data
    global DIPPR_SERAT_data, NTP_data, IARC_data, Tflash_sources
    global Tautoignition_sources, LFL_sources, UFL_sources, _safety_data_loaded
    import json
    from io import open
    file = os_path_join(folder, 'Ontario Exposure Limits.json')
    with open(file, 'r') as stream:
        Ontario_exposure_limits_dict = json.load(stream)
    NFPA_2008_data = data_source('NFPA 497 2008.tsv')
    IEC_2010_data = data_source('IS IEC 60079-20-1 2010.tsv')
    DIPPR_SERAT_data = data_source('DIPPR T_flash Serat.csv')
    NTP_data = data_source('National Toxicology Program Carcinogens.tsv')
    IARC_data = data_source('IARC Carcinogen Database.tsv')
    Tflash_sources = {IEC: IEC_2010_data,
                      NFPA: NFPA_2008_data,
                      SERAT: DIPPR_SERAT_data}
    Tautoignition_sources = {IEC: IEC_2010_data,
                             NFPA: NFPA_2008_data}
    LFL_sources = Tautoignition_sources.copy()
    UFL_sources = Tautoignition_sources.copy()
    _safety_data_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('Ontario_exposure_limits_dict', 'NFPA_2008_data', 'IEC_2010_data',
                    'DIPPR_SERAT_data'):
            _load_safety_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else: # pragma: no cover
    if can_load_data:
        _load_safety_data()

# # Used to read Ontario Expore Limits data from original file (DO NOT DELETE!)
# Ontario_exposure_limits_dict = {}
# def str_to_ppm_mgm3(line, MW):  # pragma: no cover
#     if not line:
#         return None, None
#     if 'ppm' in line:
#         _ppm = float(line.split('ppm')[0])
#         try:
#             _mgm3 = ppmv_to_mgm3(_ppm, MW)
#         except:
#             _mgm3 = None
#     elif 'mg/m3' in line:
#         _mgm3 = float(line.split('mg/m3')[0])
#         try:
#             _ppm = mgm3_to_ppmv(_mgm3, MW)
#         except:
#             _ppm = None
#     if not _ppm and not _mgm3:
#         raise Exception('failure in function')
#     return (_ppm, _mgm3)

# with open(os.path.join(folder, 'Ontario Exposure Limits.tsv'), encoding='utf-8') as f:
#     '''Read in a dict of TWAs, STELs, and Ceiling Limits. The data source
#     is the Ontario Labor Website. They have obtained their data in part from
#     their own reviews, and also from ACGIH.
#     Warning: The lowest value is taken, when multiple units or different forms
#              of a compound are listed.
#     Note that each province has a different set of values, but these serve
#     as general values.
#     '''
#     next(f)
#     for line in f:
#         values = to_num(line.strip('\n').split('\t'))

#         if values[0]:
#             if type(values[6]) == str:
#                 MWs = [float(i) if i != '' else None for i in values[6].split(';')]
#             elif values[6] is None:
#                 MWs = [None]
#             elif type(values[6]) == float:
#                 MWs = [values[6]]
#             else:
#                 MWs = [None]

#             for i, CASRN in enumerate(values[0].split(';')):
#                 try:
#                     MWi = MWs[i]
#                 except IndexError:
#                     MWi = None


#                 _ppm_TWA, _mgm3_TWA = str_to_ppm_mgm3(values[2], MWi)
#                 _ppm_STEL, _mgm3_STEL = str_to_ppm_mgm3(values[3], MWi)
#                 _ppm_C, _mgm3_C = str_to_ppm_mgm3(values[4], MWi)

#                 if values[5] == 'Skin':
#                     _skin = True
#                 else:
#                     _skin = False


#                 Ontario_exposure_limits_dict[CASRN] = {"Name": values[1],  "TWA (ppm)": _ppm_TWA,
#                 "TWA (mg/m^3)": _mgm3_TWA, "STEL (ppm)": _ppm_STEL,
#                 "STEL (mg/m^3)": _mgm3_STEL, "Ceiling (ppm)": _ppm_C,
#                 "Ceiling (mg/m^3)": _mgm3_C, "Skin":_skin, "MW": MWi}


#TODO: Add CRC exposure limits. Note that functions should be used.
#_CRCExposureLimits = {}
#with open(os.path.join(folder,'CRC Exposure Limits.csv')) as f:
#    '''Read in a dict of TWAs and STELs. The data source
#    is the CRC Handbook.. They have obtained their data from
#    NIOSH, and OSHA Chemical Information Manual, and ACGIH.
#    '''
#    f.next()
#    for line in f:
#        values = to_num(line.strip('\n').split('\t'))
#        _ppm_TWA, _mgm3_TWA = str_to_ppm_mgm3(values[2], CASRN.strip())
#        _ppm_STEL, _mgm3_STEL = str_to_ppm_mgm3(values[3], CASRN.strip())
#        _CRCExposureLimits[CASRN] = {"Name": values[1],  "TWA (ppm)": _ppm_TWA,
#        "TWA (mg/m^3)": _mgm3_TWA, "STEL (ppm)": _ppm_STEL,
#        "STEL (mg/m^3)": _mgm3_STEL}
#del _ppm_TWA, _mgm3_TWA, _ppm_STEL, _mgm3_STEL, _ppm_C, _mgm3_C, status
#print Ontario_exposure_limits_dict['109-73-9']
##{'STEL (ppm)': None, 'Name': 'n-Butylamine [109-73-9]', 'Ceiling (mg/m^3)': 14.956408997955013, 'Ceiling (ppm)': 5.0, 'TWA (mg/m^3)': None, 'STEL (mg/m^3)': None, 'TWA (ppm)': None}
#print Ontario_exposure_limits_dict['34590-94-8']
#{'STEL (ppm)': 150.0, 'Name': '(2-Methoxymethylethoxy) propanol (DPGME) [34590-94-8]', 'Ceiling (mg/m^3)': None, 'Ceiling (ppm)': None, 'TWA (mg/m^3)': None, 'STEL (mg/m^3)': None, 'TWA (ppm)': 100.0}


# TODO: Add https://www.worksafebc.com/en/law-policy/occupational-health-safety/searchable-ohs-regulation/ohs-guidelines/guidelines-part-05#ExposureLimits
### OSHA exposure limit functions

ONTARIO = 'Ontario Limits'
TWA_all_methods = (ONTARIO,)
'''Tuple of method name keys. See the :obj:`TWA` for the actual references'''

STEL_all_methods = (ONTARIO,)
'''Tuple of method name keys. See the :obj:`STEL` for the actual references'''

Ceiling_all_methods = (ONTARIO,)
'''Tuple of method name keys. See the :obj:`Ceiling` for the actual references'''

Skin_all_methods = (ONTARIO,)
'''Tuple of method name keys. See the :obj:`Skin` for the actual references'''

def TWA_methods(CASRN):
    """Return all methods available to obtain the Time-Weighted Average exposure
    limits (TWA) for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain TWA with the given inputs.

    See Also
    --------
    TWA

    Examples
    --------
    >>> TWA_methods('71-43-2')
    ['Ontario Limits']
    """
    if not _safety_data_loaded: _load_safety_data()
    if CASRN in Ontario_exposure_limits_dict:
        data = Ontario_exposure_limits_dict[CASRN]
        if (data["TWA (ppm)"] or data["TWA (mg/m^3)"]): return [ONTARIO]
    return []

def TWA(CASRN, method=None):
    """Return the Time-Weighted Average exposure
    limits (TWA) for the desired chemical if it is available.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]
    method : str
        Name of method to use, [-]

    Returns
    -------
    TWA : float
        Time-Weighted Average exposure, [ppm or mg/m^3]
    units : str
        One of ppm or mg/m^3, [-]

    Notes
    -----
    The ppm value is preferentially returned if both are available. While they
    can be converted in specific cases, it is better to work with the specified
    units of the original source.

    Examples
    --------
    >>> TWA('98-00-0')
    (10.0, 'ppm')
    >>> TWA('1303-00-0')
    (5.0742430905659505e-05, 'ppm')
    """
    if not _safety_data_loaded: _load_safety_data()
    if not method or method == ONTARIO:
        if CASRN in Ontario_exposure_limits_dict:
            data = Ontario_exposure_limits_dict[CASRN]
            value = data["TWA (ppm)"]
            if value: return value, 'ppm'
            value = data["TWA (mg/m^3)"]
            if value: return value, 'mg/m^3'
    else:
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                         method, TWA_all_methods))

def STEL_methods(CASRN):
    """Return all methods available to obtain STEL for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain STEL with the given inputs.

    See Also
    --------
    STEL
    """
    if not _safety_data_loaded: _load_safety_data()
    if CASRN in Ontario_exposure_limits_dict:
        data = Ontario_exposure_limits_dict[CASRN]
        if (data["STEL (ppm)"] or data["STEL (mg/m^3)"]): return [ONTARIO]
    return []

def STEL(CASRN, method=None):
    """This function handles the retrieval of Short-term Exposure Limit (STEL)
    on worker exposure to dangerous chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]
    method : str
        Name of method to use, [-]

    Returns
    -------
    STEL : float
        Short-term Exposure Limit, [ppm or mg/m^3]
    units : str
        One of ppm or mg/m^3, [-]

    Notes
    -----
    The ppm value is preferentially returned if both are available. While they
    can be converted in specific cases, it is better to work with the specified
    units of the original source.

    Examples
    --------
    >>> STEL('67-64-1')
    (750.0, 'ppm')
    >>> STEL('7664-38-2')
    (0.7489774978301237, 'ppm')
    >>> STEL('55720-99-5')
    (2.0, 'mg/m^3')
    """
    if not _safety_data_loaded: _load_safety_data()
    if not method or method == ONTARIO:
        if CASRN in Ontario_exposure_limits_dict:
            data = Ontario_exposure_limits_dict[CASRN]
            value = data["STEL (ppm)"]
            if value: return value, 'ppm'
            value = data["STEL (mg/m^3)"]
            if value: return value, 'mg/m^3'
    else:
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                         method, TWA_all_methods))

def Ceiling_methods(CASRN):
    """Return all methods available to obtain Ceiling limits for the desired
    chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Ceiling limits with the given inputs.

    See Also
    --------
    Ceiling
    """
    if not _safety_data_loaded: _load_safety_data()
    if CASRN in Ontario_exposure_limits_dict:
        data = Ontario_exposure_limits_dict[CASRN]
        if (data["Ceiling (ppm)"] or data["Ceiling (mg/m^3)"]): return [ONTARIO]
    return []

def Ceiling(CASRN, method=None):
    """This function handles the retrieval of ceiling limits on worker exposure
    to dangerous chemicals. Ceiling limits are not to be exceeded at any time.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]
    method : str
        Name of method to use, [-]

    Returns
    -------
    Ceiling : float
        Ceiling Limit, [ppm or mg/m^3]
    units : str
        One of ppm or mg/m^3, [-]

    Examples
    --------
    >>> Ceiling('75-07-0')
    (25.0, 'ppm')
    >>> Ceiling('1395-21-7')
    (6e-05, 'mg/m^3')
    """
    if not _safety_data_loaded: _load_safety_data()
    if not method or method == ONTARIO:
        if CASRN in Ontario_exposure_limits_dict:
            data = Ontario_exposure_limits_dict[CASRN]
            value = data["Ceiling (ppm)"]
            if value: return value, 'ppm'
            value = data["Ceiling (mg/m^3)"]
            if value: return value, 'mg/m^3'
    else:
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                         method, list(Ontario_exposure_limits_dict)))

def Skin_methods(CASRN):
    """Return all methods available to obtain whether or not a chemical can be
    absorbed through the skin.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain whether or not a chemical can
        be absorbed through the skin.

    See Also
    --------
    Skin
    """
    return [ONTARIO] if CASRN in Ontario_exposure_limits_dict else []

def Skin(CASRN, method=None):
    """This function handles the retrieval of whether or not a chemical can be
    absorbed through the skin, relevant to chemical safety calculations.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]
    method : str
        Name of method to use, [-]

    Returns
    -------
    skin : bool
        Whether or not the substance is absorbed through human skin, [-]

    Examples
    --------
    >>> Skin('108-94-1')
    True
    >>> Skin('1395-21-7')
    False
    """
    if not _safety_data_loaded: _load_safety_data()
    if not method or method == ONTARIO:
        if CASRN in Ontario_exposure_limits_dict:
            return Ontario_exposure_limits_dict[CASRN]["Skin"]
    else:
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                         method, TWA_all_methods))

### Carcinogen functions

IARC = 'International Agency for Research on Cancer'
NTP = 'National Toxicology Program 13th Report on Carcinogens'
UNLISTED = 'Unlisted'
COMBINED = 'Combined'

Carcinogen_all_methods = (IARC, NTP)
'''Tuple of method name keys. See the :obj:`Carcinogen` for the actual references'''

def Carcinogen_methods(CASRN):
    """Return all methods available to obtain Carcinogen listings for the
    desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Carcinogen listings with the given inputs.

    See Also
    --------
    Carcinogen
    """
    return list(Carcinogen_all_methods)

def Carcinogen(CASRN, method=None):
    r'''Looks up if a chemical is listed as a carcinogen or not according to
    either a specifc method or with all methods.
    Returns either the status as a string for a specified method, or the
    status of the chemical in all available data sources, in the format
    {source: status}.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    status : str or dict
        Carcinogen status information [-].

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `Carcinogen_all_methods`.

    Notes
    -----
    Supported methods are:
        * **IARC**: International Agency for Research on Cancer, [1]_. As
          extracted with a last update of  February 22, 2016. Has listing
          information of 863 chemicals with CAS numbers. Chemicals without
          CAS numbers not included here. If two listings for the same CAS
          were available, the harshest rating was used. If two
          listings were available published at different times, the latest
          value was used. All else equal, the most pessimistic value was used.
        * **NTP**: National Toxicology Program, [2]_. Has data on 228
          chemicals.

    Examples
    --------
    >>> Carcinogen('61-82-5')
    {'International Agency for Research on Cancer': 'Not classifiable as to its carcinogenicity to humans (3)', 'National Toxicology Program 13th Report on Carcinogens': 'Reasonably Anticipated'}

    References
    ----------
    .. [1] International Agency for Research on Cancer. Agents Classified by
       the IARC Monographs, Volumes 1-115. Lyon, France: IARC; 2020 Available
       from: http://monographs.iarc.fr/ENG/Classification/
    .. [2] NTP (National Toxicology Program). 2016. Report on Carcinogens,
       Fourteenth Edition. Research Triangle Park, NC: U.S. Department of
       Health and Human Services, Public Health Service.
       https://ntp.niehs.nih.gov/whatwestudy/assessments/cancer/roc/index.html
    '''
    if not _safety_data_loaded: _load_safety_data()
    if not method:
        return {
            IARC: IARC_codes[IARC_data.at[CASRN, 'group']] if CASRN in IARC_data.index else UNLISTED,
            NTP: NTP_codes[NTP_data.at[CASRN, 'Listing']] if CASRN in NTP_data.index else UNLISTED
        }
    if method == IARC:
        if CASRN in IARC_data.index:
            return IARC_codes[IARC_data.at[CASRN, 'group']]
    elif method == NTP:
        if CASRN in NTP_data.index:
            return NTP_codes[NTP_data.at[CASRN, 'Listing']]
    else:
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                       method, Carcinogen_all_methods))
    return UNLISTED


### Fire-related functions


T_flash_all_methods = (IEC, NFPA, SERAT)
'''Tuple of method name keys. See the :obj:`T_flash` for the actual references'''

def T_flash_methods(CASRN):
    """Return all methods available to obtain T_flash for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain T_flash with the given inputs.

    See Also
    --------
    T_flash
    """
    if not _safety_data_loaded: _load_safety_data()
    return list_available_methods_from_df_dict(Tflash_sources, CASRN, 'T_flash')

def T_flash(CASRN, method=None):
    r'''
    This function handles the retrieval or calculation of a chemical's
    flash point. Lookup is based on CASRNs. No predictive methods are currently
    implemented. Will automatically select a data source to use if no method
    is provided; returns None if the data is not available.

    Examples
    --------
    >>> T_flash(CASRN='64-17-5')
    285.15

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    T_flash : float
        Flash point of the chemical, [K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `T_flash_all_methods`,

    Notes
    -----
    Preferred source is 'IEC 60079-20-1 (2010)' [1]_, with the secondary source
    'NFPA 497 (2008)' [2]_ having very similar data. A third source
    'Serat DIPPR (2017)' [3]_ provides third hand experimental but evaluated
    data from the DIPPR database, version unspecified, for 870 compounds.

    The predicted values from the DIPPR databank are also available in the
    supporting material in [3]_, but are not included.

    See Also
    --------
    T_flash_methods

    References
    ----------
    .. [1] IEC. "IEC 60079-20-1:2010 Explosive atmospheres - Part 20-1:
       Material characteristics for gas and vapour classification - Test
       methods and data." https://webstore.iec.ch/publication/635. See also
       https://law.resource.org/pub/in/bis/S05/is.iec.60079.20.1.2010.pdf
    .. [2] National Fire Protection Association. NFPA 497: Recommended
       Practice for the Classification of Flammable Liquids, Gases, or Vapors
       and of Hazardous. NFPA, 2008.
    .. [3] Serat, Fatima Zohra, Ali Mustapha Benkouider, Ahmed Yahiaoui, and
       Farid Bagui. "Nonlinear Group Contribution Model for the Prediction of
       Flash Points Using Normal Boiling Points." Fluid Phase Equilibria 449
       (October 15, 2017): 52-59. doi:10.1016/j.fluid.2017.06.008.

    '''
    if not _safety_data_loaded: _load_safety_data()
    if method:
        return retrieve_from_df_dict(Tflash_sources, CASRN, 'T_flash', method)
    else:
        return retrieve_any_from_df_dict(Tflash_sources, CASRN, 'T_flash')


T_autoignition_all_methods = (IEC, NFPA)
'''Tuple of method name keys. See the :obj:`T_autoignition` for the actual references'''

def T_autoignition_methods(CASRN):
    """Return all methods available to obtain T_autoignition for the desired
    chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain T_autoignition with the given inputs.

    See Also
    --------
    T_autoignition
    """
    if not _safety_data_loaded: _load_safety_data()
    return list_available_methods_from_df_dict(Tautoignition_sources, CASRN, 'T_autoignition')

def T_autoignition(CASRN, method=None):
    r'''
    This function handles the retrieval or calculation of a chemical's
    autoifnition temperature. Lookup is based on CASRNs. No predictive methods
    are currently implemented. Will automatically select a data source to use
    if no Method is provided; returns None if the data is not available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Tautoignition : float
        Autoignition point of the chemical, [K].

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `T_autoignition_all_methods`.

    Examples
    --------
    >>> T_autoignition(CASRN='71-43-2')
    771.15

    Notes
    -----
    Preferred source is 'IEC 60079-20-1 (2010)' [1]_, with the secondary source
    'NFPA 497 (2008)' [2]_ having very similar data.

    See Also
    --------
    T_autoignition_methods

    References
    ----------
    .. [1] IEC. “IEC 60079-20-1:2010 Explosive atmospheres - Part 20-1:
       Material characteristics for gas and vapour classification - Test
       methods and data.” https://webstore.iec.ch/publication/635. See also
       https://law.resource.org/pub/in/bis/S05/is.iec.60079.20.1.2010.pdf
    .. [2] National Fire Protection Association. NFPA 497: Recommended
       Practice for the Classification of Flammable Liquids, Gases, or Vapors
       and of Hazardous. NFPA, 2008.
    '''
    if not _safety_data_loaded: _load_safety_data()
    if method:
        return retrieve_from_df_dict(Tautoignition_sources, CASRN, 'T_autoignition', method)
    else:
        return retrieve_any_from_df_dict(Tautoignition_sources, CASRN, 'T_autoignition')



LFL_all_methods = (IEC, NFPA, SUZUKI, CROWLLOUVAR)
'''Tuple of method name keys. See the :obj:`LFL` for the actual references'''

def LFL_methods(Hc=None, atoms=None, CASRN=''):
    """Return all methods available to obtain LFL for the desired chemical.

    Parameters
    ----------
    Hc : float, optional
        Heat of combustion of gas [J/mol].
    atoms : dict, optional
        Dictionary of atoms and atom counts.
    CASRN : str, optional
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain LFL with the given inputs.

    See Also
    --------
    LFL

    Examples
    --------
    Methane

    >>> LFL_methods(Hc=-890590.0, atoms={'C': 1, 'H': 4}, CASRN='74-82-8')
    ['IEC 60079-20-1 (2010)', 'NFPA 497 (2008)', 'Suzuki (1994)', 'Crowl and Louvar (2001)']
    """
    if not _safety_data_loaded: _load_safety_data()
    methods = list_available_methods_from_df_dict(LFL_sources, CASRN, 'LFL')
    if Hc is not None:
        methods.append(SUZUKI)
    if atoms is not None:
        methods.append(CROWLLOUVAR)
    return methods

def LFL(Hc=None, atoms=None, CASRN='', method=None):
    r'''This function handles the retrieval or calculation of a chemical's
    Lower Flammability Limit. Lookup is based on CASRNs. Will automatically
    select a data source to use if no Method is provided; returns None if the
    data is not available.

    Parameters
    ----------
    Hc : float, optional
        Heat of combustion of gas [J/mol].
    atoms : dict, optional
        Dictionary of atoms and atom counts.
    CASRN : str, optional
        CASRN, [-]

    Returns
    -------
    LFL : float
        Lower flammability limit of the gas in an atmosphere at STP, [mole fraction].

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `LFL_all_methods`.

    Examples
    --------
    >>> LFL(CASRN='71-43-2')
    0.012
    >>> LFL(Hc=-890590.0, atoms={'C': 1, 'H': 4}, CASRN='74-82-8')
    0.044

    Notes
    -----
    Preferred source is 'IEC 60079-20-1 (2010)' [1]_, with the secondary source
    'NFPA 497 (2008)' [2]_ having very similar data. If the heat of combustion
    is provided, the estimation method :obj:`Suzuki_LFL` can be used. If the atoms
    of the molecule are available, the method :obj:`Crowl_Louvar_LFL` can be used.

    References
    ----------
    .. [1] IEC. “IEC 60079-20-1:2010 Explosive atmospheres - Part 20-1:
       Material characteristics for gas and vapour classification - Test
       methods and data.” https://webstore.iec.ch/publication/635. See also
       https://law.resource.org/pub/in/bis/S05/is.iec.60079.20.1.2010.pdf
    .. [2] National Fire Protection Association. NFPA 497: Recommended
       Practice for the Classification of Flammable Liquids, Gases, or Vapors
       and of Hazardous. NFPA, 2008.

    '''
    if not _safety_data_loaded: _load_safety_data()
    if not method:
        LFL = retrieve_any_from_df_dict(LFL_sources, CASRN, 'LFL')
        if not LFL:
            if Hc == 0.0:
                return None
            if Hc: LFL = Suzuki_LFL(Hc)
            elif atoms: LFL = Crowl_Louvar_LFL(atoms)
        return LFL
    elif method == SUZUKI:
        return Suzuki_LFL(Hc)
    elif method == CROWLLOUVAR:
        return Crowl_Louvar_LFL(atoms)
    else:
        return retrieve_from_df_dict(LFL_sources, CASRN, 'LFL', method)

UFL_all_methods = (IEC, NFPA, SUZUKI, CROWLLOUVAR)
'''Tuple of method name keys. See the :obj:`UFL` for the actual references'''

def UFL_methods(Hc=None, atoms=None, CASRN=''):
    """Return all methods available to obtain UFL for the desired chemical.

    Parameters
    ----------
    Hc : float, optional
        Heat of combustion of gas [J/mol].
    atoms : dict, optional
        Dictionary of atoms and atom counts.
    CASRN : str, optional
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain UFL with the given inputs.

    See Also
    --------
    UFL

    Examples
    --------
    Methane

    >>> UFL_methods(Hc=-890590.0, atoms={'C': 1, 'H': 4}, CASRN='74-82-8')
    ['IEC 60079-20-1 (2010)', 'NFPA 497 (2008)', 'Suzuki (1994)', 'Crowl and Louvar (2001)']
    """
    if not _safety_data_loaded: _load_safety_data()
    methods = list_available_methods_from_df_dict(UFL_sources, CASRN, 'UFL')
    if Hc is not None:
        methods.append(SUZUKI)
    if atoms is not None:
        methods.append(CROWLLOUVAR)
    return methods

def UFL(Hc=None, atoms=None, CASRN='', method=None):
    r'''This function handles the retrieval or calculation of a chemical's
    Upper Flammability Limit. Lookup is based on CASRNs. Two predictive methods
    are currently implemented. Will automatically select a data source to use
    if no Method is provided; returns None if the data is not available.

    Examples
    --------
    >>> UFL(CASRN='71-43-2')
    0.086

    Methane

    >>> UFL(Hc=-890590.0, atoms={'C': 1, 'H': 4}, CASRN='74-82-8')
    0.17

    Parameters
    ----------
    Hc : float, optional
        Heat of combustion of gas [J/mol]
    atoms : dict, optional
        Dictionary of atoms and atom counts
    CASRN : str, optional
        CASRN [-]

    Returns
    -------
    UFL : float
        Upper flammability limit of the gas in an atmosphere at STP, [mole fraction]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `UFL_all_methods`.

    Notes
    -----
    Preferred source is 'IEC 60079-20-1 (2010)' [1]_, with the secondary source
    'NFPA 497 (2008)' [2]_ having very similar data. If the heat of combustion
    is provided, the estimation method :obj:`Suzuki_UFL` can be used. If the atoms
    of the molecule are available, the method :obj:`Crowl_Louvar_UFL` can be used.

    References
    ----------
    .. [1] IEC. “IEC 60079-20-1:2010 Explosive atmospheres - Part 20-1:
       Material characteristics for gas and vapour classification - Test
       methods and data.” https://webstore.iec.ch/publication/635. See also
       https://law.resource.org/pub/in/bis/S05/is.iec.60079.20.1.2010.pdf
    .. [2] National Fire Protection Association. NFPA 497: Recommended
       Practice for the Classification of Flammable Liquids, Gases, or Vapors
       and of Hazardous. NFPA, 2008.

    '''
    if not _safety_data_loaded: _load_safety_data()
    if not method:
        UFL = retrieve_any_from_df_dict(UFL_sources, CASRN, 'UFL')
        if not UFL:
            if Hc is not None: UFL = Suzuki_UFL(Hc)
            elif atoms is not None: UFL = Crowl_Louvar_UFL(atoms)
        return UFL
    elif method == SUZUKI:
        return Suzuki_UFL(Hc)
    elif method == CROWLLOUVAR:
        return Crowl_Louvar_UFL(atoms)
    else:
        return retrieve_from_df_dict(UFL_sources, CASRN, 'UFL', method)

ISO_10156_2017_Kks = {
    '7782-44-7': 1.0, # O2,
    '7727-37-9': 1.0, # N2
    '124-38-9': 1.5, # CO2
    '7440-59-7': 0.9, # He
    '7440-37-1': 0.55, # Ar
    '7440-01-9': 0.7, # Ne
    '7439-90-9': 0.5, # Kr
    '7440-63-3': 0.5, # Xe
    '7446-09-5': 1.5, # SO2
    '2551-62-4': 4.0, # SF6
    '75-73-0': 2.0, # CF4
    '76-19-7': 1.5, # C3F8
    '354-33-6': 3.5, # C2HF5
}

def LFL_ISO_10156_2017(zs, LFLs, CASs):
    r'''Calculate the lower flammability limit of a mixture of combustible gases
    and inert gases according to ISO 10156 (2017) [1]_.

    .. math::
        \text{LFL} = \frac{1}{\sum_{i=1}^{n_{combustible}}\frac{A_i}{\text{LFL}_i'}}

    .. math::
        \text{LFL}_i' = \frac{1 - \text{LFL}_m' - (1 - K)
        \frac{\sum_j^{n_{inert}} B_j}{\sum_j^{n_{combustible}} A_j} \text{LFL}_m'
        }
        {100 - \text{LFL}_m'}\text{LFL}_i

    .. math::
        K = \sum_i^{n_{inert}} z_i K_k

    The `B` sum is the total mole fraction of all inert gas compounds;
    and the `A` sum is the total mole fraction of all combustible compounds.
    :math:`K_k` are the looked up inert gas coefficients.
    :math:`\text{LFL}_m'` is calculated as the Le Chatelier's lower
    flammability limit if there were no inert gases in the mixture.

    Parameters
    ----------
    zs : list[float]
        Mole fractions of all components in a gas including inerts, [-]
    LFLs : list[float]
        Lower or upper flammability limits for each flammable component in a
        gas, [-]
    CASs : list[str]
        CAS numbers of each compound; required to look up inert gas factors, [-]

    Returns
    -------
    LFL : float
        Lower or flammability limit of a gas mixture, [-]

    Notes
    -----
    Inert gas parameters are available for O2, N2, CO2, He, Ar, Ne, Kr, Xe,
    SO2, SF6, CF4, C3F8, and C2HF5.

    Examples
    --------
    All the sample problems from [1]_ have been implemented as tests.

    >>> zs = [.15, .15, .3, .35+.05*.79, .05*.21]
    >>> LFLs = [.04, .044, None, None, None]
    >>> CASs = ['1333-74-0', '74-82-8', '124-38-9', '7727-37-9', '7782-44-7']
    >>> LFL_ISO_10156_2017(zs, LFLs, CASs)
    0.1427372274

    References
    ----------
    .. [1] Standardization, International Organization for. ISO 10156: 2017 :
       Gas Cylinders - Gases and Gas Mixtures - Determination of Fire Potential
       and Oxidizing Ability for the Selection of Cylinder Valve Outlets, 2017.
    '''
    N = len(zs)
    has_inerts = False
    for CAS in CASs:
        if CAS in ISO_10156_2017_Kks:
            has_inerts = True
            break

    if has_inerts:
        combustible_idxs = []
        combustible_zs = []
        for i, CAS in enumerate(CASs):
            # Combustible
            if CAS not in ISO_10156_2017_Kks:
                combustible_idxs.append(i)
                combustible_zs.append(zs[i])

        combustible_zs_norm = normalize(combustible_zs)
        Lm_prime = 0.0
        for idx, zi in zip(combustible_idxs, combustible_zs_norm):
            Lm_prime += zi/LFLs[idx]
        Lm_prime = 1.0/Lm_prime

        combustible_frac = sum(combustible_zs)
        inert_frac = 1.0 - combustible_frac
        K = 0.0
        for i, CAS in enumerate(CASs):
            if CAS in ISO_10156_2017_Kks:
                K += ISO_10156_2017_Kks[CAS]*zs[i]
        K /= inert_frac

        factor = ((1.0 - Lm_prime - (1.0 - K)*inert_frac/combustible_frac*Lm_prime)
                  /(1.0 - Lm_prime))
        Lm = 0.0
        for idx, zi in zip(combustible_idxs, combustible_zs):
            Lm += zi/(factor*LFLs[idx])
        Lm = 1.0/Lm
        return Lm

    tot = 0.0
    for i in range(len(zs)):
        tot += zs[i]/LFLs[i]
    return 1.0/tot

def fire_mixing(ys, FLs):
    '''Le Chatelier's mixing rule for lower and upper flammability limits of
    mixtures of gases.

    Parameters
    ----------
    ys : list[float]
        Normalized mole fractions of all flammable components in a gas, [-]
    FLs : list[float]
        Lower or upper flammability limits for each flammable component in a
        gas, [-]

    Returns
    -------
    FL : float
        Lower or upper flammability limit of a gas, [-]

    Notes
    -----
    This equation has a higher accuracy for lower flammability limits than
    upper flammability limits. Some sources recommend not using it for
    upper flammability limits.

    Examples
    --------
    Sample problems from [1]_ for the lower and upper flammability limit.

    >>> fire_mixing(ys=normalize([0.0024, 0.0061, 0.0015]), FLs=[.012, .053, .031])
    0.02751172136637642

    >>> fire_mixing(ys=normalize([0.0024, 0.0061, 0.0015]), FLs=[.075, .15, .32])
    0.12927551844869378

    References
    ----------
    .. [1] Crowl, Daniel A., and Joseph F. Louvar. Chemical Process Safety:
       Fundamentals with Applications. 2E. Upper Saddle River, N.J: Prentice
       Hall, 2001.
    '''
    tot = 0.0
    for i in range(len(ys)):
        tot += ys[i]/FLs[i]
    return 1.0/tot

def Suzuki_LFL(Hc):
    r'''Calculates lower flammability limit, using the Suzuki [1]_ correlation.
    Uses heat of combustion only.

    The lower flammability limit of a gas is air is:

    .. math::
        \text{LFL} = \frac{-3.42}{\Delta H_c^{\circ}} + 0.569

    .. math::
        \Delta H_c^{\circ} + 0.0538\Delta H_c^{\circ 2} + 1.80

    Parameters
    ----------
    Hc : float
        Heat of combustion of gas [J/mol]

    Returns
    -------
    LFL : float
        Lower flammability limit, mole fraction [-]

    Notes
    -----
    Fit performed with 112 compounds, r^2 was 0.977.
    LFL in percent volume in air. Hc is at standard conditions, in MJ/mol.
    11 compounds left out as they were outliers.
    Equation does not apply for molecules with halogen atoms, only hydrocarbons
    with oxygen or nitrogen or sulfur.
    No sample calculation provided with the article. However, the equation is
    straightforward.
    Limits of equations's validity are -6135596 J where it predicts a
    LFL of 0, and -48322129 J where it predicts a LFL of 1.

    Examples
    --------
    Pentane, 1.5 % LFL in literature

    >>> Suzuki_LFL(-3536600)
    0.014276107095811815

    References
    ----------
    .. [1] Suzuki, Takahiro. "Note: Empirical Relationship between Lower
       Flammability Limits and Standard Enthalpies of Combustion of Organic
       Compounds." Fire and Materials 18, no. 5 (September 1, 1994): 333-36.
       doi:10.1002/fam.810180509.
    '''
    Hc = Hc/1E6
    LFL = -3.42/Hc + 0.569*Hc + 0.0538*Hc*Hc + 1.80
    return LFL/100.


def Suzuki_UFL(Hc):
    r'''Calculates upper flammability limit, using the Suzuki [1]_ correlation.
    Uses heat of combustion only.
    The upper flammability limit of a gas is air is:

    .. math::
        \text{UFL} = 6.3\Delta H_c^\circ + 0.567\Delta H_c^{\circ 2} + 23.5

    Parameters
    ----------
    Hc : float
        Heat of combustion of gas [J/mol]

    Returns
    -------
    UFL : float
        Upper flammability limit, mole fraction

    Notes
    -----
    UFL in percent volume in air according to original equation.
    Hc is at standard conditions in the equation, in units of MJ/mol.
    AAPD = 1.2% for 95 compounds used in fit.
    Somewhat better results than the High and Danner method.
    4.9% < UFL < 23.0%
    -890.3 kJ/mol < dHc < -6380 kJ/mol
    r^2 = 0.989
    Sample calculations provided for all chemicals, both this method and
    High and Danner. Examples are from the article.
    Predicts a UFL of 1 at 7320190 J and a UFL of 0 at -5554160 J.

    Examples
    --------
    Pentane, literature 7.8% UFL

    >>> Suzuki_UFL(-3536600)
    0.0831119493052

    References
    ----------
    .. [1] Suzuki, Takahiro, and Kozo Koide. "Short Communication: Correlation
       between Upper Flammability Limits and Thermochemical Properties of
       Organic Compounds." Fire and Materials 18, no. 6 (November 1, 1994):
       393-97. doi:10.1002/fam.810180608.
    '''
    Hc = Hc/1E6
    UFL = 6.3*Hc + 0.567*Hc*Hc + 23.5
    return UFL/100.


def Crowl_Louvar_LFL(atoms):
    r'''Calculates lower flammability limit, using the Crowl-Louvar [1]_
    correlation. Uses molecular formula only.
    The lower flammability limit of a gas is air is:

    .. math::
        C_mH_xO_y + zO_2 \to mCO_2 + \frac{x}{2}H_2O

    .. math::
        \text{LFL} = \frac{0.55}{4.76m + 1.19x - 2.38y + 1}

    Parameters
    ----------
    atoms : dict
        Dictionary of atoms and atom counts

    Returns
    -------
    LFL : float
        Lower flammability limit, mole fraction

    Notes
    -----
    Coefficient of 0.55 taken from [2]_

    Examples
    --------
    Hexane, example from [1]_, lit. 1.2 %

    >>> Crowl_Louvar_LFL({'H': 14, 'C': 6})
    0.011899610558199915

    References
    ----------
    .. [1] Crowl, Daniel A., and Joseph F. Louvar. Chemical Process Safety:
       Fundamentals with Applications. 2E. Upper Saddle River, N.J:
       Prentice Hall, 2001.
    .. [2] Jones, G. W. "Inflammation Limits and Their Practical Application
       in Hazardous Industrial Operations." Chemical Reviews 22, no. 1
       (February 1, 1938): 1-26. doi:10.1021/cr60071a001
    '''
    nC, nH, nO = 0, 0, 0
    if 'C' in atoms and atoms['C']:
        nC = atoms['C']
    else:
        return None
    if 'H' in atoms:
        nH = atoms['H']
    if 'O' in atoms:
        nO = atoms['O']
    return 0.55/(4.76*nC + 1.19*nH - 2.38*nO + 1.)


def Crowl_Louvar_UFL(atoms):
    r'''Calculates upper flammability limit, using the Crowl-Louvar [1]_
    correlation. Uses molecular formula only.
    The upper flammability limit of a gas is air is:

    .. math::
        C_mH_xO_y + zO_2 \to mCO_2 + \frac{x}{2}H_2O

    .. math::
        \text{UFL} = \frac{3.5}{4.76m + 1.19x - 2.38y + 1}

    Parameters
    ----------
    atoms : dict
        Dictionary of atoms and atom counts

    Returns
    -------
    UFL : float
        Upper flammability limit, mole fraction

    Notes
    -----
    Coefficient of 3.5 taken from [2]_

    Examples
    --------
    Hexane, example from [1]_, lit. 7.5 %

    >>> Crowl_Louvar_UFL({'H': 14, 'C': 6})
    0.07572479446127219

    References
    ----------
    .. [1] Crowl, Daniel A., and Joseph F. Louvar. Chemical Process Safety:
       Fundamentals with Applications. 2E. Upper Saddle River, N.J:
       Prentice Hall, 2001.
    .. [2] Jones, G. W. "Inflammation Limits and Their Practical Application
       in Hazardous Industrial Operations." Chemical Reviews 22, no. 1
       (February 1, 1938): 1-26. doi:10.1021/cr60071a001
    '''
    nC, nH, nO = 0, 0, 0
    if 'C' in atoms and atoms['C']:
        nC = atoms['C']
    else:
        return None
    if 'H' in atoms:
        nH = atoms['H']
    if 'O' in atoms:
        nO = atoms['O']
    return 3.5/(4.76*nC + 1.19*nH - 2.38*nO + 1.)

def NFPA_30_classification(T_flash, Tb=None, Psat_100F=None):
    r'''Classify a chemical's flammability/combustibility according
    to the NFPA 30 standard Flammable and Combustible Liquids Code.

    Class IA: Flash Point < 73°F; Boiling Point < 100°F
    Class IB: Flash Point < 73°F; 100°F <= Boiling Point
    Class IC: 73°F <= Flash Point < 100°F
    Class II: 100°F <= Flash Point < 140°F
    Class IIIA: 140°F <= Flash Point < 200°F
    Class IIIB: 200°F <= Flash Point

    Class I liquids are designated as flammable; class II and II
    liquids are designated as combustible.

    Parameters
    ----------
    T_flash : float
        Flash point (closed‐cup method, adjusted for sea level), [K]
    Tb : float, optional
        Normal boiling point (needed to classify IA and IB liquids), [K]
    Psat_100F : float, optional
        Vapor pressure at 100°F (needed to classify IA and IB liquids), [K]

    Returns
    -------
    classification : str
        One of 'IA', 'IB', 'IC', 'II', 'IIIA', 'IIIB', [-]

    Notes
    -----
    Only one of `Tb` or `Psat_100F` is needed.

    Class 'IA' also includes unstable liquids.

    Examples
    --------
    Ethylene oxide

    >>> NFPA_30_classification(253.15, 283.55)
    'IA'

    Butyl alcohol

    >>> NFPA_30_classification(308.15)
    'IC'

    References
    ----------
    .. [1] NFPA (National Fire Prevention Association). NFPA 30: Flammable and
       Combustible Liquids Code, 2008. National Fire Protection Association
       (NFPA), 2007.
    '''
    F_100 = 310.92777777777777
    F_73 = 295.92777777777775
    F_140 = 333.15
    F_200 = 366.4833333333333
    Tb_above_100F = False # For numba only
    if Tb is not None:
        Tb_above_100F = Tb >= F_100
    elif Psat_100F is not None:
        Tb_above_100F = Psat_100F < 101325.0
    elif T_flash < F_73:
        raise ValueError("Tb or Psat_100F is required to classify the provided inputs")

    if T_flash < F_100:
        if T_flash < F_73 and not Tb_above_100F:
            # Also unstable flammable liquids
            # ethylene oxide, methyl chloride, pentane should go here
            return 'IA'
        elif T_flash < F_73 and Tb_above_100F:
            # acetone, benzene, ethyl alcohol, and isopropyl alcohol.
            return 'IB'
        elif F_73 <= T_flash < F_100:
            # butyl alcohol, diethyl glycol, styrene, and turpentine.
            return 'IC'
    if F_100 <= T_flash < F_140:
        return 'II'
    if F_140 <= T_flash < F_200:
        return 'IIIA'
    if F_200 <= T_flash:
        return 'IIIB'
