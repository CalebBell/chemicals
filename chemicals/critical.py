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

This module contains lookup functions for critical temperature, 
critical pressure, critical volume, and critical compressibility factors.
It also includes a few relationships between the critical properties, and a 
variety of critical mixture property estimation routines.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Critical Temperature
--------------------
.. autofunction:: chemicals.critical.Tc
.. autofunction:: chemicals.critical.Tc_methods
.. autodata:: chemicals.critical.Tc_all_methods

Critical Pressure
-----------------
.. autofunction:: chemicals.critical.Pc
.. autofunction:: chemicals.critical.Pc_methods
.. autodata:: chemicals.critical.Pc_all_methods

Critical Volume
---------------
.. autofunction:: chemicals.critical.Vc
.. autofunction:: chemicals.critical.Vc_methods
.. autodata:: chemicals.critical.Vc_all_methods
.. autofunction:: chemicals.critical.Mersmann_Kind_predictor

Critical Compressibility Factor
-------------------------------
.. autofunction:: chemicals.critical.Zc
.. autofunction:: chemicals.critical.Zc_methods
.. autodata:: chemicals.critical.Zc_all_methods

Critical Property Relationships
-------------------------------
.. autofunction:: chemicals.critical.critical_surface
.. autofunction:: chemicals.critical.critical_surface_methods
.. autodata:: chemicals.critical.critical_surface_all_methods
.. autofunction:: chemicals.critical.third_property
.. autofunction:: chemicals.critical.Ihmels
.. autofunction:: chemicals.critical.Meissner
.. autofunction:: chemicals.critical.Grigoras
.. autofunction:: chemicals.critical.Hekayati_Raeissi

Critical Temperature of Mixtures
--------------------------------
.. autofunction:: chemicals.critical.Li
.. autofunction:: chemicals.critical.Chueh_Prausnitz_Tc
.. autofunction:: chemicals.critical.Grieves_Thodos
.. autofunction:: chemicals.critical.modified_Wilson_Tc

Critical Volume of Mixtures
---------------------------
.. autofunction:: chemicals.critical.Chueh_Prausnitz_Vc
.. autofunction:: chemicals.critical.modified_Wilson_Vc
"""

__all__ = ['Tc', 'Pc', 'Vc', 'Zc',
           'Mersmann_Kind_predictor', 
           'third_property', 
           'critical_surface', 
           'Ihmels', 'Meissner', 'Grigoras', 'Hekayati_Raeissi', 'Li', 
           'Chueh_Prausnitz_Tc', 'Grieves_Thodos',
           'modified_Wilson_Tc', 'Chueh_Prausnitz_Vc', 
           'modified_Wilson_Vc',
           'Tc_methods', 'Pc_methods', 
           'Vc_methods', 'Zc_methods', 
           'critical_surface_methods',
           'Tc_all_methods', 'Pc_all_methods', 
           'Vc_all_methods', 'Zc_all_methods', 
           'critical_surface_all_methods']

__numba_additional_funcs__ = ['_assert_two_critical_properties_provided']

import os
from fluids.constants import R, R_inv, N_A
from chemicals.utils import log, PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

folder = os_path_join(source_path, 'Critical Properties')
IUPAC = 'IUPAC'
MATTHEWS = 'MATTHEWS'
CRC = 'CRC'
PSRK = 'PSRK'
PD = 'PD'
YAWS = 'YAWS'

### Register data sources and lazy load them

def _add_Zc_to_df(df):
    # Some files don't have the `Zc` column; this adds it
    # TODO: Think about adding these to the files
    import pandas as pd
    if 'Zc' in df: return
    df['Zc'] = pd.Series(df['Pc']*df['Vc']*R_inv/df['Tc'], index=df.index)

# IUPAC Organic data series
# TODO: 12E of this data http://pubsdc3.acs.org/doi/10.1021/acs.jced.5b00571
register_df_source(folder, 'IUPACOrganicCriticalProps.tsv')

# CRC Handbook from TRC Organic data section (only in 2015)
# No Inorganic table was taken, although it is already present;
# data almost all from IUPAC
register_df_source(folder, 'CRCCriticalOrganics.tsv', postload=_add_Zc_to_df)
register_df_source(folder, 'Mathews1972InorganicCriticalProps.tsv')
register_df_source(folder, 'Appendix to PSRK Revision 4.tsv', postload=_add_Zc_to_df)
register_df_source(folder, 'PassutDanner1973.tsv')
register_df_source(folder, 'Yaws Collection.tsv', postload=_add_Zc_to_df)
_critical_data_loaded = False
def _load_critical_data():
    global critical_data_IUPAC, critical_data_Matthews, critical_data_CRC
    global critical_data_PSRKR4, critical_data_Yaws, critical_data_PassutDanner
    global Tc_sources, Pc_sources, Vc_sources, Zc_sources, omega_sources
    global _critical_data_loaded
    critical_data_IUPAC = data_source('IUPACOrganicCriticalProps.tsv')
    critical_data_Matthews = data_source('Mathews1972InorganicCriticalProps.tsv')
    critical_data_CRC = data_source('CRCCriticalOrganics.tsv')
    critical_data_PSRKR4 = data_source('Appendix to PSRK Revision 4.tsv')
    critical_data_Yaws = data_source('Yaws Collection.tsv')
    critical_data_PassutDanner = data_source('PassutDanner1973.tsv')
    _critical_data_loaded = True
    Tc_sources = {
        IUPAC: critical_data_IUPAC,
        MATTHEWS: critical_data_Matthews,
        CRC: critical_data_CRC,
        PSRK: critical_data_PSRKR4,
        PD: critical_data_PassutDanner,
        YAWS: critical_data_Yaws
    }
    
    # Create copies just incase new dfs need to be added later
    Pc_sources = Tc_sources.copy()
    Vc_sources = Tc_sources.copy()
    
    # The Passut Danner tsv file doesn't have Vc, so its not included
    del Vc_sources['PD']
    Zc_sources = Vc_sources.copy()

    omega_sources = {
        PSRK: critical_data_PSRKR4,
        PD: critical_data_PassutDanner,
        YAWS: critical_data_Yaws
    }

if PY37:
    def __getattr__(name):
        if name in ('critical_data_IUPAC', 'critical_data_Matthews', 
                    'critical_data_CRC', 'critical_data_PSRKR4',
                    'critical_data_Yaws', 'critical_data_PassutDanner',
                    'Tc_sources', 'Pc_sources', 'Vc_sources', 'Zc_sources',
                    'omega_sources'):
            _load_critical_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else: # pragma: no cover
    if can_load_data:
        _load_critical_data()

### Critical point functions

Tc_all_methods = (IUPAC, MATTHEWS, CRC, PSRK, PD, YAWS)
'''Tuple of method name keys. See the `Tc` for the actual references'''

def Tc_methods(CASRN):
    """Return all methods available to obtain Tc for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Tc with the given inputs.

    See Also
    --------
    Tc
    """
    if not _critical_data_loaded: _load_critical_data()
    return list_available_methods_from_df_dict(Tc_sources, CASRN, 'Tc')
    
def Tc(CASRN, get_methods=False, method=None):
    r'''This function handles the retrieval of a chemical's critical
    temperature. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Preferred sources are 'IUPAC' for organic chemicals, and 'MATTHEWS' for 
    inorganic chemicals. Function has data for approximately 1000 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Tc : float
        Critical temperature, [K]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'IUPAC', 'MATTHEWS', 
        'CRC', 'PSRK', 'PD', and 'YAWS'. All valid values are also held  
        in the list `Tc_all_methods`.

    Notes
    -----
    A total of seven sources are available for this function. They are:

        * 'IUPAC Organic Critical Properties', a series of critically evaluated
          experimental datum for organic compounds in [1]_, [2]_, [3]_, [4]_,
          [5]_, [6]_, [7]_, [8]_, [9]_, [10]_, [11]_, and [12]_.
        * 'Matthews Inorganic Critical Properties', a series of critically
          evaluated data for inorganic compounds in [13]_.
        * 'CRC Organic Critical Properties', a compillation of critically
          evaluated data by the TRC as published in [14]_.
        * 'PSRK Revision 4 Appendix', a compillation of experimental and
          estimated data published in [15]_.
        * 'Passut Danner 1973 Critical Properties', an older compillation of
          data published in [16]_
        * 'Yaws Critical Properties', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [17]_.

    Examples
    --------
    >>> Tc(CASRN='64-17-5')
    514.0

    See Also
    --------
    Tc_methods

    References
    ----------
    .. [1] Ambrose, Douglas, and Colin L. Young. "Vapor-Liquid Critical
       Properties of Elements and Compounds. 1. An Introductory Survey."
       Journal of Chemical & Engineering Data 41, no. 1 (January 1, 1996):
       154-154. doi:10.1021/je950378q.
    .. [2] Ambrose, Douglas, and Constantine Tsonopoulos. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 2. Normal Alkanes."
       Journal of Chemical & Engineering Data 40, no. 3 (May 1, 1995): 531-46.
       doi:10.1021/je00019a001.
    .. [3] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 3. Aromatic
       Hydrocarbons." Journal of Chemical & Engineering Data 40, no. 3
       (May 1, 1995): 547-58. doi:10.1021/je00019a002.
    .. [4] Gude, Michael, and Amyn S. Teja. "Vapor-Liquid Critical Properties
       of Elements and Compounds. 4. Aliphatic Alkanols." Journal of Chemical
       & Engineering Data 40, no. 5 (September 1, 1995): 1025-36.
       doi:10.1021/je00021a001.
    .. [5] Daubert, Thomas E. "Vapor-Liquid Critical Properties of Elements
       and Compounds. 5. Branched Alkanes and Cycloalkanes." Journal of
       Chemical & Engineering Data 41, no. 3 (January 1, 1996): 365-72.
       doi:10.1021/je9501548.
    .. [6] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 6. Unsaturated Aliphatic
       Hydrocarbons." Journal of Chemical & Engineering Data 41, no. 4
       (January 1, 1996): 645-56. doi:10.1021/je9501999.
    .. [7] Kudchadker, Arvind P., Douglas Ambrose, and Constantine Tsonopoulos.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 7. Oxygen
       Compounds Other Than Alkanols and Cycloalkanols." Journal of Chemical &
       Engineering Data 46, no. 3 (May 1, 2001): 457-79. doi:10.1021/je0001680.
    .. [8] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 8. Organic Sulfur,
       Silicon, and Tin Compounds (C + H + S, Si, and Sn)." Journal of Chemical
       & Engineering Data 46, no. 3 (May 1, 2001): 480-85.
       doi:10.1021/je000210r.
    .. [9] Marsh, Kenneth N., Colin L. Young, David W. Morton, Douglas Ambrose,
       and Constantine Tsonopoulos. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 9. Organic Compounds Containing Nitrogen."
       Journal of Chemical & Engineering Data 51, no. 2 (March 1, 2006):
       305-14. doi:10.1021/je050221q.
    .. [10] Marsh, Kenneth N., Alan Abramson, Douglas Ambrose, David W. Morton,
       Eugene Nikitin, Constantine Tsonopoulos, and Colin L. Young.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 10. Organic
       Compounds Containing Halogens." Journal of Chemical & Engineering Data
       52, no. 5 (September 1, 2007): 1509-38. doi:10.1021/je700336g.
    .. [11] Ambrose, Douglas, Constantine Tsonopoulos, and Eugene D. Nikitin.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 11. Organic
       Compounds Containing B + O; Halogens + N, + O, + O + S, + S, + Si;
       N + O; and O + S, + Si." Journal of Chemical & Engineering Data 54,
       no. 3 (March 12, 2009): 669-89. doi:10.1021/je800580z.
    .. [12] Ambrose, Douglas, Constantine Tsonopoulos, Eugene D. Nikitin, David
       W. Morton, and Kenneth N. Marsh. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 12. Review of Recent Data for Hydrocarbons and
       Non-Hydrocarbons." Journal of Chemical & Engineering Data, October 5,
       2015, 151005081500002. doi:10.1021/acs.jced.5b00571.
    .. [13] Mathews, Joseph F. "Critical Constants of Inorganic Substances."
       Chemical Reviews 72, no. 1 (February 1, 1972): 71-100.
       doi:10.1021/cr60275a004.
    .. [14] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [15] Horstmann, Sven, Anna Jabłoniec, Jörg Krafczyk, Kai Fischer, and
       Jürgen Gmehling. "PSRK Group Contribution Equation of State:
       Comprehensive Revision and Extension IV, Including Critical Constants
       and Α-Function Parameters for 1000 Components." Fluid Phase Equilibria
       227, no. 2 (January 25, 2005): 157-64. doi:10.1016/j.fluid.2004.11.002.
    .. [16] Passut, Charles A., and Ronald P. Danner. "Acentric Factor. A
       Valuable Correlating Parameter for the Properties of Hydrocarbons."
       Industrial & Engineering Chemistry Process Design and Development 12,
       no. 3 (July 1, 1973): 365–68. doi:10.1021/i260047a026.
    .. [17] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    
    '''
    if not _critical_data_loaded: _load_critical_data()
    if method:
        return retrieve_from_df_dict(Tc_sources, CASRN, 'Tc', method) 
    else:
        return retrieve_any_from_df_dict(Tc_sources, CASRN, 'Tc') 

Pc_all_methods = (IUPAC, MATTHEWS, CRC, PSRK, PD, YAWS)
'''Tuple of method name keys. See the `Pc` for the actual references'''

def Pc_methods(CASRN):
    """Return all methods available to obtain Pc for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Pc with the given inputs.

    See Also
    --------
    Pc
    """
    return list_available_methods_from_df_dict(Pc_sources, CASRN, 'Pc')
    

def Pc(CASRN, get_methods=False, method=None):
    r'''This function handles the retrieval of a chemical's critical
    pressure. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Preferred sources are 'IUPAC' for organic chemicals, and 'MATTHEWS' for 
    inorganic chemicals. Function has data for approximately 7500 chemicals.

    Examples
    --------
    >>> Pc(CASRN='64-17-5')
    6137000.0

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Pc : float
        Critical pressure, [Pa]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'IUPAC', 'MATTHEWS', 
        'CRC', 'PSRK', 'PD', and 'YAWS'. All valid values are also held  
        in the list `Pc_all_methods`.

    Notes
    -----
    A total of seven sources are available for this function. They are:

        * 'IUPAC', a series of critically evaluated
          experimental datum for organic compounds in [1]_, [2]_, [3]_, [4]_,
          [5]_, [6]_, [7]_, [8]_, [9]_, [10]_, [11]_, and [12]_.
        * 'MATTHEWS', a series of critically
          evaluated data for inorganic compounds in [13]_.
        * 'CRC', a compillation of critically
          evaluated data by the TRC as published in [14]_.
        * 'PSRK', a compillation of experimental and
          estimated data published in [15]_.
        * 'PD', an older compillation of
          data published in [16]_
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [17]_.

    See Also
    --------
    Pc_methods

    References
    ----------
    .. [1] Ambrose, Douglas, and Colin L. Young. "Vapor-Liquid Critical
       Properties of Elements and Compounds. 1. An Introductory Survey."
       Journal of Chemical & Engineering Data 41, no. 1 (January 1, 1996):
       154-154. doi:10.1021/je950378q.
    .. [2] Ambrose, Douglas, and Constantine Tsonopoulos. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 2. Normal Alkanes."
       Journal of Chemical & Engineering Data 40, no. 3 (May 1, 1995): 531-46.
       doi:10.1021/je00019a001.
    .. [3] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 3. Aromatic
       Hydrocarbons." Journal of Chemical & Engineering Data 40, no. 3
       (May 1, 1995): 547-58. doi:10.1021/je00019a002.
    .. [4] Gude, Michael, and Amyn S. Teja. "Vapor-Liquid Critical Properties
       of Elements and Compounds. 4. Aliphatic Alkanols." Journal of Chemical
       & Engineering Data 40, no. 5 (September 1, 1995): 1025-36.
       doi:10.1021/je00021a001.
    .. [5] Daubert, Thomas E. "Vapor-Liquid Critical Properties of Elements
       and Compounds. 5. Branched Alkanes and Cycloalkanes." Journal of
       Chemical & Engineering Data 41, no. 3 (January 1, 1996): 365-72.
       doi:10.1021/je9501548.
    .. [6] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 6. Unsaturated Aliphatic
       Hydrocarbons." Journal of Chemical & Engineering Data 41, no. 4
       (January 1, 1996): 645-56. doi:10.1021/je9501999.
    .. [7] Kudchadker, Arvind P., Douglas Ambrose, and Constantine Tsonopoulos.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 7. Oxygen
       Compounds Other Than Alkanols and Cycloalkanols." Journal of Chemical &
       Engineering Data 46, no. 3 (May 1, 2001): 457-79. doi:10.1021/je0001680.
    .. [8] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 8. Organic Sulfur,
       Silicon, and Tin Compounds (C + H + S, Si, and Sn)." Journal of Chemical
       & Engineering Data 46, no. 3 (May 1, 2001): 480-85.
       doi:10.1021/je000210r.
    .. [9] Marsh, Kenneth N., Colin L. Young, David W. Morton, Douglas Ambrose,
       and Constantine Tsonopoulos. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 9. Organic Compounds Containing Nitrogen."
       Journal of Chemical & Engineering Data 51, no. 2 (March 1, 2006):
       305-14. doi:10.1021/je050221q.
    .. [10] Marsh, Kenneth N., Alan Abramson, Douglas Ambrose, David W. Morton,
       Eugene Nikitin, Constantine Tsonopoulos, and Colin L. Young.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 10. Organic
       Compounds Containing Halogens." Journal of Chemical & Engineering Data
       52, no. 5 (September 1, 2007): 1509-38. doi:10.1021/je700336g.
    .. [11] Ambrose, Douglas, Constantine Tsonopoulos, and Eugene D. Nikitin.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 11. Organic
       Compounds Containing B + O; Halogens + N, + O, + O + S, + S, + Si;
       N + O; and O + S, + Si." Journal of Chemical & Engineering Data 54,
       no. 3 (March 12, 2009): 669-89. doi:10.1021/je800580z.
    .. [12] Ambrose, Douglas, Constantine Tsonopoulos, Eugene D. Nikitin, David
       W. Morton, and Kenneth N. Marsh. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 12. Review of Recent Data for Hydrocarbons and
       Non-Hydrocarbons." Journal of Chemical & Engineering Data, October 5,
       2015, 151005081500002. doi:10.1021/acs.jced.5b00571.
    .. [13] Mathews, Joseph F. "Critical Constants of Inorganic Substances."
       Chemical Reviews 72, no. 1 (February 1, 1972): 71-100.
       doi:10.1021/cr60275a004.
    .. [14] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [15] Horstmann, Sven, Anna Jabłoniec, Jörg Krafczyk, Kai Fischer, and
       Jürgen Gmehling. "PSRK Group Contribution Equation of State:
       Comprehensive Revision and Extension IV, Including Critical Constants
       and Α-Function Parameters for 1000 Components." Fluid Phase Equilibria
       227, no. 2 (January 25, 2005): 157-64. doi:10.1016/j.fluid.2004.11.002.
    .. [16] Passut, Charles A., and Ronald P. Danner. "Acentric Factor. A
       Valuable Correlating Parameter for the Properties of Hydrocarbons."
       Industrial & Engineering Chemistry Process Design and Development 12,
       no. 3 (July 1, 1973): 365–68. doi:10.1021/i260047a026.
    .. [17] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if not _critical_data_loaded: _load_critical_data()
    if method:
        return retrieve_from_df_dict(Pc_sources, CASRN, 'Pc', method) 
    else:
        return retrieve_any_from_df_dict(Pc_sources, CASRN, 'Pc') 

Vc_all_methods = (IUPAC, MATTHEWS, CRC, PSRK, YAWS)
'''Tuple of method name keys. See the `Vc` for the actual references'''

def Vc_methods(CASRN):
    """Return all methods available to obtain Vc for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Vc with the given inputs.

    See Also
    --------
    Vc
    """
    if not _critical_data_loaded: _load_critical_data()
    return list_available_methods_from_df_dict(Vc_sources, CASRN, 'Vc')

def Vc(CASRN, get_methods=False, method=None):
    r'''This function handles the retrieval of a chemical's critical
    volume. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Preferred sources are 'IUPAC' for organic chemicals, and 'MATTHEWS' for 
    inorganic chemicals. Function has data for approximately 7500 chemicals.

    Examples
    --------
    >>> Vc(CASRN='64-17-5')
    0.000168

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Vc : float
        Critical volume, [m^3/mol]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'IUPAC', 'MATTHEWS', 
        'CRC', 'PSRK', and 'YAWS'. All valid values are also held  
        in the list `Vc_all_methods`.
        
    Notes
    -----
    A total of six sources are available for this function. They are:

        * 'IUPAC', a series of critically evaluated
          experimental datum for organic compounds in [1]_, [2]_, [3]_, [4]_,
          [5]_, [6]_, [7]_, [8]_, [9]_, [10]_, [11]_, and [12]_.
        * 'MATTHEWS', a series of critically
          evaluated data for inorganic compounds in [13]_.
        * 'CRC', a compillation of critically
          evaluated data by the TRC as published in [14]_.
        * 'PSRK', a compillation of experimental and
          estimated data published in [15]_.
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [16]_.

    See Also
    --------
    Vc_methods

    References
    ----------
    .. [1] Ambrose, Douglas, and Colin L. Young. "Vapor-Liquid Critical
       Properties of Elements and Compounds. 1. An Introductory Survey."
       Journal of Chemical & Engineering Data 41, no. 1 (January 1, 1996):
       154-154. doi:10.1021/je950378q.
    .. [2] Ambrose, Douglas, and Constantine Tsonopoulos. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 2. Normal Alkanes."
       Journal of Chemical & Engineering Data 40, no. 3 (May 1, 1995): 531-46.
       doi:10.1021/je00019a001.
    .. [3] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 3. Aromatic
       Hydrocarbons." Journal of Chemical & Engineering Data 40, no. 3
       (May 1, 1995): 547-58. doi:10.1021/je00019a002.
    .. [4] Gude, Michael, and Amyn S. Teja. "Vapor-Liquid Critical Properties
       of Elements and Compounds. 4. Aliphatic Alkanols." Journal of Chemical
       & Engineering Data 40, no. 5 (September 1, 1995): 1025-36.
       doi:10.1021/je00021a001.
    .. [5] Daubert, Thomas E. "Vapor-Liquid Critical Properties of Elements
       and Compounds. 5. Branched Alkanes and Cycloalkanes." Journal of
       Chemical & Engineering Data 41, no. 3 (January 1, 1996): 365-72.
       doi:10.1021/je9501548.
    .. [6] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 6. Unsaturated Aliphatic
       Hydrocarbons." Journal of Chemical & Engineering Data 41, no. 4
       (January 1, 1996): 645-56. doi:10.1021/je9501999.
    .. [7] Kudchadker, Arvind P., Douglas Ambrose, and Constantine Tsonopoulos.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 7. Oxygen
       Compounds Other Than Alkanols and Cycloalkanols." Journal of Chemical &
       Engineering Data 46, no. 3 (May 1, 2001): 457-79. doi:10.1021/je0001680.
    .. [8] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 8. Organic Sulfur,
       Silicon, and Tin Compounds (C + H + S, Si, and Sn)." Journal of Chemical
       & Engineering Data 46, no. 3 (May 1, 2001): 480-85.
       doi:10.1021/je000210r.
    .. [9] Marsh, Kenneth N., Colin L. Young, David W. Morton, Douglas Ambrose,
       and Constantine Tsonopoulos. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 9. Organic Compounds Containing Nitrogen."
       Journal of Chemical & Engineering Data 51, no. 2 (March 1, 2006):
       305-14. doi:10.1021/je050221q.
    .. [10] Marsh, Kenneth N., Alan Abramson, Douglas Ambrose, David W. Morton,
       Eugene Nikitin, Constantine Tsonopoulos, and Colin L. Young.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 10. Organic
       Compounds Containing Halogens." Journal of Chemical & Engineering Data
       52, no. 5 (September 1, 2007): 1509-38. doi:10.1021/je700336g.
    .. [11] Ambrose, Douglas, Constantine Tsonopoulos, and Eugene D. Nikitin.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 11. Organic
       Compounds Containing B + O; Halogens + N, + O, + O + S, + S, + Si;
       N + O; and O + S, + Si." Journal of Chemical & Engineering Data 54,
       no. 3 (March 12, 2009): 669-89. doi:10.1021/je800580z.
    .. [12] Ambrose, Douglas, Constantine Tsonopoulos, Eugene D. Nikitin, David
       W. Morton, and Kenneth N. Marsh. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 12. Review of Recent Data for Hydrocarbons and
       Non-Hydrocarbons." Journal of Chemical & Engineering Data, October 5,
       2015, 151005081500002. doi:10.1021/acs.jced.5b00571.
    .. [13] Mathews, Joseph F. "Critical Constants of Inorganic Substances."
       Chemical Reviews 72, no. 1 (February 1, 1972): 71-100.
       doi:10.1021/cr60275a004.
    .. [14] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [15] Horstmann, Sven, Anna Jabłoniec, Jörg Krafczyk, Kai Fischer, and
       Jürgen Gmehling. "PSRK Group Contribution Equation of State:
       Comprehensive Revision and Extension IV, Including Critical Constants
       and Α-Function Parameters for 1000 Components." Fluid Phase Equilibria
       227, no. 2 (January 25, 2005): 157-64. doi:10.1016/j.fluid.2004.11.002.
    .. [16] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if not _critical_data_loaded: _load_critical_data()
    if method:
        return retrieve_from_df_dict(Vc_sources, CASRN, 'Vc', method)
    else:
        return retrieve_any_from_df_dict(Vc_sources, CASRN, 'Vc') 

Zc_all_methods = (IUPAC, MATTHEWS, CRC, PSRK, YAWS)
'''Tuple of method name keys. See the `Zc` for the actual references'''

def Zc_methods(CASRN):
    """Return all methods available to obtain Zc for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain Zc with the given inputs.

    See Also
    --------
    Zc
    """
    if not _critical_data_loaded: _load_critical_data()
    return list_available_methods_from_df_dict(Zc_sources, CASRN, 'Zc')

def Zc(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's critical
    compressibility. Lookup is based on CASRNs. Will automatically select a
    data source to use if no method is provided; returns None if the data is
    not available.

    Preferred sources are 'IUPAC' for organic chemicals, and 'MATTHEWS' for 
    inorganic chemicals. Function has data for approximately 7500 chemicals.

    Examples
    --------
    >>> Zc(CASRN='64-17-5')
    0.24100000000000002

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Zc : float
        Critical compressibility, [-]

    Other Parameters
    ----------------
    method : string, optional
        The method name to use. Accepted methods are 'IUPAC', 'MATTHEWS', 
        'CRC', 'PSRK', and 'YAWS'. All valid values are also held  
        in `Zc_all_methods`.

    Notes
    -----
    A total of five sources are available for this function. They are:

        * 'IUPAC', a series of critically evaluated
          experimental datum for organic compounds in [1]_, [2]_, [3]_, [4]_,
          [5]_, [6]_, [7]_, [8]_, [9]_, [10]_, [11]_, and [12]_.
        * 'MATTHEWS', a series of critically
          evaluated data for inorganic compounds in [13]_.
        * 'CRC', a compillation of critically
          evaluated data by the TRC as published in [14]_.
        * 'PSRK', a compillation of experimental and
          estimated data published in [15]_.
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [16]_.

    See Also
    --------
    Zc_methods

    References
    ----------
    .. [1] Ambrose, Douglas, and Colin L. Young. "Vapor-Liquid Critical
       Properties of Elements and Compounds. 1. An Introductory Survey."
       Journal of Chemical & Engineering Data 41, no. 1 (January 1, 1996):
       154-154. doi:10.1021/je950378q.
    .. [2] Ambrose, Douglas, and Constantine Tsonopoulos. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 2. Normal Alkanes."
       Journal of Chemical & Engineering Data 40, no. 3 (May 1, 1995): 531-46.
       doi:10.1021/je00019a001.
    .. [3] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 3. Aromatic
       Hydrocarbons." Journal of Chemical & Engineering Data 40, no. 3
       (May 1, 1995): 547-58. doi:10.1021/je00019a002.
    .. [4] Gude, Michael, and Amyn S. Teja. "Vapor-Liquid Critical Properties
       of Elements and Compounds. 4. Aliphatic Alkanols." Journal of Chemical
       & Engineering Data 40, no. 5 (September 1, 1995): 1025-36.
       doi:10.1021/je00021a001.
    .. [5] Daubert, Thomas E. "Vapor-Liquid Critical Properties of Elements
       and Compounds. 5. Branched Alkanes and Cycloalkanes." Journal of
       Chemical & Engineering Data 41, no. 3 (January 1, 1996): 365-72.
       doi:10.1021/je9501548.
    .. [6] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 6. Unsaturated Aliphatic
       Hydrocarbons." Journal of Chemical & Engineering Data 41, no. 4
       (January 1, 1996): 645-56. doi:10.1021/je9501999.
    .. [7] Kudchadker, Arvind P., Douglas Ambrose, and Constantine Tsonopoulos.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 7. Oxygen
       Compounds Other Than Alkanols and Cycloalkanols." Journal of Chemical &
       Engineering Data 46, no. 3 (May 1, 2001): 457-79. doi:10.1021/je0001680.
    .. [8] Tsonopoulos, Constantine, and Douglas Ambrose. "Vapor-Liquid
       Critical Properties of Elements and Compounds. 8. Organic Sulfur,
       Silicon, and Tin Compounds (C + H + S, Si, and Sn)." Journal of Chemical
       & Engineering Data 46, no. 3 (May 1, 2001): 480-85.
       doi:10.1021/je000210r.
    .. [9] Marsh, Kenneth N., Colin L. Young, David W. Morton, Douglas Ambrose,
       and Constantine Tsonopoulos. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 9. Organic Compounds Containing Nitrogen."
       Journal of Chemical & Engineering Data 51, no. 2 (March 1, 2006):
       305-14. doi:10.1021/je050221q.
    .. [10] Marsh, Kenneth N., Alan Abramson, Douglas Ambrose, David W. Morton,
       Eugene Nikitin, Constantine Tsonopoulos, and Colin L. Young.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 10. Organic
       Compounds Containing Halogens." Journal of Chemical & Engineering Data
       52, no. 5 (September 1, 2007): 1509-38. doi:10.1021/je700336g.
    .. [11] Ambrose, Douglas, Constantine Tsonopoulos, and Eugene D. Nikitin.
       "Vapor-Liquid Critical Properties of Elements and Compounds. 11. Organic
       Compounds Containing B + O; Halogens + N, + O, + O + S, + S, + Si;
       N + O; and O + S, + Si." Journal of Chemical & Engineering Data 54,
       no. 3 (March 12, 2009): 669-89. doi:10.1021/je800580z.
    .. [12] Ambrose, Douglas, Constantine Tsonopoulos, Eugene D. Nikitin, David
       W. Morton, and Kenneth N. Marsh. "Vapor-Liquid Critical Properties of
       Elements and Compounds. 12. Review of Recent Data for Hydrocarbons and
       Non-Hydrocarbons." Journal of Chemical & Engineering Data, October 5,
       2015, 151005081500002. doi:10.1021/acs.jced.5b00571.
    .. [13] Mathews, Joseph F. "Critical Constants of Inorganic Substances."
       Chemical Reviews 72, no. 1 (February 1, 1972): 71-100.
       doi:10.1021/cr60275a004.
    .. [14] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [15] Horstmann, Sven, Anna Jabłoniec, Jörg Krafczyk, Kai Fischer, and
       Jürgen Gmehling. "PSRK Group Contribution Equation of State:
       Comprehensive Revision and Extension IV, Including Critical Constants
       and Α-Function Parameters for 1000 Components." Fluid Phase Equilibria
       227, no. 2 (January 25, 2005): 157-64. doi:10.1016/j.fluid.2004.11.002.
    .. [16] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if not _critical_data_loaded: _load_critical_data()
    if method:
        return retrieve_from_df_dict(Zc_sources, CASRN, 'Zc', method)
    else:
        return retrieve_any_from_df_dict(Zc_sources, CASRN, 'Zc') 

rcovs_Mersmann_Kind = {'C': 0.77, 'Cl': 0.99, 'I': 1.33, 'H': 0.37, 'F': 0.71, 
                       'S': 1.04, 'O': 0.6, 'N': 0.71, 'Si': 1.17, 'Br': 1.14}

rcovs_regressed =  {
    u'Nb': 0.5139380605234125,
    u'Ne': 0.7708216694154189,
    u'Al': 1.004994775098707,
    u'Re': 1.1164444694484814,
    u'Rb': 2.9910506044828837,
    u'Rn': 1.9283158156480653,
    u'Xe': 1.694221043013319,
    u'Ta': 1.1185133195453156,
    u'Bi': 1.8436438207262267,
    u'Br': 1.3081458724155532,
    u'Hf': 0.8829545460486594,
    u'Mo': 0.740396259301556,
    u'He': 0.9808144122544257,
    u'C': 0.6068586007600608,
    u'B': 0.7039677272439753,
    u'F': 0.5409105884533288,
    u'I': 1.7262432419406561,
    u'H': 0.33296601702348533,
    u'K': 0.7384112258842432,
    u'O': 0.5883254088243008,
    u'N': 0.5467979701131293,
    u'P': 1.0444655158949694,
    u'Si': 1.4181434041348049,
    u'U': 1.5530287578073485,
    u'Sn': 1.3339487990207999,
    u'W': 0.8355335838735266,
    u'V': 0.6714619384794069,
    u'Sb': 0.8840680681215854,
    u'Se': 1.5747549515496795,
    u'Ge': 1.0730584829731715,
    u'Kr': 1.393999829252709,
    u'Cl': 1.0957835025011224,
    u'S': 1.0364452121761167,
    u'Hg': 0.6750818243474633,
    u'As': 0.6750687692915264,
    u'Ar': 1.2008872952022298,
    u'Cs': 3.433699060142929,
    u'Zr': 0.9346554283483623}

def Mersmann_Kind_predictor(atoms, coeff=3.645, power=0.5, 
                            covalent_radii=rcovs_Mersmann_Kind):
    r'''Predicts the critical molar volume of a chemical based only on its
    atomic composition according to [1]_ and [2]_. This is a crude approach,
    but provides very reasonable
    estimates in practice. Optionally, the `coeff` used and the `power` in the
    fraction as well as the atomic contributions can be adjusted; this method
    is general and atomic contributions can be regressed to predict other
    properties with this routine.
    
    .. math::
        \frac{\left(\frac{V_c}{n_a N_A}\right)^{1/3}}{d_a}
        = \frac{3.645}{\left(\frac{r_a}{r_H}\right)^{1/2}}

        r_a = d_a/2
        
        d_a = 2 \frac{\sum_i (n_i r_i)}{n_a}
        
    In the above equations, :math:`n_i` is the number of atoms of species i in
    the molecule, :math:`r_i` is the covalent atomic radius of the atom, and 
    :math:`n_a` is the total number of atoms in the molecule.
    
    Parameters
    ----------
    atoms : dict
        Dictionary of atoms and their counts, [-]
    coeff : float, optional
        Coefficient used in the relationship, [m^2]
    power : float, optional
        Power applied to the relative atomic radius, [-]
    covalent_radii : dict or indexable, optional
        Object which can be indexed to atomic contrinbutions (by symbol), [-]

    Returns
    -------
    Vc : float
        Predicted critical volume of the chemical, [m^3/mol]
    
    Notes
    -----    
    Using the :obj:`chemicals.elements.periodic_table` covalent radii (from RDKit), 
    the coefficient and power should be 4.261206523632586 and 0.5597281770786228
    respectively for best results.
    
    Examples
    --------
    Prediction of critical volume of decane:
        
    >>> Mersmann_Kind_predictor({'C': 10, 'H': 22})
    0.0005851858957767497
    
    This is compared against the experimental value, 0.000624 (a 6.2% relative
    error)
    
    Using custom fitted coefficients we can do a bit better:
        
    >>> from chemicals.critical import rcovs_regressed
    >>> Mersmann_Kind_predictor({'C': 10, 'H': 22}, coeff=4.261206523632586, 
    ... power=0.5597281770786228, covalent_radii=rcovs_regressed)
    0.0005956870915974391
    
    The relative error is only 4.5% now. This is compared to an experimental 
    uncertainty of 5.6%.
    
    Evaluating 1321 critical volumes in the database, the average relative
    error is 5.0%; standard deviation 6.8%; and worst value of 79% relative
    error for phosphorus.
    
    References
    ----------
    .. [1] Mersmann, Alfons, and Matthias Kind. "Correlation for the Prediction
       of Critical Molar Volume." Industrial & Engineering Chemistry Research,
       October 16, 2017. https://doi.org/10.1021/acs.iecr.7b03171.
    .. [2] Mersmann, Alfons, and Matthias Kind. "Prediction of Mechanical and 
       Thermal Properties of Pure Liquids, of Critical Data, and of Vapor 
       Pressure." Industrial & Engineering Chemistry Research, January 31, 
       2017. https://doi.org/10.1021/acs.iecr.6b04323.
    '''
    H_RADIUS_COV = covalent_radii['H']
    tot = 0
    atom_count = 0
    for atom, count in atoms.items():
        if atom not in covalent_radii:
            raise Exception('Atom %s is not supported by the supplied dictionary' %atom)
        tot += count*covalent_radii[atom]
        atom_count += count
    da = 2.*tot/atom_count
    ra = da/2.
    da_SI = da*1e-10 # Convert from angstrom to m
    return ((coeff/(ra/H_RADIUS_COV)**power)*da_SI)**3*N_A*atom_count

### Critical Property Relationships

def _assert_two_critical_properties_provided(Tc, Pc, Vc):
    specs = 0 # numba compatibility
    if Tc is not None:
        specs += 1
    if Pc is not None:
        specs += 1
    if Vc is not None:
        specs += 1
    if specs != 2:
        raise ValueError('Two and only two of Tc, Pc, and Vc must be provided')

def Ihmels(Tc=None, Pc=None, Vc=None):
    r'''Most recent, and most recommended method of estimating critical
    properties from each other. Two of the three properties are required.
    This model uses the "critical surface", a general plot of Tc vs Pc vs Vc.
    The model used 421 organic compounds to derive equation.
    The general equation is in [1]_:

    .. math::
        P_c = -0.025 + 2.215 \frac{T_c}{V_c}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid (optional) [K]
    Pc : float
        Critical pressure of fluid (optional) [Pa]
    Vc : float
        Critical volume of fluid (optional) [m^3/mol]

    Returns
    -------
    Tc, Pc or Vc : float
        Critical property of fluid [K], [Pa], or [m^3/mol]

    Notes
    -----
    The prediction of Tc from Pc and Vc is not tested, as this is not necessary
    anywhere, but it is implemented.
    Internal units are MPa, cm^3/mol, and K. A slight error occurs when
    Pa, cm^3/mol and K are used instead, on the order of <0.2%.
    Their equation was also compared with 56 inorganic and elements.
    Devations of 20% for <200K or >1000K points.

    Examples
    --------
    Succinic acid [110-15-6]

    >>> Ihmels(Tc=851.0, Vc=0.000308)
    6095016.233766234

    References
    ----------
    .. [1] Ihmels, E. Christian. "The Critical Surface." Journal of Chemical
       & Engineering Data 55, no. 9 (September 9, 2010): 3474-80.
       doi:10.1021/je100167w.
    '''
    _assert_two_critical_properties_provided(Tc, Pc, Vc)
    if Tc is not None and Vc is not None:
        Vc = Vc*1E6  # m^3/mol to cm^3/mol
        Pc_calc = -0.025+2.215*Tc/Vc
        Pc_calc = Pc_calc*1E6  # MPa to Pa
        return Pc_calc
    elif Tc is not None and Pc is not None:
        Pc = Pc*1e-6  # Pa to MPa
        Vc_calc = 443.0*Tc/(200.0*Pc+5.0)
        Vc_calc = Vc_calc/1E6  # cm^3/mol to m^3/mol
        return Vc_calc
    else: # Pc and Vc
        Pc = Pc*1e-6  # Pa to MPa
        Vc = Vc*1E6  # m^3/mol to cm^3/mol
        Tc_calc = 5.0/443.0*(40.0*Pc*Vc + Vc)
        return Tc_calc

def Meissner(Tc=None, Pc=None, Vc=None):
    r'''Old (1942) relationship for estimating critical
    properties from each other. Two of the three properties are required.
    This model uses the "critical surface", a general plot of Tc vs Pc vs Vc.
    The model used 42 organic and inorganic compounds to derive the equation.
    The general equation is in [1]_:

    .. math::
        P_c = \frac{2.08 T_c}{V_c-8}

    Parameters
    ----------
    Tc : float, optional
        Critical temperature of fluid [K]
    Pc : float, optional
        Critical pressure of fluid [Pa]
    Vc : float, optional
        Critical volume of fluid [m^3/mol]

    Returns
    -------
    Tc, Pc or Vc : float
        Critical property of fluid [K], [Pa], or [m^3/mol]

    Notes
    -----
    The prediction of Tc from Pc and Vc is not tested, as this is not necessary
    anywhere, but it is implemented.
    Internal units are atm, cm^3/mol, and K. A slight error occurs when
    Pa, cm^3/mol and K are used instead, on the order of <0.2%.
    This equation is less accurate than that of Ihmels, but surprisingly close.
    The author also proposed means of estimated properties independently.

    Examples
    --------
    Succinic acid [110-15-6]

    >>> Meissner(Tc=851.0, Vc=0.000308)
    5978445.199999999

    References
    ----------
    .. [1] Meissner, H. P., and E. M. Redding. "Prediction of Critical
           Constants." Industrial & Engineering Chemistry 34, no. 5
           (May 1, 1942): 521-26. doi:10.1021/ie50389a003.
    '''
    _assert_two_critical_properties_provided(Tc, Pc, Vc)
    if Tc and Vc:
        Vc = Vc*1E6
        Pc = 20.8*Tc/(Vc-8)
        Pc = 101325*Pc  # atm to Pa
        return Pc
    elif Tc and Pc:
        Pc = Pc/101325.  # Pa to atm
        Vc = 104/5.0*Tc/Pc+8
        Vc = Vc/1E6  # cm^3/mol to m^3/mol
        return Vc
    else: # Pc and Vc
        Pc = Pc/101325.  # Pa to atm
        Vc = Vc*1E6  # m^3/mol to cm^3/mol
        Tc = 5./104.0*Pc*(Vc-8)
        return Tc

def Grigoras(Tc=None, Pc=None, Vc=None):
    r'''Relatively recent (1990) relationship for estimating critical
    properties from each other. Two of the three properties are required.
    This model uses the "critical surface", a general plot of Tc vs Pc vs Vc.
    The model used 137 organic and inorganic compounds to derive the equation.
    The general equation is in [1]_:

    .. math::
        P_c = 2.9 + 20.2 \frac{T_c}{V_c}

    Parameters
    ----------
    Tc : float, optional
        Critical temperature of fluid [K]
    Pc : float, optional
        Critical pressure of fluid [Pa]
    Vc : float, optional
        Critical volume of fluid [m^3/mol]

    Returns
    -------
    Tc, Pc or Vc : float
        Critical property of fluid [K], [Pa], or [m^3/mol]

    Notes
    -----
    The prediction of Tc from Pc and Vc is not tested, as this is not necessary
    anywhere, but it is implemented.
    Internal units are bar, cm^3/mol, and K. A slight error occurs when
    Pa, cm^3/mol and K are used instead, on the order of <0.2%.
    This equation is less accurate than that of Ihmels, but surprisingly close.
    The author also investigated an early QSPR model.

    Examples
    --------
    Succinic acid [110-15-6]

    >>> Grigoras(Tc=851.0, Vc=0.000308)
    5871233.766233766

    References
    ----------
    .. [1] Grigoras, Stelian. "A Structural Approach to Calculate Physical
           Properties of Pure Organic Substances: The Critical Temperature,
           Critical Volume and Related Properties." Journal of Computational
           Chemistry 11, no. 4 (May 1, 1990): 493-510.
           doi:10.1002/jcc.540110408
    '''
    _assert_two_critical_properties_provided(Tc, Pc, Vc)
    if Tc and Vc:
        Vc = Vc*1E6  # m^3/mol to cm^3/mol
        Pc = 2.9 + 20.2*Tc/Vc
        Pc = Pc*1E5  # bar to Pa
        return Pc
    elif Tc and Pc:
        Pc = Pc/1E5  # Pa to bar
        Vc = 202.0*Tc/(10*Pc-29.0)
        Vc = Vc/1E6  # cm^3/mol to m^3/mol
        return Vc
    else: # Pc and Vc
        Pc = Pc/1E5  # Pa to bar
        Vc = Vc*1E6  # m^3/mol to cm^3/mol
        Tc = 1.0/202*(10*Pc-29.0)*Vc
        return Tc

def Hekayati_Raeissi(MW, V_sat=None, Tc=None, Pc=None, Vc=None):
    r'''Estimation model for missing critical constants of a fluid
    according to [1]_. Based on the molecular weight and saturation
    molar volume of a fluid, and requires one of `Tc` or `Pc`. 
    Optionally, `Vc` can be provided to increase the accuracy of
    the prediction of `Tc` or `Pc` a little.

    Parameters
    ----------
    MW : float
        Molecular weight of fluid, [g/mol]
    V_sat : float, optional
        Molar volume of liquid at the saturation pressure of the fluid 
        at 298.15 K. Used if `Vc` is not provided. [m^3/mol]
    Tc : float, optional
        Critical temperature of fluid (optional) [K]
    Pc : float, optional
        Critical pressure of fluid (optional) [Pa]
    Vc : float, optional
        Critical volume of fluid (optional) [m^3/mol]

    Returns
    -------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    Vc : float
        Critical volume of fluid [m^3/mol]

    Notes
    -----
    Internal units are kPa, m^3/kmol, and K.

    Examples
    --------
    Toluene

    >>> Hekayati_Raeissi(MW=92.13842, V_sat=0.00010686, Pc=4108000.0)
    (599.7965819136947, 4108000.0, 0.000314909150453723)

    References
    ----------
    .. [1] Hekayati, Javad, and Sona Raeissi. "Estimation of the Critical 
       Properties of Compounds Using Volume-Based Thermodynamics." AIChE Journal
       n/a, no. n/a (n.d.): e17004. https://doi.org/10.1002/aic.17004.
    '''
    if Tc is None and Pc is None:
        raise ValueError("One of `Tc` or `Pc` are required.")
    if Pc is not None:
        Pc *= 1e-3 # Pa to kPa
    if Vc is not None:
        Vc *= 1e3 # m^3/mol to m^3/kmol
    if V_sat is not None:
        V_sat *= 1e3
    Tc_calc = Pc_calc = 0.0
    Vc_estimated = False
    
    
    if Vc is None:
        if V_sat is None:
            raise ValueError("V_sat is required when Vc is not provided")
        Vc_estimated = True
        # Direct estimation of Vc from V and MW
        # 3.2 Equations 30-33
        a3 = -2.636518E+00  
        b3 = 1.048744E-02  
        c3 = 1.222756E+00  
        d3 = 6.924349E-09  
        e3 = -1.102846E-05  
        f3 = 1.871249E-03  

        chi1 = a3
        chi2 = b3*MW + c3
        chi3 = MW*(MW*(d3*MW + e3) + f3)
        Vc = V_sat*(chi1*V_sat*V_sat + chi2) + chi3

    if Vc_estimated:
        # equqtions 21-24; slightly less accurate than Equations 25-28
        a1 = -1.434646E+02 
        b1 = 1.781802E-01 
        c1 = 4.196576E+01 
        d1 = -3.443955E-05 
        e1 = 9.812262E-03 
        f1 = 1.106396E-01 
        psi1 = a1
        psi2 = b1*MW + c1
        psi3 = d1*MW*MW + e1*MW + f1
        fact = V_sat*V_sat*(psi1*V_sat + psi2) + psi3
        if Pc is None:
            if Tc is not None:
                Pc_calc = R*Tc/fact
        elif Tc is None:
            if Pc is not None:
                Tc_calc = fact*Pc/R
    else:
        # Equations 25-28
        a2 = -1.822904E-3
        b2 = 1.947108E-5
        c2 = 3.794091E0
        d2 = -2.472981E-8
        xi1 = a2*MW
        xi2 = b2*MW*MW + c2
        xi3 = d2*MW*MW*MW
        fact = xi1*Vc*Vc + xi2*Vc + xi3
        if Tc is None:
            Tc_calc = fact*Pc/R
        elif Pc is None:
            Pc_calc = R*Tc/fact
    Tc_ans = Tc if Tc is not None else Tc_calc
    Pc_ans = Pc if Pc is not None else Pc_calc
    Pc_ans *= 1e3
    Vc *= 1e-3
    return (Tc_ans, Pc_ans, Vc)

IHMELS = 'IHMELS'
MEISSNER = 'MEISSNER'
GRIGORAS = 'GRIGORAS'
critical_surface_all_methods = (IHMELS, MEISSNER, GRIGORAS)

def critical_surface_methods(Tc=None, Pc=None, Vc=None):
    """Return all methods available to obtain the third critial property for the
    desired chemical.

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid (optional) [K].
    Pc : float
        Critical pressure of fluid (optional) [Pa].
    Vc : float
        Critical volume of fluid (optional) [m^3/mol].

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the third critical property with
        the given inputs.

    See Also
    --------
    critical_surface
    """
    if ((Tc is not None and Pc is not None) 
         or (Tc is not None and Vc is not None) 
         or (Pc is not None and Vc is not None)):
        return list(critical_surface_all_methods)
    else:
        return []

def critical_surface(Tc=None, Pc=None, Vc=None,
                     method=None):
    r'''Function for calculating a critical property of a substance from its
    other two critical properties. Calls functions Ihmels, Meissner, and
    Grigoras, each of which use a general 'Critical surface' type of equation.
    Limited accuracy is expected due to very limited theoretical backing.

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid (optional) [K].
    Pc : float
        Critical pressure of fluid (optional) [Pa].
    Vc : float
        Critical volume of fluid (optional) [m^3/mol].
    method : string
        Request calculation uses the requested method.

    Returns
    -------
    Tc, Pc or Vc : float
        Critical property of fluid [K], [Pa], or [m^3/mol].

    Examples
    --------
    Decamethyltetrasiloxane [141-62-8]

    >>> critical_surface(Tc=599.4, Pc=1.19E6, method='IHMELS')
    0.0010927333333333334
    
    See Also
    --------
    critical_surface_methods_methods
    
    '''
    if not method or method == IHMELS:
        return Ihmels(Tc=Tc, Pc=Pc, Vc=Vc)
    elif method == MEISSNER:
        return Meissner(Tc=Tc, Pc=Pc, Vc=Vc)
    elif method == GRIGORAS:
        return Grigoras(Tc=Tc, Pc=Pc, Vc=Vc)
    else:
        raise ValueError('invalid method %s'%method)

def third_property(CASRN=None, T=False, P=False, V=False):
    r'''Function for calculating a critical property of a substance from its
    other two critical properties, but retrieving the actual other critical
    values for convenient calculation.
    Calls functions Ihmels, Meissner, and
    Grigoras, each of which use a general 'Critical surface' type of equation.
    Limited accuracy is expected due to very limited theoretical backing.

    Parameters
    ----------
    CASRN : str
        The CAS number of the desired chemical
    T : bool
        Estimate critical temperature
    P : bool
        Estimate critical pressure
    V : bool
        Estimate critical volume

    Returns
    -------
    Tc, Pc or Vc : float
        Critical property of fluid [K], [Pa], or [m^3/mol]

    Examples
    --------
    Decamethyltetrasiloxane [141-62-8]
    
    >>> third_property('141-62-8', V=True)
    0.0010920041152263375

    Succinic acid [110-15-6]
    
    >>> third_property('110-15-6', P=True)
    6095016.233766234
    '''
    specs = 0
    if T:
        specs += 1
    if P:
        specs += 1
    if V:
        specs += 1
    if specs != 1:
        raise ValueError("Only one of the following arguments can be True: T, P, V")
    if V:
        Tc_value = Tc(CASRN)
        Pc_value = Pc(CASRN)
        if Tc_value is not None and Pc_value is not None:
            return critical_surface(Tc=Tc_value, Pc=Pc_value)
    elif P:
        Tc_value = Tc(CASRN)
        Vc_value = Vc(CASRN)
        if Tc_value is not None and Vc_value is not None:
            return critical_surface(Tc=Tc_value, Vc=Vc_value)
    else:
        Pc_value = Pc(CASRN)
        Vc_value = Vc(CASRN)
        if Pc_value is not None and Vc_value is not None:
            return critical_surface(Pc=Pc_value, Vc=Vc_value)
    # Fallthrough to this clause when one or both of the looked up properties fail
    raise ValueError("Could not find the required two properties")

### Crtical Temperature of Mixtures - Estimation routines

def Li(zs, Tcs, Vcs):
    r'''Calculates critical temperature of a mixture according to
    mixing rules in [1]_. Better than simple mixing rules.

    .. math::
        T_{cm} = \sum_{i=1}^n \Phi_i T_{ci}\\
        \Phi = \frac{x_i V_{ci}}{\sum_{j=1}^n x_j V_{cj}}

    Parameters
    ----------
    zs : array-like
        Mole fractions of all components
    Tcs : array-like
        Critical temperatures of all components, [K]
    Vcs : array-like
        Critical volumes of all components, [m^3/mol]

    Returns
    -------
    Tcm : float
        Critical temperatures of the mixture, [K]

    Notes
    -----
    Reviewed in many papers on critical mixture temperature.

    Second example is from Najafi (2015), for ethylene, Benzene, ethylbenzene.
    This is similar to but not identical to the result from the article. The
    experimental point is 486.9 K.

    2rd example is from Najafi (2015), for:
    butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    Its result is identical to that calculated in the article.

    Examples
    --------
    Nitrogen-Argon 50/50 mixture

    >>> Li([0.5, 0.5], [126.2, 150.8], [8.95e-05, 7.49e-05])
    137.40766423357667

    butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.

    >>> Li([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6],
    ... [0.000255, 0.000313, 0.000371])
    449.68261498555444

    References
    ----------
    .. [1] Li, C. C. "Critical Temperature Estimation for Simple Mixtures."
       The Canadian Journal of Chemical Engineering 49, no. 5
       (October 1, 1971): 709-10. doi:10.1002/cjce.5450490529.
    '''
    N = len(zs)
    denominator_inv = 0.0
    for i in range(N):
        denominator_inv += zs[i]*Vcs[i]
    denominator_inv = 1.0/denominator_inv
    Tcm = 0.0
    for i in range(N):
        Tcm += zs[i]*Vcs[i]*Tcs[i]*denominator_inv
    return Tcm

def Chueh_Prausnitz_Tc(zs, Tcs, Vcs, taus):
    r'''Calculates critical temperature of a mixture according to
    mixing rules in [1]_.

    .. math::
        T_{cm} = \sum_i^n \theta_i Tc_i + \sum_i^n\sum_j^n(\theta_i \theta_j
        \tau_{ij})T_{ref}

        \theta = \frac{x_i V_{ci}^{2/3}}{\sum_{j=1}^n x_j V_{cj}^{2/3}}

    For a binary mxiture, this simplifies to:

    .. math::
        T_{cm} = \theta_1T_{c1} + \theta_2T_{c2}  + 2\theta_1\theta_2\tau_{12}

    Parameters
    ----------
    zs : array-like
        Mole fractions of all components
    Tcs : array-like
        Critical temperatures of all components, [K]
    Vcs : array-like
        Critical volumes of all components, [m^3/mol]
    taus : array-like of shape `zs` by `zs`
        Interaction parameters, [-]

    Returns
    -------
    Tcm : float
        Critical temperatures of the mixture, [K]

    Notes
    -----
    All parameters, even if zero, must be given to this function.

    Examples
    --------
    butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.

    >>> Chueh_Prausnitz_Tc([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6],
    ... [0.000255, 0.000313, 0.000371], [[0, 1.92681, 6.80358],
    ... [1.92681, 0, 1.89312], [ 6.80358, 1.89312, 0]])
    450.122576472349

    References
    ----------
    .. [1] Chueh, P. L., and J. M. Prausnitz. "Vapor-Liquid Equilibria at High
       Pressures: Calculation of Critical Temperatures, Volumes, and Pressures
       of Nonpolar Mixtures." AIChE Journal 13, no. 6 (November 1, 1967):
       1107-13. doi:10.1002/aic.690130613.
    .. [2] Najafi, Hamidreza, Babak Maghbooli, and Mohammad Amin Sobati.
       "Prediction of True Critical Temperature of Multi-Component Mixtures:
       Extending Fast Estimation Methods." Fluid Phase Equilibria 392
       (April 25, 2015): 104-26. doi:10.1016/j.fluid.2015.02.001.
    '''
    N = len(zs)
    denominator_inv = 0.0
    zi_Vc_23s = [0.0]*N
    for i in range(N):
        v = zs[i]*Vcs[i]**(2.0/3.)
        zi_Vc_23s[i] = v
        denominator_inv += v
    denominator_inv = 1.0/denominator_inv
    Tcm = 0.0
    for i in range(N):
        Tcm += zi_Vc_23s[i]*Tcs[i]
        for j in range(N):
            Tcm += zi_Vc_23s[i]*(zi_Vc_23s[j]*denominator_inv)*taus[i][j]
    Tcm *= denominator_inv
    return Tcm

def Grieves_Thodos(zs, Tcs, Aijs):
    r'''Calculates critical temperature of a mixture according to
    mixing rules in [1]_.

    .. math::
        T_{cm} = \sum_{i} \frac{T_{ci}}{1 + (1/x_i)\sum_j A_{ij} x_j}

    For a binary mxiture, this simplifies to:

    .. math::
        T_{cm} = \frac{T_{c1}}{1 + (x_2/x_1)A_{12}} +  \frac{T_{c2}}
        {1 + (x_1/x_2)A_{21}}

    Parameters
    ----------
    zs : array-like
        Mole fractions of all components
    Tcs : array-like
        Critical temperatures of all components, [K]
    Aijs : array-like of shape `zs` by `zs`
        Interaction parameters

    Returns
    -------
    Tcm : float
        Critical temperatures of the mixture, [K]

    Notes
    -----
    All parameters, even if zero, must be given to this function.
    Giving 0s gives really bad results however.

    Examples
    --------
    butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.

    >>> Grieves_Thodos([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6], [[0, 1.2503, 1.516], [0.799807, 0, 1.23843], [0.659633, 0.807474, 0]])
    450.1839618758971

    References
    ----------
    .. [1] Grieves, Robert B., and George Thodos. "The Critical Temperatures of
       Multicomponent Hydrocarbon Systems." AIChE Journal 8, no. 4
       (September 1, 1962): 550-53. doi:10.1002/aic.690080426.
    .. [2] Najafi, Hamidreza, Babak Maghbooli, and Mohammad Amin Sobati.
       "Prediction of True Critical Temperature of Multi-Component Mixtures:
       Extending Fast Estimation Methods." Fluid Phase Equilibria 392
       (April 25, 2015): 104-26. doi:10.1016/j.fluid.2015.02.001.
    '''
    N = len(zs)
    Tcm = 0.0
    for i in range(N):
        tot = 0.0
        row = Aijs[i]
        for j in range(N):
            tot += row[j]*zs[j]
        Tcm += Tcs[i]/(1. + 1./zs[i]*tot)
    return Tcm

def modified_Wilson_Tc(zs, Tcs, Aijs):
    r'''Calculates critical temperature of a mixture according to
    mixing rules in [1]_. Equation

    .. math::
        T_{cm} = \sum_i x_i T_{ci} + C\sum_i x_i \ln \left(x_i + \sum_j x_j A_{ij}\right)T_{ref}

    For a binary mxiture, this simplifies to:

    .. math::
        T_{cm} = x_1 T_{c1} + x_2 T_{c2} + C[x_1 \ln(x_1 + x_2A_{12}) + x_2\ln(x_2 + x_1 A_{21})]

    Parameters
    ----------
    zs : float
        Mole fractions of all components
    Tcs : float
        Critical temperatures of all components, [K]
    Aijs : matrix
        Interaction parameters

    Returns
    -------
    Tcm : float
        Critical temperatures of the mixture, [K]

    Notes
    -----
    The equation and original article has been reviewed.
    [1]_ has 75 binary systems, and additional multicomponent mixture parameters.
    All parameters, even if zero, must be given to this function.

    2rd example is from [2]_, for:
    butane/pentane/hexane 0.6449/0.2359/0.1192 mixture, exp: 450.22 K.
    Its result is identical to that calculated in the article.

    Examples
    --------
    >>> modified_Wilson_Tc([0.6449, 0.2359, 0.1192], [425.12, 469.7, 507.6],
    ... [[0, 1.174450, 1.274390], [0.835914, 0, 1.21038],
    ... [0.746878, 0.80677, 0]])
    450.03059668230316

    References
    ----------
    .. [1] Teja, Amyn S., Kul B. Garg, and Richard L. Smith. "A Method for the
       Calculation of Gas-Liquid Critical Temperatures and Pressures of
       Multicomponent Mixtures." Industrial & Engineering Chemistry Process
       Design and Development 22, no. 4 (1983): 672-76.
    .. [2] Najafi, Hamidreza, Babak Maghbooli, and Mohammad Amin Sobati.
       "Prediction of True Critical Temperature of Multi-Component Mixtures:
       Extending Fast Estimation Methods." Fluid Phase Equilibria 392
       (April 25, 2015): 104-26. doi:10.1016/j.fluid.2015.02.001.
    '''
    N = len(zs)
    Tcm = 0.0
    for i in range(N):
        Tcm += zs[i]*Tcs[i]
    Tcm_add = 0.0
    for i in range(N):
        tot = 0.0
        row = Aijs[i]
        for j in range(N):
            tot += zs[j]*row[j]
        Tcm_add += zs[i]*log(zs[i] + tot)
    Tcm += -2500.0*Tcm_add
    return Tcm

### Crtical Volume of Mixtures
def Chueh_Prausnitz_Vc(zs, Vcs, nus):
    r'''Calculates critical volume of a mixture according to
    mixing rules in [1]_ with an interaction parameter.

    .. math::
        V_{cm} = \sum_i^n \theta_i V_{ci} + \sum_i^n\sum_j^n(\theta_i \theta_j \nu_{ij})V_{ref}
        \theta = \frac{x_i V_{ci}^{2/3}}{\sum_{j=1}^n x_j V_{cj}^{2/3}}

    Parameters
    ----------
    zs : float
        Mole fractions of all components
    Vcs : float
        Critical volumes of all components, [m^3/mol]
    nus : matrix
        Interaction parameters, [cm^3/mol]

    Returns
    -------
    Vcm : float
        Critical volume of the mixture, [m^3/mol]

    Notes
    -----
    All parameters, even if zero, must be given to this function.
    nu parameters are in cm^3/mol, but are converted to m^3/mol inside the function


    Examples
    --------
    1-butanol/benzene 0.4271/0.5729 mixture, Vcm = 268.096 mL/mol.

    >>> Chueh_Prausnitz_Vc([0.4271, 0.5729], [0.000273, 0.000256], [[0, 5.61847], [5.61847, 0]])
    0.00026620503424517445

    References
    ----------
    .. [1] Chueh, P. L., and J. M. Prausnitz. "Vapor-Liquid Equilibria at High
       Pressures: Calculation of Critical Temperatures, Volumes, and Pressures
       of Nonpolar Mixtures." AIChE Journal 13, no. 6 (November 1, 1967):
       1107-13. doi:10.1002/aic.690130613.
    .. [2] Najafi, Hamidreza, Babak Maghbooli, and Mohammad Amin Sobati.
       "Prediction of True Critical Volume of Multi-Component Mixtures:
       Extending Fast Estimation Methods." Fluid Phase Equilibria 386
       (January 25, 2015): 13-29. doi:10.1016/j.fluid.2014.11.008.
    '''
    N = len(zs)
    denominator_inv = 0.0
    zi_Vc_23s = [0.0]*N
    for i in range(N):
        v = zs[i]*Vcs[i]**(2.0/3.)
        zi_Vc_23s[i] = v
        denominator_inv += v
    denominator_inv = 1.0/denominator_inv
    
    Vcm = 0.0
    for i in range(N):
        Vcm += zi_Vc_23s[i]*Vcs[i]
        Vcm_tot = 0.0
        for j in range(N):
            Vcm_tot += zi_Vc_23s[i]*zi_Vc_23s[j]*nus[i][j]
        Vcm += Vcm_tot*1e-6*denominator_inv
    Vcm *= denominator_inv
    return Vcm


def modified_Wilson_Vc(zs, Vcs, Aijs):
    r'''Calculates critical volume of a mixture according to
    mixing rules in [1]_ with parameters. Equation

    .. math::
        V_{cm} = \sum_i x_i V_{ci} + C\sum_i x_i \ln \left(x_i + \sum_j x_j A_{ij}\right)V_{ref}

    For a binary mxiture, this simplifies to:

    .. math::
        V_{cm} = x_1 V_{c1} + x_2 V_{c2} + C[x_1 \ln(x_1 + x_2A_{12}) + x_2\ln(x_2 + x_1 A_{21})]

    Parameters
    ----------
    zs : float
        Mole fractions of all components
    Vcs : float
        Critical volumes of all components, [m^3/mol]
    Aijs : matrix
        Interaction parameters, [cm^3/mol]

    Returns
    -------
    Vcm : float
        Critical volume of the mixture, [m^3/mol]

    Notes
    -----
    The equation and original article has been reviewed.
    All parameters, even if zero, must be given to this function.
    C = -2500

    All parameters, even if zero, must be given to this function.
    nu parameters are in cm^3/mol, but are converted to m^3/mol inside the function


    Examples
    --------
    1-butanol/benzene 0.4271/0.5729 mixture, Vcm = 268.096 mL/mol.

    >>> modified_Wilson_Vc([0.4271, 0.5729], [0.000273, 0.000256],
    ... [[0, 0.6671250], [1.3939900, 0]])
    0.0002664335032706881

    References
    ----------
    .. [1] Teja, Amyn S., Kul B. Garg, and Richard L. Smith. "A Method for the
       Calculation of Gas-Liquid Critical Temperatures and Pressures of
       Multicomponent Mixtures." Industrial & Engineering Chemistry Process
       Design and Development 22, no. 4 (1983): 672-76.
    .. [2] Najafi, Hamidreza, Babak Maghbooli, and Mohammad Amin Sobati.
       "Prediction of True Critical Temperature of Multi-Component Mixtures:
       Extending Fast Estimation Methods." Fluid Phase Equilibria 392
       (April 25, 2015): 104-26. doi:10.1016/j.fluid.2015.02.001.
    '''
    N = len(zs)
    Vcm = 0.0
    for i in range(N):
        Vcm += zs[i]*Vcs[i]
    Vcm_add = 0.0
    for i in range(N):
        tot = 0.0
        row = Aijs[i]
        for j in range(N):
            tot += zs[j]*row[j]
        Vcm_add += zs[i]*log(zs[i] + tot)
    Vcm += -2500.0*Vcm_add*1e-6
    return Vcm
