# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
SOFTWARE.'''

from __future__ import division

__all__ = ['Tt_methods', 'Tt', 'Pt_methods', 'Pt']

import os
import numpy as np
import pandas as pd

from chemicals.utils import isnan, PY37
from chemicals.phase_change import Tm
from chemicals.data_reader import register_df_source, data_source

folder = os.path.join(os.path.dirname(__file__), 'Triple Properties')


register_df_source(folder, name='Staveley 1981.tsv')

_triple_dfs_loaded = False
def load_triple_dfs():
    global Staveley_data, _triple_dfs_loaded
    Staveley_data = data_source('Staveley 1981.tsv')
    _triple_dfs_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('Staveley_data',):
            load_triple_dfs()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    load_triple_dfs()

STAVELEY = 'STAVELEY'
MELTING = 'MELTING'
NONE = 'NONE'

Tt_methods = [STAVELEY, MELTING]


def Tt(CASRN, AvailableMethods=False, Method=None):
    r'''This function handles the retrieval of a chemical's triple temperature.
    Lookup is based on CASRNs. Will automatically select a data source to use
    if no Method is provided; returns None if the data is not available.

    Returns data from [1]_, or a chemical's melting point if available.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    Tt : float
        Triple point temperature, [K]
    methods : list, only returned if AvailableMethods == True
        List of methods which can be used to obtain Tt with the
        given inputs

    Other Parameters
    ----------------
    Method : string, optional
        A string for the method name to use, as defined by constants in
        Tt_methods
    AvailableMethods : bool, optional
        If True, function will determine which methods can be used to obtain
        the Tt for the desired chemical, and will return methods
        instead of the Tt

    Notes
    -----
    Median difference between melting points and triple points is 0.02 K.
    Accordingly, this should be more than good enough for engineering
    applications.

    Temperatures are on the ITS-68 scale.

    Examples
    --------
    Ammonia

    >>> Tt('7664-41-7')
    195.48

    References
    ----------
    .. [1] Staveley, L. A. K., L. Q. Lobo, and J. C. G. Calado. "Triple-Points
       of Low Melting Substances and Their Use in Cryogenic Work." Cryogenics
       21, no. 3 (March 1981): 131-144. doi:10.1016/0011-2275(81)90264-2.
    '''
    if not _triple_dfs_loaded:
        load_triple_dfs()
    def list_methods():
        methods = []
        if CASRN in Staveley_data.index:
            methods.append(STAVELEY)
        if Tm(CASRN):
            methods.append(MELTING)
        methods.append(NONE)
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]

    if Method == STAVELEY:
        return Staveley_data.at[CASRN, "Tt68"]
    elif Method == MELTING:
        return Tm(CASRN)
    elif Method == NONE:
        return None
    else:
        raise Exception('Failure in in function')

Pt_methods = [STAVELEY]


def Pt(CASRN, AvailableMethods=False, Method=None):
    r'''This function handles the retrieval of a chemical's triple pressure.
    Lookup is based on CASRNs. Will automatically select a data source to use
    if no Method is provided; returns None if the data is not available.

    Returns data from [1]_ only. 
    
    This function doe snot implement it but it is also possible to calculate 
    the vapor pressure at the triple temperature from a vapor pressure
    correlation, if data is available; note most Antoine-type correlations do
    not extrapolate well to this low of a pressure.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    Pt : float
        Triple point pressure, [Pa]
    methods : list, only returned if AvailableMethods == True
        List of methods which can be used to obtain Pt with the
        given inputs

    Other Parameters
    ----------------
    Method : string, optional
        A string for the method name to use, as defined by constants in
        Pt_methods
    AvailableMethods : bool, optional
        If True, function will determine which methods can be used to obtain
        the Pt for the desired chemical, and will return methods
        instead of the Pt

    Notes
    -----

    Examples
    --------
    Ammonia

    >>> Pt('7664-41-7')
    6079.5

    References
    ----------
    .. [1] Staveley, L. A. K., L. Q. Lobo, and J. C. G. Calado. "Triple-Points
       of Low Melting Substances and Their Use in Cryogenic Work." Cryogenics
       21, no. 3 (March 1981): 131-144. doi:10.1016/0011-2275(81)90264-2.
    '''
    if not _triple_dfs_loaded:
        load_triple_dfs()
    def list_methods():
        methods = []
        if CASRN in Staveley_data.index and not isnan(Staveley_data.at[CASRN, 'Pt']):
            methods.append(STAVELEY)
        methods.append(NONE)
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]

    if Method == STAVELEY:
        return Staveley_data.at[CASRN, 'Pt']
    elif Method == NONE:
        return None
    else:
        raise Exception('Failure in in function')

