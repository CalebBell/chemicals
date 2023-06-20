"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.

This module contains lookup functions for a chemical's triple temperature and
pressure. The triple temperature is the unique co-existence point between a
pure chemicals's solid, gas, and liquid state.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Triple Temperature
------------------
.. autofunction:: chemicals.triple.Tt
.. autofunction:: chemicals.triple.Tt_methods
.. autodata:: chemicals.triple.Tt_all_methods

Triple Pressure
---------------
.. autofunction:: chemicals.triple.Pt
.. autofunction:: chemicals.triple.Pt_methods
.. autodata:: chemicals.triple.Pt_all_methods

"""

__all__ = ['Tt_all_methods', 'Tt_methods', 'Tt',
           'Pt_all_methods', 'Pt_methods', 'Pt']

from chemicals import data_reader as dr
from chemicals import miscdata
from chemicals.data_reader import (
    data_source,
    database_constant_lookup,
    list_available_methods_from_df_dict,
    register_df_source,
    retrieve_any_from_df_dict,
    retrieve_from_df_dict,
)
from chemicals.phase_change import Tm
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

# Register data sources and lazy load them
folder = os_path_join(source_path, 'Triple Properties')
register_df_source(folder, 'Staveley 1981.tsv')

STAVELEY = 'STAVELEY'
MELTING = 'MELTING'

_triple_data_loaded = False
def _load_triple_data():
    global triple_data_Staveley, _triple_data_loaded, Tt_sources, Pt_sources
    triple_data_Staveley = data_source('Staveley 1981.tsv')
    _triple_data_loaded = True
    Tt_sources = {
        miscdata.HEOS: miscdata.heos_data,
        STAVELEY: triple_data_Staveley,
        miscdata.WEBBOOK: miscdata.webbook_data,
    }
    Pt_sources = Tt_sources.copy()

if PY37:
    def __getattr__(name):
        if name in ('triple_data_Staveley', 'Tt_sources', 'Pt_sources'):
            _load_triple_data()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:
    if can_load_data:
        _load_triple_data()

Tt_all_methods = (miscdata.HEOS, STAVELEY, miscdata.WEBBOOK, MELTING)
"""Tuple of method name keys. See the `Tt` for the actual references"""

@mark_numba_incompatible
def Tt_methods(CASRN):
    """Return all methods available to obtain the triple temperature for the
    desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Tt with the given
        inputs.

    See Also
    --------
    Tt
    """
    if not _triple_data_loaded: _load_triple_data()
    methods = list_available_methods_from_df_dict(Tt_sources, CASRN, 'Tt')
    if Tm(CASRN): methods.append(MELTING)
    return methods

@mark_numba_incompatible
def Tt(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's triple temperature.
    Lookup is based on CASRNs. Will automatically select a data source to use
    if no method is provided; returns None if the data is not available.

    Returns data from [1]_, [2]_ or [3]_, or a chemical's melting point if available.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Tt : float
        Triple point temperature, [K].

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `Tt_all_methods`.

    Notes
    -----
    Median difference between melting points and triple points is 0.02 K.
    Accordingly, this should be more than good enough for engineering
    applications.

    The data in [1]_ is originally on the ITS-68 temperature scale, but was
    converted to ITS-90. The numbers were rounded to 6 decimal places
    arbitrarily and the conversion was performed with this library.



    Examples
    --------
    Ammonia

    >>> Tt('7664-41-7')
    195.49

    See Also
    --------
    Tt_methods

    References
    ----------
    .. [1] Staveley, L. A. K., L. Q. Lobo, and J. C. G. Calado. "Triple-Points
       of Low Melting Substances and Their Use in Cryogenic Work." Cryogenics
       21, no. 3 (March 1981): 131-144. doi:10.1016/0011-2275(81)90264-2.
    .. [2] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    .. [3] Huber, Marcia L., Eric W. Lemmon, Ian H. Bell, and Mark O. McLinden.
       "The NIST REFPROP Database for Highly Accurate Properties of Industrially
       Important Fluids." Industrial & Engineering Chemistry Research 61, no. 42
       (October 26, 2022): 15449-72. https://doi.org/10.1021/acs.iecr.2c01427.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'Tt')
        if found: return val
    if not _triple_data_loaded: _load_triple_data()
    if method:
        if method == MELTING:
            return Tm(CASRN)
        else:
            return retrieve_from_df_dict(Tt_sources, CASRN, 'Tt', method)
    else:
        Tt = retrieve_any_from_df_dict(Tt_sources, CASRN, 'Tt')
        if Tt: return Tt
        return Tm(CASRN)

Pt_all_methods = (miscdata.HEOS, STAVELEY, miscdata.WEBBOOK)
"""Tuple of method name keys. See the `Pt` for the actual references"""

@mark_numba_incompatible
def Pt_methods(CASRN):
    """Return all methods available to obtain the Pt for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Pt with the given
        inputs.

    See Also
    --------
    Pt
    """
    if not _triple_data_loaded: _load_triple_data()
    return list_available_methods_from_df_dict(Pt_sources, CASRN, 'Pt')

@mark_numba_incompatible
def Pt(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's triple pressure.
    Lookup is based on CASRNs. Will automatically select a data source to use
    if no method is provided; returns None if the data is not available.

    Returns data from [1]_, [2]_, or [3]_.

    This function does not implement it but it is also possible to calculate
    the vapor pressure at the triple temperature from a vapor pressure
    correlation, if data is available; note most Antoine-type correlations do
    not extrapolate well to this low of a pressure.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Pt : float
        Triple point pressure, [Pa]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `Pt_all_methods`.

    Notes
    -----

    Examples
    --------
    Ammonia

    >>> Pt('7664-41-7')
    6053.386

    See Also
    --------
    Pt_methods

    References
    ----------
    .. [1] Staveley, L. A. K., L. Q. Lobo, and J. C. G. Calado. "Triple-Points
       of Low Melting Substances and Their Use in Cryogenic Work." Cryogenics
       21, no. 3 (March 1981): 131-144. doi:10.1016/0011-2275(81)90264-2.
    .. [2] Shen, V.K., Siderius, D.W., Krekelberg, W.P., and Hatch, H.W., Eds.,
       NIST WebBook, NIST, http://doi.org/10.18434/T4M88Q
    .. [3] Huber, Marcia L., Eric W. Lemmon, Ian H. Bell, and Mark O. McLinden.
       "The NIST REFPROP Database for Highly Accurate Properties of Industrially
       Important Fluids." Industrial & Engineering Chemistry Research 61, no. 42
       (October 26, 2022): 15449-72. https://doi.org/10.1021/acs.iecr.2c01427.
    '''
    if dr.USE_CONSTANTS_DATABASE and method is None:
        val, found = database_constant_lookup(CASRN, 'Pt')
        if found: return val
    if not _triple_data_loaded: _load_triple_data()
    if method:
        return retrieve_from_df_dict(Pt_sources, CASRN, 'Pt', method)
    else:
        return retrieve_any_from_df_dict(Pt_sources, CASRN, 'Pt')

