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

This module contains several tables which are common to different lookup
functions.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Temperature Dependent data
---------------------------
.. autofunction:: chemicals.miscdata.lookup_VDI_tabular_data

"""

__all__ = ['lookup_VDI_tabular_data']

import os
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import register_df_source, data_source

### Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')

### CRC Handbook general tables
register_df_source(folder, 'Physical Constants of Inorganic Compounds.csv')
register_df_source(folder, 'Physical Constants of Organic Compounds.csv')

_VDI_dict_loaded = False
def _load_VDI_saturation_dict():
    """Read in a dict of assorted chemical properties at saturation for 58
    industrially important chemicals, from:

    Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2E. Berlin : Springer, 2010.
    This listing is the successor to that in:
    Schlunder, Ernst U, and International Center for Heat and Mass Transfer.
    Heat Exchanger Design Handbook. Washington: Hemisphere Pub. Corp., 1983.
    """
    import json
    global VDI_saturation_dict, _VDI_dict_loaded

    with open(os.path.join(folder, 'VDI Saturation Compounds Data.json')) as f:
        VDI_saturation_dict = json.loads(f.read())
    _VDI_dict_loaded = True

_CRC_data_loaded = False
def _load_CRC_data():
    global CRC_inorganic_data, CRC_organic_data, _CRC_data_loaded
    CRC_inorganic_data = data_source('Physical Constants of Inorganic Compounds.csv')
    CRC_organic_data = data_source('Physical Constants of Organic Compounds.csv')
    _CRC_data_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('CRC_inorganic_data', 'CRC_organic_data'):
            _load_CRC_data()
            return globals()[name]
        elif name == 'VDI_saturation_dict':
            _load_VDI_saturation_dict()
            return VDI_saturation_dict
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_CRC_data()
        _load_VDI_saturation_dict()

### VDI Saturation

# Created with the following code. Don't delete! Updates may be necessary.
#from chemicals.utils import to_num, rho_to_Vm
#import copy
#
#emptydict = {"Name": None, "MW": None, "Tc": None, "T": [], "P": [],
#             "Density (l)": [], "Density (g)": [], "Hvap": [], "Cp (l)": [],
#            "Cp (g)": [], "Mu (l)": [], "Mu (g)": [], "K (l)": [], "K (g)": [],
#            "Pr (l)": [], "Pr (g)": [], "sigma": [], "Beta": [],
#            "Volume (l)": [], "Volume (g)": []}
#VDI_saturation_dict = {}
#with open(os.path.join(folder, 'VDI Saturation Compounds Data.csv')) as f:
#    next(f)
#    for line in f:
#        values = to_num(line.strip('\n').split('\t'))
#        (CASRN, _name, _MW, _Tc, T, P, rhol, rhog, Hvap, cpl, cpg, mul, mug, kl, kg, prl, prg, sigma, Beta) = values
#        newdict = (VDI_saturation_dict[CASRN] if CASRN in VDI_saturation_dict else copy.deepcopy(emptydict))
#        newdict["Name"] = _name
#        newdict["MW"] = _MW
#        newdict["Tc"] = _Tc
#        newdict["T"].append(T)
#        newdict["P"].append(P)
#        newdict["Density (l)"].append(rhol)
#        newdict["Density (g)"].append(rhog)  # Not actually used
#        newdict["Hvap"].append(Hvap)
#        newdict["Cp (l)"].append(cpl)  # Molar
#        newdict["Cp (g)"].append(cpg)  # Molar
#        newdict["Mu (l)"].append(mul)
#        newdict["Mu (g)"].append(mug)
#        newdict["K (l)"].append(kl)
#        newdict["K (g)"].append(kg)
#        newdict["Pr (l)"].append(prl)
#        newdict["Pr (g)"].append(prl)
#        newdict["sigma"].append(sigma)
#        newdict["Beta"].append(Beta)
#        newdict["Volume (l)"].append(rho_to_Vm(rhol, _MW))
#        newdict["Volume (g)"].append(rho_to_Vm(rhog, _MW))
#        VDI_saturation_dict[CASRN] = newdict

def lookup_VDI_tabular_data(CASRN, prop):
    r'''This function retrieves the tabular data available for a given chemical
    and a given property. Lookup is based on CASRNs. Length of data returned
    varies between chemicals. All data is at saturation condition from [1]_.

    Function has data for 58 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]
    prop : string
        Property [-]

    Returns
    -------
    Ts : list
        Temperatures where property data is available, [K]
    props : list
        Properties at each temperature, [various]

    Notes
    -----
    The available properties are 'P', 'Density (l)', 'Density (g)', 'Hvap',
    'Cp (l)', 'Cp (g)', 'Mu (l)', 'Mu (g)', 'K (l)', 'K (g)', 'Pr (l)',
    'Pr (g)', 'sigma', 'Beta', 'Volume (l)', and 'Volume (g)'.

    Data is available for all properties and all chemicals; surface tension
    data was missing for mercury, but added as estimated from the a/b
    coefficients listed in Jasper (1972) to simplify the function.

    Examples
    --------
    >>> lookup_VDI_tabular_data('67-56-1', 'Mu (g)')
    ([337.63, 360.0, 385.0, 410.0, 435.0, 460.0, 500.0], [1.11e-05, 1.18e-05, 1.27e-05, 1.36e-05, 1.46e-05, 1.59e-05, 2.04e-05])
    >>> lookup_VDI_tabular_data('7782-41-4', 'sigma')
    ([53.49, 64.0, 74.0, 85.04, 92.0, 102.0, 112.0, 122.0, 132.0, 144.41], [0.0227, 0.02, 0.0166, 0.0136, 0.0117, 0.0092, 0.0068, 0.0045, 0.0024, 0.0])

    References
    ----------
    .. [1] Gesellschaft, VDI, ed. VDI Heat Atlas. 2E. Berlin : Springer, 2010.
    '''
    if not _VDI_dict_loaded: _load_VDI_saturation_dict()
    try:
        d = VDI_saturation_dict[CASRN]
    except KeyError:
        raise LookupError('CASRN not in VDI tabulation')
    try:
        props, Ts = d[prop], d['T']
    except:
        raise ValueError('property not specified correctly')
    Ts = [T for p, T in zip(props, Ts) if p]
    props = [p for p in props if p]
    # Not all data series converege to correct values
    if prop == 'sigma':
        if Ts[-1] < d['Tc']:
            Ts.append(d['Tc'])
            props.append(0.0)
        elif Ts[-1] == d['Tc']:
            props[-1] = 0.0
    return Ts, props
