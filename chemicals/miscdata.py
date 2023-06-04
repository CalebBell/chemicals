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

from chemicals.data_reader import data_source, register_df_source
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

### Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')

### CRC Handbook general tables
register_df_source(folder, 'heos_constants.tsv')
register_df_source(folder, 'Physical Constants of Inorganic Compounds.csv')
register_df_source(folder, 'Physical Constants of Organic Compounds.csv')
register_df_source(folder, 'joback_predictions.tsv', int_CAS=True)
register_df_source(folder, 'wikidata_properties.tsv', sparsify=True, int_CAS=True)
register_df_source(folder, 'webbook_constants.tsv', sparsify=True, int_CAS=True)
register_df_source(folder, 'common_chemistry_data.tsv', sparsify=True, int_CAS=True)

# Common database IDs
JOBACK = 'JOBACK'
WIKIDATA = 'WIKIDATA'
WEBBOOK = 'WEBBOOK'
COMMON_CHEMISTRY = 'COMMON_CHEMISTRY'
PSI4_2022A = 'PSI4_2022A'
CHEMSEP = 'CHEMSEP'
JANAF = 'JANAF'
HEOS = 'HEOS'

# metadata about othe type of data being read in

EXPERIMENTAL = 'EXPERIMENTAL' # The quantity was measured directly or the standard means of measuring the quantity was performed
PREDICTED = 'PREDICTED' # A published or unpublished estimation method to guess the value
DEFINED = 'DEFINED' # e.g. enthalpy of formation of oxygen is 0
PROCESSED = 'PROCESSED' # Values that came from something that took in other values

EXPERIMENTAL_PRIMARY = 'EXPERIMENTAL_PRIMARY' # published the experimental data for the first time
EXPERIMENTAL_COMPILATION = 'EXPERIMENTAL_COMPILATION' # a paper that isn't publishing experimental data for the first time but putting multiple chemical values together
EXPERIMENTAL_REVIEW = 'EXPERIMENTAL_REVIEW' # A paper that looks at multiple single property values published in the literature for the same chemical and recommends one (potentially averaing them to get a new value)
EXPERIMENTAL_COMPILATION_SECONDARY = 'EXPERIMENTAL_COMPILATION_SECONDARY'# A paper that reproduces compilled data from another source

PREDICTED_GC = 'PREDICTED_GC'
PREDICTED_CSP = 'PREDICTED_CSP'
PREDICTED_MM = 'PREDICTED_MM' # molecular modeling
PREDICTED_QM = 'PREDICTED_QM' # quantum mechanics calc

# The definition of PROCESSED is that software was used, or a complex database (wikidata, common chemistry) involved

PROCESSED_EXPERIMENTAL = 'PROCESSED_EXPERIMENTAL' # Processed values that are all derived from experiments
PROCESSED_PREDICTED = 'PROCESSED_PREDICTED' # Processed values that are all derived from predictions
PROCESSED_EXPERIMENTAL_PREDICTED = 'PROCESSED_EXPERIMENTAL_PREDICTED' # Processed values that are derived experiments and or predictions may have no real data for a chemical

PROCESSED_EXPERIMENTAL_PREDICTED_SECONDARY = 'PROCESSED_EXPERIMENTAL_PREDICTED_SECONDARY' # same, but haven't seen original output only someone's republication
PROCESSED_EXPERIMENTAL_SECONDARY = 'PROCESSED_EXPERIMENTAL_SECONDARY' # same, but haven't seen original output only someone's republication
PROCESSED_PREDICTED_SECONDARY = 'PROCESSED_PREDICTED_SECONDARY' # same, but haven't seen original output only someone's republication


data_source_categories = [EXPERIMENTAL, PREDICTED, DEFINED, PROCESSED]
experimental_data_source_categories = [EXPERIMENTAL_PRIMARY, EXPERIMENTAL_COMPILATION, EXPERIMENTAL_REVIEW, EXPERIMENTAL_COMPILATION_SECONDARY]
predicted_data_source_categories = [PREDICTED_GC, PREDICTED_CSP, PREDICTED_MM, PREDICTED_QM]

processed_data_source_categories = [PROCESSED_EXPERIMENTAL, PROCESSED_PREDICTED, PROCESSED_EXPERIMENTAL_PREDICTED,
                                    PROCESSED_EXPERIMENTAL_PREDICTED_SECONDARY, PROCESSED_EXPERIMENTAL_SECONDARY, PROCESSED_PREDICTED_SECONDARY]


_VDI_dict_loaded = False
def _load_VDI_saturation_dict():
    """Read in a dict of assorted chemical properties at saturation for 58
    industrially important chemicals, from:

    Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2E. Berlin : Springer, 2010.
    This listing is the successor to that in:
    Schlunder, Ernst U, and International Center for Heat and Mass Transfer.
    Heat Exchanger Design Handbook. Washington: Hemisphere Pub. Corp., 1983.
    """
    import json
    global VDI_saturation_dict, _VDI_dict_loaded

    with open(os.path.join(folder, 'VDI Saturation Compounds Data.json')) as f:
        VDI_saturation_dict = json.loads(f.read())
    _VDI_dict_loaded = True

_miscdata_loaded = False
def _load_miscdata():
    global CRC_inorganic_data, CRC_organic_data, joback_predictions, wikidata_data, webbook_data, common_chemistry_data, heos_data, _miscdata_loaded
    heos_data = data_source('heos_constants.tsv')
    CRC_inorganic_data = data_source('Physical Constants of Inorganic Compounds.csv')
    CRC_organic_data = data_source('Physical Constants of Organic Compounds.csv')
    joback_predictions = data_source('joback_predictions.tsv')
    wikidata_data = data_source('wikidata_properties.tsv')
    webbook_data = data_source('webbook_constants.tsv')
    common_chemistry_data = data_source('common_chemistry_data.tsv')
    _miscdata_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('CRC_inorganic_data', 'CRC_organic_data', 'heos_data',
                    'joback_predictions', 'wikidata_data', 'webbook_data',
                    'common_chemistry_data'):
            _load_miscdata()
            return globals()[name]
        elif name == 'VDI_saturation_dict':
            _load_VDI_saturation_dict()
            return VDI_saturation_dict
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:
    if can_load_data:
        _load_miscdata()
        _load_VDI_saturation_dict()


@mark_numba_incompatible
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
    .. [1] Gesellschaft, VDI, ed. VDI Heat Atlas. 2E. Berlin : Springer, 2010.
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
    if prop == 'sigma' and props:
        if Ts[-1] < d['Tc']:
            Ts.append(d['Tc'])
            props.append(0.0)
        elif Ts[-1] == d['Tc']:
            props[-1] = 0.0
    return Ts, props
