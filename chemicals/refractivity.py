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

This module contains various refractive index lookup, calculation,
and unit conversion routines and dataframes.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:
    
Lookup Functions
----------------
.. autofunction:: chemicals.refractivity.RI
.. autofunction:: chemicals.refractivity.RI_methods
.. autodata:: chemicals.refractivity.RI_all_methods

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.refractivity.RI_IAPWS

Unit Conversions
----------------
.. autofunction:: chemicals.refractivity.brix_to_RI
.. autofunction:: chemicals.refractivity.RI_to_brix

Utility functions
-----------------
.. autofunction:: chemicals.refractivity.polarizability_from_RI
.. autofunction:: chemicals.refractivity.molar_refractivity_from_RI
.. autofunction:: chemicals.refractivity.RI_from_molar_refractivity

"""

__all__ = ['RI', 'RI_methods', 'RI_all_methods',
           'polarizability_from_RI', 'molar_refractivity_from_RI', 
           'RI_from_molar_refractivity', 'RI_IAPWS', 'RI_to_brix',
           'brix_to_RI']

from fluids.numerics import interp
from fluids.constants import pi, N_A
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)


# Register data sources and lazy load them

folder = os_path_join(source_path, 'Misc')
register_df_source(folder, 'CRC Handbook Organic RI.csv', 
                   csv_kwargs={'dtype': {'RI': float, 'RIT': float}})

CRC = 'CRC'

_RI_data_loaded = False
def _load_RI_data():
    global _RI_data_loaded, RI_data_CRC_organic, RI_sources
    RI_data_CRC_organic = data_source('CRC Handbook Organic RI.csv')
    RI_sources = {
        CRC: RI_data_CRC_organic,
    }

if PY37:
    def __getattr__(name):
        if name in ('RI_data_CRC_organic', 'RI_sources'):
            _load_RI_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_RI_data()

#  Refractive index functions

RI_all_methods = (CRC,)
'''Tuple of method name keys. See the `RI` for the actual references'''

def RI_methods(CASRN):
    """Return all methods available to obtain the RI for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the RI with the given
        inputs.

    See Also
    --------
    RI
    """
    if not _RI_data_loaded: _load_RI_data()
    return list_available_methods_from_df_dict(RI_sources, CASRN, 'RI')

def RI(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's refractive
    index. Lookup is based on CASRNs. Will automatically select a data source
    to use if no method is provided; returns None if the data is not available.

    Function has data for approximately 4500 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    RI : float
        Refractive Index on the Na D line, [-]
    T : float
        Temperature at which refractive index reading was made

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        RI_methods

    Notes
    -----
    Only one source is available in this function. It is:

        * 'CRC', a compillation of Organic RI data in [1]_.

    Examples
    --------
    >>> RI(CASRN='64-17-5')
    (1.3611, 293.15)

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    
    '''
    if not _RI_data_loaded: _load_RI_data()
    key = ('RI', 'RIT')
    if method:
        value = retrieve_from_df_dict(RI_sources, CASRN, key, method) 
    else:
        value = retrieve_any_from_df_dict(RI_sources, CASRN, key) 
    if value is None:
        value = (None, None)
    else:
        value = tuple(value)
    return value

def polarizability_from_RI(RI, Vm):
    r'''Returns the polarizability of a fluid given its molar volume and
    refractive index.

    .. math::
        \alpha = \left(\frac{3}{4\pi N_A}\right)
        \left(\frac{n^2-1}{n^2+2}\right)V_m

    Parameters
    ----------
    RI : float
        Refractive Index on Na D line, [-]
    Vm : float
        Molar volume of fluid, [m^3/mol]

    Returns
    -------
    alpha : float
        Polarizability [m^3]

    Notes
    -----
    This Lorentz-Lorentz-expression is most correct when van der Waals
    interactions dominate. Alternate conversions have been suggested.
    This is often expressed in units of cm^3 or Angstrom^3. To convert to these
    units, multiply by 1E9 or 1E30 respectively.

    Examples
    --------
    >>> polarizability_from_RI(1.3611, 5.8676E-5)
    5.147658206528923e-30

    References
    ----------
    .. [1] Panuganti, Sai R., Fei Wang, Walter G. Chapman, and Francisco M.
       Vargas. "A Simple Method for Estimation of Dielectric Constants and
       Polarizabilities of Nonpolar and Slightly Polar Hydrocarbons."
       International Journal of Thermophysics 37, no. 7 (June 6, 2016): 1-24.
       doi:10.1007/s10765-016-2075-8.
    '''
    return 3/(4*pi*N_A)*(RI**2-1)/(RI**2+2)*Vm

def molar_refractivity_from_RI(RI, Vm):
    r'''Returns the molar refractivity of a fluid given its molar volume and
    refractive index.

    .. math::
        R_m = \left(\frac{n^2-1}{n^2+2}\right)V_m

    Parameters
    ----------
    RI : float
        Refractive Index on Na D line, [-]
    Vm : float
        Molar volume of fluid, [m^3/mol]

    Returns
    -------
    Rm : float
        Molar refractivity [m^3/mol]

    Notes
    -----

    Examples
    --------
    >>> molar_refractivity_from_RI(1.3611, 5.8676E-5)
    1.2985217089649597e-05

    References
    ----------
    .. [1] Panuganti, Sai R., Fei Wang, Walter G. Chapman, and Francisco M.
       Vargas. "A Simple Method for Estimation of Dielectric Constants and
       Polarizabilities of Nonpolar and Slightly Polar Hydrocarbons."
       International Journal of Thermophysics 37, no. 7 (June 6, 2016): 1-24.
       doi:10.1007/s10765-016-2075-8.
    '''
    return (RI**2 - 1.)/(RI**2 + 2.)*Vm

def RI_from_molar_refractivity(Rm, Vm):
    r'''Returns the refractive index of a fluid given its molar volume and
    molar refractivity.

    .. math::
        RI = \sqrt{\frac{-2R_m - V_m}{R_m-V_m}}

    Parameters
    ----------
    Rm : float
        Molar refractivity [m^3/mol]
    Vm : float
        Molar volume of fluid, [m^3/mol]

    Returns
    -------
    RI : float
        Refractive Index on Na D line, [-]

    Notes
    -----

    Examples
    --------
    >>> RI_from_molar_refractivity(1.2985e-5, 5.8676E-5)
    1.3610932757685672

    References
    ----------
    .. [1] Panuganti, Sai R., Fei Wang, Walter G. Chapman, and Francisco M.
       Vargas. "A Simple Method for Estimation of Dielectric Constants and
       Polarizabilities of Nonpolar and Slightly Polar Hydrocarbons."
       International Journal of Thermophysics 37, no. 7 (June 6, 2016): 1-24.
       doi:10.1007/s10765-016-2075-8.
    '''
    Rm = ((-2*Rm - Vm)/(Rm-Vm))**0.5
    return Rm


def RI_IAPWS(T, rho, wavelength=0.5893):
    r'''Calculates the refractive index of water at a given temperature,
    density, and wavelength.

    .. math::
        n(\rho, T, \lambda) = \left(\frac{2A + 1}{1-A}\right)^{0.5}

    .. math::
        A(\delta, \theta, \Lambda) = \delta\left(a_0 + a_1\delta +
        a_2\theta + a_3\Lambda^2\theta + a_4\Lambda^{-2}
        \frac{a_5}{\Lambda^2-\Lambda_{UV}^2} + \frac{a_6}
        {\Lambda^2 - \Lambda_{IR}^2} + a_7\delta^2\right)

    .. math::
        \delta = \rho/(1000 \text{ kg/m}^3)

    .. math::
        \theta = T/273.15\text{K}

    .. math::
        \Lambda = \lambda/0.589 \mu m

    .. math::
        \Lambda_{IR} = 5.432937 

    .. math::
        \Lambda_{UV} = 0.229202

    Parameters
    ----------
    T : float
        Temperature of the water [K]
    rho : float
        Density of the water [kg/m^3]
    wavelength : float
        Wavelength of fluid [micrometers]

    Returns
    -------
    RI : float
        Refractive index of the water, [-]

    Notes
    -----
    This function is valid in the following range:
    261.15 K < T < 773.15 K
    0 < rho < 1060 kg/m^3
    0.2 < wavelength < 1.1 micrometers

    Test values are from IAPWS 2010 book.

    Examples
    --------
    >>> RI_IAPWS(298.15, 997.047435, 0.5893)
    1.3328581926471605

    References
    ----------
    .. [1] IAPWS, 1997. Release on the Refractive Index of Ordinary Water
       Substance as a Function of Wavelength, Temperature and Pressure.
    '''
    delta = rho*1e-3
    theta = T/273.15
    Lambda = wavelength/0.589

    LambdaIR = 5.432937
    LambdaUV = 0.229202
    
    Lambda2 = Lambda*Lambda

    A = delta*(0.244257733 + 0.0097463448*delta + -0.00373235*theta + 0.0002686785*Lambda2*theta + 
    0.0015892057/Lambda2 + 0.0024593426/(Lambda2 - LambdaUV*LambdaUV) + 
    0.90070492/(Lambda2 - LambdaIR*LambdaIR) - 0.0166626219*delta*delta)
    n = ((2*A + 1.)/(1. - A))**0.5
    return n

ICUMSA_1974_brix = list([float(i) for i in range(96)])
ICUMSA_1974_RIs = [1.33299, 1.33442, 1.33586, 1.33732, 1.33879, 1.34026, 1.34175, 
                   1.34325, 1.34477, 1.34629, 1.34782, 1.34937, 1.35093, 1.35250, 
                   1.35408, 1.35568, 1.35729, 1.35891, 1.36054, 1.36218, 1.36384, 
                   1.36551, 1.36720, 1.36889, 1.37060, 1.37233, 1.37406, 1.37582, 
                   1.37758, 1.37936, 1.38115, 1.38296, 1.38478, 1.38661, 1.38846, 
                   1.39032, 1.39220, 1.39409, 1.39600, 1.39792, 1.39986, 1.40181, 
                   1.40378, 1.40576, 1.40776, 1.40978, 1.41181, 1.41385, 1.41592, 
                   1.41799, 1.42009, 1.42220, 1.42432, 1.42647, 1.42863, 1.43080, 
                   1.43299, 1.43520, 1.43743, 1.43967, 1.44193, 1.44420, 1.44650, 
                   1.44881, 1.45113, 1.45348, 1.45584, 1.45822, 1.46061, 1.46303, 
                   1.46546, 1.46790, 1.47037, 1.47285, 1.47535, 1.47787, 1.48040, 
                   1.48295, 1.48552, 1.48811, 1.49071, 1.49333, 1.49597, 1.49862,
                   1.50129, 1.50398, 1.5067, 1.5094, 1.5122, 1.5149, 1.5177, 
                   1.5205, 1.5234, 1.5262, 1.5291, 1.5320]

def brix_to_RI(brix):
    """Convert a refractive index measurement on the `brix` scale to a standard
    refractive index.

    Parameters
    ----------
    brix : float
        Degrees brix to be converted, [°Bx]

    Returns
    -------
    RI : float
        Refractive index, [-]

    Notes
    -----
    The scale is officially defined from 0 to 85; but the data source contains
    values up to 95. Linear extrapolation outside of the bounds is performed;
    and a table of 96 values are linearly interpolated.
    
    The ICUMSA (International Committee of Uniform Method of Sugar Analysis)
    published a document setting out the reference values in 1974; but an 
    original data source has not been found and reviewed.

    Examples
    --------
    >>> brix_to_RI(5.8)
    1.341452
    >>> brix_to_RI(0.0)
    1.33299
    >>> brix_to_RI(95.0)
    1.532

    References
    ----------
    .. [1] "Refractometer　Data Book-Refractive Index and Brix | ATAGO CO.,
       LTD." Accessed June 13, 2020. 
       https://www.atago.net/en/databook-refractometer_relationship.php.
    """
    return interp(brix, ICUMSA_1974_brix, ICUMSA_1974_RIs, extrapolate=True)

def RI_to_brix(RI):
    """Convert a standard refractive index measurement to the `brix` scale.

    Parameters
    ----------
    RI : float
        Refractive index, [-]

    Returns
    -------
    brix : float
        Degrees brix to be converted, [°Bx]

    Notes
    -----
    The scale is officially defined from 0 to 85; but the data source contains
    values up to 95. 
    
    Linear extrapolation to values under 0 or above 95 is performed.
    
    The ICUMSA (International Committee of Uniform Method of Sugar Analysis)
    published a document setting out the reference values in 1974; but an 
    original data source has not been found and reviewed.

    Examples
    --------
    >>> RI_to_brix(1.341452)
    5.800000000000059
    >>> RI_to_brix(1.33299)
    0.0
    >>> RI_to_brix(1.532)
    95.0
    

    References
    ----------
    .. [1] "Refractometer　Data Book-Refractive Index and Brix | ATAGO CO.,
       LTD." Accessed June 13, 2020. 
       https://www.atago.net/en/databook-refractometer_relationship.php.
    """
    return interp(RI, ICUMSA_1974_RIs, ICUMSA_1974_brix, extrapolate=True)
