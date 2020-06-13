# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

__all__ = ['RI_methods', 'refractive_index', 
           'polarizability_from_RI', 'molar_refractivity_from_RI', 
           'RI_from_molar_refractivity', 'RI_IAPWS']

import os
from fluids.constants import pi, N_A
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)
from chemicals.utils import PY37


# %% Register data sources and lazy load them

folder = os.path.join(os.path.dirname(__file__), 'Misc')
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
    _load_RI_data()

# %% Refractive index functions

RI_methods = [CRC]

def refractive_index(CASRN, get_methods=False, method=None,
                     full_info=True):
    r'''This function handles the retrieval of a chemical's refractive
    index. Lookup is based on CASRNs. Will automatically select a data source
    to use if no method is provided; returns None if the data is not available.

    Function has data for approximately 4500 chemicals.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    RI : float
        Refractive Index on the Na D line, [-]
    T : float, only returned if full_info == True
        Temperature at which refractive index reading was made
    methods : list, only returned if get_methods == True
        List of methods which can be used to obtain RI with the given inputs

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by constants in
        RI_methods
    get_methods : bool, optional
        If True, function will determine which methods can be used to obtain
        RI for the desired chemical, and will return methods instead of RI
    full_info : bool, optional
        If True, function will return the temperature at which the refractive
        index reading was made

    Notes
    -----
    Only one source is available in this function. It is:

        * 'CRC', a compillation of Organic RI data in [1]_.

    Examples
    --------
    >>> refractive_index(CASRN='64-17-5')
    (1.3611, 293.15)

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    if not _RI_data_loaded: _load_RI_data()
    if get_methods:
        return list_available_methods_from_df_dict(RI_sources, CASRN, 'RI')
    else:
        key = ('RI', 'RIT') if full_info else 'RI'
        if method:
            value = retrieve_from_df_dict(RI_sources, CASRN, key, method) 
        else:
            value = retrieve_any_from_df_dict(RI_sources, CASRN, key) 
        if full_info:
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

        n(\rho, T, \lambda) = \left(\frac{2A + 1}{1-A}\right)^{0.5}\\

        A(\delta, \theta, \Lambda) = \delta\left(a_0 + a_1\delta +
        a_2\theta + a_3\Lambda^2\theta + a_4\Lambda^{-2}
        \frac{a_5}{\Lambda^2-\Lambda_{UV}^2} + \frac{a_6}
        {\Lambda^2 - \Lambda_{IR}^2} + a_7\delta^2\right)

        \delta = \rho/(1000 \text{ kg/m}^3)\\
        \theta = T/273.15\text{K}\\
        \Lambda = \lambda/0.589 \mu m

        \Lambda_{IR} = 5.432937 \\
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
