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

This module contains various vapor pressure estimation routines, dataframes
of fit coefficients, some compound-specific equations, some analytical fitting
routines, and sublimation pressure routines.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:


Fit Correlations
----------------
.. autofunction:: chemicals.vapor_pressure.Antoine
.. autofunction:: chemicals.vapor_pressure.Wagner
.. autofunction:: chemicals.vapor_pressure.Wagner_original
.. autofunction:: chemicals.vapor_pressure.TRC_Antoine_extended

Vapor Pressure Estimation Correlations
--------------------------------------
.. autofunction:: chemicals.vapor_pressure.Lee_Kesler
.. autofunction:: chemicals.vapor_pressure.Ambrose_Walton
.. autofunction:: chemicals.vapor_pressure.boiling_critical_relation
.. autofunction:: chemicals.vapor_pressure.Sanjari
.. autofunction:: chemicals.vapor_pressure.Edalat

Sublimation Pressure Estimation Correlations
--------------------------------------------
.. autofunction:: chemicals.vapor_pressure.Psub_Clapeyron

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.vapor_pressure.Psat_IAPWS
.. autofunction:: chemicals.vapor_pressure.dPsat_IAPWS_dT
.. autofunction:: chemicals.vapor_pressure.Tsat_IAPWS


Analytical Fit Equations
------------------------
.. autofunction:: chemicals.vapor_pressure.Antoine_coeffs_from_point
.. autofunction:: chemicals.vapor_pressure.Antoine_AB_coeffs_from_point
.. autofunction:: chemicals.vapor_pressure.DIPPR101_ABC_coeffs_from_point


Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: Psat_data_WagnerMcGarry

    Coefficients for the Wagner 3,6 original model equation documented in
    :obj:`Wagner_original` with data for 245 chemicals, from [1]_.

.. data:: Psat_data_WagnerPoling

    Coefficients for the Wagner 2.5, 5 model equation documented in
    :obj:`Wagner` in [2]_, with data for 104 chemicals.

.. data:: Psat_data_AntoinePoling

    Standard Antoine equation coefficients, as documented in the function
    :obj:`Antoine` and with data for 325 fluids from [2]_.
    Coefficients were altered to be in units of Pa and Celcius.

.. data:: Psat_data_AntoineExtended

    Data for 97 chemicals in [2]_ for the TRC extended Antoine model
    :obj:`TRC_Antoine_extended`.

.. data:: Psat_data_Perrys2_8

    A collection of 341 coefficient sets for :obj:`thermo.dippr.EQ101` from
    the DIPPR database published openly in [4]_.

.. data:: Psat_data_VDI_PPDS_3

    Coefficients for the Wagner equation :obj:`Wagner`, published
    openly in [3]_.

.. [1] McGarry, Jack. "Correlation and Prediction of the Vapor Pressures of
    Pure Liquids over Large Pressure Ranges." Industrial & Engineering
    Chemistry Process Design and Development 22, no. 2 (April 1, 1983):
    313-22. doi:10.1021/i200021a023.
.. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
    New York: McGraw-Hill Professional, 2000.
.. [3] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
    Berlin; New York:: Springer, 2010.
.. [4] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    Eighth Edition. McGraw-Hill Professional, 2007.

The structure of each dataframe is shown below:

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.vapor_pressure.Psat_data_WagnerMcGarry

    In [3]: chemicals.vapor_pressure.Psat_data_WagnerPoling

    In [4]: chemicals.vapor_pressure.Psat_data_AntoinePoling

    In [5]: chemicals.vapor_pressure.Psat_data_AntoineExtended

    In [6]: chemicals.vapor_pressure.Psat_data_Perrys2_8

    In [7]: chemicals.vapor_pressure.Psat_data_VDI_PPDS_3

"""

from __future__ import division

__all__ = ['Antoine', 'Wagner_original', 'Wagner', 'TRC_Antoine_extended',
           'boiling_critical_relation', 'Lee_Kesler', 'Ambrose_Walton',
           'Edalat', 'Sanjari', 'Psat_IAPWS', 'dPsat_IAPWS_dT', 'Tsat_IAPWS',
           'Psub_Clapeyron',
           'Antoine_coeffs_from_point', 'Antoine_AB_coeffs_from_point',
           'DIPPR101_ABC_coeffs_from_point']

import os
from fluids.constants import R
from fluids.numerics import numpy as np
from math import e
from chemicals.utils import log, exp, sqrt, isnan
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.dippr import EQ101
from chemicals import miscdata
from chemicals.data_reader import register_df_source, data_source

folder = os_path_join(source_path, 'Vapor Pressure')

register_df_source(folder, 'Antoine Collection Poling.tsv')
register_df_source(folder, 'Table 2-8 Vapor Pressure of Inorganic and Organic Liquids.tsv')

register_df_source(folder, 'Wagner Original McGarry.tsv', csv_kwargs={
        'dtype': {'A': float, 'B': float, 'C': float, 'D': float,
                 'Pc': float, 'Tc': float, 'Tmin': float}})

register_df_source(folder, 'Wagner Collection Poling.tsv', csv_kwargs={
        'dtype': {'A': float, 'B': float, 'C': float, 'D': float, 'Pc': float,
                  'Tc': float, 'Pc': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Antoine Extended Collection Poling.tsv', csv_kwargs={
    'dtype':{'A': float, 'B': float, 'C': float, 'Tc': float, 'to': float,
             'n': float, 'E': float, 'F': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'VDI PPDS Boiling temperatures at different pressures.tsv', csv_kwargs={
        'dtype':{'Tm': float, 'Tc': float, 'Pc': float, 'A': float,
                 'B': float, 'C': float, 'D': float}})

_vapor_pressure_dfs_loaded = False
def load_vapor_pressure_dfs():
    global Psat_data_WagnerMcGarry, Psat_values_WagnerMcGarry, Psat_data_AntoinePoling, Psat_values_AntoinePoling
    global Psat_data_WagnerPoling, Psat_values_WagnerPoling, Psat_data_AntoineExtended, Psat_values_AntoineExtended
    global Psat_data_Perrys2_8, Psat_values_Perrys2_8, Psat_data_VDI_PPDS_3, Psat_values_VDI_PPDS_3
    global _vapor_pressure_dfs_loaded

    # 57463 bytes for df; 13720 bytes for numpy
    Psat_data_WagnerMcGarry = data_source('Wagner Original McGarry.tsv')
    Psat_values_WagnerMcGarry = np.array(Psat_data_WagnerMcGarry.values[:, 1:], dtype=float)

    # 58216 bytes for df; 13000 bytes for numpy
    Psat_data_AntoinePoling = data_source('Antoine Collection Poling.tsv')
    Psat_values_AntoinePoling = np.array(Psat_data_AntoinePoling.values[:, 1:], dtype=float)

    # 20928 bytes for df; 7488 bytes for numpy
    Psat_data_WagnerPoling = data_source('Wagner Collection Poling.tsv')
    Psat_values_WagnerPoling = np.array(Psat_data_WagnerPoling.values[:, 1:], dtype=float)

    # 21388 bytes for df; 7760 bytes for numpy
    Psat_data_AntoineExtended = data_source('Antoine Extended Collection Poling.tsv')
    Psat_values_AntoineExtended = np.array(Psat_data_AntoineExtended.values[:, 1:], dtype=float)

    # 65740 bytes for df; 21760 bytes for numpy
    Psat_data_Perrys2_8 = data_source('Table 2-8 Vapor Pressure of Inorganic and Organic Liquids.tsv')
    Psat_values_Perrys2_8 = np.array(Psat_data_Perrys2_8.values[:, 1:], dtype=float)

    # 52742 bytes for df; 15400 bytes for numpy
    Psat_data_VDI_PPDS_3 = data_source('VDI PPDS Boiling temperatures at different pressures.tsv')
    Psat_values_VDI_PPDS_3 = np.array(Psat_data_VDI_PPDS_3.values[:, 1:], dtype=float)

    _vapor_pressure_dfs_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('Psat_data_WagnerMcGarry', 'Psat_values_WagnerMcGarry', 'Psat_data_AntoinePoling',
                    'Psat_values_AntoinePoling', 'Psat_data_WagnerPoling', 'Psat_values_WagnerPoling',
                    'Psat_data_AntoineExtended', 'Psat_values_AntoineExtended', 'Psat_data_Perrys2_8',
                    'Psat_values_Perrys2_8', 'Psat_data_VDI_PPDS_3', 'Psat_values_VDI_PPDS_3'):
            load_vapor_pressure_dfs()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        load_vapor_pressure_dfs()




def Antoine(T, A, B, C, base=10.0):
    r'''Calculates vapor pressure of a chemical using the Antoine equation.
    Parameters `A`, `B`, and `C` are chemical-dependent. Parameters can be
    found in numerous sources; however units of the coefficients used vary.
    Originally proposed by Antoine (1888) [2]_.

    .. math::
        \log_{\text{base}} P^{\text{sat}} = A - \frac{B}{T+C}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    base : float, optional
        Optional base of logarithm; 10 by default

    Returns
    -------
    Psat : float
        Vapor pressure calculated with coefficients [Pa]

    Notes
    -----
    Assumes coefficients are for calculating vapor pressure in Pascal.
    Coefficients should be consistent with input temperatures in Kelvin;
    however, if both the given temperature and units are specific to degrees
    Celcius, the result will still be correct.

    **Converting units in input coefficients:**

        * **ln to log10**: Divide A and B by ln(10)=2.302585 to change
          parameters for a ln equation to a log10 equation.
        * **log10 to ln**: Multiply A and B by ln(10)=2.302585 to change
          parameters for a log equation to a ln equation.
        * **mmHg to Pa**: Add log10(101325/760)= 2.1249 to A.
        * **kPa to Pa**: Add log_{base}(1000)= 6.908 to A for log(base)
        * **bar to Pa**: Add log_{base}(100000)= 11.5129254 to A for log(base)
        * **°C to K**: Subtract 273.15 from C only!

    Examples
    --------
    Methane, coefficients from [1]_, at 100 K:

    >>> Antoine(100.0, 8.7687, 395.744, -6.469)
    34478.367349639906

    Tetrafluoromethane, coefficients from [1]_, at 180 K

    >>> Antoine(180, A=8.95894, B=510.595, C=-15.95)
    702271.0518579542

    Oxygen at 94.91 K, with coefficients from [3]_ in units of °C, mmHg, log10,
    showing the conversion of coefficients A (mmHg to Pa) and C (°C to K)

    >>> Antoine(94.91, 6.83706+2.1249, 339.2095, 268.70-273.15)
    162978.88655572367

    n-hexane with Antoine coefficients from the NIST webbook in units of K and
    bar, calculating the vapor pressure in Pa at 200 K:

    >>> Antoine(T=200, A=3.45604+5, B=1044.038, C=-53.893)
    20.4329803671

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Antoine, C. 1888. Tensions des Vapeurs: Nouvelle Relation Entre les
       Tensions et les Tempé. Compt.Rend. 107:681-684.
    .. [3] Yaws, Carl L. The Yaws Handbook of Vapor Pressure: Antoine
       Coefficients. 1 edition. Houston, Tex: Gulf Publishing Company, 2007.
    '''
    return base**(A - B/(T + C))


def Antoine_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2, base=10.0):
    r'''Calculates the antoine coefficients `A`, `B`, and `C` from a known
    vapor pressure and its first and second temperature derivative.

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Psat : float
        Vapor pressure at specified `T` [Pa]
    dPsat_dT : float
        First temperature derivative of vapor pressure at specified `T` [Pa/K]
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure at specified `T` [Pa/K^2]
    Base : float, optional
        Base of logarithm; 10 by default

    Returns
    -------
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]

    Notes
    -----
    Coefficients are for calculating vapor pressure in Pascal. This is
    primarily useful for interconverting vapor pressure models, not fitting
    experimental data.

    Derived with SymPy as follows:

    >>> from sympy import * # doctest: +SKIP
    >>> base, A, B, C, T = symbols('base, A, B, C, T') # doctest: +SKIP
    >>> v = base**(A - B/(T + C)) # doctest: +SKIP
    >>> d1, d2 = diff(v, T), diff(v, T, 2) # doctest: +SKIP
    >>> vk, d1k, d2k = symbols('vk, d1k, d2k') # doctest: +SKIP
    >>> solve([Eq(v, vk), Eq(d1, d1k), Eq(d2, d2k)], [A, B, C]) # doctest: +SKIP

    Examples
    --------
    Recalculate some coefficients from a calcualted value and its derivative:


    >>> T = 178.01
    >>> A, B, C = (24.0989474955895, 4346.793091137991, -18.96968471040141)
    >>> Psat = Antoine(T, A, B, C, base=exp(1))
    >>> dPsat_dT, d2Psat_dT2 = (0.006781441203850251, 0.0010801244983894853) # precomputed
    >>> Antoine_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2, base=exp(1))
    (24.098947495155, 4346.793090994, -18.969684713118)

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    x0 = 1.0/log(base)
    x1 = Psat*d2Psat_dT2
    dPsat_dT_2 = dPsat_dT*dPsat_dT
    x3 = 1.0/(x1 - dPsat_dT_2)
    x4 = dPsat_dT_2 + dPsat_dT_2
    A = x0*log(Psat*exp(-x3*x4))
    B = 4.0*Psat*dPsat_dT*dPsat_dT_2*x0/(Psat*Psat*d2Psat_dT2*d2Psat_dT2 + dPsat_dT_2*dPsat_dT_2 - x1*x4)
    C = -x3*(2.0*Psat*dPsat_dT + T*x1 - T*dPsat_dT_2)
    return (A, B, C)

def Antoine_AB_coeffs_from_point(T, Psat, dPsat_dT, base=10.0):
    r'''Calculates the antoine coefficients `A`, `B`, with `C` set to zero to
    improve low-temperature or high-temperature extrapolation, from a known
    vapor pressure and its first temperature derivative.

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Psat : float
        Vapor pressure at specified `T` [Pa]
    dPsat_dT : float
        First temperature derivative of vapor pressure at specified `T` [Pa/K]
    Base : float, optional
        Base of logarithm; 10 by default

    Returns
    -------
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]

    Notes
    -----
    Coefficients are for calculating vapor pressure in Pascal. This is
    primarily useful for interconverting vapor pressure models, not fitting
    experimental data.

    Derived with SymPy as follows:

    >>> from sympy import * # doctest: +SKIP
    >>> base, A, B, T = symbols('base, A, B, T') # doctest: +SKIP
    >>> v = base**(A - B/T) # doctest: +SKIP
    >>> d1, d2 = diff(v, T), diff(v, T, 2) # doctest: +SKIP
    >>> vk, d1k = symbols('vk, d1k') # doctest: +SKIP
    >>> solve([Eq(v, vk), Eq(d1, d1k)], [A, B]) # doctest: +SKIP

    Examples
    --------
    Recalculate some coefficients from a calcualted value and its derivative:

    >>> T = 178.01
    >>> A, B = (27.358925161569008, 5445.569591293226)
    >>> Psat = Antoine(T, A, B, C=0, base=exp(1))
    >>> dPsat_dT = B*exp(1)**(A - B/T)*log(exp(1))/T**2
    >>> Antoine_AB_coeffs_from_point(T, Psat, dPsat_dT, base=exp(1))
    (27.35892516156901, 5445.569591293226)

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    log_base_inv = 1.0/log(base)
    Psat_inv = 1.0/Psat
    A = log(Psat*exp(T*dPsat_dT*Psat_inv))*log_base_inv
    B = T*T*dPsat_dT*log_base_inv*Psat_inv
    return (A, B)

def DIPPR101_ABC_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2):
    r'''Calculates the first three DIPPR101 coefficients `A`, `B`, and `C`
    from a known vapor pressure and its first and second temperature derivative.

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Psat : float
        Vapor pressure at specified `T` [Pa]
    dPsat_dT : float
        First temperature derivative of vapor pressure at specified `T` [Pa/K]
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure at specified `T` [Pa/K^2]

    Returns
    -------
    A : float
        DIPPR101 `A` parameter (same as Antoine `A`), [-]
    B : float
        DIPPR101 `B` parameter (same as Antoine `B`), [K]
    C: float
        DIPPR101 `C` parameter (NOT same as Antoine `C`, multiplied by log(T)),
        [-]

    Notes
    -----
    Coefficients are for calculating vapor pressure in Pascal. This is
    primarily useful for interconverting vapor pressure models, not fitting
    experimental data.

    Derived with SymPy as follows:

    >>> from sympy import * # doctest: +SKIP
    >>> base, A, B, C, T = symbols('base, A, B, C, T') # doctest: +SKIP
    >>> v = exp(A - B/T + C*log(T)) # doctest: +SKIP
    >>> d1, d2 = diff(v, T), diff(v, T, 2) # doctest: +SKIP
    >>> vk, d1k, d2k = symbols('vk, d1k, d2k') # doctest: +SKIP
    >>> solve([Eq(v, vk), Eq(d1, d1k), Eq(d2, d2k)], [A, B, C]) # doctest: +SKIP

    Examples
    --------
    Calculate the coefficients:

    >>> T = 178.01
    >>> Psat, dPsat_dT, d2Psat_dT2 = (0.03946094565666715, 0.006781441203850251, 0.0010801244983894853)
    >>> DIPPR101_ABC_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2)
    (72.47169926642, -6744.620564969, -7.2976291987890)
    '''
    x0 = Psat*Psat
    x1 = 1.0/x0
    x2 = Psat*dPsat_dT
    x3 = -T*dPsat_dT*dPsat_dT
    x4 = T*d2Psat_dT2
    x5 = T*(Psat*(2.0*dPsat_dT + x4) + x3)
    A = x1*(T*x2 + x0*log(Psat) - x5*log(T) - x5)
    B = T*T*x1*(Psat*x4 + x2 + x3)
    C = x1*x5
    return (A, B, C)

def TRC_Antoine_extended(T, Tc, to, A, B, C, n, E, F):
    r'''Calculates vapor pressure of a chemical using the TRC Extended Antoine
    equation. Parameters are chemical dependent, and said to be from the
    Thermodynamics Research Center (TRC) at Texas A&M. Coefficients for various
    chemicals can be found in [1]_.

    .. math::
        \log_{10} P^{sat} = A - \frac{B}{T + C} + 0.43429x^n + Ex^8 + Fx^{12}

    .. math::
        x = \max \left(\frac{T-t_o-273.15}{T_c}, 0 \right)

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    A, B, C, n, E, F : floats
        Regressed coefficients for the Antoine Extended (TRC) equation,
        specific for each chemical, [-]

    Returns
    -------
    Psat : float
        Vapor pressure calculated with coefficients [Pa]

    Notes
    -----
    Assumes coefficients are for calculating vapor pressure in Pascal.
    Coefficients should be consistent with input temperatures in Kelvin;

    Examples
    --------
    Tetrafluoromethane, coefficients from [1]_, at 180 K:

    >>> TRC_Antoine_extended(180.0, 227.51, -120., 8.95894, 510.595, -15.95,
    ... 2.41377, -93.74, 7425.9)
    706317.0898414153

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    x = (T - to - 273.15)/Tc
    if x < 0.0:
        x = 0.0
    x4 = x*x*x*x
    return 10.**(A - B/(T+C) + 0.43429*x**n + x4*x4*(E + F*x4))


def Wagner_original(T, Tc, Pc, a, b, c, d):
    r'''Calculates vapor pressure using the Wagner equation (3, 6 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \ln P^{sat}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^3 + d\tau^6}
        {T_r}

    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a, b, c, d : floats
        Parameters for wagner equation. Specific to each chemical. [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Warning: Pc is often treated as adjustable constant.

    Examples
    --------
    Methane, coefficients from [2]_, at 100 K.

    >>> Wagner_original(100.0, 190.53, 4596420., a=-6.00435, b=1.1885,
    ... c=-0.834082, d=-1.22833)
    34520.44601450499

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] McGarry, Jack. "Correlation and Prediction of the Vapor Pressures of
       Pure Liquids over Large Pressure Ranges." Industrial & Engineering
       Chemistry Process Design and Development 22, no. 2 (April 1, 1983):
       313-22. doi:10.1021/i200021a023.
    '''
    Tr = T/Tc
    tau = 1.0 - Tr
    tau2 = tau*tau
    tau_Tr = tau/Tr
    return Pc*exp(((d*tau2*tau + c)*tau2 + a + b*sqrt(tau))*tau_Tr)


def Wagner(T, Tc, Pc, a, b, c, d):
    r'''Calculates vapor pressure using the Wagner equation (2.5, 5 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \ln P^{sat}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^{2.5}
        + d\tau^5} {T_r}

    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a, b, c, d : floats
        Parameters for wagner equation. Specific to each chemical. [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Warning: Pc is often treated as adjustable constant.

    Examples
    --------
    Methane, coefficients from [2]_, at 100 K.

    >>> Wagner(100., 190.551, 4599200, -6.02242, 1.26652, -0.5707, -1.366)
    34415.00476263708

    References
    ----------
    .. [1] Wagner, W. "New Vapour Pressure Measurements for Argon and Nitrogen and
       a New Method for Establishing Rational Vapour Pressure Equations."
       Cryogenics 13, no. 8 (August 1973): 470-82. doi:10.1016/0011-2275(73)90003-9
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    tau = 1.0 - T/Tc
    return Pc*exp((a*tau + b*tau**1.5 + c*tau**2.5 + d*tau**5)/Tr)


def Psat_IAPWS(T):
    r'''Calculates vapor pressure of water using the IAPWS explicit equation.

    .. math::
        P^{sat} = 10^6 \left[ \frac{2C}{-B + \sqrt{B^2 - 4AC}}  \right]^4

    .. math::
        A = \nu^2 + n_1 \nu + n_2

    .. math::
        B = n_3 \nu^2 + n_4\nu + n_5

    .. math::
        C = n_6\nu^2 + n_7\nu + n_8

    .. math::
        \nu = T + \frac{n_9}{T - n_{10}}

    Parameters
    ----------
    T : float
        Temperature of water, [K]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    This formulation is quite efficient, and can also be solved backward.
    The range of validity of this equation is 273.15 K < T < 647.096 K, the
    IAPWS critical point.

    Extrapolation to lower temperatures is very poor. The function continues to
    decrease until a pressure of 5.7 mPa is reached at 159.77353993926621 K;
    under that pressure the vapor pressure increases, which is obviously wrong.

    Examples
    --------
    >>> Psat_IAPWS(300.)
    3536.58941301301

    References
    ----------
    .. [1] Kretzschmar, Hans-Joachim, and Wolfgang Wagner. International Steam
       Tables: Properties of Water and Steam Based on the Industrial
       Formulation IAPWS-IF97. Springer, 2019.
    '''
    v = T - 0.23855557567849/(T - 0.65017534844798E3)
    v2 = v*v
    A = v2 + 0.11670521452767E4*v - 0.72421316703206E6
    B = -0.17073846940092E2*v2 + 0.12020824702470E5*v - 0.32325550322333E7
    C = 0.14915108613530E2*v2 - 0.48232657361591E4*v + 0.40511340542057E6
    x = ((C + C)/(sqrt(B*B - 4.0*A*C) - B))
    x2 = x*x
    P = 1E6*x2*x2
    return P


def dPsat_IAPWS_dT(T):
    r'''Calculates the first temperature dervative of vapor pressure of water
    using the IAPWS explicit equation. This was derived with SymPy, using the
    CSE method.

    Parameters
    ----------
    T : float
        Temperature of water, [K]

    Returns
    -------
    dPsat_dT : float
        Temperature dervative of vapor pressure at T [Pa/K]

    Notes
    -----
    The derivative of this is useful when solving for water dew point.

    Examples
    --------
    >>> dPsat_IAPWS_dT(300.)
    207.88388134164282

    References
    ----------
    .. [1] Kretzschmar, Hans-Joachim, and Wolfgang Wagner. International Steam
       Tables: Properties of Water and Steam Based on the Industrial
       Formulation IAPWS-IF97. Springer, 2019.
    '''
    # 1 power, four divisions
    x0 = 1.0/(T - 650.175348447980014)
    x1 = T - 0.238555575678490006*x0
    x2 = x1*x1
    x3 = -0.0119059642846224972*T + 0.00284023416392566131*x0 + 0.0000368171193891888733*x2 + 1.0
    x4 = 0.00371867596455583913*T
    x5 = 0.000887110885486382334*x0
    x6 = 5.28184261979789455e-6*x2
    x7 = x4 - x5 - x6 - 1.0
    x8 = 4668.20858110680001*T - 1113.62718545319967*x0 + 4.0*x2 - 2896852.66812823992
    x9 = sqrt(x7*x7 - 9.56991643658934639e-14*x8*(-4823.26573615910002*T + 1150.61693433977007*x0 + 14.9151086135300002*x2 + 405113.405420569994))
    x10 = -x4 + x5 + x6 + x9 + 1.0
    x10_inv = 1.0/x10
    x11 = 0.00153804662447919488*T - 1.0
    x11 = 1.0/(x11*x11)
    x12 = x1*(1.12864813714895499e-6*x11 + 2.0)
    x20 = x10_inv*x10_inv*x10_inv*x10_inv
    return (-3946.77829948379076*x3*x3*x3*(2.68752888216023447e-8*x11 - 0.000147268477556755493*x12
                + 0.0476238571384899889 + x3*(-8.39415340011308276e-9*x11 + 0.0000211273704791915782*x12
                - 0.0148747038582233565 + 2.0*(x7*(4.19707670005654138e-9*x11
            - 0.0000105636852395957891*x12 + 0.00743735192911167825) + x8*(2.60482114645230015e-16*x11
            - 1.42736343074136086e-12*x12 + 4.6158250046507185e-10) + (0.00263438245944447809*x11
            + 4.0*x12 + 4668.20858110680001)*(4.6158250046507185e-10*T - 1.10113079121562103e-10*x0
            - 1.42736343074136086e-12*x2 - 3.8769014372169964e-8))/x9)*x10_inv)*x20)

def Tsat_IAPWS(P):
    r'''Calculates the saturation temperature of water using the IAPWS explicit
    equation.

    .. math::
        T_s = \frac{n_{10} + D - \left[(n_{10}+D)^2 - 4(n_9 + n_{10}D) \right]^{0.5}}{2}

    .. math:
        D = \frac{2G}{-F - (F^2 - 4EG)^{0.5}}

    .. math::
        E = \beta^2 + n_3 \beta + n_6

    .. math::
        F = n_1 \beta^2 + n_4\beta + n_7

    .. math::
        G = n_2\beta^2 + n_5\beta + n_8

    .. math::
        \beta = \left(P_{sat} \right)^{0.25}


    Parameters
    ----------
    Psat : float
        Vapor pressure at T [Pa]

    Returns
    -------
    T : float
        Temperature of water along the saturation curve at `Psat`, [K]

    Notes
    -----
    The range of validity of this equation is 273.15 K < T < 647.096 K, the
    IAPWS critical point.

    The coefficients `n1` to `n10` are (0.11670521452767E4, -0.72421316703206E6,
    -0.17073846940092E2, 0.12020824702470E5, -0.32325550322333E7, 0.14915108613530E2,
    -0.48232657361591E4, 0.40511340542057E6, -0.23855557567849, 0.65017534844798E3)

    Examples
    --------
    >>> Tsat_IAPWS(1E5)
    372.75591861133773

    References
    ----------
    .. [1] Kretzschmar, Hans-Joachim, and Wolfgang Wagner. International Steam
       Tables: Properties of Water and Steam Based on the Industrial
       Formulation IAPWS-IF97. Springer, 2019.
    '''
    B = sqrt(sqrt(P*1E-6))
    E = B*(B + -0.17073846940092E2) + 0.14915108613530E2
    F = B*(0.11670521452767E4*B + 0.12020824702470E5) + -0.48232657361591E4
    G = B*(-0.72421316703206E6*B + -0.32325550322333E7) + 0.40511340542057E6
    D = 2.0*G/(-F - sqrt(F*F - 4.0*E*G))
    n10 = 0.65017534844798E3
    x0 = (n10 + D)
    T = (n10 + D - sqrt(x0*x0 - 4.0*(-0.23855557567849+n10*D)))*0.5
    return T

### CSP Methods

def boiling_critical_relation(T, Tb, Tc, Pc):
    r'''Calculates vapor pressure of a fluid at arbitrary temperatures using a
    CSP relationship as in [1]_; requires a chemical's critical temperature
    and pressure as well as boiling point.

    The vapor pressure is given by:

    .. math::
        \ln P^{sat}_r = h\left( 1 - \frac{1}{T_r}\right)

    .. math::
        h = T_{br} \frac{\ln(P_c/101325)}{1-T_{br}}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Units are Pa. Formulation makes intuitive sense; a logarithmic form of
    interpolation.

    Examples
    --------
    Example as in [1]_ for ethylbenzene

    >>> boiling_critical_relation(347.2, 409.3, 617.1, 36E5)
    15209.467273093938

    References
    ----------
    .. [1] Reid, Robert C..; Prausnitz, John M.;; Poling, Bruce E.
       The Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Tbr = Tb/Tc
    Tr = T/Tc
    h = Tbr*log(Pc/101325.)/(1 - Tbr)
    return exp(h*(1-1/Tr))*Pc


def Lee_Kesler(T, Tc, Pc, omega):
    r'''Calculates vapor pressure of a fluid at arbitrary temperatures using a
    CSP relationship by [1]_; requires a chemical's critical temperature and
    acentric factor.

    The vapor pressure is given by:

    .. math::
        \ln P^{sat}_r = f^{(0)} + \omega f^{(1)}

    .. math::
        f^{(0)} = 5.92714-\frac{6.09648}{T_r}-1.28862\ln T_r + 0.169347T_r^6

    .. math::
        f^{(1)} = 15.2518-\frac{15.6875}{T_r} - 13.4721 \ln T_r + 0.43577T_r^6

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    This equation appears in [1]_ in expanded form.
    The reduced pressure form of the equation ensures predicted vapor pressure
    cannot surpass the critical pressure.

    Examples
    --------
    Example from [2]_; ethylbenzene at 347.2 K.

    >>> Lee_Kesler(347.2, 617.1, 36E5, 0.299)
    13078.694162949312

    References
    ----------
    .. [1] Lee, Byung Ik, and Michael G. Kesler. "A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States." AIChE Journal
       21, no. 3 (1975): 510-527. doi:10.1002/aic.690210313.
    .. [2] Reid, Robert C..; Prausnitz, John M.;; Poling, Bruce E.
       The Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Tr = T/Tc
    logTr = log(Tr)
    Tr6 = Tr*Tr
    Tr6 *= Tr6*Tr6
    f0 = 5.92714 - 6.09648/Tr - 1.28862*logTr+ 0.169347*Tr6
    f1 = 15.2518 - 15.6875/Tr - 13.4721*logTr + 0.43577*Tr6
    return exp(f0 + omega*f1)*Pc


def Ambrose_Walton(T, Tc, Pc, omega):
    r'''Calculates vapor pressure of a fluid at arbitrary temperatures using a
    CSP relationship by [1]_; requires a chemical's critical temperature and
    acentric factor.

    The vapor pressure is given by:

    .. math::
        \ln P_r=f^{(0)}+\omega f^{(1)}+\omega^2f^{(2)}

    .. math::
        f^{(0)}=\frac{-5.97616\tau + 1.29874\tau^{1.5}- 0.60394\tau^{2.5}
        -1.06841\tau^5}{T_r}

    .. math::
        f^{(1)}=\frac{-5.03365\tau + 1.11505\tau^{1.5}- 5.41217\tau^{2.5}
        -7.46628\tau^5}{T_r}

    .. math::
        f^{(2)}=\frac{-0.64771\tau + 2.41539\tau^{1.5}- 4.26979\tau^{2.5}
        +3.25259\tau^5}{T_r}

    .. math::
        \tau = 1-T_{r}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Somewhat more accurate than the :obj:`Lee_Kesler` formulation.

    Examples
    --------
    Example from [2]_; ethylbenzene at 347.25 K.

    >>> Ambrose_Walton(347.25, 617.15, 36.09E5, 0.304)
    13278.878504306222

    References
    ----------
    .. [1] Ambrose, D., and J. Walton. "Vapour Pressures up to Their Critical
       Temperatures of Normal Alkanes and 1-Alkanols." Pure and Applied
       Chemistry 61, no. 8 (1989): 1395-1403. doi:10.1351/pac198961081395.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    tau = 1.0 - Tr
    tau15 = tau**1.5
    tau25 = tau*tau15
    tau5 = tau25*tau25
    f0 = (-5.97616*tau + 1.29874*tau15 - 0.60394*tau25 - 1.06841*tau5)
    f1 = (-5.03365*tau + 1.11505*tau15 - 5.41217*tau25 - 7.46628*tau5)
    f2 = (-0.64771*tau + 2.41539*tau15 - 4.26979*tau25 + 3.25259*tau5)
    return Pc*exp((f0 + omega*(f1 + f2*omega))/Tr)


def Sanjari(T, Tc, Pc, omega):
    r'''Calculates vapor pressure of a fluid at arbitrary temperatures using a
    CSP relationship by [1]_. Requires a chemical's critical temperature,
    pressure, and acentric factor. Although developed for refrigerants,
    this model should have some general predictive ability.

    The vapor pressure of a chemical at `T` is given by:

    .. math::
        P^{sat} = P_c\exp(f^{(0)} + \omega f^{(1)} + \omega^2 f^{(2)})

    .. math::
        f^{(0)} = a_1 + \frac{a_2}{T_r} + a_3\ln T_r + a_4 T_r^{1.9}

    .. math::
        f^{(1)} = a_5 + \frac{a_6}{T_r} + a_7\ln T_r + a_8 T_r^{1.9}

    .. math::
        f^{(2)} = a_9 + \frac{a_{10}}{T_r} + a_{11}\ln T_r + a_{12} T_r^{1.9}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Psat : float
        Vapor pressure, [Pa]

    Notes
    -----
    a[1-12] are as follows:
    6.83377, -5.76051, 0.90654, -1.16906,
    5.32034, -28.1460, -58.0352, 23.57466,
    18.19967, 16.33839, 65.6995, -35.9739.

    For a claimed fluid not included in the regression, R128, the claimed AARD
    was 0.428%. A re-calculation using 200 data points from 125.45 K to
    343.90225 K evenly spaced by 1.09775 K as generated by NIST Webbook April
    2016 produced an AARD of 0.644%. It is likely that the author's regression
    used more precision in its coefficients than was shown here. Nevertheless,
    the function is reproduced as shown in [1]_.

    For Tc=808 K, Pc=1100000 Pa, omega=1.1571, this function actually declines
    after 770 K.

    Examples
    --------
    >>> Sanjari(347.2, 617.1, 36E5, 0.299)
    13651.916109552523

    References
    ----------
    .. [1] Sanjari, Ehsan, Mehrdad Honarmand, Hamidreza Badihi, and Ali
       Ghaheri. "An Accurate Generalized Model for Predict Vapor Pressure of
       Refrigerants." International Journal of Refrigeration 36, no. 4
       (June 2013): 1327-32. doi:10.1016/j.ijrefrig.2013.01.007.
    '''
    Tr = T/Tc
    Tr_inv = 1.0/Tr
    log_Tr = log(Tr)
    Tr_19 = Tr**1.9
    f0 = 6.83377 + -5.76051*Tr_inv + 0.90654*log_Tr + -1.16906*Tr_19
    f1 = 5.32034 + -28.1460*Tr_inv + -58.0352*log_Tr + 23.57466*Tr_19
    f2 = 18.19967 + 16.33839*Tr_inv + 65.6995*log_Tr + -35.9739*Tr_19
    return Pc*exp(f0 + omega*f1 + omega*omega*f2)


def Edalat(T, Tc, Pc, omega):
    r'''Calculates vapor pressure of a fluid at arbitrary temperatures using a
    CSP relationship by [1]_. Requires a chemical's critical temperature,
    pressure, and acentric factor. Claimed to have a higher accuracy than the
    Lee-Kesler CSP relationship.

    The vapor pressure of a chemical at `T` is given by:

    .. math::
        \ln(P^{sat}/P_c) = \frac{a\tau + b\tau^{1.5} + c\tau^3 + d\tau^6}
        {1-\tau}

    .. math::
        a = -6.1559 - 4.0855\omega

    .. math::
        b = 1.5737 - 1.0540\omega - 4.4365\times 10^{-3} d

    .. math::
        c = -0.8747 - 7.8874\omega

    .. math::
        d = \frac{1}{-0.4893 - 0.9912\omega + 3.1551\omega^2}

    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Psat : float
        Vapor pressure, [Pa]

    Notes
    -----
    [1]_ found an average error of 6.06% on 94 compounds and 1106 data points.

    Examples
    --------
    >>> Edalat(347.2, 617.1, 36E5, 0.299)
    13461.273080743307

    References
    ----------
    .. [1] Edalat, M., R. B. Bozar-Jomehri, and G. A. Mansoori. "Generalized
       Equation Predicts Vapor Pressure of Hydrocarbons." Oil and Gas Journal;
       91:5 (February 1, 1993).
    '''
    tau = 1. - T/Tc
    a = -6.1559 - 4.0855*omega
    c = -0.8747 - 7.8874*omega
    d = 1./(-0.4893 - 0.9912*omega + 3.1551*omega*omega)
    b = 1.5737 - 1.0540*omega - 4.4365E-3*d
    tau_15 = tau**1.5
    tau3 = tau_15*tau_15
    lnPr = (a*tau + b*tau_15 + c*tau3 + d*tau3*tau3)/(1.-tau)
    return exp(lnPr)*Pc


### Sublimation Pressure

def Psub_Clapeyron(T, Tt, Pt, Hsub_t):
    r'''Calculates sublimation pressure of a solid at arbitrary temperatures
    using an approximate themodynamic identity - the Clapeyron equation as
    described in [1]_ and [2]_.
    Requires a chemical's triple temperature, triple pressure, and triple
    enthalpy of sublimation.

    The sublimation pressure of a chemical at `T` is given by:

    .. math::
        \ln \frac{P}{P_{tp}} = -\frac{\Delta H_{sub}}{R}
        \left(\frac{1}{T}-\frac{1}{T_{tp}} \right)

    Parameters
    ----------
    T : float
        Temperature of solid [K]
    Tt : float
        Triple temperature of solid [K]
    Pt : float
        Truple pressure of solid [Pa]
    Hsub_t : float
        Enthalpy of fusion at the triple point of the chemical, [J/mol]

    Returns
    -------
    Psub : float
        Sublimation pressure, [Pa]

    Notes
    -----
    Does not seem to capture the decrease in sublimation pressure quickly
    enough.

    Examples
    --------
    >>> Psub_Clapeyron(250, Tt=273.15, Pt=611.0, Hsub_t=51100.0)
    76.06457150831804
    >>> Psub_Clapeyron(300, Tt=273.15, Pt=611.0, Hsub_t=51100.0)
    4577.282832876156

    References
    ----------
    .. [1] Goodman, B. T., W. V. Wilding, J. L. Oscarson, and R. L. Rowley.
       "Use of the DIPPR Database for the Development of QSPR Correlations:
       Solid Vapor Pressure and Heat of Sublimation of Organic Compounds."
       International Journal of Thermophysics 25, no. 2 (March 1, 2004):
       337-50. https://doi.org/10.1023/B:IJOT.0000028471.77933.80.
    .. [2] Feistel, Rainer, and Wolfgang Wagner. "Sublimation Pressure and
       Sublimation Enthalpy of H2O Ice Ih between 0 and 273.16K." Geochimica et
       Cosmochimica Acta 71, no. 1 (January 1, 2007): 36-45.
       https://doi.org/10.1016/j.gca.2006.08.034.
    '''
    return Pt*exp(Hsub_t*(T - Tt)/(R*T*Tt))


