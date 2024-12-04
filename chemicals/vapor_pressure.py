"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021 Caleb Bell
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
.. autofunction:: chemicals.vapor_pressure.Yaws_Psat
.. autofunction:: chemicals.vapor_pressure.TDE_PVExpansion
.. autofunction:: chemicals.vapor_pressure.Arrhenius_extrapolation

Fit Correlation Derivatives
---------------------------
.. autofunction:: chemicals.vapor_pressure.dAntoine_dT
.. autofunction:: chemicals.vapor_pressure.d2Antoine_dT2
.. autofunction:: chemicals.vapor_pressure.dWagner_dT
.. autofunction:: chemicals.vapor_pressure.d2Wagner_dT2
.. autofunction:: chemicals.vapor_pressure.dWagner_original_dT
.. autofunction:: chemicals.vapor_pressure.d2Wagner_original_dT2
.. autofunction:: chemicals.vapor_pressure.dTRC_Antoine_extended_dT
.. autofunction:: chemicals.vapor_pressure.d2TRC_Antoine_extended_dT2
.. autofunction:: chemicals.vapor_pressure.dYaws_Psat_dT
.. autofunction:: chemicals.vapor_pressure.d2Yaws_Psat_dT2
.. autofunction:: chemicals.vapor_pressure.dArrhenius_extrapolation_dT
.. autofunction:: chemicals.vapor_pressure.d2Arrhenius_extrapolation_dT2
.. autofunction:: chemicals.vapor_pressure.d3Arrhenius_extrapolation_dT3


Jacobians (for fitting)
-----------------------
.. autofunction:: chemicals.vapor_pressure.Wagner_fitting_jacobian
.. autofunction:: chemicals.vapor_pressure.Wagner_original_fitting_jacobian
.. autofunction:: chemicals.vapor_pressure.Antoine_fitting_jacobian
.. autofunction:: chemicals.vapor_pressure.Yaws_Psat_fitting_jacobian
.. autofunction:: chemicals.vapor_pressure.TRC_Antoine_extended_fitting_jacobian

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
.. autofunction:: chemicals.vapor_pressure.Arrhenius_parameters


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

.. data:: Psat_data_Alcock_elements

    Coefficients for the DIPPR 101 equation :obj:`chemicals.dippr.EQ101`,
    published in [5]_ and converted to provide base SI units (and use the
    natural logarithm). The conversions are as follows (original paper -> chemicals's `EQ101`):

    - A -> A * ln(10)
    - B -> B * ln(10)
    - C -> C (logarithm base-10 conversion cancels out, so C remains unchanged)
    - D -> 1e-3 * D * ln(10)
    - E is set to 1

.. data:: Psub_data_Alcock_elements

    Coefficients for the DIPPR 101 equation :obj:`chemicals.dippr.EQ101`,
    published in [5]_ and converted to provide base SI units (and use the
    natural logarithm). Note this is a sublimation pressure data set.
    Note that the `E` parameter in the :obj:`chemicals.dippr.EQ101` is 1
    for all chemicals, not the default of that function which is 0.0 and means
    the `D` parameter is not used.

.. data:: Psub_data_Landolt_Antoine

    Standard Antoine equation coefficients for sublimation pressure,
    as documented in the function
    :obj:`Antoine` and with data for ~1000 solids from [6]_, [7]_, and [8]_.
    Coefficients were altered to be in units of Pa and Kelvin with the
    exponential instead of base-10 power.

.. data:: Psat_data_Landolt_Antoine

    Standard Antoine equation coefficients for vapor pressure,
    as documented in the function
    :obj:`Antoine` and with data for ~6000 liquids from [6]_, [7]_, and [8]_.
    Coefficients were altered to be in units of Pa and Kelvin with the
    exponential instead of base-10 power.

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
.. [5] Alcock, C. B., V. P. Itkin, and M. K. Horrigan. "Vapour Pressure
    Equations for the Metallic Elements: 298-2500K." Canadian Metallurgical
    Quarterly 23, no. 3 (July 1, 1984): 309-13.
    https://doi.org/10.1179/cmq.1984.23.3.309.
.. [6] Hall, K. R. Vapor Pressure and Antoine Constants for Hydrocarbons,
   and S, Se, Te, and Halogen Containing Organic Compounds. Springer, 1999.
.. [7] Dykyj, J., and K. R. Hall. "Vapor Pressure and Antoine Constants for
    Oxygen Containing Organic Compounds". 2000.
.. [8] Hall, K. R. Vapor Pressure and Antoine Constants for Nitrogen
   Containing Organic Compounds. Springer, 2001.



The structure of each dataframe is shown below:

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.vapor_pressure.Psat_data_WagnerMcGarry

    In [3]: chemicals.vapor_pressure.Psat_data_WagnerPoling

    In [4]: chemicals.vapor_pressure.Psat_data_AntoinePoling

    In [5]: chemicals.vapor_pressure.Psat_data_AntoineExtended

    In [6]: chemicals.vapor_pressure.Psat_data_Perrys2_8

    In [7]: chemicals.vapor_pressure.Psat_data_VDI_PPDS_3

    In [8]: chemicals.vapor_pressure.Psat_data_Alcock_elements

    In [9]: chemicals.vapor_pressure.Psub_data_Alcock_elements

    In [10]: chemicals.vapor_pressure.Psub_data_Landolt_Antoine

    In [11]: chemicals.vapor_pressure.Psat_data_Landolt_Antoine
"""


__all__ = ['Antoine','dAntoine_dT', 'd2Antoine_dT2',
           'Wagner_original',  'dWagner_original_dT', 'd2Wagner_original_dT2',
           'Wagner', 'dWagner_dT', 'd2Wagner_dT2',
           'TRC_Antoine_extended', 'dTRC_Antoine_extended_dT',
           'd2TRC_Antoine_extended_dT2', 'dYaws_Psat_dT',
           'boiling_critical_relation', 'Lee_Kesler', 'Ambrose_Walton',
           'Edalat', 'Sanjari', 'Psat_IAPWS', 'dPsat_IAPWS_dT', 'Tsat_IAPWS',
           'Psub_Clapeyron', 'Yaws_Psat', 'd2Yaws_Psat_dT2',
           'Antoine_coeffs_from_point', 'Antoine_AB_coeffs_from_point',
           'DIPPR101_ABC_coeffs_from_point', 'Wagner_original_fitting_jacobian',
           'Wagner_fitting_jacobian', 'Yaws_Psat_fitting_jacobian',
           'Antoine_fitting_jacobian', 'TRC_Antoine_extended_fitting_jacobian',
           'TDE_PVExpansion', 'Arrhenius_extrapolation', 'Arrhenius_parameters']

from math import isinf, log10

from fluids.constants import R
from fluids.numerics import exp, log, sqrt, trunc_exp
from fluids.numerics import numpy as np

from chemicals.data_reader import data_source, register_df_source
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

folder = os_path_join(source_path, 'Vapor Pressure')

register_df_source(folder, 'Antoine Collection Poling.tsv')
register_df_source(folder, 'Table 2-8 Vapor Pressure of Inorganic and Organic Liquids.tsv')

register_df_source(folder, 'Alcock_Itkin_Horrigan_metalic_elements.tsv',
                   csv_kwargs={'dtype':{'Tmin': float, 'Tmax': float, 'A': float,
                            'B': float, 'C': float, 'D': float, 'E': float}})

register_df_source(folder, 'Landolt_antoine_sublimation_V20.tsv',
                   csv_kwargs={'dtype':{'Tmin': float, 'Tmax': float, 'A': float,
                            'B': float, 'C': float}})

register_df_source(folder, 'Landolt_antoine_V20.tsv',
                   csv_kwargs={'dtype':{'Tmin': float, 'Tmax': float, 'A': float,
                            'B': float, 'C': float}})

register_df_source(folder, 'Alcock_Itkin_Horrigan_metalic_elements_sublimation.tsv',
                   csv_kwargs={'dtype':{'Tmin': float, 'Tmax': float, 'A': float,
                            'B': float, 'C': float, 'D': float}}, index_col=1)

register_df_source(folder, 'Wagner Original McGarry.tsv', csv_kwargs={
        'dtype': {'A': float, 'B': float, 'C': float, 'D': float,
                 'Pc': float, 'Tc': float, 'Tmin': float}})

register_df_source(folder, 'Wagner Collection Poling.tsv', csv_kwargs={
        'dtype': {'A': float, 'B': float, 'C': float, 'D': float, 'Pc': float,
                  'Tc': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Antoine Extended Collection Poling.tsv', csv_kwargs={
    'dtype':{'A': float, 'B': float, 'C': float, 'Tc': float, 'to': float,
             'n': float, 'E': float, 'F': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'VDI PPDS Boiling temperatures at different pressures.tsv', csv_kwargs={
        'dtype':{'Tm': float, 'Tc': float, 'Pc': float, 'A': float,
                 'B': float, 'C': float, 'D': float}})

_vapor_pressure_dfs_loaded = False
@mark_numba_incompatible
def load_vapor_pressure_dfs():
    global Psat_data_WagnerMcGarry, Psat_values_WagnerMcGarry, Psat_data_AntoinePoling, Psat_values_AntoinePoling
    global Psat_data_WagnerPoling, Psat_values_WagnerPoling, Psat_data_AntoineExtended, Psat_values_AntoineExtended
    global Psat_data_Perrys2_8, Psat_values_Perrys2_8, Psat_data_VDI_PPDS_3, Psat_values_VDI_PPDS_3
    global Psat_data_Alcock_elements, Psat_values_Alcock_elements, Psub_data_Alcock_elements, Psub_values_Alcock_elements
    global Psub_values_Landolt_Antoine, Psub_data_Landolt_Antoine
    global Psat_values_Landolt_Antoine, Psat_data_Landolt_Antoine
    global _vapor_pressure_dfs_loaded
    if _vapor_pressure_dfs_loaded:
        return

    # 57463 bytes for df; 13720 bytes for numpy
    Psat_data_WagnerMcGarry = data_source('Wagner Original McGarry.tsv')
    Psat_values_WagnerMcGarry = np.array(Psat_data_WagnerMcGarry.values[:, 1:], dtype=float)

    # 58216 bytes for df; 13000 bytes for numpy
    Psat_data_AntoinePoling = data_source('Antoine Collection Poling.tsv')
    Psat_values_AntoinePoling = np.array(Psat_data_AntoinePoling.values[:, 1:], dtype=float)

    Psat_data_Alcock_elements = data_source('Alcock_Itkin_Horrigan_metalic_elements.tsv')
    Psat_values_Alcock_elements = np.array(Psat_data_Alcock_elements.values[:, 1:], dtype=float)

    Psub_data_Alcock_elements = data_source('Alcock_Itkin_Horrigan_metalic_elements_sublimation.tsv')
    Psub_values_Alcock_elements = np.array(Psub_data_Alcock_elements.values[:, 1:], dtype=float)

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

    Psub_data_Landolt_Antoine = data_source('Landolt_antoine_sublimation_V20.tsv')
    Psub_values_Landolt_Antoine = np.array(Psub_data_Landolt_Antoine.values[:, 1:], dtype=float)

    Psat_data_Landolt_Antoine = data_source('Landolt_antoine_V20.tsv')
    Psat_values_Landolt_Antoine = np.array(Psat_data_Landolt_Antoine.values[:, 1:], dtype=float)

    _vapor_pressure_dfs_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('Psat_data_WagnerMcGarry', 'Psat_values_WagnerMcGarry', 'Psat_data_AntoinePoling',
                    'Psat_values_AntoinePoling', 'Psat_data_WagnerPoling', 'Psat_values_WagnerPoling',
                    'Psat_data_AntoineExtended', 'Psat_values_AntoineExtended', 'Psat_data_Perrys2_8',
                    'Psat_values_Perrys2_8', 'Psat_data_VDI_PPDS_3', 'Psat_values_VDI_PPDS_3',
                    'Psat_data_Alcock_elements', 'Psat_values_Alcock_elements',
                    'Psub_data_Alcock_elements', 'Psub_values_Alcock_elements',
                    'Landolt_data_sublimation_Antoine', 'Psub_values_Landolt_Antoine',
                    'Psat_data_Landolt_Antoine', 'Psat_values_Landolt_Antoine'):
            load_vapor_pressure_dfs()
            return globals()[name]
        raise AttributeError(f"module {__name__} has no attribute {name}")
else:
    if can_load_data:
        load_vapor_pressure_dfs()


def TDE_PVExpansion(T, a1, a2, a3, a4=0.0, a5=0.0, a6=0.0, a7=0.0, a8=0.0):
    r'''Calculates vapor pressure or sublimation pressure of a chemical using
    the PVExpansion equation for vapor pressure or sublimation pressure.
    Parameters `a1`, `a2`, `a3`, `a4`, `a5`, `a6`, `a7`, and `a8`
    are chemical-dependent. Parameters
    can be found in various sources; however units of the coefficients used
    vary.

    .. math::
        \log P^{\text{sat}} = a_1 + \frac{a_2}{T} + a_3\ln(T) + a_4T + a_5T^2
        + \frac{a_6}{T^2} + a_7 T^6 + \frac{a_8}{T^4}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    a1 : float
        Regression parameter, [-]
    a2 : float
        Regression parameter, [-]
    a3 : float
        Regression parameter, [-]
    a4 : float
        Regression parameter, [-]
    a5 : float
        Regression parameter, [-]
    a6 : float
        Regression parameter, [-]
    a7 : float
        Regression parameter, [-]
    a8 : float
        Regression parameter, [-]

    Returns
    -------
    Psat : float
        Vapor pressure calculated with coefficients [Pa]

    Notes
    -----
    Coefficients in [1]_ produce a vapor pressure in kPa; add log(1000) to
    `a1` to make them produce vapor pressure in Pa.

    Examples
    --------
    Coefficients for sublimation pressure from [1]_:

    >>> TDE_PVExpansion(T=273.16, a1=23.7969+log(1000), a2=-11422, a3=0.177978)
    4.06220657398e-05

    References
    ----------
    .. [1] "ThermoData Engine (TDE103a V10.1) User`s Guide."
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-PhaseBoundaryLG/PVExpansion.htm
    '''
    T2 = T*T
    return trunc_exp(a1 + a2/T + a3*log(T) + a4*T + a5*T2 + a6/T2 + a7*T2*T2*T2 + a8/(T2*T2))

def Yaws_Psat(T, A, B, C, D, E):
    r'''Calculates vapor pressure of a chemical using the Yaws equation for
    vapor pressure.
    Parameters `A`, `B`, `C`, `D`, and `E` are chemical-dependent. Parameters
    can be found in numerous sources; however units of the coefficients used
    vary.

    .. math::
        \log_{10} P^{\text{sat}} = A + \frac{B}{T} + C\log_{10}(T) + DT + ET^2

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    A : float
        `A` parameter, [-]
    B : float
        `B` parameter, [K]
    C : float
        `C` parameter, [-]
    D : float
        `D` parameter, [1/K]
    E : float
        `E` parameter, [1/K^2]

    Returns
    -------
    Psat : float
        Vapor pressure calculated with coefficients [Pa]

    Notes
    -----
    Assumes coefficients are for calculating vapor pressure in Pascal.
    Coefficients should be consistent with input temperatures in Kelvin;

    **Converting units in input coefficients:**

        * **mmHg to Pa**: Add log10(101325/760)= 2.1249 to A.
        * **kPa to Pa**: Add log_{10}(1000)= 3 to A
        * **bar to Pa**: Add log_{10}(100000)= 5 to A

    Examples
    --------
    Acetone, coefficients from [1]_, at 400 K and with the conversion of `A`
    to obtain a result in Pa:

    >>> Yaws_Psat(T=400.0, A=28.588 + log10(101325/760), B=-2469, C=-7.351, D=2.8025E-10, E=2.7361E-6)
    708657.089106

    Coefficients for benzene from [2]_ at 400 K; that source outputs vapor
    pressure in kPa. That style of coefficients can be converted to `Pa`
    by adding `3` to `A`.

    >>> Yaws_Psat(T=400.0, A=39.7918+3, B=-2965.83, C=-12.073, D=0.0033269, E=1.58609e-6)
    352443.191026

    References
    ----------
    .. [1] Yaws, Carl L. Chemical Properties Handbook: Physical, Thermodynamic,
       Environmental, Transport, Safety, and Health Related Properties for
       Organic and Inorganic Chemicals. McGraw-Hill, 2001.
    .. [2] "ThermoData Engine (TDE103a V10.1) User`s Guide."
       https://trc.nist.gov/TDE/Help/TDE103a/Eqns-Pure-PhaseBoundaryLG/Yaws-VaporPressure.htm.
    '''
    exponent = (A + B/T + C*log10(T) + T*(D + E*T))
    if exponent > 308.0:
        return 1e308
    return 10.0**exponent

def dYaws_Psat_dT(T, A, B, C, D, E):
    r'''Calculates the first temperature derivative of vapor pressure of a
    chemical using the Yaws equation for vapor pressure.
    Parameters `A`, `B`, `C`, `D`, and `E` are chemical-dependent. Parameters
    can be found in numerous sources; however units of the coefficients used
    vary.

    .. math::
        \frac{\partial  P^{\text{sat}} }{\partial T} = 10^{A + \frac{B}{T}
        + \frac{C \log{\left(T \right)}}{\log{\left(10 \right)}} + D T
        + E T^{2}} \left(- \frac{B}{T^{2}} + \frac{C}{T \log{\left(10 \right)}}
        + D + 2 E T\right) \log{\left(10 \right)}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    A : float
        `A` parameter, [-]
    B : float
        `B` parameter, [K]
    C : float
        `C` parameter, [-]
    D : float
        `D` parameter, [1/K]
    E : float
        `E` parameter, [1/K^2]

    Returns
    -------
    dPsat_dT : float
        First temperature derivative of vapor pressure calculated with
        coefficients [Pa/K]

    Examples
    --------
    Benzene:

    >>> dYaws_Psat_dT(T=400.0, A=42.7918, B=-2965.83, C=-12.073, D=0.0033269, E=1.58609e-6)
    8134.87548930
    '''
    x0 = 2.302585092994046
    x1 = T*T
    x2 = 1.0/T
    x3 = C/x0
    return 10.0**(A + B*x2 + D*T + E*x1 + x3*log(T))*x0*(-B/x1 + D + 2.0*E*T + x2*x3)

def d2Yaws_Psat_dT2(T, A, B, C, D, E):
    r'''Calculates the second temperature derivative of vapor pressure of a
    chemical using the Yaws equation for vapor pressure.
    Parameters `A`, `B`, `C`, `D`, and `E` are chemical-dependent. Parameters
    can be found in numerous sources; however units of the coefficients used
    vary.

    .. math::
        \frac{\partial^2 P^{\text{sat}} }{\partial T^2} = 10^{A + \frac{B}{T}
        + \frac{C \log{\left(T \right)}}{\log{\left(10 \right)}} + D T
        + E T^{2}} \left(\frac{2 B}{T^{3}} - \frac{C}{T^{2} \log{\left(10
        \right)}} + 2 E + \left(- \frac{B}{T^{2}} + \frac{C}{T \log{\left(10
        \right)}} + D + 2 E T\right)^{2} \log{\left(10 \right)}\right)
        \log{\left(10 \right)}

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    A : float
        `A` parameter, [-]
    B : float
        `B` parameter, [K]
    C : float
        `C` parameter, [-]
    D : float
        `D` parameter, [1/K]
    E : float
        `E` parameter, [1/K^2]

    Returns
    -------
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure calculated with
        coefficients [Pa/K^2]

    Examples
    --------
    Benzene:

    >>> d2Yaws_Psat_dT2(T=400.0, A=42.7918, B=-2965.83, C=-12.073, D=0.0033269, E=1.58609e-6)
    141.7181045862
    '''
    x0 = 2.302585092994046
    x1 = 1.0/T
    x2 = T*T
    x3 = C/x0
    x4 = 2.0*E
    x5 = 1.0/x2
    x6 = (-B*x5 + D + T*x4 + x1*x3)
    return 10.0**(A + B*x1 + D*T + E*x2 + x3*log(T))*x0*(2.0*B*x1*x1*x1 + x0*x6*x6 - x3*x5 + x4)

def Yaws_Psat_fitting_jacobian(Ts, A, B, C, D, E):
    r'''Compute and return the Jacobian of the property predicted by
    the Yaws vapor pressure equation with respect to all the coefficients. This is
    used in fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A : float
        `A` parameter, [-]
    B : float
        `B` parameter, [K]
    C : float
        `C` parameter, [-]
    D : float
        `D` parameter, [1/K]
    E : float
        `E` parameter, [1/K^2]

    Returns
    -------
    jac : list[list[float, 5], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]

    '''
    x0 = 2.302585092994046
    N = len(Ts)
#    out = np.zeros((N, 5)) # numba: uncomment
    out = [[0.0]*5 for _ in range(N)] # numba: delete
    for i in range(N):
        T = Ts[i]
        r = out[i]
        x1 = 1.0/T
        x2 = T*T
        x3 = log(T)
        x4 = 10.0**(A + B*x1 + C*x3/x0 + D*T + E*x2)
        x5 = x0*x4
        r[0] = x5
        r[1] = x1*x5
        r[2] = x3*x4
        r[3] = T*x5
        r[4] = x2*x5
    return out

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


    Note that if `C` is negative and `T` is less than `C`, the predicted vapor
    pressure would be high and positive at those temperatures under `C`; and
    a singularity would occur at `T` == `C`. This implementation is corrected
    to return zero for the case of `T + C < 0.0`, which matches the intention
    of the Antoine equation.

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
    T_C = T + C
    if T_C <= 0.0:
        return 0.0
    return base**(A - B/(T_C))

def dAntoine_dT(T, A, B, C, base=10.0):
    r'''Calculates the first temperature derivative of vapor pressure of a
    chemical using the Antoine equation.
    Parameters `A`, `B`, and `C` are chemical-dependent.

    .. math::
        \frac{\partial  P^{\text{sat}} }{\partial T} =
        \frac{B \text{base}^{A - \frac{B}{C + T}} \log{\left(\text{base} \right)}}
        {\left(C + T\right)^{2}}

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
    dPsat_dT : float
        First temperature derivative of vapor pressure calculated with
        coefficients [Pa/K]

    Notes
    -----

    Examples
    --------
    Methane at 100 K:

    >>> dAntoine_dT(100.0, 8.7687, 395.744, -6.469)
    3591.4147747481
    '''
    T_C = T + C
    if T_C <= 0.0:
        return 0.0
    den = 1.0/(T_C)
    return B*base**(A - B*den)*log(base)*den*den

def d2Antoine_dT2(T, A, B, C, base=10.0):
    r'''Calculates the second temperature derivative of vapor pressure of a
    chemical using the Antoine equation.
    Parameters `A`, `B`, and `C` are chemical-dependent.

    .. math::
        \frac{\partial^2  P^{\text{sat}} }{\partial T^2} =
        \frac{B \text{base}^{A - \frac{B}{C + T}} \left(\frac{B \log{\left(
        \text{base} \right)}}{C + T} - 2\right) \log{\left(\text{base}
        \right)}}{\left(C + T\right)^{3}}

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
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure calculated with
        coefficients [Pa/K^2]

    Notes
    -----

    Examples
    --------
    Methane at 100 K:

    >>> d2Antoine_dT2(100.0, 8.7687, 395.744, -6.469)
    297.30093799054
    '''
    T_C = T + C
    if T_C <= 0.0:
        return 0.0
    den = 1.0/(T_C)
    log_base = log(base)
    x0 = B*den*log_base
    return x0*base**(A - B*den)*(x0 - 2.0)*den*den

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
    base : float, optional
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
    base : float, optional
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
    # The expression from SymPy is as follows
    #A = log(Psat*exp(T*dPsat_dT*Psat_inv))*log_base_inv
    # However, that expression has overflows a lot
    # Mathematical manipulation yields the following, which avoids overflows
    A = (T*dPsat_dT*Psat_inv + log(Psat))*log_base_inv

    B = T*T*dPsat_dT*log_base_inv*Psat_inv
    return (A, B)

def DIPPR101_ABC_coeffs_from_point(T, Psat, dPsat_dT, d2Psat_dT2):
    r'''Calculates the first three DIPPR101 coefficients `A`, `B`, and `C`
    from a known vapor pressure and its first and second temperature derivative.

    If the second derivative is infinity as is the case in some vapor pressure
    models at the critical point, only the `A` and `C` coefficients are fit,
    using the first derivative an the actual value of vapor pressure.

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
    >>> v = exp(A + B/T + C*log(T)) # doctest: +SKIP
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
    if isinf(d2Psat_dT2):
        # We cannot match the second derivative, so there will definitely be
        # a discontinuity
        # This means only two parameters can be obtained
        # If the A and B parameters are matched, the vapor pressure will go down
        # so the A and C paramters have to be matched.
        return (-T*dPsat_dT*log(T)/Psat + log(Psat), 0.0, T*dPsat_dT/Psat)
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
    Tc : float
        Critical temperature of fluid, [K]
    to : float
        Fit temperature-transition parameter, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    n : float
        Fit parameter, [-]
    E : float
        Fit parameter, [-]
    F : float
        Fit parameter, [-]

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

    >>> TRC_Antoine_extended(T=180.0, Tc=227.51, to=-120., A=8.95894,
    ... B=510.595, C=-15.95, n=2.41377, E=-93.74, F=7425.9)
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
    T_C = T + C
    if T_C <= 0.0:
        return 0.0
    return 10.**(A - B/T_C + 0.43429*x**n + x4*x4*(E + F*x4))

def dTRC_Antoine_extended_dT(T, Tc, to, A, B, C, n, E, F):
    r'''Calculates the first temperature derivative of vapor pressure of a
    chemical using the TRC Extended Antoine equation.

    .. math::
        \frac{\partial  P^{\text{sat}} }{\partial T} =
        10^{A - \frac{B}{C + T} + \frac{E \left(T - T_{ref} - to\right)^{8}}
        {T_{c}^{8}} + \frac{F \left(T - T_{ref} - to\right)^{12}}{T_{c}^{12}}
        + f \left(\frac{T - T_{ref} - to}{T_{c}}\right)^{n}}
        \left(\frac{B}{\left(C + T\right)^{2}} + \frac{8 E \left(T - T_{ref}
        - to\right)^{7}}{T_{c}^{8}} + \frac{12 F \left(T - T_{ref} - to
        \right)^{11}}{T_{c}^{12}} + \frac{f n \left(\frac{T - T_{ref} - to}
        {T_{c}}\right)^{n}}{T - T_{ref} - to}\right) \log{\left(10 \right)}

    .. math::
        x = \max \left(\frac{T-t_o-273.15}{T_c}, 0 \right)

    .. math::
        T_{ref} = 273.15 \text{ K}

    .. math::
        f = 0.43429

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Tc : float
        Critical temperature of fluid, [K]
    to : float
        Fit temperature-transition parameter, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    n : float
        Fit parameter, [-]
    E : float
        Fit parameter, [-]
    F : float
        Fit parameter, [-]

    Returns
    -------
    dPsat_dT : float
        First temperature derivative of vapor pressure calculated with
        coefficients [Pa/K]

    Notes
    -----

    Examples
    --------
    Tetrafluoromethane at 180 K:

    >>> dTRC_Antoine_extended_dT(T=180.0, Tc=227.51, to=-120., A=8.95894,
    ... B=510.595, C=-15.95, n=2.41377, E=-93.74, F=7425.9)
    31219.6061263
    '''
    T_ref = 273.15
    f = 0.43429
    ln10 = 2.302585092994046
    x = (T - to - 273.15)/Tc
    if x < 0.0:
        return dAntoine_dT(T, A, B, C, base=10.0)
    x0 = C + T
    x1 = -T + T_ref + to
    x2 = E/Tc**8
    x3 = F/Tc**12
    x4 = f*(-x1/Tc)**n
    return (-10**(A - B/x0 + x1**12*x3 + x1**8*x2 + x4)*(-B/x0**2 + n*x4/x1
                  + 12.0*x1**11*x3 + 8.0*x1**7*x2)*ln10)

def d2TRC_Antoine_extended_dT2(T, Tc, to, A, B, C, n, E, F):
    r'''Calculates the second temperature derivative of vapor pressure of a
    chemical using the TRC Extended Antoine equation.

    .. math::
        \frac{\partial^2  P^{\text{sat}} }{\partial T^2} =
        10^{A - \frac{B}{C + T} + \frac{E \left(- T + T_{ref} + to\right)^{8}}
        {T_{c}^{8}} + \frac{F \left(- T + T_{ref} + to\right)^{12}}{T_{c}^{12}}
        + f \left(- \frac{- T + T_{ref} + to}{T_{c}}\right)^{n}} \left(
        - \frac{2 B}{\left(C + T\right)^{3}} + \frac{56 E \left(- T + T_{ref}
        + to\right)^{6}}{T_{c}^{8}} + \frac{132 F \left(- T + T_{ref}
        + to\right)^{10}}{T_{c}^{12}} + \frac{f n^{2} \left(- \frac{- T
        + T_{ref} + to}{T_{c}}\right)^{n}}{\left(- T + T_{ref} + to\right)^{2}}
        - \frac{f n \left(- \frac{- T + T_{ref} + to}{T_{c}}\right)^{n}}
        {\left(- T + T_{ref} + to\right)^{2}} + \left(- \frac{B}{\left(C
        + T\right)^{2}} + \frac{8 E \left(- T + T_{ref} + to\right)^{7}}
        {T_{c}^{8}} + \frac{12 F \left(- T + T_{ref} + to\right)^{11}}
        {T_{c}^{12}} + \frac{f n \left(- \frac{- T + T_{ref} + to}
        {T_{c}}\right)^{n}}{- T + T_{ref} + to}\right)^{2}
        \log{\left(10 \right)}\right) \log{\left(10 \right)}

    .. math::
        x = \max \left(\frac{T-t_o-273.15}{T_c}, 0 \right)

    .. math::
        T_{ref} = 273.15 \text{ K}

    .. math::
        f = 0.43429

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    Tc : float
        Critical temperature of fluid, [K]
    to : float
        Fit temperature-transition parameter, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    n : float
        Fit parameter, [-]
    E : float
        Fit parameter, [-]
    F : float
        Fit parameter, [-]

    Returns
    -------
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure calculated with
        coefficients [Pa/K]

    Notes
    -----

    Examples
    --------
    Tetrafluoromethane at 180 K:

    >>> d2TRC_Antoine_extended_dT2(T=180.0, Tc=227.51, to=-120., A=8.95894,
    ... B=510.595, C=-15.95, n=2.41377, E=-93.74, F=7425.9)
    1022.550368944
    '''
    T_ref = 273.15
    f = 0.43429
    x = (T - to - 273.15)/Tc
    if x < 0.0:
        return d2Antoine_dT2(T, A, B, C, base=10.0)
    x0 = 2.302585092994046
    x1 = C + T
    x2 = -T + T_ref + to
    x3 = E/Tc**8
    x4 = F/Tc**12
    x5 = f*(-x2/Tc)**n
    x6 = x2**(-2)
    x7 = n*x5
    return (10**(A - B/x1 + x2**12*x4 + x2**8*x3 + x5)*x0*(-2.0*B/x1**3 + n**2*x5*x6
            + x0*(-B/x1**2 + 12.0*x2**11*x4 + 8.0*x2**7*x3 + x7/x2)**2
            + 132.0*x2**10*x4 + 56.0*x2**6*x3 - x6*x7))


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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        Cubic power tau coefficient, [-]
    d : float
        6 power tau coefficient, [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Warning: Pc is often treated as adjustable constant.
    This is also called the PPDS1 equation [3]_.

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
    .. [3] PPDS2 Temperature-Dependent Equation Forms. National Engineering
       Laboratory, 2004
       https://web.archive.org/web/20050510061545/http://www.ppds.co.uk/library/pdf/PPDS_EquationForms.pdf
    '''
    Tr = T/Tc
    if Tr > 1.0:
        Tr = 1.0
    tau = 1.0 - Tr
    tau2 = tau*tau
    try:
        tau_Tr = tau/Tr
    except:
        # T = 0; Tr = 0
        return 0.0
    return Pc*exp(((d*tau2*tau + c)*tau2 + a + b*sqrt(tau))*tau_Tr)

def dWagner_original_dT(T, Tc, Pc, a, b, c, d):
    r'''Calculates first temperature derivative of vapor pressure using the
    Wagner equation (3, 6 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \frac{\partial  P^{\text{sat}} }{\partial T} =
        P_{c} \left(\frac{T_{c} \left(- \frac{a}{T_{c}} - \frac{1.5 b
        \tau^{0.5}}{T_{c}} - \frac{3 c \tau^{2}}{T_{c}} - \frac{6 d \tau^{5}}
        {T_{c}}\right)}{T} - \frac{T_{c} \left(a \tau + b \tau^{1.5}
        + c \tau^{3} + d \tau^{6}\right)}{T^{2}}\right) e^{\frac{T_{c} \left(a
        \tau + b \tau^{1.5} + c \tau^{3} + d \tau^{6}\right)}{T}}

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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        Cubic power tau coefficient, [-]
    d : float
        6 power tau coefficient, [-]

    Returns
    -------
    dPsat_dT : float
        First temperature derivative of vapor pressure at T [Pa/K]

    Notes
    -----

    Examples
    --------
    Methane at 100 K.

    >>> dWagner_original_dT(100.0, 190.53, 4596420., a=-6.00435, b=1.1885,
    ... c=-0.834082, d=-1.22833)
    3593.70783283
    '''
    Tr = T/Tc
    tau = 1.0 - Tr
    tau2 = tau*tau
    tau3 = tau2*tau
    try:
        T_inv = 1.0/T
        Tr_inv = Tc*T_inv
    except:
        # T = 0; Tr = 0
        return 0.0
    tau_rt = sqrt(tau)
    y0 = a*tau + b*tau*tau_rt + tau3*(c + d*tau3)
    exp_term = exp(Tr_inv*y0)
    if exp_term == 0.0:
        # Avoid underflowing to nan
        return 0.0
    dPsat_dT = Pc*(T_inv*((-a - 1.5*b*tau_rt - tau2*(3.0*c + 6.0*d*tau3))
                   - Tc*y0*T_inv)
                    *exp_term)
    return dPsat_dT

def d2Wagner_original_dT2(T, Tc, Pc, a, b, c, d):
    r'''Calculates second temperature derivative of vapor pressure using the
    Wagner equation (3, 6 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \frac{\partial^2  P^{\text{sat}} }{\partial T^2} =
        \frac{P_{c} \left(\frac{\frac{0.75 b}{\tau^{0.5}} + 6 c \tau + 30 d
        \tau^{4}}{T_{c}} + \frac{2 \left(a + 1.5 b \tau^{0.5} + 3 c \tau^{2}
        + 6 d \tau^{5}\right)}{T} + \frac{36 \left(\frac{a}{6} + 0.25 b
        \tau^{0.5} + \frac{c \tau^{2}}{2} + d \tau^{5} - \frac{T_{c} \left(
        - a \tau - b \tau^{1.5} - c \tau^{3} - d \tau^{6}\right)}{6 T}
        \right)^{2}}{T} - \frac{2 T_{c} \left(- a \tau - b \tau^{1.5}
        - c \tau^{3} - d \tau^{6}\right)}{T^{2}}\right) e^{- \frac{T_{c}
        \left(- a \tau - b \tau^{1.5} - c \tau^{3} - d \tau^{6}\right)}{T}}}{T}

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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        Cubic power tau coefficient, [-]
    d : float
        6 power tau coefficient, [-]

    Returns
    -------
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure at T [Pa/K^2]

    Notes
    -----
    This second derivative is infinity at T == Tc.

    Examples
    --------
    Methane at 100 K.

    >>> d2Wagner_original_dT2(100.0, 190.53, 4596420., a=-6.00435, b=1.1885,
    ... c=-0.834082, d=-1.22833)
    296.87593368224
    '''
    Tr = T/Tc
    tau = 1.0 - Tr
    tau2 = tau*tau
    tau3 = tau2*tau
    try:
        T_inv = 1.0/T
#        Tr_inv = Tc*T_inv
    except:
        # T = 0; Tr = 0
        return 0.0
#    y0 = a*tau + b*tau*tau_rt + tau3*(c + d*tau3)
#
#    tau_rt = sqrt(tau)
    tau_rt = sqrt(tau)
    x1 = Tc*tau*(a + b*tau_rt + tau2*(c + d*tau3))
    x2 = T_inv*x1
    exp_term = exp(x2)
    if exp_term == 0.0:
        # Avoid underflowing to nan
        return 0.0
    x4 = b*tau_rt
    x5 = c*tau2
    x6 = d*tau2*tau3
    x7 = (a*(1.0/6.0) + x2*(1.0/6.0) + 0.25*x4 + 0.5*x5 + x6)
    return (Pc*T_inv*T_inv*(2.0*(a + 1.5*x4 + 3.0*x5 + 6.0*x6) + 36.0*x7*x7
                   + (0.75*b/tau_rt + 6.0*c*tau + 30.0*d*tau2*tau2)*Tr
                   + 2.0*x1*T_inv)*exp_term)

def Wagner_original_fitting_jacobian(Ts, Tc, Pc, a, b, c, d):
    r'''Calculates the jacobian of the Wagner (3, 6) vapor pressure equation
    for use in fitting these parameters when experimental values are known.

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of fluid data points, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        Cubic power tau coefficient, [-]
    d : float
        6 power tau coefficient, [-]

    Returns
    -------
    jac : list[list[float, 4], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]
    '''
    N = len(Ts)
#    out = np.zeros((N, 4)) # numba: uncomment
    out = [[0.0]*4 for _ in range(N)] # numba: delete
    for i in range(N):
        Tr = Ts[i]/Tc
        tau = 1.0 - Tr
        x2 = tau*sqrt(tau)
        x3 = x2*x2
        x4 = x3*x3
        x1 = 1.0/Tr
        x5 = Pc*x1*exp(x1*(a*tau + b*x2 + c*x3 + d*x4))
        row = out[i]
        row[0] = tau*x5
        row[1] = x2*x5
        row[2] = x3*x5
        row[3] = x4*x5
    return out

def Wagner_fitting_jacobian(Ts, Tc, Pc, a, b, c, d):
    r'''Calculates the jacobian of the Wagner (2.5, 5) vapor pressure equation
    for use in fitting these parameters when experimental values are known.

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of fluid data points, [K]
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        2.5 power tau coefficient, [-]
    d : float
        5 power tau coefficient, [-]

    Returns
    -------
    jac : list[list[float, 4], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]
    '''
    N = len(Ts)
#    out = np.zeros((N, 4)) # numba: uncomment
    out = [[0.0]*4 for _ in range(N)] # numba: delete
    for i in range(N):
        Tr = Ts[i]/Tc
        tau = 1.0 - Tr
        x2 = tau*sqrt(tau)
        x3 = tau*x2
        x4 = x3*x3
        x1 = 1.0/Tr
        x5 = Pc*x1*exp(x1*(a*tau + b*x2 + c*x3 + d*x4))
        row = out[i]
        row[0] = tau*x5
        row[1] = x2*x5
        row[2] = x3*x5
        row[3] = x4*x5
    return out

def Antoine_fitting_jacobian(Ts, A, B, C, base=10.0):
    r'''Calculates the jacobian of the Antoine vapor pressure equation
    for use in fitting these parameters when experimental values are known.

    Requires three coefficients specific to each chemical.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of fluid data points, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    base : float, optional
        Optional base of logarithm; 10 by default, [-]

    Returns
    -------
    jac : list[list[float, 3], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]
    '''
    N = len(Ts)
#    out = np.zeros((N, 3)) # numba: uncomment
    out = [[0.0]*3 for _ in range(N)] # numba: delete
    ln_base = log(base)
    for i in range(N):
        row = out[i]
        x0 = C + Ts[i]
        if x0 <= 0.0:
            row[0] = 0.0
            row[1] = 0.0
            row[2] = 0.0
        else:
            x1 = 1.0/x0
            x2 = base**(A - B*x1)*ln_base
            row[0] = x2
            row[1] = -x1*x2
            row[2] =  B*x2*x1*x1
    return out

def TRC_Antoine_extended_fitting_jacobian(Ts, Tc, to, A, B, C, n, E, F):
    r'''Calculates the jacobian of the TRC Antoine extended vapor pressure
    equation for use in fitting these parameters when experimental values are
    known.

    Requires 7 coefficients specific to each chemical.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of fluid data points, [K]
    Tc : float
        Critical temperature of fluid, [K]
    to : float
        Fit temperature-transition parameter, [K]
    A : float
        Antoine `A` parameter, [-]
    B : float
        Antoine `B` parameter, [K]
    C : float
        Antoine `C` parameter, [K]
    n : float
        Fit parameter, [-]
    E : float
        Fit parameter, [-]
    F : float
        Fit parameter, [-]


    Returns
    -------
    jac : list[list[float, 7], len(Ts)]
        Matrix of derivatives of the equation with respect to the fitting
        parameters, [various]
    '''
    N = len(Ts)
#    out = np.zeros((N, 7)) # numba: uncomment
    out = [[0.0]*7 for _ in range(N)] # numba: delete
    ln_base = 2.302585092994046 #log(10)
    c0 = 273.15
    c1 = 0.43429
    Tc_inv = 1.0/Tc
    Tc_inv2 = Tc_inv*Tc_inv
    Tc_inv4 = Tc_inv2*Tc_inv2
    Tc_n8 = Tc_inv4*Tc_inv4
    Tc_n12 = Tc_n8*Tc_inv4

    v0 = - to - c0

    for i in range(N):
        row = out[i]
        T = Ts[i]
        x = (T + v0)*Tc_inv
        if x < 0.0:
            x0 = C + T
            if x0 > 0.0:
                x1 = 1.0/x0
                x2 = 10.0**(A - B*x1)*ln_base
                x1x2 = x1*x2
                row[1] = x2
                row[2] = -x1x2
                row[3] = B*x1x2*x1
        else:
            x1 = -T + c0 + to
            x1_2 = x1*x1
            x1_4 = x1_2*x1_2
            x1_8 = x1_4*x1_4
            x2 = -x1*Tc_inv
            x3 = c1*x2**n
            x5 = Tc_n8*(E + F*x1_4*Tc_inv4)
            x6 = C + T
            x7 = 1.0/x6
            x9 = 10.0**(A - B*x7 + x3 + x5*x1_8)*ln_base
            row[0] = x9*(4.0*F*Tc_n12*x1_8*x1_2*x1 + n*x3/x1 + 8.0*x1_4*x1_2*x1*x5)
            row[1] = x9
            row[2] = -x7*x9
            row[3] = B*x9*x7*x7
            row[4] = x3*x9*log(x2)
            row[5] = Tc_n8*x1_8*x9
            row[6] = Tc_n12*x1_8*x1_4*x9
    return out

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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        2.5 power tau coefficient, [-]
    d : float
        5 power tau coefficient, [-]

    Returns
    -------
    Psat : float
        Vapor pressure at T [Pa]

    Notes
    -----
    Warning: Pc is often treated as adjustable constant.
    This is also called the PPDS16 equation [3]_.

    Examples
    --------
    Methane, coefficients from [2]_, at 100 K.

    >>> Wagner(100., 190.551, 4599200, -6.02242, 1.26652, -0.5707, -1.366)
    34415.004762637

    References
    ----------
    .. [1] Wagner, W. "New Vapour Pressure Measurements for Argon and Nitrogen and
       a New Method for Establishing Rational Vapour Pressure Equations."
       Cryogenics 13, no. 8 (August 1973): 470-82. doi:10.1016/0011-2275(73)90003-9
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [3] PPDS2 Temperature-Dependent Equation Forms. National Engineering
       Laboratory, 2004
       https://web.archive.org/web/20050510061545/http://www.ppds.co.uk/library/pdf/PPDS_EquationForms.pdf
    '''
    Tr = T/Tc
    if Tr > 1.0:
        Tr = 1.0
    tau = 1.0 - Tr
    tau_rt = sqrt(tau)
    tau15 = tau*tau_rt
    tau25 = tau*tau15
    return Pc*exp((a + b*tau_rt + tau15*(c + d*tau25))*tau/Tr)

def dWagner_dT(T, Tc, Pc, a, b, c, d):
    r'''Calculates the first temperature derivative of vapor pressure using the
    Wagner equation (2.5, 5 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \frac{\partial  P^{\text{sat}} }{\partial T} =
        P_{c} \left(\frac{T_{c} \left(- \frac{a}{T_{c}} - \frac{1.5 b
        \tau^{0.5}}{T_{c}} - \frac{2.5 c \tau^{1.5}}{T_{c}} - \frac{5 d
        \tau^{4}}{T_{c}}\right)}{T} - \frac{T_{c} \left(a \tau + b \tau^{1.5}
        + c \tau^{2.5} + d \tau^{5}\right)}{T^{2}}\right) e^{\frac{T_{c}
        \left(a \tau + b \tau^{1.5} + c \tau^{2.5} + d \tau^{5}\right)}{T}}

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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        2.5 power tau coefficient, [-]
    d : float
        5 power tau coefficient, [-]

    Returns
    -------
    dPsat_dT : float
        First temperature derivative of vapor pressure at T [Pa/K]

    Notes
    -----

    Examples
    --------
    Methane at 100 K.

    >>> dWagner_dT(100., 190.551, 4599200, -6.02242, 1.26652, -0.5707, -1.366)
    3587.2910498076

    '''
    tau = 1.0 - T/Tc
    tau2 = tau*tau
    try:
        x0 = 1.0/T
    except:
        # T = 0
        return 0.0
    tau_rt = sqrt(tau)
    x1 = tau*tau_rt
    x2 = Tc*x0*(a*tau + b*x1 + tau2*(c*tau_rt + d*tau2*tau))
    exp_term = exp(x2)
    if exp_term == 0.0:
        # Avoid nan issues
        return 0.0
    return -Pc*x0*(a + 1.5*b*tau_rt + 2.5*c*x1 + 5.0*d*tau2*tau2 + x2)*exp_term

def d2Wagner_dT2(T, Tc, Pc, a, b, c, d):
    r'''Calculates the second temperature derivative of vapor pressure using the
    Wagner equation (2.5, 5 form).

    Requires critical temperature and pressure as well as four coefficients
    specific to each chemical.

    .. math::
        \frac{\partial^2  P^{\text{sat}} }{\partial T^2} =
        \frac{P_{c} \left(\frac{\frac{0.75 b}{\tau^{0.5}} + 3.75 c \tau^{0.5}
        + 20 d \tau^{3}}{T_{c}} + \frac{2 \left(a + 1.5 b \tau^{0.5}
        + 2.5 c \tau^{1.5} + 5 d \tau^{4}\right)}{T} + \frac{25 \left(
        \frac{a}{5} + 0.3 b \tau^{0.5} + 0.5 c \tau^{1.5} + d \tau^{4}
        - \frac{T_{c} \left(- a \tau - b \tau^{1.5} - c \tau^{2.5}
        - d \tau^{5}\right)}{5 T}\right)^{2}}{T} - \frac{2 T_{c} \left(- a
        \tau - b \tau^{1.5} - c \tau^{2.5} - d \tau^{5}\right)}{T^{2}}\right)
        e^{- \frac{T_{c} \left(- a \tau - b \tau^{1.5} - c \tau^{2.5}
        - d \tau^{5}\right)}{T}}}{T}

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
    a : float
        Linear tau coefficient, [-]
    b : float
        1.5 power tau coefficient, [-]
    c : float
        2.5 power tau coefficient, [-]
    d : float
        5 power tau coefficient, [-]

    Returns
    -------
    d2Psat_dT2 : float
        Second temperature derivative of vapor pressure at T [Pa/K^2]

    Notes
    -----
    This second derivative is infinity at T == Tc.

    Examples
    --------
    Methane at 100 K.

    >>> d2Wagner_dT2(100., 190.551, 4599200, -6.02242, 1.26652, -0.5707, -1.366)
    296.7091513877
    '''
    tau = 1.0 - T/Tc
    tau_rt = sqrt(tau)
    tau2 = tau*tau
    try:
        T_inv = 1.0/T
    except:
        # T = 0
        return 0.0
    x1 = tau*tau_rt
    x2 = Tc*(a*tau + b*x1 + tau*(c*x1 + d*tau2*tau2))
    x3 = T_inv*x2
    x5 = b*tau_rt
    x6 = c*x1
    x7 = d*tau2*tau2
    exp_term = exp(x3)
    if exp_term == 0.0:
        # Avoid nan issues
        return 0.0
    return (Pc*T_inv*(2.0*T_inv*(a + 1.5*x5 + 2.5*x6 + 5.0*x7) + 25.0*T_inv*(a*(1.0/5.0) + x3*(1.0/5.0)
            + 0.3*x5 + 0.5*x6 + x7)**2 + (0.75*b/tau_rt + 3.75*c*tau_rt
            + 20.0*d*tau*tau2)/Tc + 2.0*x2*T_inv*T_inv)*exp_term)

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
    P : float
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
    return exp(h*(1.0 - 1.0/Tr))*Pc


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
    if T > Tc:
        T = Tc
    Tr = T/Tc
    tau = 1.0 - Tr
    tau15 = tau*sqrt(tau)
    tau25 = tau*tau15
    tau5 = tau25*tau25
    if omega < 0.0:
        # The omega term is based on positive omegas only, prevent it from
        # going negative
        omega = 0.0
    f0 = (-5.97616*tau + 1.29874*tau15 - 0.60394*tau25 - 1.06841*tau5)
    f1 = (-5.03365*tau + 1.11505*tau15 - 5.41217*tau25 - 7.46628*tau5)
    f2 = (-0.64771*tau + 2.41539*tau15 - 4.26979*tau25 + 3.25259*tau5)
    return Pc*trunc_exp((f0 + omega*(f1 + f2*omega))/Tr)


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
    tau_15 = sqrt(tau)*tau
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
    ans = Pt*exp(Hsub_t*(T - Tt)/(R*T*Tt))
    # Truncation to avoid issues in later calculations
    return max(ans, 1e-200)


def Arrhenius_parameters(T, P, dP_dT):
    r'''Calculates parameters for Arrhenius-style vapor pressure extrapolation. This
    converts a vapor pressure and its temperature derivative into a slope suitable
    for extrapolation in ln(P) vs 1/T coordinates.
    
    .. math::
        \text{slope} = -T^2\frac{d \ln P}{dT} = -\frac{T^2}{P}\frac{dP}{dT}
        
    Parameters
    ----------
    T : float
        Temperature point for extrapolation [K]
    P : float
        Pressure at the temperature point [Pa] 
    dP_dT : float
        Temperature derivative of pressure at the temperature point [Pa/K]
        
    Returns
    -------
    T : float
        Temperature point [K]
    P : float
        Pressure at temperature point [Pa] 
    slope : float
        Slope d(ln P)/d(1/T) at the temperature point
    
    Notes
    -----
    Useful for extrapolating vapor pressures via the Clausius-Clapeyron relation.
    The slope represents the negative of enthalpy of vaporization divided by the 
    gas constant.
    
    Examples
    --------
    >>> Arrhenius_parameters(400.0, 1E5, 1E3)
    (400.0, 100000.0, -1600.0)
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return (T, P, -T*T*dP_dT/P)

def Arrhenius_extrapolation(T, T_ref, P_ref, slope):
    r'''Calculates extrapolated vapor pressure using Arrhenius-style extrapolation
    in ln(P) vs 1/T coordinates. This form of extrapolation is appropriate for 
    vapor pressures following the Clausius-Clapeyron relation.
    
    .. math::
        \ln P = \ln P_{ref} + \text{slope}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)

    Parameters
    ----------
    T : float
        Temperature to extrapolate to [K]
    T_ref : float
        Reference temperature point [K]
    P_ref : float
        Pressure at reference point [Pa]
    slope : float
        Slope d(ln P)/d(1/T) at reference point
        
    Returns
    -------
    P : float
        Extrapolated vapor pressure [Pa]
    
    Notes
    -----
    
    Examples
    --------
    >>> Arrhenius_extrapolation(300.0, 400.0, 1E5, -1600)
    26359.713811
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return P_ref*exp(slope*(1.0/T - 1.0/T_ref))



def dArrhenius_extrapolation_dT(T, T_ref, P_ref, slope):
    r'''Calculates the first temperature derivative of vapor pressure using the 
    Arrhenius-style vapor pressure extrapolation in ln(P) vs 1/T coordinates.

    .. math::
        \frac{\partial P}{\partial T} = -\frac{P_{ref} \cdot \text{slope} \cdot 
        e^{\text{slope}(1/T - 1/T_{ref})}}{T^2}

    Parameters
    ----------
    T : float
        Temperature to evaluate derivative at [K]
    T_ref : float
        Reference temperature point [K]
    P_ref : float
        Pressure at reference point [Pa]
    slope : float
        Slope d(ln P)/d(1/T) at reference point

    Returns
    -------
    dP_dT : float
        First temperature derivative of vapor pressure [Pa/K]

    Notes
    -----
    This derivative is useful for numerical solutions and for evaluating the 
    rate of change of vapor pressure with temperature.

    Examples
    --------
    >>> dArrhenius_extrapolation_dT(300.0, 400.0, 1E5, -1600)
    468.617134427
    '''    
    return -P_ref*slope*exp(slope*(1.0/T - 1.0/T_ref))/(T*T)

def d2Arrhenius_extrapolation_dT2(T, T_ref, P_ref, slope):
    r'''Calculates the second temperature derivative of vapor pressure using the
    Arrhenius-style vapor pressure extrapolation in ln(P) vs 1/T coordinates.

    .. math::
        \frac{\partial^2 P}{\partial T^2} = \frac{P_{ref} \cdot \text{slope} \cdot 
        e^{\text{slope}(1/T - 1/T_{ref})}(\text{slope} + 2T)}{T^4}

    Parameters
    ----------
    T : float
        Temperature to evaluate derivative at [K]
    T_ref : float
        Reference temperature point [K]
    P_ref : float
        Pressure at reference point [Pa]
    slope : float
        Slope d(ln P)/d(1/T) at reference point

    Returns
    -------
    d2P_dT2 : float
        Second temperature derivative of vapor pressure [Pa/K^2]

    Notes
    -----
    This derivative is useful for numerical solutions requiring higher-order 
    derivatives and for analyzing the curvature of vapor pressure with temperature.

    Examples
    --------
    >>> d2Arrhenius_extrapolation_dT2(300.0, 400.0, 1E5, -1600)
    5.206857049199541
    '''
    e = exp(slope*(1.0/T - 1.0/T_ref))
    return P_ref*slope*e*(slope + 2.0*T)/(T*T*T*T)

def d3Arrhenius_extrapolation_dT3(T, T_ref, P_ref, slope):
    r'''Calculates the third temperature derivative of vapor pressure using the
    Arrhenius-style vapor pressure extrapolation in ln(P) vs 1/T coordinates.

    .. math::
        \frac{\partial^3 P}{\partial T^3} = -\frac{P_{ref} \cdot \text{slope} \cdot 
        e^{\text{slope}(1/T - 1/T_{ref})}(\text{slope}^2 + 6\cdot\text{slope}T + 6T^2)}{T^6}

    Parameters
    ----------
    T : float
        Temperature to evaluate derivative at [K]
    T_ref : float
        Reference temperature point [K]
    P_ref : float
        Pressure at reference point [Pa]
    slope : float
        Slope d(ln P)/d(1/T) at reference point

    Returns
    -------
    d3P_dT3 : float
        Third temperature derivative of vapor pressure [Pa/K^3]

    Notes
    -----
    This derivative provides additional detail for numerical solutions requiring 
    higher-order derivatives and for analyzing the rate of change of vapor pressure 
    curvature.

    Examples
    --------
    >>> d3Arrhenius_extrapolation_dT3(300.0, 400.0, 1E5, -1600)
    0.012727872786932212
    '''
    e = exp(slope*(1.0/T - 1.0/T_ref))
    return -P_ref*slope*e*(slope*slope + 6.0*slope*T + 6.0*T*T)/(T*T*T*T*T*T)