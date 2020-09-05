# -*- coding: utf-8 -*-
r"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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

This module contains lookup functions for melting and boiling point, heat of 
fusion, various enthalpy of vaporization estimation routines, and dataframes
of fit coefficients.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Boiling Point
-------------
.. autofunction:: chemicals.phase_change.Tb
.. autofunction:: chemicals.phase_change.Tb_methods
.. autodata:: chemicals.phase_change.Tb_all_methods

Melting Point
-------------
.. autofunction:: chemicals.phase_change.Tm
.. autofunction:: chemicals.phase_change.Tm_methods
.. autodata:: chemicals.phase_change.Tm_all_methods

Heat of Fusion
--------------
Heat of fusion does not strongly depend on temperature or pressure. This is the
standard value, at 1 atm and the normal melting point.

.. autofunction:: chemicals.phase_change.Hfus
.. autofunction:: chemicals.phase_change.Hfus_methods
.. autodata:: chemicals.phase_change.Hfus_all_methods

Heat of Vaporization at Tb Correlations 
---------------------------------------
.. autofunction:: chemicals.phase_change.Riedel
.. autofunction:: chemicals.phase_change.Chen
.. autofunction:: chemicals.phase_change.Liu
.. autofunction:: chemicals.phase_change.Vetere

Heat of Vaporization at T Correlations 
--------------------------------------
.. autofunction:: chemicals.phase_change.Pitzer
.. autofunction:: chemicals.phase_change.SMK
.. autofunction:: chemicals.phase_change.MK
.. autofunction:: chemicals.phase_change.Velasco
.. autofunction:: chemicals.phase_change.Pitzer
.. autofunction:: chemicals.phase_change.Clapeyron
.. autofunction:: chemicals.phase_change.Watson

Heat of Vaporization at T Model Equations 
-----------------------------------------
.. autofunction:: chemicals.phase_change.Alibakhshi
.. autofunction:: chemicals.phase_change.PPDS12

Heat of Sublimation 
-------------------
No specific correlation is provided. This value is fairly strongly temperature
dependent; the dependency comes almost entirely from the vaporization
enthalpy's dependence. To calculate heat of sublimation at any temperature, use
the equation :math:`H_{sub} = H_{fus} + H_{vap}`.

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: phase_change_data_Perrys2_150

    A collection of 344 coefficient sets from the DIPPR database published
    openly in [1]_. Provides temperature limits for all its fluids. 
    See :obj:`chemicals.dippr.EQ106` for the model equation.

.. data:: phase_change_data_VDI_PPDS_4

    Coefficients for a equation form developed by the PPDS, published 
    openly in [2]_. Extrapolates poorly at low temperatures. See :obj:`PPDS12` 
    for the model equation.

.. data:: phase_change_data_Alibakhshi_Cs

    One-constant limited temperature range regression coefficients presented
    in [3]_, with constants for ~2000 chemicals from the DIPPR database.
    Valid up to 100 K below the critical point, and 50 K under the boiling
    point. See :obj:`Alibakhshi` for the model equation.

.. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    8E. McGraw-Hill Professional, 2007.
.. [2] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
    Berlin; New York:: Springer, 2010.
.. [3] Alibakhshi, Amin. "Enthalpy of Vaporization, Its Temperature
    Dependence and Correlation with Surface Tension: A Theoretical Approach."
    Fluid Phase Equilibria 432 (January 25, 2017): 62-69.
    https://doi.org/10.1016/j.fluid.2016.10.013.

The structure of each dataframe is shown below:
    

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.phase_change.phase_change_data_Perrys2_150

    In [3]: chemicals.phase_change.phase_change_data_VDI_PPDS_4

    In [4]: chemicals.phase_change.phase_change_data_Alibakhshi_Cs

"""

__all__ = ['Tb_methods', 'Tb', 'Tm_methods', 'Tm', 
           'Clapeyron', 'Pitzer', 'SMK', 'MK', 'Velasco', 'Riedel', 'Chen', 
           'Liu', 'Vetere', 'Alibakhshi','PPDS12', 'Watson', 'Hfus', 'Hfus_methods']

import os
from fluids.numerics import numpy as np
from fluids.constants import R, N_A, pi
from chemicals.utils import log
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals import miscdata
from chemicals.data_reader import (register_df_source,
                                   data_source,
                                   retrieve_from_df_dict,
                                   retrieve_any_from_df_dict,
                                   list_available_methods_from_df_dict)

###  Register data sources and lazy load them

folder = os_path_join(source_path, 'Phase Change')
register_df_source(folder, 'Yaws Boiling Points.tsv')
register_df_source(folder, 'OpenNotebook Melting Points.tsv')
register_df_source(folder, 'Ghazerati Appendix Vaporization Enthalpy.tsv',
                   csv_kwargs={'dtype': {'Hvap298': float}})
register_df_source(folder, 'CRC Handbook Heat of Vaporization.tsv')
register_df_source(folder, 'CRC Handbook Heat of Fusion.tsv')
register_df_source(folder, 'Ghazerati Appendix Sublimation Enthalpy.tsv')
register_df_source(folder, 'Table 2-150 Heats of Vaporization of Inorganic and Organic Liquids.tsv')
register_df_source(folder, 'VDI PPDS Enthalpies of vaporization.tsv')
register_df_source(folder, 'Alibakhshi one-coefficient enthalpy of vaporization.tsv')

CRC_ORG = 'CRC_ORG'
CRC_INORG = 'CRC_INORG'
YAWS = 'YAWS'
OPEN_NTBKM = 'OPEN_NTBKM'
CRC = 'CRC'

_phase_change_const_loaded = False
def _load_phase_change_constants():
    global Tb_data_Yaws, Tm_ON_data, Hvap_data_Gharagheizi, Hvap_data_CRC
    global Hfus_data_CRC, Hsub_data_Gharagheizi, _phase_change_const_loaded
    global Tb_sources, Tm_sources, Hfus_sources
    Tb_data_Yaws = data_source('Yaws Boiling Points.tsv')
    Tm_ON_data = data_source('OpenNotebook Melting Points.tsv')
    Hvap_data_Gharagheizi = data_source('Ghazerati Appendix Vaporization Enthalpy.tsv')
    Hvap_data_CRC = data_source('CRC Handbook Heat of Vaporization.tsv')
    Hfus_data_CRC = data_source('CRC Handbook Heat of Fusion.tsv')
    Hsub_data_Gharagheizi = data_source('Ghazerati Appendix Sublimation Enthalpy.tsv')
    _phase_change_const_loaded = True
    Tb_sources = {
        CRC_ORG: miscdata.CRC_organic_data,
        CRC_INORG: miscdata.CRC_inorganic_data,
        YAWS: Tb_data_Yaws,
    }
    Tm_sources = {
        OPEN_NTBKM: Tm_ON_data,
        CRC_INORG: miscdata.CRC_inorganic_data,
        CRC_ORG: miscdata.CRC_organic_data,
    }
    Hfus_sources = {
        CRC: Hfus_data_CRC,
    }
    

_phase_change_corrs_loaded = False
def _load_phase_change_correlations():
    global phase_change_data_Perrys2_150, phase_change_values_Perrys2_150
    global phase_change_data_VDI_PPDS_4, phase_change_values_VDI_PPDS_4
    global phase_change_data_Alibakhshi_Cs, _phase_change_corrs_loaded

    # 66554 for pandas; 19264 bytes for numpy
    phase_change_data_Perrys2_150 = data_source('Table 2-150 Heats of Vaporization of Inorganic and Organic Liquids.tsv')
    phase_change_values_Perrys2_150 = np.array(phase_change_data_Perrys2_150.values[:, 1:], dtype=float)
    
    # 52187 bytes for pandas, 13056 bytes for numpy
    phase_change_data_VDI_PPDS_4 = data_source('VDI PPDS Enthalpies of vaporization.tsv')
    phase_change_values_VDI_PPDS_4 = np.array(phase_change_data_VDI_PPDS_4.values[:, 2:], dtype=float)
    
    phase_change_data_Alibakhshi_Cs = data_source('Alibakhshi one-coefficient enthalpy of vaporization.tsv')
    _phase_change_corrs_loaded = True
    
if PY37:
    def __getattr__(name):
        if name in ('Tb_data_Yaws', 'Tm_ON_data', 'Hvap_data_Gharagheizi',
                    'Hvap_data_CRC', 'Hfus_data_CRC', 'Hsub_data_Gharagheizi',
                    'Tb_sources', 'Tm_sources', 'Hfus_sources'):
            _load_phase_change_constants()
            return globals()[name]
        elif name in ('phase_change_data_Perrys2_150',
                      'phase_change_values_Perrys2_150', 
                      'phase_change_data_VDI_PPDS_4', 
                      'phase_change_values_VDI_PPDS_4', 
                      'phase_change_data_Alibakhshi_Cs'):
            _load_phase_change_correlations()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_phase_change_constants()
        _load_phase_change_correlations()

### Phase change functions

### Boiling Point at 1 atm

Tb_all_methods = (CRC_INORG, CRC_ORG, YAWS)
'''Tuple of method name keys. See the `Tb` for the actual references'''

def Tb_methods(CASRN):
    """Return all methods available to obtain the Tb for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Tb with the given inputs.

    See Also
    --------
    Tb
    """
    if not _phase_change_const_loaded: _load_phase_change_constants()
    return list_available_methods_from_df_dict(Tb_sources, CASRN, 'Tb')

def Tb(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's boiling
    point. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Preferred sources are 'CRC Physical Constants, organic' for organic
    chemicals, and 'CRC Physical Constants, inorganic' for inorganic
    chemicals. Function has data for approximately 13000 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Tb : float
        Boiling temperature, [K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined in the variable,
        `Tb_all_methods`.

    Notes
    -----
    A total of three methods are available for this function. They are:

        * 'CRC_ORG', a compillation of data on organics
          as published in [1]_.
        * 'CRC_INORG', a compillation of data on
          inorganic as published in [1]_.
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [2]_.

    Examples
    --------
    >>> Tb('7732-18-5')
    373.124

    See Also
    --------
    Tb_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [2] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    if not _phase_change_const_loaded: _load_phase_change_constants()
    if method:
        return retrieve_from_df_dict(Tb_sources, CASRN, 'Tb', method) 
    else:
        return retrieve_any_from_df_dict(Tb_sources, CASRN, 'Tb') 

### Melting Point

Tm_all_methods = (OPEN_NTBKM, CRC_INORG, CRC_ORG)
'''Tuple of method name keys. See the `Tm` for the actual references'''

def Tm_methods(CASRN):
    """Return all methods available to obtain the Tm for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Tm with the given inputs.

    See Also
    --------
    Tm
    """
    if not _phase_change_const_loaded: _load_phase_change_constants()
    return list_available_methods_from_df_dict(Tm_sources, CASRN, 'Tm')

def Tm(CASRN, method=None):
    r'''This function handles the retrieval of a chemical's melting
    point. Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    Preferred sources are 'Open Notebook Melting Points', with backup sources
    'CRC Physical Constants, organic' for organic chemicals, and
    'CRC Physical Constants, inorganic' for inorganic chemicals. Function has
    data for approximately 14000 chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Tm : float
        Melting temperature, [K]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by the vairable
        `Tm_all_methods`.

    Notes
    -----
    A total of three sources are available for this function. They are:

        * 'OPEN_NTBKM, a compillation of data on organics
          as published in [1]_ as Open Notebook Melting Points; Averaged 
          (median) values were used when
          multiple points were available. For more information on this
          invaluable and excellent collection, see
          http://onswebservices.wikispaces.com/meltingpoint.
        * 'CRC_ORG', a compillation of data on organics
          as published in [2]_.
        * 'CRC_INORG', a compillation of data on
          inorganic as published in [2]_.

    Examples
    --------
    >>> Tm(CASRN='7732-18-5')
    273.15

    See Also
    --------
    Tm_methods

    References
    ----------
    .. [1] Bradley, Jean-Claude, Antony Williams, and Andrew Lang.
       "Jean-Claude Bradley Open Melting Point Dataset", May 20, 2014.
       https://figshare.com/articles/Jean_Claude_Bradley_Open_Melting_Point_Datset/1031637.
    .. [2] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    if not _phase_change_const_loaded: _load_phase_change_constants()
    elif method:
        return retrieve_from_df_dict(Tm_sources, CASRN, 'Tm', method) 
    else:
        return retrieve_any_from_df_dict(Tm_sources, CASRN, 'Tm') 

### Enthalpy of Vaporization at T

def Clapeyron(T, Tc, Pc, dZ=1, Psat=101325):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using the
    Clapeyron equation.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} = RT \Delta Z \frac{\ln (P_c/Psat)}{(1-T_{r})}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    dZ : float
        Change in compressibility factor between liquid and gas, []
    Psat : float
        Saturation pressure of fluid [Pa], optional

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    No original source is available for this equation.
    [1]_ claims this equation overpredicts enthalpy by several percent.
    Under Tr = 0.8, dZ = 1 is a reasonable assumption.
    This equation is most accurate at the normal boiling point.

    Internal units are bar.

    WARNING: I believe it possible that the adjustment for pressure may be incorrect

    Examples
    --------
    Problem from Perry's examples.

    >>> Clapeyron(T=294.0, Tc=466.0, Pc=5.55E6)
    26512.36357131963

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    return R*T*dZ*log(Pc/Psat)/(1. - Tr)

def Pitzer(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    fit by [2]_ to the work of Pitzer [1]_; requires a chemical's critical
    temperature and acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \frac{\Delta_{vap} H}{RT_c}=7.08(1-T_r)^{0.354}+10.95\omega(1-T_r)^{0.456}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    This equation is listed in [3]_, page 2-487 as method #2 for estimating
    Hvap. This cites [2]_.

    The recommended range is 0.6 to 1 Tr. Users should expect up to 5% error.
    T must be under Tc, or an exception is raised.

    The original article has been reviewed and found to have a set of tabulated
    values which could be used instead of the fit function to provide additional
    accuracy.

    Examples
    --------
    Example as in [3]_, p2-487; exp: 37.51 kJ/mol

    >>> Pitzer(452, 645.6, 0.35017)
    36696.749078320056

    References
    ----------
    .. [1] Pitzer, Kenneth S. "The Volumetric and Thermodynamic Properties of
       Fluids. I. Theoretical Basis and Virial Coefficients."
       Journal of the American Chemical Society 77, no. 13 (July 1, 1955):
       3427-33. doi:10.1021/ja01618a001
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [3] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Tr = T/Tc
    return R*Tc * (7.08*(1. - Tr)**0.354 + 10.95*omega*(1. - Tr)**0.456)

def SMK(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
         \frac{\Delta H_{vap}} {RT_c} =
         \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)} + \left(
         \frac{\omega - \omega^{(R1)}} {\omega^{(R2)} - \omega^{(R1)}} \right)
         \left[\left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R2)} - \left(
         \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)} \right]

    .. math::
        \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)}
        = 6.537 \tau^{1/3} - 2.467 \tau^{5/6} - 77.251 \tau^{1.208} +
        59.634 \tau + 36.009 \tau^2 - 14.606 \tau^3

    .. math::
        \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R2)} - \left(
        \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)}=-0.133 \tau^{1/3} - 28.215
        \tau^{5/6} - 82.958 \tau^{1.208} + 99.00 \tau  + 19.105 \tau^2 -2.796 \tau^3

    .. math::
        \tau = 1-T/T_c

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed and found to have coefficients with
    slightly more precision. Additionally, the form of the equation is slightly
    different, but numerically equivalent.

    The refence fluids are:

    :math:`\omega_0` = benzene = 0.212

    :math:`\omega_1` = carbazole = 0.461

    A sample problem in the article has been verified. The numerical result
    presented by the author requires high numerical accuracy to obtain.

    Examples
    --------
    Problem in [1]_:

    >>> SMK(553.15, 751.35, 0.302)
    39866.18999046229

    References
    ----------
    .. [1] Sivaraman, Alwarappa, Joe W. Magee, and Riki Kobayashi. "Generalized
       Correlation of Latent Heats of Vaporization of Coal-Liquid Model Compounds
       between Their Freezing Points and Critical Points." Industrial &
       Engineering Chemistry Fundamentals 23, no. 1 (February 1, 1984): 97-100.
       doi:10.1021/i100013a017.
    '''
    if T > Tc:
        return 0.0
    omegaR1, omegaR2 = 0.212, 0.461
    A10 = 6.536924
    A20 = -2.466698
    A30 = -77.52141
    B10 = 59.63435
    B20 = 36.09887
    B30 = -14.60567

    A11 = -0.132584
    A21 = -28.21525
    A31 = -82.95820
    B11 = 99.00008
    B21 = 19.10458
    B31 = -2.795660
    # 1 branch, two powers, 1 division
    tau = 1. - T/Tc
    tau_16 = tau**(1.0/6.0)
    tau_13 = tau_16*tau_16
    tau_56 = tau_13*tau_13*tau_16
    tau_other = tau**1.2083333333333333 #(1-1/8. + 1/3.))
    L0 = (A10*tau_13 + A20*tau_56 + A30*tau_other +
          tau*(B10 + tau*(B20 + B30*tau)))

    L1 = (A11*tau_13 + A21*tau_56 + A31*tau_other + 
           tau*(B11 + tau*(B21 + B31*tau)))

    domega = 4.016064257028112*(omega - omegaR1) # 4.016... = 1.0/(omegaR2 - omegaR1
#    domega = (omega - omegaR1)/(omegaR2 - omegaR1)
    return R*Tc*(L0 + domega*L1)

def MK(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} =  \Delta H_{vap}^{(0)} + \omega \Delta H_{vap}^{(1)} + \omega^2 \Delta H_{vap}^{(2)}

    .. math::
        \frac{\Delta H_{vap}^{(i)}}{RT_c} = b^{(j)} \tau^{1/3} + b_2^{(j)} \tau^{5/6}
        + b_3^{(j)} \tau^{1.2083} + b_4^{(j)}\tau + b_5^{(j)} \tau^2 + b_6^{(j)} \tau^3

    .. math::
        \tau = 1-T/T_c

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed. A total of 18 coefficients are used:

    WARNING: The correlation has been implemented as described in the article,
    but its results seem different and with some error.
    Its results match with other functions however.

    Has poor behavior for low-temperature use.

    Examples
    --------
    Problem in article for SMK function.

    >>> MK(553.15, 751.35, 0.302)
    38728.00667307733

    References
    ----------
    .. [1] Morgan, David L., and Riki Kobayashi. "Extension of Pitzer CSP
       Models for Vapor Pressures and Heats of Vaporization to Long-Chain
       Hydrocarbons." Fluid Phase Equilibria 94 (March 15, 1994): 51-87.
       doi:10.1016/0378-3812(94)87051-9.
    '''
    bs0 = [5.2804, 0.080022, 7.2543]
    bs1 = [12.8650, 273.23, -346.45]
    bs2 = [1.1710, 465.08, -610.48]
    bs3 = [-13.1160, -638.51, 839.89]
    bs4 = [0.4858, -145.12, 160.05]
    bs5 = [-1.0880, 74.049, -50.711]

    tau = 1. - T/Tc
    tau_third = tau**0.3333
    tau_83 = tau**0.8333
    tau_other = tau**1.2083
    tau2 = tau*tau
    tau3 = tau2*tau
    
    H0 = (bs0[0]*tau_third + bs1[0]*tau_83 + bs2[0]*tau_other +
    bs3[0]*tau + bs4[0]*tau2 + bs5[0]*tau3)

    H1 = (bs0[1]*tau_third + bs1[1]*tau_83 + bs2[1]*tau_other +
    bs3[1]*tau + bs4[1]*tau2 + bs5[1]*tau3)

    H2 = (bs0[2]*tau_third + bs1[2]*tau_83 + bs2[2]*tau_other +
    bs3[2]*tau + bs4[2]*tau2 + bs5[2]*tau3)

    return (H0 + omega*(H1 + omega*H2))*R*Tc

def Velasco(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta_{vap} H = RT_c(7.2729 + 10.4962\omega + 0.6061\omega^2)(1-T_r)^{0.38}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed. It is regressed from enthalpy of
    vaporization values at 0.7Tr, from 121 fluids in REFPROP 9.1.
    A value in the article was read to be similar, but slightly too low from
    that calculated here.

    Examples
    --------
    From graph, in [1]_ for perfluoro-n-heptane.

    >>> Velasco(333.2, 476.0, 0.5559)
    33299.428636069264

    References
    ----------
    .. [1] Velasco, S., M. J. Santos, and J. A. White. "Extended Corresponding
       States Expressions for the Changes in Enthalpy, Compressibility Factor
       and Constant-Volume Heat Capacity at Vaporization." The Journal of
       Chemical Thermodynamics 85 (June 2015): 68-76.
       doi:10.1016/j.jct.2015.01.011.
    '''
    return (7.2729 + 10.4962*omega + 0.6061*omega**2)*(1-T/Tc)**0.38*R*Tc

### Enthalpy of Vaporization at Normal Boiling Point.

def Riedel(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization at the boiling point, using the
    Ridel [1]_ CSP method. Required information are critical temperature
    and pressure, and boiling point. Equation taken from [2]_ and [3]_.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta_{vap} H=1.093 T_b R\frac{\ln P_c-1.013}{0.930-T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at the normal boiling point, [J/mol]

    Notes
    -----
    This equation has no example calculation in any source. The source has not
    been verified. It is equation 4-144 in Perry's. Perry's also claims that
    errors seldom surpass 5%.

    [2]_ is the source of example work here, showing a calculation at 0.0%
    error.

    Internal units of pressure are bar.

    Examples
    --------
    Pyridine, 0.0% err vs. exp: 35090 J/mol; from Poling [2]_.

    >>> Riedel(388.4, 620.0, 56.3E5)
    35089.80179000598

    References
    ----------
    .. [1] Riedel, L. "Eine Neue Universelle Dampfdruckformel Untersuchungen
       Uber Eine Erweiterung Des Theorems Der Ubereinstimmenden Zustande. Teil
       I." Chemie Ingenieur Technik 26, no. 2 (February 1, 1954): 83-89.
       doi:10.1002/cite.330260206.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [3] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Pc = Pc/1E5  # Pa to bar
    Tbr = Tb/Tc
    return 1.093*Tb*R*(log(Pc) - 1.013)/(0.93 - Tbr)

def Chen(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization using the Chen [1]_ correlation
    and a chemical's critical temperature, pressure and boiling point.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vb} = RT_b \frac{3.978 T_r - 3.958 + 1.555 \ln P_c}{1.07 - T_r}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The formulation presented in the original article is similar, but uses
    units of atm and calorie instead. The form in [2]_ has adjusted for this.
    A method for estimating enthalpy of vaporization at other conditions
    has also been developed, but the article is unclear on its implementation.
    Based on the Pitzer correlation.

    Internal units: bar and K

    Examples
    --------
    Same problem as in Perry's examples.

    >>> Chen(294.0, 466.0, 5.55E6)
    26705.902558030946

    References
    ----------
    .. [1] Chen, N. H. "Generalized Correlation for Latent Heat of Vaporization."
       Journal of Chemical & Engineering Data 10, no. 2 (April 1, 1965): 207-10.
       doi:10.1021/je60025a047
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tbr = Tb/Tc
    Pc = Pc/1E5  # Pa to bar
    return R*Tb*(3.978*Tbr - 3.958 + 1.555*log(Pc))/(1.07 - Tbr)

def Liu(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization at the normal boiling point using
    the Liu [1]_ correlation, and a chemical's critical temperature, pressure
    and boiling point.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} = RT_b \left[ \frac{T_b}{220}\right]^{0.0627} \frac{
        (1-T_{br})^{0.38} \ln(P_c/P_A)}{1-T_{br} + 0.38 T_{br} \ln T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    This formulation can be adjusted for lower boiling points, due to the use
    of a rationalized pressure relationship. The formulation is taken from
    the original article.

    A correction for alcohols and organic acids based on carbon number,
    which only modifies the boiling point, is available but not implemented.

    No sample calculations are available in the article.

    Internal units: Pa and K

    Examples
    --------
    Same problem as in Perry's examples

    >>> Liu(294.0, 466.0, 5.55E6)
    26378.575260517395

    References
    ----------
    .. [1] LIU, ZHI-YONG. "Estimation of Heat of Vaporization of Pure Liquid at
       Its Normal Boiling Temperature." Chemical Engineering Communications
       184, no. 1 (February 1, 2001): 221-28. doi:10.1080/00986440108912849.
    '''
    Tbr = Tb/Tc
    return R*Tb*(Tb/220.)**0.0627*(1. - Tbr)**0.38*log(Pc/101325.) \
        / (1 - Tbr + 0.38*Tbr*log(Tbr))

def Vetere(Tb, Tc, Pc, F=1):
    r'''Calculates enthalpy of vaporization at the boiling point, using the
    Vetere [1]_ CSP method. Required information are critical temperature
    and pressure, and boiling point. Equation taken from [2]_.

    The enthalpy of vaporization is given by:

    .. math::
        \frac {\Delta H_{vap}}{RT_b} = \frac{\tau_b^{0.38}
        \left[ \ln P_c - 0.513 + \frac{0.5066}{P_cT_{br}^2}\right]}
        {\tau_b + F(1-\tau_b^{0.38})\ln T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    F : float, optional
        Constant for a fluid, [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at the boiling point, [J/mol]

    Notes
    -----
    The equation cannot be found in the original source. It is believed that a
    second article is its source, or that DIPPR staff have altered the formulation.

    Internal units of pressure are bar.

    Examples
    --------
    Example as in [2]_, p2-487; exp: 25.73

    >>> Vetere(294.0, 466.0, 5.55E6)
    26363.43895706672

    References
    ----------
    .. [1] Vetere, Alessandro. "Methods to Predict the Vaporization Enthalpies
       at the Normal Boiling Temperature of Pure Compounds Revisited."
       Fluid Phase Equilibria 106, no. 1-2 (May 1, 1995): 1–10.
       doi:10.1016/0378-3812(94)02627-D.
    .. [2] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Tbr = Tb/Tc
    taub = 1-Tb/Tc
    Pc = Pc/1E5
    term = taub**0.38*(log(Pc)-0.513 + 0.5066/Pc/Tbr**2) / (taub + F*(1-taub**0.38)*log(Tbr))
    return R*Tb*term

### Enthalpy of Vaporization adjusted for T

def Watson(T, Hvap_ref, T_ref, Tc, exponent=0.38):
    r'''Calculates enthalpy of vaporization of a chemical at a temperature 
    using the known heat of vaporization at another temperature according to
    the Watson [1]_ [2]_ correlation. This is an application of the 
    corresponding-states principle, with an emperical temperature dependence.

    .. math::
        \frac{\Delta H_{vap}^{T2}}{\Delta H_{vap}^{T1}}  = \left(
        \frac{1-T_{r,1}}{1-T_{r,2}} \right)^{0.38}

    Parameters
    ----------
    T : float
        Temperature for which to calculate heat of vaporization, [K]
    Hvap_ref : float
        Enthalpy of vaporization at the known temperature point, [J/mol]
    T_ref : float
        Reference temperature; ideally as close to `T` as posible, [K]
    Tc : float
        Critical temperature of fluid [K]
    exponent : float, optional
        A fit exponent can optionally be used instead of the Watson 0.38 
        exponent, [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at `T`, [J/mol]

    Notes
    -----

    Examples
    --------
    Predict the enthalpy of vaporization of water at 320 K from a point at
    300 K:
        
    >>> Watson(T=320, Hvap_ref=43908, T_ref=300.0, Tc=647.14)
    42928.990094915454
    
    The error is 0.38% compared to the correct value of 43048 J/mol.

    References
    ----------
    .. [1] Watson, KM. "Thermodynamics of the Liquid State." Industrial & 
       Engineering Chemistry 35, no. 4 (1943): 398–406.
    .. [2] Martin, Joseph J., and John B. Edwards. "Correlation of Latent Heats
       of Vaporization.” AIChE Journal 11, no. 2 (1965): 331-33. 
       https://doi.org/10.1002/aic.690110226.
    '''
    Tr = T/Tc
    Trefr = T_ref/Tc
    H2 = Hvap_ref*((1.0 - Tr)/(1.0 - Trefr))**exponent
    return H2

### Enthalpy of Vaporization model equations
    
def Alibakhshi(T, Tc, C):
    r'''Calculates enthalpy of vaporization of a chemical at a temperature 
    using a theoretically-derived single-coefficient fit equation developed in
    [1]_. This model falls apart at ~0.8 Tc.

    .. math::
        \Delta H_{vap} = \left(4.5\pi N_A\right)^{1/3.}4.2\times 10^{-7}
        (T_c - 6) - 0.5RT\log(T) + CT

    Parameters
    ----------
    T : float
        Temperature for which to calculate heat of vaporization, [K]
    Tc : float
        Critical temperature of fluid [K]
    C : float
        Alibakhshi fit coefficient, [J/mol/K]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at `T`, [J/mol]

    Notes
    -----
    The authors of [1]_ evaluated their model on 1890 compounds for a 
    temperature range of 50 K under `Tb` to 100 K below `Tc`, and obtained an
    average absolute relative error of 4.5%.

    Examples
    --------
    Predict the enthalpy of vaporization of water at 320 K:
        
    >>> Alibakhshi(T=320.0, Tc=647.14, C=-16.7171)
    41961.30490225752
    
    The error is 2.5% compared to the correct value of 43048 J/mol.

    References
    ----------
    .. [1] Alibakhshi, Amin. "Enthalpy of Vaporization, Its Temperature
       Dependence and Correlation with Surface Tension: A Theoretical Approach."
       Fluid Phase Equilibria 432 (January 25, 2017): 62-69.
       https://doi.org/10.1016/j.fluid.2016.10.013.
    '''
    return (4.5*pi*N_A)**(1/3.)*4.2E-7*(Tc-6.) - R/2.*T*log(T) + C*T

def PPDS12(T, Tc, A, B, C, D, E):
    r'''Calculate the enthalpy of vaporization of a fluid using the 5-term 
    power fit developed by the PPDS and named PPDS equation 12.
    
    .. math::
       H_{vap} = RT_c \left(A\tau^{1/3} + B\tau^{2/3} + C\tau + D\tau^2 
       + E\tau^6\right)
    
    .. math::
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [-]
    C : float
        Coefficient, [-]
    D : float
        Coefficient, [-]
    E : float
        Coefficient, [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at `T`, [J/mol]

    Notes
    -----
    Coefficients can be found in [1]_, but no other source for these 
    coefficients has been found.

    Examples
    --------
    >>> PPDS12(300.0, 591.75, 4.60584, 13.97224, -10.592315, 2.120205, 4.277128)
    37948.76862035925
    
    References
    ----------
    .. [1] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
       Berlin; New York:: Springer, 2010.
    '''
    tau = 1. - T/Tc
    tau_cbrt = tau**(1/3.)
    tau2 = tau*tau
    Hvap = R*Tc*(tau_cbrt*(A + B*tau_cbrt) + C*tau
                               + tau2*(D + E*tau2*tau2))
    return Hvap

### Heat of Fusion

Hfus_all_methods = (CRC,)
'''Tuple of method name keys. See the `Hfus` for the actual references'''

def Hfus_methods(CASRN):
    """Return all methods available to obtain the Hfus for the desired chemical.

    Parameters
    ----------
    CASRN : str
        CASRN, [-]

    Returns
    -------
    methods : list[str]
        Methods which can be used to obtain the Hfus with the given inputs.

    See Also
    --------
    Hfus
    """
    if not _phase_change_const_loaded: _load_phase_change_constants()
    return list_available_methods_from_df_dict(Hfus_sources, CASRN, 'Hfus')

def Hfus(CASRN, method=None): 
    r'''This function handles the retrieval of a chemical's heat of fusion.
    Lookup is based on CASRNs. Will automatically select a data
    source to use if no method is provided; returns None if the data is not
    available.

    The Preferred source is 'CRC'. Function has data for approximately 1100 
    chemicals.

    Parameters
    ----------
    CASRN : str
        CASRN [-]

    Returns
    -------
    Hfus : float
        Molar enthalpy of fusion at normal melting point, [J/mol]

    Other Parameters
    ----------------
    method : string, optional
        A string for the method name to use, as defined by the variable,
        `Hfus_all_methods`.

    Notes
    -----
    A total of one method is available for this function. They are:

        * 'CRC', a compillation of data on organics and inorganics as published 
        in [1]_.

    Examples
    --------
    >>> Hfus('7732-18-5')
    6010.0

    See Also
    --------
    Hfus_methods

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    if not _phase_change_const_loaded: _load_phase_change_constants()
    if method:
        return retrieve_from_df_dict(Hfus_sources, CASRN, 'Hfus', method) 
    else:
        return retrieve_any_from_df_dict(Hfus_sources, CASRN, 'Hfus') 


