# -*- coding: utf-8 -*-
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

This module contains various surface tension estimation routines, dataframes
of fit coefficients, fitting model equations, mixing rules, and
water-hydrocarbon interfacial tension estimation routines.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Pure Component Correlations
---------------------------
.. autofunction:: chemicals.interface.Brock_Bird
.. autofunction:: chemicals.interface.Pitzer_sigma
.. autofunction:: chemicals.interface.Sastri_Rao
.. autofunction:: chemicals.interface.Zuo_Stenby
.. autofunction:: chemicals.interface.Hakim_Steinberg_Stiel
.. autofunction:: chemicals.interface.Miqueu
.. autofunction:: chemicals.interface.Aleem
.. autofunction:: chemicals.interface.Mersmann_Kind_sigma

Mixing Rules
------------
.. autofunction:: chemicals.interface.Winterfeld_Scriven_Davis
.. autofunction:: chemicals.interface.Weinaug_Katz
.. autofunction:: chemicals.interface.Diguilio_Teja

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.interface.sigma_IAPWS

Petroleum Correlations
----------------------
.. autofunction:: chemicals.interface.API10A32

Oil-Water Interfacial Tension Correlations
------------------------------------------
.. autofunction:: chemicals.interface.Meybodi_Daryasafar_Karimi

Fit Correlations
----------------
.. autofunction:: chemicals.interface.REFPROP_sigma
.. autofunction:: chemicals.interface.Somayajulu
.. autofunction:: chemicals.interface.Jasper
.. autofunction:: chemicals.interface.PPDS14
.. autofunction:: chemicals.interface.Watson_sigma
.. autofunction:: chemicals.interface.ISTExpansion

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: sigma_data_Mulero_Cachadina

    Data from [5]_ with :obj:`REFPROP_sigma` coefficients.

.. data:: sigma_data_Jasper_Lange

    Data as shown in [4]_ but originally in [3]_ with :obj:`Jasper` coefficients.

.. data:: sigma_data_Somayajulu

    Data from [1]_ with :obj:`Somayajulu` coefficients.

.. data:: sigma_data_Somayajulu2

    Data from [2]_ with :obj:`Somayajulu` coefficients. These should be
    preferred over the original coefficients.

.. data:: sigma_data_VDI_PPDS_11

    Data from [6]_ with :obj:`chemicals.dippr.EQ106` coefficients.

.. [1] Somayajulu, G. R. "A Generalized Equation for Surface Tension from
   the Triple Point to the Critical Point." International Journal of
   Thermophysics 9, no. 4 (July 1988): 559-66. doi:10.1007/BF00503154.
.. [2] Mulero, A., M. I. Parra, and I. Cachadina. "The Somayajulu
   Correlation for the Surface Tension Revisited." Fluid Phase
   Equilibria 339 (February 15, 2013): 81-88.
   doi:10.1016/j.fluid.2012.11.038.
.. [3] Jasper, Joseph J. "The Surface Tension of Pure Liquid Compounds."
   Journal of Physical and Chemical Reference Data 1, no. 4
   (October 1, 1972): 841-1010. doi:10.1063/1.3253106.
.. [4] Speight, James. Lange's Handbook of Chemistry. 16 edition.
   McGraw-Hill Professional, 2005.
.. [5] Mulero, A., I. Cachadiña, and M. I. Parra. “Recommended
   Correlations for the Surface Tension of Common Fluids.” Journal of
   Physical and Chemical Reference Data 41, no. 4 (December 1, 2012):
   043105. doi:10.1063/1.4768782.
.. [6] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
   Berlin; New York:: Springer, 2010.

The structure of each dataframe is shown below:


.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.interface.sigma_data_Mulero_Cachadina

    In [3]: chemicals.interface.sigma_data_Jasper_Lange

    In [4]: chemicals.interface.sigma_data_Somayajulu

    In [5]: chemicals.interface.sigma_data_Somayajulu2

    In [6]: chemicals.interface.sigma_data_VDI_PPDS_11
"""


from __future__ import division

__all__ = ['REFPROP_sigma', 'Somayajulu', 'Jasper',
           'Brock_Bird', 'Pitzer_sigma', 'Sastri_Rao', 'Zuo_Stenby',
           'sigma_IAPWS', 'PPDS14', 'Watson_sigma',
           'Mersmann_Kind_sigma', 'API10A32',
           'Hakim_Steinberg_Stiel', 'Miqueu', 'Aleem',
           'Winterfeld_Scriven_Davis', 'Diguilio_Teja', 'Weinaug_Katz',
           'Meybodi_Daryasafar_Karimi', 'ISTExpansion']

import os
from fluids.numerics import numpy as np
from fluids.constants import N_A, k
from chemicals.utils import log, exp, sqrt
from chemicals.utils import mixing_simple, PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import register_df_source, data_source

folder = os_path_join(source_path, 'Interface')


register_df_source(folder, 'MuleroCachadinaParameters.tsv')
register_df_source(folder, 'Jasper-Lange.tsv')
register_df_source(folder, 'Somayajulu.tsv')
register_df_source(folder, 'SomayajuluRevised.tsv')
register_df_source(folder, 'VDI PPDS surface tensions.tsv')

_interface_dfs_loaded = False
def load_interface_dfs():
    global _interface_dfs_loaded, sigma_data_Mulero_Cachadina, sigma_values_Mulero_Cachadina
    global sigma_data_Jasper_Lange, sigma_values_Jasper_Lange
    global sigma_data_Somayajulu, sigma_values_Somayajulu, sigma_data_Somayajulu2
    global sigma_values_Somayajulu2, sigma_data_VDI_PPDS_11, sigma_values_VDI_PPDS_11

    sigma_data_Mulero_Cachadina = data_source('MuleroCachadinaParameters.tsv')
    sigma_values_Mulero_Cachadina = np.array(sigma_data_Mulero_Cachadina.values[:, 1:], dtype=float)

    sigma_data_Jasper_Lange = data_source('Jasper-Lange.tsv')
    sigma_values_Jasper_Lange = np.array(sigma_data_Jasper_Lange.values[:, 1:], dtype=float)

    sigma_data_Somayajulu = data_source('Somayajulu.tsv')
    sigma_values_Somayajulu = np.array(sigma_data_Somayajulu.values[:, 1:], dtype=float)

    sigma_data_Somayajulu2 = data_source('SomayajuluRevised.tsv')
    sigma_values_Somayajulu2 = np.array(sigma_data_Somayajulu2.values[:, 1:], dtype=float)

    sigma_data_VDI_PPDS_11 = data_source('VDI PPDS surface tensions.tsv')
    sigma_values_VDI_PPDS_11 = np.array(sigma_data_VDI_PPDS_11.values[:, 1:], dtype=float)

if PY37:
    def __getattr__(name):
        if name in ('sigma_data_Mulero_Cachadina', 'sigma_values_Mulero_Cachadina',
                    'sigma_data_Jasper_Lange', 'sigma_values_Jasper_Lange',
                    'sigma_data_Somayajulu', 'sigma_values_Somayajulu', 'sigma_data_Somayajulu2',
                    'sigma_values_Somayajulu2', 'sigma_data_VDI_PPDS_11', 'sigma_values_VDI_PPDS_11'
                    ):
            load_interface_dfs()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        load_interface_dfs()



def sigma_IAPWS(T):
    r'''Calculate the surface tension of pure water as a function of .
    temperature. Assumes the 2011 IAPWS [1]_ formulation.

    .. math::
        \sigma = B\tau^\mu(1+b\tau)\\

    .. math::
        \tau = 1-T/T_c\\

    .. math::
        B = 0.2358 \text{N/m}\\

    .. math::
        b = -0.625\\

    .. math::
        \mu = 1.256

    Parameters
    ----------
    T : float
        Temperature of liquid [K]

    Returns
    -------
    sigma : float
        Air-water surface tension, [N/m]

    Notes
    -----
    This function is valid from the triple temperature to the critical
    temperature. No effects for pressure are included in the formulation.
    Test values are from IAPWS 2010 book.

    Examples
    --------
    >>> sigma_IAPWS(300.)
    0.0716859625271
    >>> sigma_IAPWS(450.)
    0.0428914991565
    >>> sigma_IAPWS(600.)
    0.0083756108728

    References
    ----------
    .. [1] IAPWS. 2014. Revised Release on Surface Tension of Ordinary Water
       Substance
    '''
    tau = 1. - T*(1.0/647.096)
    return 0.2358*tau**1.256*(1.0 - 0.625*tau)

### Regressed coefficient-based functions

def REFPROP_sigma(T, Tc, sigma0, n0, sigma1=0.0, n1=0.0, sigma2=0.0, n2=0.0):
    r'''Calculates air-liquid surface tension  using the REFPROP_sigma [1]_
    regression-based method. Relatively recent, and most accurate.

    .. math::
        \sigma(T)=\sigma_0\left(1-\frac{T}{T_c}\right)^{n_0}+
        \sigma_1\left(1-\frac{T}{T_c}\right)^{n_1}+
        \sigma_2\left(1-\frac{T}{T_c}\right)^{n_2}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    sigma0 : float
        First emperical coefficient of a fluid
    n0 : float
        First emperical exponent of a fluid
    sigma1 : float, optional
        Second emperical coefficient of a fluid.
    n1 : float, optional
        Second emperical exponent of a fluid.
    sigma1 : float, optional
        Third emperical coefficient of a fluid.
    n2 : float, optional
        Third emperical exponent of a fluid.

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Notes
    -----
    Function as implemented in [1]_. No example necessary; results match
    literature values perfectly.
    Form of function returns imaginary results when T > Tc; None is returned
    if this is the case.


    Examples
    --------
    Parameters for water at 298.15 K

    >>> REFPROP_sigma(298.15, 647.096, -0.1306, 2.471, 0.2151, 1.233)
    0.07205503890847453

    References
    ----------
    .. [1] Diky, Vladimir, Robert D. Chirico, Chris D. Muzny, Andrei F.
       Kazakov, Kenneth Kroenlein, Joseph W. Magee, Ilmutdin Abdulagatov, and
       Michael Frenkel. "ThermoData Engine (TDE): Software Implementation of
       the Dynamic Data Evaluation Concept." Journal of Chemical Information
       and Modeling 53, no. 12 (2013): 3418-30. doi:10.1021/ci4005699.
    '''
    Tr = T/Tc
    one_minus_Tr = 1.0 - Tr
    sigma = sigma0*(one_minus_Tr)**n0 + sigma1*(one_minus_Tr)**n1 + sigma2*(one_minus_Tr)**n2
    return sigma


def PPDS14(T, Tc, a0, a1, a2):
    r'''Calculates air-water surface tension  using the [1]_
    emperical (parameter-regressed) method, called the PPDS 14 equation for
    surface tension.
    
    .. math::
        \sigma = a_0 \tau^{a_1}(1 + a_2 \tau)

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    a0 : float
        Regression parameter, [N/m]
    a1 : float
        Regression parameter, [-]
    a2 : float
        Regression parameter, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Notes
    -----

    Examples
    --------
    Benzene at 280 K from [1]_

    >>> PPDS14(T=280, Tc=562.05, a0=0.0786269, a1=1.28646, a2=-0.112304)
    0.030559764256249854

    References
    ----------
    .. [1] "ThermoData Engine (TDE103b V10.1) User’s Guide." 
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-SurfaceTension/PPDS14.htm.
    .. [2] Frenkel, Michael, Robert D. Chirico, Vladimir Diky, Xinjian Yan, 
       Qian Dong, and Chris Muzny. "ThermoData Engine (TDE):  Software
       Implementation of the Dynamic Data Evaluation Concept." Journal of 
       Chemical Information and Modeling 45, no. 4 (July 1, 2005): 816-38. 
       https://doi.org/10.1021/ci050067b.
    '''
    tau = 1.0 - T/Tc
    return a0*tau**a1*(1.0 + a2*tau)

def Watson_sigma(T, Tc, a1, a2, a3=0.0, a4=0.0, a5=0.0):
    r'''Calculates air-water surface tension using the Watson [1]_
    emperical (parameter-regressed) method developed by NIST.
    
    .. math::
        \sigma = \exp\left[a_{1} + \ln(1 - T_r)\left(
        a_2 + a_3T_r + a_4T_r^2 + a_5T_r^3 \right)\right]

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
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

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Notes
    -----
    This expression is also used for enthalpy of vaporization in [1]_.
    The coefficients from NIST TDE for enthalpy of vaporization are kJ/mol.

    Examples
    --------
    Isooctane at 350 K from [1]_:

    >>> Watson_sigma(T=350.0, Tc=543.836, a1=-3.02417, a2=1.21792, a3=-5.26877e-9, a4=5.62659e-9, a5=-2.27553e-9)
    0.0138340926605649

    References
    ----------
    .. [1] "ThermoData Engine (TDE103b V10.1) User’s Guide." 
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-SurfaceTension/HVPExpansion-SurfaceTension.htm
    '''
    Tr = T/Tc
    l = log(1.0 - Tr)
    return exp(a1 + l*(a2 + Tr*(a3 + Tr*(a4 + a5*Tr))))

def ISTExpansion(T, Tc, a1, a2, a3=0.0, a4=0.0, a5=0.0):
    r'''Calculates air-water surface tension using the IST expansion [1]_
    emperical (parameter-regressed) method developed by NIST.
    
    .. math::
        \sigma = \sum_i a_i\left(1 - \frac{T}{T_c} \right)^i

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
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

    Returns
    -------
    sigma : float
        Liquid surface tension, [N/m]

    Notes
    -----

    Examples
    --------
    Diethyl phthalate at 400 K from [1]_:

    >>> ISTExpansion(T=400.0, Tc=776.0, a1=0.037545, a2=0.0363288)
    0.02672100905515996

    References
    ----------
    .. [1] "ThermoData Engine (TDE103b V10.1) User’s Guide." 
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-SurfaceTension/ISTExpansion-SurfaceTension.htm
    '''
    tau = 1.0 - T/Tc
    return tau*(a1 + tau*(a2 + tau*(a3 + tau*(a4 + a5*tau))))

def Somayajulu(T, Tc, A, B, C):
    r'''Calculates air-water surface tension  using the [1]_
    emperical (parameter-regressed) method. Well regressed, no recent data.

    .. math::
        \sigma=aX^{5/4}+bX^{9/4}+cX^{13/4}

    .. math::
        X=(T_c-T)/T_c

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    A : float
        Regression parameter
    B : float
        Regression parameter
    C : float
        Regression parameter

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    Presently untested, but matches expected values. Internal units are mN/m.
    Form of function returns imaginary results when T > Tc; None is returned
    if this is the case. Function is claimed valid from the triple to the
    critical point. Results can be evaluated beneath the triple point.

    Examples
    --------
    Water at 300 K

    >>> Somayajulu(300, 647.126, 232.713514, -140.18645, -4.890098)
    0.07166386387996758

    References
    ----------
    .. [1] Somayajulu, G. R. "A Generalized Equation for Surface Tension from
       the Triple Point to the Critical Point." International Journal of
       Thermophysics 9, no. 4 (July 1988): 559-66. doi:10.1007/BF00503154.
    '''
    X = (Tc-T)/Tc
    return X*sqrt(sqrt(X))*(A + X*(B + C*X))*1e-3


def Jasper(T, a, b):
    r'''Calculates surface tension of a fluid given two parameters, a linear
    fit in Celcius from [1]_ with data reprinted in [2]_.

    .. math::
        \sigma = a - bT

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    a : float
        Parameter for equation. Chemical specific.
    b : float
        Parameter for equation. Chemical specific.

    Returns
    -------
    sigma : float
        Surface tension [N/m]

    Notes
    -----
    Internal units are mN/m, and degrees Celcius.
    This function has been checked against several references.

    Examples
    --------
    >>> Jasper(298.15, 24, 0.0773)
    0.0220675

    References
    ----------
    .. [1] Jasper, Joseph J. "The Surface Tension of Pure Liquid Compounds."
       Journal of Physical and Chemical Reference Data 1, no. 4
       (October 1, 1972): 841-1010. doi:10.1063/1.3253106.
    .. [2] Speight, James. Lange's Handbook of Chemistry. 16 edition.
       McGraw-Hill Professional, 2005.
    '''
    sigma = (a - b*(T-273.15))*1e-3
    return sigma


### CSP methods


def Brock_Bird(T, Tb, Tc, Pc):
    r'''Calculates air-water surface tension  using the [1]_
    emperical method. Old and tested.

    .. math::
        \sigma = P_c^{2/3}T_c^{1/3}Q(1-T_r)^{11/9}

    .. math::
        Q = 0.1196 \left[ 1 + \frac{T_{br}\ln (P_c/1.01325)}{1-T_{br}}\right]-0.279

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    Numerous arrangements of this equation are available.
    This is DIPPR Procedure 7A: Method for the Surface Tension of Pure,
    Nonpolar, Nonhydrocarbon Liquids
    The exact equation is not in the original paper.
    If the equation yields a negative result, return None.

    Examples
    --------
    p-dichloribenzene at 412.15 K, from DIPPR; value differs due to a slight
    difference in method.

    >>> Brock_Bird(412.15, 447.3, 685, 3.952E6)
    0.02208448325192495

    Chlorobenzene from Poling, as compared with a % error value at 293 K.

    >>> Brock_Bird(293.15, 404.75, 633.0, 4530000.0)
    0.032985686413713036

    References
    ----------
    .. [1] Brock, James R., and R. Byron Bird. "Surface Tension and the
       Principle of Corresponding States." AIChE Journal 1, no. 2
       (June 1, 1955): 174-77. doi:10.1002/aic.690010208
    '''
    Tc_inv = 1.0/Tc
    Tbr = Tb*Tc_inv
    Tr = T*Tc_inv
    Pc = Pc*1e-5  # Convert to bar
    Q = 0.1196*(1.0 + Tbr*log(Pc*(1.0/1.01325))/(1.0 - Tbr)) - 0.279
    sigma = (Pc)**(2.0/3.0)*Tc**(1.0/3.0)*Q*(1.0 - Tr)**(11.0/9.0)
    sigma = sigma*1e-3  # convert to N/m
    return sigma


def Pitzer_sigma(T, Tc, Pc, omega):
    r'''Calculates air-water surface tension using the correlation derived
    by [1]_ from the works of [2]_ and [3]_. Based on critical property CSP
    methods.

    .. math::
        \sigma = P_c^{2/3}T_c^{1/3}\frac{1.86 + 1.18\omega}{19.05}
        \left[ \frac{3.75 + 0.91 \omega}{0.291 - 0.08 \omega}\right]^{2/3} (1-T_r)^{11/9}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    The source of this equation has not been reviewed.
    Internal units of presure are bar, surface tension of mN/m.

    Examples
    --------
    Chlorobenzene from Poling, as compared with a % error value at 293 K.

    >>> Pitzer_sigma(293., 633.0, 4530000.0, 0.249)
    0.03458453513446388

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Curl, R. F., and Kenneth Pitzer. "Volumetric and Thermodynamic
       Properties of Fluids-Enthalpy, Free Energy, and Entropy." Industrial &
       Engineering Chemistry 50, no. 2 (February 1, 1958): 265-74.
       doi:10.1021/ie50578a047
    .. [3] Pitzer, K. S.: Thermodynamics, 3d ed., New York, McGraw-Hill,
       1995, p. 521.
    '''
    Tr = T/Tc
    Pc = Pc*1e-5  # Convert to bar
    sigma = Pc**(2.0/3.0)*Tc**(1.0/3.0)*(1.86 + 1.18*omega)*(1.0/19.05)*(
        (3.75 + 0.91*omega)/(0.291 - 0.08*omega))**(2.0/3.0)*(1.0 - Tr)**(11.0/9.0)
    return sigma*1e-3  # N/m, please


def Sastri_Rao(T, Tb, Tc, Pc, chemicaltype=None):
    r'''Calculates air-water surface tension using the correlation derived by
    [1]_ based on critical property CSP methods and chemical classes.

    .. math::
        \sigma = K P_c^xT_b^y T_c^z\left[\frac{1-T_r}{1-T_{br}}\right]^m

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    The source of this equation has not been reviewed.
    Internal units of presure are bar, surface tension of mN/m.

    Examples
    --------
    Chlorobenzene from Poling, as compared with a % error value at 293 K.

    >>> Sastri_Rao(293.15, 404.75, 633.0, 4530000.0)
    0.03234567739694441

    References
    ----------
    .. [1] Sastri, S. R. S., and K. K. Rao. "A Simple Method to Predict
       Surface Tension of Organic Liquids." The Chemical Engineering Journal
       and the Biochemical Engineering Journal 59, no. 2 (October 1995): 181-86.
       doi:10.1016/0923-0467(94)02946-6.
    '''
    if chemicaltype == 'alcohol':
        k, x, y, z, m = 2.28, 0.25, 0.175, 0, 0.8
    elif chemicaltype == 'acid':
        k, x, y, z, m = 0.125, 0.50, -1.5, 1.85, 11.0/9.0
    else:
        k, x, y, z, m = 0.158, 0.50, -1.5, 1.85, 11.0/9.0
    Tr = T/Tc
    Tbr = Tb/Tc
    Pc = Pc*1E-5  # Convert to bar
    sigma = k*Pc**x*Tb**y*Tc**z*((1.0 - Tr)/(1.0 - Tbr))**m
    sigma = sigma*1e-3  # N/m
    return sigma


def Zuo_Stenby(T, Tc, Pc, omega):
    r'''Calculates air-water surface tension using the reference fluids
    methods of [1]_.

    .. math::
        \sigma^{(1)} = 40.520(1-T_r)^{1.287}

    .. math::
        \sigma^{(2)} = 52.095(1-T_r)^{1.21548}

    .. math::
        \sigma_r = \sigma_r^{(1)}+ \frac{\omega - \omega^{(1)}}
        {\omega^{(2)}-\omega^{(1)}} (\sigma_r^{(2)}-\sigma_r^{(1)})

    .. math::
        \sigma = T_c^{1/3}P_c^{2/3}[\exp{(\sigma_r)} -1]

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    Presently untested. Have not personally checked the sources.
    I strongly believe it is broken.
    The reference values for methane and n-octane are from the DIPPR database.

    Examples
    --------
    Chlorobenzene

    >>> Zuo_Stenby(293., 633.0, 4530000.0, 0.249)
    0.03345569011871088

    References
    ----------
    .. [1] Zuo, You-Xiang, and Erling H. Stenby. "Corresponding-States and
       Parachor Models for the Calculation of Interfacial Tensions." The
       Canadian Journal of Chemical Engineering 75, no. 6 (December 1, 1997):
       1130-37. doi:10.1002/cjce.5450750617
    '''
    Tc_1, Pc_1, omega_1 = 190.56, 4599000.0*1e-5, 0.012
    Tc_2, Pc_2, omega_2 = 568.7, 2490000.0*1e-5, 0.4
    Pc = Pc*1e-5
    Tr = T/Tc

    ST_1 = 40.520*(1.0 - Tr)**1.287  # Methane
    ST_2 = 52.095*(1.0 - Tr)**1.21548  # n-octane

    ST_r_1 = log(1.0 + 0.013537770442486932*ST_1) # Constant from 1/(Tc_1**(1.0/3.0)*Pc_1**(2.0/3.0))
#    ST_r_1 = log(1.0 + ST_1/(Tc_1**(1.0/3.0)*Pc_1**(2.0/3.0)))
    ST_r_2 = log(1.0 + 0.014154874587259097*ST_2) # Constant from /(Tc_2**(1.0/3.0)*Pc_2**(2.0/3.0))
    sigma_r = ST_r_1 + (omega-omega_1)*(ST_r_2-ST_r_1)*2.5773195876288657
#    sigma_r = ST_r_1 + (omega-omega_1)/(omega_2 - omega_1)*(ST_r_2-ST_r_1)
    sigma = Tc**(1.0/3.0)*Pc**(2.0/3.0)*(exp(sigma_r) - 1.0)
    sigma = sigma*1e-3  # N/m, please
    return sigma


def Hakim_Steinberg_Stiel(T, Tc, Pc, omega, StielPolar=0.0):
    r'''Calculates air-water surface tension using the reference fluids methods
    of [1]_.

    .. math::
        \sigma = 4.60104\times 10^{-7} P_c^{2/3}T_c^{1/3}Q_p \left(\frac{1-T_r}{0.4}\right)^m

    .. math::
        Q_p = 0.1574+0.359\omega-1.769\chi-13.69\chi^2-0.51\omega^2+1.298\omega\chi

    .. math::
        m = 1.21+0.5385\omega-14.61\chi-32.07\chi^2-1.65\omega^2+22.03\omega\chi

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]
    StielPolar : float, optional
        Stiel Polar Factor, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    Original equation for m and Q are used. Internal units are atm and mN/m.

    Examples
    --------
    1-butanol, as compared to value in CRC Handbook of 0.02493.

    >>> Hakim_Steinberg_Stiel(298.15, 563.0, 4414000.0, 0.59, StielPolar=-0.07872)
    0.02190790257519

    References
    ----------
    .. [1] Hakim, D. I., David Steinberg, and L. I. Stiel. "Generalized
       Relationship for the Surface Tension of Polar Fluids." Industrial &
       Engineering Chemistry Fundamentals 10, no. 1 (February 1, 1971): 174-75.
       doi:10.1021/i160037a032.
    '''
    omega2 = omega*omega
    StielPolar2 = StielPolar*StielPolar
    Q = (0.1574 + 0.359*omega - 1.769*StielPolar - 13.69*StielPolar2
        - 0.510*omega2 + 1.298*StielPolar*omega)
    m = (1.210 + 0.5385*omega - 14.61*StielPolar - 32.07*StielPolar2
        - 1.656*omega2 + 22.03*StielPolar*omega)
    Tr = T/Tc
    Pc = Pc*(1.0/101325.0)
    sigma = Pc**(2.0/3.)*Tc**(1.0/3.0)*Q*(2.5*(1.0 - Tr))**m
    sigma = sigma*1e-3  # convert to N/m
    return sigma


def Miqueu(T, Tc, Vc, omega):
    r'''Calculates air-water surface tension using the methods of [1]_.

    .. math::
        \sigma = k T_c \left( \frac{N_a}{V_c}\right)^{2/3}
        (4.35 + 4.14 \omega)t^{1.26}(1+0.19t^{0.5} - 0.487t)

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol]
    omega : float
        Acentric factor for fluid, [-]

    Returns
    -------
    sigma : float
        Liquid surface tension, N/m

    Notes
    -----
    Uses Avogadro's constant and the Boltsman constant.
    Internal units of volume are mL/mol and mN/m. However, either a typo
    is in the article or author's work, or my value of k is off by 10; this is
    corrected nonetheless.
    Created with 31 normal fluids, none polar or hydrogen bonded. Has an
    AARD of 3.5%.

    Examples
    --------
    Bromotrifluoromethane, 2.45 mN/m

    >>> Miqueu(300., 340.1, 0.000199, 0.1687)
    0.003474100774091376

    References
    ----------
    .. [1] Miqueu, C, D Broseta, J Satherley, B Mendiboure, J Lachaise, and
       A Graciaa. "An Extended Scaled Equation for the Temperature Dependence
       of the Surface Tension of Pure Compounds Inferred from an Analysis of
       Experimental Data." Fluid Phase Equilibria 172, no. 2 (July 5, 2000):
       169-82. doi:10.1016/S0378-3812(00)00384-8.
    '''
    Vc = Vc*1E6
    t = 1. - T/Tc
    sigma = k*Tc*(N_A/Vc)**(2.0/3.0)*(4.35 + 4.14*omega)*t**1.26*(1.0 + 0.19*sqrt(t) - 0.25*t)*10000.0
    return sigma


def Aleem(T, MW, Tb, rhol, Hvap_Tb, Cpl):
    r'''Calculates vapor-liquid surface tension using the correlation derived by
    [1]_ based on critical property CSP methods.

    .. math::
        \sigma = \phi \frac{MW^{1/3}} {6N_A^{1/3}}\rho_l^{2/3}\left[H_{vap}
        + C_{p,l}(T_b-T)\right]

    .. math::
        \phi = 1 - 0.0047MW + 6.8\times 10^{-6} MW^2

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    MW : float
        Molecular weight [g/mol]
    Tb : float
        Boiling temperature of the fluid [K]
    rhol : float
        Liquid density at T and P [kg/m^3]
    Hvap_Tb : float
        Mass enthalpy of vaporization at the normal boiling point [kg/m^3]
    Cpl : float
        Liquid heat capacity of the chemical at T [J/kg/K]

    Returns
    -------
    sigma : float
        Liquid-vapor surface tension [N/m]

    Notes
    -----
    Internal units of molecuar weight are kg/mol. This model is dimensionally
    consistent.

    This model does not use the critical temperature. After it predicts a
    surface tension of 0 at a sufficiently high temperature, it returns
    negative results. The temperature at which this occurs (the "predicted"
    critical temperature) can be calculated as follows:

    .. math::
        \sigma = 0 \to T_{c,predicted} \text{ at } T_b + \frac{H_{vap}}{Cp_l}

    Because of its dependence on density, it has the potential to model the
    effect of pressure on surface tension.

    Claims AAD of 4.3%. Developed for normal alkanes. Total of 472 data points.
    Behaves worse for higher alkanes. Behaves very poorly overall.

    Examples
    --------
    Methane at 90 K

    >>> Aleem(T=90, MW=16.04246, Tb=111.6, rhol=458.7, Hvap_Tb=510870.,
    ... Cpl=2465.)
    0.01669970230131523

    References
    ----------
    .. [1] Aleem, W., N. Mellon, S. Sufian, M. I. A. Mutalib, and D. Subbarao.
       "A Model for the Estimation of Surface Tension of Pure Hydrocarbon
       Liquids." Petroleum Science and Technology 33, no. 23-24 (December 17,
       2015): 1908-15. doi:10.1080/10916466.2015.1110593.
    '''
    MW = MW*1e-3 # Use kg/mol for consistency with the other units
    sphericity = 1. - MW*(0.0047 - 6.8E-6*MW)
    return sphericity*MW**(1.0/3.0)/(6.*N_A**(1.0/3.0))*rhol**(2.0/3.)*(Hvap_Tb + Cpl*(Tb-T))


def Mersmann_Kind_sigma(T, Tm, Tb, Tc, Pc, n_associated=1):
    r'''Estimates the surface tension of organic liquid substances
    according to the method of [1]_.

    .. math::
        \sigma^* = \frac{\sigma n_{ass}^{1/3}} {(kT_c)^{1/3} T_{rm}P_c^{2/3}}

    .. math::
        \sigma^* = \left(\frac{T_b - T_m}{T_m}\right)^{1/3}
        \left[6.25(1-T_r) + 31.3(1-T_r)^{4/3}\right]

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tm : float
        Melting temperature [K]
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    n_associated : float
        Number of associated molecules in a cluster (2 for alcohols, 1
        otherwise), [-]

    Returns
    -------
    sigma : float
        Liquid-vapor surface tension [N/m]

    Notes
    -----
    In the equation, all quantities must be in SI units. `k` is the boltzman
    constant.

    Examples
    --------
    MTBE at STP (the actual value is 0.0181):

    >>> Mersmann_Kind_sigma(298.15, 164.15, 328.25, 497.1, 3430000.0)
    0.016744311449290426

    References
    ----------
    .. [1] Mersmann, Alfons, and Matthias Kind. "Prediction of Mechanical and
       Thermal Properties of Pure Liquids, of Critical Data, and of Vapor
       Pressure." Industrial & Engineering Chemistry Research, January 31,
       2017. https://doi.org/10.1021/acs.iecr.6b04323.
    '''
    Tr = T/Tc
    sigma_star = ((Tb - Tm)/Tm)**(1.0/3.)*(6.25*(1. - Tr) + 31.3*(1. - Tr)**(4.0/3.))
    sigma = sigma_star*(k*Tc)**(1.0/3.0)*(Tm/Tc)*Pc**(2.0/3.0)*n_associated**(-1.0/3.0)
    return sigma


def API10A32(T, Tc, K_W):
    r'''Calculates the interfacial tension between
    a liquid petroleum fraction and air, using the oil's pseudocritical
    temperature and Watson K Characterization factor.

    .. math::
        \sigma = \frac{673.7\left[\frac{\left(T_c - T\right)}{T_c}\right]^{1.232}}{K_W}

    Parameters
    ----------
    T : float
        Liquid temperature, [K]
    Tc : float
        Pseudocritical temperature (or critical temperature if using
        the equation with a pure component), [K]
    K_W : float
        Watson characterization factor

    Returns
    -------
    sigma : float
        Air-water surface tension, [N/m]

    Notes
    -----
    [1]_ cautions that this should not be applied to coal liquids,
    and that it will give higher errors at pressures above 500 psi.
    [1]_ claims this has an average error of 10.7%.

    This function converges to zero at `Tc`; do not use it above that
    temperature!

    Examples
    --------
    Sample problem in Comments on Procedure 10A3.2.1 of [1]_;

    >>> from fluids.core import F2K, R2K
    >>> API10A32(T=F2K(60), Tc=R2K(1334), K_W=12.4)
    29.577333312096968

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    '''
    return 673.7*((Tc-T)/Tc)**1.232/K_W

### Surface Tension Mixtures

def Winterfeld_Scriven_Davis(xs, sigmas, rhoms):
    r'''Calculates surface tension of a liquid mixture according to
    mixing rules in [1]_ and also in [2]_.

    .. math::
        \sigma_M = \sum_i \sum_j \frac{1}{V_L^{L2}}\left(x_i V_i \right)
        \left( x_jV_j\right)\sqrt{\sigma_i\cdot \sigma_j}

    Parameters
    ----------
    xs : array-like
        Mole fractions of all components, [-]
    sigmas : array-like
        Surface tensions of all components, [N/m]
    rhoms : array-like
        Molar densities of all components, [mol/m^3]

    Returns
    -------
    sigma : float
        Air-liquid surface tension of mixture, [N/m]

    Notes
    -----
    DIPPR Procedure 7C: Method for the Surface Tension of Nonaqueous Liquid
    Mixtures

    Becomes less accurate as liquid-liquid critical solution temperature is
    approached. DIPPR Evaluation:  3-4% AARD, from 107 nonaqueous binary
    systems, 1284 points. Internally, densities are converted to kmol/m^3. The
    Amgat function is used to obtain liquid mixture density in this equation.

    Raises a ZeroDivisionError if either molar volume are zero, and a
    ValueError if a surface tensions of a pure component is negative.

    Examples
    --------
    >>> Winterfeld_Scriven_Davis([0.1606, 0.8394], [0.01547, 0.02877],
    ... [8610., 15530.])
    0.02496738845043982

    References
    ----------
    .. [1] Winterfeld, P. H., L. E. Scriven, and H. T. Davis. "An Approximate
       Theory of Interfacial Tensions of Multicomponent Systems: Applications
       to Binary Liquid-Vapor Tensions." AIChE Journal 24, no. 6
       (November 1, 1978): 1010-14. doi:10.1002/aic.690240610.
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    N = len(xs)
    Vms = [0.0]*N
    rho = 0.0
    for i in range(N):
        Vms[i] = 1e3/rhoms[i]
        rho += xs[i]*Vms[i]
#    rho = 1./rho
    rho = 1.4142135623730951/rho # factor out rt2
    # For speed, transform the Vms array to contain
#    xs[i]*Vms[i]*sigmas_05[i]*rho
    tot = 0.0
    for i in range(N):
        val = sqrt(sigmas[i])*xs[i]*rho*Vms[i]
        Vms[i] = val
        tot += val*val
    tot *= 0.5
    for i in range(N):
        # Symmetric - can be slightly optimized
        temp = 0.0
        for j in range(i):
            temp += Vms[j]
        tot += Vms[i]*temp

    return tot


def Diguilio_Teja(T, xs, sigmas_Tb, Tbs, Tcs):
    r'''Calculates surface tension of a liquid mixture according to
    mixing rules in [1]_.

    .. math::
        \sigma = 1.002855(T^*)^{1.118091} \frac{T}{T_b} \sigma_r

    .. math::
        T^* = \frac{(T_c/T)-1}{(T_c/T_b)-1}

    .. math::
        \sigma_r = \sum x_i \sigma_i

    .. math::
        T_b = \sum x_i T_{b,i}

    .. math::
        T_c = \sum x_i T_{c,i}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    xs : array-like
        Mole fractions of all components
    sigmas_Tb : array-like
        Surface tensions of all components at the boiling point, [N/m]
    Tbs : array-like
        Boiling temperatures of all components, [K]
    Tcs : array-like
        Critical temperatures of all components, [K]

    Returns
    -------
    sigma : float
        Air-liquid surface tension of mixture, [N/m]

    Notes
    -----
    Simple model, however it has 0 citations. Gives similar results to the
    `Winterfeld_Scriven_Davis` model.

    Raises a ValueError if temperature is greater than the mixture's critical
    temperature or if the given temperature is negative, or if the mixture's
    boiling temperature is higher than its critical temperature.

    [1]_ claims a 4.63 percent average absolute error on 21 binary and 4
    ternary non-aqueous systems. [1]_ also considered Van der Waals mixing
    rules for `Tc`, but found it provided a higher error of 5.58%

    Examples
    --------
    >>> Diguilio_Teja(T=298.15, xs=[0.1606, 0.8394],
    ... sigmas_Tb=[0.01424, 0.02530], Tbs=[309.21, 312.95], Tcs=[469.7, 508.0])
    0.025716823875045505

    References
    ----------
    .. [1] Diguilio, Ralph, and Amyn S. Teja. "Correlation and Prediction of
       the Surface Tensions of Mixtures." The Chemical Engineering Journal 38,
       no. 3 (July 1988): 205-8. doi:10.1016/0300-9467(88)80079-0.
    '''
    Tc, Tb, sigmar = 0.0, 0.0, 0.0
    for i in range(len(xs)):
        Tc += Tcs[i]*xs[i]
        Tb += Tbs[i]*xs[i]
        sigmar += sigmas_Tb[i]*xs[i]
    if T > Tc:
        raise ValueError('T > Tc according to Kays rule - model is not valid in this range.')
    Tst = (Tc/T - 1.)/(Tc/Tb - 1.0)
    return 1.002855*Tst**1.118091*(T/Tb)*sigmar


def Weinaug_Katz(parachors, Vml, Vmg, xs, ys):
    r'''Calculates surface tension of a liquid mixture according to
    mixing rules in [1]_ and also in [2]_. This is based on the
    Parachor concept. This is called the Macleod-Sugden model in some places.

    .. math::
        \sigma_M = \left[\sum_i P_i\left( \frac{x_i}{V_{m,l}}
         - \frac{y_i}{V_{m,g}}\right) \right]^4

    Parameters
    ----------
    parachors : list[float]
        Parachors of each component, [N^0.25*m^2.75/mol]
    Vml : float
        Liquid mixture molar volume, [m^3/mol]
    Vmg : float
        Gas mixture molar volume; this can be set to zero at
        low pressures, [m^3/mol]
    xs : list[float]
        Mole fractions of all components in liquid phase, [-]
    xs : list[float]
        Mole fractions of all components in gas phase, [-]

    Returns
    -------
    sigma : float
        Air-liquid surface tension of mixture, [N/m]

    Notes
    -----
    This expression is efficient and does not require pure component
    surface tensions. Its accuracy is dubious.

    Examples
    --------
    >>> Weinaug_Katz([5.1e-5, 7.2e-5], Vml=0.000125, Vmg=0.02011, xs=[.4, .6], ys=[.6, .4])
    0.06547479150776776

    Neglect the vapor phase density by setting `Vmg` to a high value:

    >>> Weinaug_Katz([5.1e-5, 7.2e-5], Vml=0.000125, Vmg=1e100, xs=[.4, .6], ys=[.6, .4])
    0.06701752894095361

    References
    ----------
    .. [1] Weinaug, Charles F., and Donald L. Katz. "Surface Tensions of
       Methane-Propane Mixtures." Industrial & Engineering Chemistry 35,
       no. 2 (February 1, 1943): 239-246. https://doi.org/10.1021/ie50398a028.
    .. [2] Pedersen, Karen Schou, Aage Fredenslund, and Per Thomassen.
       Properties of Oils and Natural Gases. Vol. 5. Gulf Pub Co, 1989.
    '''
    tot = 0.0
    rhoml = 1.0/Vml
    rhomg = 1.0/Vmg
    for i in range(len(parachors)):
        tot += parachors[i]*(xs[i]*rhoml - ys[i]*rhomg)
    tot *= tot
    tot *= tot # fourth power it
    return tot

### Water-hydrocarbon interfacial tensions


def Meybodi_Daryasafar_Karimi(rho_water, rho_oil, T, Tc):
    r'''Calculates the interfacial tension between water and a hydrocabon
    liquid according to the correlation of [1]_.

    .. math::
        \gamma_{hw} = \left(\frac{A_1 + A_2 \Delta \rho + A_3\Delta\rho^2
        + A_4\Delta\rho^3} {A_5 + A_6\frac{T^{A_7}}{T_{c,h}} + A_8T^{A_9}}
        \right)^{A_{10}}

    Parameters
    ----------
    rho_water : float
        The density of the aqueous phase, [kg/m^3]
    rho_oil : float
        The density of the hydrocarbon phase, [kg/m^3]
    T : float
        Temperature of the fluid, [K]
    Tc : float
        Critical temperature of the hydrocarbon mixture, [K]

    Returns
    -------
    sigma : float
        Hydrocarbon-water surface tension [N/m]

    Notes
    -----
    Internal units of the equation are g/mL and mN/m.

    Examples
    --------
    >>> Meybodi_Daryasafar_Karimi(980, 760, 580, 914)
    0.02893598143089256

    References
    ----------
    .. [1] Kalantari Meybodi, Mahdi, Amin Daryasafar, and Masoud Karimi.
       "Determination of Hydrocarbon-Water Interfacial Tension Using a New
       Empirical Correlation."  Fluid Phase Equilibria 415 (May 15, 2016):
       42-50. doi:10.1016/j.fluid.2016.01.037.
    '''
    A1 = -1.3687340042E-1
    A2 = -3.0391828884E-1
    A3 = 5.6225871072E-1
    A4 = -3.3074367079E-1
    A5 = -3.0050179309E0
    A6 = 5.8914210205E-5
    A7 = -4.1388901263E0
    A8 = 3.0084299030E0
    A9 = -3.8203072876E-3
#    A10 = 3.5000000000E0
    drho = abs(rho_water - rho_oil)*1e-3 # Correlation in units of g/mL
    sigma = ((A1 + drho*(A2 + drho*(A3 + A4*drho)))
             /(A5 + A6*T**A7/Tc + A8*T**A9))
    return sigma*sigma*sigma*sqrt(sigma)*1e-3 # mN/m to N/m
