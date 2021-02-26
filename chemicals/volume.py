# -*- coding: utf-8 -*-
r"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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

This module contains various volume/density estimation routines, dataframes
of fit coefficients, and mixing rules.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Pure Low Pressure Liquid Correlations
-------------------------------------
.. autofunction:: chemicals.volume.volume_VDI_PPDS
.. autofunction:: chemicals.volume.Rackett
.. autofunction:: chemicals.volume.COSTALD
.. autofunction:: chemicals.volume.Yen_Woods_saturation
.. autofunction:: chemicals.volume.Yamada_Gunn
.. autofunction:: chemicals.volume.Townsend_Hales
.. autofunction:: chemicals.volume.Bhirud_normal
.. autofunction:: chemicals.volume.Campbell_Thodos
.. autofunction:: chemicals.volume.SNM0

Pure High Pressure Liquid Correlations
--------------------------------------
.. autofunction:: chemicals.volume.COSTALD_compressed

Liquid Mixing Rules
-------------------
.. autofunction:: chemicals.volume.Amgat
.. autofunction:: chemicals.volume.Rackett_mixture
.. autofunction:: chemicals.volume.COSTALD_mixture

Gas Correlations
----------------
Gas volumes are predicted with one of:

1) An equation of state
2) A virial coefficient model
3) The ideal gas law

Equations of state do much more than predict volume however. An implementation
of many of them can be found in `thermo <https://github.com/CalebBell/thermo>`_.

Virial functions are implemented in :obj:`chemicals.virial`.

.. autofunction:: chemicals.volume.ideal_gas

Pure Solid Correlations
-----------------------
Solid density does not depend on pressure significantly, and unless operating
in the geochemical or astronomical domain is normally neglected.

.. autofunction:: chemicals.volume.Goodman

Pure Component Solid Fit Correlations
-------------------------------------
.. autofunction:: chemicals.volume.CRC_inorganic

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: rho_data_COSTALD

    Coefficients for the :obj:`COSTALD` method from [3]_; 192 fluids have
    coefficients published.

.. data:: rho_data_SNM0

    Coefficients for the :obj:`SNM0` method for 73 fluids from [2]_.

.. data:: rho_data_Perry_8E_105_l

    Coefficients for :obj:`chemicals.dippr.EQ105` from [1]_ for 344 fluds. Note
    this is in terms of molar density; to obtain molar volume, invert the result!

.. data:: rho_data_VDI_PPDS_2

    Coefficients in [5]_ developed by the PPDS using :obj:`chemicals.dippr.EQ116`
    but in terms of mass density [kg/m^3]; Valid up to the critical temperature,
    and extrapolates to very low temperatures well.

.. data:: rho_data_CRC_inorg_l

    Single-temperature coefficient linear model in terms of mass density
    for the density of inorganic liquids. Data is available for 177 fluids
    normally valid over a
    narrow range above the melting point, from [4]_; described in
    :obj:`CRC_inorganic`.

.. data:: rho_data_CRC_inorg_l_const

    Constant inorganic liquid molar volumes published in [4]_.

.. data:: rho_data_CRC_inorg_s_const

    Constant solid densities molar volumes published in [4]_.

.. data:: rho_data_CRC_virial

    Coefficients for a tempereture polynomial (T in Kelvin) for the second `B`
    virial coefficient published in [4]_. The form of the equation is
    :math:`B = (a_1 + t(a_2 + t(a_3 + t(a_4 + a_5 t)))) \times 10^{-6}` with
    :math:`t = \frac{298.15}{T} - 1` and then `B` will be in units of m^3/mol.


.. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    8E. McGraw-Hill Professional, 2007.
.. [2] Mchaweh, A., A. Alsaygh, Kh. Nasrifar, and M. Moshfeghian.
    "A Simplified Method for Calculating Saturated Liquid Densities."
    Fluid Phase Equilibria 224, no. 2 (October 1, 2004): 157-67.
    doi:10.1016/j.fluid.2004.06.054
.. [3] Hankinson, Risdon W., and George H. Thomson. "A New Correlation for
    Saturated Densities of Liquids and Their Mixtures." AIChE Journal
    25, no. 4 (1979): 653-663. doi:10.1002/aic.690250412
.. [4] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
    Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
.. [5] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
    Berlin; New York:: Springer, 2010.

The structure of each dataframe is shown below:

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.volume.rho_data_COSTALD

    In [3]: chemicals.volume.rho_data_SNM0

    In [4]: chemicals.volume.rho_data_Perry_8E_105_l

    In [5]: chemicals.volume.rho_data_VDI_PPDS_2

    In [6]: chemicals.volume.rho_data_CRC_inorg_l

    In [7]: chemicals.volume.rho_data_CRC_inorg_l_const

    In [8]: chemicals.volume.rho_data_CRC_inorg_s_const

    In [9]: chemicals.volume.rho_data_CRC_virial

"""

from __future__ import division

__all__ = ['volume_VDI_PPDS', 'Yen_Woods_saturation', 'Rackett', 'Yamada_Gunn', 'Townsend_Hales',
'Bhirud_normal', 'COSTALD', 'Campbell_Thodos', 'SNM0', 'CRC_inorganic',
'COSTALD_compressed', 'Amgat', 'Rackett_mixture', 'COSTALD_mixture',
'ideal_gas', 'Goodman',
]

import os
from fluids.numerics import np, splev, implementation_optimize_tck
from fluids.constants import R, atm_inv
from chemicals.utils import log, exp, isnan
from chemicals.utils import Vm_to_rho, mixing_simple
from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import data_source, register_df_source

folder = os_path_join(source_path, 'Density')

register_df_source(folder, 'COSTALD Parameters.tsv')
register_df_source(folder, 'Mchaweh SN0 deltas.tsv')
register_df_source(folder, 'Perry Parameters 105.tsv')
register_df_source(folder, 'CRC Liquid Inorganic Constant Densities.tsv')
register_df_source(folder, 'CRC Solid Inorganic Constant Densities.tsv')

register_df_source(folder, 'VDI PPDS Density of Saturated Liquids.tsv', csv_kwargs={
        'dtype':{'rhoc': float}})
register_df_source(folder, 'CRC Inorganics densties of molten compounds and salts.tsv', csv_kwargs={
        'dtype':{'rho': float}})
register_df_source(folder, 'CRC Virial polynomials.tsv', csv_kwargs={
        'dtype':{'a1': float, 'a2': float, 'a3': float, 'a4': float, 'a5': float}})


_rho_data_loaded = False
def _load_rho_data():
    global _rho_data_loaded, rho_data_COSTALD, rho_data_SNM0
    global rho_data_Perry_8E_105_l, rho_values_Perry_8E_105_l
    global rho_data_VDI_PPDS_2, rho_values_VDI_PPDS_2
    global rho_data_CRC_inorg_l, rho_values_CRC_inorg_l
    global rho_data_CRC_inorg_l_const, rho_data_CRC_inorg_s_const
    global rho_data_CRC_virial, rho_values_CRC_virial

    rho_data_COSTALD = data_source('COSTALD Parameters.tsv')
    rho_data_SNM0 = data_source('Mchaweh SN0 deltas.tsv')
    rho_data_Perry_8E_105_l = data_source('Perry Parameters 105.tsv')
    rho_values_Perry_8E_105_l = np.array(rho_data_Perry_8E_105_l.values[:, 1:], dtype=float)

    rho_data_VDI_PPDS_2 = data_source('VDI PPDS Density of Saturated Liquids.tsv')
    rho_values_VDI_PPDS_2 = np.array(rho_data_VDI_PPDS_2.values[:, 1:], dtype=float)

    rho_data_CRC_inorg_l = data_source('CRC Inorganics densties of molten compounds and salts.tsv')
    rho_values_CRC_inorg_l = np.array(rho_data_CRC_inorg_l.values[:, 1:], dtype=float)

    rho_data_CRC_inorg_l_const = data_source('CRC Liquid Inorganic Constant Densities.tsv')
    rho_data_CRC_inorg_s_const = data_source('CRC Solid Inorganic Constant Densities.tsv')

    rho_data_CRC_virial = data_source('CRC Virial polynomials.tsv')
    rho_values_CRC_virial = np.array(rho_data_CRC_virial.values[:, 1:], dtype=float)

if PY37:
    def __getattr__(name):
        if name in ('rho_data_COSTALD', 'rho_data_SNM0', 'rho_data_Perry_8E_105_l',
                    'rho_values_Perry_8E_105_l', 'rho_data_VDI_PPDS_2',
                    'rho_values_VDI_PPDS_2', 'rho_data_CRC_inorg_l',
                    'rho_values_CRC_inorg_l', 'rho_data_CRC_inorg_l_const',
                    'rho_data_CRC_inorg_s_const', 'rho_data_CRC_virial',
                    'rho_values_CRC_virial'):
            _load_rho_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_rho_data()

def volume_VDI_PPDS(T, Tc, rhoc, a, b, c, d, MW=None):
    r'''Calculates saturation liquid volume, using the critical properties
    and fitted coefficients from [1]_.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    rhoc : float
        Critical density of fluid [kg/m^3]
    a,b,c,d : float
        Fitted coefficients [-]
    MW : float, optional
        Molecular weight of chemical [g/mol]

    Returns
    -------
    Vs : float
        Saturation liquid molar volume or density, [m^3/mol if MW given; kg/m^3 otherwise]

    Examples
    --------
    Calculate density of nitrogen in kg/m3 at 300 K:

    >>> volume_VDI_PPDS(300, 126.19, 313, 470.922, 493.251, -560.469, 389.611)
    313.0

    Calculate molar volume of nitrogen in m3/mol at 300 K:

    >>> volume_VDI_PPDS(300, 126.19, 313, 470.922, 493.251, -560.469, 389.611, 28.01)
    8.9488817891e-05

    References
    ----------
    .. [1] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
           Berlin; New York:: Springer, 2010.
    '''
    tau = 1. - T/Tc if T < Tc else 0.
    rho = rhoc + a*tau**0.35 + b*tau**(2/3.) + c*tau + d*tau**(4/3.)
    return rho if MW is None else 0.001 * MW / rho

### Critical-properties based

def Yen_Woods_saturation(T, Tc, Vc, Zc):
    r'''Calculates saturation liquid volume, using the Yen and Woods [1]_ CSP
    method and a chemical's critical properties.

    The molar volume of a liquid is given by:

    .. math::
        Vc/Vs = 1 + A(1-T_r)^{1/3} + B(1-T_r)^{2/3} + D(1-T_r)^{4/3}

    .. math::
        D = 0.93-B

    .. math::
        A = 17.4425 - 214.578Z_c + 989.625Z_c^2 - 1522.06Z_c^3

    .. math::
        B = -3.28257 + 13.6377Z_c + 107.4844Z_c^2-384.211Z_c^3
        \text{ if } Zc \le 0.26

    .. math::
        B = 60.2091 - 402.063Z_c + 501.0 Z_c^2 + 641.0 Z_c^3
        \text{ if } Zc \ge 0.26


    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol]
    Zc : float
        Critical compressibility of fluid, [-]

    Returns
    -------
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Notes
    -----
    Original equation was in terms of density, but it is converted here.

    No example has been found, nor are there points in the article. However,
    it is believed correct. For compressed liquids with the Yen-Woods method,
    see the `YenWoods_compressed` function.

    Examples
    --------
    >>> Yen_Woods_saturation(300, 647.14, 55.45E-6, 0.245)
    1.769533076529574e-05

    References
    ----------
    .. [1] Yen, Lewis C., and S. S. Woods. "A Generalized Equation for Computer
       Calculation of Liquid Densities." AIChE Journal 12, no. 1 (1966):
       95-99. doi:10.1002/aic.690120119
    '''
    Tr = T/Tc
    A = Zc*(Zc*(989.625 - 1522.06*Zc) - 214.578) + 17.4425
    if Zc <= 0.26:
        B = Zc*(Zc*(107.4844 - 384.211*Zc) + 13.6377) - 3.28257
    else:
        B = Zc*(Zc*(641.0*Zc + 501.0) - 402.063) + 60.2091
    D = 0.93 - B
    tau_cbrt = (1.0 - Tr)**(1/3.)
    Vm = Vc/(tau_cbrt*(A + tau_cbrt*(B + D*tau_cbrt*tau_cbrt)) + 1.0)
    return Vm


def Rackett(T, Tc, Pc, Zc):
    r'''Calculates saturation liquid volume, using Rackett CSP method and
    critical properties.

    The molar volume of a liquid is given by:

    .. math::
        V_s = \frac{RT_c}{P_c}{Z_c}^{[1+(1-{T/T_c})^{2/7} ]}

    Units are all currently in m^3/mol - this can be changed to kg/m^3

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    Zc : float
        Critical compressibility of fluid, [-]

    Returns
    -------
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Notes
    -----
    According to Reid et. al, underpredicts volume for compounds with Zc < 0.22

    Examples
    --------
    Propane, example from the API Handbook

    >>> Vm_to_rho(Rackett(272.03889, 369.83, 4248000.0, 0.2763), 44.09562)
    531.3221411755724

    References
    ----------
    .. [1] Rackett, Harold G. "Equation of State for Saturated Liquids."
       Journal of Chemical & Engineering Data 15, no. 4 (1970): 514-517.
       doi:10.1021/je60047a012
    '''
    return R*Tc/Pc*Zc**(1.0 + (1.0 - T/Tc)**(2.0/7.))


def Yamada_Gunn(T, Tc, Pc, omega):
    r'''Calculates saturation liquid volume, using Yamada and Gunn CSP method
    and a chemical's critical properties and acentric factor.

    The molar volume of a liquid is given by:

    .. math::
        V_s = \frac{RT_c}{P_c}{(0.29056-0.08775\omega)}^{[1+(1-{T/T_c})^{2/7}]}

    Units are in m^3/mol.

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
    Vs : float
        saturation liquid volume, [m^3/mol]

    Notes
    -----
    This equation is an improvement on the Rackett equation.
    This is often presented as the Rackett equation.
    The acentric factor is used here, instead of the critical compressibility
    A variant using a reference fluid also exists

    Examples
    --------
    >>> Yamada_Gunn(300, 647.14, 22048320.0, 0.245)
    2.188284384699659e-05

    References
    ----------
    .. [1] Gunn, R. D., and Tomoyoshi Yamada. "A Corresponding States
        Correlation of Saturated Liquid Volumes." AIChE Journal 17, no. 6
        (1971): 1341-45. doi:10.1002/aic.690170613
    .. [2] Yamada, Tomoyoshi, and Robert D. Gunn. "Saturated Liquid Molar
        Volumes. Rackett Equation." Journal of Chemical & Engineering Data 18,
        no. 2 (1973): 234-36. doi:10.1021/je60057a006
    '''
    return R*Tc/Pc*(0.29056 - 0.08775*omega)**(1.0 + (1.0 - T/Tc)**(2.0/7.))


def Townsend_Hales(T, Tc, Vc, omega):
    r'''Calculates saturation liquid density, using the Townsend and Hales
    CSP method as modified from the original Riedel equation. Uses
    chemical critical volume and temperature, as well as acentric factor

    The density of a liquid is given by:

    .. math::
        Vs = V_c/\left(1+0.85(1-T_r)+(1.692+0.986\omega)(1-T_r)^{1/3}\right)

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
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Notes
    -----
    The requirement for critical volume and acentric factor requires all data.

    Examples
    --------
    >>> Townsend_Hales(300, 647.14, 55.95E-6, 0.3449)
    1.8007361992619923e-05

    References
    ----------
    .. [1] Hales, J. L, and R Townsend. "Liquid Densities from 293 to 490 K of
       Nine Aromatic Hydrocarbons." The Journal of Chemical Thermodynamics
       4, no. 5 (1972): 763-72. doi:10.1016/0021-9614(72)90050-X
    '''
    Tr = T/Tc
    return Vc/(1.0 + 0.85*(1.0-Tr) + (1.692 + 0.986*omega)*(1.0-Tr)**(1.0/3.))

Bhirud_normal_Trs = [0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994,
            0.996, 0.998, 0.999, 1]
Bhirud_normal_lnU0s = [-1.6198, -1.604, -1.59, -1.578, -1.564, -1.548, -1.533,
              -1.515, -1.489, -1.454, -1.425, -1.243]
Bhirud_normal_lnU1 = [-0.4626, -0.459, -0.451, -0.441, -0.428, -0.412, -0.392,
              -0.367, -0.337, -0.302, -0.283, -0.2629]

Bhirud_normal_lnU0_tck = implementation_optimize_tck([[0.98, 0.98, 0.98, 0.98, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.996, 0.998, 1.0, 1.0, 1.0, 1.0],
                                                     [-1.6198000000000001, -1.6090032590995371, -1.5930934818009272, -1.5785097772986094, -1.5645133351302174, -1.547436882180519, -1.5337391361477044, -1.5156065732286643, -1.4938345709376377, -1.4430551430207872, -1.4529816189930713, -1.243, 0.0, 0.0, 0.0, 0.0],
                                                     3])
Bhirud_normal_lnU1_tck = implementation_optimize_tck([[0.98, 0.98, 0.98, 0.98, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.996, 0.998, 1.0, 1.0, 1.0, 1.0],
                                                      [-0.4626000000000001, -0.4624223420145324, -0.4543553159709354, -0.4416003593769303, -0.4284319856698448, -0.41267169794369013, -0.39288122255539465, -0.3678034118347314, -0.3379051301056794, -0.30257606774255136, -0.2767190885302606, -0.2629, 0.0, 0.0, 0.0, 0.0],
                                                     3])


def Bhirud_normal(T, Tc, Pc, omega):
    r'''Calculates saturation liquid density using the Bhirud [1]_ CSP method.
    Uses Critical temperature and pressure and acentric factor.

    The density of a liquid is given by:

    .. math::
        \ln \frac{P_c}{\rho RT} = \ln U^{(0)} + \omega\ln U^{(1)}

    .. math::
        \ln U^{(0)} = 1.396 44 - 24.076T_r+ 102.615T_r^2
        -255.719T_r^3+355.805T_r^4-256.671T_r^5 + 75.1088T_r^6

    .. math::
        \ln U^{(1)} = 13.4412 - 135.7437 T_r + 533.380T_r^2-
        1091.453T_r^3+1231.43T_r^4 - 728.227T_r^5 + 176.737T_r^6

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
    Vm : float
        Saturated liquid molar volume, [mol/m^3]

    Notes
    -----
    Claimed inadequate by others.

    An interpolation table for ln U values are used from Tr = 0.98 - 1.000.
    Has terrible behavior at low reduced temperatures.

    Examples
    --------
    Pentane

    >>> Bhirud_normal(280.0, 469.7, 33.7E5, 0.252)
    0.00011249657842514176

    References
    ----------
    .. [1] Bhirud, Vasant L. "Saturated Liquid Densities of Normal Fluids."
       AIChE Journal 24, no. 6 (November 1, 1978): 1127-31.
       doi:10.1002/aic.690240630
    '''
    Tr = T/Tc
    if Tr <= 0.98:
        lnU0 = Tr*(Tr*(Tr*(Tr*(Tr*(75.1088*Tr - 256.671) + 355.805) - 255.719)
               + 102.615) - 24.076) + 1.39644
        lnU1 = Tr*(Tr*(Tr*(Tr*(Tr*(176.737*Tr - 728.227) + 1231.43) - 1091.453)
               + 533.38) - 135.7437) + 13.4412
    elif Tr > 1.0:
        raise ValueError('Critical phase, correlation does not apply')
    else:
        lnU0 = float(splev(Tr, Bhirud_normal_lnU0_tck))
        lnU1 = float(splev(Tr, Bhirud_normal_lnU1_tck))

    Unonpolar = exp(lnU0 + omega*lnU1)
    Vm = Unonpolar*R*T/Pc
    return Vm


def COSTALD(T, Tc, Vc, omega):
    r'''Calculate saturation liquid density using the COSTALD CSP method.

    A popular and accurate estimation method. If possible, fit parameters are
    used; alternatively critical properties work well.

    The density of a liquid is given by:

    .. math::
        V_s=V^*V^{(0)}[1-\omega_{SRK}V^{(\delta)}]

    .. math::
        V^{(0)}=1-1.52816(1-T_r)^{1/3}+1.43907(1-T_r)^{2/3}
        - 0.81446(1-T_r)+0.190454(1-T_r)^{4/3}

    .. math::
        V^{(\delta)}=\frac{-0.296123+0.386914T_r-0.0427258T_r^2-0.0480645T_r^3}
        {T_r-1.00001}

    Units are that of critical or fit constant volume.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol].
        This parameter is alternatively a fit parameter
    omega : float
        (ideally SRK) Acentric factor for fluid, [-]
        This parameter is alternatively a fit parameter.

    Returns
    -------
    Vs : float
        Saturation liquid volume

    Notes
    -----
    196 constants are fit to this function in [1]_.
    Range: 0.25 < Tr < 0.95, often said to be to 1.0

    This function has been checked with the API handbook example problem.

    Examples
    --------
    Propane, from an example in the API Handbook:

    >>> Vm_to_rho(COSTALD(272.03889, 369.83333, 0.20008161E-3, 0.1532), 44.097)
    530.3009967969844

    References
    ----------
    .. [1] Hankinson, Risdon W., and George H. Thomson. "A New Correlation for
       Saturated Densities of Liquids and Their Mixtures." AIChE Journal
       25, no. 4 (1979): 653-663. doi:10.1002/aic.690250412
    '''
    if T > Tc:
        T = Tc
    Tr = T/Tc
    tau = 1.0 - Tr
    tau_cbrt = (tau)**(1.0/3.)
    V_delta = (-0.296123 + Tr*(Tr*(-0.0480645*Tr - 0.0427258) + 0.386914))/(Tr - 1.00001)
    V_0 = tau_cbrt*(tau_cbrt*(tau_cbrt*(0.190454*tau_cbrt - 0.81446) + 1.43907) - 1.52816) + 1.0
    return Vc*V_0*(1.0 - omega*V_delta)


def Campbell_Thodos(T, Tb, Tc, Pc, MW, dipole=0.0, has_hydroxyl=False):
    r'''Calculate saturation liquid density using the Campbell-Thodos [1]_
    CSP method.

    An old and uncommon estimation method.

    .. math::
        V_s = \frac{RT_c}{P_c}{Z_{RA}}^{[1+(1-T_r)^{2/7}]}

    .. math::
        Z_{RA} = \alpha + \beta(1-T_r)

    .. math::
        \alpha = 0.3883-0.0179s

    .. math::
        s = T_{br} \frac{\ln P_c}{(1-T_{br})}

    .. math::
        \beta = 0.00318s-0.0211+0.625\Lambda^{1.35}

    .. math::
        \Lambda = \frac{P_c^{1/3}} { MW^{1/2} T_c^{5/6}}

    For polar compounds:

    .. math::
        \theta = P_c \mu^2/T_c^2

    .. math::
        \alpha = 0.3883 - 0.0179s - 130540\theta^{2.41}

    .. math::
        \beta = 0.00318s - 0.0211 + 0.625\Lambda^{1.35} + 9.74\times
        10^6 \theta^{3.38}

    Polar Combounds with hydroxyl groups (water, alcohols)

    .. math::
        \alpha = \left[0.690T_{br} -0.3342 + \frac{5.79\times 10^{-10}}
        {T_{br}^{32.75}}\right] P_c^{0.145}

    .. math::
        \beta = 0.00318s - 0.0211 + 0.625 \Lambda^{1.35} + 5.90\Theta^{0.835}

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
    MW : float
        Molecular weight of the fluid [g/mol]
    dipole : float, optional
        Dipole moment of the fluid [debye]
    has_hydroxyl : bool, optional
        Swith to use the hydroxyl variant for polar fluids

    Returns
    -------
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Notes
    -----
    If a dipole is provided, the polar chemical method is used.
    The paper is an excellent read.
    Pc is internally converted to atm.

    Examples
    --------
    Ammonia, from [1]_.

    >>> Campbell_Thodos(T=405.45, Tb=239.82, Tc=405.45, Pc=111.7*101325, MW=17.03, dipole=1.47)
    7.347366126245e-05

    References
    ----------
    .. [1] Campbell, Scott W., and George Thodos. "Prediction of Saturated
       Liquid Densities and Critical Volumes for Polar and Nonpolar
       Substances." Journal of Chemical & Engineering Data 30, no. 1
       (January 1, 1985): 102-11. doi:10.1021/je00039a032.
    '''
    Tc_inv = 1.0/Tc
    Tr = T * Tc_inv
    Tbr = Tb * Tc_inv
    Pc = Pc/101325.
    s = Tbr*log(Pc)/(1.0 - Tbr)
    Lambda = Pc**(1.0/3.)/(MW**0.5*Tc**(5/6.))
    beta = 0.00318*s - 0.0211 + 0.625*Lambda**(1.35)
    if dipole is None:
        alpha = 0.3883 - 0.0179*s
    else:
        theta = Pc*dipole*dipole/(Tc*Tc)
        beta += 9.74E6 * theta**3.38
        if has_hydroxyl:
            beta += 5.90*theta**0.835
            alpha = (0.69*Tbr - 0.3342 + 5.79E-10*Tbr**-32.75)*Pc**0.145
        else:
            alpha = 0.3883 - 0.0179*s - 130540 * theta**2.41
    Zra = alpha + beta*(1.0 - Tr)
    p = 1.0 if T == Tc else (1.0 + (1.0 - Tr)**(2.0/7.))
    Vs = R*Tc/(Pc*101325.0)*Zra**p
    return Vs


def SNM0(T, Tc, Vc, omega, delta_SRK=None):
    r'''Calculates saturated liquid density using the Mchaweh, Moshfeghian
    model [1]_. Designed for simple calculations.

    .. math::
        V_s = V_c/(1+1.169\tau^{1/3}+1.818\tau^{2/3}-2.658\tau+2.161\tau^{4/3}

    .. math::
        \tau = 1-\frac{(T/T_c)}{\alpha_{SRK}}

    .. math::
        \alpha_{SRK} = [1 + m(1-\sqrt{T/T_C}]^2

    .. math::
        m = 0.480+1.574\omega-0.176\omega^2

    If the fit parameter `delta_SRK` is provided, the following is used:

    .. math::
        V_s = V_C/(1+1.169\tau^{1/3}+1.818\tau^{2/3}-2.658\tau+2.161\tau^{4/3})
        /\left[1+\delta_{SRK}(\alpha_{SRK}-1)^{1/3}\right]

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
    delta_SRK : float, optional
        Fitting parameter [-]

    Returns
    -------
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Notes
    -----
    73 fit parameters have been gathered from the article.

    Examples
    --------
    Argon, without the fit parameter and with it. Tabulated result in Perry's
    is 3.4613e-05. The fit increases the error on this occasion.

    >>> SNM0(121, 150.8, 7.49e-05, -0.004)
    3.4402256402733416e-05
    >>> SNM0(121, 150.8, 7.49e-05, -0.004, -0.03259620)
    3.493288100008123e-05

    References
    ----------
    .. [1] Mchaweh, A., A. Alsaygh, Kh. Nasrifar, and M. Moshfeghian.
       "A Simplified Method for Calculating Saturated Liquid Densities."
       Fluid Phase Equilibria 224, no. 2 (October 1, 2004): 157-67.
       doi:10.1016/j.fluid.2004.06.054
    '''
    Tr = T/Tc
    m = 0.480 + 1.574*omega - 0.176*omega*omega
    alpha_SRK = (1.0 + m*(1.0 - Tr**0.5))**2
    tau = 1. - Tr/alpha_SRK

    rho0 = 1. + 1.169*tau**(1/3.) + 1.818*tau**(2/3.) - 2.658*tau + 2.161*tau**(4/3.)
    V0 = 1./rho0

    if delta_SRK is None:
        return Vc*V0
    else:
        return Vc*V0/(1. + delta_SRK*(alpha_SRK - 1.0)**(1.0/3.0))


def CRC_inorganic(T, rho0, k, Tm, MW=None):
    r'''Calculates liquid density of a molten element or salt at temperature
    above the melting point. Some coefficients are given nearly up to the
    boiling point.

    The mass density of the inorganic liquid is given by:

    .. math::
        \rho = \rho_{0} - k(T-T_m)

    Parameters
    ----------
    T : float
        Temperature of the liquid, [K]
    rho0 : float
        Mass density of the liquid at Tm, [kg/m^3]
    k : float
        Linear temperature dependence of the mass density, [kg/m^3/K]
    Tm : float
        The normal melting point, used in the correlation [K]
    MW : float, optional
        Molecular weight of chemical [g/mol]

    Returns
    -------
    rho : float
        Mass density of molten metal or salt, [m^3/mol if MW given; kg/m^3 otherwise]

    Notes
    -----
    [1]_ has units of g/mL. While the individual densities could have been
    converted to molar units, the temperature coefficient could only be
    converted by refitting to calculated data. To maintain compatibility with
    the form of the equations, this was not performed.

    This linear form is useful only in small temperature ranges.
    Coefficients for one compound could be used to predict the temperature
    dependence of density of a similar compound.

    Examples
    --------
    >>> CRC_inorganic(300, 2370.0, 2.687, 239.08)
    2206.30796

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
        Chemistry and Physics, 95E. [Boca Raton, FL]: CRC press, 2014.
    '''
    rho = rho0 - k*(T-Tm)
    return rho if MW is None else 0.001 * MW / rho


def COSTALD_compressed(T, P, Psat, Tc, Pc, omega, Vs):
    r'''Calculates compressed-liquid volume, using the COSTALD [1]_ CSP
    method and a chemical's critical properties.

    The molar volume of a liquid is given by:

    .. math::
        V = V_s\left( 1 - C \ln \frac{B + P}{B + P^{sat}}\right)

    .. math::
        \frac{B}{P_c} = -1 + a\tau^{1/3} + b\tau^{2/3} + d\tau + e\tau^{4/3}

    .. math::
        e = \exp(f + g\omega_{SRK} + h \omega_{SRK}^2)

    .. math::
        C = j + k \omega_{SRK}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    P : float
        Pressure of fluid [Pa]
    Psat : float
        Saturation pressure of the fluid [Pa]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        (ideally SRK) Acentric factor for fluid, [-]
        This parameter is alternatively a fit parameter.
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Returns
    -------
    V_dense : float
        High-pressure liquid volume, [m^3/mol]

    Notes
    -----
    Original equation was in terms of density, but it is converted here.

    The example is from DIPPR, and exactly correct.
    This is DIPPR Procedure 4C: Method for Estimating the Density of Pure
    Organic Liquids under Pressure.

    Examples
    --------
    >>> COSTALD_compressed(303., 9.8E7, 85857.9, 466.7, 3640000.0, 0.281, 0.000105047)
    9.287482879788505e-05

    References
    ----------
    .. [1]  Thomson, G. H., K. R. Brobst, and R. W. Hankinson. "An Improved
       Correlation for Densities of Compressed Liquids and Liquid Mixtures."
       AIChE Journal 28, no. 4 (July 1, 1982): 671-76. doi:10.1002/aic.690280420
    '''
    a = -9.070217
    b = 62.45326
    d = -135.1102
    f = 4.79594
    g = 0.250047
    h = 1.14188
    j = 0.0861488
    k = 0.0344483
    e = exp(f + omega*(g + h*omega))
    C = j + k*omega
    tau = 1.0 - T/Tc
    tau13 = tau**(1.0/3.0)
    B = Pc*(-1.0 + a*tau13 + b*tau13*tau13 + d*tau + e*tau*tau13)
    return Vs*(1.0 - C*log((B + P)/(B + Psat)))


### Liquid Mixtures

def Amgat(xs, Vms):
    r'''Calculate mixture liquid density using the Amgat mixing rule.
    Highly inacurate, but easy to use. Assumes idea liquids with
    no excess volume. Average molecular weight should be used with it to obtain
    density.

    .. math::
        V_{mix} = \sum_i x_i V_i

    or in terms of density:

    .. math::

        \rho_{mix} = \sum\frac{x_i}{\rho_i}

    Parameters
    ----------
    xs : array
        Mole fractions of each component, []
    Vms : array
        Molar volumes of each fluids at conditions [m^3/mol]

    Returns
    -------
    Vm : float
        Mixture liquid volume [m^3/mol]

    Notes
    -----
    Units are that of the given volumes.
    It has been suggested to use this equation with weight fractions,
    but the results have been less accurate.

    Examples
    --------
    >>> Amgat([0.5, 0.5], [4.057e-05, 5.861e-05])
    4.9590000000000005e-05
    '''
    return mixing_simple(xs, Vms)


def Rackett_mixture(T, xs, MWs, Tcs, Pcs, Zrs):
    r'''Calculate mixture liquid density using the Rackett-derived mixing rule
    as shown in [2]_.

    .. math::
        V_m = \sum_i\frac{x_i T_{ci}}{MW_i P_{ci}} Z_{R,m}^{(1 + (1 - T_r)^{2/7})} R \sum_i x_i MW_i

    Parameters
    ----------
    T : float
        Temperature of liquid [K]
    xs: list
        Mole fractions of each component, []
    MWs : list
        Molecular weights of each component [g/mol]
    Tcs : list
        Critical temperatures of each component [K]
    Pcs : list
        Critical pressures of each component [Pa]
    Zrs : list
        Rackett parameters of each component []

    Returns
    -------
    Vm : float
        Mixture liquid volume [m^3/mol]

    Notes
    -----
    Model for pure compounds in [1]_ forms the basis for this model, shown in
    [2]_. Molecular weights are used as weighing by such has been found to
    provide higher accuracy in [2]_. The model can also be used without
    molecular weights, but results are somewhat different.

    As with the Rackett model, critical compressibilities may be used if
    Rackett parameters have not been regressed.

    Critical mixture temperature, and compressibility are all obtained with
    simple mixing rules.

    Examples
    --------
    Calculation in [2]_ for methanol and water mixture. Result matches example.

    >>> Rackett_mixture(T=298., xs=[0.4576, 0.5424], MWs=[32.04, 18.01], Tcs=[512.58, 647.29], Pcs=[8.096E6, 2.209E7], Zrs=[0.2332, 0.2374])
    2.6252894930056885e-05

    References
    ----------
    .. [1] Rackett, Harold G. "Equation of State for Saturated Liquids."
       Journal of Chemical & Engineering Data 15, no. 4 (1970): 514-517.
       doi:10.1021/je60047a012
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    bigsum, Tc, Zr, MW = 0.0, 0.0, 0.0, 0.0

    # Fastest for numba and PyPy and CPython
    for i in range(len(xs)):
        x0 = Tcs[i]*xs[i]
        Tc += x0
        Zr += Zrs[i]*xs[i]
        MW += MWs[i]*xs[i]
        bigsum += x0/(Pcs[i]*MWs[i])
    Tr = T/Tc
    return (R*bigsum*Zr**(1.0 + (1.0 - Tr)**(2.0/7.0)))*MW


def COSTALD_mixture(xs, T, Tcs, Vcs, omegas):
    r'''Calculate mixture liquid density using the COSTALD CSP method.

    A popular and accurate estimation method. If possible, fit parameters are
    used; alternatively critical properties work well.

    The mixing rules giving parameters for the pure component COSTALD
    equation are:

    .. math::
        T_{cm} = \frac{\sum_i\sum_j x_i x_j (V_{ij}T_{cij})}{V_m}

    .. math::
        V_m = 0.25\left[ \sum x_i V_i + 3(\sum x_i V_i^{2/3})(\sum_i x_i V_i^{1/3})\right]

    .. math::
        V_{ij}T_{cij} = (V_iT_{ci}V_{j}T_{cj})^{0.5}

    .. math::
        \omega = \sum_i z_i \omega_i

    Parameters
    ----------
    xs: list
        Mole fractions of each component
    T : float
        Temperature of fluid [K]
    Tcs : list
        Critical temperature of fluids [K]
    Vcs : list
        Critical volumes of fluids [m^3/mol].
        This parameter is alternatively a fit parameter
    omegas : list
        (ideally SRK) Acentric factor of all fluids, [-]
        This parameter is alternatively a fit parameter.

    Returns
    -------
    Vs : float
        Saturation liquid mixture volume

    Notes
    -----
    Range: 0.25 < Tr < 0.95, often said to be to 1.0
    No example has been found.
    Units are that of critical or fit constant volume.

    Examples
    --------
    >>> COSTALD_mixture([0.4576, 0.5424], 298.,  [512.58, 647.29], [0.000117, 5.6e-05], [0.559,0.344])
    2.7065887732713534e-05

    References
    ----------
    .. [1] Hankinson, Risdon W., and George H. Thomson. "A New Correlation for
       Saturated Densities of Liquids and Their Mixtures." AIChE Journal
       25, no. 4 (1979): 653-663. doi:10.1002/aic.690250412
    '''
    N = len(xs)
    sum1, sum2, sum3, omega = 0.0, 0.0, 0.0, 0.0
    for i in range(N):
        sum1 += xs[i]*Vcs[i]
        p = Vcs[i]**(1.0/3.)
        v = xs[i]*p
        sum2 += v
        sum3 += v*p
        omega += xs[i]*omegas[i]

    root_two = 1.4142135623730951 #2.0**0.5
    Vm = 0.25*(sum1 + 3.0*sum2*sum3)
    Vm_inv_root = root_two*(Vm)**-0.5
    vec = [0.0]*N
    for i in range(N):
        vec[i] = (Tcs[i]*Vcs[i])**0.5*xs[i]*Vm_inv_root

    Tcm = 0.0
    for i in range(N):
        for j in range(i):
            Tcm += vec[i]*vec[j]
        Tcm += 0.5*vec[i]*vec[i]
    return COSTALD(T, Tcm, Vm, omega)


### Gases


def ideal_gas(T, P):
    r'''Calculates ideal gas molar volume.
    The molar volume of an ideal gas is given by:

    .. math::
        V = \frac{RT}{P}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    P : float
        Pressure of fluid [Pa]

    Returns
    -------
    V : float
        Gas volume, [m^3/mol]

    Examples
    --------
    >>> ideal_gas(298.15, 101325.)
    0.024465403697038125
    '''
    return R*T/P


### Solids

def Goodman(T, Tt, Vml):
    r'''Calculates solid density at T using the simple relationship
    by a member of the DIPPR.

    The molar volume of a solid is given by:

    .. math::
        \frac{1}{V_m} = \left( 1.28 - 0.16 \frac{T}{T_t}\right)
        \frac{1}{{Vm}_L(T_t)}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tt : float
        Triple temperature of fluid [K]
    Vml : float
        Liquid molar volume of the organic liquid at the triple point,
        [m^3/mol]

    Returns
    -------
    Vms : float
        Solid molar volume, [m^3/mol]

    Notes
    -----
    Works to the next solid transition temperature or to approximately 0.3Tt.

    Examples
    --------
    Decane at 200 K:

    >>> Goodman(200, 243.225, 0.00023585)
    0.0002053665090860923

    References
    ----------
    .. [1] Goodman, Benjamin T., W. Vincent Wilding, John L. Oscarson, and
       Richard L. Rowley. "A Note on the Relationship between Organic Solid
       Density and Liquid Density at the Triple Point." Journal of Chemical &
       Engineering Data 49, no. 6 (2004): 1512-14. doi:10.1021/je034220e.
    '''
    return Vml/(1.28 - 0.16*(T/Tt))


