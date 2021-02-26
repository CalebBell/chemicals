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


This module contains various thermal conductivity estimation routines, dataframes
of fit coefficients, and mixing rules.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Pure Low Pressure Liquid Correlations
-------------------------------------
.. autofunction:: chemicals.thermal_conductivity.Sheffy_Johnson
.. autofunction:: chemicals.thermal_conductivity.Sato_Riedel
.. autofunction:: chemicals.thermal_conductivity.Lakshmi_Prasad
.. autofunction:: chemicals.thermal_conductivity.Gharagheizi_liquid
.. autofunction:: chemicals.thermal_conductivity.Nicola_original
.. autofunction:: chemicals.thermal_conductivity.Nicola
.. autofunction:: chemicals.thermal_conductivity.Bahadori_liquid
.. autofunction:: chemicals.thermal_conductivity.kl_Mersmann_Kind

Pure High Pressure Liquid Correlations
--------------------------------------
.. autofunction:: chemicals.thermal_conductivity.DIPPR9G
.. autofunction:: chemicals.thermal_conductivity.Missenard

Liquid Mixing Rules
-------------------
.. autofunction:: chemicals.thermal_conductivity.DIPPR9H
.. autofunction:: chemicals.thermal_conductivity.DIPPR9I
.. autofunction:: chemicals.thermal_conductivity.Filippov

Pure Low Pressure Gas Correlations
----------------------------------
.. autofunction:: chemicals.thermal_conductivity.Eucken
.. autofunction:: chemicals.thermal_conductivity.Eucken_modified
.. autofunction:: chemicals.thermal_conductivity.DIPPR9B
.. autofunction:: chemicals.thermal_conductivity.Chung
.. autofunction:: chemicals.thermal_conductivity.Eli_Hanley
.. autofunction:: chemicals.thermal_conductivity.Gharagheizi_gas
.. autofunction:: chemicals.thermal_conductivity.Bahadori_gas

Pure High Pressure Gas Correlations
-----------------------------------
.. autofunction:: chemicals.thermal_conductivity.Stiel_Thodos_dense
.. autofunction:: chemicals.thermal_conductivity.Eli_Hanley_dense
.. autofunction:: chemicals.thermal_conductivity.Chung_dense

Gas Mixing Rules
----------------
.. autofunction:: chemicals.thermal_conductivity.Lindsay_Bromley
.. autofunction:: chemicals.thermal_conductivity.Wassiljewa_Herning_Zipperer

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.thermal_conductivity.k_IAPWS
.. autofunction:: chemicals.thermal_conductivity.k_air_lemmon


Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: k_data_Perrys_8E_2_315

    Data from [1]_ with :obj:`chemicals.dippr.EQ100` coefficients for liquids.

.. data:: k_data_Perrys_8E_2_314

    Data from [1]_ with :obj:`chemicals.dippr.EQ102` coefficients for gases.

.. data:: k_data_VDI_PPDS_9

    Data from [2]_ with polynomial coefficients for liquids.

.. data:: k_data_VDI_PPDS_10

    Data from [2]_ with polynomial coefficients for gases.

.. [1] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
   Eighth Edition. McGraw-Hill Professional, 2007.
.. [2] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
   Berlin; New York:: Springer, 2010.

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.thermal_conductivity.k_data_Perrys_8E_2_315

    In [3]: chemicals.thermal_conductivity.k_data_Perrys_8E_2_314

    In [4]: chemicals.thermal_conductivity.k_data_VDI_PPDS_9

    In [5]: chemicals.thermal_conductivity.k_data_VDI_PPDS_10

"""

from __future__ import division

__all__ = ['Sheffy_Johnson', 'Sato_Riedel', 'Lakshmi_Prasad',
'Gharagheizi_liquid', 'Nicola_original', 'Nicola', 'Bahadori_liquid',
'kl_Mersmann_Kind', 'DIPPR9G', 'DIPPR9I','k_IAPWS',
'Missenard', 'DIPPR9H', 'Filippov', 'Eucken', 'Eucken_modified', 'DIPPR9B',
'Chung', 'Eli_Hanley', 'Gharagheizi_gas', 'Bahadori_gas',
'Stiel_Thodos_dense', 'Eli_Hanley_dense', 'Chung_dense', 'Lindsay_Bromley',
'Wassiljewa_Herning_Zipperer', 'k_air_lemmon']

from fluids.numerics import bisplev, implementation_optimize_tck, numpy as np
from fluids.constants import R, R_inv, N_A, k, pi
from chemicals.utils import log, exp, sqrt, atan, PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import register_df_source, data_source
from chemicals.viscosity import Herning_Zipperer

folder = os_path_join(source_path, 'Thermal Conductivity')


register_df_source(folder, 'Table 2-314 Vapor Thermal Conductivity of Inorganic and Organic Substances.tsv')
register_df_source(folder, 'Table 2-315 Thermal Conductivity of Inorganic and Organic Liquids.tsv')
register_df_source(folder, 'VDI PPDS Thermal conductivity of saturated liquids.tsv', csv_kwargs={'float_precision': 'legacy'})
register_df_source(folder, 'VDI PPDS Thermal conductivity of gases.tsv', csv_kwargs={'float_precision': 'legacy'})

_k_data_loaded = False
def _load_k_data():
    global _k_data_loaded, k_data_Perrys_8E_2_314, k_values_Perrys_8E_2_314
    global k_data_Perrys_8E_2_315, k_values_Perrys_8E_2_315, k_data_VDI_PPDS_9
    global k_values_VDI_PPDS_9, k_data_VDI_PPDS_10, k_values_VDI_PPDS_10

    k_data_Perrys_8E_2_314 = data_source('Table 2-314 Vapor Thermal Conductivity of Inorganic and Organic Substances.tsv')
    k_values_Perrys_8E_2_314 = np.array(k_data_Perrys_8E_2_314.values[:, 1:], dtype=float)

    k_data_Perrys_8E_2_315 = data_source('Table 2-315 Thermal Conductivity of Inorganic and Organic Liquids.tsv')
    k_values_Perrys_8E_2_315 = np.array(k_data_Perrys_8E_2_315.values[:, 1:], dtype=float)

    k_data_VDI_PPDS_9 = data_source('VDI PPDS Thermal conductivity of saturated liquids.tsv')
    k_values_VDI_PPDS_9 = np.array(k_data_VDI_PPDS_9.values[:, 1:], dtype=float)

    k_data_VDI_PPDS_10 = data_source('VDI PPDS Thermal conductivity of gases.tsv')
    k_values_VDI_PPDS_10 = np.array(k_data_VDI_PPDS_10.values[:, 1:], dtype=float)

if PY37:
    def __getattr__(name):
        if name in ('k_data_Perrys_8E_2_314', 'k_values_Perrys_8E_2_314', 'k_data_Perrys_8E_2_315',
                    'k_values_Perrys_8E_2_315', 'k_data_VDI_PPDS_9', 'k_values_VDI_PPDS_9', 'k_data_VDI_PPDS_10',
                    'k_values_VDI_PPDS_10'):
            _load_k_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_k_data()

pi_inv = 1.0/pi # todo move to fluids.constants

def k_IAPWS(T, rho, Cp=None, Cv=None, mu=None, drho_dP=None, drho_dP_Tr=None):
    r'''Calculate the thermal conductivity of water or steam according to the
    2011 IAPWS [1]_ formulation. Critical enhancement is ignored unless
    parameters for it are provided.

    .. math::
        \bar\lambda = \bar\lambda_0\times \bar\lambda_1(\bar T, \bar \rho)
        + \bar\lambda_2(\bar T, \bar\rho)

    .. math::
        \bar\lambda_0 = \frac{\sqrt{\bar T}}
        {\sum_{k=0}^4 \frac{L_k}{\bar T^k}}

    .. math::
        \bar \lambda_1(\bar T, \bar \rho) = \exp\left[ \bar\rho \sum_{i=0}^4
        \left(\left(\frac{1}{\bar T} - 1 \right)^i
        \sum_{j=0}^5 L_{ij}(\bar\rho - 1)^j\right)\right]

    .. math::
        \bar\lambda_2 = \Gamma\frac{\bar\rho \bar c_p \bar T}{\bar \mu} Z(y)

    .. math::
        Z(y) = \frac{2}{\pi y} \left\{\left[(1 - \kappa^{-1})\arctan(y)
        + \kappa^{-1}y\right] - \left[1 - \exp\left(\frac{-1}{y^{-1}
        + y^{-2}/3\bar\rho^2}\right)\right]\right\}

    .. math::
        y = \bar q_D \xi(\bar T, \bar \rho)

    .. math::
        \xi = \xi_0 \left(\frac{\Delta \bar\chi}{\Gamma_0}\right)^{\nu/\gamma}

    .. math::
        \Delta \bar\chi(\bar T, \bar \rho) = \bar\rho\left[
        \zeta(\bar T, \bar \rho) - \zeta(\bar T_R, \bar \rho)\frac{\bar T_R}{\bar T}
        \right]

    .. math::
        \zeta = \left(\frac{\partial \bar \rho}{\partial \bar p}\right)_{\bar T}

    Parameters
    ----------
    T : float
        Temperature water [K]
    rho : float
        Density of water [kg/m^3]
    Cp : float, optional
        Constant pressure heat capacity of water, [J/kg/K]
    Cv : float, optional
        Constant volume heat capacity of water, [J/kg/K]
    mu : float, optional
        Viscosity of water, [Pa*s]
    drho_dP : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature, [kg/m^3/Pa]
    drho_dP_Tr : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature (at the reference temperature (970.644 K) and the actual
        density of water); will be calculated from the industrial formulation
        fit if omitted, [kg/m^3/Pa]

    Returns
    -------
    k : float
        Thermal condiuctivity, [W/m/K]

    Notes
    -----
    Gamma = 177.8514;

    qd = 0.4E-9;

    nu = 0.630;

    gamma = 1.239;

    zeta0 = 0.13E-9;

    Gamma0 = 0.06;

    TRC = 1.5

    The formulation uses the industrial variant of the critical enhancement.
    It matches to 5E-6 relative tolerance at the check temperature, and should
    match even closer outside it.

    Examples
    --------
    >>> k_IAPWS(647.35, 750.)
    0.5976194153179502

    Region 1, test 1, from MPEI, exact match:

    >>> k_IAPWS(T=620., rho=613.227777440324, Cp=7634.337046792,
    ... Cv=3037.934412104, mu=70.905106751524E-6, drho_dP=5.209378197916E-6)
    0.48148519510200044

    Full scientific calculation:

    >>> from chemicals.iapws import iapws95_properties, iapws95_P, iapws95_Tc
    >>> from chemicals.viscosity import mu_IAPWS
    >>> T, P = 298.15, 1e5
    >>> rho, _, _, _, Cv, Cp, _, _, _, _, drho_dP = iapws95_properties(T, P)
    >>> P_ref = iapws95_P(1.5*iapws95_Tc, rho)
    >>> _, _, _, _, _, _, _, _, _, _, drho_dP_Tr = iapws95_properties(1.5*iapws95_Tc, P_ref)
    >>> mu = mu_IAPWS(T, rho, drho_dP, drho_dP_Tr)
    >>> k_IAPWS(T, rho, Cp, Cv, mu, drho_dP, drho_dP_Tr)
    0.60651532815

    References
    ----------
    .. [1] Huber, M. L., R. A. Perkins, D. G. Friend, J. V. Sengers, M. J.
       Assael, I. N. Metaxa, K. Miyagawa, R. Hellmann, and E. Vogel. "New
       International Formulation for the Thermal Conductivity of H2O."
       Journal of Physical and Chemical Reference Data 41, no. 3 (September 1,
       2012): 033102. doi:10.1063/1.4738955.
    '''
    rhor = rho*0.003105590062111801#1/322.0
    Tr = T*0.0015453657571674064 # 1/647.096
    Tr_inv = 1.0/Tr

#     Lijs = [[1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258],
#             [2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245],
#             [2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816],
#             [-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0, 0],
#             [-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842]]

#     Aijs = [[6.53786807199516, 6.52717759281799, 5.35500529896124, 1.55225959906681, 1.11999926419994],
#             [-5.61149954923348, -6.30816983387575, -3.96415689925446, 0.464621290821181, 0.595748562571649],
#             [3.39624167361325, 8.08379285492595, 8.91990208918795, 8.93237374861479, 9.8895256507892],
#             [-2.27492629730878, -9.82240510197603, -12.033872950579, -11.0321960061126, -10.325505114704],
#             [10.2631854662709, 12.1358413791395, 9.19494865194302, 6.1678099993336, 4.66861294457414],
#             [1.97815050331519, -5.54349664571295, -2.16866274479712, -0.965458722086812, -0.503243546373828]]
    '''Unoptimized (but editable) code; the below is generated with sympy
    Ls = [2.443221E-3, 1.323095E-2, 6.770357E-3, -3.454586E-3, 4.096266E-4]
    lambda0 = 0
    for i, L in enumerate(Ls):
        lambda0 += L/Tr**i
    lambda0 = Tr**0.5/lambda0
    '''
    lambda0 = sqrt(Tr)/(Tr_inv*(Tr_inv*(Tr_inv*(0.0004096266*Tr_inv - 0.003454586)
                        + 0.006770357) + 0.01323095) + 0.002443221)

    '''Unoptimized (but editable) code; the below is generated with sympy
    tot1 = 0
    for i, Ljs in enumerate(Lijs):
        tot2 = 0
        for j, L in enumerate(Ljs):
            tot2 += L*(rhor - 1.)**j
        tot1 += (1./Tr -1.)**i*tot2
    '''
    x0 = rhor - 1.0
    x1 = (Tr_inv - 1.0)
    x12 = x1*x1
    tot1 = (x0*(x0*(x0*(x0*(x0*(x1*(x1*(0.012913842*x12 - 0.0205938816) - 0.00719201245) + 0.00609859258)
                            + x1*(x1*(0.275418278 - 0.19268305*x12) + 0.0832827019) - 0.0504123634)
                        + x1*(x1*(x1*(1.1168348*x1 + 0.0716373224) - 1.40944978) - 0.463045512) + 0.102997357)
                    + x1*(x1*(x1*(-3.18369245*x1 - 0.621178141) + 3.55777244) + 1.53616167) + 0.111443906)
                + x1*(x1*(x1*(4.57586331*x1 + 1.60812989) - 4.54580785) - 2.78843778) - 0.646013523)
            + x1*(x1*(x1*(-2.720337*x1 - 1.21051378) + 2.19650529) + 2.33771842) + 1.60397357)
    lambda1 = exp(rhor*tot1)

    if Cp is not None and Cv is not None and mu is not None and drho_dP is not None:
        Cpr = Cp*0.0021667624917378636 #1/461.51805 # J/kg/K
        if Cpr < 0.0 or Cpr > 1E13:
            Cpr = 1E13
            Cp = Cpr*461.51805 # This is correct
        mur = mu*1E6
        kappa_inv = Cv/Cp

        Gamma = 177.8514
        qd = 2500000000.0#(0.4E-9)**-1
        # nu = 0.630
        xi0 = 0.13E-9
        # gamma = 1.239
        # Gamma0 = 0.06
        TRC = 1.5

        zeta_drho_dP = drho_dP*68521.73913043478#22.064E6/322.0
        if drho_dP_Tr is None:
            if rhor <= 0.310559006:
                tot1 = (rhor*(rhor*(rhor*(rhor*(1.97815050331519*rhor + 10.2631854662709) - 2.27492629730878)
                            + 3.39624167361325) - 5.61149954923348) + 6.53786807199516)
            elif rhor <= 0.776397516:
                tot1 = (rhor*(rhor*(rhor*(rhor*(12.1358413791395 - 5.54349664571295*rhor) - 9.82240510197603)
                            + 8.08379285492595) - 6.30816983387575) + 6.52717759281799)
            elif rhor <= 1.242236025:
                tot1 = (rhor*(rhor*(rhor*(rhor*(9.19494865194302 - 2.16866274479712*rhor) - 12.033872950579)
                            + 8.91990208918795) - 3.96415689925446) + 5.35500529896124)
            elif rhor <= 1.863354037:
                tot1 = (rhor*(rhor*(rhor*(rhor*(6.1678099993336 - 0.965458722086812*rhor) - 11.0321960061126)
                            + 8.93237374861479) + 0.464621290821181) + 1.55225959906681)
            else:
                tot1 = (rhor*(rhor*(rhor*(rhor*(4.66861294457414 - 0.503243546373828*rhor) - 10.325505114704)
                            + 9.8895256507892) + 0.595748562571649) + 1.11999926419994)
            '''Original code:
            zeta_drho_dP_Tr = 1./sum([Aijs[i][j]*rhor**i for i in range(6)])
            '''
            zeta_drho_dP_Tr = 1./tot1
        else:
            zeta_drho_dP_Tr = drho_dP_Tr*68521.73913043478#22.064E6/322.0
        dchi = rhor*(zeta_drho_dP - zeta_drho_dP_Tr*TRC*Tr_inv)
        if dchi < 0.0:
            xi = 0.0
        else:
            # 16.666666666666668 = 1.0/Gamma0
            xi = xi0*(dchi*16.666666666666668)**0.5084745762711864#(nu/gamma)
        if xi < 0.0 or xi > 1E4:
            xi = 1E4

        y = qd*xi
        if y < 1.2E-7:
            Z_y = 0.0
        else:
            y_inv = 1.0 / y
            Z_y = 2.*pi_inv*y_inv*(((1.0 - kappa_inv)*atan(y) + kappa_inv*y)
                               - (1.0 - exp(-1.0/(y_inv + y*y/(3.0*rhor*rhor)))))
        lambda2 = Gamma*rhor*Cpr*Tr/mur*Z_y

    else:
        lambda2 = 0.0

    k = (lambda0*lambda1 + lambda2)*1e-3
    return k

def k_air_lemmon(T, rho, Cp=None, Cv=None, drho_dP=None, drho_dP_Tr=None, mu=None):
    r'''Calculate the thermal conductivity of air using the Lemmon and Jacobsen
    (2004) [1]_ formulation. The critical enhancement term is ignored unless
    all the rquired parameters for it are provided.

    .. math::
        \lambda = \lambda^0(T) + \lambda^r(\tau, \delta) + \lambda^c(\tau, \delta)

    .. math::
        \lambda^0 = N_1\left[\frac{\eta^0(T)}{1 \mu \text{Pa}\cdot \text{s}}
        \right] + N_2\tau^{t_2} + N_3\tau^{t_3}

    .. math::
        \lambda^r = \sum_{i=4}^n N_i \tau^{t_i} \delta^{d_i} \exp(-\gamma_i
        \delta^{l_i})

    .. math::
        \lambda^c = \rho C_p \frac{kR_0 T}{6\pi\xi\cdot \eta(T, \rho)}\left(
        \tilde \Omega -\tilde \Omega_0\right)

    .. math::
        \tilde \Omega = \frac{2}{\pi}\left[
        \left(\frac{C_p - C_v}{C_p}\right)\tan^{-1} (\xi/q_D) + \frac{C_v}
        {C_p}(\xi/q_D) \right]

    .. math::
        \tilde \Omega_0 = \frac{2}{\pi}\left\{1 - \exp\left[\frac{-1}{q_D/\xi
        + 1/3(\xi/q_D)^2(\rho_c/\rho)^2} \right] \right\}

    .. math::
        \xi = \xi_0 \left[\frac{\tilde \chi(T, \rho) - \tilde \chi(T_{ref},
        \rho)\frac{T_{ref}}{T}}{\Gamma}  \right]^{\nu/\gamma}

    .. math::
        \tilde \chi(T, \rho) = \frac{P_c \rho}{\rho_c^2} \left(\frac{\partial
        \rho}{\partial P} \right)_{T}

    Parameters
    ----------
    T : float
        Temperature air [K]
    rho : float
        Molar density of air [mol/m^3]
    Cp : float, optional
        Molar constant pressure heat capacity of air, [J/mol/K]
    Cv : float, optional
        Molar constant volume heat capacity of air, [J/mol/K]
    mu : float, optional
        Viscosity of air, [Pa*s]
    drho_dP : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature, [mol/m^3/Pa]
    drho_dP_Tr : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature (at the reference temperature (265.262 K) and the actual
        density of air), [mol/m^3/Pa]

    Returns
    -------
    k : float
        Thermal condiuctivity of air, [W/m/K]

    Notes
    -----
    The constnts are as follows:

    Ni = [1.308, 1.405, -1.036, 8.743, 14.76, -16.62, 3.793, -6.142, -0.3778]

    ti = [None, -1.1, -0.3, 0.1, 0.0, 0.5, 2.7, 0.3, 1.3]

    di = [None, None, None, 1, 2, 3, 7, 7, 11]

    li  = [None, None, None, 0, 0, 2, 2, 2, 2]

    gammai = [None, None, None, 0, 0, 1, 1, 1, 1]

    R0 = 1.01; Pc = 3.78502E6 Pa; xi0 = 0.11E-9 nm; qd = 0.31E-9 nm;
    Tc = 132.6312 K (actually the maxcondentherm); T_ref = 265.262 (2Tc
    rounded differently); rhoc = 10447.7 mol/m^3 (actually the maxcondentherm);
    k = 1.380658E-23 J/K; nu = 0.63 and gamma = 1.2415, sigma = 0.36,
    MW = 28.9586 g/mol.


    Examples
    --------
    Basic calculation at 300 K and approximately 1 bar:

    >>> k_air_lemmon(300, 40.0)
    0.0263839695044

    Calculation near critical point:

    >>> k_air_lemmon(132.64, 10400, 2137.078854678728, 35.24316159996235, 0.07417878614315769, 0.00035919027241528256, 1.7762253265868595e-05)
    0.07562307234760

    References
    ----------
    .. [1] Lemmon, E. W., and R. T. Jacobsen. "Viscosity and Thermal
       Conductivity Equations for Nitrogen, Oxygen, Argon, and Air."
       International Journal of Thermophysics 25, no. 1 (January 1, 2004):
       21-69. https://doi.org/10.1023/B:IJOT.0000022327.04529.f3.
    '''
    R0 = 1.01
    Pc = 3.78502E6
    xi0 = 0.11E-9
    qd = 0.31E-9
#     Gamma = 0.055
    Tc = 132.6312 # K, maxcondentherm actually
    T_ref = 265.262 # Tc*2 but rounded differently
    rhoc = 10447.7
    rhoc2 = rhoc*rhoc

    qd_inv = 3225806451.612903 # 10.31E-9
    gamma_inv = 18.181818181818183 # 1/.055

    tau = Tc/T
    tau_10 = tau**0.1
    tau2_10 = tau_10*tau_10
    tau3_10 = tau_10*tau2_10
    tau12_10 = tau3_10*tau3_10
    tau12_10 *= tau12_10
    tau24_10 = tau12_10*tau12_10

    delta = rho*9.571484632981421e-05 # 9.57...E-5 = 1/10447.7

    Ts = T*0.00968054211035818 # 1/e_k
    lnTs = log(Ts)
    Omega_inv = exp(-0.431 - lnTs*(lnTs*(lnTs*(0.005341 - 0.00331*lnTs) + 0.08406) - 0.4623))

    #12.7658... = 0.0266958*sqrt(28.9586)/(0.360*0.360)*sqrt(132.6312)
    eta0 = 12.765845058845755*Omega_inv/(tau2_10*tau3_10)

    k0 = 1.308*eta0 + 1.405/(tau*tau_10) - 1.036/tau3_10

#     kr = 0.0
#     for i in range(3, 9):
#         kr += Ni[i]*tau**ti[i]*delta**di[i]*exp(-gammai[i]*delta**li[i])

    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta*delta3
    kr = (8.743*delta*tau_10 + 14.76*delta2 -  exp(-delta2)*(16.62*delta3*tau3_10*tau2_10
          + delta3*delta4*(0.3778*delta4*tau_10*tau12_10 + tau3_10*(6.142 - 3.793*tau24_10))))

    if Cp is not None and Cv is not None and mu is not None and drho_dP is not None and drho_dP_Tr is not None:
        x2 = Pc*rho/rhoc2
        xi_bar = x2*drho_dP
        xi_bar_ref = x2*drho_dP_Tr

        xi = xi0*((xi_bar - xi_bar_ref*T_ref/T)*gamma_inv)**0.5074506645187273# .50745... = (0.63/1.2415)
        if xi < 0.0:
            kc = 0.0
        else:
            xi_qd = xi*qd_inv

            term0 = qd/xi + (1.0/3.0)*xi_qd*xi_qd*(rhoc*rhoc/(rho*rho))
            Omega_bar0 = 2.0*pi_inv*(1.0 - exp(-1.0/term0))

            Omega_bar = 2.0*pi_inv*((Cp - Cv)/Cp*atan(xi_qd) + Cv/Cp*xi_qd)
            k = 1.380658E-23 # J/K

            # Mu should still be in Pa*s
            kc = rho*Cp*k*R0*T/(6.0*pi*xi*mu)*(Omega_bar - Omega_bar0)
            kc *= 1e3 # Convert to mW/m/K, same as others
        return (k0 + kr + kc)*1e-3
    return (k0 + kr)*1e-3

### Purely CSP Methods - Liquids

def Sheffy_Johnson(T, MW, Tm):
    r'''Calculate the thermal conductivity of a liquid as a function of
    temperature using the Sheffy-Johnson (1961) method. Requires
    Temperature, molecular weight, and melting point.

    .. math::
        k = 1.951 \frac{1-0.00126(T-T_m)}{T_m^{0.216}MW^{0.3}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tm : float
        Melting point of the fluid [K]

    Returns
    -------
    kl : float
        Thermal conductivity of the fluid, W/m/k

    Notes
    -----
    The origin of this equation has been challenging to trace. It is
    presently unknown, and untested.

    Examples
    --------
    >>> Sheffy_Johnson(300, 47, 280)
    0.17740150413112193

    References
    ----------
    .. [1] Scheffy, W. J., and E. F. Johnson. "Thermal Conductivities of
       Liquids at High Temperatures." Journal of Chemical & Engineering Data
       6, no. 2 (April 1, 1961): 245-49. doi:10.1021/je60010a019
    '''
    return 1.951*(1.0 - 0.00126*(T - Tm))*Tm**-0.216*MW**-0.3


def Sato_Riedel(T, MW, Tb, Tc):
    r'''Calculate the thermal conductivity of a liquid as a function of
    temperature using the CSP method of Sato-Riedel [1]_, [2]_, published in
    Reid [3]_. Requires temperature, molecular weight, and boiling and critical
    temperatures.

    .. math::
        k = \frac{1.1053}{\sqrt{MW}}\frac{3+20(1-T_r)^{2/3}}
        {3+20(1-T_{br})^{2/3}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    This equation has a complicated history. It is proposed by Reid [3]_.
    Limited accuracy should be expected. Uncheecked.

    Examples
    --------
    >>> Sato_Riedel(300, 47, 390, 520)
    0.21037692461337687

    References
    ----------
    .. [1] Riedel, L.: Chem. Ing. Tech., 21, 349 (1949); 23: 59, 321, 465 (1951)
    .. [2] Maejima, T., private communication, 1973
    .. [3] Properties of Gases and Liquids", 3rd Ed., McGraw-Hill, 1977
    '''
    Tr = T/Tc
    Tbr = Tb/Tc
    return 1.1053*(3. + 20.*(1.0 - Tr)**(2.0/3.0))*MW**-0.5/(3.0 + 20.0*(1.0 - Tbr)**(2.0/3.))


def Lakshmi_Prasad(T, MW):
    r'''Estimates thermal conductivity of pure liquids as a function of
    temperature using a reference fluid approach. Low accuracy but quick.
    Developed using several organic fluids.

    .. math::
        \lambda = 0.0655-0.0005T + \frac{1.3855-0.00197T}{MW^{0.5}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    This equation returns negative numbers at high T sometimes.
    This equation is one of those implemented by DDBST.
    If this results in a negative thermal conductivity, no value is returned.

    Examples
    --------
    >>> Lakshmi_Prasad(273.15, 100)
    0.013664450000000009

    References
    ----------
    .. [1] Lakshmi, D. S., and D. H. L. Prasad. "A Rapid Estimation Method for
       Thermal Conductivity of Pure Liquids." The Chemical Engineering Journal
       48, no. 3 (April 1992): 211-14. doi:10.1016/0300-9467(92)80037-B
    '''
    return 0.0655 - 0.0005*T + (1.3855 - 0.00197*T)*MW**-0.5


def Gharagheizi_liquid(T, MW, Tb, Pc, omega):
    r'''Estimates the thermal conductivity of a liquid as a function of
    temperature using the CSP method of Gharagheizi [1]_. A  convoluted
    method claiming high-accuracy and using only statistically significant
    variable following analalysis.

    Requires temperature, molecular weight, boiling temperature and critical
    pressure and acentric factor.

    .. math::
        k = 10^{-4}\left[10\omega + 2P_c-2T+4+1.908(T_b+\frac{1.009B^2}{MW^2})
        +\frac{3.9287MW^4}{B^4}+\frac{A}{B^8}\right]

    .. math::
        A = 3.8588MW^8(1.0045B+6.5152MW-8.9756)

    .. math::
        B = 16.0407MW+2T_b-27.9074

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tb : float
        Boiling temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor of the fluid [-]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    Pressure is internally converted into bar, as used in the original equation.

    This equation was derived with 19000 points representing 1640 unique compounds.

    Examples
    --------
    >>> Gharagheizi_liquid(300, 40, 350, 1E6, 0.27)
    0.2171113029534838

    References
    ----------
    .. [1] Gharagheizi, Farhad, Poorandokht Ilani-Kashkouli, Mehdi Sattari,
        Amir H. Mohammadi, Deresh Ramjugernath, and Dominique Richon.
        "Development of a General Model for Determination of Thermal
        Conductivity of Liquid Chemical Compounds at Atmospheric Pressure."
        AIChE Journal 59, no. 5 (May 1, 2013): 1702-8. doi:10.1002/aic.13938
    '''
    M2 = MW*MW
    M4 = M2*M2
    Pc = Pc*1E-5
    B = 16.0407*MW + 2.*Tb - 27.9074
    A = 3.8588*M4*M4*(1.0045*B + 6.5152*MW - 8.9756)
    B_inv4 = 1.0/(B*B*B*B)
    kl = 1E-4*(10.*omega + 2.*(Pc - T) + 4. + 1.908*(Tb + 1.009*B*B/(M2))
        + 3.9287*M4*B_inv4 + A*B_inv4*B_inv4)
    return kl


def Nicola_original(T, MW, Tc, omega, Hfus):
    r'''Estimates the thermal conductivity of a liquid as a function of
    temperature using the CSP method of Nicola [1]_. A  simpler but long
    method claiming high-accuracy and using only statistically significant
    variable following analalysis.

    Requires temperature, molecular weight, critical temperature, acentric
    factor and the heat of vaporization.

    .. math::
        \frac{\lambda}{1 \text{Wm/K}}=-0.5694-0.1436T_r+5.4893\times10^{-10}
        \frac{\Delta_{fus}H}{\text{kmol/J}}+0.0508\omega +
        \left(\frac{1 \text{kg/kmol}}{MW}\right)^{0.0622}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tc : float
        Critical temperature of the fluid [K]
    omega : float
        Acentric factor of the fluid [-]
    Hfus : float
        Heat of fusion of the fluid [J/mol]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    A weird statistical correlation. Recent and yet to be reviewed.
    This correlation has been superceded by the author's later work.
    Hfus is internally converted to be in J/kmol.

    Examples
    --------
    >>> Nicola_original(300, 142.3, 611.7, 0.49, 201853)
    0.2305018632230984

    References
    ----------
    .. [1] Nicola, Giovanni Di, Eleonora Ciarrocchi, Mariano Pierantozzi, and
        Roman Stryjek. "A New Equation for the Thermal Conductivity of Organic
        Compounds." Journal of Thermal Analysis and Calorimetry 116, no. 1
        (April 1, 2014): 135-40. doi:10.1007/s10973-013-3422-7
    '''
    Tr = T/Tc
    Hfus = Hfus*1000.0
    return -0.5694 - 0.1436*Tr + 5.4893E-10*Hfus + 0.0508*omega + MW**-0.0622


def Nicola(T, MW, Tc, Pc, omega):
    r'''Estimates the thermal conductivity of a liquid as a function of
    temperature using the CSP method of [1]_. A statistically derived
    equation using any correlated terms.

    Requires temperature, molecular weight, critical temperature and pressure,
    and acentric factor.

    .. math::
        \frac{\lambda}{0.5147 W/m/K} = -0.2537T_r+\frac{0.0017Pc}{\text{bar}}
        +0.1501 \omega + \left(\frac{1}{MW}\right)^{-0.2999}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor of the fluid [-]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    A statistical correlation. A revision of an original correlation.

    Examples
    --------
    >>> Nicola(300, 142.3, 611.7, 2110000.0, 0.49)
    0.10863821554584034

    References
    ----------
    .. [1] Di Nicola, Giovanni, Eleonora Ciarrocchi, Gianluca Coccia, and
       Mariano Pierantozzi. "Correlations of Thermal Conductivity for
       Liquid Refrigerants at Atmospheric Pressure or near Saturation."
       International Journal of Refrigeration. 2014.
       doi:10.1016/j.ijrefrig.2014.06.003
    '''
    return 0.5147*(-0.2537*T/Tc + 0.0017E-5*Pc + 0.1501*omega + MW**-0.2999)


def Bahadori_liquid(T, MW):
    r'''Estimates the thermal conductivity of parafin liquid hydrocarbons.
    Fits their data well, and is useful as only MW is required.
    X is the Molecular weight, and Y the temperature.

    .. math::
        K = a + bY + CY^2 + dY^3

    .. math::
        a = A_1 + B_1 X + C_1 X^2 + D_1 X^3

    .. math::
        b = A_2 + B_2 X + C_2 X^2 + D_2 X^3

    .. math::
        c = A_3 + B_3 X + C_3 X^2 + D_3 X^3

    .. math::
        d = A_4 + B_4 X + C_4 X^2 + D_4 X^3

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    The accuracy of this equation has not been reviewed.

    Examples
    --------
    Data point from [1]_.

    >>> Bahadori_liquid(273.15, 170)
    0.1427427810827268

    References
    ----------
    .. [1] Bahadori, Alireza, and Saeid Mokhatab. "Estimating Thermal
       Conductivity of Hydrocarbons." Chemical Engineering 115, no. 13
       (December 2008): 52-54
    '''
    A = (-6.48326E-2, 2.715015E-3, -1.08580E-5, 9.853917E-9)
    B = (1.565612E-2, -1.55833E-4, 5.051114E-7, -4.68030E-10)
    C = (-1.80304E-4, 1.758693E-6, -5.55224E-9, 5.201365E-12)
    D = (5.880443E-7, -5.65898E-9, 1.764384E-11, -1.65944E-14)
    X, Y = MW, T

    a = A[0] + X*(B[0] + X*(C[0] + D[0]*X))
    b = A[1] + X*(B[1] + X*(C[1] + D[1]*X))
    c = A[2] + X*(B[2] + X*(C[2] + D[2]*X))
    d = A[3] + X*(B[3] + X*(C[3] + D[3]*X))
    return a + Y*(b + Y*(c + d*Y))


def kl_Mersmann_Kind(T, MW, Tc, Vc, na):
    r'''Estimates the thermal conductivity of organic liquid substances
    according to the method of [1]_.

    .. math::
        \lambda^* = \frac{\lambda\cdot V_c^{2/3}\cdot T_c\cdot \text{MW}^{0.5}}
        {(k\cdot T_c)^{1.5}\cdot N_A^{7/6}}

    .. math::
        \lambda^* = \frac{2}{3}\left(n_a + 40\sqrt{1-T_r}\right)

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW : float
        Molecular weight of the fluid [g/mol]
    Tc : float
        Critical temperature of the fluid [K]
    Vc : float
        Critical volume of the fluid [m^3/mol]
    na : float
        Number of atoms in the molecule, [-]

    Returns
    -------
    kl : float
        Estimated liquid thermal conductivity [W/m/k]

    Notes
    -----
    In the equation, all quantities must be in SI units but N_A is in a kmol
    basis and Vc is in units of (m^3/kmol); this is converted internally.

    Examples
    --------
    Dodecane at 400 K:

    >>> kl_Mersmann_Kind(400, 170.33484, 658.0,
    ... 0.000754, 38)
    0.0895271829899285

    References
    ----------
    .. [1] Mersmann, Alfons, and Matthias Kind. "Prediction of Mechanical and
       Thermal Properties of Pure Liquids, of Critical Data, and of Vapor
       Pressure." Industrial & Engineering Chemistry Research, January 31,
       2017. https://doi.org/10.1021/acs.iecr.6b04323.
    '''
    lambda_star = (2/3.)*(na + 40.*sqrt(1. - T/Tc))
    Vc = Vc*1000.0 # m^3/mol to m^3/kmol
    N_A2 = N_A*1000.0 # Their avogadro's constant is per kmol
    kTc = k*Tc
    kl = lambda_star*kTc*sqrt(kTc/MW)*N_A2**(7.0/6.)*Vc**(-2.0/3.)/Tc
    return kl




### Thermal Conductivity of Dense Liquids

def DIPPR9G(T, P, Tc, Pc, kl):
    r'''Adjustes for pressure the thermal conductivity of a liquid using an
    emperical formula based on [1]_, but as given in [2]_.

    .. math::
        k = k^* \left[ 0.98 + 0.0079 P_r T_r^{1.4} + 0.63 T_r^{1.2}
        \left( \frac{P_r}{30 + P_r}\right)\right]

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    P : float
        Pressure of fluid [Pa]
    Tc: float
        Critical point of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    kl : float
        Thermal conductivity of liquid at 1 atm or saturation, [W/m/K]

    Returns
    -------
    kl_dense : float
        Thermal conductivity of liquid at P, [W/m/K]

    Notes
    -----
    This equation is entrely dimensionless; all dimensions cancel.
    The original source has not been reviewed.

    This is DIPPR Procedure 9G: Method for the Thermal Conductivity of Pure
    Nonhydrocarbon Liquids at High Pressures

    Examples
    --------
    From [2]_, for butyl acetate.

    >>> DIPPR9G(515.05, 3.92E7, 579.15, 3.212E6, 7.085E-2)
    0.0864419738671184

    References
    ----------
    .. [1] Missenard, F. A., Thermal Conductivity of Organic Liquids of a
       Series or a Group of Liquids , Rev. Gen.Thermodyn., 101 649 (1970).
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    Tr = T/Tc
    Pr = P/Pc
    Tr_2_10 = Tr**(0.2)
    Tr_12_10 = Tr_2_10*Tr_2_10*Tr_2_10
    Tr_12_10 *= Tr_12_10
    return kl*(0.98 + Tr_12_10*(0.0079*Pr*Tr_2_10 + 0.63*(Pr/(30. + Pr))))


Trs_Missenard = [0.8, 0.7, 0.6, 0.5]
Prs_Missenard = [1, 5, 10, 50, 100, 200]
Qs_Missenard = [[0.036, 0.038, 0.038, 0.038, 0.038, 0.038],
                [0.018, 0.025, 0.027, 0.031, 0.032, 0.032],
                [0.015, 0.020, 0.022, 0.024, 0.025, 0.025],
                [0.012, 0.0165, 0.017, 0.019, 0.020, 0.020]]
# tck obtained with interp1d's regrid_smth
Missenard_tck = implementation_optimize_tck([[1.0, 1.0, 5.0, 10.0, 50.0, 100.0, 200.0, 200.0],
                                             [0.5, 0.5, 0.6, 0.7, 0.8, 0.8],
                                             [0.012, 0.015, 0.018, 0.036, 0.0165, 0.02,
                                              0.025, 0.038, 0.017, 0.022, 0.027, 0.038,
                                              0.019, 0.024, 0.031, 0.038, 0.02, 0.025,
                                              0.032, 0.038, 0.02, 0.025, 0.032, 0.038],
                                              1,1])

def Missenard(T, P, Tc, Pc, kl):
    r'''Adjustes for pressure the thermal conductivity of a liquid using an
    emperical formula based on [1]_, but as given in [2]_.

    .. math::
        \frac{k}{k^*} = 1 + Q P_r^{0.7}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    P : float
        Pressure of fluid [Pa]
    Tc: float
        Critical point of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    kl : float
        Thermal conductivity of liquid at 1 atm or saturation, [W/m/K]

    Returns
    -------
    kl_dense : float
        Thermal conductivity of liquid at P, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    An interpolation routine is used here from tabulated values of Q.
    The original source has not been reviewed.

    Examples
    --------
    Example from [2]_, toluene; matches.

    >>> Missenard(304., 6330E5, 591.8, 41E5, 0.129)
    0.2198375777069657

    References
    ----------
    .. [1] Missenard, F. A., Thermal Conductivity of Organic Liquids of a
       Series or a Group of Liquids , Rev. Gen.Thermodyn., 101 649 (1970).
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    Pr = P/Pc
    Q = float(bisplev(Pr, Tr, Missenard_tck))
    return kl*(1. + Q*Pr**0.7)

### Thermal conductivity of liquid mixtures


def DIPPR9H(ws, ks):
    r'''Calculates thermal conductivity of a liquid mixture according to
    mixing rules in [1]_ and also in [2]_.

    .. math::
        \lambda_m = \left( \sum_i w_i \lambda_i^{-2}\right)^{-1/2}

    This is also called the Vredeveld (1973) equation. A review in [3]_ finds
    this the best model on average. However, they did caution that in some
    cases a linear mole-fraction mixing rule performs better. This equation
    according to Poling [1]_ should not be used if some components have
    thermal conductivities more than twice other components. They also say this
    should not be used with water.

    Parameters
    ----------
    ws : float
        Mass fractions of components
    ks : float
        Liquid thermal conductivites of all components, [W/m/K]

    Returns
    -------
    kl : float
        Thermal conductivity of liquid mixture, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The example is from [2]_; all results agree.
    The original source has not been reviewed.

    DIPPR Procedure 9H: Method for the Thermal Conductivity of Nonaqueous Liquid Mixtures

    Average deviations of 3%. for 118 nonaqueous systems with 817 data points.
    Max deviation 20%. According to DIPPR.

    In some sources, this equation is given with the molecular weights included:

    .. math::
        \lambda_m^{-2} = \frac{\sum_i z_i {MW}_i \lambda_i^{-2}}
        {\sum_i z_i {MW}_i}

    Examples
    --------
    >>> DIPPR9H([0.258, 0.742], [0.1692, 0.1528])
    0.15657104706719646

    References
    ----------
    .. [1] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E. The
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    .. [3] Focke, Walter W. "Correlating Thermal-Conductivity Data for Ternary
       Liquid Mixtures." International Journal of Thermophysics 29, no. 4
       (August 1, 2008): 1342-60. https://doi.org/10.1007/s10765-008-0465-2.
    '''
    kl = 0.0
    for i in range(len(ws)):
        kl += ws[i]/(ks[i]*ks[i])
    return 1.0/sqrt(kl)

def DIPPR9I(zs, Vms, ks):
    r'''Calculates thermal conductivity of a liquid mixture according to
    mixing rules in [1]_. This is recommended in [2]_ for aqueous and
    nonaqueous systems.

    .. math::
        k_{mix} = \sum_{i}\sum_j \phi_i\phi_j k_{i,j}

    .. math::
        k_{i,j} = \frac{2}{\frac{1}{k_i} + \frac{1}{k_j}}

    .. math::
        \phi_i = \frac{z_i V_{m,i}}{\sum_j^n z_j V_{m,j}}

    Parameters
    ----------
    zs : list[float]
        Mole fractions of components, [-]
    Vms : list[float]
        Molar volumes of each component, [m^3/mol]
    ks : float
        Liquid thermal conductivites of all components, [W/m/K]

    Returns
    -------
    kl : float
        Thermal conductivity of liquid mixture, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The example is from [2]_; all results agree.

    [2]_ found average deviations of 4-6% for 118 nonaqueous systems
    and 15 aqueous systems at atmospheric pressure, with a maximum deviation of
    33%.

    The computational complexity here is N^2, with a division present in the
    inner loop.

    Examples
    --------
    >>> DIPPR9I(zs=[.682, .318], Vms=[1.723e-2, 7.338e-2], ks=[.6037, .1628])
    0.25397430656658937

    References
    ----------
    .. [1] Li, C. C. "Thermal Conductivity of Liquid Mixtures." AIChE Journal
       22, no. 5 (1976): 927–30. https://doi.org/10.1002/aic.690220520.
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    N = len(zs)
    k = 0.0
    # Precomputation
    ks_inv = [0.0]*N
    phis = [0.0]*N
    tot = 0.0
    for i in range(N):
        val = zs[i]*Vms[i]
        phis[i] = val
        tot += val
    tot = 1.0/tot
    for i in range(N):
        phis[i] *= tot

    # Compute the diagonal and store ks_inv
    for i in range(N):
        k_inv = 1.0/ks[i]
        k += phis[i]*phis[i]*ks[i]
        ks_inv[i] = k_inv

    # Main loop
    main_k_sum = 0.0
    for i in range(N):
        tot = 0.0
        for j in range(i):
            tot += phis[j]/(ks_inv[i] + ks_inv[j])
        main_k_sum += tot*phis[i]

    # factored out 4 - 2 from inner loop, two from symmetry
    k += 4.0*main_k_sum
    return k

def Filippov(ws, ks):
    r'''Calculates thermal conductivity of a binary liquid mixture according to
    mixing rules in [2]_ as found in [1]_.

    .. math::
        \lambda_m = w_1 \lambda_1 + w_2\lambda_2
        - 0.72 w_1 w_2(\lambda_2-\lambda_1)

    Parameters
    ----------
    ws : float
        Mass fractions of components
    ks : float
        Liquid thermal conductivites of all components, [W/m/K]

    Returns
    -------
    kl : float
        Thermal conductivity of liquid mixture, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The original source has not been reviewed.
    Only useful for binary mixtures.

    Examples
    --------
    >>> Filippov([0.258, 0.742], [0.1692, 0.1528])
    0.15929167628799998

    References
    ----------
    .. [1] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E. The
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    .. [2] Filippov, L. P.: Vest. Mosk. Univ., Ser. Fiz. Mat. Estestv. Nauk,
       (8I0E): 67-69A955); Chem. Abstr., 50: 8276 A956).
       Filippov, L. P., and N. S. Novoselova: Vestn. Mosk. Univ., Ser. F
       iz. Mat. Estestv.Nauk, CI0B): 37-40A955); Chem. Abstr., 49: 11366 A955).
    '''
    len_ks = len(ks)
    if len_ks != len(ws) or len_ks != 2:
        raise ValueError("Filippov method is only defined for mixtures of two components")
    return ws[0]*ks[0] + ws[1]*ks[1] - 0.72*ws[0]*ws[1]*(ks[1] - ks[0])



### Thermal Conductivity of Gases

def Eucken(MW, Cvm, mu):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the CSP method of Eucken [1]_.

    .. math::
        \frac{\lambda MW}{\eta C_v} = 1 + \frac{9/4}{C_v/R}

    Parameters
    ----------
    MW : float
        Molecular weight of the gas [g/mol]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]
    mu : float
        Gas viscosity [Pa*s]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    Temperature dependence is introduced via heat capacity and viscosity.
    A theoretical equation. No original author located.
    MW internally converted to kg/g-mol.

    Examples
    --------
    2-methylbutane at low pressure, 373.15 K. Mathes calculation in [1]_.

    >>> Eucken(MW=72.151, Cvm=135.9, mu=8.77E-6)
    0.018792645058456698

    References
    ----------
    .. [1] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    MW = MW*1e-3
    return (1. + 9.0/4.0*R/Cvm)*mu*Cvm/MW


def Eucken_modified(MW, Cvm, mu):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the Modified CSP method of Eucken [1]_.

    .. math::
        \frac{\lambda MW}{\eta C_v} = 1.32 + \frac{1.77}{C_v/R}

    Parameters
    ----------
    MW : float
        Molecular weight of the gas [g/mol]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]
    mu : float
        Gas viscosity [Pa*s]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    Temperature dependence is introduced via heat capacity and viscosity.
    A theoretical equation. No original author located.
    MW internally converted to kg/g-mol.

    Examples
    --------
    2-methylbutane at low pressure, 373.15 K. Mathes calculation in [1]_.

    >>> Eucken_modified(MW=72.151, Cvm=135.9, mu=8.77E-6)
    0.02359353760551249

    References
    ----------
    .. [1] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    MW = MW*1e-3
    return (1.32 + 1.77*R/Cvm)*mu*Cvm/MW


def DIPPR9B(T, MW, Cvm, mu, Tc=None, chemtype=None):
    r'''Calculates the thermal conductivity of a gas using one of several
    emperical equations developed in [1]_, [2]_, and presented in [3]_.

    For monoatomic gases:

    .. math::
        k = 2.5 \frac{\eta C_v}{MW}

    For linear molecules:

    .. math::
        k = \frac{\eta}{MW} \left( 1.30 C_v + 14644.00 - \frac{2928.80}{T_r}\right)

    For nonlinear molecules:

    .. math::
        k = \frac{\eta}{MW}(1.15C_v + 16903.36)

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    MW : float
        Molwcular weight of fluid [g/mol]
    Cvm : float
        Molar heat capacity at constant volume of fluid, [J/mol/K]
    mu : float
        Viscosity of gas, [Pa*s]

    Returns
    -------
    k_g : float
        Thermal conductivity of gas, [W/m/k]

    Notes
    -----
    Tested with DIPPR values.
    Cvm is internally converted to J/kmol/K.

    Examples
    --------
    CO:

    >>> DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear')
    0.01813208676438415

    References
    ----------
    .. [1] Bromley, LeRoy A., Berkeley. University of California, and U.S.
       Atomic Energy Commission. Thermal Conductivity of Gases at Moderate
       Pressures. UCRL;1852. Berkeley, CA: University of California Radiation
       Laboratory, 1952.
    .. [2] Stiel, Leonard I., and George Thodos. "The Thermal Conductivity of
       Nonpolar Substances in the Dense Gaseous and Liquid Regions." AIChE
       Journal 10, no. 1 (January 1, 1964): 26-30. doi:10.1002/aic.690100114
    .. [3] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    Cvm = Cvm*1000.  # J/g/K to J/kmol/K
    if chemtype == 'monoatomic':
        return 2.5*mu*Cvm/MW
    elif chemtype == 'linear' or chemtype is None:
        Tr = T/Tc
        return mu/MW*(1.30*Cvm + 14644 - 2928.80/Tr)
    elif chemtype == 'nonlinear':
        return mu/MW*(1.15*Cvm + 16903.36)
    else:
        raise ValueError('Specified chemical type is not an option')


def Chung(T, MW, Tc, omega, Cvm, mu):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the CSP method of Chung [1]_.

    .. math::
        \frac{\lambda MW}{\eta C_v} = \frac{3.75 \Psi}{C_v/R}

    .. math::
        \Psi = 1 + \alpha \left\{[0.215+0.28288\alpha-1.061\beta+0.26665Z]/
        [0.6366+\beta Z + 1.061 \alpha \beta]\right\}

    .. math::
        \alpha = \frac{C_v}{R}-1.5

    .. math::
        \beta = 0.7862-0.7109\omega + 1.3168\omega^2

    .. math::
        Z=2+10.5T_r^2

    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]
    Tc : float
        Critical temperature of the gas [K]
    omega : float
        Acentric factor of the gas [-]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]
    mu : float
        Gas viscosity [Pa*s]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    MW internally converted to kg/g-mol.

    Examples
    --------
    2-methylbutane at low pressure, 373.15 K. Mathes calculation in [2]_.

    >>> Chung(T=373.15, MW=72.151, Tc=460.4, omega=0.227, Cvm=135.9, mu=8.77E-6)
    0.023015653797111124

    References
    ----------
    .. [1] Chung, Ting Horng, Lloyd L. Lee, and Kenneth E. Starling.
       "Applications of Kinetic Gas Theories and Multiparameter Correlation for
       Prediction of Dilute Gas Viscosity and Thermal Conductivity."
       Industrial & Engineering Chemistry Fundamentals 23, no. 1
       (February 1, 1984): 8-13. doi:10.1021/i100013a002
    .. [2] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    MW = MW*1e-3
    alpha = Cvm*R_inv - 1.5
    beta = 0.7862 - 0.7109*omega + 1.3168*omega*omega
    Z = 2.0 + 10.5*T*T/(Tc*Tc)
    psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)
                      /(0.6366 + beta*Z + 1.061*alpha*beta))
    return 3.75*psi/(Cvm*MW)*R*mu*Cvm


def Eli_Hanley(T, MW, Tc, Vc, Zc, omega, Cvm):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the reference fluid method of Eli and Hanley [1]_ as
    shown in [2]_.

    .. math::
        \lambda = \lambda^* + \frac{\eta^*}{MW}(1.32)\left(C_v - \frac{3R}{2}\right)

    .. math::
        Tr = \text{min}(Tr, 2)

    .. math::
        \theta = 1 + (\omega-0.011)\left(0.56553 - 0.86276\ln Tr - \frac{0.69852}{Tr}\right)

    .. math::
        \psi = [1 + (\omega - 0.011)(0.38560 - 1.1617\ln Tr)]\frac{0.288}{Z_c}

    .. math::
        f = \frac{T_c}{190.4}\theta

    .. math::
        h = \frac{V_c}{9.92E-5}\psi

    .. math::
        T_0 = T/f

    .. math::
        \eta_0^*(T_0)= \sum_{n=1}^9 C_n T_0^{(n-4)/3}

    .. math::
        \theta_0 = 1944 \eta_0

    .. math::
        \lambda^* = \lambda_0 H

    .. math::
        \eta^* = \eta^*_0 H \frac{MW}{16.04}

    .. math::
        H = \left(\frac{16.04}{MW}\right)^{0.5}f^{0.5}/h^{2/3}

    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]
    Tc : float
        Critical temperature of the gas [K]
    Vc : float
        Critical volume of the gas [m^3/mol]
    Zc : float
        Critical compressibility of the gas []
    omega : float
        Acentric factor of the gas [-]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    Reference fluid is Methane.
    MW internally converted to kg/g-mol.

    Examples
    --------
    2-methylbutane at low pressure, 373.15 K. Matches calculation in [2]_.

    >>> Eli_Hanley(T=373.15, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267,
    ... omega=0.227, Cvm=135.9)
    0.02247951724513664

    References
    ----------
    .. [1] Ely, James F., and H. J. M. Hanley. "Prediction of Transport
       Properties. 2. Thermal Conductivity of Pure Fluids and Mixtures."
       Industrial & Engineering Chemistry Fundamentals 22, no. 1 (February 1,
       1983): 90-97. doi:10.1021/i100009a016.
    .. [2] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''

    Tr = T/Tc
    if Tr > 2.0:
        Tr = 2.0
    logTr = log(Tr)
    theta = 1.0 + (omega - 0.011)*(0.56553 - 0.86276*logTr - 0.69852/Tr)
    psi = (1.0 + (omega - 0.011)*(0.38560 - 1.1617*logTr))*0.288/Zc
    f = Tc/190.4*theta
    h = Vc/9.92E-5*psi
    T0 = T/f

    T0_third = T0**(1.0/3.0)
    T0_moving = 1.0/T0
    tot = (2907741.307*T0_moving + T0_third*(-3312874.033*T0_moving
            + T0_third*(1608101.838*T0_moving + T0_third*(-433190.4871*T0_moving
            + T0_third*(70624.8133*T0_moving + T0_third*(-7116.62075*T0_moving
         + T0_third*(432.51744*T0_moving + T0_third*(0.2037119479*T0_moving*T0_third
        - 14.4591121*T0_moving))))))))

#    Cs = [2.907741307E6, -3.312874033E6, 1.608101838E6, -4.331904871E5,
#          7.062481330E4, -7.116620750E3, 4.325174400E2, -1.445911210E1, 2.037119479E-1]
#    tot = 0.0
#    for i in range(9):
#        tot += Cs[i]*T0_moving
#        T0_moving *= T0_third
    eta0 = 1e-7*tot
    k0 = 1944.0*eta0

    H = sqrt(f*16.04/MW)*h**(-2.0/3.)
    etas = eta0*H*MW*0.06234413965087282 # /16.04
    ks = k0*H
    return ks + 1320.0*etas/MW*(Cvm - 1.5*R)


def Gharagheizi_gas(T, MW, Tb, Pc, omega):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the CSP method of Gharagheizi [1]_. A  convoluted
    method claiming high-accuracy and using only statistically significant
    variable following analalysis.

    Requires temperature, molecular weight, boiling temperature and critical
    pressure and acentric factor.

    .. math::
        k = 7.9505\times 10^{-4} + 3.989\times 10^{-5} T
        -5.419\times 10^-5 MW + 3.989\times 10^{-5} A

    .. math::
       A = \frac{\left(2\omega + T - \frac{(2\omega + 3.2825)T}{T_b} + 3.2825\right)}{0.1MP_cT}
        \times (3.9752\omega + 0.1 P_c + 1.9876B + 6.5243)^2


    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    MW: float
        Molecular weight of the fluid [g/mol]
    Tb : float
        Boiling temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor of the fluid [-]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    Pressure is internally converted into 10*kPa but author used correlation with
    kPa; overall, errors have been corrected in the presentation of the formula.

    This equation was derived with 15927 points and 1574 compounds.
    Example value from [1]_ is the first point in the supportinf info, for CH4.

    Examples
    --------
    >>> Gharagheizi_gas(580., 16.04246, 111.66, 4599000.0, 0.0115478000)
    0.09594861261873211

    References
    ----------
    .. [1] Gharagheizi, Farhad, Poorandokht Ilani-Kashkouli, Mehdi Sattari,
       Amir H. Mohammadi, Deresh Ramjugernath, and Dominique Richon.
       "Development of a General Model for Determination of Thermal
       Conductivity of Liquid Chemical Compounds at Atmospheric Pressure."
       AIChE Journal 59, no. 5 (May 1, 2013): 1702-8. doi:10.1002/aic.13938
    '''
    Pc = Pc*1e-4
    Tb_inv = 1.0/Tb
    B = (T + (2.*omega + 2.*T - 2.*T*(2.*omega + 3.2825)*Tb_inv + 3.2825)
         /(2.0*omega + T - T*(2.0*omega + 3.2825)*Tb_inv + 3.2825)
         - T*(2.0*omega + 3.2825)*Tb_inv)

    x0 = (3.9752*omega + 0.1*Pc + 1.9876*B + 6.5243)
    A = (2.0*omega + T - T*(2.0*omega + 3.2825)*Tb_inv + 3.2825)/(0.1*MW*Pc*T) * x0*x0
    return 7.9505E-4 + 3.989E-5*T - 5.419E-5*MW + 3.989E-5*A


def Bahadori_gas(T, MW):
    r'''Estimates the thermal conductivity of hydrocarbons gases at low P.
    Fits their data well, and is useful as only MW is required.
    Y is the Molecular weight, and X the temperature.

    .. math::
        K = a + bY + CY^2 + dY^3

    .. math::
        a = A_1 + B_1 X + C_1 X^2 + D_1 X^3

    .. math::
        b = A_2 + B_2 X + C_2 X^2 + D_2 X^3

    .. math::
        c = A_3 + B_3 X + C_3 X^2 + D_3 X^3

    .. math::
        d = A_4 + B_4 X + C_4 X^2 + D_4 X^3

    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]

    Returns
    -------
    kg : float
        Estimated gas thermal conductivity [W/m/k]

    Notes
    -----
    The accuracy of this equation has not been reviewed.

    Examples
    --------
    >>> Bahadori_gas(40+273.15, 20.0) # Point from article
    0.03196816533787329

    References
    ----------
    .. [1] Bahadori, Alireza, and Saeid Mokhatab. "Estimating Thermal
       Conductivity of Hydrocarbons." Chemical Engineering 115, no. 13
       (December 2008): 52-54
    '''
    A = (4.3931323468E-1, -3.88001122207E-2, 9.28616040136E-4, -6.57828995724E-6)
    B = (-2.9624238519E-3, 2.67956145820E-4, -6.40171884139E-6, 4.48579040207E-8)
    C = (7.54249790107E-6, -6.46636219509E-7, 1.5124510261E-8, -1.0376480449E-10)
    D = (-6.0988433456E-9, 5.20752132076E-10, -1.19425545729E-11, 8.0136464085E-14)
    X, Y = T, MW
    a = A[0] + X*(B[0] + X*(C[0] + D[0]*X))
    b = A[1] + X*(B[1] + X*(C[1] + D[1]*X))
    c = A[2] + X*(B[2] + X*(C[2] + D[2]*X))
    d = A[3] + X*(B[3] + X*(C[3] + D[3]*X))
    return a + Y*(b + Y*(c + d*Y))



### Thermal Conductivity of dense gases

def Stiel_Thodos_dense(T, MW, Tc, Pc, Vc, Zc, Vm, kg):
    r'''Estimates the thermal conductivity of a gas at high pressure as a
    function of temperature using difference method of Stiel and Thodos [1]_
    as shown in [2]_.

    if :math:`\rho_r < 0.5`:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

    if :math:`0.5 < \rho_r < 2.0`:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

    if :math:`2 < \rho_r < 2.8`:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

    .. math::
        \Gamma = 210 \left(\frac{T_cMW^3}{P_c^4}\right)^{1/6}

    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]
    Tc : float
        Critical temperature of the gas [K]
    Pc : float
        Critical pressure of the gas [Pa]
    Vc : float
        Critical volume of the gas [m^3/mol]
    Zc : float
        Critical compressibility of the gas [-]
    Vm : float
        Molar volume of the gas at T and P [m^3/mol]
    kg : float
        Low-pressure gas thermal conductivity [W/m/k]

    Returns
    -------
    kg : float
        Estimated dense gas thermal conductivity [W/m/k]

    Notes
    -----
    Pc is internally converted to bar.

    Examples
    --------
    >>> Stiel_Thodos_dense(T=378.15, MW=44.013, Tc=309.6, Pc=72.4E5,
    ... Vc=97.4E-6, Zc=0.274, Vm=144E-6, kg=2.34E-2)
    0.041245574404863684

    References
    ----------
    .. [1] Stiel, Leonard I., and George Thodos. "The Thermal Conductivity of
       Nonpolar Substances in the Dense Gaseous and Liquid Regions." AIChE
       Journal 10, no. 1 (January 1, 1964): 26-30. doi:10.1002/aic.690100114.
    .. [2] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    gamma = 210.0*(Tc*MW*MW*MW*(Pc*1e-5)**-4.0)**(1.0/6.0)
    rhor = Vc/Vm
    if rhor < 0.5:
        term = 1.22E-2*(exp(0.535*rhor) - 1.)
    elif rhor < 2.0:
        term = 1.14E-2*(exp(0.67*rhor) - 1.069)
    else:
        # Technically only up to 2.8
        term = 2.60E-3*(exp(1.155*rhor) + 2.016)
    diff = term*Zc**-5/gamma
    kg = kg + diff
    return kg


def Eli_Hanley_dense(T, MW, Tc, Vc, Zc, omega, Cvm, Vm):
    r'''Estimates the thermal conductivity of a gas at high pressure as a
    function of temperature using the reference fluid method of Eli and
    Hanley [1]_ as shown in [2]_.

    .. math::
        Tr = min(Tr, 2)

    .. math::
        Vr = min(Vr, 2)

    .. math::
        f = \frac{T_c}{190.4}\theta

    .. math::
        h = \frac{V_c}{9.92E-5}\psi

    .. math::
        T_0 = T/f

    .. math::
        \rho_0 = \frac{16.04}{V}h

    .. math::
        \theta = 1 + (\omega-0.011)\left(0.09057 - 0.86276\ln Tr + \left(
        0.31664 - \frac{0.46568}{Tr}\right) (V_r - 0.5)\right)

    .. math::
        \psi = [1 + (\omega - 0.011)(0.39490(V_r - 1.02355) - 0.93281(V_r -
        0.75464)\ln T_r]\frac{0.288}{Z_c}

    .. math::
        \lambda_1 = 1944 \eta_0

    .. math::
        \lambda_2 = \left\{b_1 + b_2\left[b_3 - \ln \left(\frac{T_0}{b_4}
        \right)\right]^2\right\}\rho_0

    .. math::
        \lambda_3 = \exp\left(a_1 + \frac{a_2}{T_0}\right)\left\{\exp[(a_3 +
        \frac{a_4}{T_0^{1.5}})\rho_0^{0.1} + (\frac{\rho_0}{0.1617} - 1)
        \rho_0^{0.5}(a_5 + \frac{a_6}{T_0} + \frac{a_7}{T_0^2})] - 1\right\}

    .. math::
        \lambda^{**} = [\lambda_1 + \lambda_2 + \lambda_3]H

    .. math::
        H = \left(\frac{16.04}{MW}\right)^{0.5}f^{0.5}/h^{2/3}

    .. math::
        X = \left\{\left[1 - \frac{T}{f}\left(\frac{df}{dT}\right)_v \right]
        \frac{0.288}{Z_c}\right\}^{1.5}

    .. math::
        \left(\frac{df}{dT}\right)_v = \frac{T_c}{190.4}\left(\frac{d\theta}
        {d T}\right)_v

    .. math::
        \left(\frac{d\theta}{d T}\right)_v = (\omega-0.011)\left[
        \frac{-0.86276}{T} + (V_r-0.5)\frac{0.46568T_c}{T^2}\right]

    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]
    Tc : float
        Critical temperature of the gas [K]
    Vc : float
        Critical volume of the gas [m^3/mol]
    Zc : float
        Critical compressibility of the gas []
    omega : float
        Acentric factor of the gas [-]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]
    Vm : float
        Volume of the gas at T and P [m^3/mol]

    Returns
    -------
    kg : float
        Estimated dense gas thermal conductivity [W/m/k]

    Notes
    -----
    Reference fluid is Methane.
    MW internally converted to kg/g-mol.

    Examples
    --------
    >>> Eli_Hanley_dense(T=473., MW=42.081, Tc=364.9, Vc=1.81E-4, Zc=0.274,
    ... omega=0.144, Cvm=82.70, Vm=1.721E-4)
    0.06038475754109959

    References
    ----------
    .. [1] Ely, James F., and H. J. M. Hanley. "Prediction of Transport
       Properties. 2. Thermal Conductivity of Pure Fluids and Mixtures."
       Industrial & Engineering Chemistry Fundamentals 22, no. 1 (February 1,
       1983): 90-97. doi:10.1021/i100009a016.
    .. [2] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Cs = [2.907741307E6, -3.312874033E6, 1.608101838E6, -4.331904871E5,
          7.062481330E4, -7.116620750E3, 4.325174400E2, -1.445911210E1,
          2.037119479E-1]

    Tr = T/Tc
    if Tr > 2.0:
        Tr = 2.0
    Vr = Vm/Vc
    if Vr > 2.0:
        Vr = 2.0
    logTr = log(Tr)
    theta = 1.0 + (omega - 0.011)*(0.09057 - 0.86276*logTr + (0.31664 - 0.46568/Tr)*(Vr-0.5))
    psi = (1.0 + (omega- 0.011)*(0.39490*(Vr-1.02355) - 0.93281*(Vr-0.75464)*logTr))*0.288/Zc
    f = Tc/190.4*theta
    h = Vc/9.92E-5*psi
    T0 = T/f
    rho0 = 16.04/(Vm*1E6)*h  # Vm must be in cm^3/mol here.

    T0_third = T0**(1.0/3.0)
    T0_moving = 1.0/T0
    tot = 0.0
    for i in range(9):
        tot += Cs[i]*T0_moving
        T0_moving *= T0_third

    eta0 = 1E-7*tot
    k1 = 1944*eta0
    b1 = -0.25276920E0
    b2 = 0.334328590E0
    b3 = 1.12
    b4 = 0.1680E3
    x0 = (b3 - log(T0/b4))
    k2 = (b1 + b2*x0*x0)*1e-3*rho0

    a1 = -7.19771
    a2 = 85.67822
    a3 = 12.47183
    a4 = -984.6252
    a5 = 0.3594685
    a6 = 69.79841
    a7 = -872.8833

    k3 = exp(a1 + a2/T0)*(exp((a3 + a4/T0**1.5)*rho0**0.1 + (rho0/0.1617 - 1.0)*rho0**0.5*(a5 + a6/T0 + a7/T0**2)) - 1)/1000.

    if T/Tc > 2.0:
        dtheta = 0.0
    else:
        dtheta = (omega - 0.011)*(-0.86276/T + (Vr-0.5)*0.46568*Tc/T**2)
    dfdT = Tc/190.4*dtheta
    X = ((1.0 - T/f*dfdT)*0.288/Zc)**1.5

    H = (16.04/MW)**0.5*f**0.5/h**(2/3.)
    ks = (k1*X + k2 + k3)*H

    ### Uses calculations similar to those for pure species here
    theta = 1.0 + (omega - 0.011)*(0.56553 - 0.86276*logTr - 0.69852/Tr)
    psi = (1.0 + (omega - 0.011)*(0.38560 - 1.1617*logTr))*0.288/Zc
    f = Tc/190.4*theta
    h = Vc/9.92E-5*psi
    T0 = T/f

    T0_third = T0**(1.0/3.0)
    T0_moving = 1.0/T0
    tot = 0.0
    for i in range(9):
        tot += Cs[i]*T0_moving
        T0_moving *= T0_third

    eta0 = 1E-7*tot
    H = (16.04*f/MW)**0.5*h**(-2.0/3.)
    etas = eta0*H*MW/16.04
    k = ks + etas*1e3/(MW)*1.32*(Cvm - 1.5*R)
    return k


def Chung_dense(T, MW, Tc, Vc, omega, Cvm, Vm, mu, dipole, association=0.0):
    r'''Estimates the thermal conductivity of a gas at high pressure as a
    function of temperature using the reference fluid method of
    Chung [1]_ as shown in [2]_.

    .. math::
        \lambda = \frac{31.2 \eta^\circ \Psi}{M'}(G_2^{-1} + B_6 y)+qB_7y^2T_r^{1/2}G_2

    .. math::
        \Psi = 1 + \alpha \left\{[0.215+0.28288\alpha-1.061\beta+0.26665Z]/
        [0.6366+\beta Z + 1.061 \alpha \beta]\right\}

    .. math::
        \alpha = \frac{C_v}{R}-1.5

    .. math::
        \beta = 0.7862-0.7109\omega + 1.3168\omega^2

    .. math::
        Z=2+10.5T_r^2

    .. math::
        q = 3.586\times 10^{-3} (T_c/M')^{1/2}/V_c^{2/3}

    .. math::
        y = \frac{V_c}{6V}

    .. math::
        G_1 = \frac{1-0.5y}{(1-y)^3}

    .. math::
        G_2 = \frac{(B_1/y)[1-\exp(-B_4y)]+ B_2G_1\exp(B_5y) + B_3G_1}
        {B_1B_4 + B_2 + B_3}

    .. math::
        B_i = a_i + b_i \omega + c_i \mu_r^4 + d_i \kappa


    Parameters
    ----------
    T : float
        Temperature of the gas [K]
    MW : float
        Molecular weight of the gas [g/mol]
    Tc : float
        Critical temperature of the gas [K]
    Vc : float
        Critical volume of the gas [m^3/mol]
    omega : float
        Acentric factor of the gas [-]
    Cvm : float
        Molar contant volume heat capacity of the gas [J/mol/K]
    Vm : float
        Molar volume of the gas at T and P [m^3/mol]
    mu : float
        Low-pressure gas viscosity [Pa*s]
    dipole : float
        Dipole moment [debye]
    association : float, optional
        Association factor [-]

    Returns
    -------
    kg : float
        Estimated dense gas thermal conductivity [W/m/k]

    Notes
    -----
    MW internally converted to kg/g-mol.
    Vm internally converted to mL/mol.
    [1]_ is not the latest form as presented in [1]_.
    Association factor is assumed 0. Relates to the polarity of the gas.

    Coefficients as follows:

    ais = [2.4166E+0, -5.0924E-1, 6.6107E+0, 1.4543E+1, 7.9274E-1, -5.8634E+0, 9.1089E+1]

    bis = [7.4824E-1, -1.5094E+0, 5.6207E+0, -8.9139E+0, 8.2019E-1, 1.2801E+1, 1.2811E+2]

    cis = [-9.1858E-1, -4.9991E+1, 6.4760E+1, -5.6379E+0, -6.9369E-1, 9.5893E+0, -5.4217E+1]

    dis = [1.2172E+2, 6.9983E+1, 2.7039E+1, 7.4344E+1, 6.3173E+0, 6.5529E+1, 5.2381E+2]


    Examples
    --------
    >>> Chung_dense(T=473., MW=42.081, Tc=364.9, Vc=184.6E-6, omega=0.142,
    ... Cvm=82.67, Vm=172.1E-6, mu=134E-7, dipole=0.4)
    0.06160569232570781

    References
    ----------
    .. [1] Chung, Ting Horng, Mohammad Ajlan, Lloyd L. Lee, and Kenneth E.
       Starling. "Generalized Multiparameter Correlation for Nonpolar and Polar
       Fluid Transport Properties." Industrial & Engineering Chemistry Research
       27, no. 4 (April 1, 1988): 671-79. doi:10.1021/ie00076a024.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    mur = 131.3*dipole*(Vc*1E6*Tc)**-0.5
    mur4 = mur*mur
    mur4 *= mur4

    # From Chung Method
    alpha = Cvm*R_inv - 1.5
    beta = 0.7862 - 0.7109*omega + 1.3168*omega*omega
    Z = 2.0 + 10.5*Tr*Tr
    psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta))

    y = Vc/(6.0*Vm)
    B1 = 2.4166 + 0.74824*omega + -0.91858*mur4 + 121.72*association
    B2 = -0.50924 + -1.5094*omega + -49.991*mur4 + 69.983*association
    B3 = 6.6107 + 5.6207*omega + 64.76*mur4 + 27.039*association
    B4 = 14.543 + -8.9139*omega + -5.6379*mur4 + 74.344*association
    B5 = 0.79274 + 0.82019*omega + -0.69369*mur4 + 6.3173*association
    B6 = -5.8634 + 12.801*omega + 9.5893*mur4 + 65.529*association
    B7 = 91.089 + 128.11*omega + -54.217*mur4 + 523.81*association

    x0 = 1.0 - y
    G1 = (1.0 - 0.5*y)/(x0*x0*x0)
    G2 = (B1/y*(1.0 - exp(-B4*y)) + B2*G1*exp(B5*y) + B3*G1)/(B1*B4 + B2 + B3)
    MW_kg_inv = 1000.0/MW
    q = 3.586E-3*(Tc*MW_kg_inv)**0.5*(Vc*1E6)**(-2.0/3.)
    return 31.2*mu*psi*MW_kg_inv*(1.0/G2 + B6*y) + q*B7*y*y*Tr**0.5*G2


### Thermal conductivity of gas mixtures


def Lindsay_Bromley(T, ys, ks, mus, Tbs, MWs):
    r'''Calculates thermal conductivity of a gas mixture according to
    mixing rules in [1]_ and also in [2]_. It is significantly more complicated
    than other kinetic theory models.

    .. math::
        k = \sum_i \frac{y_i k_i}{\sum_j y_i A_{ij}}

    .. math::
        A_{ij} = \frac{1}{4} \left\{ 1 + \left[\frac{\eta_i}{\eta_j}
        \left(\frac{MW_j}{MW_i}\right)^{0.75} \left( \frac{T+S_i}{T+S_j}\right)
        \right]^{0.5} \right\}^2 \left( \frac{T+S_{ij}}{T+S_i}\right)

    .. math::
        S_{ij} = S_{ji} = (S_i S_j)^{0.5}

    .. math::
        S_i = 1.5 T_b

    Parameters
    ----------
    T : float
        Temperature of gas [K]
    ys : float
        Mole fractions of gas components
    ks : float
        Gas thermal conductivites of all components, [W/m/K]
    mus : float
        Gas viscosities of all components, [Pa*s]
    Tbs : float
        Boiling points of all components, [K]
    MWs : float
        Molecular weights of all components, [g/mol]

    Returns
    -------
    kg : float
        Thermal conductivity of gas mixture, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The example is from [2]_; all results agree.
    The original source has not been reviewed.

    DIPPR Procedure 9D: Method for the Thermal Conductivity of Gas Mixtures

    Average deviations of 4-5% for 77 binary mixtures reviewed in [2]_, from
    1342 points; also six ternary mixtures (70  points); max deviation observed
    was 40%. (DIPPR)

    Examples
    --------
    >>> Lindsay_Bromley(323.15, [0.23, 0.77], [1.939E-2, 1.231E-2], [1.002E-5, 1.015E-5], [248.31, 248.93], [46.07, 50.49])
    0.013902644179693132

    References
    ----------
    .. [1] Lindsay, Alexander L., and LeRoy A. Bromley. "Thermal Conductivity
       of Gas Mixtures." Industrial & Engineering Chemistry 42, no. 8
       (August 1, 1950): 1508-11. doi:10.1021/ie50488a017.
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    .. [3] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    N = len(ys)
    S_roots = [0.0]*N
    bigis = [0.0]*N
    Ss_invT = [0.0]*N
    TSrootSsinv = [0.0]*N
    bigis_inv = [0.0]*N

    for i in range(N):
        Si = 1.5*Tbs[i]
        S_roots[i] = sqrt(Si)
        T_Si = T + Si
        S_inv = 1.0/T_Si
        Ss_invT[i] = T*S_inv
        TSrootSsinv[i] = S_roots[i]*S_inv
        rt05MW = sqrt(MWs[i])
        rt25MW = sqrt(rt05MW)
        bigis[i] = sqrt(T_Si*mus[i]/(rt05MW*rt25MW))# correct and clever - compute MW^0.375
        bigis_inv[i] = 1.0/bigis[i]

    k = 0.0
    for i in range(N):
        den = 0.0
        S_rooti = S_roots[i]
        for j in range(N):
            # 1 multiply, 3 indexes into different arrays
            x0 = Ss_invT[i] + TSrootSsinv[i]*S_roots[j]
            # 2 multiplies, 3 indexes
#            x0 = (T + S_rooti*S_roots[j])*Ss_inv[i]
            big = 1.0 + bigis[i]*bigis_inv[j]
            Aij = big*big*x0
            den += ys[j]*Aij
        k += ys[i]*ks[i]/den
    k *= 4.0 # constant
    return k

    # Original, unoptimized implementation
#    cmps = range(len(ys))
#    Ss = [1.5*Tb for Tb in Tbs]
#    Sij = [[(Si*Sj)**0.5 for Sj in Ss] for Si in Ss]
#
#    Aij = [[0.25*(1. + (mus[i]/mus[j]*(MWs[j]/MWs[i])**0.75
#            *(T+Ss[i])/(T+Ss[j]))**0.5 )**2 *(T+Sij[i][j])/(T+Ss[i])
#            for j in cmps] for i in cmps]
#
#    return sum([ys[i]*ks[i]/sum(ys[j]*Aij[i][j] for j in cmps) for i in cmps])


def Wassiljewa_Herning_Zipperer(zs, ks, MWs, MW_roots=None):
    r'''Calculates thermal conductivity of a gas mixture according to
    the kinetic theory expression of Wassiljewa with the interaction
    term from the Herning-Zipperer expression. This is also used for
    the prediction of gas mixture viscosity.

    .. math::
        k = \sum \frac{y_i k_i}{\sum y_i A_{ij}}

    .. math::
        A_{ij} = \left(\frac{MW_j}{MW_i}\right)^{0.5}

    Parameters
    ----------
    zs : float
        Mole fractions of gas components, [-]
    ks : float
        gas thermal conductivites of all components, [W/m/K]
    MWs : float
        Molecular weights of all components, [g/mol]
    MW_roots : float, optional
        Square roots of molecular weights of all components;
        speeds up the calculation if provided, [g^0.5/mol^0.5]

    Returns
    -------
    kg : float
        Thermal conductivity of gas mixture, [W/m/K]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.

    Examples
    --------
    >>> Wassiljewa_Herning_Zipperer(zs=[.1, .4, .5], ks=[1.002E-5, 1.15E-5, 2e-5], MWs=[40.0, 50.0, 60.0])
    1.5861181979916883e-05

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return Herning_Zipperer(zs, ks, MWs, MW_roots)
