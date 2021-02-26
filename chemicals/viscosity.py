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

This module contains various viscosity estimation routines, dataframes
of fit coefficients, and mixing rules.

For reporting bugs, adding feature requests, or submitting pull requests,
please use the `GitHub issue tracker <https://github.com/CalebBell/chemicals/>`_.

.. contents:: :local:

Pure Low Pressure Liquid Correlations
-------------------------------------
.. autofunction:: chemicals.viscosity.Letsou_Stiel
.. autofunction:: chemicals.viscosity.Przedziecki_Sridhar

Pure High Pressure Liquid Correlations
--------------------------------------
.. autofunction:: chemicals.viscosity.Lucas

Liquid Mixing Rules
-------------------
No specific correlations are implemented but
:obj:`chemicals.utils.mixing_logarithmic` with weight fractions is the
recommended form.

Pure Low Pressure Gas Correlations
----------------------------------
.. autofunction:: chemicals.viscosity.Yoon_Thodos
.. autofunction:: chemicals.viscosity.Stiel_Thodos
.. autofunction:: chemicals.viscosity.Lucas_gas
.. autofunction:: chemicals.viscosity.viscosity_gas_Gharagheizi

Pure High Pressure Gas Correlations
-----------------------------------
No correlations are implemented yet.

Gas Mixing Rules
----------------
.. autofunction:: chemicals.viscosity.Herning_Zipperer
.. autofunction:: chemicals.viscosity.Brokaw
.. autofunction:: chemicals.viscosity.Wilke
.. autofunction:: chemicals.viscosity.Wilke_prefactors
.. autofunction:: chemicals.viscosity.Wilke_prefactored
.. autofunction:: chemicals.viscosity.Wilke_large

Correlations for Specific Substances
------------------------------------
.. autofunction:: chemicals.viscosity.mu_IAPWS
.. autofunction:: chemicals.viscosity.mu_air_lemmon

Petroleum Correlations
----------------------
.. autofunction:: chemicals.viscosity.Twu_1985
.. autofunction:: chemicals.viscosity.Lorentz_Bray_Clarke

Fit Correlations
----------------
.. autofunction:: chemicals.viscosity.PPDS9
.. autofunction:: chemicals.viscosity.dPPDS9_dT
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_2
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_2_exponential
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_3


Conversion functions
--------------------
.. autofunction:: chemicals.viscosity.viscosity_converter
.. autofunction:: chemicals.viscosity.viscosity_index

Fit Coefficients
----------------
All of these coefficients are lazy-loaded, so they must be accessed as an
attribute of this module.

.. data:: mu_data_Dutt_Prasad

    Coefficient sfor :obj:`chemicals.viscosity.Viswanath_Natarajan_3` from [1]_
    for 100 fluids.

.. data:: mu_data_VN3

    Coefficients for :obj:`chemicals.viscosity.Viswanath_Natarajan_3` from [1]_
    with data for 432 fluids.

.. data:: mu_data_VN2

    Coefficients for :obj:`chemicals.viscosity.Viswanath_Natarajan_2` from [1]_
    with data for 135 fluids.

.. data:: mu_data_VN2E

    Coefficients for :obj:`chemicals.viscosity.Viswanath_Natarajan_2_exponential`
    from [1]_ with data for 14 fluids.

.. data:: mu_data_Perrys_8E_2_313

    A collection of 337 coefficient sets for :obj:`chemicals.dippr.EQ101` from the
    DIPPR database published openly in [3]_.

.. data:: mu_data_Perrys_8E_2_312

    A collection of 345 coefficient sets for :obj:`chemicals.dippr.EQ102` from the
    DIPPR database published openly in [3]_.

.. data:: mu_data_VDI_PPDS_7

    Coefficients for the model equation :obj:`PPDS9`, published openly in [2]_.
    Provides no temperature limits, but has been designed
    for extrapolation. Extrapolated to low temperatures it provides a
    smooth exponential increase. However, for some chemicals such as
    glycerol, extrapolated to higher temperatures viscosity is predicted
    to increase above a certain point.

.. data:: mu_data_VDI_PPDS_8

    Coefficients for a tempereture polynomial (T in Kelvin) developed by the
    PPDS, published openly in [2]_. :math:`\mu = A + BT + CT^2 + DT^3 + ET^4`.

.. [1] Viswanath, Dabir S., and G. Natarajan. Databook On The Viscosity Of
   Liquids. New York: Taylor & Francis, 1989
.. [2] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
   Berlin; New York:: Springer, 2010.
.. [3] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
   Eighth Edition. McGraw-Hill Professional, 2007.

The structure of each dataframe is shown below:

.. ipython::

    In [1]: import chemicals

    In [2]: chemicals.viscosity.mu_data_Dutt_Prasad

    In [3]: chemicals.viscosity.mu_data_VN3

    In [4]: chemicals.viscosity.mu_data_VN2

    In [5]: chemicals.viscosity.mu_data_VN2E

    In [6]: chemicals.viscosity.mu_data_Perrys_8E_2_313

    In [7]: chemicals.viscosity.mu_data_Perrys_8E_2_312

    In [8]: chemicals.viscosity.mu_data_VDI_PPDS_7

    In [9]: chemicals.viscosity.mu_data_VDI_PPDS_8

"""

from __future__ import division

__all__ = ['Viswanath_Natarajan_3','Letsou_Stiel', 'Przedziecki_Sridhar', 'PPDS9', 'dPPDS9_dT',
'Viswanath_Natarajan_2', 'Viswanath_Natarajan_2_exponential', 'Lucas', 'Brokaw',
'Yoon_Thodos', 'Stiel_Thodos', 'Lucas_gas', 'viscosity_gas_Gharagheizi', 'Herning_Zipperer',
'Wilke', 'Wilke_prefactors', 'Wilke_prefactored', 'Wilke_large',
'viscosity_index', 'viscosity_converter', 'Lorentz_Bray_Clarke', 'Twu_1985', 'mu_IAPWS', 'mu_air_lemmon']

from fluids.numerics import secant, interp, numpy as np, trunc_exp
from chemicals.utils import log, exp, sqrt, atan, tan, sin, acos

from chemicals.utils import PY37, source_path, os_path_join, can_load_data
from chemicals.data_reader import register_df_source, data_source


folder = os_path_join(source_path, 'Viscosity')

register_df_source(folder, 'Dutt Prasad 3 term.tsv', csv_kwargs={
            'dtype':{'A': float, 'B': float, 'C': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Viswanath Natarajan Dynamic 3 term.tsv', csv_kwargs={
        'dtype':{'A': float, 'B': float, 'C': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Viswanath Natarajan Dynamic 2 term.tsv', csv_kwargs={
    'dtype':{'A': float, 'B': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Viswanath Natarajan Dynamic 2 term Exponential.tsv', csv_kwargs={
    'dtype':{'C': float, 'D': float, 'Tmin': float, 'Tmax': float}})

register_df_source(folder, 'Table 2-313 Viscosity of Inorganic and Organic Liquids.tsv')

register_df_source(folder, 'Table 2-312 Vapor Viscosity of Inorganic and Organic Substances.tsv', csv_kwargs={
    'dtype':{'C1': float, 'C2': float, 'C3': float, 'C4': float, 'Tmin': float, 'Tmax': float}})


register_df_source(folder, 'VDI PPDS Dynamic viscosity of saturated liquids polynomials.tsv', csv_kwargs={'float_precision': 'legacy'})
register_df_source(folder, 'VDI PPDS Dynamic viscosity of gases polynomials.tsv', csv_kwargs={'float_precision': 'legacy'})



_mu_data_loaded = False
def _load_mu_data():
    global _mu_data_loaded, mu_data_Dutt_Prasad, mu_values_Dutt_Prasad
    global mu_data_VN3, mu_values_VN3, mu_data_VN2, mu_values_VN2
    global mu_data_VN2E, mu_values_VN2E, mu_data_Perrys_8E_2_313, mu_values_Perrys_8E_2_313
    global mu_data_Perrys_8E_2_312, mu_values_Perrys_8E_2_312
    global mu_data_VDI_PPDS_7, mu_values_PPDS_7, mu_data_VDI_PPDS_8, mu_values_PPDS_8

    mu_data_Dutt_Prasad = data_source('Dutt Prasad 3 term.tsv')
    mu_values_Dutt_Prasad = np.array(mu_data_Dutt_Prasad.values[:, 1:], dtype=float)

    mu_data_VN3 = data_source('Viswanath Natarajan Dynamic 3 term.tsv')
    mu_values_VN3 = np.array(mu_data_VN3.values[:, 2:], dtype=float)

    mu_data_VN2 = data_source('Viswanath Natarajan Dynamic 2 term.tsv')
    mu_values_VN2 = np.array(mu_data_VN2.values[:, 2:], dtype=float)

    mu_data_VN2E = data_source('Viswanath Natarajan Dynamic 2 term Exponential.tsv')
    mu_values_VN2E = np.array(mu_data_VN2E.values[:, 2:], dtype=float)

    mu_data_Perrys_8E_2_313 = data_source('Table 2-313 Viscosity of Inorganic and Organic Liquids.tsv')
    mu_values_Perrys_8E_2_313 = np.array(mu_data_Perrys_8E_2_313.values[:, 1:], dtype=float)

    mu_data_Perrys_8E_2_312 = data_source('Table 2-312 Vapor Viscosity of Inorganic and Organic Substances.tsv')
    mu_values_Perrys_8E_2_312 = np.array(mu_data_Perrys_8E_2_312.values[:, 1:], dtype=float)

    mu_data_VDI_PPDS_7 = data_source('VDI PPDS Dynamic viscosity of saturated liquids polynomials.tsv')
    mu_values_PPDS_7 = np.array(mu_data_VDI_PPDS_7.values[:, 2:], dtype=float)

    mu_data_VDI_PPDS_8 = data_source('VDI PPDS Dynamic viscosity of gases polynomials.tsv')
    mu_values_PPDS_8 = np.array(mu_data_VDI_PPDS_8.values[:, 1:], dtype=float)

    _mu_data_loaded = True

if PY37:
    def __getattr__(name):
        if name in ('mu_data_Dutt_Prasad', 'mu_values_Dutt_Prasad', 'mu_data_VN3',
                    'mu_values_VN3', 'mu_data_VN2', 'mu_values_VN2', 'mu_data_VN2E',
                    'mu_values_VN2E', 'mu_data_Perrys_8E_2_313', 'mu_values_Perrys_8E_2_313',
                    'mu_data_Perrys_8E_2_312', 'mu_values_Perrys_8E_2_312', 'mu_data_VDI_PPDS_7',
                    'mu_values_PPDS_7', 'mu_data_VDI_PPDS_8', 'mu_values_PPDS_8'):
            _load_mu_data()
            return globals()[name]
        raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    if can_load_data:
        _load_mu_data()



def mu_IAPWS(T, rho, drho_dP=None, drho_dP_Tr=None):
    r'''Calculates and returns the viscosity of water according to the IAPWS
    (2008) release.

    Viscosity is calculated as a function of three terms;
    the first is the dilute-gas limit; the second is the contribution due to
    finite density; and the third and most complex is a critical enhancement
    term.

    .. math::
        \mu = \mu_0 \cdot \mu_1(T, \rho)
        \cdot \mu_2(T, \rho)

    .. math::
        \mu_0(T) = \frac{100\sqrt{T}}{\sum_{i=0}^3 \frac{H_i}{T^i}}

    .. math::
        \mu_1(T, \rho) = \exp\left[\rho \sum_{i=0}^5
        \left(\left(\frac{1}{T} - 1 \right)^i
        \sum_{j=0}^6 H_{ij}(\rho - 1)^j\right)\right]

    .. math::
        \text{if }\xi < 0.3817016416 \text{ nm:}

    .. math::
        Y = 0.2 q_c \xi(q_D \xi)^5 \left(1 - q_c\xi + (q_c\xi)^2 -
        \frac{765}{504}(q_D\xi)^2\right)

    .. math::
        \text{else:}

    .. math::
        Y = \frac{1}{12}\sin(3\psi_D) - \frac{1}{4q_c \xi}\sin(2\psi_D) +
        \frac{1}{(q_c\xi)^2}\left[1 - 1.25(q_c\xi)^2\right]\sin(\psi_D)
        - \frac{1}{(q_c\xi)^3}\left\{\left[1 - 1.5(q_c\xi)^2\right]\psi_D
        - \left|(q_c\xi)^2 - 1\right|^{1.5}L(w)\right\}

    .. math::
        w = \left| \frac{q_c \xi -1}{q_c \xi +1}\right|^{0.5} \tan\left(
        \frac{\psi_D}{2}\right)

    .. math::
        L(w) = \ln\frac{1 + w}{1 - w} \text{  if  }q_c \xi > 1

    .. math::
        L(w) = 2\arctan|w| \text{  if  }q_c \xi \le 1

    .. math::
        \psi_D = \arccos\left[\left(1 + q_D^2 \xi^2\right)^{-0.5}\right]

    .. math::
        \Delta \bar\chi(\bar T, \bar \rho) = \bar\rho\left[\zeta(\bar T, \bar
        \rho) - \zeta(\bar T_R, \bar \rho)\frac{\bar T_R}{\bar T}\right]

    .. math::
        \xi = \xi_0 \left(\frac{\Delta \bar\chi}{\Gamma_0}\right)^{\nu/\gamma}

    .. math::
        \zeta = \left(\frac{\partial\bar\rho}{\partial \bar p}\right)_{\bar T}

    Parameters
    ----------
    T : float
        Temperature of water [K]
    rho : float
        Density of water [kg/m^3]
    drho_dP : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature (at the temperature and density of water), [kg/m^3/Pa]
    drho_dP_Tr : float, optional
        Partial derivative of density with respect to pressure at constant
        temperature (at the reference temperature (970.644 K) and the actual
        density of water), [kg/m^3/Pa]

    Returns
    -------
    mu : float
        Viscosity, [Pa*s]

    Notes
    -----
    There are three ways to use this formulation.

    1) Compute the Industrial formulation value which does not include the
       critical enhacement, by leaving `drho_dP` and `drho_dP_Tr` None.
    2) Compute the Scientific formulation value by accurately computing and
       providing `drho_dP` and `drho_dP_Tr`, both with IAPWS-95.
    3) Get a non-standard but 8 decimal place matching result by providing
       `drho_dP` computed with either IAPWS-95 or IAPWS-97, but not providing
       `drho_dP_Tr`; which is calculated internally. There is a formulation
       for that term in the thermal conductivity IAPWS equation which is used.

    xmu = 0.068

    qc = (1.9E-9)**-1

    qd = (1.1E-9)**-1

    nu = 0.630

    gamma = 1.239

    xi0 = 0.13E-9

    Gamma0 = 0.06

    TRC = 1.5

    This forulation is highly optimized, spending most of its time in the
    logarithm, power, and square root.

    Examples
    --------
    >>> mu_IAPWS(298.15, 998.)
    0.000889735100149808

    >>> mu_IAPWS(1173.15, 400.)
    6.415460784836147e-05

    Point 4 of formulation, compared with MPEI and IAPWS, matches.

    >>> mu_IAPWS(T=647.35, rho=322., drho_dP=1.213641949033E-2)
    4.2961578738287e-05

    Full scientific calculation:

    >>> from chemicals.iapws import iapws95_properties, iapws95_P, iapws95_Tc
    >>> T, P = 298.15, 1e5
    >>> rho, _, _, _, _, _, _, _, _, _, drho_dP = iapws95_properties(T, P)
    >>> P_ref = iapws95_P(1.5*iapws95_Tc, rho)
    >>> _, _, _, _, _, _, _, _, _, _, drho_dP_Tr = iapws95_properties(1.5*iapws95_Tc, P_ref)
    >>> mu_IAPWS(T, rho, drho_dP, drho_dP_Tr)
    0.00089002267377


    References
    ----------
    .. [1] Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
       Sengers, M. J. Assael, I. N. Metaxa, E. Vogel, R. Mares, and
       K. Miyagawa. "New International Formulation for the Viscosity of H2O."
       Journal of Physical and Chemical Reference Data 38, no. 2
       (June 1, 2009): 101-25. doi:10.1063/1.3088050.
    '''
    Tr = T*0.0015453657571674064 #/647.096
    Tr_inv = 1.0/Tr
    rhor = rho*0.003105590062111801 #1/322.
    x0 = rhor - 1.
    x1 = Tr_inv - 1.
#    His = [1.67752, 2.20462, 0.6366564, -0.241605]
#    mu0 = 0
#    for i in range(4):
#        mu0 += His[i]/Tr**i
    mu0 = 100.0*sqrt(Tr)/(Tr_inv*(Tr_inv*(0.6366564 - 0.241605*Tr_inv) + 2.20462) + 1.67752)

#    i_coefs = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6]
#    j_coef = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5]
#    Hijs = [0.520094, .0850895, -1.08374, -0.289555, 0.222531, 0.999115,
#          1.88797, 1.26613, 0.120573, -0.281378, -0.906851, -0.772479,
#          -0.489837, -0.257040, 0.161913, 0.257399, -0.0325372, 0.0698452,
#          0.00872102, -0.00435673, -0.000593264]
#    tot = 0
#    for i in range(21):
#        tot += Hijs[i]*(rhor - 1.)**i_coefs[i]*(Tr_inv - 1.)**j_coef[i]
    x02 = x0*x0
    tot = (x0*(x0*(x0*(0.161913 - 0.0325372*x0) - 0.281378) + 0.222531)
           + x1*(x0*(x0*(0.257399*x0 - 0.906851) + 0.999115) + x1*(x0*(1.88797 - 0.772479*x0)
            + x1*(x0*(x0*(x02*(0.0698452 - 0.00435673*x02) - 0.489837) + 1.26613)
            + x1*(x02*(0.00872102*x0*x02 - 0.25704) + x0*x1*(0.120573 - 0.000593264*x02*x02*x0)) - 0.289555)
            - 1.08374) + 0.0850895) + 0.520094)
    mu1 = exp(rhor*tot)

    if drho_dP is not None:
        xmu = 0.068
        qc = 526315789.4736842#(1.9E-9)**-1
        qD = 909090909.0909091#(1.1E-9)**-1
#        nu = 0.630
#        gamma = 1.239
        xi0 = 0.13E-9
#        Gamma0 = 0.06
        TRC = 1.5

        # Not a perfect match because of
        zeta_drho_dP = drho_dP*68521.73913043478 #22.064E6/322.0
        if drho_dP_Tr is None:
            # Brach needed to allow scientific points to work
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
            drho_dP_Tr2 = 1./tot1
        else:
            drho_dP_Tr2 = drho_dP_Tr*68521.73913043478 #22.064E6/322.0
        dchi = rhor*(zeta_drho_dP - drho_dP_Tr2*TRC*Tr_inv)
        if dchi < 0.0:
            # By definition
            return mu0*mu1*1e-6

        # 16.666 = 1/Gamma0
        xi = xi0*(dchi*16.666666666666668)**0.5084745762711864 #(nu/gamma)
        qD2 = qD*qD
        xi2 = xi*xi
        x2 = qD2*xi2
        psiD = acos(1.0/sqrt(1.0 + x2))
        qcxi = qc*xi
        qcxi2 = qcxi*qcxi
        qcxi_inv = 1.0/qcxi

        w = sqrt(abs((qcxi - 1.0)/(qcxi + 1.0)))*tan(psiD*0.5)
        if qc*xi > 1.0:
            Lw = log((1.0 + w)/(1.0 - w))
        else:
            Lw = 2.0*atan(w)

        if xi <= 0.381706416E-9:
            # 1.5178571428571428 = 765./504
            Y = 0.2*qcxi*x2*x2*qD*xi*(1.0 - qcxi + qcxi2 - 1.5178571428571428*x2)
        else:
            # sin(ax) = 2cos(ax/2)*sin(ax/2)
            # It would be possible to compute the sin(2psiD) and sin(psiD) together with a sincos
            # operation, but not the sin(3psid)
            x3 = (abs(qcxi2 - 1.0))
            Y = (1/12.*sin(3.0*psiD) + (-0.25*sin(2.0*psiD)
                 +((1.0 - 1.25*qcxi2)*sin(psiD)
                 -qcxi_inv*((1.0 - 1.5*qcxi2)*psiD - x3*sqrt(x3)*Lw))*qcxi_inv)*qcxi_inv )

        mu2 = exp(xmu*Y)
    else:
        mu2 = 1.0
    mu = mu0*mu1*mu2*1e-6
    return mu


def mu_air_lemmon(T, rho):
    r'''Calculates and returns the viscosity of air according to Lemmon
    and Jacobsen (2003) [1]_.

    Viscosity is calculated as a function of two terms;
    the first is the dilute-gas limit; the second is the contribution due to
    finite density.

    .. math::
        \mu = \mu^0(T) + \mu^r(T, \rho)

    .. math::
        \mu^0(T) = \frac{0.9266958\sqrt{MT}}{\sigma^2 \Omega(T^*)}

    .. math::
        \Omega(T^*) = \exp\left( \sum_{i=0}^4 b_i [\ln(T^*)]^i \right)

    .. math::
        \mu^r = \sum_{i=1}^n N_i \tau^{t_i} \delta^{d_i} \exp\left(
        -\gamma_i \delta^{l_i}\right)

    Parameters
    ----------
    T : float
        Temperature of air [K]
    rho : float
        Molar density of air [mol/m^3]

    Returns
    -------
    mu : float
        Viscosity of air, [Pa*s]

    Notes
    -----

    The coefficients are:

    Ni = [10.72, 1.122, 0.002019, -8.876, -0.02916]

    ti = [0.2, 0.05, 2.4, 0.6, 3.6]

    di = [1, 4, 9, 1, 8]

    gammai = Ii = [0, 0, 0, 1, 1]

    bi = [.431, -0.4623, 0.08406, 0.005341, -0.00331]

    The reducing parameters are :math:`T_c = 132.6312` K and
    :math:`\rho_c = 10447.7` mol/m^3. Additional parameters used are
    :math:`\sigma = 0.36` nm,
    :math:`M = 28.9586` g/mol and :math:`\frac{e}{k} = 103.3` K.

    This is an implementation optimized for speed, spending its time
    in the calclulation of 1 log; 2 exp; 1 power; and 2 divisions.

    Examples
    --------
    Viscosity at 300 K and 1 bar:

    >>> mu_air_lemmon(300.0, 40.10292351061862)
    1.85371518556e-05

    Calculate the density in-place:

    >>> from chemicals.air import lemmon2000_rho
    >>> mu_air_lemmon(300.0, lemmon2000_rho(300.0, 1e5))
    1.85371518556e-05

    References
    ----------
    .. [1] Lemmon, E. W., and R. T. Jacobsen. "Viscosity and Thermal
       Conductivity Equations for Nitrogen, Oxygen, Argon, and Air."
       International Journal of Thermophysics 25, no. 1 (January 1, 2004):
       21-69. https://doi.org/10.1023/B:IJOT.0000022327.04529.f3.
    '''
    # Cost: 1 log; 2 exp; 1 power; 2 divisions
#     sigma = 0.360 # nm
#     M = 28.9586 # g/mol
#     rhoc = 10447.7 # mol/m^3, maxcondentherm actually
    Tc = 132.6312 # K, maxcondentherm actually
    tau = Tc/T
    delta = rho*9.571484632981421e-05 # 9.57...E-5 = 1/10447.7

    delta2 = delta*delta
    delta4 = delta2*delta2
    delta8 = delta4*delta4
    tau_20 = tau**0.05
    tau2_20 = tau_20*tau_20
    tau4_20 = tau2_20*tau2_20
    tau8_20 = tau4_20*tau4_20

    tau12_20 = tau4_20*tau8_20
    tau24_20 = tau12_20*tau12_20
    tau48_20 = tau24_20*tau24_20
    x0 = exp(-delta)
    etar = (delta*(-8.876e-6*tau12_20*x0 + 0.002019e-6*tau48_20*delta8
                  + 10.72e-6*tau4_20) + 1.122e-6*delta4*tau_20
            - 0.02916e-6*delta8*tau24_20*x0*tau48_20)

#     e_k = 103.3 # K
#     Ts = T/e_k
    Ts = T*0.00968054211035818 # 1/e_k
    lnTs = log(Ts)
#     tot = 0.0
#     for i in range(5):
#         tot += CIs[i]*lnTs**i
    Omega = exp(lnTs*(lnTs*(lnTs*(0.005341 - 0.00331*lnTs) + 0.08406) - 0.4623) + 0.431)

    # 0.0266958*sqrt(28.9586)/(0.360*0.360)*sqrt(132.6312) = 12.76...
    eta0 = 12.765845058845755e-6/(Omega*tau8_20*tau2_20)
#     eta0 = 0.0266958*sqrt(T*M)/(sigma*sigma*Omega)
#     etar = 0.0
#     for i in range(5):
#         etar += Ni[i]*tau**ti[i]*delta**di[i]*exp(-gammai[i]*delta**Ii[i])
    return (eta0 + etar)


def Viswanath_Natarajan_2(T, A, B):
    r'''Calculate the viscosity of a liquid using the 2-term form
    representation developed in [1]_. Requires input coefficients. The `A`
    coefficient is assumed to yield coefficients in Pa*s; if it yields
    values in 1E-3 Pa*s, remove log(100) for A.

    .. math::
        \mu = \exp\left(A + \frac{B}{T}\right)

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----
    No other source for these coefficients than [1]_ has been found.

    Examples
    --------
    DDBST has 0.0004580 as a value at this temperature for 1-Butanol.

    >>> Viswanath_Natarajan_2(348.15, -5.9719-log(100), 1007.0)
    0.000459836869568295

    References
    ----------
    .. [1] Viswanath, Dabir S., and G. Natarajan. Databook On The Viscosity Of
       Liquids. New York: Taylor & Francis, 1989
    '''
    return exp(A + B/T)


def Viswanath_Natarajan_2_exponential(T, C, D):
    r'''Calculate the viscosity of a liquid using the 2-term exponential form
    representation developed in [1]_. Requires input coefficients. The `A`
    coefficient is assumed to yield coefficients in Pa*s, as all
    coefficients found so far have been.

    .. math::
        \mu = C T^D

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    C : float
        Linear coefficient, [Pa*s]
    D : float
        Exponential coefficient, [-]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----
    No other source for these coefficients has been found.

    Examples
    --------
    >>> Ts = [283.15, 288.15, 303.15, 349.65]
    >>> mus = [2.2173, 2.1530, 1.741, 1.0091] # in cP
    >>> Viswanath_Natarajan_2_exponential(288.15, 4900800, -3.8075)
    0.002114798866203873

    Calculation of the AARD of the fit (1% is the value stated in [1]_.:

    >>> mu_calc = [Viswanath_Natarajan_2_exponential(T, 4900800, -3.8075) for T in Ts]
    >>> np.mean([abs((mu - mu_i*1000)/mu) for mu, mu_i in zip(mus, mu_calc)])
    0.010467928813061298

    References
    ----------
    .. [1] Viswanath, Dabir S., and G. Natarajan. Databook On The Viscosity Of
       Liquids. New York: Taylor & Francis, 1989
    '''
    return C*T**D


def Viswanath_Natarajan_3(T, A, B, C):
    r'''Calculate the viscosity of a liquid using the 3-term Antoine form
    representation developed in [1]_. Requires input coefficients. If the
    coefficients do not yield viscosity in Pa*s, but rather cP, remove
    log10(1000) from `A`.

    .. math::
        \log_{10} \mu = A + B/(T + C)

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]
    C : float
        Coefficient, [K]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----
    No other source for these coefficients has been found.

    Examples
    --------
    >>> from math import log10
    >>> Viswanath_Natarajan_3(298.15, -2.7173-log10(1000), -1071.18, -129.51)
    0.0006129806445142113

    References
    ----------
    .. [1] Viswanath, Dabir S., and G. Natarajan. Databook On The Viscosity Of
       Liquids. New York: Taylor & Francis, 1989
    '''
    return 10.0**(A + B/(C - T))

def PPDS9(T, A, B, C, D, E):
    r'''Calculate the viscosity of a liquid using the 5-term exponential power
    fit developed by the PPDS and named PPDS equation 9.

    .. math::
       \mu = E \exp\left[A \left(\frac{C-T}{T-D}\right)^{1/3}
        + B \left(\frac{C-T}{T-D}\right)^{4/3}  \right]

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [-]
    C : float
        Coefficient, [K]
    D : float
        Coefficient, [K]
    E : float
        Coefficient, [Pa*s]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----
    No other source for these coefficients has been found.

    There can be a singularity in this equation when `T` approaches `C` or
    `D`; it may be helpful to take as a limit to this equation `D` + 5 K.

    Examples
    --------
    >>> PPDS9(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    0.00035091137378230684

    References
    ----------
    .. [1] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
       Berlin; New York:: Springer, 2010.
    '''
    term = (C - T)/(T-D)
    if term < 0.0:
        term1 = -((T - C)/(T-D))**(1/3.)
    else:
        term1 = term**(1/3.)
    term2 = term*term1
    mu = E*trunc_exp(A*term1 + B*term2)
    return mu

def dPPDS9_dT(T, A, B, C, D, E):
    r'''Calculate the temperature derivative of  viscosity of a liquid using
    the 5-term exponential power fit developed by the PPDS and named PPDS
    equation 9.

    Normally, the temperature derivative is:

    .. math::
        \frac{\partial \mu}{\partial T} = E \left(\frac{A \sqrt[3]{\frac{C - T}
        {- D + T}} \left(- D + T\right) \left(- \frac{C - T}{3 \left(- D
        + T\right)^{2}} - \frac{1}{3 \left(- D + T\right)}\right)}{C - T}
        - \frac{B \sqrt[3]{\frac{C - T}{- D + T}} \left(C - T\right)}{\left(
        - D + T\right)^{2}} + B \sqrt[3]{\frac{C - T}{- D + T}} \left(- \frac{
        C - T}{3 \left(- D + T\right)^{2}} - \frac{1}{3 \left(- D + T\right)}
        \right) - \frac{B \sqrt[3]{\frac{C - T}{- D + T}}}{- D + T}\right)
        e^{A \sqrt[3]{\frac{C - T}{- D + T}} + \frac{B \sqrt[3]{\frac{C - T}
        {- D + T}} \left(C - T\right)}{- D + T}}

    For the low-temperature region:

    .. math::
        \frac{\partial \mu}{\partial T} = E \left(- \frac{A \sqrt[3]{\frac{
        - C + T}{- D + T}} \left(- D + T\right) \left(- \frac{- C + T}{3
        \left(- D + T\right)^{2}} + \frac{1}{3 \left(- D + T\right)}\right)
        }{- C + T} + \frac{B \sqrt[3]{\frac{- C + T}{- D + T}} \left(C
        - T\right)}{\left(- D + T\right)^{2}} + \frac{B \sqrt[3]{\frac{
        - C + T}{- D + T}}}{- D + T} - \frac{B \sqrt[3]{\frac{- C + T}{
        - D + T}} \left(C - T\right) \left(- \frac{- C + T}{3 \left(- D
        + T\right)^{2}} + \frac{1}{3 \left(- D + T\right)}\right)}{- C
        + T}\right) e^{- A \sqrt[3]{\frac{- C + T}{- D + T}} - \frac{B
        \sqrt[3]{\frac{- C + T}{- D + T}} \left(C - T\right)}{- D + T}}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [-]
    C : float
        Coefficient, [K]
    D : float
        Coefficient, [K]
    E : float
        Coefficient, [Pa*s]

    Returns
    -------
    dmu_dT : float
        First temperature derivative of liquid viscosity, [Pa*s]
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----

    Examples
    --------
    >>> dPPDS9_dT(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    (-3.186540635882627e-06, 0.00035091137378230684)

    References
    ----------
    .. [1] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
       Berlin; New York:: Springer, 2010.
    '''
    term = (C - T)/(T-D)
    if term < 0.0:
        x0 = 1.0/(-D + T)
        x1 = x0*(-C + T)
        x2 = -T
        x3 = C + x2
        x4 = B*x3
        mu = E*trunc_exp(-x1**(1.0/3.0)*(A + x0*x4))
        x6 = D + x2
        x7 = 1.0/x6
        x8 = x0*(x1 - 1.0)/3.0
        dmu_dT = -mu*(x3*x7)**(1.0/3.0)*(-A*x6*x8/x3 + B*x7 + B*x8 - x4*x7*x7)
    else:
        x0 = -T
        x1 = C + x0
        x2 = D + x0
        x3 = 1.0/x2
        x4 = x1*x3
        x5 = (-x4)**(1.0/3.0)
        mu = E*trunc_exp(x5*(A - B*x4))
        x7 = 1.0/(-D + T)
        x8 = x7*(x1*x7 + 1.0)*(1.0/3.0)
        dmu_dT = -x5*mu*(-A*x2*x8/x1 + B*x1*x3*x3 - B*x3 + B*x8)
    return (dmu_dT, mu)


def Letsou_Stiel(T, MW, Tc, Pc, omega):
    r'''Calculates the viscosity of a liquid using an emperical model
    developed in [1]_. However. the fitting parameters for tabulated values
    in the original article are found in ChemSep.

    .. math::
        \xi = \frac{2173.424 T_c^{1/6}}{\sqrt{MW} P_c^{2/3}}

    .. math::
        \xi^{(0)} = (1.5174 - 2.135T_r + 0.75T_r^2)\cdot 10^{-5}

    .. math::
        \xi^{(1)} = (4.2552 - 7.674 T_r + 3.4 T_r^2)\cdot 10^{-5}

    .. math::
        \mu = (\xi^{(0)} + \omega \xi^{(1)})/\xi

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    MW : float
        Molwcular weight of fluid [g/mol]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    omega : float
        Acentric factor of compound

    Returns
    -------
    mu_l : float
        Viscosity of liquid, [Pa*s]

    Notes
    -----
    The form of this equation is a polynomial fit to tabulated data.
    The fitting was performed by the DIPPR. This is DIPPR Procedure 8G: Method
    for the viscosity of pure, nonhydrocarbon liquids at high temperatures
    internal units are SI standard. [1]_'s units were different.
    DIPPR test value for ethanol is used.

    Average error 34%. Range of applicability is 0.76 < Tr < 0.98.

    Examples
    --------
    >>> Letsou_Stiel(400., 46.07, 516.25, 6.383E6, 0.6371)
    0.0002036150875308

    References
    ----------
    .. [1] Letsou, Athena, and Leonard I. Stiel. "Viscosity of Saturated
       Nonpolar Liquids at Elevated Pressures." AIChE Journal 19, no. 2 (1973):
       409-11. doi:10.1002/aic.690190241.
    '''
    Tr = T/Tc
    xi0 = (1.5174 - Tr*(2.135 - 0.75*Tr))*1E-5
    xi1 = (4.2552 - Tr*(7.674 - 3.4*Tr))*1E-5
    xi = 2173.424*Tc**(1.0/6.)/sqrt(MW)*Pc**(-2.0/3.)
    return (xi0 + omega*xi1)/xi


def Przedziecki_Sridhar(T, Tm, Tc, Pc, Vc, Vm, omega, MW):
    r'''Calculates the viscosity of a liquid using an emperical formula
    developed in [1]_.

    .. math::
        \mu=\frac{V_o}{E(V-V_o)}

    .. math::
        E=-1.12+\frac{V_c}{12.94+0.10MW-0.23P_c+0.0424T_{m}-11.58(T_{m}/T_c)}

    .. math::
        V_o = 0.0085\omega T_c-2.02+\frac{V_{m}}{0.342(T_m/T_c)+0.894}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tm : float
        Melting point of fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    Vc : float
        Critical volume of the fluid [m^3/mol]
    Vm : float
        Molar volume of the fluid at temperature [K]
    omega : float
        Acentric factor of compound
    MW : float
        Molwcular weight of fluid [g/mol]

    Returns
    -------
    mu_l : float
        Viscosity of liquid, [Pa*s]

    Notes
    -----
    A test by Reid (1983) is used, but only mostly correct.
    This function is not recommended.
    Internal units are bar and mL/mol.

    Examples
    --------
    >>> Przedziecki_Sridhar(383., 178., 591.8, 41E5, 316E-6, 95E-6, .263, 92.14)
    0.00021981479956033846

    References
    ----------
    .. [1] Przedziecki, J. W., and T. Sridhar. "Prediction of Liquid
       Viscosities." AIChE Journal 31, no. 2 (February 1, 1985): 333-35.
       doi:10.1002/aic.690310225.
    '''
    Pc = Pc*1e-5  # Pa to atm
    Vm, Vc = Vm*1E6, Vc*1E6  # m^3/mol to mL/mol
    Tc_inv = 1.0/Tc
    Tr = T*Tc_inv

    Tr2 = Tr*Tr
    Gamma = 0.29607 - 0.09045*Tr - 0.04842*Tr2
    VrT = 0.33593 - 0.33953*Tr + 1.51941*Tr2 - 2.02512*Tr*Tr2 + 1.11422*Tr2*Tr2
    V = VrT*(1.0 - omega*Gamma)*Vc

    Vo = 0.0085*omega*Tc - 2.02 + Vm/(0.342*(Tm*Tc_inv) + 0.894)  # checked
    E = -1.12 + Vc/(12.94 + 0.1*MW - 0.23*Pc + 0.0424*Tm - 11.58*(Tm*Tc_inv))
    return Vo/(E*(V-Vo))*1e-3

### Viscosity of Dense Liquids


def Lucas(T, P, Tc, Pc, omega, Psat, mu_l):
    r'''Adjustes for pressure the viscosity of a liquid using an emperical
    formula developed in [1]_, but as discussed in [2]_ as the original source
    is in German.

    .. math::
        \frac{\mu}{\mu_{sat}}=\frac{1+D(\Delta P_r/2.118)^A}{1+C\omega \Delta P_r}

    .. math::
        \Delta P_r = \frac{P-P^{sat}}{P_c}

    .. math::
        A=0.9991-\frac{4.674\times 10^{-4}}{1.0523T_r^{-0.03877}-1.0513}

    .. math::
        D = \frac{0.3257}{(1.0039-T_r^{2.573})^{0.2906}}-0.2086

    .. math::
        C = -0.07921+2.1616T_r-13.4040T_r^2+44.1706T_r^3-84.8291T_r^4+
        96.1209T_r^5-59.8127T_r^6+15.6719T_r^7

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
    omega : float
        Acentric factor of compound
    Psat : float
        Saturation pressure of the fluid [Pa]
    mu_l : float
        Viscosity of liquid at 1 atm or saturation, [Pa*s]

    Returns
    -------
    mu_l_dense : float
        Viscosity of liquid, [Pa*s]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The example is from Reid (1987); all results agree.
    Above several thousand bar, this equation does not represent true behavior.
    If Psat is larger than P, the fluid may not be liquid; dPr is set to 0.

    Examples
    --------
    >>> Lucas(300., 500E5, 572.2, 34.7E5, 0.236, 0, 0.00068) # methylcyclohexane
    0.0010683738499316494

    References
    ----------
    .. [1] Lucas, Klaus. "Ein Einfaches Verfahren Zur Berechnung Der
       Viskositat von Gasen Und Gasgemischen." Chemie Ingenieur Technik 46, no. 4
       (February 1, 1974): 157-157. doi:10.1002/cite.330460413.
    .. [2] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Tr = min(T/Tc, 1.0)
    C = Tr*(Tr*(Tr*(Tr*(Tr*(Tr*(15.6719*Tr - 59.8127) + 96.1209) - 84.8291) + 44.1706) - 13.404) + 2.1616) - 0.07921
    D = 0.3257*(1.0039-Tr**2.573)**-0.2906 - 0.2086
    A = 0.9991 - 4.674E-4/(1.0523*Tr**-0.03877 - 1.0513)
    dPr = (P-Psat)/Pc
    if dPr < 0.0:
        dPr = 0.0
    return (1. + D*(dPr*(1.0/2.118))**A)/(1. + C*omega*dPr)*mu_l

### Viscosity of liquid mixtures
### Viscosity of Gases - low pressure

def Yoon_Thodos(T, Tc, Pc, MW):
    r'''Calculates the viscosity of a gas using an emperical formula
    developed in [1]_.

    .. math::
        \eta \xi \times 10^8 = 46.10 T_r^{0.618} - 20.40 \exp(-0.449T_r) + 1
        9.40\exp(-4.058T_r)+1

    .. math::
        \xi = 2173.424 T_c^{1/6} MW^{-1/2} P_c^{-2/3}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    MW : float
        Molwcular weight of fluid [g/mol]

    Returns
    -------
    mu_g : float
        Viscosity of gas, [Pa*s]

    Notes
    -----
    This equation has been tested. The equation uses SI units only internally.
    The constant 2173.424 is an adjustment factor for units.
    Average deviation within 3% for most compounds.
    Greatest accuracy with dipole moments close to 0.
    Hydrogen and helium have different coefficients, not implemented.
    This is DIPPR Procedure 8B: Method for the Viscosity of Pure,
    non hydrocarbon, nonpolar gases at low pressures

    Examples
    --------
    >>> Yoon_Thodos(300., 556.35, 4.5596E6, 153.8)
    1.019488572777e-05

    References
    ----------
    .. [1] Yoon, Poong, and George Thodos. "Viscosity of Nonpolar Gaseous
       Mixtures at Normal Pressures." AIChE Journal 16, no. 2 (1970): 300-304.
       doi:10.1002/aic.690160225.
    '''
    Tr = T/Tc
    xi = 2173.4241*Tc**(1/6.)/sqrt(MW)*Pc**(-2.0/3.)
    a = 46.1
    b = 0.618
    c = 20.4
    d = -0.449
    e = 19.4
    f = -4.058
    return (1. + a*Tr**b - c * exp(d*Tr) + e*exp(f*Tr))/(1E8*xi)


def Stiel_Thodos(T, Tc, Pc, MW):
    r'''Calculates the viscosity of a gas using an emperical formula
    developed in [1]_.

    if :math:`T_r > 1.5`:

    .. math::
        \mu_g = 17.78\times 10^{-5} (4.58T_r - 1.67)^{0.625}/\xi

    else:

    .. math::
        \mu_g = 34\times 10^{-5} T_r^{0.94}/\xi

    .. math::
        \xi = \frac{T_c^{(1/6)}}{\sqrt{MW} P_c^{2/3}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    MW : float
        Molwcular weight of fluid [g/mol]

    Returns
    -------
    mu_g : float
        Viscosity of gas, [Pa*s]

    Notes
    -----
    Claimed applicability from 0.2 to 5 atm.
    Developed with data from 52 nonpolar, and 53 polar gases.
    internal units are poise and atm.
    Seems to give reasonable results.

    Examples
    --------
    >>> Stiel_Thodos(300., 556.35, 4.5596E6, 153.8) #CCl4
    1.040892622360e-05

    References
    ----------
    .. [1] Stiel, Leonard I., and George Thodos. "The Viscosity of Nonpolar
       Gases at Normal Pressures." AIChE Journal 7, no. 4 (1961): 611-15.
       doi:10.1002/aic.690070416.
    '''
    Pc = Pc*(1.0/101325.)
    Tr = T/Tc
    xi = Tc**(1/6.)/(sqrt(MW)*Pc**(2/3.))
    if Tr > 1.5:
        mu_g = 17.78E-5*(4.58*Tr-1.67)**0.625/xi
    else:
        mu_g = 34E-5*Tr**0.94/xi
    return mu_g*1e-3


def Lucas_gas(T, Tc, Pc, Zc, MW, dipole=0.0, CASRN=None):
    r'''Estimate the viscosity of a gas using an emperical
    formula developed in several sources, but as discussed in [1]_ as the
    original sources are in German or merely personal communications with the
    authors of [1]_.

    .. math::
        \eta  = \left[0.807T_r^{0.618}-0.357\exp(-0.449T_r) + 0.340\exp(-4.058
        T_r) + 0.018\right]F_p^\circ F_Q^\circ /\xi

    .. math::
        F_p^\circ=1, 0 \le \mu_{r} < 0.022

    .. math::
        F_p^\circ = 1+30.55(0.292-Z_c)^{1.72}, 0.022 \le \mu_{r} < 0.075

    .. math::
        F_p^\circ = 1+30.55(0.292-Z_c)^{1.72}|0.96+0.1(T_r-0.7)| 0.075 < \mu_{r}

    .. math::
        F_Q^\circ = 1.22Q^{0.15}\left\{ 1+0.00385[(T_r-12)^2]^{1/M}\text{sign}
        (T_r-12)\right\}

    .. math::
        \mu_r = 52.46 \frac{\mu^2 P_c}{T_c^2}

    .. math::
        \xi=0.176\left(\frac{T_c}{MW^3 P_c^4}\right)^{1/6}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc: float
        Critical point of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    Zc : float
        Critical compressibility of the fluid [Pa]
    dipole : float
        Dipole moment of fluid [debye]
    CASRN : str, optional
        CAS of the fluid

    Returns
    -------
    mu_g : float
        Viscosity of gas, [Pa*s]

    Notes
    -----
    The example is from [1]_; all results agree.
    Viscosity is calculated in micropoise, and converted to SI internally (1E-7).
    Q for He = 1.38; Q for H2 = 0.76; Q for D2 = 0.52.

    Examples
    --------
    >>> Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)
    1.7822676912698925e-05

    References
    ----------
    .. [1] Reid, Robert C.; Prausnitz, John M.; Poling, Bruce E.
       Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Tc_inv = 1.0/Tc
    Tr = T*Tc_inv
    MW_inv = 1.0/MW
    Pc_bar = Pc*1e-5
    xi = 0.176*(Tc*MW_inv*MW_inv*MW_inv/(Pc_bar*Pc_bar*Pc_bar*Pc_bar))**(1.0/6.0)  # bar arrording to example in Poling
    if dipole is None:
        dipole = 0.0
    dipoler = 52.46*dipole*dipole*Pc_bar*Tc_inv*Tc_inv  # bar arrording to example in Poling
    if dipoler < 0.022:
        Fp = 1.0
    elif 0.022 <= dipoler < 0.075:
        Fp = 1.0 + 30.55*(0.292 - Zc)**1.72
    else:
        Fp = 1.0 + 30.55*(0.292 - Zc)**1.72*abs(0.96 + 0.1*(Tr - 0.7))
    FQ = 1.0
    if CASRN is not None:
        Q = 0.0
        if CASRN == '7440-59-7':
            Q = 1.38
        elif CASRN == '1333-74-0':
            Q = 0.76
        elif CASRN == '7782-39-0':
            Q = 0.52
        if Q != 0.0:
            if Tr - 12.0 > 0.0:
                value = 1.0
            else:
                value = -1.0
            x0 = (Tr-12.0)
            FQ = 1.22*Q**0.15*(1.0 + 0.00385*(x0*x0)**(MW_inv)*value)
    eta = (0.807*Tr**0.618 - 0.357*exp(-0.449*Tr) + 0.340*exp(-4.058*Tr) + 0.018)*Fp*FQ/xi
    return eta*1E-7


def viscosity_gas_Gharagheizi(T, Tc, Pc, MW):
    r'''Calculates the viscosity of a gas using an emperical formula
    developed in [1]_.

    .. math::
        \mu = 10^{-7} | 10^{-5} P_cT_r + \left(0.091-\frac{0.477}{M}\right)T +
        M \left(10^{-5}P_c-\frac{8M^2}{T^2}\right)
        \left(\frac{10.7639}{T_c}-\frac{4.1929}{T}\right)|

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    Tc : float
        Critical temperature of the fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    MW : float
        Molwcular weight of fluid [g/mol]

    Returns
    -------
    mu_g : float
        Viscosity of gas, [Pa*s]

    Notes
    -----
    Example is first point in supporting information of article, for methane.
    This is the prefered function for gas viscosity.
    7% average relative deviation. Deviation should never be above 30%.
    Developed with the DIPPR database. It is believed theoretically predicted values
    are included in the correlation.

    Under 0.2Tc, this correlation has been modified to provide values at the
    limit.

    Examples
    --------
    >>> viscosity_gas_Gharagheizi(120., 190.564, 45.99E5, 16.04246)
    5.215761625399613e-06

    References
    ----------
    .. [1] Gharagheizi, Farhad, Ali Eslamimanesh, Mehdi Sattari, Amir H.
       Mohammadi, and Dominique Richon. "Corresponding States Method for
       Determination of the Viscosity of Gases at Atmospheric Pressure."
       Industrial & Engineering Chemistry Research 51, no. 7
       (February 22, 2012): 3179-85. doi:10.1021/ie202591f.
    '''
    Tr = T/Tc
    if Tr < 0.2:
        Tr = 0.2
        T = 0.2*Tc

    mu_g = 1E-5*Pc*Tr + (0.091 - 0.477/MW)*T + MW*(1E-5*Pc - 8.0*MW*MW/(T*T))*(10.7639/Tc - 4.1929/T)
    mu_g = 1e-7*mu_g
    return mu_g

### Viscosity of gas mixtures


def Herning_Zipperer(zs, mus, MWs, MW_roots=None):
    r'''Calculates viscosity of a gas mixture according to
    mixing rules in [1]_.

    .. math::
        \mu = \frac{\sum x_i \mu_i \sqrt{MW_i}}
        {\sum x_i \sqrt{MW_i}}

    Parameters
    ----------
    zs : float
        Mole fractions of components, [-]
    mus : float
        Gas viscosities of all components, [Pa*s]
    MWs : float
        Molecular weights of all components, [g/mol]
    MW_roots : float, optional
        Square roots of molecular weights of all components, [g^0.5/mol^0.5]

    Returns
    -------
    mug : float
        Viscosity of gas mixture, [Pa*s]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The original source has not been reviewed.

    Adding the square roots can speed up the calculation.

    Examples
    --------
    >>> Herning_Zipperer([0.5, 0.25, 0.25], [1.78e-05, 1.12e-05, 9.35e-06], [28.0134, 16.043, 30.07])
    1.4174908599465168e-05

    References
    ----------
    .. [1] Herning, F. and Zipperer, L,: "Calculation of the Viscosity of
       Technical Gas Mixtures from the Viscosity of Individual Gases, german",
       Gas u. Wasserfach (1936) 79, No. 49, 69.
    '''
    N = len(zs)
    if MW_roots is None:
        MW_roots = [0.0]*N
        for i in range(N):
            MW_roots[i] = sqrt(MWs[i])
    denominator = k = 0.0
    for i in range(N):
        v = zs[i]*MW_roots[i]
        k += v*mus[i]
        denominator += v
    return k/denominator

def Wilke(ys, mus, MWs):
    r'''Calculates viscosity of a gas mixture according to
    mixing rules in [1]_.

    .. math::
        \eta_{mix} = \sum_{i=1}^n \frac{y_i \eta_i}{\sum_{j=1}^n y_j \phi_{ij}}

    .. math::
        \phi_{ij} = \frac{(1 + \sqrt{\eta_i/\eta_j}(MW_j/MW_i)^{0.25})^2}
        {\sqrt{8(1+MW_i/MW_j)}}

    Parameters
    ----------
    ys : float
        Mole fractions of gas components, [-]
    mus : float
        Gas viscosities of all components, [Pa*s]
    MWs : float
        Molecular weights of all components, [g/mol]

    Returns
    -------
    mug : float
        Viscosity of gas mixture, [Pa*s]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The original source has not been reviewed or found.

    See Also
    --------
    Wilke_prefactors
    Wilke_prefactored
    Wilke_large

    Examples
    --------
    >>> Wilke([0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07])
    9.701614885866193e-06

    References
    ----------
    .. [1] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." The Journal of
       Chemical Physics 18, no. 4 (April 1, 1950): 517-19.
       https://doi.org/10.1063/1.1747673.
    '''
    cmps = range(len(ys))
    phis = [[(1.0 + (mus[i]/mus[j])**0.5*(MWs[j]/MWs[i])**0.25)**2.0/(8.0*(1.0 + MWs[i]/MWs[j]))**0.5
                    for j in cmps] for i in cmps]
    # Some publications show the denominator sum should not consider i ==j and have only the
    # mole fraction  but this reduces to that as phi[i][i] == 1
    return sum([ys[i]*mus[i]/sum([ys[j]*phis[i][j] for j in cmps]) for i in cmps])


def Wilke_prefactors(MWs):
    r'''The :obj:`Wilke` gas viscosity method can be sped up by precomputing several
    matrices. The memory used is proportional to N^2, so it can be significant,
    but is still a substantial performance increase even when they are so large
    they cannot fit into cached memory. These matrices are functions of
    molecular weights only. These are used by the :obj:`Wilke_prefactored` function.

    .. math::
        t0_{i,j} = \frac{ \sqrt{\frac{MW_{j}}{MW_{i}}}}{
        \sqrt{\frac{8 MW_{i}}{MW_{j}} + 8}}

    .. math::
        t1_{i,j} = \frac{2 \sqrt[4]{\frac{MW_{j}}{MW_{i}}}
        }{\sqrt{\frac{8 MW_{i}}{MW_{j}} + 8}}

    .. math::
        t2_{i,j} = \frac{1}{\sqrt{\frac{8 MW_{i}}{MW_{j}} + 8}}

    Parameters
    ----------
    MWs : list[float]
        Molecular weights of all components, [g/mol]

    Returns
    -------
    t0s : list[list[float]]
        First terms, [-]
    t1s : list[list[float]]
        Second terms, [-]
    t2s : list[list[float]]
        Third terms, [-]

    Notes
    -----
    These terms are derived as follows using SymPy. The viscosity terms are not
    known before hand so they are not included in the factors, but otherwise
    these parameters simplify the computation of the :math:`\phi_{ij}` term
    to the following:

    .. math::
        \phi_{ij} = \frac{\mu_i}{\mu_j}t0_{i,j} + \sqrt{\frac{\mu_i}{\mu_j}}t1_{i,j} + t2_{i,j}

    >>> from sympy import * # doctest: +SKIP
    >>> MWi, MWj, mui, muj = symbols('MW_i, MW_j, mu_i, mu_j') # doctest: +SKIP
    >>> f = (1 + sqrt(mui/muj)*(MWj/MWi)**Rational(1,4))**2 # doctest: +SKIP
    >>> denom = sqrt(8*(1+MWi/MWj)) # doctest: +SKIP
    >>> (expand(simplify(expand(f))/denom)) # doctest: +SKIP
    mu_i*sqrt(MW_j/MW_i)/(mu_j*sqrt(8*MW_i/MW_j + 8)) + 2*(MW_j/MW_i)**(1/4)*sqrt(mu_i/mu_j)/sqrt(8*MW_i/MW_j + 8) + 1/sqrt(8*MW_i/MW_j + 8) # doctest: +SKIP

    Examples
    --------
    >>> Wilke_prefactors([64.06, 46.07])
    ([[0.25, 0.19392193320396522], [0.3179655106303118, 0.25]], [[0.5, 0.421161930934918], [0.5856226024677849, 0.5]], [[0.25, 0.22867110638055677], [0.2696470380083788, 0.25]])
    >>> Wilke_prefactored([0.05, 0.95], [1.34E-5, 9.5029E-6], *Wilke_prefactors([64.06, 46.07]))
    9.701614885866193e-06
    '''
    cmps = range(len(MWs))
    MWs_inv = [1.0/MWi for MWi in MWs]
    phi_fact_invs = [[1.0/(8.0*(1.0 + MWs[i]*MWs_inv[j]))**0.5
                    for j in cmps] for i in cmps]

    t0s = [[(MWs[j]*MWs_inv[i])**0.5*phi_fact_invs[i][j]
                    for j in cmps] for i in cmps]

    t1s = [[2.0*(MWs[j]*MWs_inv[i])**0.25*phi_fact_invs[i][j]
                    for j in cmps] for i in cmps]
    return t0s, t1s, phi_fact_invs

def Wilke_prefactored(ys, mus, t0s, t1s, t2s):
    r'''Calculates viscosity of a gas mixture according to
    mixing rules in [1]_, using precomputed parameters.

    .. math::
        \eta_{mix} = \sum_{i=1}^n \frac{y_i \eta_i}{\sum_{j=1}^n y_j \phi_{ij}}

    .. math::
        \phi_{ij} = \frac{\mu_i}{\mu_j}t0_{i,j} + \sqrt{\frac{\mu_i}{\mu_j}}
        t1_{i,j} + t2_{i,j}

    Parameters
    ----------
    ys : float
        Mole fractions of gas components, [-]
    mus : float
        Gas viscosities of all components, [Pa*s]
    t0s : list[list[float]]
        First terms, [-]
    t1s : list[list[float]]
        Second terms, [-]
    t2s : list[list[float]]
        Third terms, [-]

    Returns
    -------
    mug : float
        Viscosity of gas mixture, [Pa*s]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.

    See Also
    --------
    Wilke_prefactors
    Wilke
    Wilke_large

    Examples
    --------
    >>> Wilke_prefactored([0.05, 0.95], [1.34E-5, 9.5029E-6], *Wilke_prefactors([64.06, 46.07]))
    9.701614885866193e-06

    References
    ----------
    .. [1] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." The Journal of
       Chemical Physics 18, no. 4 (April 1, 1950): 517-19.
       https://doi.org/10.1063/1.1747673.
    '''
    N = len(ys)
    mu_root_invs = [0.0]*N
    mu_roots = [0.0]*N
    mus_inv = [0.0]*N
    for i in range(N):
        # 1/sqrt(mus)
        mu_root_invs[i] = muirtinv = 1.0/sqrt(mus[i])
        # sqrt(mus)
        mu_roots[i] = muirtinv*mus[i]
        # 1/mus
        mus_inv[i] = muirtinv*muirtinv

    mu = 0.0
    for i in range(N): # numba's p range does not help here
        tot = 0.0
        for j in range(N):
            phiij = mus[i]*mus_inv[j]*t0s[i][j] + mu_roots[i]*mu_root_invs[j]*t1s[i][j] + t2s[i][j]
            tot += ys[j]*phiij
        mu += ys[i]*mus[i]/tot
    return mu

def Wilke_large(ys, mus, MWs):
    r'''Calculates viscosity of a gas mixture according to
    mixing rules in [1]_.

    This function is a slightly faster version of :obj:`Wilke`. It achieves its
    extra speed by avoiding some checks, some powers, and by allocating less
    memory during the computation. For very large component vectors, this
    function should be called instead.

    Parameters
    ----------
    ys : float
        Mole fractions of gas components, [-]
    mus : float
        Gas viscosities of all components, [Pa*s]
    MWs : float
        Molecular weights of all components, [g/mol]

    Returns
    -------
    mug : float
        Viscosity of gas mixture, [Pa*s]

    See Also
    --------
    Wilke_prefactors
    Wilke_prefactored
    Wilke

    Examples
    --------
    >>> Wilke_large([0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07])
    9.701614885866193e-06

    References
    ----------
    .. [1] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." The Journal of
       Chemical Physics 18, no. 4 (April 1, 1950): 517-19.
       https://doi.org/10.1063/1.1747673.
    '''
    # For the cases where memory is sparse or not desired to be consumed
    N = len(MWs)

#   Compute the MW and assorted power vectors
    MW_invs = [0.0]*N
    MW_inv_mus = [0.0]*N
    mu_roots = [0.0]*N
    mus_inv_MW_roots = [0.0]*N
    mu_root_invs_MW_25s = [0.0]*N

    for i in range(N):
        MW_root = sqrt(MWs[i])
        MW_root_inv = 1.0/MW_root
        MW_25_inv = sqrt(MW_root_inv)
        mu_root_inv = 1.0/sqrt(mus[i])
        x0 = mu_root_inv*MW_root
        # Stored values
        mu_roots[i] = 2.0*mu_root_inv*mus[i]*MW_25_inv
        MW_invs[i] = 8.0*MW_root_inv*MW_root_inv
        MW_inv_mus[i] = mus[i]*MW_root_inv
        mus_inv_MW_roots[i] = mu_root_inv*x0
        mu_root_invs_MW_25s[i] = x0*MW_25_inv

    mu = 0.0
    for i in range(N):
        # numba's p range does help here but only when large, when small it hinders
        tot = 0.0
        MWi = MWs[i]
        MWs_root_invi = MW_inv_mus[i]
        MW_25_invi = mu_roots[i]
        # Not a symmetric matrix unfortunately
        for j in range(N):
            # sqrt call is important for PyPy to make this fast
            # Numba sees as 25% performance increase by making this an pow(x, -0.5)
            phii_denom = ys[j]/sqrt(8.0 + MWi*MW_invs[j])
            tot += phii_denom + phii_denom*(mus_inv_MW_roots[j]*MWs_root_invi
                                + mu_root_invs_MW_25s[j]*MW_25_invi)
        mu += ys[i]*mus[i]/tot
    return mu



def Brokaw(T, ys, mus, MWs, molecular_diameters, Stockmayers):
    r'''Calculates viscosity of a gas mixture according to
    mixing rules in [1]_.

    .. math::
        \eta_{mix} = \sum_{i=1}^n \frac{y_i \eta_i}{\sum_{j=1}^n y_j \phi_{ij}}

    .. math::
        \phi_{ij} = \left( \frac{\eta_i}{\eta_j} \right)^{0.5} S_{ij} A_{ij}

    .. math::
        A_{ij} = m_{ij} M_{ij}^{-0.5} \left[1 +
        \frac{M_{ij} - M_{ij}^{0.45}}
        {2(1+M_{ij}) + \frac{(1 + M_{ij}^{0.45}) m_{ij}^{-0.5}}{1 + m_{ij}}} \right]

    .. math::
        m_{ij} = \left[ \frac{4}{(1+M_{ij}^{-1})(1+M_{ij})}\right]^{0.25}

    .. math::
        M_{ij} = \frac{M_i}{M_j}

    .. math::
        S_{ij} = \frac{1 + (T_i^* T_j^*)^{0.5} + (\delta_i \delta_j/4)}
        {[1+T_i^* + (\delta_i^2/4)]^{0.5}[1+T_j^*+(\delta_j^2/4)]^{0.5}}

    .. math::
        T^* = kT/\epsilon

    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    ys : float
        Mole fractions of gas components, [-]
    mus : float
        Gas viscosities of all components, [Pa*s]
    MWs : float
        Molecular weights of all components, [g/mol]
    molecular_diameters : float
        L-J molecular diameter  of all components, [angstroms]
    Stockmayers : float
        L-J Stockmayer energy parameters of all components, []

    Returns
    -------
    mug : float
        Viscosity of gas mixture, [Pa*s]

    Notes
    -----
    This equation is entirely dimensionless; all dimensions cancel.
    The original source has not been reviewed.

    This is DIPPR Procedure 8D: Method for the Viscosity of Nonhydrocarbon
    Vapor Mixtures at Low Pressure (Polar and Nonpolar)

    Examples
    --------
    >>> Brokaw(308.2, [0.05, 0.95], [1.34E-5, 9.5029E-6], [64.06, 46.07], [0.42, 0.19], [347, 432])
    9.699085099801568e-06

    References
    ----------
    .. [1] Brokaw, R. S. "Predicting Transport Properties of Dilute Gases."
       Industrial & Engineering Chemistry Process Design and Development
       8, no. 2 (April 1, 1969): 240-53. doi:10.1021/i260030a015.
    .. [2] Brokaw, R. S. Viscosity of Gas Mixtures, NASA-TN-D-4496, 1968.
    .. [3] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    N = len(ys)
    cmps = range(len(ys))
    MDs = molecular_diameters
    Tsts = [T/Stockmayer_i for Stockmayer_i in Stockmayers]
    Tstrs = [i**0.5 for i in Tsts]
    Aij = [[0.0]*N for j in cmps]
    phiij =[[0.0]*N for j in cmps]

    for i in cmps:
        for j in cmps:
            Sij = (1.0 +Tstrs[i]*Tstrs[j] + (MDs[i]*MDs[j])/4.)/(
                    1.0 + Tsts[i] + (0.25*MDs[i]*MDs[i]))**0.5/(1.0 + Tsts[j]
                    + (0.25*MDs[j]*MDs[j]))**0.5
            if MDs[i] <= 0.1 and MDs[j] <= 0.1:
                Sij = 1.0
            Mij = MWs[i]/MWs[j]
            Mij45 = Mij**0.45

            mij = (4./((1.0 + 1.0/Mij)*(1.0 + Mij)))**0.25

            Aij[i][j] = mij*Mij**-0.5*(1.0 + (Mij - Mij45)/(2.0*(1.0 + Mij)
                + (1.0 + Mij45)*mij**-0.5/(1.0 + mij)))

            phiij[i][j] = (mus[i]/mus[j])**0.5*Sij*Aij[i][j]

    return sum([ys[i]*mus[i]/sum([ys[j]*phiij[i][j] for j in cmps]) for i in cmps])

### Petroleum liquids

def Twu_1985_internal(T, Tb, SG):
    Tb2 = Tb*Tb
    Tb10 = Tb2*Tb2
    Tb10 *= Tb10*Tb2 # compute Tb^-10
    Tb_inv = 1.0/Tb
    Tb_sqrt_inv = 1.0/sqrt(Tb)

    # equation 15
    Tc0 = Tb/(0.533272 + 0.191017e-3*Tb + 0.779681e-7*Tb2 - 0.284376e-10*Tb2*Tb
             + 0.959468e28/(Tb10*Tb2*Tb))
    alpha = 1.0 - Tb/Tc0
    alpha3 = alpha*alpha*alpha

    SG0 = 0.843593  -0.128624*alpha - 3.36159*alpha3
    alpha6 = alpha3*alpha3
    SG0 -= 13749.5*alpha6*alpha6
    dSG = SG - SG0
    nu20 = (exp(4.73227 - 27.0975*alpha + alpha*(49.4491*alpha
             - 50.4706*alpha3)) - 1.5)


    nu10 = exp(0.801621 + 1.37179*log(nu20))

    x = abs(1.99873 - 56.7394*Tb_sqrt_inv)
    f1 = 1.33932*x*dSG - 21.1141*dSG*dSG*Tb_sqrt_inv
    f2 = x*dSG - 21.1141*dSG*dSG*Tb_sqrt_inv

    square_term2 = (1.0 + f2 + f2)/(1.0 - f2 - f2)
    square_term2 *= square_term2

    square_term1 = (1.0 + f1 + f1)/(1.0 - f1 - f1)
    square_term1 *= square_term1

    x0 = 450.0*Tb_inv
    nu1 = exp(log(nu10 + x0)*square_term1) - x0
    nu2 = exp(log(nu20 + x0)*square_term2) - x0

    # T1 = 559.67 # 100 deg F
    # T2 = 669.67 # 210 deg F
    logT1 = 6.3273473243178415 # log(559.67)
    # logT2 = 6.506785053735233 # log(669.67)

    Z1 = nu1 + 0.7 + exp(-1.47 - nu1*(1.84 + 0.51*nu1))
    Z2 = nu2 + 0.7 + exp(-1.47 - nu2*(1.84 + 0.51*nu2))

    loglogZ1 = log(log(Z1))
    try:
        B = (loglogZ1 - log(log(Z2)))*-5.572963964974682 #/(logT1 - logT2)
    except:
        B = 0.0
    try:
        Z = exp(exp(loglogZ1 + B*(log(T) - logT1)))
    except:
        Z = 1.0

    # cSt
    x0 = Z - 0.7
    nu = x0 - exp(-0.7487 + x0*(x0*(0.6119 - 0.3193*x0) - 3.295))
    return nu

def Twu_1985(T, Tb, rho):
    r'''Calculate the viscosity of a petroleum liquid using the
    Twu (1985) correlation
    developed in [1]_. Based on a fit to n-alkanes that used as a
    reference. Requires the boiling point and density of
    the system.

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tb : float
        Normal boiling point, [K]
    rho : float
        Liquid density liquid as measured at 60 deg F, [kg/m^3]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----
    The formulas are as follows:

    .. math::
        T_{c}^{\circ}=T_{b}\left(0.533272+0.191017 \times 10^{-3} T_{b}
        +0.779681 \times 10^{-7} T_{b}^{2}
        -0.284376 \times 10^{-10} T_{b}^{3}+0.959468
        \times 10^{28}/T_{b}^{13}\right)^{-1}

    .. math::
        \alpha=1-T_{b} / T_{c}^{\circ}

    .. math::
        \ln \left(\nu_2^{\circ}+1.5\right)=4.73227-27.0975 \alpha
        +49.4491 \alpha^{2}-50.4706 \alpha^{4}

    .. math::
        \ln \left(\nu_1^{\circ}\right)=0.801621+1.37179 \ln \left(\nu_2^{\circ}\right)

    .. math::
        {SG}^{\circ}=0.843593-0.128624 \alpha-3.36159 \alpha^{3}-13749.5 \alpha^{12}

    .. math::
        \Delta {SG} = {SG} - {SG}^\circ

    .. math::
        |x|=\left|1.99873-56.7394 / \sqrt{T_{b}}\right|

    .. math::
        f_{1}=1.33932|x| \Delta {SG} - 21.1141 \Delta {SG}^{2} / \sqrt{T_{b}}

    .. math::
        f_{2}=|x| \Delta {SG}-21.1141 \Delta {SG}^{2} / \sqrt{T_{b}}

    .. math::
        \ln \left(\nu_{1}+\frac{450}{T_{b}}\right)=\ln \left(\nu_{1}^{\circ}
        +\frac{450}{T_{b}}\right)\left(\frac{1+2 f_{1}}{1-2 f_{1}}\right)^{2}

    .. math::
        \ln \left(\nu_{2}+\frac{450}{T_{b}}\right)=\ln \left(\nu_{2}^{\circ}
        +\frac{450}{T_{b}}\right)\left(\frac{1+2 f_{2}}{1-2 f_{2}}\right)^{2}

    .. math::
        Z = \nu+0.7+\exp \left(-1.47-1.84 \nu-0.51 \nu^{2}\right)

    .. math::
        B=\frac{\ln \ln Z_{1}-\ln \ln Z_{2}}{\ln T_1-\ln T_2}

    .. math::
        \ln \ln Z=\ln \ln Z_{1}+B(\ln T-\ln T_1)

    .. math::
        \nu=(Z-0.7)-\exp \left(-0.7487-3.295Z-0.7)+0.6119Z-0.7)^{2}-0.3193Z-0.7)^{3}\right)


    Examples
    --------
    Sample point from article:

    >>> Twu_1985(T=338.7055, Tb=672.3166, rho=895.5189)
    0.008235009644854494

    References
    ----------
    .. [1] Twu, Chorng H. "Internally Consistent Correlation for Predicting
       Liquid Viscosities of Petroleum Fractions." Industrial & Engineering
       Chemistry Process Design and Development 24, no. 4 (October 1, 1985):
       1287-93. https://doi.org/10.1021/i200031a064.
    '''
    SG = rho*0.00100098388466972 #1/999.0170824078306
    nu = Twu_1985_internal(T*1.8, Tb*1.8, SG)
    nu = nu*1e-6 # to m^2/s
    rho = SG*999.0170824078306 # calculate density from SG
    mu = nu*rho
    return mu



### Viscosity for Liquids or Gases

def Lorentz_Bray_Clarke(T, P, Vm, zs, MWs, Tcs, Pcs, Vcs):
    r'''Calculates the viscosity of a gas or a liquid using the method of
    Lorentz, Bray, and Clarke [1]_. This method is not quite the same as the
    original, but rather the form commonly presented and used today. The
    original had a different formula for pressure correction for gases which
    was tabular and not presented entirely in [1]_. However using that
    distinction introduces a discontinuity between the liquid and gas viscosity,
    so it is not normally used.

    .. math::
        \mu [\text{centipoise}] = \mu_{\text{P low, Stiel-hThodos}} [\text{centipoise}]
        + \frac{\text{poly}^4 - 0.0001}{\xi}

    .. math::
        \text{poly} = (0.1023 + 0.023364 \rho_r + 0.058533\rho_r^2
            - 0.040758\rho_r^3 + 0.0093724\rho_r^4)

    .. math::
        \xi = T_c^{1/6} MW^{-1/2} (P_c\text{[atm]})^{-2/3}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    P : float
        Pressure of the fluid [Pa]
    Vm : float
        Molar volume of the fluid at the actual conditions, [m^3/mol]
    zs : list[float]
        Mole fractions of chemicals in the fluid, [-]
    MWs : list[float]
        Molwcular weights of chemicals in the fluid [g/mol]
    Tcs : float
        Critical temperatures of chemicals in the fluid [K]
    Pcs : float
        Critical pressures of chemicals in the fluid [Pa]
    Vcs : float
        Critical molar volumes of chemicals in the fluid; these are often used
        as tuning parameters, fit to match a pure component experimental
        viscosity value [m^3/mol]

    Returns
    -------
    mu : float
        Viscosity of phase at actual conditions , [Pa*s]

    Notes
    -----
    An example from [2]_ was implemented and checked for validation. Somewhat
    different rounding is used in [2]_.

    The mixing of the pure component Stiel-Thodos viscosities happens with the
    Herning-Zipperer mixing rule:

    .. math::
        \mu = \frac{\sum x_i \mu_i \sqrt{MW_i}}{\sum x_i \sqrt{MW_i}}

    Examples
    --------
    >>> Lorentz_Bray_Clarke(T=300.0, P=1e6, Vm=0.0023025, zs=[.4, .3, .3],
    ... MWs=[16.04246, 30.06904, 44.09562], Tcs=[190.564, 305.32, 369.83],
    ... Pcs=[4599000.0, 4872000.0, 4248000.0], Vcs=[9.86e-05, 0.0001455, 0.0002])
    9.925488160761484e-06

    References
    ----------
    .. [1] Lohrenz, John, Bruce G. Bray, and Charles R. Clark. "Calculating
       Viscosities of Reservoir Fluids From Their Compositions." Journal of
       Petroleum Technology 16, no. 10 (October 1, 1964): 1,171-1,176.
       https://doi.org/10.2118/915-PA.
    .. [2] Whitson, Curtis H., and Michael R. Brulé. Phase Behavior. Henry L.
       Doherty Memorial Fund of AIME, Society of Petroleum Engineers, 2000.
    '''
    Tc, Pc, Vc, MW = 0.0, 0.0, 0.0, 0.0
    N = len(zs)
    for i in range(N):
        Tc += Tcs[i]*zs[i]
        Pc += Pcs[i]*zs[i]
        Vc += Vcs[i]*zs[i]
        MW += MWs[i]*zs[i]

    Pc = Pc/101325. # The `xi` parameter is defined using P in atmospheres
    xi = Tc**(1.0/6.0)*MW**-0.5*Pc**(-2.0/3.0)
    rhoc = 1.0/Vc # Molar pseudocritical density
    rhom = 1.0/Vm

    rhor = rhom/rhoc
    # mu star is computed here
    mus_low_gas = [0.0]*N
    for i in range(N):
        mus_low_gas[i] = Stiel_Thodos(T, Tcs[i], Pcs[i], MWs[i])
    mu_low_gas = Herning_Zipperer(zs, mus_low_gas, MWs)

    # Polynomial - in horner form, validated
    poly = rhor*(rhor*(rhor*(0.0093724*rhor - 0.040758) + 0.058533) + 0.023364) + 0.1023

    mu_low_gas *= 1e3 # Convert low-pressure viscosity to cP

    poly2 = poly*poly
    mu = (mu_low_gas*xi + poly2*poly2 - 0.0001)/xi
    return mu*1e-3 # Convert back from cP to Pa



### Misc functions


def _round_whole_even(i):
    r'''Round a number to the nearest whole number. If the number is exactly
    between two numbers, round to the even whole number. Used by
    `viscosity_index`.

    Parameters
    ----------
    i : float
        Number, [-]

    Returns
    -------
    i : int
        Rounded number, [-]

    Notes
    -----
    Should never run with inputs from a practical function, as numbers on
    computers aren't really normally exactly between two numbers.

    Examples
    --------
    _round_whole_even(116.5)
    116
    '''
    if i % .5 == 0.0:
        if (i + 0.5) % 2 == 0.0:
            i = i + 0.5
        else:
            i = i - 0.5
    else:
        i = round(i, 0)
    return int(i)


VI_nus = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
    3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7,
    4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1,
    6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5,
    7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9,
    9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2,
    10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4,
    11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
    12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8,
    13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0,
    15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0, 16.1, 16.2,
    16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4,
    17.5, 17.6, 17.7, 17.8, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6,
    18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8,
    19.9, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2, 21.4, 21.6, 21.8, 22.0,
    22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6, 23.8, 24.0, 24.2, 24.4,
    24.6, 24.8, 25.0, 25.2, 25.4, 25.6, 25.8, 26.0, 26.2, 26.4, 26.6, 26.8,
    27.0, 27.2, 27.4, 27.6, 27.8, 28.0, 28.2, 28.4, 28.6, 28.8, 29.0, 29.2,
    29.4, 29.6, 29.8, 30.0, 30.5, 31.0, 31.5, 32.0, 32.5, 33.0, 33.5, 34.0,
    34.5, 35.0, 35.5, 36.0, 36.5, 37.0, 37.5, 38.0, 38.5, 39.0, 39.5, 40.0,
    40.5, 41.0, 41.5, 42.0, 42.5, 43.0, 43.5, 44.0, 44.5, 45.0, 45.5, 46.0,
    46.5, 47.0, 47.5, 48.0, 48.5, 49.0, 49.5, 50.0, 50.5, 51.0, 51.5, 52.0,
    52.5, 53.0, 53.5, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 57.0, 57.5, 58.0,
    58.5, 59.0, 59.5, 60.0, 60.5, 61.0, 61.5, 62.0, 62.5, 63.0, 63.5, 64.0,
    64.5, 65.0, 65.5, 66.0, 66.5, 67.0, 67.5, 68.0, 68.5, 69.0, 69.5, 70.0
]
VI_Ls = [7.994, 8.64, 9.309, 10.0, 10.71, 11.45, 12.21, 13.0, 13.8, 14.63,
    15.49, 16.36, 17.26, 18.18, 19.12, 20.09, 21.08, 22.09, 23.13, 24.19,
    25.32, 26.5, 27.75, 29.07, 30.48, 31.96, 33.52, 35.13, 36.79, 38.5,
    40.23, 41.99, 43.76, 45.53, 47.31, 49.09, 50.87, 52.64, 54.42, 56.2,
    57.97, 59.74, 61.52, 63.32, 65.18, 67.12, 69.16, 71.29, 73.48, 75.72,
    78.0, 80.25, 82.39, 84.53, 86.66, 88.85, 91.04, 93.2, 95.43, 97.72,
    100.0, 102.3, 104.6, 106.9, 109.2, 111.5, 113.9, 116.2, 118.5, 120.9,
    123.3, 125.7, 128.0, 130.4, 132.8, 135.3, 137.7, 140.1, 142.7, 145.2,
    147.7, 150.3, 152.9, 155.4, 158.0, 160.6, 163.2, 165.8, 168.5, 171.2,
    173.9, 176.6, 179.4, 182.1, 184.9, 187.6, 190.4, 193.3, 196.2, 199.0,
    201.9, 204.8, 207.8, 210.7, 213.6, 216.6, 219.6, 222.6, 225.7, 228.8,
    231.9, 235.0, 238.1, 241.2, 244.3, 247.4, 250.6, 253.8, 257.0, 260.1,
    263.3, 266.6, 269.8, 273.0, 276.3, 279.6, 283.0, 286.4, 289.7, 293.0,
    296.5, 300.0, 303.4, 306.9, 310.3, 313.9, 317.5, 321.1, 324.6, 328.3,
    331.9, 335.5, 339.2, 342.9, 346.6, 350.3, 354.1, 358.0, 361.7, 365.6,
    369.4, 373.3, 377.1, 381.0, 384.9, 388.9, 392.7, 396.7, 400.7, 404.6,
    408.6, 412.6, 416.7, 420.7, 424.9, 429.0, 433.2, 437.3, 441.5, 445.7,
    449.9, 454.2, 458.4, 462.7, 467.0, 471.3, 475.7, 479.7, 483.9, 488.6,
    493.2, 501.5, 510.8, 519.9, 528.8, 538.4, 547.5, 556.7, 566.4, 575.6,
    585.2, 595.0, 604.3, 614.2, 624.1, 633.6, 643.4, 653.8, 663.3, 673.7,
    683.9, 694.5, 704.2, 714.9, 725.7, 736.5, 747.2, 758.2, 769.3, 779.7,
    790.4, 801.6, 812.8, 824.1, 835.5, 847.0, 857.5, 869.0, 880.6, 892.3,
    904.1, 915.8, 927.6, 938.6, 951.2, 963.4, 975.4, 987.1, 998.9, 1011.0,
    1023.0, 1055.0, 1086.0, 1119.0, 1151.0, 1184.0, 1217.0, 1251.0, 1286.0,
    1321.0, 1356.0, 1391.0, 1427.0, 1464.0, 1501.0, 1538.0, 1575.0, 1613.0,
    1651.0, 1691.0, 1730.0, 1770.0, 1810.0, 1851.0, 1892.0, 1935.0, 1978.0,
    2021.0, 2064.0, 2108.0, 2152.0, 2197.0, 2243.0, 2288.0, 2333.0, 2380.0,
    2426.0, 2473.0, 2521.0, 2570.0, 2618.0, 2667.0, 2717.0, 2767.0, 2817.0,
    2867.0, 2918.0, 2969.0, 3020.0, 3073.0, 3126.0, 3180.0, 3233.0, 3286.0,
    3340.0, 3396.0, 3452.0, 3507.0, 3563.0, 3619.0, 3676.0, 3734.0, 3792.0,
    3850.0, 3908.0, 3966.0, 4026.0, 4087.0, 4147.0, 4207.0, 4268.0, 4329.0,
    4392.0, 4455.0, 4517.0, 4580.0, 4645.0, 4709.0, 4773.0, 4839.0, 4905.0
]

VI_Hs = [6.394, 6.894, 7.41, 7.944, 8.496, 9.063, 9.647, 10.25, 10.87, 11.5,
    12.15, 12.82, 13.51, 14.21, 14.93, 15.66, 16.42, 17.19, 17.97, 18.77,
    19.56, 20.37, 21.21, 22.05, 22.92, 23.81, 24.71, 25.63, 26.57, 27.53,
    28.49, 29.46, 30.43, 31.4, 32.37, 33.34, 34.32, 35.29, 36.26, 37.23,
    38.19, 39.17, 40.15, 41.13, 42.14, 43.18, 44.24, 45.33, 46.44, 47.51,
    48.57, 49.61, 50.69, 51.78, 52.88, 53.98, 55.09, 56.2, 57.31, 58.45,
    59.6, 60.74, 61.89, 63.05, 64.18, 65.32, 66.48, 67.64, 68.79, 69.94,
    71.1, 72.27, 73.42, 74.57, 75.73, 76.91, 78.08, 79.27, 80.46, 81.67,
    82.87, 84.08, 85.3, 86.51, 87.72, 88.95, 90.19, 91.4, 92.65, 93.92,
    95.19, 96.45, 97.71, 98.97, 100.2, 101.5, 102.8, 104.1, 105.4, 106.7,
    108.0, 109.4, 110.7, 112.0, 113.3, 114.7, 116.0, 117.4, 118.7, 120.1,
    121.5, 122.9, 124.2, 125.6, 127.0, 128.4, 129.8, 131.2, 132.6, 134.0,
    135.4, 136.8, 138.2, 139.6, 141.0, 142.4, 143.9, 145.3, 146.8, 148.2,
    149.7, 151.2, 152.6, 154.1, 155.6, 157.0, 158.6, 160.1, 161.6, 163.1,
    164.6, 166.1, 167.7, 169.2, 170.7, 172.3, 173.8, 175.4, 177.0, 178.6,
    180.2, 181.7, 183.3, 184.9, 186.5, 188.1, 189.7, 191.3, 192.9, 194.6,
    196.2, 197.8, 199.4, 201.0, 202.6, 204.3, 205.9, 207.6, 209.3, 211.0,
    212.7, 214.4, 216.1, 217.7, 219.4, 221.1, 222.8, 224.5, 226.2, 227.7,
    229.5, 233.0, 236.4, 240.1, 243.5, 247.1, 250.7, 254.2, 257.8, 261.5,
    264.9, 268.6, 272.3, 275.8, 279.6, 283.3, 286.8, 290.5, 294.4, 297.9,
    301.8, 305.6, 309.4, 313.0, 317.0, 320.9, 324.9, 328.8, 332.7, 336.7,
    340.5, 344.4, 348.4, 352.3, 356.4, 360.5, 364.6, 368.3, 372.3, 376.4,
    380.6, 384.6, 388.8, 393.0, 396.6, 401.1, 405.3, 409.5, 413.5, 417.6,
    421.7, 432.4, 443.2, 454.0, 464.9, 475.9, 487.0, 498.1, 509.6, 521.1,
    532.5, 544.0, 555.6, 567.1, 579.3, 591.3, 603.1, 615.0, 627.1, 639.2,
    651.8, 664.2, 676.6, 689.1, 701.9, 714.9, 728.2, 741.3, 754.4, 767.6,
    780.9, 794.5, 808.2, 821.9, 835.5, 849.2, 863.0, 876.9, 890.9, 905.3,
    919.6, 933.6, 948.2, 962.9, 977.5, 992.1, 1007.0, 1021.0, 1036.0,
    1051.0, 1066.0, 1082.0, 1097.0, 1112.0, 1127.0, 1143.0, 1159.0, 1175.0,
    1190.0, 1206.0, 1222.0, 1238.0, 1254.0, 1270.0, 1286.0, 1303.0, 1319.0,
    1336.0, 1352.0, 1369.0, 1386.0, 1402.0, 1419.0, 1436.0, 1454.0, 1471.0,
    1488.0, 1506.0, 1523.0, 1541.0, 1558.0
]

def viscosity_index(nu_40, nu_100, rounding=False):
    r'''Calculates the viscosity index of a liquid. Requires dynamic viscosity
    of a liquid at 40°C and 100°C. Value may either be returned with or
    without rounding. Rounding is performed per the standard.

    if nu_100 < 70:

    .. math::
        L, H = \text{interp}(nu_100)

    else:

    .. math::
        L = 0.8353\nu_{100}^2 + 14.67\nu_{100} - 216

    .. math::
        H = 0.1684\nu_{100}^2 + 11.85\nu_{100} - 97

    if nu_40 > H:

    .. math::
        VI = \frac{L-nu_{40}}{L-H}\cdot 100

    else:

    .. math::
        N = \frac{\ln(H) - \ln(\nu_{40})}{\ln (\nu_{100})}

    .. math::
        VI = \frac{10^N-1}{0.00715} + 100

    Parameters
    ----------
    nu_40 : float
        Dynamic viscosity of fluid at 40°C, [m^2/s]
    nu_100 : float
        Dynamic viscosity of fluid at 100°C, [m^2/s]
    rounding : bool, optional
        Whether to round the value or not.

    Returns
    -------
    VI: float
        Viscosity index [-]

    Notes
    -----
    VI is undefined for nu_100 under 2 mm^2/s. None is returned if this is the
    case. Internal units are mm^2/s. Higher values of viscosity index suggest
    a lesser decrease in kinematic viscosity as temperature increases.

    Note that viscosity is a pressure-dependent property, and that the
    viscosity index is defined for a fluid at whatever pressure it is at.
    The viscosity index is thus also a function of pressure.

    Examples
    --------
    >>> viscosity_index(73.3E-6, 8.86E-6, rounding=True)
    92

    References
    ----------
    .. [1] ASTM D2270-10(2016) Standard Practice for Calculating Viscosity
       Index from Kinematic Viscosity at 40 °C and 100 °C, ASTM International,
       West Conshohocken, PA, 2016, http://dx.doi.org/10.1520/D2270-10R16
    '''
    nu_40, nu_100 = nu_40*1E6, nu_100*1E6  # m^2/s to mm^2/s
    if nu_100 < 2.0:
        return None  # Not defined for under this
    elif nu_100 < 70.0:
        L = interp(nu_100, VI_nus, VI_Ls)
        H = interp(nu_100, VI_nus, VI_Hs)
    else:
        L = (0.8353*nu_100 + 14.67)*nu_100 - 216.0
        H = (0.1684*nu_100 + 11.85)*nu_100 - 97.0
    if nu_40 > H:
        VI = (L-nu_40)/(L-H)*100.0
    else:
        N = (log(H) - log(nu_40))/log(nu_100)
        VI = (10**N - 1.0)/0.00715 + 100.0
    if rounding:
        VI = _round_whole_even(VI)
    return VI





# All results in units of seconds, except engler and barbey which are degrees
# Data from Hydraulic Institute Handbook

viscosity_scales = {}

SSU_SSU = [31.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 15000.0, 20000.0]
SSU_nu = [1.0, 2.56, 4.3, 7.4, 10.3, 13.1, 15.7, 18.2, 20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['saybolt universal'] = (SSU_SSU, SSU_nu)

SSF_SSF = [12.95, 13.7, 14.44, 15.24, 19.3, 23.5, 28.0, 32.5, 41.9, 51.6, 61.4, 71.1, 81.0, 91.0, 100.7, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1500.0, 2000.0]
SSF_nu = [13.1, 15.7, 18.2, 20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['saybolt furol'] = (SSF_SSF, SSF_nu)

SRS_SRS = [29.0, 32.1, 36.2, 44.3, 52.3, 60.9, 69.2, 77.6, 85.6, 128.0, 170.0, 212.0, 254.0, 338.0, 423.0, 508.0, 592.0, 677.0, 762.0, 896.0, 1270.0, 1690.0, 2120.0, 2540.0, 3380.0, 4230.0, 5080.0, 5920.0, 6770.0, 7620.0, 8460.0, 13700.0, 18400.0]
SRS_nu = SSU_nu
viscosity_scales['redwood standard'] = (SRS_SRS, SRS_nu)

SRA_SRA = [5.1, 5.83, 6.77, 7.6, 8.44, 9.3, 10.12, 14.48, 18.9, 23.45, 28.0, 37.1, 46.2, 55.4, 64.6, 73.8, 83.0, 92.1, 138.2, 184.2, 230.0, 276.0, 368.0, 461.0, 553.0, 645.0, 737.0, 829.0, 921.0]
SRA_nu = [4.3, 7.4, 10.3, 13.1, 15.7, 18.2, 20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0]
viscosity_scales['redwood admiralty'] = (SRA_SRA, SRA_nu)

Engler_degrees = [1.0, 1.16, 1.31, 1.58, 1.88, 2.17, 2.45, 2.73, 3.02, 4.48, 5.92, 7.35, 8.79, 11.7, 14.6, 17.5, 20.45, 23.35, 26.3, 29.2, 43.8, 58.4, 73.0, 87.6, 117.0, 146.0, 175.0, 204.5, 233.5, 263.0, 292.0, 438.0, 584.0]
Engler_nu = SSU_nu
viscosity_scales['engler'] = (Engler_degrees, Engler_nu)

# Note: Barbey is decreasing not increasing
Barbey_degrees = [6200.0, 2420.0, 1440.0, 838.0, 618.0, 483.0, 404.0, 348.0, 307.0, 195.0, 144.0, 114.0, 95.0, 70.8, 56.4, 47.0, 40.3, 35.2, 31.3, 28.2, 18.7, 14.1, 11.3, 9.4, 7.05, 5.64, 4.7, 4.03, 3.52, 3.13, 2.82, 2.5, 1.4]
Barbey_nu = SSU_nu
viscosity_scales['barbey'] = (Barbey_degrees, Barbey_nu)
#
PC7_PC7 = [40.0, 46.0, 52.5, 66.0, 79.0, 92.0, 106.0, 120.0, 135.0, 149.0]
PC7_nu = [43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0]
viscosity_scales['parlin cup #7'] = (PC7_PC7, PC7_nu)

PC10_PC10 = [15.0, 21.0, 25.0, 30.0, 35.0, 39.0, 41.0, 43.0, 65.0, 86.0, 108.0, 129.0, 172.0, 215.0, 258.0, 300.0, 344.0, 387.0, 430.0, 650.0, 860.0]
PC10_nu = [65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['parlin cup #10'] = (PC10_PC10, PC10_nu)

PC15_PC15 = [6.0, 7.2, 7.8, 8.5, 9.0, 9.8, 10.7, 11.5, 15.2, 19.5, 24.0, 28.5, 37.0, 47.0, 57.0, 67.0, 76.0, 86.0, 96.0, 147.0, 203.0]
PC15_nu = PC10_nu
viscosity_scales['parlin cup #15'] = (PC15_PC15, PC15_nu)

PC20_PC20 = [3.0, 3.2, 3.4, 3.6, 3.9, 4.1, 4.3, 4.5, 6.3, 7.5, 9.0, 11.0, 14.0, 18.0, 22.0, 25.0, 29.0, 32.0, 35.0, 53.0, 70.0]
PC20_nu = PC10_nu
viscosity_scales['parlin cup #20'] = (PC20_PC20, PC20_nu)

FC3_FC3 = [30.0, 42.0, 50.0, 58.0, 67.0, 74.0, 82.0, 90.0, 132.0, 172.0, 218.0, 258.0, 337.0, 425.0, 520.0, 600.0, 680.0, 780.0, 850.0, 1280.0, 1715.0]
FC3_nu = PC10_nu
viscosity_scales['ford cup #3'] = (FC3_FC3, FC3_nu)

FC4_FC4 = [20.0, 28.0, 34.0, 40.0, 45.0, 50.0, 57.0, 62.0, 90.0, 118.0, 147.0, 172.0, 230.0, 290.0, 350.0, 410.0, 465.0, 520.0, 575.0, 860.0, 1150.0]
FC4_nu = PC10_nu
viscosity_scales['ford cup #4'] = (FC4_FC4, FC4_nu)

MM_MM = [125.0, 145.0, 165.0, 198.0, 225.0, 270.0, 320.0, 370.0, 420.0, 470.0, 515.0, 570.0, 805.0, 1070.0, 1325.0, 1690.0, 2110.0, 2635.0, 3145.0, 3670.0, 4170.0, 4700.0, 5220.0, 7720.0, 10500.0]
MM_nu = [20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['mac michael'] = (MM_MM, MM_nu)

ZC1_ZC1 = [38.0, 47.0, 54.0, 62.0, 73.0, 90.0]
ZC1_nu = [20.6, 32.1, 43.2, 54.0, 65.0, 87.6]
viscosity_scales['zahn cup #1'] = (ZC1_ZC1, ZC1_nu)

ZC2_ZC2 = [18.0, 20.0, 23.0, 26.0, 29.0, 37.0, 46.0, 55.0, 63.0, 72.0, 80.0, 88.0]
ZC2_nu = [20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0]
viscosity_scales['zahn cup #2'] = (ZC2_ZC2, ZC2_nu)

ZC3_ZC3 = [22.5, 24.5, 27.0, 29.0, 40.0, 51.0, 63.0, 75.0]
ZC3_nu = [154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0]
viscosity_scales['zahn cup #3'] = (ZC3_ZC3, ZC3_nu)

ZC4_ZC4 = [18.0, 20.0, 28.0, 34.0, 41.0, 48.0, 63.0, 77.0]
ZC4_nu = [198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0]
viscosity_scales['zahn cup #4'] = (ZC4_ZC4, ZC4_nu)

ZC5_ZC5 = [13.0, 18.0, 24.0, 29.0, 33.0, 43.0, 50.0, 65.0, 75.0, 86.0, 96.0]
ZC5_nu = [220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0]
viscosity_scales['zahn cup #5'] = (ZC5_ZC5, ZC5_nu)

D1_D1 = [1.3, 2.3, 3.2, 4.1, 4.9, 5.7, 6.5, 10.0, 13.5, 16.9, 20.4, 27.4, 34.5, 41.0, 48.0, 55.0, 62.0, 69.0, 103.0, 137.0, 172.0, 206.0, 275.0, 344.0, 413.0, 481.0, 550.0, 620.0, 690.0, 1030.0, 1370.0]
D1_nu = [4.3, 7.4, 10.3, 13.1, 15.7, 18.2, 20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['demmier #1'] = (D1_D1, D1_nu)

D10_D10 = [1.0, 1.4, 1.7, 2.0, 2.7, 3.5, 4.1, 4.8, 5.5, 6.2, 6.9, 10.3, 13.7, 17.2, 20.6, 27.5, 34.4, 41.3, 48.0, 55.0, 62.0, 69.0, 103.0, 137.0]
D10_nu = [32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['demmier #10'] = (D10_D10, D10_nu)

S100_S100 = [2.6, 3.6, 4.6, 5.5, 6.4, 7.3, 11.3, 15.2, 19.0, 23.0, 31.0, 39.0, 46.0, 54.0, 62.0, 70.0, 77.0, 116.0, 154.0, 193.0, 232.0, 308.0, 385.0, 462.0, 540.0, 618.0, 695.0, 770.0, 1160.0, 1540.0]
S100_nu = [7.4, 10.3, 13.1, 15.7, 18.2, 20.6, 32.1, 43.2, 54.0, 65.0, 87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['stormer 100g load'] = (S100_S100, S100_nu)

PLF_PLF = [7.0, 8.0, 9.0, 9.5, 10.8, 11.9, 12.4, 16.8, 22.0, 27.6, 33.7, 45.0, 55.8, 65.5, 77.0, 89.0, 102.0, 113.0, 172.0, 234.0]
PLF_nu = [87.6, 110.0, 132.0, 154.0, 176.0, 198.0, 220.0, 330.0, 440.0, 550.0, 660.0, 880.0, 1100.0, 1320.0, 1540.0, 1760.0, 1980.0, 2200.0, 3300.0, 4400.0]
viscosity_scales['pratt lambert f'] = (PLF_PLF, PLF_nu)

viscosity_scales['kinematic viscosity'] = (SSU_nu, SSU_nu)


viscosity_converters_to_nu = {}
viscosity_converters_from_nu = {}
viscosity_converter_limits = {}

_created_viscosity_converters = False
def _create_viscosity_converters():
    global _created_viscosity_converters
    from scipy.interpolate import UnivariateSpline


    for key, val in viscosity_scales.items():
        if key == 'barbey':
            continue
        values, nus = val
        viscosity_converter_limits[key] = (values[0], values[-1], nus[0], nus[-1])
        values, nus = np.log(values), np.log(nus)
        viscosity_converters_to_nu[key] = UnivariateSpline(values, nus, k=3, s=0)
        viscosity_converters_from_nu[key] = UnivariateSpline(nus, values, k=3, s=0)

    # Barbey gets special treatment because of its reversed values
    viscosity_converter_limits['barbey'] = (Barbey_degrees[-1], Barbey_degrees[0], Barbey_nu[0], Barbey_nu[-1])
    barbey_values, barbey_nus = np.log(list(reversed(Barbey_degrees))), np.log(list(reversed(Barbey_nu)))
    viscosity_converters_to_nu['barbey'] = UnivariateSpline(barbey_values, barbey_nus, k=3, s=0)
    viscosity_converters_from_nu['barbey'] = UnivariateSpline(np.log(Barbey_nu), np.log(Barbey_degrees), k=3, s=0)

    _created_viscosity_converters = True

# originally from  Euverard, M. R., "The Efflux Type Viscosity Cup," National
# Paint, Varnish, and Lacquer Association, 9 April 1948.
# actually found in the Paint Testing Manual
# stored are (coefficient, and minimum time (seconds))
# some of these overlap with the tabulated values; those are used in preference

# Note: Engler can also be reported in units of time? Would be good to have a reference.

viscosity_scales_linear = {
    'american can': (3.5, 35),
    'astm 0.07': (1.4, 60),
    'astm 0.10': (4.8, 25),
    'astm 0.15': (21, 9),
    'astm 0.20': (61, 5),
    'astm 0.25': (140, 4),
    'a&w b': (18.5, 10),
    'a&w crucible': (11.7, 12),
    'caspers tin plate': (3.6, 39),
    'continental can': (3.3, 12),
    'crown cork and seal': (3.3, 12),
    'engler': (7.3, 18),
    'ford cup #3': (2.4, 34),
    'ford cup #4': (3.7, 23),
    'murphy varnish': (3.1, 24),
    'parlin cup #7': (1.3, 60),
    'parlin cup #10': (4.8, 21),
    'parlin cup #15': (21.5, 10),
    'parlin cup #20': (60, 5),
    'parlin cup #25': (140, 15),
    'parlin cup #30': (260, 10),
    'pratt lambert a': (0.61, 70),
    'pratt lambert b': (1.22, 60),
    'pratt lambert c': (2.43, 40),
    'pratt lambert d': (4.87, 25),
    'pratt lambert e': (9.75, 15),
    'pratt lambert f': (19.5, 9),
    'pratt lambert g': (38, 7),
    'pratt lambert h': (76, 5),
    'pratt lambert i': (152, 4),
    'redwood standard': (0.23, 320),
    'saybolt furol': (2.1, 17),
    'saybolt universal': (0.21, 70),
    'scott': (1.6, 20),
    'westinghouse': (3.4, 30),
    'zahn cup #1': (0.75, 50),
    'zahn cup #2': (3.1, 30),
    'zahn cup #3': (9.8, 25),
    'zahn cup #4': (12.5, 14),
    'zahn cup #5': (23.6, 12)
}


def Saybolt_universal_eq(nu):
    return (4.6324*nu + (1E5 + 3264.*nu)/(nu*(nu*(1.646*nu + 23.97)
                                          + 262.7) + 3930.2))


def viscosity_converter(val, old_scale, new_scale, extrapolate=False):
    r'''Converts kinematic viscosity values from different scales which have
    historically been used. Though they may not be in use much, some standards
    still specify values in these scales.

    Parameters
    ----------
    val : float
        Viscosity value in the specified scale; [m^2/s] if
        'kinematic viscosity'; [degrees] if Engler or Barbey; [s] for the other
        scales.
    old_scale : str
        String representing the scale that `val` is in originally.
    new_scale : str
        String representing the scale that `val` should be converted to.
    extrapolate : bool
        If True, a conversion will be performed even if outside the limits of
        either scale; if False, and either value is outside a limit, an
        exception will be raised.

    Returns
    -------
    result : float
        Viscosity value in the specified scale; [m^2/s] if
        'kinematic viscosity'; [degrees] if Engler or Barbey; [s] for the other
        scales

    Notes
    -----
    The valid scales for this function are any of the following:

    ['a&w b', 'a&w crucible', 'american can', 'astm 0.07', 'astm 0.10',
    'astm 0.15', 'astm 0.20', 'astm 0.25', 'barbey', 'caspers tin plate',
    'continental can', 'crown cork and seal', 'demmier #1', 'demmier #10',
    'engler', 'ford cup #3', 'ford cup #4', 'kinematic viscosity',
    'mac michael', 'murphy varnish', 'parlin cup #10', 'parlin cup #15',
    'parlin cup #20', 'parlin cup #25', 'parlin cup #30', 'parlin cup #7',
    'pratt lambert a', 'pratt lambert b', 'pratt lambert c', 'pratt lambert d',
    'pratt lambert e', 'pratt lambert f', 'pratt lambert g', 'pratt lambert h',
    'pratt lambert i', 'redwood admiralty', 'redwood standard',
    'saybolt furol', 'saybolt universal', 'scott', 'stormer 100g load',
    'westinghouse', 'zahn cup #1', 'zahn cup #2', 'zahn cup #3', 'zahn cup #4',
    'zahn cup #5']

    Some of those scales are converted linearly; the rest use tabulated data
    and splines.

    Because the conversion is performed by spline functions, a re-conversion
    of a value will not yield exactly the original value. However, it is quite
    close.

    The method 'Saybolt universal' has a special formula implemented for its
    conversion, from [4]_. It is designed for maximum backwards compatibility
    with prior experimental data. It is solved by newton's method when
    kinematic viscosity is desired as an output.

    .. math::
        SUS_{eq} = 4.6324\nu_t + \frac{[1.0 + 0.03264\nu_t]}
        {[(3930.2 + 262.7\nu_t + 23.97\nu_t^2 + 1.646\nu_t^3)\times10^{-5})]}

    Examples
    --------
    >>> viscosity_converter(8.79, 'engler', 'parlin cup #7')
    52.5
    >>> viscosity_converter(700, 'Saybolt Universal Seconds', 'kinematic viscosity')
    0.00015108914751515542

    References
    ----------
    .. [1] Hydraulic Institute. Hydraulic Institute Engineering Data Book.
       Cleveland, Ohio: Hydraulic Institute, 1990.
    .. [2] Gardner/Sward. Paint Testing Manual. Physical and Chemical
       Examination of Paints, Varnishes, Lacquers, and Colors. 13th Edition.
       ASTM, 1972.
    .. [3] Euverard, M. R., The Efflux Type Viscosity Cup. National Paint,
       Varnish, and Lacquer Association, 1948.
    .. [4] API Technical Data Book: General Properties & Characterization.
       American Petroleum Institute, 7E, 2005.
    .. [5] ASTM. Standard Practice for Conversion of Kinematic Viscosity to
       Saybolt Universal Viscosity or to Saybolt Furol Viscosity. D 2161 - 93.
    '''
    if not _created_viscosity_converters:
        _create_viscosity_converters()
    def range_check(visc, scale):
        scale_min, scale_max, nu_min, nu_max = viscosity_converter_limits[scale]

        if visc < scale_min*(1.-1E-7) or visc > scale_max*(1.+1E-7):
            raise ValueError('Viscosity conversion is outside the limits of the '
                            '%s scale; given value is %s, but the range of the '
                            'scale is from %s to %s. Set `extrapolate` to True '
                            'to perform the conversion anyway.' %(scale, visc, scale_min, scale_max))

    def range_check_linear(val, c, tmin, scale):
        if val < tmin:
            raise ValueError('Viscosity conversion is outside the limits of the '
                            '%s scale; given value is %s, but the minimum time '
                            'for this scale is %s s. Set `extrapolate` to True '
                            'to perform the conversion anyway.' %(scale, val, tmin))

    old_scale = old_scale.lower().replace('degrees', '').replace('seconds', '').strip()
    new_scale = new_scale.lower().replace('degrees', '').replace('seconds', '').strip()

    # Convert to kinematic viscosity
    if old_scale == 'kinematic viscosity':
        val = 1E6*val # convert to centistokes, the basis of the functions
    elif old_scale == 'saybolt universal':
        if not extrapolate:
            range_check(val, old_scale)
        to_solve = lambda nu: Saybolt_universal_eq(nu) - val
        val = secant(to_solve, 1)
    elif old_scale in viscosity_converters_to_nu:
        if not extrapolate:
            range_check(val, old_scale)
        val = exp(viscosity_converters_to_nu[old_scale](log(val)))
    elif old_scale in viscosity_scales_linear:
        c, tmin = viscosity_scales_linear[old_scale]
        if not extrapolate:
            range_check_linear(val, c, tmin, old_scale)
        val = c*val # convert from seconds to centistokes
    else:
        keys = sorted(set(list(viscosity_scales.keys()) + list(viscosity_scales_linear.keys())))
        raise ValueError('Scale "%s" not recognized - allowable values are any of %s.' %(old_scale, keys))

    # Convert to desired scale
    if new_scale == 'kinematic viscosity':
        val = 1E-6*val # convert to m^2/s
    elif new_scale == 'saybolt universal':
        val = Saybolt_universal_eq(val)
    elif new_scale in viscosity_converters_from_nu:
        val = exp(viscosity_converters_from_nu[new_scale](log(val)))
        if not extrapolate:
            range_check(val, new_scale)
    elif new_scale in viscosity_scales_linear:
        c, tmin = viscosity_scales_linear[new_scale]
        val = val/c # convert from centistokes to seconds
        if not extrapolate:
            range_check_linear(val, c, tmin, new_scale)
    else:
        keys = sorted(set(list(viscosity_scales.keys()) + list(viscosity_scales_linear.keys())))
        raise ValueError('Scale "%s" not recognized - allowable values are any of %s.' %(new_scale, keys))
    return float(val)


