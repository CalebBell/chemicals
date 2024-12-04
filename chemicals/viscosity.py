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
.. autofunction:: chemicals.viscosity.PPDS5
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_2
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_2_exponential
.. autofunction:: chemicals.viscosity.Viswanath_Natarajan_3
.. autofunction:: chemicals.viscosity.mu_Yaws
.. autofunction:: chemicals.viscosity.dmu_Yaws_dT
.. autofunction:: chemicals.viscosity.mu_Yaws_fitting_jacobian
.. autofunction:: chemicals.viscosity.mu_TDE

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


__all__ = ['Viswanath_Natarajan_3','Letsou_Stiel', 'Przedziecki_Sridhar', 'PPDS9', 'dPPDS9_dT',
'Viswanath_Natarajan_2', 'Viswanath_Natarajan_2_exponential', 'Lucas', 'Brokaw',
'mu_TDE',
'Yoon_Thodos', 'Stiel_Thodos', 'Lucas_gas', 'viscosity_gas_Gharagheizi', 'Herning_Zipperer',
'Wilke', 'Wilke_prefactors', 'Wilke_prefactored', 'Wilke_large', 'mu_Yaws', 'dmu_Yaws_dT', 'mu_Yaws_fitting_jacobian',
'viscosity_index', 'viscosity_converter', 'Lorentz_Bray_Clarke', 'Twu_1985', 'mu_IAPWS', 'mu_air_lemmon',
'PPDS5']

from math import acos, atan, tan

from fluids.numerics import exp, interp, log, secant, sin, sqrt, trunc_exp, implementation_optimize_tck, splev
from fluids.numerics import numpy as np

from chemicals.data_reader import data_source, register_df_source
from chemicals.utils import PY37, can_load_data, mark_numba_incompatible, os_path_join, source_path

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
@mark_numba_incompatible
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
        raise AttributeError(f"module {__name__} has no attribute {name}")
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
    >>> float(np.mean([abs((mu - mu_i*1000)/mu) for mu, mu_i in zip(mus, mu_calc)]))
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

def mu_Yaws(T, A, B, C=0.0, D=0.0):
    r'''Calculate the viscosity of a liquid using the 4-term Yaws polynomial
    form. Requires input coefficients. If the
    coefficients do not yield viscosity in Pa*s, but rather cP, remove
    log10(1000) from `A`; this is required for the coefficients in [1]_.

    .. math::
        \log_{10} \mu = A + B/T + CT + DT^2

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]
    C : float
        Coefficient, [1/K]
    D : float
        Coefficient, [1/K^2]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----

    Examples
    --------
    >>> from math import log10
    >>> mu_Yaws(300.0, -6.4406-log10(1000), 1117.6, 0.0137, -0.000015465)
    0.0010066612081

    References
    ----------
    .. [1] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. 2 edition. Amsterdam Boston: Gulf
       Professional Publishing, 2014.
    '''
    exponent = (A + B/T + T*(C + D*T))
    if exponent > 308.0:
        return 1e308
    return 10.0**exponent

def dmu_Yaws_dT(T, A, B, C=0.0, D=0.0):
    r'''Calculate the temperature derivative of the viscosity of a liquid using
    the 4-term Yaws polynomial form. Requires input coefficients.

    .. math::
        \frac{\partial \mu}{\partial T} = 10^{A + \frac{B}{T} + T \left(C
        + D T\right)} \left(- \frac{B}{T^{2}} + C + 2 D T\right)
        \log{\left(10 \right)}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]
    C : float
        Coefficient, [1/K]
    D : float
        Coefficient, [1/K^2]

    Returns
    -------
    dmu_dT : float
        First temperature derivative of liquid viscosity, [Pa*s/K]

    Notes
    -----

    Examples
    --------
    >>> dmu_Yaws_dT(300.0, -9.4406, 1117.6, 0.0137, -0.000015465)
    -1.853591586963e-05
    '''
    x0 = D*T
    B_T = B/T
    return 10.0**(A + B_T + T*(C + x0))*(-B_T/T + C + 2.0*x0)*2.302585092994046

def mu_Yaws_fitting_jacobian(Ts, A, B, C, D):
    r'''Compute and return the Jacobian of the property predicted by
    the Yaws viscosity equation with respect to all the coefficients. This is
    used in fitting parameters for chemicals.

    Parameters
    ----------
    Ts : list[float]
        Temperatures of the experimental data points, [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]
    C : float
        Coefficient, [1/K]
    D : float
        Coefficient, [1/K^2]

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
        T = Ts[i]
        r = out[i]
        x0 = 1.0/T
        x1 = 10.0**(A + B*x0 + T*(C + D*T))*2.302585092994046
        r[0] = x1
        r[1] = x0*x1
        r[2] = T*x1
        r[3] = T*T*x1
    return out

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
        x8 = x0*(x1 - 1.0)*(1.0/3.0)
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

def mu_TDE(T, A, B, C, D):
    r'''Calculate the viscosity of a liquid using the 4-term exponential
    inverse-temperature fit equation used in NIST's TDE.

    .. math::
       \mu = \exp\left[A + \frac{B}{T} + \frac{C}{T^2} + \frac{D}{T^3}\right]

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    A : float
        Coefficient, [-]
    B : float
        Coefficient, [K]
    C : float
        Coefficient, [K^2]
    D : float
        Coefficient, [K^3]

    Returns
    -------
    mu : float
        Liquid viscosity, [Pa*s]

    Notes
    -----

    Examples
    --------
    Coefficients for isooctane at 400 K, as shown in [1]_.

    >>> mu_TDE(400.0, -14.0878, 3500.26, -678132.0, 6.17706e7)
    0.0001822175281438

    References
    ----------
    .. [1] "ThermoData Engine (TDE103b V10.1) User`s Guide."
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-ViscositySatL/ViscosityL.htm.
    '''
    T_inv = 1.0/T
    expr = A + T_inv*(B + T_inv*(C + D*T_inv))
    return trunc_exp(expr)

def PPDS5(T, Tc, a0, a1, a2):
    r'''Calculate the viscosity of a low-pressure gas using the 3-term
    exponential power fit developed by the PPDS and named PPDS equation 5.

    .. math::
       \mu = \frac{a_0 T_r}{\left( 1 + a_1 T_r^{a_2}(T_r - 1) \right)^{1/6}}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    a0 : float
        Coefficient, [-]
    a1 : float
        Coefficient, [-]
    a2 : float
        Coefficient, [-]

    Returns
    -------
    mu : float
        Low pressure gas viscosity, [Pa*s]

    Notes
    -----

    Examples
    --------
    Sample coefficients for n-pentane in [1]_, at 350 K:

    >>> PPDS5(T=350.0, Tc=470.008, a0=1.08003e-5, a1=0.19583, a2=0.811897)
    8.096643275836e-06

    References
    ----------
    .. [1] "ThermoData Engine (TDE103b V10.1) User`s Guide."
       https://trc.nist.gov/TDE/Help/TDE103b/Eqns-Pure-ViscosityG/PPDS5-ViscosityGas.htm.
    '''
    Tr = T/Tc
    return a0*Tr/(1.0 + a1*(Tr - 1.0)*Tr**a2)**(1.0/6.0)

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
        Molecular weight of fluid [g/mol]
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
        Molar volume of the fluid at temperature [m^3]
    omega : float
        Acentric factor of compound
    MW : float
        Molecular weight of fluid [g/mol]

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
    Tc : float
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
        Molecular weight of fluid [g/mol]

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
        Molecular weight of fluid [g/mol]

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
    Tc : float
        Critical point of fluid [K]
    Pc : float
        Critical pressure of the fluid [Pa]
    Zc : float
        Critical compressibility of the fluid [Pa]
    MW : float
        Molecular weight of fluid [g/mol]
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
        Molecular weight of fluid [g/mol]

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
    """ # Alternate variant which may be able to be faster in parallel
    N = len(ys)
    mu_root_invs = [0.0]*N
    mu_roots = [0.0]*N
    mus_inv = [0.0]*N
    tots = [0.0]*N
    for i in range(N):
        # 1/sqrt(mus)
        mu_root_invs[i] = muirtinv = 1.0/sqrt(mus[i])
        # sqrt(mus)
        mu_roots[i] = muirtinv*mus[i]
        # 1/mus
        mus_inv[i] = muirtinv*muirtinv*ys[i]
        mu_root_invs[i] *= ys[i]

    mu = 0.0
    for i in range(N):
        tot = 0.0
        # Not a symmetric matrix unfortunately
        for j in range(N):
            tot += ys[j]*t2s[i][j]
        tots[i] += tot
    for i in range(N):
        tot1 = 0.0
        for j in range(N):
            tot1 += mus_inv[j]*t0s[i][j]
        tots[i] += tot1*mus[i]
    for i in range(N):
        tot2 = 0.0
        for j in range(N):
            tot2 += mu_root_invs[j]*t1s[i][j]
        tots[i] += tot2*mu_roots[i]
    for i in range(N):
        mu += ys[i]*mus[i]/tots[i]
    return mu
    """

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
    Tstrs = [sqrt(i) for i in Tsts]
    rt_mus = [sqrt(mu) for mu in mus]
    MW_45 = [MWi**0.45 for MWi in MWs]
    term0s = [1.0/sqrt(1.0 + Tsts[j] + 0.25*MDs[j]*MDs[j]) for j in range(N)]

    tot = 0.0
    for i in cmps:
        inner_sum = 0.0
        x1 = 1.0/(sqrt(1.0 + Tsts[i] + (0.25*MDs[i]*MDs[i])))
        for j in cmps:
            # possible to optimize 3 of the divs out but still left with 6 and 4 sqrt
            if MDs[i] <= 0.1 and MDs[j] <= 0.1:
                Sij = 1.0
            else:
                Sij = (1.0 + Tstrs[i]*Tstrs[j] + 0.25*MDs[i]*MDs[j])*x1*term0s[j]
            Mij = MWs[i]/MWs[j]
            Mij45 = MW_45[i]/MW_45[j]
            mij = sqrt(sqrt(4./((1.0 + 1.0/Mij)*(1.0 + Mij))))
            Aij = mij/sqrt(Mij)*(1.0 + (Mij - Mij45)/(2.0*(1.0 + Mij)
                + (1.0 + Mij45)/(sqrt(mij)*(1.0 + mij))))
            phiij = (rt_mus[i]/rt_mus[j])*Sij*Aij
            inner_sum += ys[j]*phiij
        tot += ys[i]*mus[i]/inner_sum
    return tot

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
        Molecular weights of chemicals in the fluid [g/mol]
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
       Index from Kinematic Viscosity at 40 °C and 100 °C, ASTM International,
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
        N = (log(H/nu_40))/log(nu_100)
        VI = (10**N - 1.0)*(1.0/0.00715) + 100.0
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

viscosity_converter_limits = {
    scale: (values[0], values[-1], nus[0], nus[-1])
    for scale, (values, nus) in viscosity_scales.items()
}
# special handling
viscosity_converter_limits['barbey'] = (Barbey_degrees[-1], Barbey_degrees[0], Barbey_nu[0], Barbey_nu[-1])

"""
def generate_viscosity_dicts():
    print("# Auto-generated viscosity conversion splines")
    # Process each scale except barbey first
    for key, val in viscosity_scales.items():
        if key == 'barbey' or key == 'kinematic viscosity':
            continue
            
        values, nus = val
        values_log, nus_log = np.log(values), np.log(nus)
        
        to_nu_tck = splrep(values_log, nus_log, k=3, s=0)
        from_nu_tck = splrep(nus_log, values_log, k=3, s=0)
        
        print(f"viscosity_converters_to_nu['{key}'] = implementation_optimize_tck([")
        print(f"    {to_nu_tck[0].tolist()!r},")
        print(f"    {to_nu_tck[1].tolist()!r},")
        print(f"    {to_nu_tck[2]}")
        print("])")
        print(f"viscosity_converters_from_nu['{key}'] = implementation_optimize_tck([")
        print(f"    {from_nu_tck[0].tolist()!r},")
        print(f"    {from_nu_tck[1].tolist()!r},")
        print(f"    {from_nu_tck[2]}")
        print("])")
        print()
    
    # Handle barbey scale separately
    barbey_values, barbey_nus = list(reversed(Barbey_degrees)), list(reversed(Barbey_nu))
    barbey_values_log, barbey_nus_log = np.log(barbey_values), np.log(barbey_nus)
    
    to_nu_tck = splrep(barbey_values_log, barbey_nus_log, k=3, s=0)
    from_nu_tck = splrep(np.log(Barbey_nu), np.log(Barbey_degrees), k=3, s=0)
    
    print("    viscosity_converters_to_nu['barbey'] = implementation_optimize_tck([")
    print(f"        {to_nu_tck[0].tolist()!r},")
    print(f"        {to_nu_tck[1].tolist()!r},")
    print(f"        {to_nu_tck[2]}")
    print("    ])")
    print("    viscosity_converters_from_nu['barbey'] = implementation_optimize_tck([")
    print(f"        {from_nu_tck[0].tolist()!r},")
    print(f"        {from_nu_tck[1].tolist()!r},")
    print(f"        {from_nu_tck[2]}")
    print("    ])")
    print()
    print("create_viscosity_converter_dicts()")
"""
viscosity_converters_to_nu = {}
viscosity_converters_from_nu = {}

viscosity_converters_to_nu['saybolt universal'] = implementation_optimize_tck([
    [3.4339872044851463, 3.4339872044851463, 3.4339872044851463, 3.4339872044851463, 3.6888794541139363, 3.912023005428146, 4.0943445622221, 4.248495242049359, 4.382026634673881, 4.499809670330265, 4.605170185988092, 5.0106352940962555, 5.298317366548036, 5.521460917862246, 5.703782474656201, 5.991464547107982, 6.214608098422191, 6.396929655216146, 6.551080335043404, 6.684611727667927, 6.802394763324311, 6.907755278982137, 7.313220387090301, 7.600902459542082, 7.824046010856292, 8.006367567650246, 8.294049640102028, 8.517193191416238, 8.699514748210191, 8.85366542803745, 8.987196820661973, 9.104979856318357, 9.210340371976184, 9.903487552536127, 9.903487552536127, 9.903487552536127, 9.903487552536127],
    [-2.575601309673699e-18, 0.8982062175689627, 1.4796172487516794, 1.9931137257692506, 2.3217449048160703, 2.5685211362061713, 2.7483500452967573, 2.8981247808410164, 3.142680014381249, 3.4299774843937145, 3.746135121917541, 3.9743336288866606, 4.209367192735221, 4.451711048011456, 4.687660784671951, 4.873250999168711, 5.030113681170718, 5.1652270656408446, 5.284127795272641, 5.4936605984622675, 5.759832137976685, 6.065261800997406, 6.2963109627307094, 6.527360000337889, 6.758409068683301, 6.9894581269751255, 7.175996723320295, 7.332664599660219, 7.467819635712368, 7.586711283688449, 7.892141527648389, 8.158310759718983, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['saybolt universal'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 1.4586150226995167, 2.0014800002101243, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [3.4339872044851454, 3.4770122406396635, 3.5499978345728866, 3.866476068856097, 4.072959316830347, 4.230716131588708, 4.372175807676851, 4.491832893375579, 4.695665316073382, 4.961810515135175, 5.271921068674815, 5.509918917935071, 5.7414955220589965, 5.967553884029658, 6.1987002440732715, 6.387697257042421, 6.544173647107132, 6.679369625020088, 6.798252346087102, 7.007791924297624, 7.273958888244225, 7.579389703241187, 7.810438663611398, 8.041487745992798, 8.272536798159697, 8.503585860208874, 8.690124455826915, 8.846792332316191, 8.98194736833642, 9.100839016319444, 9.406269260275643, 9.672438492350183, 9.903487552536127, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['saybolt furol'] = implementation_optimize_tck([
    [2.5610957881455465, 2.5610957881455465, 2.5610957881455465, 2.5610957881455465, 2.67000213346468, 2.7239235502585, 2.9601050959108397, 3.1570004211501135, 3.332204510175204, 3.481240089335692, 3.735285826928092, 3.9435216724875173, 4.117409835153096, 4.264087336809195, 4.394449154672439, 4.51085950651685, 4.612145799724517, 5.0106352940962555, 5.298317366548036, 5.521460917862246, 5.703782474656201, 5.991464547107982, 6.214608098422191, 6.396929655216146, 6.551080335043404, 6.684611727667927, 6.802394763324311, 6.907755278982137, 7.600902459542082, 7.600902459542082, 7.600902459542082, 7.600902459542082],
    [2.5726122302071066, 2.693958194199299, 2.8665941879570487, 3.161076872192932, 3.4587143215270197, 3.7649798217372896, 3.977044772380672, 4.220248931358957, 4.458502089358454, 4.69015119303059, 4.872555263761827, 5.032212011814911, 5.166268190856598, 5.282140863641715, 5.498077785250693, 5.762602127463929, 6.065184312201249, 6.296325467420656, 6.527354049755852, 6.758410509314402, 6.989457851050498, 7.17599677958435, 7.332664587690123, 7.467819638330826, 7.586711283128077, 7.892141528812376, 8.15831075906347, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['saybolt furol'] = implementation_optimize_tck([
    [2.5726122302071057, 2.5726122302071057, 2.5726122302071057, 2.5726122302071057, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [2.5610957881455465, 2.5955176376740865, 2.637849385380515, 2.7703990007554493, 2.9240008937754385, 3.1317252362242067, 3.3227410264313173, 3.509762868371335, 3.7116986808655983, 3.9275052472031002, 4.109154250219636, 4.256607525913555, 4.388775885164053, 4.507754014184948, 4.706943671155109, 4.970898952984426, 5.2768865961958396, 5.507838224130827, 5.738908948929746, 5.969950180925843, 6.201001059153068, 6.387539303303461, 6.544207251986937, 6.679362272571954, 6.798253923918287, 7.10368416605005, 7.369853400049709, 7.600902459542082, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['redwood standard'] = implementation_optimize_tck([
    [3.367295829986474, 3.367295829986474, 3.367295829986474, 3.367295829986474, 3.5890591188317256, 3.7909846770510898, 3.9569963710708773, 4.109233174715851, 4.237000862623624, 4.351567427189173, 4.449685283147696, 4.852030263919617, 5.135798437050262, 5.356586274672012, 5.537334267018537, 5.823045895483019, 6.0473721790462776, 6.230481447578482, 6.3835066348840055, 6.517671272912275, 6.635946555686647, 6.79794041297493, 7.146772179452637, 7.432483807917119, 7.659171367666058, 7.8399193600125825, 8.125630988477065, 8.349957272040324, 8.533066540572527, 8.68609172787805, 8.820256365906321, 8.938531648680692, 9.04310445260027, 9.820105943597078, 9.820105943597078, 9.820105943597078, 9.820105943597078],
    [4.282841451291465e-18, 0.947037676439903, 1.477251010205093, 1.9944132211219097, 2.3359307229445525, 2.5626774509345287, 2.751329437909557, 2.8939892962699667, 3.1570848844933703, 3.425197413541664, 3.7472394803637106, 3.9744318887640464, 4.209503949548924, 4.453538574390683, 4.687826096235816, 4.871468740026616, 5.032526465958892, 5.16127260644653, 5.317080486672467, 5.397590053158708, 5.808831075970844, 6.063381839195508, 6.294129060227128, 6.527602371787728, 6.760205843032018, 6.989473697437789, 7.174666175341072, 7.333745087287264, 7.467974952604229, 7.585333136497559, 7.9301678574581995, 8.088976161950345, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['redwood standard'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 1.4586150226995167, 2.0014800002101243, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [3.3672958299864737, 3.400434306614668, 3.460603811995803, 3.752138785880567, 3.932074643902959, 4.094246535336774, 4.226012267827194, 4.345701404190787, 4.530788286407653, 4.806653210921934, 5.109115574154829, 5.3451446199103785, 5.574777317040083, 5.798471702527616, 6.03102484318245, 6.222700954612007, 6.373865463382127, 6.519045916863379, 6.611498649709145, 6.989076379744764, 7.069058817296409, 7.415486958966173, 7.6458745591303465, 7.874729352199147, 8.103355246560241, 8.336210343280028, 8.524358493490668, 8.678746964907798, 8.814741785965277, 8.935541033869432, 9.227940102424007, 9.639786192479114, 9.820105943597078, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['redwood admiralty'] = implementation_optimize_tck([
    [1.62924053973028, 1.62924053973028, 1.62924053973028, 1.62924053973028, 1.9125010869241836, 2.028148247292285, 2.1329823086078656, 2.2300144001592104, 2.3145136638593193, 2.6727683869575705, 2.9391619220655967, 3.1548704948922883, 3.332204510175204, 3.6136169696133895, 3.832979798087693, 4.014579593753238, 4.168214410788556, 4.301358731606427, 4.418840607796598, 4.522874943261261, 4.92870191133357, 5.216022123821206, 5.438079308923196, 5.62040086571715, 5.908082938168931, 6.133398042996649, 6.315358001522335, 6.4692503167957724, 6.602587892189336, 6.825460036255307, 6.825460036255307, 6.825460036255307, 6.825460036255307],
    [1.4586150226995171, 1.99245547818402, 2.198645351736692, 2.5763836703653205, 2.752292613488804, 2.8954053090461627, 3.1625950736882227, 3.434585678625406, 3.7529639955991314, 3.9744587249425765, 4.210317953537125, 4.452063844926447, 4.689274258752756, 4.873187490282762, 5.030131949628652, 5.1653692478957085, 5.283387927421693, 5.496215927041389, 5.758987858475604, 6.0646762410300905, 6.296991371437176, 6.527051338901946, 6.760290518752237, 6.987855244018296, 7.176173074639375, 7.332648915198676, 7.502977730497019, 7.621834230268609, 7.696212639346407, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['redwood admiralty'] = implementation_optimize_tck([
    [1.4586150226995167, 1.4586150226995167, 1.4586150226995167, 1.4586150226995167, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.696212639346407, 7.696212639346407, 7.696212639346407, 7.696212639346407],
    [1.6292405397302805, 1.6319773022180983, 1.817654349717493, 2.0133354226652287, 2.1235649294247927, 2.224712463645648, 2.3841932793299736, 2.6283858383869063, 2.911412080325757, 3.1438472334063343, 3.3685671133919604, 3.5902734096354734, 3.816522814754836, 4.005469480368168, 4.161334122990485, 4.296004503218459, 4.41509117349568, 4.62094487520665, 4.890037248239102, 5.194852236817404, 5.424154014131652, 5.65583129326589, 5.88540376188756, 6.120554382333644, 6.305823619357463, 6.462414898741196, 6.632394626184034, 6.751166651395489, 6.825460036255307, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['engler'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 0.2700271372130602, 0.4574248470388755, 0.6312717768418578, 0.7747271675523681, 0.8960880245566357, 1.0043016091968684, 1.1052568313867783, 1.4996230464268938, 1.7783364488959144, 1.9947003132247452, 2.1736147116970854, 2.4595888418037104, 2.681021528714291, 2.8622008809294686, 3.017982882488811, 3.150596984114906, 3.269568939183719, 3.374168709274236, 3.7796338173824005, 4.067315889834181, 4.290459441148391, 4.472780997942346, 4.762173934797756, 4.983606621708336, 5.1647859739235145, 5.320567975482857, 5.453182077108952, 5.572154032177765, 5.676753802268282, 6.369900982828227, 6.369900982828227, 6.369900982828227, 6.369900982828227],
    [2.1594153113692496e-17, 0.7037275545363345, 1.4308578642142395, 2.0240892855650605, 2.3148847620025803, 2.566472239097515, 2.749062627900356, 2.9013846175406637, 3.1438575205043926, 3.428976409714024, 3.7461303364216247, 3.9756254444891543, 4.211084523228673, 4.451340173494784, 4.687463832818195, 4.875335812879517, 5.028124457195314, 5.166916625693337, 5.28276475320344, 5.495581312016972, 5.759354077441313, 6.065390229772868, 6.296165964372802, 6.528422737952099, 6.7562300018197465, 6.98938683539475, 7.178058496345786, 7.330681085096422, 7.469506660198741, 7.585353117283889, 7.895447960631904, 8.156588918877056, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['engler'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 1.4586150226995167, 2.0014800002101243, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [0.0, 0.059332625496393246, 0.16886365696695008, 0.4084889360836577, 0.6125425274276474, 0.7589962438515636, 0.887107192852041, 0.9957434380856023, 1.192512459576446, 1.4527460137348578, 1.7528371474791022, 1.982961773267287, 2.2102322897395914, 2.436034716546568, 2.6655403224011383, 2.851796029543637, 3.012217775025054, 3.1443639230386924, 3.2661410158203377, 3.4725380675102833, 3.7408499307017538, 4.04567463885504, 4.276997482383534, 4.507408788850714, 4.74169563074714, 4.970266628489593, 5.1542513521928495, 5.3148300272908555, 5.446944441665345, 5.568722811658543, 5.8696284711234075, 6.140571969833711, 6.369900982828227, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['parlin cup #7'] = implementation_optimize_tck([
    [3.6888794541139363, 3.6888794541139363, 3.6888794541139363, 3.6888794541139363, 3.960813169597578, 4.189654742026425, 4.3694478524670215, 4.5217885770490405, 4.663439094112067, 4.787491742782046, 5.003946305945459, 5.003946305945459, 5.003946305945459, 5.003946305945459],
    [3.7658404952500653, 3.9233280575512803, 4.160734697228943, 4.45214587299842, 4.69074301735995, 4.882817316303614, 5.029261804025008, 5.207909888410892, 5.310044699429509, 5.393627546352362, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['parlin cup #7'] = implementation_optimize_tck([
    [3.765840495250065, 3.765840495250065, 3.765840495250065, 3.765840495250065, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.393627546352362, 5.393627546352362, 5.393627546352362, 5.393627546352362],
    [3.688879454113936, 3.765900816511903, 3.92856617129852, 4.1709431993242765, 4.35623876083105, 4.509783729225598, 4.658337101103849, 4.811041442198458, 4.940771350207313, 5.003946305945459, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['parlin cup #10'] = implementation_optimize_tck([
    [2.70805020110221, 2.70805020110221, 2.70805020110221, 2.70805020110221, 3.2188758248682006, 3.4011973816621555, 3.5553480614894135, 3.6635616461296463, 3.713572066704308, 3.7612001156935624, 4.174387269895637, 4.454347296253507, 4.68213122712422, 4.859812404361672, 5.147494476813453, 5.3706380281276624, 5.552959584921617, 5.703782474656201, 5.840641657373398, 5.958424693029782, 6.063785208687608, 6.756932389247553, 6.756932389247553, 6.756932389247553, 6.756932389247553],
    [4.174387269895635, 4.177868605279104, 4.622888024731438, 4.866822485352646, 5.0302665415023196, 5.109906760441441, 5.289605560506511, 5.667364571070318, 5.7029903115726714, 6.082652956001828, 6.288016568001154, 6.533015667965044, 6.757142200848558, 6.990355208478638, 7.172880278163543, 7.3373236076608706, 7.465487531312411, 7.586875920994592, 7.892705481130982, 8.141543260362633, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['parlin cup #10'] = implementation_optimize_tck([
    [4.174387269895637, 4.174387269895637, 4.174387269895637, 4.174387269895637, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [2.7080502011022096, 3.0050746948658005, 3.084113298472601, 3.3931966213361364, 3.5516524090096575, 3.672367384479932, 3.710718293097653, 3.793657747253398, 4.161319062622642, 4.424188181146925, 4.672857873926469, 4.891207941699611, 5.127143618258207, 5.356161383763761, 5.545566574329714, 5.6944535245951435, 5.836627005131794, 5.954114617700045, 6.259152912628183, 6.542772412646688, 6.756932389247553, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['parlin cup #15'] = implementation_optimize_tck([
    [1.791759469228055, 1.791759469228055, 1.791759469228055, 1.791759469228055, 2.0541237336955462, 2.1400661634962708, 2.1972245773362196, 2.2823823856765264, 2.3702437414678603, 2.4423470353692043, 2.7212954278522306, 2.970414465569701, 3.1780538303479458, 3.349904087274605, 3.6109179126442243, 3.8501476017100584, 4.04305126783455, 4.204692619390966, 4.330733340286331, 4.454347296253507, 4.564348191467836, 5.313205979041787, 5.313205979041787, 5.313205979041787, 5.313205979041787],
    [4.174387269895636, 4.069411219997944, 4.635941726808014, 4.839118622854668, 5.082895259664981, 5.167895904788395, 5.279390094845426, 5.49457694855415, 5.807506325678317, 6.07144265735525, 6.299044565451546, 6.519201243080843, 6.786529971237586, 6.984500950950229, 7.178001184972424, 7.321785222638846, 7.476997641984404, 7.585341856635592, 7.902476308077364, 8.17692436468151, 8.389359819906355, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['parlin cup #15'] = implementation_optimize_tck([
    [4.174387269895637, 4.174387269895637, 4.174387269895637, 4.174387269895637, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [1.7917594692280552, 1.970830229617362, 1.9791186272895795, 2.1459249074862647, 2.1830257440343495, 2.2767416042446986, 2.3692390446678893, 2.5098427507224996, 2.6748536525858753, 2.9495897244648175, 3.1637250489772666, 3.387554813535431, 3.5758290532766868, 3.839225148663546, 4.030963246826631, 4.203587113976053, 4.32075904756381, 4.4510953972340355, 4.766950417856681, 5.043008453764842, 5.313205979041787, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['parlin cup #20'] = implementation_optimize_tck([
    [1.0986122886681098, 1.0986122886681098, 1.0986122886681098, 1.0986122886681098, 1.2237754316221157, 1.2809338454620642, 1.3609765531356006, 1.410986973710262, 1.4586150226995167, 1.5040773967762742, 1.840549633397487, 2.0149030205422647, 2.1972245773362196, 2.3978952727983707, 2.6390573296152584, 2.8903717578961645, 3.091042453358316, 3.2188758248682006, 3.367295829986474, 3.4657359027997265, 3.5553480614894135, 4.248495242049359, 4.248495242049359, 4.248495242049359, 4.248495242049359],
    [4.174387269895636, 4.39660473046765, 4.621951004757847, 4.929363588260073, 4.993461117294169, 5.174073014206359, 5.285955320978621, 5.632708492823754, 5.652999406401367, 6.1166972563540964, 6.328730773353977, 6.481094525964685, 6.810237508549038, 6.989742039396471, 7.140761691334212, 7.36479254067722, 7.442085614485727, 7.5894464026039214, 7.933298535525166, 8.118111141334893, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['parlin cup #20'] = implementation_optimize_tck([
    [4.174387269895637, 4.174387269895637, 4.174387269895637, 4.174387269895637, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [1.0986122886681102, 1.123675515426659, 1.19860175331427, 1.2651364479210463, 1.3662080607592082, 1.4059546827393599, 1.4584559265920132, 1.5341047757161053, 1.8539096701151676, 1.9870353450750486, 2.1668514375991434, 2.458329365130117, 2.5917386078007136, 2.874448418835185, 3.0993850182758695, 3.194743066616606, 3.374934434419905, 3.4598590976024552, 3.721050274617513, 4.054812225158822, 4.248495242049359, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['ford cup #3'] = implementation_optimize_tck([
    [3.4011973816621555, 3.4011973816621555, 3.4011973816621555, 3.4011973816621555, 3.912023005428146, 4.060443010546419, 4.204692619390966, 4.30406509320417, 4.406719247264253, 4.499809670330265, 4.882801922586371, 5.147494476813453, 5.384495062789089, 5.552959584921617, 5.820082930352362, 6.052089168924417, 6.253828811575473, 6.396929655216146, 6.522092798170152, 6.659293919683638, 6.745236349484363, 7.44716835960004, 7.44716835960004, 7.44716835960004, 7.44716835960004],
    [4.174387269895635, 4.235010951495693, 4.546306401368109, 4.891301119016985, 5.006035087987643, 5.179475622998687, 5.282927421563667, 5.5071623985516815, 5.744909093520779, 6.092188693108466, 6.27711122625197, 6.530168441864454, 6.773058874431685, 6.999063720043492, 7.158310121353729, 7.331424013637432, 7.488770336766041, 7.559595844806968, 7.976158876764658, 8.115360502604224, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['ford cup #3'] = implementation_optimize_tck([
    [4.174387269895637, 4.174387269895637, 4.174387269895637, 4.174387269895637, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [3.4011973816621555, 3.6658150644295513, 3.8217904719932014, 4.0438126549887965, 4.209411716855145, 4.293444443486365, 4.404349254525439, 4.584716740249135, 4.858012747496513, 5.1116516587850995, 5.380437068564622, 5.58277904895064, 5.794305549281043, 6.031853923159892, 6.25209862156876, 6.3925603941058675, 6.504933849634002, 6.668104522488961, 6.873672983280552, 7.249970267314259, 7.44716835960004, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['ford cup #4'] = implementation_optimize_tck([
    [2.995732273553991, 2.995732273553991, 2.995732273553991, 2.995732273553991, 3.5263605246161616, 3.6888794541139363, 3.8066624897703196, 3.912023005428146, 4.04305126783455, 4.127134385045092, 4.499809670330265, 4.770684624465665, 4.990432586778736, 5.147494476813453, 5.438079308923196, 5.66988092298052, 5.857933154483459, 6.016157159698354, 6.142037405587356, 6.253828811575473, 6.354370040797351, 7.047517221357296, 7.047517221357296, 7.047517221357296, 7.047517221357296],
    [4.174387269895637, 4.251488580015067, 4.578802297622766, 4.857630648725153, 5.03030241586516, 5.195074682944892, 5.257846026150228, 5.529297583846938, 5.752866695891904, 6.077911699309654, 6.2767900479705325, 6.554004012227638, 6.758509674974357, 6.988689456935371, 7.176894254362655, 7.324981145497996, 7.46897293490252, 7.586769337379531, 7.904669560066052, 8.158272029516167, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['ford cup #4'] = implementation_optimize_tck([
    [4.174387269895637, 4.174387269895637, 4.174387269895637, 4.174387269895637, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [2.9957322735539913, 3.2543421273633903, 3.4129326031173814, 3.6865473452369435, 3.80421925615069, 3.8923855836695638, 4.051448296413546, 4.190864134323703, 4.4719077976775665, 4.7420473255424485, 4.987515192397558, 5.16785281223535, 5.417319039066321, 5.656229034075097, 5.847207555041156, 6.012893732297337, 6.136167475065811, 6.250046666279471, 6.539886402559759, 6.816102294644018, 7.047517221357296, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['mac michael'] = implementation_optimize_tck([
    [4.8283137373023015, 4.8283137373023015, 4.8283137373023015, 4.8283137373023015, 5.10594547390058, 5.288267030694535, 5.41610040220442, 5.598421958998375, 5.768320995793772, 5.91350300563827, 6.040254711277414, 6.152732694704104, 6.244166900663736, 6.345636360828596, 6.690842277418564, 6.975413927455952, 7.1891677384203225, 7.432483807917119, 7.6544432264701125, 7.876638460975463, 8.05356916913454, 8.207946941048617, 8.335671314792847, 8.45531778769815, 8.560252680876685, 9.259130536145614, 9.259130536145614, 9.259130536145614, 9.259130536145614],
    [3.0252910757955362, 3.309048660225974, 3.8015578068795866, 3.9477102091193084, 4.197562932422006, 4.4818256053756516, 4.6899230339151625, 4.875816005626788, 5.033739813098702, 5.156529510812838, 5.300152694461719, 5.4642538160845735, 5.8047350194763325, 6.046487048467108, 6.350275147569952, 6.4365377098076495, 6.810164293656892, 6.97909139000128, 7.18182393923967, 7.327133812599172, 7.473172812291389, 7.584998019355741, 7.895273291164656, 8.194667189042956, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['mac michael'] = implementation_optimize_tck([
    [3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [4.8283137373023015, 4.977572196004898, 4.920899182300126, 5.288415669109919, 5.443641391461551, 5.572517443890477, 5.755665691788839, 5.905643458046299, 6.032351381772725, 6.153169546264002, 6.2336628658740185, 6.453386619988654, 6.627709633315948, 6.9723057914990045, 7.142439342344769, 7.5197673498103255, 7.603779301481185, 7.87108111536658, 8.04047575929093, 8.20446815597734, 8.32770260556312, 8.452105814008403, 8.755175884170438, 8.990090610140427, 9.259130536145614, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['zahn cup #1'] = implementation_optimize_tck([
    [3.6375861597263857, 3.6375861597263857, 3.6375861597263857, 3.6375861597263857, 3.9889840465642745, 4.127134385045092, 4.499809670330265, 4.499809670330265, 4.499809670330265, 4.499809670330265],
    [3.025291075795535, 3.222644075465119, 3.666419096064154, 4.127083853487198, 4.2250148808030294, 4.472780997942346, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['zahn cup #1'] = implementation_optimize_tck([
    [3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.765840495250065, 3.9889840465642745, 4.472780997942346, 4.472780997942346, 4.472780997942346, 4.472780997942346],
    [3.6375861597263857, 3.7769696282041414, 3.8873423262901667, 4.1609287874314695, 4.449826287955393, 4.499809670330265, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['zahn cup #2'] = implementation_optimize_tck([
    [2.8903717578961645, 2.8903717578961645, 2.8903717578961645, 2.8903717578961645, 3.1354942159291497, 3.258096538021482, 3.367295829986474, 3.6109179126442243, 3.828641396489095, 4.007333185232471, 4.143134726391533, 4.276666119016055, 4.477336814478207, 4.477336814478207, 4.477336814478207, 4.477336814478207],
    [3.0252910757955345, 3.4939623726886793, 3.6846822351087916, 3.983448284883306, 4.2603024753818, 4.463633546652626, 4.691320918802126, 4.859619373592127, 5.04374341823447, 5.183883618647157, 5.324994691713382, 5.393627546352363, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['zahn cup #2'] = implementation_optimize_tck([
    [3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.0252910757955354, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.393627546352362, 5.393627546352362, 5.393627546352362, 5.393627546352362],
    [2.890371757896164, 2.9019929388211776, 3.0355921974217805, 3.251781703198547, 3.3802279222783604, 3.58558219226917, 3.810689745569852, 4.0057220423981725, 4.129534879273454, 4.3146159960105335, 4.404734261587551, 4.477336814478207, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['zahn cup #3'] = implementation_optimize_tck([
    [3.1135153092103742, 3.1135153092103742, 3.1135153092103742, 3.1135153092103742, 3.295836866004329, 3.367295829986474, 3.6888794541139363, 3.9318256327243257, 4.31748811353631, 4.31748811353631, 4.31748811353631, 4.31748811353631],
    [5.0369526024136295, 5.159791792082915, 5.224416524098557, 5.527395600207853, 5.761664633167972, 6.1555615441667895, 6.351712744533946, 6.492239835020471, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['zahn cup #3'] = implementation_optimize_tck([
    [5.0369526024136295, 5.0369526024136295, 5.0369526024136295, 5.0369526024136295, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.492239835020471, 6.492239835020471, 6.492239835020471, 6.492239835020471],
    [3.1135153092103742, 3.143178103410941, 3.2706425864829805, 3.4269498134078895, 3.6603264183165023, 3.954219559147537, 4.192558405324681, 4.31748811353631, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['zahn cup #4'] = implementation_optimize_tck([
    [2.8903717578961645, 2.8903717578961645, 2.8903717578961645, 2.8903717578961645, 3.332204510175204, 3.5263605246161616, 3.713572066704308, 3.871201010907891, 4.343805421853684, 4.343805421853684, 4.343805421853684, 4.343805421853684],
    [5.288267030694534, 5.432615247755282, 5.651928263819595, 6.101421209819679, 6.295187294274101, 6.621999291047122, 6.8126881329071205, 7.003065458786462, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['zahn cup #4'] = implementation_optimize_tck([
    [5.288267030694535, 5.288267030694535, 5.288267030694535, 5.288267030694535, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 7.003065458786462, 7.003065458786462, 7.003065458786462, 7.003065458786462],
    [2.8903717578961654, 3.066256960033695, 3.2992505316778855, 3.4968229248625007, 3.703841744777038, 3.9591375071077697, 4.204162454398758, 4.343805421853684, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['zahn cup #5'] = implementation_optimize_tck([
    [2.5649493574615367, 2.5649493574615367, 2.5649493574615367, 2.5649493574615367, 3.1780538303479458, 3.367295829986474, 3.4965075614664802, 3.7612001156935624, 3.912023005428146, 4.174387269895637, 4.31748811353631, 4.564348191467836, 4.564348191467836, 4.564348191467836, 4.564348191467836],
    [5.393627546352365, 5.710312512118053, 5.925474983703362, 6.27000301923224, 6.577527820379961, 6.686902295630258, 7.105058919635149, 7.112863665222863, 7.394253732954623, 7.488520782028579, 7.590852123688581, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['zahn cup #5'] = implementation_optimize_tck([
    [5.393627546352362, 5.393627546352362, 5.393627546352362, 5.393627546352362, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.590852123688581, 7.590852123688581, 7.590852123688581, 7.590852123688581],
    [2.564949357461537, 2.6822088440741987, 3.046146719461773, 3.3677420173958623, 3.5006795420821426, 3.7907683089722517, 3.838797898224285, 4.198368376327341, 4.32566686002696, 4.502481713423904, 4.564348191467836, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['demmier #1'] = implementation_optimize_tck([
    [0.26236426446749106, 0.26236426446749106, 0.26236426446749106, 0.26236426446749106, 1.1631508098056809, 1.410986973710262, 1.589235205116581, 1.7404661748405046, 1.8718021769015913, 2.302585092994046, 2.6026896854443837, 2.8273136219290276, 3.0155349008501706, 3.3105430133940246, 3.5409593240373143, 3.713572066704308, 3.871201010907891, 4.007333185232471, 4.127134385045092, 4.23410650459726, 4.634728988229636, 4.919980925828125, 5.147494476813453, 5.327876168789581, 5.616771097666572, 5.840641657373398, 6.023447592961033, 6.175867270105761, 6.309918278226516, 6.429719478039138, 6.536691597591305, 7.222566018822171, 7.222566018822171, 7.222566018822171, 7.222566018822171],
    [1.4586150226995171, 1.71073020413693, 2.130557177301486, 2.5459071200692494, 2.746335602097296, 2.8965003266930434, 3.115192670969138, 3.4329157094051497, 3.7390508458630056, 3.9782988851980874, 4.20697132189436, 4.456141285762206, 4.674944448477957, 4.882253102916237, 5.028896980521317, 5.16527336142792, 5.2840697109140615, 5.489424075705094, 5.76014617385828, 6.0702812571485065, 6.292044976707557, 6.530315846838931, 6.757650191781952, 6.989819784532628, 7.174344986163944, 7.3339959706348585, 7.468668358430634, 7.586630009023315, 7.885398972076964, 8.16285981672486, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['demmier #1'] = implementation_optimize_tck([
    [1.4586150226995167, 1.4586150226995167, 1.4586150226995167, 1.4586150226995167, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [0.26236426446749106, 0.6010699462118122, 0.9363299696306459, 1.3947077517127915, 1.5764812081499955, 1.7306738789371188, 1.9890056748728902, 2.2459961971737625, 2.579693991216237, 2.813292095638698, 3.0562074724399007, 3.282189336952426, 3.5318885806104703, 3.6998539305408715, 3.865145982650516, 4.001868007130916, 4.122916610774629, 4.336272582432516, 4.5959904225597015, 4.89567561549399, 5.136056991938173, 5.361088907432062, 5.595850897703561, 5.8265983368689955, 6.0149611478018645, 6.168424236702356, 6.304302416540273, 6.425530993946575, 6.736467367503558, 6.989354337693505, 7.222566018822171, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['demmier #10'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 0.5306282510621704, 0.6931471805599453, 0.9932517730102834, 1.252762968495368, 1.410986973710262, 1.5686159179138452, 1.7047480922384253, 1.824549292051046, 1.9315214116032138, 2.33214389523559, 2.617395832834079, 2.8449093838194073, 3.0252910757955354, 3.3141860046725258, 3.5380565643793527, 3.720862499966987, 3.871201010907891, 4.007333185232471, 4.127134385045092, 4.23410650459726, 4.919980925828125, 4.919980925828125, 4.919980925828125, 4.919980925828125],
    [3.4688560301359708, 3.569770497469311, 3.8483584159193303, 4.229655585481154, 4.475874515937343, 4.64509586335484, 4.893176395602696, 5.027436022534686, 5.1655917774463616, 5.284001533060033, 5.4894988836778875, 5.7601244001096354, 6.070288063445428, 6.292035849980337, 6.5303467651576215, 6.757474200127011, 6.990271203658143, 7.172404210492034, 7.336919929874899, 7.467181958740508, 7.586802071296843, 7.885050970603688, 8.163055382426146, 8.389359819906355, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['demmier #10'] = implementation_optimize_tck([
    [3.4688560301359703, 3.4688560301359703, 3.4688560301359703, 3.4688560301359703, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [-2.5940108011161694e-18, 0.24439543096273444, 0.43218802462619954, 0.7238737291530682, 0.9516479228540552, 1.2578479050323754, 1.391504805145228, 1.5637871969499424, 1.6990146050637553, 1.8203891091709945, 2.033622370617629, 2.2934244853383956, 2.5930843570290696, 2.8334805868344586, 3.0584733811479787, 3.293436233953942, 3.523575372206134, 3.713605683315231, 3.862220007233949, 4.002509008119753, 4.1227764758193075, 4.434234201765468, 4.686571049547071, 4.919980925828125, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['stormer 100g load'] = implementation_optimize_tck([
    [0.9555114450274363, 0.9555114450274363, 0.9555114450274363, 0.9555114450274363, 1.5260563034950492, 1.7047480922384253, 1.8562979903656263, 1.9878743481543455, 2.424802725718295, 2.7212954278522306, 2.9444389791664403, 3.1354942159291497, 3.4339872044851463, 3.6635616461296463, 3.828641396489095, 3.9889840465642745, 4.127134385045092, 4.248495242049359, 4.343805421853684, 4.7535901911063645, 5.0369526024136295, 5.262690188904886, 5.44673737166631, 5.730099782973574, 5.953243334287785, 6.135564891081739, 6.29156913955832, 6.42648845745769, 6.543911845564792, 6.646390514847729, 7.3395376954076745, 7.3395376954076745, 7.3395376954076745, 7.3395376954076745],
    [2.0014800002101247, 2.2081982869987318, 2.4359447967016497, 2.7461554921163596, 2.896398061351864, 3.117564806356968, 3.426312248721431, 3.739877931108494, 3.9809047078493927, 4.205991082632947, 4.455558910023335, 4.6686049571745585, 4.8892765924161745, 5.027328255028008, 5.167059609893703, 5.274710399165016, 5.516870605111536, 5.747956242836144, 6.070687236497746, 6.2958893029864935, 6.523788392473561, 6.7613656433125335, 6.98895758543227, 7.177304852304456, 7.332463277016886, 7.467132482765323, 7.584968899386042, 7.901591060586778, 8.144361746774743, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['stormer 100g load'] = implementation_optimize_tck([
    [2.0014800002101243, 2.0014800002101243, 2.0014800002101243, 2.0014800002101243, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [0.955511445027437, 1.129055253782241, 1.4021249764841106, 1.692140023209894, 1.8465608724801996, 2.1046550526727743, 2.371478655042039, 2.6982483378714024, 2.9291968044688286, 3.177359581499901, 3.4050767523313756, 3.6581125927944678, 3.8116229477343673, 3.9840076714494184, 4.120026127241616, 4.248853234980905, 4.426912765186266, 4.722920627303147, 5.012396641859333, 5.249203877271702, 5.48341394310556, 5.707079426964617, 5.940135790488285, 6.125478770623342, 6.28473919227368, 6.421360662788953, 6.540644363371554, 6.834008583321535, 7.122439894929551, 7.3395376954076745, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['pratt lambert f'] = implementation_optimize_tck([
    [1.9459101490553132, 1.9459101490553132, 1.9459101490553132, 1.9459101490553132, 2.1972245773362196, 2.2512917986064953, 2.379546134130174, 2.4765384001174837, 2.517696472610991, 2.8213788864092133, 3.091042453358316, 3.3178157727231046, 3.517497837358316, 3.8066624897703196, 4.021773869387265, 4.182050142641207, 4.343805421853684, 4.48863636973214, 4.624972813284271, 4.727387818712341, 5.455321115357702, 5.455321115357702, 5.455321115357702, 5.455321115357702],
    [4.472780997942348, 4.735723623279403, 4.643157968700142, 5.139532436075622, 5.151621937826048, 5.227025496128334, 5.646901991476376, 5.756030986163222, 6.082273367390175, 6.304393363254469, 6.514223938630912, 6.75734568699049, 6.974838232501322, 7.195432554445979, 7.331926380253746, 7.475243201027783, 7.574048611098696, 7.923204111231037, 8.153345682904884, 8.389359819906353, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['pratt lambert f'] = implementation_optimize_tck([
    [4.472780997942346, 4.472780997942346, 4.472780997942346, 4.472780997942346, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [1.945910149055314, 1.970233575484605, 2.212617992316628, 2.218944817704878, 2.380488654268661, 2.4863130470989145, 2.5367199460437266, 2.785993314533494, 3.0666917243801786, 3.2986259656253227, 3.5617751404077778, 3.782999363183989, 4.017336918420929, 4.163858930396986, 4.338371609655286, 4.478013030014296, 4.626990641459872, 4.903498060635997, 5.2171271309421146, 5.455321115357702, 0.0, 0.0, 0.0, 0.0],
    3
])

viscosity_converters_to_nu['barbey'] = implementation_optimize_tck([
    [0.3364722366212129, 0.3364722366212129, 0.3364722366212129, 0.3364722366212129, 1.0367368849500223, 1.1410330045520618, 1.2584609896100056, 1.3937663759585917, 1.547562508716013, 1.7298840655099674, 1.9530276168241774, 2.2407096892759584, 2.424802725718295, 2.6461747973841225, 2.928523523860541, 3.339321977944068, 3.4436180975461075, 3.5610460826040513, 3.696351468952637, 3.8501476017100584, 4.032469158504013, 4.259859000699674, 4.553876891600541, 4.736198448394496, 4.969813299576001, 5.272999558563747, 5.726847747587197, 5.8522024797744745, 6.0014148779611505, 6.180016653652572, 6.42648845745769, 6.731018100482083, 7.272398392570047, 8.732304571033183, 8.732304571033183, 8.732304571033183, 8.732304571033183],
    [8.389359819906355, 11.016039308886484, 7.714822572010306, 7.602753323347697, 7.461128819509047, 7.335996764118376, 7.174879950521344, 6.989767099286427, 6.758533281838772, 6.525700312748302, 6.298388483905108, 6.067101933617808, 5.7517671769019465, 5.497843241516832, 5.283655849227734, 5.163953580434897, 5.031653546990693, 4.872885650123412, 4.6855823490804, 4.450958400202121, 4.213557960341426, 3.9682272779345245, 3.747010272216809, 3.417210801373589, 3.1342657975409898, 2.8932365145450674, 2.7455867304919956, 2.545778447103095, 2.3232392010369707, 1.9028183731782218, 1.1732532186088065, 0.46813938214393835, 0.0, 0.0, 0.0, 0.0, 0.0],
    3
])
viscosity_converters_from_nu['barbey'] = implementation_optimize_tck([
    [0.0, 0.0, 0.0, 0.0, 1.4586150226995167, 2.0014800002101243, 2.33214389523559, 2.5726122302071057, 2.7536607123542622, 2.9014215940827497, 3.0252910757955354, 3.4688560301359703, 3.765840495250065, 3.9889840465642745, 4.174387269895637, 4.472780997942346, 4.700480365792417, 4.882801922586371, 5.0369526024136295, 5.170483995038151, 5.288267030694535, 5.393627546352362, 5.799092654460526, 6.0867747269123065, 6.309918278226516, 6.492239835020471, 6.779921907472252, 7.003065458786462, 7.1853870155804165, 7.3395376954076745, 7.473069088032197, 7.590852123688581, 7.696212639346407, 8.389359819906353, 8.389359819906353, 8.389359819906353, 8.389359819906353],
    [8.73230457103318, 8.22931707401472, 7.596562094846855, 6.783238261660232, 6.466667813349762, 6.195757103615377, 6.014212277460382, 5.8599436146689365, 5.6197783228924045, 5.32037565213582, 4.998898870253533, 4.745014146796411, 4.518261786803321, 4.283771553970103, 4.047728697974775, 3.8591241028036105, 3.70409882799983, 3.565715518826268, 3.447530238307399, 3.241302582823786, 2.9633181827865713, 2.668317377842094, 2.439320840853108, 2.204501278713815, 1.9747300870228615, 1.7436385592434631, 1.556406462440578, 1.4019214017072539, 1.2618866366700296, 1.1488238145304066, 0.8015807047280459, 1.2114951712296043, 0.3364722366212129, 0.0, 0.0, 0.0, 0.0],
    3
])


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

def errf(T_other, initial_T, backward_calculator):
    return backward_calculator(T_other) - initial_T

def polish_conversion(T_initial, forward_calculator, backward_calculator):
    # Get initial guess from direct conversion
    initial_conversion = forward_calculator(T_initial)
    T_polished = secant(
        func=errf,
        x0=initial_conversion,
        xtol=1e-12,
        maxiter=10,
        require_eval=False,
        require_xtol=False,
        x1=initial_conversion * (1 + 1e-6),
        args=(T_initial, backward_calculator)
    )
    return T_polished


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

    The conversion uses cubic spline interpolation and is reversible to
    nearly machine precision.

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
    def range_check(visc, scale):
        scale_min, scale_max, kinematic_viscosity_min, kinematic_viscosity_max = viscosity_converter_limits[scale]

        if visc < scale_min*(1.-1E-7) or visc > scale_max*(1.+1E-7):
            raise ValueError('Viscosity conversion is outside the limits of the '
                            f'{scale} scale; given value is {visc}, but the range of the '
                            f'scale is from {scale_min} to {scale_max}. Set `extrapolate` to True '
                            'to perform the conversion anyway.')

    def range_check_linear(val, c, tmin, scale):
        if val < tmin:
            raise ValueError('Viscosity conversion is outside the limits of the '
                            f'{scale} scale; given value is {val}, but the minimum time '
                            f'for this scale is {tmin} s. Set `extrapolate` to True '
                            'to perform the conversion anyway.')

    old_scale = old_scale.lower().replace('degrees', '').replace('seconds', '').strip()
    new_scale = new_scale.lower().replace('degrees', '').replace('seconds', '').strip()

    # Convert to kinematic viscosity
    if old_scale == 'kinematic viscosity':
        kinematic_viscosity = 1E6*val # convert to centistokes, the basis of the functions
    elif old_scale == 'saybolt universal':
        if not extrapolate:
            range_check(val, old_scale)
        to_solve = lambda nu: Saybolt_universal_eq(nu) - val
        kinematic_viscosity = secant(to_solve, 1)
    elif old_scale in viscosity_converters_to_nu:
        if not extrapolate:
            range_check(val, old_scale)
        kinematic_viscosity = exp(splev(log(val), viscosity_converters_to_nu[old_scale]))
    elif old_scale in viscosity_scales_linear:
        c, tmin = viscosity_scales_linear[old_scale]
        if not extrapolate:
            range_check_linear(val, c, tmin, old_scale)
        kinematic_viscosity = c*val # convert from seconds to centistokes
    else:
        keys = sorted(set(list(viscosity_scales.keys()) + list(viscosity_scales_linear.keys())))
        raise ValueError(f'Scale "{old_scale}" not recognized - allowable values are any of {keys}.')

    # Convert to desired scale
    if new_scale == 'kinematic viscosity':
        val = 1E-6*kinematic_viscosity # convert to m^2/s
    elif new_scale == 'saybolt universal':
        val = Saybolt_universal_eq(kinematic_viscosity)
    elif new_scale in viscosity_converters_from_nu:
        forward = lambda x: exp(splev(log(x), viscosity_converters_from_nu[new_scale]))
        backward = lambda x: exp(splev(log(x), viscosity_converters_to_nu[new_scale]))
        val = polish_conversion(kinematic_viscosity, forward, backward)
        if not extrapolate:
            range_check(val, new_scale)
    elif new_scale in viscosity_scales_linear:
        c, tmin = viscosity_scales_linear[new_scale]
        val = kinematic_viscosity/c # convert from centistokes to seconds
        if not extrapolate:
            range_check_linear(val, c, tmin, new_scale)
    else:
        keys = sorted(set(list(viscosity_scales.keys()) + list(viscosity_scales_linear.keys())))
        raise ValueError(f'Scale "{new_scale}" not recognized - allowable values are any of {keys}.')
    return float(val)


