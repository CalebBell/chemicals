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

from __future__ import division

__all__ = ['Sheffy_Johnson', 'Sato_Riedel', 'Lakshmi_Prasad', 
'Gharagheizi_liquid', 'Nicola_original', 'Nicola', 'Bahadori_liquid', 
'Mersmann_Kind_thermal_conductivity_liquid', 'DIPPR9G',
'Missenard', 'DIPPR9H', 'Filippov', 'Eucken', 'Eucken_modified', 'DIPPR9B',
'Chung', 'Eli_Hanley', 'Gharagheizi_gas', 'Bahadori_gas', 
'Stiel_Thodos_dense', 'Eli_Hanley_dense', 'Chung_dense', 'Lindsay_Bromley']

import os
import numpy as np
import pandas as pd

from fluids.numerics import horner, bisplev, implementation_optimize_tck
from fluids.constants import R, R_inv, N_A, k
from chemicals.utils import log, exp, PY37
from chemicals.data_reader import register_df_source, data_source


folder = os.path.join(os.path.dirname(__file__), 'Thermal Conductivity')


register_df_source(folder, 'Table 2-314 Vapor Thermal Conductivity of Inorganic and Organic Substances.tsv')
register_df_source(folder, 'Table 2-315 Thermal Conductivity of Inorganic and Organic Liquids.tsv')
register_df_source(folder, 'VDI PPDS Thermal conductivity of saturated liquids.tsv')
register_df_source(folder, 'VDI PPDS Thermal conductivity of gases.tsv')

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
    _load_k_data()

### Purely CSP Methods - Liquids

def Sheffy_Johnson(T, M, Tm):
    r'''Calculate the thermal conductivity of a liquid as a function of
    temperature using the Sheffy-Johnson (1961) method. Requires
    Temperature, molecular weight, and melting point.

    .. math::
        k = 1.951 \frac{1-0.00126(T-T_m)}{T_m^{0.216}MW^{0.3}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    M : float
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
    return 1.951*(1.0 - 0.00126*(T - Tm))*Tm**-0.216*M**-0.3


def Sato_Riedel(T, M, Tb, Tc):
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
    M : float
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
    return 1.1053*(3. + 20.*(1.0 - Tr)**(2.0/3.0))*M**-0.5/(3.0 + 20.0*(1.0 - Tbr)**(2.0/3.))


def Lakshmi_Prasad(T, M):
    r'''Estimates thermal conductivity of pure liquids as a function of
    temperature using a reference fluid approach. Low accuracy but quick.
    Developed using several organic fluids.

    .. math::
        \lambda = 0.0655-0.0005T + \frac{1.3855-0.00197T}{M^{0.5}}

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    M : float
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
    return 0.0655 - 0.0005*T + (1.3855 - 0.00197*T)*M**-0.5


def Gharagheizi_liquid(T, M, Tb, Pc, omega):
    r'''Estimates the thermal conductivity of a liquid as a function of
    temperature using the CSP method of Gharagheizi [1]_. A  convoluted
    method claiming high-accuracy and using only statistically significant
    variable following analalysis.

    Requires temperature, molecular weight, boiling temperature and critical
    pressure and acentric factor.

    .. math::
        &k = 10^{-4}\left[10\omega + 2P_c-2T+4+1.908(T_b+\frac{1.009B^2}{MW^2})
        +\frac{3.9287MW^4}{B^4}+\frac{A}{B^8}\right]

        &A = 3.8588MW^8(1.0045B+6.5152MW-8.9756)

        &B = 16.0407MW+2T_b-27.9074

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    M : float
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
    Pc = Pc*1E-5
    B = 16.0407*M + 2.*Tb - 27.9074
    A = 3.8588*M**8*(1.0045*B + 6.5152*M - 8.9756)
    M2 = M*M
    B_inv4 = B**-4.0
    kl = 1E-4*(10.*omega + 2.*Pc - 2.*T + 4. + 1.908*(Tb + 1.009*B*B/(M2))
        + 3.9287*M2*M2*B_inv4 + A*B_inv4*B_inv4)
    return kl


def Nicola_original(T, M, Tc, omega, Hfus):
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
    M : float
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
    Hfus = Hfus*1000
    return -0.5694 - 0.1436*Tr + 5.4893E-10*Hfus + 0.0508*omega + M**-0.0622


def Nicola(T, M, Tc, Pc, omega):
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
    M : float
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
    Tr = T/Tc
    Pc = Pc*1E-5
    return 0.5147*(-0.2537*Tr + 0.0017*Pc + 0.1501*omega + M**-0.2999)


def Bahadori_liquid(T, M):
    r'''Estimates the thermal conductivity of parafin liquid hydrocarbons.
    Fits their data well, and is useful as only MW is required.
    X is the Molecular weight, and Y the temperature.

    .. math::
        K = a + bY + CY^2 + dY^3

        a = A_1 + B_1 X + C_1 X^2 + D_1 X^3

        b = A_2 + B_2 X + C_2 X^2 + D_2 X^3

        c = A_3 + B_3 X + C_3 X^2 + D_3 X^3

        d = A_4 + B_4 X + C_4 X^2 + D_4 X^3

    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    M : float
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
    X, Y = M, T
    
    a = A[0] + X*(B[0] + X*(C[0] + D[0]*X))
    b = A[1] + X*(B[1] + X*(C[1] + D[1]*X))
    c = A[2] + X*(B[2] + X*(C[2] + D[2]*X))
    d = A[3] + X*(B[3] + X*(C[3] + D[3]*X))
    return a + Y*(b + Y*(c + d*Y))


def Mersmann_Kind_thermal_conductivity_liquid(T, MW, Tc, Vc, atoms):
    r'''Estimates the thermal conductivity of organic liquid substances
    according to the method of [1]_.

    .. math::
        \lambda^* = \frac{\lambda\cdot V_c^{2/3}\cdot T_c\cdot \text{MW}^{0.5}}
        {(k\cdot T_c)^{1.5}\cdot N_A^{7/6}}

        \lambda^* = \frac{2}{3}\left(n_a + 40\sqrt{1-T_r}\right)
        
    Parameters
    ----------
    T : float
        Temperature of the fluid [K]
    M : float
        Molecular weight of the fluid [g/mol]
    Tc : float
        Critical temperature of the fluid [K]
    Vc : float
        Critical volume of the fluid [m^3/mol]
    atoms : dict
        Dictionary of atoms and their counts, [-]

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
        
    >>> Mersmann_Kind_thermal_conductivity_liquid(400, 170.33484, 658.0, 
    ... 0.000754, {'C': 12, 'H': 26})
    0.0895271829899285

    References
    ----------
    .. [1] Mersmann, Alfons, and Matthias Kind. "Prediction of Mechanical and 
       Thermal Properties of Pure Liquids, of Critical Data, and of Vapor 
       Pressure." Industrial & Engineering Chemistry Research, January 31, 
       2017. https://doi.org/10.1021/acs.iecr.6b04323.
    '''
    na = sum(atoms.values())
    lambda_star = 2/3.*(na + 40.*(1. - T/Tc)**0.5)
    Vc = Vc*1000.0 # m^3/mol to m^3/kmol
    N_A2 = N_A*1000.0 # Their avogadro's constant is per kmol
    kl = lambda_star*(k*Tc)**1.5*N_A2**(7.0/6.)*Vc**(-2.0/3.)/Tc*MW**-0.5
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
    return kl*(0.98 + 0.0079*Pr*Tr**1.4 + 0.63*Tr**1.2*(Pr/(30. + Pr)))


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
    '''
    kl = 0.0
    for i in range(len(ws)):
        kl += ws[i]/(ks[i]*ks[i])
    return kl**-0.5


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
        \frac{\lambda M}{\eta C_v} = 1 + \frac{9/4}{C_v/R}

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
        \frac{\lambda M}{\eta C_v} = 1.32 + \frac{1.77}{C_v/R}

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
        \frac{\lambda M}{\eta C_v} = \frac{3.75 \Psi}{C_v/R}

        \Psi = 1 + \alpha \left\{[0.215+0.28288\alpha-1.061\beta+0.26665Z]/
        [0.6366+\beta Z + 1.061 \alpha \beta]\right\}

        \alpha = \frac{C_v}{R}-1.5

        \beta = 0.7862-0.7109\omega + 1.3168\omega^2

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

        Tr = \text{min}(Tr, 2)

        \theta = 1 + (\omega-0.011)\left(0.56553 - 0.86276\ln Tr - \frac{0.69852}{Tr}\right)

        \psi = [1 + (\omega - 0.011)(0.38560 - 1.1617\ln Tr)]\frac{0.288}{Z_c}

        f = \frac{T_c}{190.4}\theta

        h = \frac{V_c}{9.92E-5}\psi

        T_0 = T/f

        \eta_0^*(T_0)= \sum_{n=1}^9 C_n T_0^{(n-4)/3}

        \theta_0 = 1944 \eta_0

        \lambda^* = \lambda_0 H

        \eta^* = \eta^*_0 H \frac{MW}{16.04}

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
    2-methylbutane at low pressure, 373.15 K. Mathes calculation in [2]_.

    >>> Eli_Hanley(T=373.15, MW=72.151, Tc=460.4, Vc=3.06E-4, Zc=0.267,
    ... omega=0.227, Cvm=135.9)
    0.02247951724514062

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
          7.062481330E4, -7.116620750E3, 4.325174400E2, -1.445911210E1, 2.037119479E-1]

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
    tot = 0.0
    for i in range(9):
        tot += Cs[i]*T0_moving
        T0_moving *= T0_third
        
    eta0 = 1e-7*tot
    k0 = 1944*eta0

    H = (f*16.04/MW)**0.5*h**(-2.0/3.)
    etas = eta0*H*MW/16.04
    ks = k0*H
    return ks + etas/(MW*1e-3)*1.32*(Cvm - 1.5*R)


def Gharagheizi_gas(T, MW, Tb, Pc, omega):
    r'''Estimates the thermal conductivity of a gas as a function of
    temperature using the CSP method of Gharagheizi [1]_. A  convoluted
    method claiming high-accuracy and using only statistically significant
    variable following analalysis.

    Requires temperature, molecular weight, boiling temperature and critical
    pressure and acentric factor.

    .. math::
        k = 7.9505\times 10^{-4} + 3.989\times 10^{-5} T
        -5.419\times 10^-5 M + 3.989\times 10^{-5} A

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
    B = T + (2.*omega + 2.*T - 2.*T*(2.*omega + 3.2825)*Tb_inv + 3.2825)/(2.0*omega + T - T*(2.0*omega + 3.2825)*Tb_inv + 3.2825) - T*(2.0*omega + 3.2825)*Tb_inv
    
    x0 = (3.9752*omega + 0.1*Pc + 1.9876*B + 6.5243)
    A = (2.0*omega + T - T*(2.0*omega + 3.2825)*Tb_inv + 3.2825)/(0.1*MW*Pc*T) * x0*x0
    return 7.9505E-4 + 3.989E-5*T - 5.419E-5*MW + 3.989E-5*A


def Bahadori_gas(T, MW):
    r'''Estimates the thermal conductivity of hydrocarbons gases at low P.
    Fits their data well, and is useful as only MW is required.
    Y is the Molecular weight, and X the temperature.

    .. math::
        K = a + bY + CY^2 + dY^3

        a = A_1 + B_1 X + C_1 X^2 + D_1 X^3

        b = A_2 + B_2 X + C_2 X^2 + D_2 X^3

        c = A_3 + B_3 X + C_3 X^2 + D_3 X^3

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
    >>> Bahadori_gas(40+273.15, 20) # Point from article
    0.031968165337873326

    References
    ----------
    .. [1] Bahadori, Alireza, and Saeid Mokhatab. "Estimating Thermal
       Conductivity of Hydrocarbons." Chemical Engineering 115, no. 13
       (December 2008): 52-54
    '''
    A = [4.3931323468E-1, -3.88001122207E-2, 9.28616040136E-4, -6.57828995724E-6]
    B = [-2.9624238519E-3, 2.67956145820E-4, -6.40171884139E-6, 4.48579040207E-8]
    C = [7.54249790107E-6, -6.46636219509E-7, 1.5124510261E-8, -1.0376480449E-10]
    D = [-6.0988433456E-9, 5.20752132076E-10, -1.19425545729E-11, 8.0136464085E-14]
    X, Y = T, MW
    a = A[0] + B[0]*X + C[0]*X**2 + D[0]*X**3
    b = A[1] + B[1]*X + C[1]*X**2 + D[1]*X**3
    c = A[2] + B[2]*X + C[2]*X**2 + D[2]*X**3
    d = A[3] + B[3]*X + C[3]*X**2 + D[3]*X**3
    return a + b*Y + c*Y**2 + d*Y**3



### Thermal Conductivity of dense gases

def Stiel_Thodos_dense(T, MW, Tc, Pc, Vc, Zc, Vm, kg):
    r'''Estimates the thermal conductivity of a gas at high pressure as a
    function of temperature using difference method of Stiel and Thodos [1]_
    as shown in [2]_.

    if \rho_r < 0.5:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

    if 0.5 < \rho_r < 2.0:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

    if 2 < \rho_r < 2.8:

    .. math::
        (\lambda-\lambda^\circ)\Gamma Z_c^5=1.22\times 10^{-2} [\exp(0.535 \rho_r)-1]

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

        Vr = min(Vr, 2)

        f = \frac{T_c}{190.4}\theta

        h = \frac{V_c}{9.92E-5}\psi

        T_0 = T/f

        \rho_0 = \frac{16.04}{V}h

        \theta = 1 + (\omega-0.011)\left(0.09057 - 0.86276\ln Tr + \left(
        0.31664 - \frac{0.46568}{Tr}\right) (V_r - 0.5)\right)

        \psi = [1 + (\omega - 0.011)(0.39490(V_r - 1.02355) - 0.93281(V_r -
        0.75464)\ln T_r]\frac{0.288}{Z_c}

        \lambda_1 = 1944 \eta_0

        \lambda_2 = \left\{b_1 + b_2\left[b_3 - \ln \left(\frac{T_0}{b_4}
        \right)\right]^2\right\}\rho_0

        \lambda_3 = \exp\left(a_1 + \frac{a_2}{T_0}\right)\left\{\exp[(a_3 +
        \frac{a_4}{T_0^{1.5}})\rho_0^{0.1} + (\frac{\rho_0}{0.1617} - 1)
        \rho_0^{0.5}(a_5 + \frac{a_6}{T_0} + \frac{a_7}{T_0^2})] - 1\right\}

        \lambda^{**} = [\lambda_1 + \lambda_2 + \lambda_3]H

        H = \left(\frac{16.04}{MW}\right)^{0.5}f^{0.5}/h^{2/3}

        X = \left\{\left[1 - \frac{T}{f}\left(\frac{df}{dT}\right)_v \right]
        \frac{0.288}{Z_c}\right\}^{1.5}

        \left(\frac{df}{dT}\right)_v = \frac{T_c}{190.4}\left(\frac{d\theta}
        {d T}\right)_v

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

        \Psi = 1 + \alpha \left\{[0.215+0.28288\alpha-1.061\beta+0.26665Z]/
        [0.6366+\beta Z + 1.061 \alpha \beta]\right\}

        \alpha = \frac{C_v}{R}-1.5

        \beta = 0.7862-0.7109\omega + 1.3168\omega^2

        Z=2+10.5T_r^2

        q = 3.586\times 10^{-3} (T_c/M')^{1/2}/V_c^{2/3}

        y = \frac{V_c}{6V}

        G_1 = \frac{1-0.5y}{(1-y)^3}

        G_2 = \frac{(B_1/y)[1-\exp(-B_4y)]+ B_2G_1\exp(B_5y) + B_3G_1}
        {B_1B_4 + B_2 + B_3}

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
    ais = [2.4166E+0, -5.0924E-1, 6.6107E+0, 1.4543E+1, 7.9274E-1, -5.8634E+0, 9.1089E+1]
    bis = [7.4824E-1, -1.5094E+0, 5.6207E+0, -8.9139E+0, 8.2019E-1, 1.2801E+1, 1.2811E+2]
    cis = [-9.1858E-1, -4.9991E+1, 6.4760E+1, -5.6379E+0, -6.9369E-1, 9.5893E+0, -5.4217E+1]
    dis = [1.2172E+2, 6.9983E+1, 2.7039E+1, 7.4344E+1, 6.3173E+0, 6.5529E+1, 5.2381E+2]
    Tr = T/Tc
    mur = 131.3*dipole/(Vc*1E6*Tc)**0.5
    mur4 = mur*mur
    mur4 *= mur4 

    # From Chung Method
    alpha = Cvm/R - 1.5
    beta = 0.7862 - 0.7109*omega + 1.3168*omega**2
    Z = 2 + 10.5*(T/Tc)**2
    psi = 1 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta))

    y = Vc/(6*Vm)
    B1, B2, B3, B4, B5, B6, B7 = [ais[i] + bis[i]*omega + cis[i]*mur4 + dis[i]*association for i in range(7)]
    G1 = (1 - 0.5*y)/(1. - y)**3
    G2 = (B1/y*(1 - exp(-B4*y)) + B2*G1*exp(B5*y) + B3*G1)/(B1*B4 + B2 + B3)
    q = 3.586E-3*(Tc/(MW/1000.))**0.5/(Vc*1E6)**(2/3.)
    return 31.2*mu*psi/(MW/1000.)*(G2**-1 + B6*y) + q*B7*y**2*Tr**0.5*G2


### Thermal conductivity of gas mixtures

def Lindsay_Bromley(T, ys, ks, mus, Tbs, MWs):
    r'''Calculates thermal conductivity of a gas mixture according to
    mixing rules in [1]_ and also in [2]_.

    .. math::
        k = \sum \frac{y_i k_i}{\sum y_i A_{ij}}

        A_{ij} = \frac{1}{4} \left\{ 1 + \left[\frac{\eta_i}{\eta_j}
        \left(\frac{MW_j}{MW_i}\right)^{0.75} \left( \frac{T+S_i}{T+S_j}\right)
        \right]^{0.5} \right\}^2 \left( \frac{T+S_{ij}}{T+S_i}\right)

        S_{ij} = S_{ji} = (S_i S_j)^{0.5}
        
        S_i = 1.5 T_b

    Parameters
    ----------
    T : float
        Temperature of gas [K]
    ys : float
        Mole fractions of gas components
    ks : float
        Liquid thermal conductivites of all components, [W/m/K]
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
    0.01390264417969313

    References
    ----------
    .. [1] Lindsay, Alexander L., and LeRoy A. Bromley. "Thermal Conductivity
       of Gas Mixtures." Industrial & Engineering Chemistry 42, no. 8
       (August 1, 1950): 1508-11. doi:10.1021/ie50488a017.
    .. [2] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    cmps = range(len(ys))
    Ss = [1.5*Tb for Tb in Tbs]
    Sij = [[(Si*Sj)**0.5 for Sj in Ss] for Si in Ss]

    Aij = [[0.25*(1. + (mus[i]/mus[j]*(MWs[j]/MWs[i])**0.75
            *(T+Ss[i])/(T+Ss[j]))**0.5 )**2 *(T+Sij[i][j])/(T+Ss[i])
            for j in cmps] for i in cmps]
            
    return sum([ys[i]*ks[i]/sum(ys[j]*Aij[i][j] for j in cmps) for i in cmps])
